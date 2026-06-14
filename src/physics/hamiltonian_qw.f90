module hamiltonian_qw

  ! ==============================================================================
  ! QW Hamiltonian builder in CSR format.
  !
  ! Builds the 8N x 8N quantum well Hamiltonian using COO assembly (same
  ! helpers as hamiltonian_wire) but with 1D kpterms(N,N,10) instead of
  ! 2D CSR kpterms arrays.  This enables sparse FEAST solver on QW geometry.
  !
  ! Key insight: the QW has only 1D spatial coupling (z-direction FD stencils).
  ! The kpterms(nz,nz,10) matrices have finite-bandwidth structure from FD
  ! stencils, but the bandwidth depends on FD order and variable coefficients.
  ! We build each kp-term CSR block by scanning the dense kpterms arrays and
  ! extracting non-zero entries.
  !
  ! The COO insertion helpers from hamiltonian_wire are grid-dimensionality-
  ! agnostic and are reused directly.
  ! ==============================================================================

  use definitions, only: IU, RQS2, SQR3, UM, ZERO, bir_pikus_blocks, &
    dp, simulation_config, wavevector
  use sparse_matrices
  use hamiltonian_blocks, only: kp_entry, get_kp_block_table, &
    KP_Q, KP_T, KP_S, KP_SC, KP_R, KP_RC, &
    KP_PP, KP_PM, KP_PZ, KP_A, KP_DIFF, KP_HALF_SUM
  use strain_solver, only: bir_pikus_blocks_free, lookup_bp_field, &
    strain_entry, get_strain_table
  use magnetic_field, only: zeeman_entry, get_zeeman_table
  use hamiltonian_wire, only: wire_coo_cache, wire_coo_cache_free, &
    insert_main_blocks, insert_profile_diagonal, insert_strain_coo, &
    insert_zeeman_coo, finalize_coo_to_csr

  implicit none

  private

  public :: qw_workspace, qw_workspace_free
  public :: ZB8bandQW_csr

  ! ------------------------------------------------------------------
  ! Workspace for QW CSR Hamiltonian k-sweep: pre-allocated CSR blocks
  ! and COO buffers.  Structure is fixed from the first call; only
  ! values are updated for subsequent k-points.
  ! ------------------------------------------------------------------
  type :: qw_workspace
    ! Pre-allocated kp-term CSR blocks
    type(csr_matrix) :: blk_Q, blk_T, blk_S, blk_SC
    type(csr_matrix) :: blk_R, blk_RC, blk_PZ, blk_PP, blk_PM, blk_A
    type(csr_matrix) :: blk_diff, blk_temp

    ! COO→CSR sort cache (populated on slow path, reused on fast path)
    type(wire_coo_cache) :: coo_cache

    ! Pre-allocated COO buffers
    integer, allocatable          :: coo_rows(:), coo_cols(:)
    complex(kind=dp), allocatable :: coo_vals(:)
    integer                       :: coo_capacity = 0

    ! COO-to-CSR sort cache
    integer, allocatable          :: coo_to_csr(:)
    integer                       :: coo_nnz_in = 0
    logical                       :: coo_cache_valid = .false.

    logical :: initialized = .false.
    logical :: was_freed = .false.
  contains
    final :: qw_workspace_finalize
  end type qw_workspace

contains

  ! ==================================================================
  ! Free the QW workspace: all CSR blocks and COO buffers.
  ! ==================================================================
  subroutine qw_workspace_free(ws)
    type(qw_workspace), intent(inout) :: ws
    if (ws%was_freed) return
    ws%was_freed = .true.
    call csr_free(ws%blk_Q)
    call csr_free(ws%blk_T)
    call csr_free(ws%blk_S)
    call csr_free(ws%blk_SC)
    call csr_free(ws%blk_R)
    call csr_free(ws%blk_RC)
    call csr_free(ws%blk_PZ)
    call csr_free(ws%blk_PP)
    call csr_free(ws%blk_PM)
    call csr_free(ws%blk_A)
    call csr_free(ws%blk_diff)
    call csr_free(ws%blk_temp)
    if (allocated(ws%coo_rows)) deallocate(ws%coo_rows)
    if (allocated(ws%coo_cols)) deallocate(ws%coo_cols)
    if (allocated(ws%coo_vals)) deallocate(ws%coo_vals)
    if (allocated(ws%coo_to_csr)) deallocate(ws%coo_to_csr)
    call wire_coo_cache_free(ws%coo_cache)
    ws%coo_nnz_in = 0
    ws%coo_cache_valid = .false.
    ws%coo_capacity = 0
    ws%initialized = .false.
  end subroutine qw_workspace_free

  ! ==================================================================
  ! Finalizer: automatically called when qw_workspace goes out of scope.
  ! ==================================================================
  subroutine qw_workspace_finalize(ws)
    type(qw_workspace), intent(inout) :: ws
    call qw_workspace_free(ws)
  end subroutine qw_workspace_finalize

  ! ==================================================================
  ! Helper: convert a complex dense matrix to COO (skipping zeros).
  ! Returns arrays suitable for csr_build_from_coo.
  ! ==================================================================
  subroutine dense_to_coo(mat, N, tol, rows, cols, vals, nnz)
    integer, intent(in) :: N
    complex(kind=dp), intent(in) :: mat(N, N)
    real(kind=dp), intent(in) :: tol
    integer, allocatable, intent(out) :: rows(:), cols(:)
    complex(kind=dp), allocatable, intent(out) :: vals(:)
    integer, intent(out) :: nnz

    integer :: ii, jj

    ! Count nonzeros first
    nnz = 0
    do jj = 1, N
      do ii = 1, N
        if (abs(mat(ii, jj)) > tol) nnz = nnz + 1
      end do
    end do

    allocate(rows(nnz), cols(nnz), vals(nnz))
    nnz = 0
    do jj = 1, N
      do ii = 1, N
        if (abs(mat(ii, jj)) > tol) then
          nnz = nnz + 1
          rows(nnz) = ii
          cols(nnz) = jj
          vals(nnz) = mat(ii, jj)
        end if
      end do
    end do
  end subroutine dense_to_coo

  ! ==================================================================
  ! Helper: build a CSR block from a dense complex matrix.
  ! ==================================================================
  subroutine dense_to_csr_block(mat, N, tol, blk)
    integer, intent(in) :: N
    complex(kind=dp), intent(in) :: mat(N, N)
    real(kind=dp), intent(in) :: tol
    type(csr_matrix), intent(out) :: blk

    integer, allocatable :: rows(:), cols(:)
    complex(kind=dp), allocatable :: vals(:)
    integer :: nnz

    call dense_to_coo(mat, N, tol, rows, cols, vals, nnz)
    if (nnz > 0) then
      call csr_build_from_coo(blk, N, N, nnz, rows, cols, vals)
    else
      call csr_init(blk, N, N)
    end if
    deallocate(rows, cols, vals)
  end subroutine dense_to_csr_block

  ! ==================================================================
  ! Compute the 10 dense kp-term blocks for a given (kx, ky).
  ! Extracted from ZB8bandQW_csr so the workspace-init path can build
  ! the cached CSR structure at a sentinel k that exposes the FULL
  ! k-independent nonzero pattern (at k=0, the k-prefactor-dependent
  ! blocks S, SC, R, RC, PP, PM collapse to zero and their structure
  ! would be lost). Values are recomputed by update_kp_term_values on
  ! the fast path; only the structure matters here.
  ! ==================================================================
  subroutine compute_dense_kp_blocks(kpterms, N, kx, ky, &
      Q, T, S, SC, R, RC, PZ, PP, PM, A)
    integer, intent(in) :: N
    real(kind=dp), intent(in) :: kpterms(N, N, 10)
    real(kind=dp), intent(in) :: kx, ky
    complex(kind=dp), intent(out) :: Q(N,N), T(N,N), S(N,N), SC(N,N)
    complex(kind=dp), intent(out) :: R(N,N), RC(N,N), PZ(N,N), PP(N,N), PM(N,N), A(N,N)

    integer :: ii, jj
    real(kind=dp) :: kx2, ky2, k2, kxky
    complex(kind=dp) :: kminus, kplus

    kx2 = kx**2; ky2 = ky**2; k2 = kx2 + ky2
    kxky = kx * ky
    kplus = kx + IU*ky; kminus = kx - IU*ky

    Q = ZERO; T = ZERO; S = ZERO; SC = ZERO
    R = ZERO; RC = ZERO; PZ = ZERO; PP = ZERO; PM = ZERO; A = ZERO

    do jj = 1, N
      do ii = 1, N
        Q(ii,jj) = -((kpterms(ii,jj,1) + kpterms(ii,jj,2))*k2 + kpterms(ii,jj,7))
        T(ii,jj) = -((kpterms(ii,jj,1) - kpterms(ii,jj,2))*k2 + kpterms(ii,jj,8))
        S(ii,jj)  =  2.0_dp * SQR3 * kminus * kpterms(ii,jj,9)
        SC(jj,ii) =  2.0_dp * SQR3 * kplus  * kpterms(ii,jj,9)
        PZ(ii,jj) = kpterms(ii,jj,6) * (-IU)
        A(ii,jj)  = cmplx(kpterms(ii,jj,5) + k2*kpterms(ii,jj,10), 0.0_dp, kind=dp)
      end do
    end do

    do ii = 1, N
      R(ii,ii)  = -SQR3 * (kpterms(ii,ii,2)*(kx2 - ky2) - 2.0_dp*IU*kpterms(ii,ii,3)*kxky)
      RC(ii,ii) = -SQR3 * (kpterms(ii,ii,2)*(kx2 - ky2) + 2.0_dp*IU*kpterms(ii,ii,3)*kxky)
      PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
      PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
    end do
  end subroutine compute_dense_kp_blocks

  ! ==================================================================
  ! Fast-path value update for one cached kp-term CSR block.
  ! Walks blk%rowptr/colind (structure fixed) and recomputes each value
  ! from kpterms(i,j,·) and the current kx,ky. Formulas mirror the
  ! dense build (ZB8bandQW_csr slow path) exactly — verified by the
  ! verify_qw_sparse_solver.py fast-vs-slow equivalence test.
  ! ==================================================================
  subroutine update_kp_term_values(blk, kpterms, kx, ky, term_id)
    type(csr_matrix), intent(inout) :: blk
    real(kind=dp), intent(in)       :: kpterms(:,:,:)
    real(kind=dp), intent(in)       :: kx, ky
    integer, intent(in)             :: term_id

    integer :: i, j, p, n
    real(kind=dp) :: kx2, ky2, k2, kxky
    complex(kind=dp) :: kplus, kminus, val

    n = blk%nrows
    kx2 = kx**2; ky2 = ky**2; k2 = kx2 + ky2; kxky = kx*ky
    kplus  = cmplx(kx,  ky, kind=dp)
    kminus = cmplx(kx, -ky, kind=dp)

    do i = 1, n
      do p = blk%rowptr(i), blk%rowptr(i+1) - 1
        j = blk%colind(p)
        select case (term_id)
        case (KP_Q)
          val = -((kpterms(i,j,1) + kpterms(i,j,2))*k2 + kpterms(i,j,7))
        case (KP_T)
          val = -((kpterms(i,j,1) - kpterms(i,j,2))*k2 + kpterms(i,j,8))
        case (KP_S)
          val = 2.0_dp*SQR3*kminus*kpterms(i,j,9)
        case (KP_SC)
          ! SC is the conjugate transpose of S: stored at (j,i) of S's formula
          val = 2.0_dp*SQR3*kplus*kpterms(j,i,9)
        case (KP_R)
          val = -SQR3*(kpterms(i,j,2)*(kx2 - ky2) - 2.0_dp*IU*kpterms(i,j,3)*kxky)
        case (KP_RC)
          val = -SQR3*(kpterms(i,j,2)*(kx2 - ky2) + 2.0_dp*IU*kpterms(i,j,3)*kxky)
        case (KP_PZ)
          val = kpterms(i,j,6)*(-IU)
        case (KP_PP)
          val = kpterms(i,j,4)*kplus*RQS2
        case (KP_PM)
          val = kpterms(i,j,4)*kminus*RQS2
        case (KP_A)
          val = cmplx(kpterms(i,j,5) + k2*kpterms(i,j,10), 0.0_dp, kind=dp)
        case default
          error stop 'update_kp_term_values: unknown term_id'
        end select
        blk%values(p) = val
      end do
    end do
  end subroutine update_kp_term_values

  ! ==================================================================
  ! Main entry point: build QW Hamiltonian in CSR format.
  !
  ! HT_csr  : output CSR matrix (8N x 8N)
  ! wv      : wavevector (kx, ky, kz)
  ! profile : band-offset profile (N, 3) [EV, EV-DSO, EC]
  ! kpterms : k.p term matrices (N, N, 10)
  ! cfg     : simulation configuration
  ! coo_cache : optional COO sort cache for repeated builds
  ! ws      : optional workspace for fast-path reuse
  ! ==================================================================
  subroutine ZB8bandQW_csr(HT_csr, wv, profile, kpterms, cfg, coo_cache, ws)

    type(csr_matrix), intent(inout)         :: HT_csr
    type(wavevector), intent(in)            :: wv
    real(kind=dp), intent(in), contiguous   :: profile(:,:)
    real(kind=dp), intent(in), contiguous   :: kpterms(:,:,:)
    type(simulation_config), intent(in)     :: cfg
    type(wire_coo_cache), intent(inout), optional :: coo_cache
    type(qw_workspace), intent(inout), optional :: ws

    integer :: N, Ntot, ii, jj
    real(kind=dp) :: kx, ky, kx2, ky2, k2, kxky
    complex(kind=dp) :: kminus, kplus

    ! Dense kp-term blocks (N x N complex)
    complex(kind=dp), allocatable :: Q(:,:), T(:,:), S(:,:), SC(:,:)
    complex(kind=dp), allocatable :: R(:,:), RC(:,:), PZ(:,:), PP(:,:), PM(:,:), A(:,:)

    ! kp-term CSR blocks
    type(csr_matrix) :: blk_Q, blk_T, blk_S, blk_SC
    type(csr_matrix) :: blk_R, blk_RC, blk_PZ, blk_PP, blk_PM, blk_A
    type(csr_matrix) :: blk_diff, blk_temp

    ! COO assembly arrays
    integer, allocatable :: coo_rows(:), coo_cols(:)
    complex(kind=dp), allocatable :: coo_vals(:)
    integer :: coo_idx, coo_capacity, nnz_est

    real(kind=dp), parameter :: sparse_tol = 1.0e-300_dp

    N = size(kpterms, 1)
    Ntot = 8 * N

    kx = wv%kx; ky = wv%ky
    kx2 = kx**2; ky2 = ky**2; k2 = kx2 + ky2
    kxky = kx * ky
    kplus = kx + IU*ky; kminus = kx - IU*ky

    ! ================================================================
    ! FAST PATH (subsequent k-points): structure is fixed from the
    ! first (slow) call. Update kp-term block values in place, re-scatter
    ! into the cached COO buffers, and rebuild HT_csr via the cached
    ! COO→CSR sort map. No dense matrices, no dense_to_csr_block O(N^2)
    ! scans.
    !
    ! COO scatter order is byte-identical to the slow path (driven by
    ! the unchanged block sparsity), so the cached coo_to_csr sort map
    ! stays valid across k-points. The map is reused every iteration
    ! (the real O(NNZ) vs O(NNZ log NNZ) win); the rowptr/colind are
    ! rebuilt unconditionally because the active QW FEAST sweep frees
    ! HT_csr each iteration (HT_csr%nnz == 0 on entry). A value-only
    ! path is intentionally omitted: it would be unsafe for a caller
    ! reusing a persistent HT_csr across a Γ→off-Γ transition (Γ-point
    ! structure is sparser than the sentinel-recorded map).
    ! ================================================================
    if (present(ws)) then
      if (ws%initialized) then
        ! Update the 10 cached kp-term blocks in place (O(NNZ) each)
        call update_kp_term_values(ws%blk_Q,  kpterms, kx, ky, KP_Q)
        call update_kp_term_values(ws%blk_T,  kpterms, kx, ky, KP_T)
        call update_kp_term_values(ws%blk_S,  kpterms, kx, ky, KP_S)
        call update_kp_term_values(ws%blk_SC, kpterms, kx, ky, KP_SC)
        call update_kp_term_values(ws%blk_R,  kpterms, kx, ky, KP_R)
        call update_kp_term_values(ws%blk_RC, kpterms, kx, ky, KP_RC)
        call update_kp_term_values(ws%blk_PZ, kpterms, kx, ky, KP_PZ)
        call update_kp_term_values(ws%blk_PP, kpterms, kx, ky, KP_PP)
        call update_kp_term_values(ws%blk_PM, kpterms, kx, ky, KP_PM)
        call update_kp_term_values(ws%blk_A,  kpterms, kx, ky, KP_A)

        ! blk_diff = Q - T, blk_temp = 0.5*(Q + T).  csr_add reallocates
        ! its output (intent(out)); ws%blk_diff/blk_temp are distinct
        ! from blk_Q/blk_T so no aliasing.  Two tiny tridiagonal blocks
        ! per call — negligible vs the 10 dense matrices eliminated.
        call csr_add(ws%blk_Q, ws%blk_T, ws%blk_diff, UM, &
                     cmplx(-1.0_dp, 0.0_dp, kind=dp))
        call csr_add(ws%blk_Q, ws%blk_T, ws%blk_temp, &
                     cmplx(0.5_dp, 0.0_dp, kind=dp), &
                     cmplx(0.5_dp, 0.0_dp, kind=dp))

        ! Re-scatter into cached COO buffers (reset index, overwrite in
        ! place; same insertion order as the slow path).
        coo_idx = 0
        call insert_main_blocks(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
          ws%coo_capacity, coo_idx, ws%blk_Q, ws%blk_T, ws%blk_S, ws%blk_SC, &
          ws%blk_R, ws%blk_RC, ws%blk_PZ, ws%blk_PP, ws%blk_PM, ws%blk_A, &
          ws%blk_diff, ws%blk_temp, N)
        call insert_profile_diagonal(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
          ws%coo_capacity, coo_idx, profile, N)
        if (allocated(cfg%strain_blocks%delta_Ec)) then
          call insert_strain_coo(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
            ws%coo_capacity, coo_idx, cfg%strain_blocks, N)
        end if
        if (any(abs(cfg%bdg%B_vec) > 1.0e-12_dp)) then
          call insert_zeeman_coo(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
            ws%coo_capacity, coo_idx, cfg%bdg%B_vec, cfg%bdg%g_factor, &
            cfg%grid, N)
        end if

        ! Rebuild HT_csr from the re-scattered COO via the cached COO→CSR
        ! sort map (single rebuild path). The sort map is recorded once
        ! (slow path / workspace init) and reused every iteration: this is
        ! the O(NNZ) win over a fresh O(NNZ log NNZ) sort. We rebuild the
        ! rowptr/colind unconditionally rather than value-only, because the
        ! active QW FEAST sweep frees HT_csr each iteration (HT_csr%nnz == 0
        ! on entry); a value-only path would be a latent footgun for any
        ! future caller reusing a persistent HT_csr across a Γ→off-Γ
        ! transition (Γ-point structure is sparser than the sentinel-recorded
        ! map). The map is deterministic so repopulation is harmless.
        call csr_build_from_coo_cached(HT_csr, Ntot, Ntot, coo_idx, &
          ws%coo_rows(1:coo_idx), ws%coo_cols(1:coo_idx), &
          ws%coo_vals(1:coo_idx), ws%coo_cache%coo_to_csr)
        ws%coo_cache%coo_nnz_in = coo_idx
        ws%coo_cache%initialized = .true.
        return
      end if
    end if

    ! ================================================================
    ! Build dense kp-term matrices (exact same formulas as ZB8bandQW)
    ! ================================================================
    allocate(Q(N,N), T(N,N), S(N,N), SC(N,N))
    allocate(R(N,N), RC(N,N), PZ(N,N), PP(N,N), PM(N,N), A(N,N))
    Q = ZERO; T = ZERO; S = ZERO; SC = ZERO
    R = ZERO; RC = ZERO; PZ = ZERO; PP = ZERO; PM = ZERO; A = ZERO

    do jj = 1, N
      do ii = 1, N
        Q(ii,jj) = -((kpterms(ii,jj,1) + kpterms(ii,jj,2))*k2 + kpterms(ii,jj,7))
        T(ii,jj) = -((kpterms(ii,jj,1) - kpterms(ii,jj,2))*k2 + kpterms(ii,jj,8))
        S(ii,jj)  =  2.0_dp * SQR3 * kminus * kpterms(ii,jj,9)
        SC(jj,ii) =  2.0_dp * SQR3 * kplus  * kpterms(ii,jj,9)
        PZ(ii,jj) = kpterms(ii,jj,6) * (-IU)
        A(ii,jj)  = cmplx(kpterms(ii,jj,5) + k2*kpterms(ii,jj,10), 0.0_dp, kind=dp)
      end do
    end do

    do ii = 1, N
      R(ii,ii)  = -SQR3 * (kpterms(ii,ii,2)*(kx2 - ky2) - 2.0_dp*IU*kpterms(ii,ii,3)*kxky)
      RC(ii,ii) = -SQR3 * (kpterms(ii,ii,2)*(kx2 - ky2) + 2.0_dp*IU*kpterms(ii,ii,3)*kxky)
      PP(ii,ii) = kpterms(ii,ii,4) * kplus  * RQS2
      PM(ii,ii) = kpterms(ii,ii,4) * kminus * RQS2
    end do

    ! ================================================================
    ! Convert dense kp-term blocks to CSR
    ! ================================================================
    call dense_to_csr_block(Q, N, sparse_tol, blk_Q)
    call dense_to_csr_block(T, N, sparse_tol, blk_T)
    call dense_to_csr_block(S, N, sparse_tol, blk_S)
    call dense_to_csr_block(SC, N, sparse_tol, blk_SC)
    call dense_to_csr_block(R, N, sparse_tol, blk_R)
    call dense_to_csr_block(RC, N, sparse_tol, blk_RC)
    call dense_to_csr_block(PZ, N, sparse_tol, blk_PZ)
    call dense_to_csr_block(PP, N, sparse_tol, blk_PP)
    call dense_to_csr_block(PM, N, sparse_tol, blk_PM)
    call dense_to_csr_block(A, N, sparse_tol, blk_A)

    deallocate(Q, T, S, SC, R, RC, PZ, PP, PM, A)

    ! blk_diff = Q - T, blk_temp = 0.5*(Q + T)
    call csr_add(blk_Q, blk_T, blk_diff, UM, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    call csr_add(blk_Q, blk_T, blk_temp, cmplx(0.5_dp, 0.0_dp, kind=dp), &
      cmplx(0.5_dp, 0.0_dp, kind=dp))

    ! ================================================================
    ! COO capacity estimation
    ! ================================================================
    nnz_est = 2*blk_Q%nnz + 2*blk_T%nnz + 6*blk_S%nnz + 6*blk_SC%nnz
    nnz_est = nnz_est + 5*blk_R%nnz + 3*blk_RC%nnz + 8*blk_PZ%nnz
    nnz_est = nnz_est + 6*blk_PP%nnz + 5*blk_PM%nnz + 2*blk_A%nnz
    nnz_est = nnz_est + 4*(blk_Q%nnz + blk_T%nnz)
    nnz_est = nnz_est + 2*(blk_Q%nnz + blk_T%nnz)
    nnz_est = nnz_est + 8 * N
    if (allocated(cfg%strain_blocks%delta_Ec)) then
      nnz_est = nnz_est + 32 * N
    end if
    if (any(abs(cfg%bdg%B_vec) > 1.0e-12_dp)) then
      nnz_est = nnz_est + 8 * N
    end if
    coo_capacity = nnz_est + nnz_est / 5

    allocate(coo_rows(coo_capacity))
    allocate(coo_cols(coo_capacity))
    allocate(coo_vals(coo_capacity))
    coo_idx = 0

    ! ================================================================
    ! Insert 8x8 blocks into COO arrays
    ! ================================================================
    call insert_main_blocks(coo_rows, coo_cols, coo_vals, coo_capacity, &
      coo_idx, blk_Q, blk_T, blk_S, blk_SC, blk_R, blk_RC, blk_PZ, &
      blk_PP, blk_PM, blk_A, blk_diff, blk_temp, N)

    ! Band-offset profile
    call insert_profile_diagonal(coo_rows, coo_cols, coo_vals, coo_capacity, &
      coo_idx, profile, N)

    ! Bir-Pikus strain corrections
    if (allocated(cfg%strain_blocks%delta_Ec)) then
      call insert_strain_coo(coo_rows, coo_cols, coo_vals, coo_capacity, &
        coo_idx, cfg%strain_blocks, N)
    end if

    ! Zeeman splitting
    if (any(abs(cfg%bdg%B_vec) > 1.0e-12_dp)) then
      call insert_zeeman_coo(coo_rows, coo_cols, coo_vals, coo_capacity, &
        coo_idx, cfg%bdg%B_vec, cfg%bdg%g_factor, cfg%grid, N)
    end if

    ! ================================================================
    ! Build final CSR from COO
    ! Route through the workspace's COO→CSR sort cache when ws is present
    ! so the sort map is recorded on the first (slow) call and reused on
    ! subsequent calls. The values-only path (csr_set_values_from_coo) is
    ! only safe when HT_csr already holds the matching structure; callers
    ! that free HT_csr between calls (HT_csr%nnz == 0) must rebuild from
    ! scratch. This mirrors the wire path's HT_csr%nnz > 0 guard.
    ! ================================================================
    if (present(ws)) then
      if (ws%coo_cache%initialized .and. HT_csr%nnz > 0) then
        call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
          coo_idx, coo_cache=ws%coo_cache)
      else if (.not. ws%coo_cache%initialized) then
        call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
          coo_idx, coo_cache=ws%coo_cache)
      else
        ! Cache initialized but HT_csr was freed between calls: rebuild
        ! from scratch (the cached map stays valid for future reuse).
        call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
          coo_idx)
      end if
    else if (present(coo_cache)) then
      call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
        coo_idx, coo_cache=coo_cache)
    else
      call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
        coo_idx)
    end if

    ! ================================================================
    ! Workspace initialization: capture the k-INDEPENDENT CSR structure
    ! for all cached blocks. Built at a sentinel k=(1,0) so the
    ! k-prefactor-dependent blocks (S, SC, R, RC, PP, PM) expose their
    ! full nonzero pattern — at k=0 these blocks are identically zero
    ! and dense_to_csr_block would collapse them to empty structure,
    ! silently dropping their contributions at k!=0 on the fast path.
    ! Values are recomputed by update_kp_term_values; only structure
    ! matters here.  blk_diff/blk_temp derive their structure from the
    ! sentinel-k Q and T (same union pattern as any real k).
    ! ================================================================
    if (present(ws)) then
      if (.not. ws%initialized) then
        block
          ! Sentinel-k dense temporaries (structure capture only)
          complex(kind=dp), allocatable :: sQ(:,:), sT(:,:), sS(:,:), sSC(:,:)
          complex(kind=dp), allocatable :: sR(:,:), sRC(:,:), sPZ(:,:), sPP(:,:), sPM(:,:), sA(:,:)
          type(csr_matrix) :: sblk_Q, sblk_T, sblk_S, sblk_SC
          type(csr_matrix) :: sblk_R, sblk_RC, sblk_PZ, sblk_PP, sblk_PM, sblk_A
          type(csr_matrix) :: sblk_diff, sblk_temp

          ! kx=1, ky=0: makes k2=1, kplus=1, kminus=1, kx2-ky2=1, kxky=0.
          ! All k-prefactors are nonzero, so every structurally-nonzero
          ! entry of each block is exposed (no sparsity collapse).
          allocate(sQ(N,N), sT(N,N), sS(N,N), sSC(N,N))
          allocate(sR(N,N), sRC(N,N), sPZ(N,N), sPP(N,N), sPM(N,N), sA(N,N))
          call compute_dense_kp_blocks(kpterms, N, 1.0_dp, 0.0_dp, &
            sQ, sT, sS, sSC, sR, sRC, sPZ, sPP, sPM, sA)

          call dense_to_csr_block(sQ, N, sparse_tol, sblk_Q)
          call dense_to_csr_block(sT, N, sparse_tol, sblk_T)
          call dense_to_csr_block(sS, N, sparse_tol, sblk_S)
          call dense_to_csr_block(sSC, N, sparse_tol, sblk_SC)
          call dense_to_csr_block(sR, N, sparse_tol, sblk_R)
          call dense_to_csr_block(sRC, N, sparse_tol, sblk_RC)
          call dense_to_csr_block(sPZ, N, sparse_tol, sblk_PZ)
          call dense_to_csr_block(sPP, N, sparse_tol, sblk_PP)
          call dense_to_csr_block(sPM, N, sparse_tol, sblk_PM)
          call dense_to_csr_block(sA, N, sparse_tol, sblk_A)
          deallocate(sQ, sT, sS, sSC, sR, sRC, sPZ, sPP, sPM, sA)

          ! Derived blocks: structure = union of Q and T patterns.
          call csr_add(sblk_Q, sblk_T, sblk_diff, UM, &
                       cmplx(-1.0_dp, 0.0_dp, kind=dp))
          call csr_add(sblk_Q, sblk_T, sblk_temp, &
                       cmplx(0.5_dp, 0.0_dp, kind=dp), &
                       cmplx(0.5_dp, 0.0_dp, kind=dp))

          call csr_clone_structure(sblk_Q, ws%blk_Q)
          call csr_clone_structure(sblk_T, ws%blk_T)
          call csr_clone_structure(sblk_S, ws%blk_S)
          call csr_clone_structure(sblk_SC, ws%blk_SC)
          call csr_clone_structure(sblk_R, ws%blk_R)
          call csr_clone_structure(sblk_RC, ws%blk_RC)
          call csr_clone_structure(sblk_PZ, ws%blk_PZ)
          call csr_clone_structure(sblk_PP, ws%blk_PP)
          call csr_clone_structure(sblk_PM, ws%blk_PM)
          call csr_clone_structure(sblk_A, ws%blk_A)
          call csr_clone_structure(sblk_diff, ws%blk_diff)
          call csr_clone_structure(sblk_temp, ws%blk_temp)

          call csr_free(sblk_Q); call csr_free(sblk_T)
          call csr_free(sblk_S); call csr_free(sblk_SC)
          call csr_free(sblk_R); call csr_free(sblk_RC)
          call csr_free(sblk_PZ); call csr_free(sblk_PP)
          call csr_free(sblk_PM); call csr_free(sblk_A)
          call csr_free(sblk_diff); call csr_free(sblk_temp)
        end block

        ! Size the cached COO buffers from the SENTINEL-k block nnz
        ! (>= the real-k=0 nnz, since sentinel exposes the full k-
        ! independent sparsity).  The slow path's local coo_rows/cols/
        ! vals were sized for the k=0 pattern and may be too small for
        ! the fast path's denser scatter, so we discard them and
        ! allocate fresh at the sentinel capacity.
        nnz_est = 2*ws%blk_Q%nnz + 2*ws%blk_T%nnz + 6*ws%blk_S%nnz + 6*ws%blk_SC%nnz
        nnz_est = nnz_est + 5*ws%blk_R%nnz + 3*ws%blk_RC%nnz + 8*ws%blk_PZ%nnz
        nnz_est = nnz_est + 6*ws%blk_PP%nnz + 5*ws%blk_PM%nnz + 2*ws%blk_A%nnz
        nnz_est = nnz_est + 4*(ws%blk_Q%nnz + ws%blk_T%nnz)
        nnz_est = nnz_est + 2*(ws%blk_Q%nnz + ws%blk_T%nnz)
        nnz_est = nnz_est + 8 * N
        if (allocated(cfg%strain_blocks%delta_Ec)) then
          nnz_est = nnz_est + 32 * N
        end if
        if (any(abs(cfg%bdg%B_vec) > 1.0e-12_dp)) then
          nnz_est = nnz_est + 8 * N
        end if
        ws%coo_capacity = nnz_est + nnz_est / 5
        allocate(ws%coo_rows(ws%coo_capacity))
        allocate(ws%coo_cols(ws%coo_capacity))
        allocate(ws%coo_vals(ws%coo_capacity))

        ! Record the COO→CSR sort map from the sentinel-k structure so
        ! the fast path's re-scatter (same sentinel structure → same
        ! (row,col) sequence) can reuse it via csr_build_from_coo_cached
        ! without re-sorting.  The cached blocks currently hold zero
        ! values (from csr_clone_structure); only the (row,col) pattern
        ! drives the map, so values are irrelevant here.  We scatter
        ! into the just-allocated ws COO buffers and build a throwaway
        ! CSR to capture the map.
        block
          type(csr_matrix) :: map_scratch
          integer :: map_idx

          map_idx = 0
          call insert_main_blocks(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
            ws%coo_capacity, map_idx, ws%blk_Q, ws%blk_T, ws%blk_S, ws%blk_SC, &
            ws%blk_R, ws%blk_RC, ws%blk_PZ, ws%blk_PP, ws%blk_PM, ws%blk_A, &
            ws%blk_diff, ws%blk_temp, N)
          call insert_profile_diagonal(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
            ws%coo_capacity, map_idx, profile, N)
          if (allocated(cfg%strain_blocks%delta_Ec)) then
            call insert_strain_coo(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
              ws%coo_capacity, map_idx, cfg%strain_blocks, N)
          end if
          if (any(abs(cfg%bdg%B_vec) > 1.0e-12_dp)) then
            call insert_zeeman_coo(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
              ws%coo_capacity, map_idx, cfg%bdg%B_vec, cfg%bdg%g_factor, &
              cfg%grid, N)
          end if

          call csr_build_from_coo_cached(map_scratch, Ntot, Ntot, map_idx, &
            ws%coo_rows(1:map_idx), ws%coo_cols(1:map_idx), &
            ws%coo_vals(1:map_idx), ws%coo_cache%coo_to_csr)
          ws%coo_cache%coo_nnz_in = map_idx
          ws%coo_cache%initialized = .true.
          call csr_free(map_scratch)
        end block

        ws%initialized = .true.
      end if
    end if

    ! Cleanup
    call csr_free(blk_Q)
    call csr_free(blk_T)
    call csr_free(blk_S)
    call csr_free(blk_SC)
    call csr_free(blk_R)
    call csr_free(blk_RC)
    call csr_free(blk_PZ)
    call csr_free(blk_PP)
    call csr_free(blk_PM)
    call csr_free(blk_A)
    call csr_free(blk_diff)
    call csr_free(blk_temp)
    if (allocated(coo_rows)) deallocate(coo_rows)
    if (allocated(coo_cols)) deallocate(coo_cols)
    if (allocated(coo_vals)) deallocate(coo_vals)

  end subroutine ZB8bandQW_csr

end module hamiltonian_qw
