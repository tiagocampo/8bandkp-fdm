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
    ! ================================================================
    call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
      coo_idx, coo_cache=coo_cache)

    ! ================================================================
    ! Workspace initialization
    ! ================================================================
    if (present(ws)) then
      if (.not. ws%initialized) then
        call csr_clone_structure(blk_Q, ws%blk_Q)
        call csr_clone_structure(blk_T, ws%blk_T)
        call csr_clone_structure(blk_S, ws%blk_S)
        call csr_clone_structure(blk_SC, ws%blk_SC)
        call csr_clone_structure(blk_R, ws%blk_R)
        call csr_clone_structure(blk_RC, ws%blk_RC)
        call csr_clone_structure(blk_PZ, ws%blk_PZ)
        call csr_clone_structure(blk_PP, ws%blk_PP)
        call csr_clone_structure(blk_PM, ws%blk_PM)
        call csr_clone_structure(blk_A, ws%blk_A)
        call csr_clone_structure(blk_diff, ws%blk_diff)
        call csr_clone_structure(blk_temp, ws%blk_temp)

        ws%coo_capacity = coo_capacity
        call move_alloc(coo_rows, ws%coo_rows)
        call move_alloc(coo_cols, ws%coo_cols)
        call move_alloc(coo_vals, ws%coo_vals)
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
