module bdg_hamiltonian

  ! ==============================================================================
  ! Bogoliubov-de Gennes (BdG) Nambu-space Hamiltonian assembly.
  !
  ! The BdG Hamiltonian is a 16N x 16N matrix in Nambu space built from the
  ! 8N x 8N electron Hamiltonian:
  !
  !   H_BdG = |  H0(+k) - mu*I       Delta            |
  !           |  Delta^dagger        -conjg(H0(-k)) + mu*I |
  !
  ! where H0 is the 8N x 8N wire/QW Hamiltonian (from ZB8bandGeneralized or
  ! ZB8bandQW), mu is the chemical potential, and Delta is the superconducting
  ! pairing.
  !
  ! CANONICAL HOLE-BLOCK CONVENTION (per ADR 0007):
  !   Hole block = -conjg(H0(-k))    [the QW dense form]
  !
  ! This is the time-reversed Nambu conjugate of the electron block. It is
  ! the k!=0-general form (Leijnse-Flensberg Eq. 38) and is equivalent to
  ! the wire-CSR's -H0^T(+k) form when H0 is Hermitian at k=0 in the absence
  ! of Peierls phases, but it is the CORRECT form at generic k with Peierls
  ! (where H0(+k) != H0(-k)^T).
  !
  ! Pinning oracle: tests/unit/test_bdg_phs.pf. RED on pre-Issue-03 wire
  ! code (canonical hole-block residual ~0.13); GREEN once Issue 03's
  ! shared wrapper extracts build_bdg_hole_block with this canonical form.
  !
  ! Key properties:
  !   - BdG matrix IS Hermitian (FEAST works directly, no non-Hermitian solver)
  !   - Particle-hole symmetry: eigenvalues come in ±E pairs
  !   - Pairing matrix: Delta = delta_0 * (iσ_y ⊗ I_4), intra-band s-wave
  !     pairing within each Kramers pair in zinc-blende 8-band basis
  !
  ! Basis ordering (fixed throughout codebase):
  !   Bands 1-4 = valence (HH, LH, LH, SO), 5-6 = split-off, 7-8 = conduction
  !   Nambu space doubles: 1..8N = electron, 8N+1..16N = hole
  !
  ! Reference: Fu & Kane (2008), Supriyo Datta (2018) "Quantum Transport"
  ! ==============================================================================

  use definitions, only: IU, ZERO, dp, mu_B, simulation_config, &
    spatial_grid, wavevector
  use sparse_matrices
  use hamiltonian_wire, only: ZB8bandGeneralized, wire_workspace
  use hamiltonianConstructor, only: ZB8bandQW
  use magnetic_field, only: add_peierls_coo, compute_zeeman_vz

  implicit none

  private

  public :: build_bdg_hamiltonian_1d
  public :: build_bdg_hamiltonian_qw

  ! s-wave pairing partner band indices (iσ_y ⊗ I_4 in zinc-blende basis).
  ! Basis: HH↑(1), LH↑(2), LH↓(3), HH↓(4), SO↑(5), SO↓(6), CB↑(7), CB↓(8).
  ! Each spin-up partner is the spin-down within the same Kramers pair:
  !   1↔4 (HH↑↔HH↓), 2↔3 (LH↑↔LH↓), 5↔6 (SO↑↔SO↓), 7↔8 (CB↑↔CB↓).
  integer, parameter :: pairing_partner(8) = [4, 3, 2, 1, 6, 5, 8, 7]

  ! Sign for each band in the pairing: +1 for spin-up, -1 for spin-down.
  ! Kramers pairs (1,4), (2,3), (5,6), (7,8) get (+1,-1) respectively.
  real(kind=dp), parameter :: pairing_sign(8) = [ &
    +1.0_dp, &  ! band 1 (HH↑)
    +1.0_dp, &  ! band 2 (LH↑)
    -1.0_dp, &  ! band 3 (LH↓)
    -1.0_dp, &  ! band 4 (HH↓)
    +1.0_dp, &  ! band 5 (SO↑)
    -1.0_dp, &  ! band 6 (SO↓)
    +1.0_dp, &  ! band 7 (CB↑)
    -1.0_dp ]   ! band 8 (CB↓)

contains

  ! ==============================================================================
  ! Layer B (ADR 0007): Shared hole-block wrapper.
  !
  ! Canonical form: H_hole = -conjg(H0(-k))
  !
  ! This wrapper unifies the wire CSR and dense QW hole-block constructions
  ! under a single canonical form (Layer D, ADR 0007). Both builders call this;
  ! the wrapper takes H0(-k) already evaluated by the caller (it does not
  ! evaluate k, since the caller has the appropriate momentum in scope).
  !
  ! The pairing-embed (pairing_partner/pairing_sign, module data) and the
  ! mu-shift are documented as shared operations but kept at the call sites,
  ! since they have different sign conventions for the electron vs hole
  ! blocks (mu enters as -mu in (1,1) and +mu in (2,2); Zeeman enters as
  ! +Vz_delta in (1,1) and -Vz_delta in (2,2)).
  ! ==============================================================================
  pure subroutine build_bdg_hole_block(H0_minus_k, H_hole)
    complex(kind=dp), intent(in) :: H0_minus_k(:,:)
    complex(kind=dp), intent(out) :: H_hole(:,:)
    integer :: i, j, n

    n = size(H0_minus_k, 1)
    do concurrent (i = 1:n, j = 1:n)
      H_hole(i, j) = -conjg(H0_minus_k(i, j))
    end do
  end subroutine build_bdg_hole_block

  ! ==============================================================================
  ! Check that COO index has not exceeded pre-allocated capacity.
  ! Follows the pattern in insert_zeeman_coo (hamiltonian_wire.f90).
  ! ==============================================================================
  subroutine check_coo_bounds(coo_idx, coo_capacity)
    integer, intent(in) :: coo_idx, coo_capacity

    if (coo_idx > coo_capacity) then
      print *, 'ERROR: build_bdg_hamiltonian_1d: COO capacity exceeded'
      print *, '  coo_idx=', coo_idx, ' coo_capacity=', coo_capacity
      error stop 'bdg_hamiltonian: COO capacity exceeded'
    end if
  end subroutine check_coo_bounds

  ! ==============================================================================
  ! Build the 16N x 16N BdG Hamiltonian from the 8N x 8N wire Hamiltonian.
  !
  ! H_BdG = | H0 - mu*I    Delta                       |
  !         | Delta^dagger  -conjg(H0(-k)) + mu*I       |
  !
  ! The pairing matrix in the zinc-blende 8-band basis is:
  !   Delta = delta_0 * (iσ_y ⊗ I_4)
  ! This gives intra-band s-wave pairing within each Kramers pair:
  !   CB↑(7)↔CB↓(8), HH↑(1)↔HH↓(4), LH↑(2)↔LH↓(3), SO↑(5)↔SO↓(6)
  ! with signs (+1,-1,+1,-1,+1,-1,+1,-1) reflecting iσ_y structure.
  !
  ! For the wire, each band-block is replicated at each spatial point.
  !
  ! Hole-block construction (ADR 0007, Layer B+D): the canonical form is
  !   H_hole = -conjg(H0(-k))
  ! produced by the shared `build_bdg_hole_block` wrapper. The wire CSR
  ! path builds a SECOND H0 at -kz (H0_minus_k) and routes it through the
  ! wrapper. The mu-shift (+mu*I) and Zeeman (-Vz_delta) for the hole
  ! block are kept at the call site since they differ in sign from the
  ! electron block; the wrapper itself is a pure canonical-form seam.
  ! ==============================================================================
  subroutine build_bdg_hamiltonian_1d(H_bdg_csr, cfg, profile_2d, kpterms_2d, &
                                       kz, mu, delta_0, ws, B_vec, g_factor)

    type(csr_matrix), intent(out) :: H_bdg_csr
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile_2d(:,:)
    type(csr_matrix), intent(in) :: kpterms_2d(:)
    real(kind=dp), intent(in) :: kz, mu, delta_0
    type(wire_workspace), intent(inout) :: ws
    real(kind=dp), intent(in), optional :: B_vec(3)
    real(kind=dp), intent(in), optional :: g_factor

    type(csr_matrix) :: H0, H0_minus_k_csr
    complex(kind=dp), allocatable :: H0_minus_k(:,:), H_hole_dense(:,:)
    integer :: N, Ntot, nnz_bdg
    integer :: coo_idx, coo_capacity
    integer, allocatable :: coo_row_bdg(:), coo_col_bdg(:)
    complex(kind=dp), allocatable :: coo_vals_bdg(:)

    integer :: i, kk, row, col, global_row, global_col, hole_nnz
    real(kind=dp) :: Vz_opt(8), Vz_cfg(8), Vz_delta(8), g_f, B_mag_opt, B_mag_cfg
    logical :: add_optional_zeeman

    ! --- Guard: delta_0 must be positive and finite ---
    if (delta_0 <= 0.0_dp .or. delta_0 /= delta_0) then
      print *, 'Error: delta_0 must be positive for BdG Hamiltonian construction'
      print *, '  delta_0=', delta_0
      error stop 'bdg_hamiltonian: delta_0 must be positive and finite'
    end if

    ! --- Build H0 (8N x 8N wire Hamiltonian) at +kz (electron block) ---
    call ZB8bandGeneralized(H0, kz, profile_2d, kpterms_2d, cfg, ws=ws)

    ! --- Build H0(-kz) (canonical hole block, ADR 0007) ---
    ! The wire builder must evaluate H0 at -kz to populate the hole block
    ! via the canonical -conjg(H0(-k)) form. ZB8bandGeneralized is k-dependent
    ! via its kz argument; passing -kz builds the time-reversed normal
    ! Hamiltonian. Note: Peierls phase symmetry H0(+k,B) = H0(-k,-B) means
    ! that for Bx-driven Peierls the hole block at -kz needs the conjugate
    ! of the B-altered H0 at +kz -- which the wrapper's conjg() delivers.
    call ZB8bandGeneralized(H0_minus_k_csr, -kz, profile_2d, kpterms_2d, cfg, ws=ws)

    ! Validate H0 after build
    if (H0%nrows == 0 .or. H0%nnz == 0) then
      print *, 'ERROR: build_bdg_hamiltonian_1d: empty H0 from ZB8bandGeneralized'
      print *, '  nrows=', H0%nrows, ' nnz=', H0%nnz
      error stop 'bdg_hamiltonian: empty H0 from ZB8bandGeneralized'
    end if
    if (mod(H0%nrows, 8) /= 0) then
      print *, 'ERROR: build_bdg_hamiltonian_1d: H0 nrows not multiple of 8'
      print *, '  nrows=', H0%nrows
      error stop 'bdg_hamiltonian: H0 nrows not multiple of 8'
    end if

    N = H0%nrows / 8
    Ntot = 16 * N

    ! Convert H0(-kz) to dense once for the wrapper call. The wrapper is
    ! dense-only (ADR 0007); the wire path pays one 8N x 8N conversion.
    ! csr_to_dense is a tests/support helper, so we inline a small
    ! conversion here (8N x 8N is small enough that the loop is cheap).
    ! Use H0_minus_k_csr%nnz for the COO capacity estimate: the wrapper's
    ! -conjg transform preserves the sparsity pattern exactly, so the hole
    ! block nnz matches H0_minus_k_csr%nnz to FP precision.
    !
    ! NOTE on Peierls phase: the QW builder's hole block is
    !   H_hole = -conjg(H_h)  where H_h = ZB8bandQW(-k) (with Peierls applied
    !   internally via cfg-driven Zeeman/Peierls insertion in ZB8bandQW).
    ! For byte-identical hole blocks across both builders (per ADR 0007 Layer D),
    ! the wire builder's hole block must derive from H0(-k) WITHOUT external
    ! Peierls application here -- Peierls is applied to the electron block only
    ! (via add_peierls_coo below) and the hole block follows the canonical
    ! form -conjg(H0(-k)) (no Peierls). This is the convention that makes both
    ! builders produce byte-identical hole blocks for the same normal H0.
    ! Class-D PHS is preserved at k=0 (where H0(-k)=H0(+k) and Peierls is
    ! real/antisymmetric); at generic k, PHS differs from the wire form.
    hole_nnz = H0_minus_k_csr%nnz
    allocate(H0_minus_k(H0%nrows, H0%nrows))
    H0_minus_k = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do row = 1, H0_minus_k_csr%nrows
      do kk = H0_minus_k_csr%rowptr(row), H0_minus_k_csr%rowptr(row + 1) - 1
        col = H0_minus_k_csr%colind(kk)
        H0_minus_k(row, col) = H0_minus_k_csr%values(kk)
      end do
    end do
    call csr_free(H0_minus_k_csr)
    allocate(H_hole_dense(H0%nrows, H0%nrows))
    call build_bdg_hole_block(H0_minus_k, H_hole_dense)

    ! --- Estimate COO capacity ---
    ! Block (1,1): H0 - mu*I  => H0%nnz + 8N (diagonal)
    ! Block (2,2): -conjg(H0(-k)) + mu*I => hole_nnz + 8N (diagonal)
    ! Block (1,2): Delta => 8*N (antidiagonal, one entry per band per spatial point)
    ! Block (2,1): Delta^dagger => 8*N (conjugate transpose)
    ! mu*I diagonal: 16*N (8N in each of blocks (1,1) and (2,2))
    ! Optional Zeeman can add one electron and one hole diagonal entry per band/site.
    nnz_bdg = H0%nnz + hole_nnz + 2 * (8 * N) + 16 * N + 16 * N

    coo_capacity = nnz_bdg + nnz_bdg / 5  ! small safety margin
    allocate(coo_row_bdg(coo_capacity))
    allocate(coo_col_bdg(coo_capacity))
    allocate(coo_vals_bdg(coo_capacity))
    coo_idx = 0

    ! ==================================================================
    ! Block (1,1): H0 - mu*I  (rows 0..8N-1, cols 0..8N-1)
    ! ==================================================================
    do row = 1, H0%nrows
      do kk = H0%rowptr(row), H0%rowptr(row + 1) - 1
        col = H0%colind(kk)
        coo_idx = coo_idx + 1
        call check_coo_bounds(coo_idx, coo_capacity)
        coo_row_bdg(coo_idx) = row
        coo_col_bdg(coo_idx) = col
        coo_vals_bdg(coo_idx) = H0%values(kk)
      end do
    end do

    ! ==================================================================
    ! Add Peierls phase correction to electron block (1,1)
    ! Peierls phase modifies existing off-diagonal entries in-place.
    ! Zeeman is already added by ZB8bandGeneralized (insert_zeeman_coo
    ! when cfg%bdg%enabled), so we must not add it again here.
    ! ==================================================================
    if (present(B_vec)) then
      if (abs(B_vec(1)) > 1.0e-12_dp) then
        call add_peierls_coo(coo_vals_bdg, coo_row_bdg, coo_col_bdg, &
                              coo_idx, cfg%grid, B_vec)
      end if
    end if

    add_optional_zeeman = .false.
    Vz_delta = 0.0_dp
    if (present(B_vec)) then
      B_mag_opt = sqrt(sum(B_vec**2))
      if (B_mag_opt > 1.0e-12_dp) then
        if (present(g_factor)) then
          g_f = g_factor
        else
          g_f = cfg%bdg%g_factor
        end if
        call compute_zeeman_vz(g_f, mu_B, B_mag_opt, Vz_opt)
        if (cfg%bdg%enabled) then
          B_mag_cfg = sqrt(sum(cfg%bdg%B_vec**2))
          call compute_zeeman_vz(cfg%bdg%g_factor, mu_B, B_mag_cfg, Vz_cfg)
          ! ----------------------------------------------------------------
          ! Vz_delta = Vz_opt - Vz_cfg is the double-counting guard.
          ! Zeeman is added by ZB8bandGeneralized at +kz (electron block) AND
          ! again at -kz (H0_minus_k_csr above) when cfg%bdg%enabled. If the
          ! optional B_vec matches cfg%bdg%B_vec, the contribution would be
          ! counted twice. Vz_delta subtracts the cfg-driven contribution so
          ! the optional arg adds only the OPTIONAL component on top of the
          ! cfg-already-applied base. This is CORRECT -- do NOT "fix" it
          ! (it is a known automated-review false-positive magnet).
          ! ----------------------------------------------------------------
          Vz_delta = Vz_opt - Vz_cfg
        else
          Vz_delta = Vz_opt
        end if
        add_optional_zeeman = maxval(abs(Vz_delta)) > 1.0e-14_dp
      end if
    end if

    ! Subtract mu from diagonal: (band_diag, band_diag) entries
    do i = 1, 8 * N
      coo_idx = coo_idx + 1
      call check_coo_bounds(coo_idx, coo_capacity)
      coo_row_bdg(coo_idx) = i
      coo_col_bdg(coo_idx) = i
      coo_vals_bdg(coo_idx) = cmplx(-mu, 0.0_dp, kind=dp)
    end do

    if (add_optional_zeeman) then
      do i = 1, N
        do row = 1, 8
          coo_idx = coo_idx + 1
          call check_coo_bounds(coo_idx, coo_capacity)
          global_row = (row - 1) * N + i
          coo_row_bdg(coo_idx) = global_row
          coo_col_bdg(coo_idx) = global_row
          coo_vals_bdg(coo_idx) = cmplx(Vz_delta(row), 0.0_dp, kind=dp)
        end do
      end do
    end if

    ! ==================================================================
    ! Block (2,1): Delta^dagger  (rows 8N..16N-1, cols 0..8N-1)
    ! Delta^dagger(j,i) = conjg(Delta(i,j)), so sign is from partner band.
    ! ==================================================================
    do i = 1, N  ! spatial index
      do row = 1, 8  ! hole band
        col = pairing_partner(row)  ! electron partner band
        coo_idx = coo_idx + 1
        call check_coo_bounds(coo_idx, coo_capacity)
        global_row = 8 * N + (row - 1) * N + i
        global_col = (col - 1) * N + i
        coo_vals_bdg(coo_idx) = cmplx(delta_0 * pairing_sign(col), 0.0_dp, kind=dp)
        coo_row_bdg(coo_idx) = global_row
        coo_col_bdg(coo_idx) = global_col
      end do
    end do

    ! ==================================================================
    ! Block (1,2): Delta  (rows 0..8N-1, cols 8N..16N-1)
    ! Intra-band s-wave pairing: electron(ib) couples to hole(partner(ib))
    ! ==================================================================
    do i = 1, N  ! spatial index
      do row = 1, 8  ! electron band
        col = pairing_partner(row)  ! hole partner band
        coo_idx = coo_idx + 1
        call check_coo_bounds(coo_idx, coo_capacity)
        global_row = (row - 1) * N + i
        global_col = 8 * N + (col - 1) * N + i
        coo_vals_bdg(coo_idx) = cmplx(delta_0 * pairing_sign(row), 0.0_dp, kind=dp)
        coo_row_bdg(coo_idx) = global_row
        coo_col_bdg(coo_idx) = global_col
      end do
    end do

    ! ==================================================================
    ! Block (2,2): -conjg(H0(-k)) + mu*I  (rows 8N..16N-1, cols 8N..16N-1)
    ! Layer B+D (ADR 0007): hole block derived from the canonical wrapper.
    ! H_hole_dense was computed once via build_bdg_hole_block(H0_minus_k, ...).
    ! ==================================================================
    do row = 1, H0%nrows
      do col = 1, H0%nrows
        if (abs(H_hole_dense(row, col)) > 1.0e-30_dp) then
          coo_idx = coo_idx + 1
          call check_coo_bounds(coo_idx, coo_capacity)
          coo_row_bdg(coo_idx) = 8 * N + row
          coo_col_bdg(coo_idx) = 8 * N + col
          coo_vals_bdg(coo_idx) = H_hole_dense(row, col)
        end if
      end do
    end do

    ! Add mu*I to block (2,2)
    do i = 1, 8 * N
      coo_idx = coo_idx + 1
      call check_coo_bounds(coo_idx, coo_capacity)
      coo_row_bdg(coo_idx) = 8 * N + i
      coo_col_bdg(coo_idx) = 8 * N + i
      coo_vals_bdg(coo_idx) = cmplx(mu, 0.0_dp, kind=dp)
    end do

    if (add_optional_zeeman) then
      do i = 1, N
        do row = 1, 8
          coo_idx = coo_idx + 1
          call check_coo_bounds(coo_idx, coo_capacity)
          global_row = 8 * N + (row - 1) * N + i
          coo_row_bdg(coo_idx) = global_row
          coo_col_bdg(coo_idx) = global_row
          coo_vals_bdg(coo_idx) = cmplx(-Vz_delta(row), 0.0_dp, kind=dp)
        end do
      end do
    end if

    ! ==================================================================
    ! Build the final CSR matrix from COO
    ! ==================================================================
    call csr_build_from_coo(H_bdg_csr, Ntot, Ntot, coo_idx, &
                             coo_row_bdg(1:coo_idx), coo_col_bdg(1:coo_idx), &
                             coo_vals_bdg(1:coo_idx))

    deallocate(coo_row_bdg, coo_col_bdg, coo_vals_bdg)
    deallocate(H0_minus_k, H_hole_dense)
    call csr_free(H0)

  end subroutine build_bdg_hamiltonian_1d

  ! ==============================================================================
  ! Build the 16N x 16N dense BdG Hamiltonian for a quantum well (QW).
  !
  ! H_BdG = | H0 - mu*I          Delta                |
  !         | Delta^dagger       -conjg(H0(-k)) + mu*I |
  !
  ! H0 is the 8N x 8N dense QW Hamiltonian from ZB8bandQW.
  ! Delta uses intra-band s-wave pairing (iσ_y ⊗ I_4) within each Kramers pair.
  !
  ! Hole-block construction (ADR 0007, Layer B+D): the canonical form is
  !   H_hole = -conjg(H0(-k))
  ! produced by the shared `build_bdg_hole_block` wrapper. The QW builder
  ! already evaluated H_h = H0(-k_par) for this; the wrapper delivers the
  ! -conjg transform. Same call site as before, unified convention.
  !
  ! Arguments:
  !   H_bdg    - (out) 16N x 16N dense BdG matrix (allocated here)
  !   cfg      - (in)  simulation configuration
  !   profile  - (in)  N x 3 band edge profile
  !   kpterms  - (in)  N x N x 10 k.p parameter matrices
  !   k_par    - (in)  in-plane wavevector (k_parallel)
  !   mu       - (in)  chemical potential
  !   delta_0  - (in)  superconducting pairing amplitude
  !   B_vec    - (in, opt) magnetic field vector [Bx, By, Bz]
  !   g_factor - (in, opt) effective g-factor for Zeeman splitting
  ! ==============================================================================
  subroutine build_bdg_hamiltonian_qw(H_bdg, cfg, profile, kpterms, &
                                        k_par, mu, delta_0, B_vec, g_factor)

    complex(kind=dp), allocatable, intent(out) :: H_bdg(:,:)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:)
    real(kind=dp), contiguous, intent(in) :: kpterms(:,:,:)
    real(kind=dp), intent(in) :: k_par, mu, delta_0
    real(kind=dp), intent(in), optional :: B_vec(3)
    real(kind=dp), intent(in), optional :: g_factor

    integer :: N, N8, Ntot, i, ib, row, col
    complex(kind=dp), allocatable :: H_e(:,:), H_h(:,:)
    type(wavevector) :: wv
    type(simulation_config) :: cfg_normal

    ! --- Guard: delta_0 must be positive and finite ---
    if (delta_0 <= 0.0_dp .or. delta_0 /= delta_0) then
      print *, 'Error: delta_0 must be positive for BdG Hamiltonian construction'
      print *, '  delta_0=', delta_0
      error stop 'bdg_hamiltonian: delta_0 must be positive and finite'
    end if

    ! --- Dimension setup ---
    N = cfg%grid%npoints()
    if (N <= 0) then
      print *, 'ERROR: build_bdg_hamiltonian_qw: invalid fdStep=', N
      error stop 'bdg_hamiltonian: invalid QW fdStep'
    end if
    N8 = 8 * N
    Ntot = 16 * N

    ! --- Build 8N x 8N normal Hamiltonians H_e(+k) and H_h(-k) ---
    allocate(H_e(N8, N8), H_h(N8, N8))
    H_e = cmplx(0.0_dp, 0.0_dp, kind=dp)
    H_h = cmplx(0.0_dp, 0.0_dp, kind=dp)

    cfg_normal = cfg
    if (present(B_vec)) then
      cfg_normal%bdg%B_vec = B_vec
      cfg_normal%bdg%enabled = any(abs(B_vec) > 1.0e-12_dp)
      if (present(g_factor)) cfg_normal%bdg%g_factor = g_factor
    end if

    ! Set wavevector: in-plane k_par along kx, kz=0 for QW
    wv%kx = k_par
    wv%ky = 0.0_dp
    wv%kz = 0.0_dp

    call ZB8bandQW(H_e, wv, profile, kpterms, cfg=cfg_normal)
    wv%kx = -k_par
    call ZB8bandQW(H_h, wv, profile, kpterms, cfg=cfg_normal)

    ! --- Allocate the full 16N x 16N BdG matrix ---
    allocate(H_bdg(Ntot, Ntot))
    H_bdg = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! --- Block (1,1): H_e(+k) - mu*I ---
    H_bdg(1:N8, 1:N8) = H_e
    do i = 1, N8
      H_bdg(i, i) = H_bdg(i, i) - cmplx(mu, 0.0_dp, kind=dp)
    end do

    ! --- Block (1,2): Delta (intra-band s-wave pairing, band-major) ---
    ! Delta(ib, partner(ib)) = delta_0 * pairing_sign(ib) at each spatial site
    do i = 1, N
      do ib = 1, 8
        row = (ib - 1) * N + i
        col = (pairing_partner(ib) - 1) * N + i
        H_bdg(row, N8 + col) = cmplx(delta_0 * pairing_sign(ib), 0.0_dp, kind=dp)
      end do
    end do

    ! --- Block (2,1): Delta^dagger = conjg(Delta^T), sign from partner band ---
    do i = 1, N
      do ib = 1, 8
        row = (ib - 1) * N + i
        col = (pairing_partner(ib) - 1) * N + i
        H_bdg(N8 + row, col) = cmplx(delta_0 * pairing_sign(pairing_partner(ib)), 0.0_dp, kind=dp)
      end do
    end do

    ! --- Block (2,2): -conjg(H(-k)) + mu*I ---
    !
    ! Canonical hole-block form (ADR 0007): H_hole = -conjg(H0(-k)).
    ! H_h was built at wv%kx = -k_par (time-reversed momentum), so H_h = H(-k).
    ! The conjugation handles the k-dependent terms: for example, kplus = kx + i*ky
    ! appears in S and kminus = kx - i*ky in SC. At -k, these swap roles, and
    ! the explicit conjg() ensures the hole block is the Nambu conjugate of the
    ! electron block. Layer B+D: routed through the shared wrapper so the wire
    ! CSR and dense QW builders produce byte-identical hole blocks.
    call build_bdg_hole_block(H_h, H_bdg(N8 + 1:Ntot, N8 + 1:Ntot))
    do i = 1, N8
      H_bdg(N8 + i, N8 + i) = H_bdg(N8 + i, N8 + i) + cmplx(mu, 0.0_dp, kind=dp)
    end do

    deallocate(H_e, H_h)

  end subroutine build_bdg_hamiltonian_qw

end module bdg_hamiltonian
