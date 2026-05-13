module bdg_hamiltonian

  ! ==============================================================================
  ! Bogoliubov-de Gennes (BdG) Nambu-space Hamiltonian assembly.
  !
  ! The BdG Hamiltonian is a 16N x 16N matrix in Nambu space built from the
  ! 8N x 8N electron Hamiltonian:
  !
  !   H_BdG = |  H0 - mu*I    Delta   |
  !           |  Delta^dagger  -H0^T+mu*I |
  !
  ! where H0 is the 8N x 8N wire Hamiltonian (from ZB8bandGeneralized),
  ! mu is the chemical potential, and Delta is the superconducting pairing.
  !
  ! Key properties:
  !   - BdG matrix IS Hermitian (FEAST works directly, no non-Hermitian solver)
  !   - Particle-hole symmetry: eigenvalues come in ±E pairs
  !   - Pairing matrix: Xi = delta_0 * antidiag(+1,-1,+1,-1,-1,+1,+1,-1)
  !     in zinc-blende 8-band basis
  !
  ! Basis ordering (fixed throughout codebase):
  !   Bands 1-4 = valence (HH, LH, LH, SO), 5-6 = split-off, 7-8 = conduction
  !   Nambu space doubles: 1..8N = electron, 8N+1..16N = hole
  !
  ! Reference: Fu & Kane (2008), Supriyo Datta (2018) "Quantum Transport"
  ! ==============================================================================

  use definitions
  use sparse_matrices
  use hamiltonian_wire, only: ZB8bandGeneralized, wire_workspace
  use hamiltonianConstructor, only: ZB8bandQW
  use magnetic_field, only: add_peierls_coo, compute_zeeman_vz

  implicit none

  private

  public :: build_bdg_hamiltonian_1d
  public :: build_bdg_hamiltonian_qw

  ! Pairing sign antidiagonal pattern for zinc-blende 8-band basis.
  ! Xi = delta_0 * antidiag(+1,-1,+1,-1,-1,+1,+1,-1) in bands 1..8.
  real(kind=dp), parameter :: pairing_sign_xi(8) = [ &
    +1.0_dp, &  ! band 1 (HH)
    -1.0_dp, &  ! band 2 (LH)
    +1.0_dp, &  ! band 3 (LH)
    -1.0_dp, &  ! band 4 (HH)
    -1.0_dp, &  ! band 5 (SO)
    +1.0_dp, &  ! band 6 (SO)
    +1.0_dp, &  ! band 7 (CB)
    -1.0_dp ]   ! band 8 (CB)

contains

  ! ==============================================================================
  ! Check that COO index has not exceeded pre-allocated capacity.
  ! Follows the pattern in insert_zeeman_coo (hamiltonian_wire.f90).
  ! ==============================================================================
  subroutine check_coo_bounds(coo_idx, coo_capacity)
    integer, intent(in) :: coo_idx, coo_capacity

    if (coo_idx > coo_capacity) then
      print *, 'ERROR: build_bdg_hamiltonian_1d: COO capacity exceeded'
      print *, '  coo_idx=', coo_idx, ' coo_capacity=', coo_capacity
      stop 1
    end if
  end subroutine check_coo_bounds

  ! ==============================================================================
  ! Build the 16N x 16N BdG Hamiltonian from the 8N x 8N wire Hamiltonian.
  !
  ! H_BdG = | H0 - mu*I    Delta      |
  !         | Delta^dagger  -H0^T+mu*I |
  !
  ! The pairing matrix in the zinc-blende 8-band basis is:
  !   Xi = delta_0 * antidiag(+1,-1,+1,-1,-1,+1,+1,-1)
  !
  ! The 8x8 block pattern (+1,-1,+1,-1,-1,+1,+1,-1) reflects the spinor
  ! structure: CB states get (+1,-1) for spin-up/spin-down pairing, and
  ! VB states get paired with opposite signs according to their character.
  !
  ! For the wire, each band-block is replicated at each spatial point.
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

    type(csr_matrix) :: H0
    integer :: N, Ntot, nnz_bdg
    integer :: coo_idx, coo_capacity
    integer, allocatable :: coo_row_bdg(:), coo_col_bdg(:)
    complex(kind=dp), allocatable :: coo_vals_bdg(:)

    ! Temporary for -H0^T + mu*I
    integer :: i, kk, row, col, global_row, global_col, normal_nnz
    real(kind=dp) :: Vz_opt(8), Vz_cfg(8), Vz_delta(8), g_f, B_mag_opt, B_mag_cfg
    logical :: add_optional_zeeman

    ! --- Guard: delta_0 must be positive ---
    if (delta_0 <= 0.0_dp) then
      print *, 'Error: delta_0 must be positive for BdG Hamiltonian construction'
      print *, '  delta_0=', delta_0
      stop 1
    end if

    ! --- Build H0 (8N x 8N wire Hamiltonian) ---
    call ZB8bandGeneralized(H0, kz, profile_2d, kpterms_2d, cfg, ws=ws)

    ! Validate H0 after build
    if (H0%nrows == 0 .or. H0%nnz == 0) then
      print *, 'ERROR: build_bdg_hamiltonian_1d: empty H0 from ZB8bandGeneralized'
      print *, '  nrows=', H0%nrows, ' nnz=', H0%nnz
      stop 1
    end if
    if (mod(H0%nrows, 8) /= 0) then
      print *, 'ERROR: build_bdg_hamiltonian_1d: H0 nrows not multiple of 8'
      print *, '  nrows=', H0%nrows
      stop 1
    end if

    N = H0%nrows / 8
    Ntot = 16 * N

    ! --- Estimate COO capacity ---
    ! H0 has H0%nnz nonzeros
    ! Block (1,1): H0 - mu*I  => H0%nnz + 8N (diagonal)
    ! Block (2,2): -H0^T + mu*I => H0%nnz + 8N (diagonal, transpose)
    ! Block (1,2): Delta => 8*N (antidiagonal, one entry per band per spatial point)
    ! Block (2,1): Delta^dagger => 8*N (conjugate transpose)
    ! mu*I diagonal: 16*N (8N in each of blocks (1,1) and (2,2))
    ! Optional Zeeman can add one electron and one hole diagonal entry per band/site.
    nnz_bdg = 2 * H0%nnz + 2 * (8 * N) + 16 * N + 16 * N

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
    normal_nnz = coo_idx

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
    ! Antidiagonal pairing in band-major: row = (band-1)*N + site
    ! ==================================================================
    do i = 1, N  ! spatial index
      do row = 1, 8  ! hole band
        col = 9 - row  ! antidiagonal electron partner
        coo_idx = coo_idx + 1
        call check_coo_bounds(coo_idx, coo_capacity)
        global_row = 8 * N + (row - 1) * N + i
        global_col = (col - 1) * N + i
        coo_row_bdg(coo_idx) = global_row
        coo_col_bdg(coo_idx) = global_col
        coo_vals_bdg(coo_idx) = cmplx(delta_0 * pairing_sign_xi(col), 0.0_dp, kind=dp)
      end do
    end do

    ! ==================================================================
    ! Block (1,2): Delta  (rows 0..8N-1, cols 8N..16N-1)
    ! Antidiagonal pairing in band-major: row = (band-1)*N + site
    ! ==================================================================
    do i = 1, N  ! spatial index
      do row = 1, 8  ! electron band
        col = 9 - row  ! antidiagonal hole partner
        coo_idx = coo_idx + 1
        call check_coo_bounds(coo_idx, coo_capacity)
        global_row = (row - 1) * N + i
        global_col = 8 * N + (col - 1) * N + i
        coo_row_bdg(coo_idx) = global_row
        coo_col_bdg(coo_idx) = global_col
        coo_vals_bdg(coo_idx) = cmplx(delta_0 * pairing_sign_xi(row), 0.0_dp, kind=dp)
      end do
    end do

    ! ==================================================================
    ! Block (2,2): -H0^T + mu*I  (rows 8N..16N-1, cols 8N..16N-1)
    ! -H0^T flips the transpose relationship
    ! ==================================================================
    do kk = 1, normal_nnz
      row = coo_row_bdg(kk)
      col = coo_col_bdg(kk)
      coo_idx = coo_idx + 1
      call check_coo_bounds(coo_idx, coo_capacity)
      ! Hole block is -H_e^T. Do not conjugate here: complex Peierls
      ! phases must transpose into the hole sector with their phase intact.
      coo_row_bdg(coo_idx) = 8 * N + col
      coo_col_bdg(coo_idx) = 8 * N + row
      coo_vals_bdg(coo_idx) = -coo_vals_bdg(kk)
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
    call csr_free(H0)

  end subroutine build_bdg_hamiltonian_1d

  ! ==============================================================================
  ! Build the 16N x 16N dense BdG Hamiltonian for a quantum well (QW).
  !
  ! H_BdG = | H0 - mu*I    Delta       |
  !         | Delta^dagger  -H0^T+mu*I  |
  !
  ! H0 is the 8N x 8N dense QW Hamiltonian from ZB8bandQW.
  ! Delta is diagonal with antidiagonal band pairing within each spatial site.
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

    ! --- Guard: delta_0 must be positive ---
    if (delta_0 <= 0.0_dp) then
      print *, 'Error: delta_0 must be positive for BdG Hamiltonian construction'
      print *, '  delta_0=', delta_0
      stop 1
    end if

    ! --- Dimension setup ---
    N = cfg%fdStep
    if (N <= 0) then
      print *, 'ERROR: build_bdg_hamiltonian_qw: invalid fdStep=', N
      stop 1
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

    ! --- Block (1,2): Delta (antidiagonal pairing, band-major) ---
    ! Delta(b, 9-b) = delta_0 * pairing_sign(b) at each spatial site
    do i = 1, N
      do ib = 1, 8
        row = (ib - 1) * N + i
        col = (9 - ib - 1) * N + i
        H_bdg(row, N8 + col) = cmplx(delta_0 * pairing_sign_xi(ib), 0.0_dp, kind=dp)
      end do
    end do

    ! --- Block (2,1): Delta^dagger = Delta^T (real-valued, symmetric within site) ---
    ! Since Delta is real and antidiagonal: Delta^dagger(b',b) = Delta(b,b')
    do i = 1, N
      do ib = 1, 8
        row = (ib - 1) * N + i
        col = (9 - ib - 1) * N + i
        H_bdg(N8 + row, col) = cmplx(delta_0 * pairing_sign_xi(9 - ib), 0.0_dp, kind=dp)
      end do
    end do

    ! --- Block (2,2): -conjg(H(-k)) + mu*I ---
    !
    ! This is the correct BdG hole block for all k-values. H_h was built at
    ! wv%kx = -k_par (time-reversed momentum), so H_h = H(-k). The conjugation
    ! handles the k-dependent terms: for example, kplus = kx + i*ky appears in
    ! S and kminus = kx - i*ky in SC. At -k, these swap roles, and the explicit
    ! conjg() ensures the hole block is the Nambu conjugate of the electron block.
    ! This is NOT a shortcut — it is the full, k-correct construction.
    do row = 1, N8
      do col = 1, N8
        H_bdg(N8 + row, N8 + col) = -conjg(H_h(row, col))
      end do
    end do
    do i = 1, N8
      H_bdg(N8 + i, N8 + i) = H_bdg(N8 + i, N8 + i) + cmplx(mu, 0.0_dp, kind=dp)
    end do

    deallocate(H_e, H_h)

  end subroutine build_bdg_hamiltonian_qw

end module bdg_hamiltonian
