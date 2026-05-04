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
  use magnetic_field, only: add_zeeman_coo, add_peierls_coo

  implicit none

  private

  public :: build_bdg_hamiltonian_1d

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
    real(kind=dp), intent(in) :: profile_2d(:,:)
    type(csr_matrix), intent(in) :: kpterms_2d(:)
    real(kind=dp), intent(in) :: kz, mu, delta_0
    type(wire_workspace), intent(inout) :: ws
    real(kind=dp), intent(in), optional :: B_vec(3)
    real(kind=dp), intent(in), optional :: g_factor

    type(csr_matrix) :: H0
    integer :: N, Ntot, nnz_bdg, nnz_zeeman
    integer :: coo_idx, coo_capacity
    integer, allocatable :: coo_row_bdg(:), coo_col_bdg(:)
    complex(kind=dp), allocatable :: coo_vals_bdg(:)
    ! Separate REAL COO for Zeeman/Peierls (merged during assembly)
    integer, allocatable :: coo_row_zeeman(:), coo_col_zeeman(:)
    real(kind=dp), allocatable :: coo_vals_zeeman(:)

    ! Pairing sign pattern for 8-band zinc-blende basis
    ! CB1=+1, CB2=-1, VB signs alternate, then negate for particle-hole
    real(kind=dp) :: pairing_sign(8)

    ! Temporary for -H0^T + mu*I
    integer :: i, kk, row, col, global_row, global_col, nnz_z
    real(kind=dp) :: g_f

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
    ! Block (1,2): Delta => 64*N (antidiagonal, 8x8 band matrix per spatial point)
    ! Block (2,1): Delta^dagger => 64*N (conjugate transpose, 8x8 band matrix)
    ! mu*I diagonal: 16*N (8N in each of blocks (1,1) and (2,2))
    ! Zeeman/Peierls: 8*N additional diagonal entries (added to (1,1) when B != 0)
    nnz_bdg = 2 * H0%nnz + 2 * (64 * N) + 16 * N + 8 * N

    coo_capacity = nnz_bdg + nnz_bdg / 5  ! small safety margin
    allocate(coo_row_bdg(coo_capacity))
    allocate(coo_col_bdg(coo_capacity))
    allocate(coo_vals_bdg(coo_capacity))
    nnz_zeeman = 8 * N  ! maximum Zeeman diagonal entries
    allocate(coo_row_zeeman(nnz_zeeman))
    allocate(coo_col_zeeman(nnz_zeeman))
    allocate(coo_vals_zeeman(nnz_zeeman))
    coo_idx = 0

    ! --- Pairing sign pattern: antidiag(+1,-1,+1,-1,-1,+1,+1,-1) ---
    pairing_sign = [ &
      +1.0_dp, &  ! CB1 (band 7)
      -1.0_dp, &  ! CB2 (band 8)
      +1.0_dp, &  ! HH1 (band 1)
      -1.0_dp, &  ! HH2 (band 4)
      -1.0_dp, &  ! SO1 (band 5)
      +1.0_dp, &  ! SO2 (band 6)
      +1.0_dp, &  ! LH1 (band 2)
      -1.0_dp &   ! LH2 (band 3)
    ]

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
    ! Add Peierls and Zeeman corrections to electron block (1,1)
    ! Peierls phase modifies existing off-diagonal entries in-place.
    ! Zeeman splitting adds new diagonal entries (merged during CSR assembly).
    ! ==================================================================
    if (present(B_vec) .and. (abs(B_vec(1)) > 1e-12_dp .or. &
                              abs(B_vec(2)) > 1e-12_dp .or. &
                              abs(B_vec(3)) > 1e-12_dp)) then
      if (present(g_factor)) then
        g_f = g_factor
      else
        g_f = cfg%bdg%g_factor
      end if
      ! Apply Peierls phase to off-diagonal entries of the already-populated
      ! electron block (1,1) in the main complex COO array.
      call add_peierls_coo(coo_vals_bdg, coo_row_bdg, coo_col_bdg, &
                            coo_idx, cfg%grid, B_vec)
      ! Add Zeeman splitting to separate real COO array (diagonal-only)
      nnz_z = 0
      call add_zeeman_coo(coo_vals_zeeman, coo_row_zeeman, coo_col_zeeman, &
                          nnz_z, cfg%grid, B_vec, g_f)
      ! Copy Zeeman entries into the main complex COO
      do i = 1, nnz_z
        coo_idx = coo_idx + 1
        call check_coo_bounds(coo_idx, coo_capacity)
        coo_row_bdg(coo_idx) = coo_row_zeeman(i)
        coo_col_bdg(coo_idx) = coo_col_zeeman(i)
        coo_vals_bdg(coo_idx) = cmplx(coo_vals_zeeman(i), 0.0_dp, kind=dp)
      end do
    end if

    ! Subtract mu from diagonal: (band_diag, band_diag) entries
    do i = 1, 8 * N
      coo_idx = coo_idx + 1
      call check_coo_bounds(coo_idx, coo_capacity)
      coo_row_bdg(coo_idx) = i
      coo_col_bdg(coo_idx) = i
      coo_vals_bdg(coo_idx) = cmplx(-mu, 0.0_dp, kind=dp)
    end do

    ! ==================================================================
    ! Block (2,1): Delta^dagger  (rows 8N..16N-1, cols 0..8N-1)
    ! Delta^dagger(i,j) = conjg(Delta(j,i))
    ! For wire: Delta is block-antidiagonal per spatial point.
    ! Delta^dagger has the pairing signs on the hole side.
    ! ==================================================================
    do i = 1, N  ! spatial index
      do row = 1, 8  ! band index (electron)
        do col = 1, 8  ! band index (hole)
          if (abs(pairing_sign(row)) < 0.5_dp) cycle  ! skip zero entries
          coo_idx = coo_idx + 1
          call check_coo_bounds(coo_idx, coo_capacity)
          ! Hole row = 8N + (i-1)*8 + row
          ! Electron col = (j-1)*8 + col
          ! But Delta is antidiagonal in band space, so hole at spatial i
          ! couples to electron at spatial i with sign pairing_sign(row)
          global_row = 8 * N + (i - 1) * 8 + row
          global_col = (i - 1) * 8 + col
          coo_row_bdg(coo_idx) = global_row
          coo_col_bdg(coo_idx) = global_col
          ! Delta^dagger(i,j) = conjg(Delta(j,i))
          ! Delta(j,i) = delta_0 * pairing_sign(col) * delta_{ij}
          ! So Delta^dagger(i,j) = conjg(delta_0 * pairing_sign(j) * delta_{ji})
          !                     = delta_0 * pairing_sign(j) * delta_{ij}
          coo_vals_bdg(coo_idx) = cmplx(delta_0 * pairing_sign(col), 0.0_dp, kind=dp)
        end do
      end do
    end do

    ! ==================================================================
    ! Block (1,2): Delta  (rows 0..8N-1, cols 8N..16N-1)
    ! Delta(i,j) = delta_0 * pairing_sign(i) * delta_{ij} (antidiagonal in band)
    ! ==================================================================
    do i = 1, N  ! spatial index
      do row = 1, 8  ! band index (electron)
        do col = 1, 8  ! band index (hole)
          if (abs(pairing_sign(row)) < 0.5_dp) cycle  ! skip zero entries
          coo_idx = coo_idx + 1
          call check_coo_bounds(coo_idx, coo_capacity)
          global_row = (i - 1) * 8 + row
          global_col = 8 * N + (i - 1) * 8 + col
          coo_row_bdg(coo_idx) = global_row
          coo_col_bdg(coo_idx) = global_col
          coo_vals_bdg(coo_idx) = cmplx(delta_0 * pairing_sign(row), 0.0_dp, kind=dp)
        end do
      end do
    end do

    ! ==================================================================
    ! Block (2,2): -H0^T + mu*I  (rows 8N..16N-1, cols 8N..16N-1)
    ! -H0^T flips the transpose relationship
    ! ==================================================================
    do row = 1, H0%nrows
      do kk = H0%rowptr(row), H0%rowptr(row + 1) - 1
        col = H0%colind(kk)
        coo_idx = coo_idx + 1
        call check_coo_bounds(coo_idx, coo_capacity)
        ! Transpose: row in H0^T becomes col, col becomes row
        ! Offset by 8N for the hole block
        coo_row_bdg(coo_idx) = 8 * N + col
        coo_col_bdg(coo_idx) = 8 * N + row
        coo_vals_bdg(coo_idx) = -conjg(H0%values(kk))  ! -H0^T: negate and conjugate
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

    ! ==================================================================
    ! Build the final CSR matrix from COO
    ! ==================================================================
    call csr_build_from_coo(H_bdg_csr, Ntot, Ntot, coo_idx, &
                             coo_row_bdg(1:coo_idx), coo_col_bdg(1:coo_idx), &
                             coo_vals_bdg(1:coo_idx))

    deallocate(coo_row_bdg, coo_col_bdg, coo_vals_bdg)
    deallocate(coo_row_zeeman, coo_col_zeeman, coo_vals_zeeman)
    call csr_free(H0)

  end subroutine build_bdg_hamiltonian_1d

end module bdg_hamiltonian