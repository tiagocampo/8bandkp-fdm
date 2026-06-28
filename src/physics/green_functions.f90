module green_functions

  ! ==============================================================================
  ! Green function and spectral function utilities.
  !
  ! LDOS(r, E) = -(1/pi) * Im[G(r, r, E + i*eta)]  (requires PARDISO)
  ! Spectral function A(k, E) = sum_n Lorentzian(E - E_n, eta)
  !
  ! LDOS uses MKL PARDISO complex solver via pardiso_c to invert
  ! the shifted matrix A = E + i*eta - H for each energy point.
  ! ==============================================================================

  use definitions, only: dp, pi_dp, simulation_config, wavevector, grid_ngrid
  use sparse_matrices
  use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t

  use linalg, only: pardiso_c
  use hamiltonianConstructor, only: ZB8bandBulk, ZB8bandQW
  use hamiltonian_wire, only: ZB8bandGeneralized
  use eigensolver, only: eigensolver_base, eigensolver_config, eigensolver_result, &
    & eigensolver_result_free, make_eigensolver, auto_compute_energy_window, &
    & asw_evals, &
    & EIGEN_MODE_ENERGY, EIGEN_MODE_FULL
  use wire_setup_mod, only: wire_setup, wire_setup_init, wire_setup_free

  implicit none
  private

  public :: compute_spectral_function_bulk
  public :: compute_spectral_function_qw
  public :: compute_spectral_function_wire
  public :: compute_landauer_transmission_1d
  public :: compute_ldos_csr

contains

  elemental real(kind=dp) function lorentzian_broadening(E, En, eta) result(val)
    real(kind=dp), intent(in) :: E, En, eta

    val = eta / (pi_dp * ((E - En)**2 + eta**2))
  end function lorentzian_broadening

  elemental function spectral_wavevector(cfg, kval) result(wv)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in) :: kval
    type(wavevector) :: wv

    wv%kx = 0.0_dp
    wv%ky = 0.0_dp
    wv%kz = 0.0_dp
    select case (trim(cfg%wave_vector%mode))
    case ('ky')
      wv%ky = kval
    case ('kz')
      wv%kz = kval
    case ('kxky')
      wv%kx = kval
      wv%ky = kval
    case ('kxkz')
      wv%kx = kval
      wv%kz = kval
    case ('kykz')
      wv%ky = kval
      wv%kz = kval
    case default
      wv%kx = kval
    end select
  end function spectral_wavevector

  ! Conductance is returned in units of the single-channel quantum e^2/h.
  function compute_landauer_transmission_1d(E, onsite, hopping, eta) result(T)
    real(kind=dp), intent(in) :: E, onsite, hopping, eta
    real(kind=dp) :: T
    complex(kind=dp) :: z, root_term, sigma_l, sigma_r, G
    real(kind=dp) :: gamma_l, gamma_r

    if (eta <= 0.0_dp .or. abs(hopping) < 1.0e-30_dp) then
      T = 0.0_dp
      return
    end if

    z = cmplx(E - onsite, eta, kind=dp)
    root_term = sqrt(z*z - cmplx(4.0_dp * hopping*hopping, 0.0_dp, kind=dp))
    sigma_l = 0.5_dp * (z - root_term)
    sigma_r = sigma_l
    gamma_l = -2.0_dp * aimag(sigma_l)
    gamma_r = -2.0_dp * aimag(sigma_r)
    G = 1.0_dp / (z - sigma_l - sigma_r)
    T = min(1.0_dp, max(0.0_dp, gamma_l * gamma_r * abs(G)**2))
  end function compute_landauer_transmission_1d

  ! ==============================================================================
  ! Compute spectral function A(k, E) for a QW system.
  !
  ! A(k, E) = sum_n (1/pi) * eta / ((E - E_n(k))^2 + eta^2)
  !
  ! For each k-point the 8N x 8N QW Hamiltonian is built and diagonalized
  ! with zheev. The resulting eigenvalues are broadened with a Lorentzian
  ! of half-width eta.
  !
  ! Args:
  !   cfg      -- simulation configuration (confinement, grid, etc.)
  !   profile  -- material profile array (N x 3)
  !   kpterms  -- k.p parameter terms (N x N x 10)
  !   k_arr    -- k-points at which to compute A (nk)
  !   E_arr    -- energy grid (nE)
  !   eta      -- Lorentzian broadening half-width (eV)
  !   A_kE     -- output spectral function (nk x nE)
  ! ==============================================================================
  subroutine compute_spectral_function_qw(cfg, profile, kpterms, k_arr, E_arr, &
                                           eta, A_kE)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: profile(:,:), kpterms(:,:,:)
    real(kind=dp), contiguous, intent(in) :: k_arr(:), E_arr(:)
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: A_kE(:,:)

    integer :: Ngrid, nk, nE, ik, iE, i
    integer :: dim
    complex(kind=dp), allocatable :: H(:,:)
    real(kind=dp), allocatable :: evals(:)
    real(kind=dp) :: lorentz
    type(wavevector) :: wv
    type(eigensolver_config) :: spec_cfg
    class(eigensolver_base), allocatable :: spec_solver
    type(eigensolver_result) :: spec_result

    Ngrid = size(profile, 1)
    dim = 8 * Ngrid
    nk = size(k_arr)
    nE = size(E_arr)

    if (eta <= 0.0_dp) then
      allocate(A_kE(0, 0))
      return
    end if

    allocate(A_kE(nk, nE))
    A_kE = 0.0_dp

    allocate(H(dim, dim))
    allocate(evals(dim))

    ! Eigensolver: DENSE + FULL (need all eigenvalues for spectral function)
    spec_cfg%method = 'DENSE'
    spec_cfg%mode = EIGEN_MODE_FULL
    spec_cfg%nev = dim
    spec_solver = make_eigensolver(spec_cfg)

    do ik = 1, nk
      wv%kx = k_arr(ik)
      wv%ky = 0.0_dp
      wv%kz = 0.0_dp

      H = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call ZB8bandQW(H, wv, profile, kpterms, cfg=cfg)

      ! Diagonalize
      call spec_solver%solve_dense(H, spec_cfg, spec_result)
      if (.not. spec_result%converged) then
        A_kE = 0.0_dp
        call eigensolver_result_free(spec_result)
        return
      end if
      evals = spec_result%eigenvalues
      call eigensolver_result_free(spec_result)

      ! Compute A(k, E) = sum_n delta_eta(E - E_n)
      do iE = 1, nE
        do i = 1, dim
          lorentz = lorentzian_broadening(E_arr(iE), evals(i), eta)
          A_kE(ik, iE) = A_kE(ik, iE) + lorentz
        end do
      end do
    end do

    deallocate(H, evals)
    if (allocated(spec_solver)) deallocate(spec_solver)
  end subroutine compute_spectral_function_qw

  subroutine compute_spectral_function_bulk(cfg, k_arr, E_arr, eta, A_kE)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: k_arr(:), E_arr(:)
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: A_kE(:,:)

    integer :: nk, nE, ik, iE, i
    complex(kind=dp) :: H(8, 8)
    real(kind=dp) :: evals(8)
    type(wavevector) :: wv
    type(eigensolver_config) :: bulk_spec_cfg
    class(eigensolver_base), allocatable :: bulk_spec_solver
    type(eigensolver_result) :: bulk_spec_result

    nk = size(k_arr)
    nE = size(E_arr)

    if (eta <= 0.0_dp) then
      allocate(A_kE(0, 0))
      return
    end if

    allocate(A_kE(nk, nE))
    A_kE = 0.0_dp

    ! Eigensolver: DENSE + FULL
    bulk_spec_cfg%method = 'DENSE'
    bulk_spec_cfg%mode = EIGEN_MODE_FULL
    bulk_spec_cfg%nev = 8
    bulk_spec_solver = make_eigensolver(bulk_spec_cfg)

    do ik = 1, nk
      wv = spectral_wavevector(cfg, k_arr(ik))
      H = cmplx(0.0_dp, 0.0_dp, kind=dp)
      call ZB8bandBulk(H, wv, cfg%params(1:1), cfg=cfg)
      call bulk_spec_solver%solve_dense(H, bulk_spec_cfg, bulk_spec_result)
      if (.not. bulk_spec_result%converged) then
        A_kE = 0.0_dp
        call eigensolver_result_free(bulk_spec_result)
        return
      end if
      evals = bulk_spec_result%eigenvalues
      call eigensolver_result_free(bulk_spec_result)

      do iE = 1, nE
        do i = 1, 8
          A_kE(ik, iE) = A_kE(ik, iE) + lorentzian_broadening(E_arr(iE), evals(i), eta)
        end do
      end do
    end do

    if (allocated(bulk_spec_solver)) deallocate(bulk_spec_solver)
  end subroutine compute_spectral_function_bulk

  subroutine compute_spectral_function_wire(cfg, k_arr, E_arr, eta, A_kE)
    type(simulation_config), intent(in) :: cfg
    real(kind=dp), contiguous, intent(in) :: k_arr(:), E_arr(:)
    real(kind=dp), intent(in) :: eta
    real(kind=dp), allocatable, intent(out) :: A_kE(:,:)

    ! Local mutable copy: wire_setup_init populates cfg_local%strain_blocks
    ! when strain is enabled (the step the pre-#04 path omitted).
    type(simulation_config) :: cfg_local
    type(wire_setup) :: wsetup
    type(csr_matrix) :: H_csr
    type(eigensolver_config) :: eigen_cfg
    class(eigensolver_base), allocatable :: solver
    type(eigensolver_result) :: eigen_res
    integer :: nk, nE, ik, iE, i, ngrid, matrix_dim
    real(kind=dp) :: emin_auto, emax_auto

    nk = size(k_arr)
    nE = size(E_arr)

    if (eta <= 0.0_dp) then
      allocate(A_kE(0, 0))
      return
    end if

    allocate(A_kE(nk, nE))
    A_kE = 0.0_dp

    ! Strain-aware wire init (Issue #04): run confinementInitialization_2d
    ! + compute_strain + compute_bir_pikus_blocks so the wire spectral
    ! function honors strain. Pre-#04 this site called
    ! confinementInitialization_2d directly and silently dropped strain.
    cfg_local = cfg
    call wire_setup_init(wsetup, cfg_local)

    ngrid = grid_ngrid(cfg_local%grid)
    matrix_dim = 8 * ngrid
    eigen_cfg%method = 'FEAST'
    eigen_cfg%mode = EIGEN_MODE_ENERGY
    eigen_cfg%nev = max(1, min(matrix_dim, cfg_local%bands%num_cb + cfg_local%bands%num_vb))
    eigen_cfg%max_iter = 200
    eigen_cfg%tol = 1.0e-10_dp
    eigen_cfg%m0 = min(matrix_dim, max(cfg_local%solver%m0, &
      & min(matrix_dim, max(4 * eigen_cfg%nev, 128))))
    solver = make_eigensolver(eigen_cfg)

    do ik = 1, nk
      call ZB8bandGeneralized(H_csr, k_arr(ik), wsetup%profile_2d, wsetup%kpterms_2d, &
        & cfg_local, ws=wsetup%ws)
      ! KTD6 closure (ADR 0005): route the wire spectral window through
      ! apply_solver_window. The override path (cfg_local%solver%emin/emax
      ! nonzero) is honored verbatim. Otherwise, the E-grid envelope
      ! is fed through asw_evals so the Gershgorin margin convention
      ! applies uniformly. The legacy auto_compute_energy_window
      ! fallback is retained only as the narrow-case degenerate-window
      ! backstop (eigen_cfg%emin >= emax).
      call asw_evals(E_arr, cfg_local%solver%emin, cfg_local%solver%emax, &
        & eigen_cfg%emin, eigen_cfg%emax)
      if (eigen_cfg%emin >= eigen_cfg%emax) then
        call auto_compute_energy_window(H_csr, emin_auto, emax_auto)
        eigen_cfg%emin = emin_auto
        eigen_cfg%emax = emax_auto
      end if

      call solver%solve_sparse(H_csr, eigen_cfg, eigen_res)
      if (.not. eigen_res%converged) then
        print *, 'ERROR: compute_spectral_function_wire: ', solver%backend_name(), &
          ' failed at k index ', ik
        call fail_wire_spectral()
        return
      end if
      if (eigen_res%nev_found <= 0) then
        print *, 'ERROR: compute_spectral_function_wire: ', solver%backend_name(), &
          ' found no eigenvalues at k index ', ik
        call fail_wire_spectral()
        return
      end if
      if (eigen_res%m0_used > 0 .and. eigen_res%nev_found >= eigen_res%m0_used .and. &
          eigen_res%m0_used < matrix_dim) then
        print *, 'ERROR: compute_spectral_function_wire: eigensolver returned ', &
          & eigen_res%nev_found, ' states, reaching m0=', eigen_res%m0_used
        print *, '  Increase m0 or narrow the spectral energy window.'
        call fail_wire_spectral()
        return
      end if
      do iE = 1, nE
        do i = 1, eigen_res%nev_found
          A_kE(ik, iE) = A_kE(ik, iE) + &
            & lorentzian_broadening(E_arr(iE), eigen_res%eigenvalues(i), eta)
        end do
      end do
      call eigensolver_result_free(eigen_res)
      call csr_free(H_csr)
    end do

    if (allocated(solver)) deallocate(solver)
    call wire_setup_free(wsetup)

  contains

    subroutine cleanup_wire_spectral()
      call eigensolver_result_free(eigen_res)
      call csr_free(H_csr)
      if (allocated(solver)) deallocate(solver)
      call wire_setup_free(wsetup)
    end subroutine cleanup_wire_spectral

    subroutine fail_wire_spectral()
      call cleanup_wire_spectral()
      if (allocated(A_kE)) deallocate(A_kE)
      allocate(A_kE(0, 0))
    end subroutine fail_wire_spectral

  end subroutine compute_spectral_function_wire

  ! ==============================================================================
  ! Compute LDOS at each grid point for a given energy E.
  !
  ! LDOS(r, E) = -(1/pi) * Im[G(r, r, E + i*eta)]
  ! where G = (E + i*eta - H)^-1
  !
  ! Args:
  !   H          -- CSR Hamiltonian matrix (complex)
  !   E          -- energy (real, eV)
  !   eta        -- broadening (real, eV)
  !   ldos       -- output LDOS at each grid point (real, 1/eV)
  !
  ! The shifted system is: A * x = b  with  A = E + i*eta - H
  ! For diagonal G(r,r) we solve A * e_r = e_r where e_r is unit vector at row r,
  ! then G(r,r) = x(r).  We extract diagonal by solving with unit vectors
  ! at each diagonal position.
  ! ==============================================================================
  subroutine compute_ldos_csr(H, E, eta, ldos)
    type(csr_matrix), intent(in) :: H
    real(kind=dp), intent(in) :: E
    real(kind=dp), intent(in) :: eta
    real(kind=dp), intent(out) :: ldos(:)

    integer :: N, nnz, i, r
    integer :: maxfct, mnum, mtype, phase, nrhs, msglvl, error_loc
    integer(kind=c_intptr_t) :: pt(64)
    integer :: iparm(64)
    complex(kind=dp), allocatable :: a_val(:), rhs(:), sol(:)
    integer, allocatable :: ia(:), ja(:), perm(:)
    complex(kind=dp) :: shift
    logical :: found_diag

    N = H%nrows

    if (N /= size(ldos)) then
      print *, 'ERROR: compute_ldos_csr: ldos size mismatch (N=', N, ', ldos=', size(ldos), ')'
      stop 1
    end if

    shift = cmplx(E, eta, kind=dp)

    ! Copy CSR structure
    nnz = H%nnz
    allocate(a_val(nnz), ia(N+1), ja(nnz))
    a_val = -H%values
    ia = H%rowptr
    ja = H%colind
    do r = 1, N
      found_diag = .false.
      do i = ia(r), ia(r + 1) - 1
        if (ja(i) == r) then
          a_val(i) = a_val(i) + shift
          found_diag = .true.
        end if
      end do
      if (.not. found_diag) then
        print *, 'ERROR: compute_ldos_csr: missing structural diagonal at row ', r
        stop 1
      end if
    end do

    ! PARDISO setup
    pt = 0
    iparm = 0
    maxfct = 1
    mnum = 1
    mtype = 3  ! complex structurally symmetric
    nrhs = 1
    msglvl = 0
    error_loc = 0

    iparm(1) = 1   ! no solver default
    iparm(2) = 2   ! OpenMP nested dissection
    iparm(8) = 2   ! max iterative refinement steps
    iparm(10) = 13 ! perturb pivots with 1E-13
    iparm(11) = 1  ! scaling
    iparm(13) = 1  ! matching
    iparm(18) = -1 ! report number of nonzeros
    iparm(40) = 1  ! distributed matrix input (CSR)

    allocate(perm(N), rhs(N), sol(N))

    ! Phase 11: Symbolic factorization
    phase = 11
    call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                   a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                   rhs, sol, error_loc)
    if (error_loc /= 0) then
      print *, 'ERROR: compute_ldos_csr: symbolic factorization failed (error=', error_loc, ')'
      phase = -1
      call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                     a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                     rhs, sol, error_loc)
      stop 1
    end if

    ! Phase 22: Numeric factorization
    phase = 22
    call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                   a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                   rhs, sol, error_loc)
    if (error_loc /= 0) then
      print *, 'ERROR: compute_ldos_csr: numeric factorization failed (error=', error_loc, ')'
      phase = -1
      call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                     a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                     rhs, sol, error_loc)
      stop 1
    end if

    ! Extract diagonal of G by solving A * x = e_r for each r
    do r = 1, N
      rhs = cmplx(0.0_dp, 0.0_dp, kind=dp)
      rhs(r) = cmplx(1.0_dp, 0.0_dp, kind=dp)

      ! Phase 33: Solve
      phase = 33
      call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                     a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                     rhs, sol, error_loc)
      if (error_loc /= 0) then
        print *, 'ERROR: compute_ldos_csr: solve failed at r=', r, ' (error=', error_loc, ')'
        phase = -1
        call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                       a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                       rhs, sol, error_loc)
        stop 1
      end if

      ! LDOS = -(1/pi) * Im[G(r,r)] = -(1/pi) * Im[sol(r)]
      ldos(r) = -aimag(sol(r)) / pi_dp
    end do

    ! Phase -1: Release memory
    phase = -1
    call pardiso_c(pt, maxfct, mnum, mtype, phase, N, &
                   a_val, ia, ja, perm, nrhs, iparm, msglvl, &
                   rhs, sol, error_loc)

    deallocate(a_val, ia, ja, perm, rhs, sol)
  end subroutine compute_ldos_csr

end module green_functions
