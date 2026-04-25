module exciton_solver

  use definitions
  use parameters
  use outputFunctions, only: ensure_output_dir, get_unit
  implicit none

  private

  public :: compute_exciton_binding
  public :: sommerfeld_2d
  public :: apply_excitonic_corrections

  ! Coulomb constant: e^2/(4*pi*eps0) = hbar*c/alpha_fine = 14.40 eV*AA
  real(kind=dp), parameter :: ALPHA_FINE = 137.035999084_dp
  real(kind=dp), parameter :: COULOMB_CONST = hbar * c / ALPHA_FINE  ! eV*AA

  ! Golden-section ratio
  real(kind=dp), parameter :: GR = 0.5_dp * (sqrt(5.0_dp) - 1.0_dp)

  ! Search bounds for lambda (AA)
  real(kind=dp), parameter :: LAMBDA_LO = 1.0_dp
  real(kind=dp), parameter :: LAMBDA_HI = 500.0_dp
  real(kind=dp), parameter :: LAMBDA_TOL = 0.01_dp

  ! Radial quadrature (odd for Simpson's 1/3 rule)
  integer, parameter :: NR_RADIAL = 401

contains

  ! ------------------------------------------------------------------
  ! Variational exciton binding energy for a QW (Bastard, PRB 1982).
  !
  ! Trial function:
  !   psi(r,z_e,z_h) = phi_e(z_e)*phi_h(z_h)
  !                     * sqrt(2/(pi*lambda^2)) * exp(-r/lambda)
  !
  ! E(lambda) = hbar^2/(2*mu*lambda^2)
  !           - (e^2/(4*pi*eps_r*eps0)) * I_C(lambda)
  !
  ! Minimised over lambda via golden-section search.
  ! ------------------------------------------------------------------
  subroutine compute_exciton_binding(eigvals, eigvecs, z_grid, dz, &
    & nlayers, params, numcb, numvb, fdstep, E_binding, lambda_opt, &
    & material_id)

    real(kind=dp), intent(in)    :: eigvals(:)
    complex(kind=dp), intent(in) :: eigvecs(:,:)
    real(kind=dp), intent(in)    :: z_grid(:)
    real(kind=dp), intent(in)    :: dz
    integer, intent(in)          :: nlayers
    type(paramStruct), intent(in):: params(nlayers)
    integer, intent(in)          :: numcb, numvb, fdstep
    real(kind=dp), intent(out)   :: E_binding          ! meV
    real(kind=dp), intent(out)   :: lambda_opt         ! AA
    integer, intent(in), optional :: material_id(:)    ! (fdstep)

    real(kind=dp), allocatable :: phi_e(:), phi_h(:)
    real(kind=dp) :: eps_r, mu_eff
    real(kind=dp) :: Eb_x1, Eb_x2, lam_lo, lam_hi, lam_x1, lam_x2
    integer :: cb_state, vb_state, iter
    integer, parameter :: MAX_ITER = 200
    integer, allocatable :: mat_id(:)

    cb_state = numvb + 1  ! first CB state
    vb_state = numvb      ! highest VB state (HH1)

    if (cb_state > size(eigvals)) then
      print '(A)', 'exciton_solver: no CB state found, skipping.'
      E_binding = 0.0_dp
      lambda_opt = 0.0_dp
      return
    end if

    ! Set up material_id
    allocate(mat_id(fdstep))
    if (present(material_id)) then
      mat_id = material_id(1:fdstep)
    else
      mat_id = 1
    end if

    ! 1. Extract envelope densities
    allocate(phi_e(fdstep), phi_h(fdstep))
    call extract_envelope(eigvecs, fdstep, cb_state, dz, phi_e)
    call extract_envelope(eigvecs, fdstep, vb_state, dz, phi_h)

    ! 2. Material parameters
    eps_r = avg_dielectric(phi_e, phi_h, params, nlayers, mat_id)
    mu_eff = reduced_mass(phi_e, phi_h, params, nlayers, mat_id)

    if (mu_eff <= 0.0_dp) then
      print '(A)', 'Warning: exciton reduced mass <= 0, using 0.059*m0'
      mu_eff = 0.059_dp * m0
    end if

    ! 3. Golden-section search
    lam_lo = LAMBDA_LO
    lam_hi = LAMBDA_HI
    lam_x1 = lam_hi - GR * (lam_hi - lam_lo)
    lam_x2 = lam_lo + GR * (lam_hi - lam_lo)

    Eb_x1 = total_energy(lam_x1, phi_e, phi_h, dz, fdstep, eps_r, mu_eff)
    Eb_x2 = total_energy(lam_x2, phi_e, phi_h, dz, fdstep, eps_r, mu_eff)

    do iter = 1, MAX_ITER
      if ((lam_hi - lam_lo) < LAMBDA_TOL) exit

      if (Eb_x1 < Eb_x2) then
        lam_hi = lam_x2
        lam_x2 = lam_x1
        Eb_x2 = Eb_x1
        lam_x1 = lam_hi - GR * (lam_hi - lam_lo)
        Eb_x1 = total_energy(lam_x1, phi_e, phi_h, dz, fdstep, eps_r, mu_eff)
      else
        lam_lo = lam_x1
        lam_x1 = lam_x2
        Eb_x1 = Eb_x2
        lam_x2 = lam_lo + GR * (lam_hi - lam_lo)
        Eb_x2 = total_energy(lam_x2, phi_e, phi_h, dz, fdstep, eps_r, mu_eff)
      end if
    end do

    lambda_opt = 0.5_dp * (lam_lo + lam_hi)
    E_binding = -total_energy(lambda_opt, phi_e, phi_h, dz, fdstep, eps_r, mu_eff)
    E_binding = E_binding * 1000.0_dp  ! eV -> meV

    ! 4. Output
    print '(A)', ''
    print '(A)', '=== Exciton binding energy (variational, Bastard 1982) ==='
    print '(A,ES12.4,A)', '  lambda_opt  = ', lambda_opt, ' AA'
    print '(A,ES12.4,A)', '  E_binding   = ', E_binding, ' meV'
    print '(A,ES12.4)',   '  mu/m0       = ', mu_eff / m0
    print '(A,ES12.4)',   '  eps_r       = ', eps_r
    print '(A)', ''

    call write_exciton_output(lambda_opt, E_binding, mu_eff, eps_r)

    deallocate(phi_e, phi_h, mat_id)

  end subroutine compute_exciton_binding


  ! ------------------------------------------------------------------
  ! Extract envelope density: phi(n) = sum_b |eigvecs(idx,state)|^2
  ! Band-outer ordering: idx = (b-1)*fdstep + n
  ! Normalised so that sum(phi)*dz = 1.
  ! ------------------------------------------------------------------
  subroutine extract_envelope(eigvecs, fdstep, state_idx, dz, density)
    complex(kind=dp), intent(in)  :: eigvecs(:,:)
    integer, intent(in)           :: fdstep, state_idx
    real(kind=dp), intent(in)     :: dz
    real(kind=dp), intent(out)    :: density(:)

    integer :: n, b, idx
    real(kind=dp) :: norm

    do n = 1, fdstep
      density(n) = 0.0_dp
      do b = 1, 8
        idx = (b - 1) * fdstep + n
        density(n) = density(n) + real(eigvecs(idx, state_idx) &
          & * conjg(eigvecs(idx, state_idx)), kind=dp)
      end do
    end do

    norm = sum(density) * dz
    if (norm > 0.0_dp) density = density / norm

  end subroutine extract_envelope


  ! ------------------------------------------------------------------
  ! Density-weighted average dielectric constant.
  ! ------------------------------------------------------------------
  function avg_dielectric(phi_e, phi_h, params, nlayers, mat_id) result(eps_r)
    real(kind=dp), intent(in)    :: phi_e(:), phi_h(:)
    type(paramStruct), intent(in):: params(nlayers)
    integer, intent(in)          :: nlayers, mat_id(:)
    real(kind=dp)                :: eps_r

    integer :: n, layer
    real(kind=dp) :: w, wtotal

    eps_r = 0.0_dp
    wtotal = 0.0_dp
    do n = 1, size(phi_e)
      layer = mat_id(n)
      if (layer < 1 .or. layer > nlayers) cycle
      if (params(layer)%eps0 <= 0.0_dp) cycle
      w = phi_e(n) + phi_h(n)
      eps_r = eps_r + w * params(layer)%eps0
      wtotal = wtotal + w
    end do

    if (wtotal > 0.0_dp) then
      eps_r = eps_r / wtotal
    else
      eps_r = 0.0_dp
      do layer = 1, nlayers
        eps_r = eps_r + params(layer)%eps0
      end do
      eps_r = eps_r / real(nlayers, kind=dp)
    end if

    if (eps_r < 1.0_dp) then
      print '(A,ES12.4)', 'Warning: eps_r unexpectedly small, using 12.9 (GaAs)'
      eps_r = 12.9_dp
    end if

  end function avg_dielectric


  ! ------------------------------------------------------------------
  ! Density-weighted reduced mass.
  ! 1/mu = 1/m_e + 1/m_hh
  ! m_e: meff (density-weighted over electron envelope)
  ! m_hh: 1/(gamma1 - 2*gamma2) for [001] growth
  ! ------------------------------------------------------------------
  function reduced_mass(phi_e, phi_h, params, nlayers, mat_id) result(mu_eff)
    real(kind=dp), intent(in)    :: phi_e(:), phi_h(:)
    type(paramStruct), intent(in):: params(nlayers)
    integer, intent(in)          :: nlayers, mat_id(:)
    real(kind=dp)                :: mu_eff

    integer :: n, layer
    real(kind=dp) :: m_e_ratio, m_hh_inv, we, wh, gamma1, gamma2

    ! Electron: average meff/m0
    m_e_ratio = 0.0_dp
    we = 0.0_dp
    do n = 1, size(phi_e)
      layer = mat_id(n)
      if (layer < 1 .or. layer > nlayers) cycle
      m_e_ratio = m_e_ratio + phi_e(n) * params(layer)%meff
      we = we + phi_e(n)
    end do
    if (we > 0.0_dp) m_e_ratio = m_e_ratio / we

    ! HH: average 1/(gamma1-2*gamma2)
    m_hh_inv = 0.0_dp
    wh = 0.0_dp
    do n = 1, size(phi_h)
      layer = mat_id(n)
      if (layer < 1 .or. layer > nlayers) cycle
      gamma1 = params(layer)%gamma1
      gamma2 = params(layer)%gamma2
      if (gamma1 > 2.0_dp * abs(gamma2)) then
        m_hh_inv = m_hh_inv + phi_h(n) / (gamma1 - 2.0_dp * gamma2)
      else
        m_hh_inv = m_hh_inv + phi_h(n) * 2.0_dp  ! ~ 0.5 m0
      end if
      wh = wh + phi_h(n)
    end do
    if (wh > 0.0_dp) m_hh_inv = m_hh_inv / wh

    ! 1/mu [1/m0] = 1/m_e_ratio + 1/m_hh_inv
    ! mu [m0] = 1/(1/m_e_ratio + 1/m_hh_inv)
    ! mu [eV*s^2/AA^2] = mu/m0 * m0
    mu_eff = 1.0_dp / (1.0_dp / m_e_ratio + 1.0_dp / m_hh_inv) * m0

  end function reduced_mass


  ! ------------------------------------------------------------------
  ! E(lambda) = T + V_C
  ! ------------------------------------------------------------------
  function total_energy(lambda, phi_e, phi_h, dz, fdstep, eps_r, mu_eff) &
    & result(E)

    real(kind=dp), intent(in) :: lambda, phi_e(:), phi_h(:), dz
    integer, intent(in)       :: fdstep
    real(kind=dp), intent(in) :: eps_r, mu_eff
    real(kind=dp)             :: E

    real(kind=dp) :: T_kinetic, V_coulomb, I_C

    T_kinetic = hbar**2 / (2.0_dp * mu_eff * lambda**2)
    I_C = coulomb_integral(lambda, phi_e, phi_h, dz, fdstep)
    V_coulomb = -(COULOMB_CONST / eps_r) * I_C

    E = T_kinetic + V_coulomb

  end function total_energy


  ! ------------------------------------------------------------------
  ! Coulomb integral I_C = <1/r> for the 1s trial function.
  !
  ! I_C = (4/lambda^2) * sum_{n,m} phi_e(n)*phi_h(m)
  !       * I_r(|z_n-z_m|, lambda) * dz^2
  ! ------------------------------------------------------------------
  function coulomb_integral(lambda, phi_e, phi_h, dz, fdstep) result(I_C)
    real(kind=dp), intent(in) :: lambda, phi_e(:), phi_h(:), dz
    integer, intent(in)       :: fdstep
    real(kind=dp)             :: I_C

    integer :: n, m
    real(kind=dp) :: dz_nm, Ir, prefactor

    prefactor = 4.0_dp / lambda**2

    I_C = 0.0_dp
    do n = 1, fdstep
      do m = 1, fdstep
        dz_nm = abs(real(n - m, kind=dp) * dz)
        Ir = radial_integral(dz_nm, lambda)
        I_C = I_C + phi_e(n) * phi_h(m) * Ir
      end do
    end do

    I_C = prefactor * I_C * dz * dz

  end function coulomb_integral


  ! ------------------------------------------------------------------
  ! Radial Coulomb integral:
  ! I_r(dz, lambda) = (lambda/2) * int_0^inf u*exp(-u)/sqrt(u^2+alpha^2) du
  ! where alpha = 2*dz/lambda.  Evaluated by Simpson's rule.
  ! ------------------------------------------------------------------
  function radial_integral(dz_val, lambda) result(Ir)
    real(kind=dp), intent(in) :: dz_val, lambda
    real(kind=dp)             :: Ir

    integer :: i
    real(kind=dp) :: alpha, u, du, u_max, h, wsum
    real(kind=dp) :: f(NR_RADIAL)

    if (lambda < 1.0e-10_dp) then
      Ir = 0.0_dp
      return
    end if

    alpha = 2.0_dp * dz_val / lambda
    u_max = 30.0_dp
    du = u_max / real(NR_RADIAL - 1, kind=dp)

    do i = 1, NR_RADIAL
      u = real(i - 1, kind=dp) * du
      if (u < 1.0e-30_dp) then
        if (alpha < 1.0e-30_dp) then
          f(i) = 1.0_dp
        else
          f(i) = 0.0_dp
        end if
      else
        f(i) = u * exp(-u) / sqrt(u**2 + alpha**2)
      end if
    end do

    h = du / 3.0_dp
    wsum = f(1) + f(NR_RADIAL)
    do i = 2, NR_RADIAL - 1, 2
      wsum = wsum + 4.0_dp * f(i)
    end do
    do i = 3, NR_RADIAL - 1, 2
      wsum = wsum + 2.0_dp * f(i)
    end do

    Ir = (lambda / 2.0_dp) * h * wsum

  end function radial_integral


  ! ------------------------------------------------------------------
  ! Write exciton results to output/exciton.dat
  ! ------------------------------------------------------------------
  subroutine write_exciton_output(lambda_opt, E_binding, mu_eff, eps_r)
    real(kind=dp), intent(in) :: lambda_opt    ! AA
    real(kind=dp), intent(in) :: E_binding     ! meV
    real(kind=dp), intent(in) :: mu_eff        ! eV*s^2/AA^2
    real(kind=dp), intent(in) :: eps_r

    integer :: iounit

    call ensure_output_dir()
    call get_unit(iounit)

    open(unit=iounit, file='output/exciton.dat', status='replace', &
      & action='write')
    write(iounit, '(A)') '# Exciton binding energy (variational, Bastard 1982)'
    write(iounit, '(A)') '# lambda_opt(AA)  E_binding(meV)  mu/m0  eps_r'
    write(iounit, '(4(ES14.6,2X))') lambda_opt, E_binding, mu_eff / m0, eps_r
    close(iounit)

    print '(A)', '  Exciton results written to output/exciton.dat'

  end subroutine write_exciton_output

  ! ------------------------------------------------------------------
  ! 2D Sommerfeld enhancement factor for the absorption continuum.
  ! ------------------------------------------------------------------
  function sommerfeld_2d(E_excess, E_binding) result(S)
    real(kind=dp), intent(in) :: E_excess    ! energy above band edge (eV)
    real(kind=dp), intent(in) :: E_binding    ! exciton binding energy (eV, positive)
    real(kind=dp) :: S
    real(kind=dp) :: D, x

    if (E_binding < 1.0e-10_dp .or. E_excess < 0.0_dp) then
      S = 1.0_dp
      return
    end if

    D = E_excess / (E_binding / 4.0_dp)

    if (D < 1.0e-10_dp) then
      S = 2.0_dp
      return
    end if

    x = pi_dp / sqrt(D)

    if (x > 500.0_dp) then
      S = 2.0_dp
      return
    end if

    S = exp(x) / cosh(x)

  end function sommerfeld_2d


  ! ------------------------------------------------------------------
  ! Local pseudo-Voigt lineshape.
  ! ------------------------------------------------------------------
  function lineshape_voigt_local(E, E0, gamma_l, gamma_g) result(V)
    real(kind=dp), intent(in) :: E, E0, gamma_l, gamma_g
    real(kind=dp) :: V

    real(kind=dp) :: fwhm_l, fwhm_g, eta
    real(kind=dp) :: lorentz, gaussian
    real(kind=dp) :: x, sigma

    fwhm_l = 2.0_dp * gamma_l
    fwhm_g = 2.0_dp * gamma_g

    if (fwhm_l + fwhm_g > 0.0_dp) then
      eta = fwhm_l / (fwhm_l + fwhm_g)
    else
      eta = 0.5_dp
    end if

    if (gamma_l > 0.0_dp) then
      lorentz = gamma_l / (pi_dp * ((E - E0)**2 + gamma_l**2))
    else
      lorentz = 0.0_dp
    end if

    if (gamma_g > 0.0_dp) then
      sigma = fwhm_g / (2.0_dp * sqrt(2.0_dp * log(2.0_dp)))
      x = (E - E0) / sigma
      gaussian = exp(-0.5_dp * x**2) / (sigma * sqrt(2.0_dp * pi_dp))
    else
      gaussian = 0.0_dp
    end if

    V = eta * lorentz + (1.0_dp - eta) * gaussian

  end function lineshape_voigt_local


  ! ------------------------------------------------------------------
  ! Apply excitonic corrections to the absorption spectrum.
  ! ------------------------------------------------------------------
  subroutine apply_excitonic_corrections(E_grid, alpha_te, alpha_tm, &
    & E_gap, E_binding, optcfg)
    real(kind=dp), intent(in)    :: E_grid(:)
    real(kind=dp), intent(inout) :: alpha_te(:), alpha_tm(:)
    real(kind=dp), intent(in)    :: E_gap       ! band gap (eV)
    real(kind=dp), intent(in)    :: E_binding    ! binding energy (meV, positive)
    type(optics_config), intent(in) :: optcfg

    integer :: ie, npts
    real(kind=dp) :: E_excess, S, E_exciton
    real(kind=dp) :: gamma_l, gamma_g
    real(kind=dp) :: A_exciton, peak_norm
    real(kind=dp) :: Eb_eV

    Eb_eV = E_binding * 1.0e-3_dp

    if (Eb_eV < 1.0e-10_dp) return

    npts = size(E_grid)
    if (npts == 0) return
    if (size(alpha_te) /= npts .or. size(alpha_tm) /= npts) return

    gamma_l = optcfg%linewidth_lorentzian / 2.0_dp
    gamma_g = optcfg%linewidth_gaussian / 2.0_dp
    E_exciton = E_gap - Eb_eV

    A_exciton = 0.0_dp
    peak_norm = 0.0_dp
    do ie = 1, npts
      if (E_grid(ie) > E_gap .and. E_grid(ie) < E_gap + 5.0_dp * Eb_eV) then
        A_exciton = A_exciton + alpha_te(ie)
        peak_norm = peak_norm + 1.0_dp
      end if
    end do

    if (peak_norm > 0.0_dp .and. gamma_l > 0.0_dp) then
      A_exciton = (A_exciton / peak_norm) / gamma_l
    else if (gamma_l > 0.0_dp) then
      A_exciton = 0.0_dp
    else
      A_exciton = 0.0_dp
    end if

    do ie = 1, npts
      if (E_grid(ie) >= E_gap) then
        E_excess = E_grid(ie) - E_gap
        S = sommerfeld_2d(E_excess, Eb_eV)
        alpha_te(ie) = alpha_te(ie) * S
        alpha_tm(ie) = alpha_tm(ie) * S
      end if

      if (A_exciton > 0.0_dp) then
        alpha_te(ie) = alpha_te(ie) + A_exciton &
          & * lineshape_voigt_local(E_grid(ie), E_exciton, gamma_l, gamma_g)
        alpha_tm(ie) = alpha_tm(ie) + A_exciton &
          & * lineshape_voigt_local(E_grid(ie), E_exciton, gamma_l, gamma_g)
      end if
    end do

  end subroutine apply_excitonic_corrections

end module exciton_solver
