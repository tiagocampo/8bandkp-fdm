module scattering_solver
  ! LO-phonon Froehlich scattering rates for intersubband transitions
  ! in a quantum well.  Follows Ferreira & Bastard, PRB 40, 1074 (1989).
  !
  ! For each pair of conduction subbands (i, j) at k_par = 0:
  !   1. Extract envelope functions psi_i(z), psi_j(z)
  !   2. Compute the form factor F_ij = sum_n psi_i(z_n) psi_j(z_n) dz
  !   3. Compute Froehlich coupling constant
  !   4. Evaluate Bose-Einstein phonon occupation
  !   5. Compute emission and absorption rates
  !
  ! Temperature is taken from cfg%sc%temperature when SC is active,
  ! otherwise defaults to 300 K.

  use definitions
  use parameters
  use outputFunctions, only: ensure_output_dir, get_unit
  implicit none

  private

  public :: compute_phonon_scattering

  ! Number of q_z points for the form-factor integration
  integer, parameter :: NQZ = 256

  ! Integration cutoff for q_z (1/AA)
  real(kind=dp), parameter :: QZ_MAX = 2.0_dp

  ! Default temperature when SC is not active (K)
  real(kind=dp), parameter :: T_DEFAULT = 300.0_dp

contains

  ! ------------------------------------------------------------------
  ! Main entry point: compute LO-phonon scattering rates between all
  ! CB-CB subband pairs.
  !
  ! eigvals(k, state) -- eigenvalues at k_par = 0 (eV)
  ! eigvecs(:, state) -- eigenvectors at k_par = 0
  ! z_grid(:)         -- spatial grid (AA)
  ! params(:)         -- material parameters per layer
  ! dz                -- grid spacing (AA)
  ! numcb, numvb      -- number of conduction / valence bands
  ! fdstep            -- number of spatial grid points
  ! ------------------------------------------------------------------
  subroutine compute_phonon_scattering(cfg, eigvals, eigvecs, z_grid, &
    & params, dz, numcb, numvb, fdstep)

    type(simulation_config), intent(in) :: cfg
    real(kind=dp), intent(in)    :: eigvals(:)
    complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)
    real(kind=dp), intent(in)    :: z_grid(:)
    type(paramStruct), intent(in):: params(:)
    real(kind=dp), intent(in)    :: dz
    integer, intent(in)          :: numcb, numvb, fdstep

    integer :: i, j, pair_count
    integer :: cb_lo, cb_hi, ncb, ncb_raw
    integer, allocatable :: cb_state_idx(:)
    real(kind=dp) :: temperature, omega_lo, eps_inf, eps_0
    real(kind=dp), parameter :: e0_AA = e0 / 10.0_dp   ! C/(V*AA)
    real(kind=dp) :: coupling_eVAA, m0_over_hbar2_AA, N_lo
    real(kind=dp), allocatable :: phi_i(:), phi_j(:)
    real(kind=dp), allocatable :: rate_em(:,:), rate_ab(:,:)
    real(kind=dp) :: E_i, E_j, dE, r_em, r_ab
    integer :: iounit, ios
    real(kind=dp), parameter :: KRAMERS_TOL_EV = 1.0e-6_dp

    ! First CB state index in the eigenvalue array
    cb_lo = numvb + 1
    ncb_raw = numcb
    cb_hi = cb_lo + ncb_raw - 1

    if (cb_hi > size(eigvals)) then
      print '(A)', 'scattering_solver: not enough eigenvalues for CB states.'
      return
    end if

    if (ncb_raw < 2) then
      print '(A)', 'scattering_solver: only one CB subband, no intersubband scattering.'
      return
    end if

    ! Collapse consecutive Kramers-degenerate partners into one physical
    ! subband representative so intersubband scattering does not include
    ! zero-energy spin-partner transitions.
    allocate(cb_state_idx(ncb_raw))
    ncb = 0
    do i = cb_lo, cb_hi
      if (ncb == 0 .or. abs(eigvals(i) - eigvals(cb_state_idx(ncb))) > KRAMERS_TOL_EV) then
        ncb = ncb + 1
        cb_state_idx(ncb) = i
      end if
    end do

    if (ncb < 2) then
      print '(A)', 'scattering_solver: fewer than two unique CB subbands after'
      print '(A)', '  collapsing Kramers partners. No intersubband scattering.'
      deallocate(cb_state_idx)
      return
    end if

    ! Material parameters from scattering config
    omega_lo = cfg%scattering%phonon_energy / hbar  ! rad/s
    eps_inf  = cfg%scattering%eps_inf
    eps_0    = cfg%scattering%eps_0

    ! Temperature: use SC temperature if available, else default
    if (cfg%sc%enabled == 1 .and. cfg%sc%temperature > 0.0_dp) then
      temperature = cfg%sc%temperature
    else
      temperature = T_DEFAULT
    end if

    ! Froehlich coupling in eV*AA:
    !   e^2/(4*pi*eps0) = 14.4 eV*AA  (standard Coulomb constant)
    ! Coupling includes the dielectric factor (1/eps_inf - 1/eps_0).
    ! Rate formula (Ferreira & Bastard 1989):
    !   Gamma = coupling * omega_LO * (m0/hbar^2) * (N+1) * J_ij
    ! where J_ij = (L_z/2pi) * integral over q_z of |I_ij(q_z)|^2 * H(q_z)
    ! Units: [eV*AA] * [1/s] * [1/(eV*AA^2)] * [AA] = [1/s]
    coupling_eVAA = e**2 / (4.0_dp * pi_dp * e0_AA) / e &
      & * (1.0_dp / eps_inf - 1.0_dp / eps_0)
    ! m0/hbar^2 in code units where hbar2O2m0 acts as eV*AA^2.
    ! When hbar2O2m0 is multiplied by k^2 in 1/AA^2, the result is in eV
    ! (the 1e-20 m^2/AA^2 and 1e20 AA^2/m^2 cancel). So we use the same
    ! convention: m0/hbar^2 = 1/(2*hbar2O2m0) gives 1/(eV*AA^2) in code units.
    m0_over_hbar2_AA = 1.0_dp / (2.0_dp * hbar2O2m0)

    ! Bose-Einstein occupation N_LO = 1/(exp(hbar*omega/(kB*T)) - 1)
    N_lo = bose_einstein(cfg%scattering%phonon_energy, temperature)

    print '(A)', ''
    print '(A)', '=== LO-phonon Froehlich scattering (Ferreira & Bastard 1989) ==='
    print '(A,ES12.4,A)', '  hbar*omega_LO = ', cfg%scattering%phonon_energy * 1000.0_dp, ' meV'
    print '(A,ES12.4,A)', '  eps_inf       = ', eps_inf
    print '(A,ES12.4,A)', '  eps_0         = ', eps_0
    print '(A,ES12.4,A)', '  temperature   = ', temperature, ' K'
    print '(A,ES12.4)',   '  N_LO          = ', N_lo
    print '(A,I0)',       '  CB states     = ', ncb_raw
    print '(A,I0)',       '  unique CB subbands = ', ncb

    ! Allocate work arrays
    allocate(phi_i(fdstep), phi_j(fdstep))
    allocate(rate_em(ncb, ncb), rate_ab(ncb, ncb))
    rate_em = 0.0_dp
    rate_ab = 0.0_dp

    pair_count = 0

    do i = 1, ncb
      do j = i + 1, ncb
        ! Extract envelopes
        call extract_cb_envelope(eigvecs, fdstep, cb_state_idx(i), dz, phi_i)
        call extract_cb_envelope(eigvecs, fdstep, cb_state_idx(j), dz, phi_j)

        E_i = eigvals(cb_state_idx(i))
        E_j = eigvals(cb_state_idx(j))
        dE  = E_j - E_i  ! always positive since j > i => E_j > E_i for CB

        ! Compute scattering rates
        call compute_pair_rates(phi_i, phi_j, z_grid, dz, fdstep, &
          & dE, coupling_eVAA, m0_over_hbar2_AA, omega_lo, N_lo, &
          & cfg%scattering%phonon_energy, &
          & r_em, r_ab)

        rate_em(i, j) = r_em
        rate_ab(i, j) = r_ab
        ! Reverse direction: i <- j
        rate_em(j, i) = r_ab  ! emission from j to i uses absorption-like kinematics
        rate_ab(j, i) = r_em  ! absorption from j to i requires E_i - E_j (negative, suppressed)

        pair_count = pair_count + 1

        print '(A,I0,A,I0,A,ES12.4,A)', '  CB', i, ' -> CB', j, &
          & ': dE = ', dE * 1000.0_dp, ' meV'
        if (r_em > 0.0_dp) then
          print '(A,ES12.4,A)', '    emission:     1/tau = ', r_em, ' 1/s'
          print '(A,ES12.4,A)', '    tau_emission:         ', 1.0_dp / r_em * 1.0e12_dp, ' ps'
        else
          print '(A)', '    emission: forbidden (energy below phonon threshold)'
        end if
        if (r_ab > 0.0_dp) then
          print '(A,ES12.4,A)', '    absorption:   1/tau = ', r_ab, ' 1/s'
          print '(A,ES12.4,A)', '    tau_absorption:       ', 1.0_dp / r_ab * 1.0e12_dp, ' ps'
        else
          print '(A)', '    absorption: rate = 0'
        end if

      end do
    end do

    ! Write output file
    call write_scattering_output(ncb, cb_state_idx(1:ncb), eigvals, rate_em, rate_ab)

    print '(A)', ''
    print '(A,I0,A)', '  Computed ', pair_count, ' subband pairs.'
    print '(A)', '  Results written to output/scattering_rates.dat'
    print '(A)', ''

    deallocate(phi_i, phi_j, rate_em, rate_ab, cb_state_idx)

  end subroutine compute_phonon_scattering


  ! ------------------------------------------------------------------
  ! Extract the CB envelope density for one state.
  ! phi(n) = sum_{b=7,8} |eigvecs((n-1)*8 + b, state)|^2
  ! Normalised so that sum(phi)*dz = 1.
  ! Bands 7 and 8 are the two CB states in the 8-band basis.
  ! ------------------------------------------------------------------
  subroutine extract_cb_envelope(eigvecs, fdstep, state_idx, dz, density)
    complex(kind=dp), intent(in), contiguous  :: eigvecs(:,:)
    integer, intent(in)           :: fdstep, state_idx
    real(kind=dp), intent(in)     :: dz
    real(kind=dp), intent(out)    :: density(:)

    integer :: n, b, idx
    real(kind=dp) :: norm

    do n = 1, fdstep
      density(n) = 0.0_dp
      do b = 7, 8
        idx = (b - 1) * fdstep + n
        density(n) = density(n) + real(eigvecs(idx, state_idx) &
          & * conjg(eigvecs(idx, state_idx)), kind=dp)
      end do
    end do

    norm = sum(density) * dz
    if (norm > 0.0_dp) density = density / norm

  end subroutine extract_cb_envelope


  ! ------------------------------------------------------------------
  ! Compute emission and absorption rates for a single subband pair.
  !
  ! Following Ferreira & Bastard (PRB 40, 1074, 1989):
  !
  !   Gamma = coupling * omega_LO * (m0/hbar^2) * (N +/- 1) * J_ij
  !
  ! where J_ij = (L_z / 2*pi) * integral dq_z |I_ij(q_z)|^2 * H(q_z)
  !
  ! and H(q_z) = ln((q_z^2 + k_plus^2)/(q_z^2 + k_minus^2))
  !              / sqrt(q_z^2 + k_minus^2)
  !
  ! with k^2(delta) = 2*m0*delta/hbar^2 in 1/AA^2.
  ! All wavevectors consistently in 1/AA.
  ! ------------------------------------------------------------------
  subroutine compute_pair_rates(phi_i, phi_j, z_grid, dz, fdstep, &
    & dE, coupling_eVAA, m0_over_hbar2_AA, omega_lo, N_lo, phonon_energy, &
    & rate_em, rate_ab)

    real(kind=dp), intent(in)  :: phi_i(:), phi_j(:), z_grid(:)
    real(kind=dp), intent(in)  :: dz, dE, coupling_eVAA, m0_over_hbar2_AA
    real(kind=dp), intent(in)  :: omega_lo, N_lo, phonon_energy
    integer, intent(in)        :: fdstep
    real(kind=dp), intent(out) :: rate_em, rate_ab

    real(kind=dp) :: qz, dqz, I_ij, I_ij2, J_em, J_ab
    real(kind=dp) :: k_lo2, k_em2, k_ab2, L_z
    real(kind=dp) :: log_arg, g_factor
    real(kind=dp) :: delta_em, delta_ab
    integer :: iq

    rate_em = 0.0_dp
    rate_ab = 0.0_dp

    L_z = real(fdstep, kind=dp) * dz

    ! LO phonon wavevector squared in 1/AA^2
    k_lo2 = 2.0_dp * m0_over_hbar2_AA * phonon_energy

    ! Energy differences after phonon exchange
    delta_em = dE - phonon_energy   ! must be > 0 for emission
    delta_ab = dE + phonon_energy   ! always > 0

    ! Additional wavevector from excess kinetic energy (1/AA^2)
    k_em2 = 2.0_dp * m0_over_hbar2_AA * delta_em
    k_ab2 = 2.0_dp * m0_over_hbar2_AA * delta_ab

    ! q_z integration using Simpson's rule
    dqz = QZ_MAX / real(NQZ - 1, kind=dp)

    J_em = 0.0_dp
    J_ab = 0.0_dp

    do iq = 1, NQZ
      qz = real(iq - 1, kind=dp) * dqz

      I_ij = form_factor(phi_i, phi_j, z_grid, dz, fdstep, qz)
      I_ij2 = I_ij * I_ij

      ! Emission kinematic factor (only if delta_em > 0)
      if (delta_em > 0.0_dp) then
        log_arg = (qz**2 + k_lo2 + k_em2) / (qz**2 + k_lo2)
        g_factor = log(log_arg) / sqrt(qz**2 + k_lo2)
        J_em = J_em + simpson_weight(iq, NQZ, dqz) * I_ij2 * g_factor
      end if

      ! Absorption kinematic factor
      log_arg = (qz**2 + k_lo2 + k_ab2) / (qz**2 + k_lo2)
      g_factor = log(log_arg) / sqrt(qz**2 + k_lo2)
      J_ab = J_ab + simpson_weight(iq, NQZ, dqz) * I_ij2 * g_factor

    end do

    ! Convert sum to integral via (L_z/2*pi) normalisation for discrete q_z
    J_em = (L_z / (2.0_dp * pi_dp)) * J_em
    J_ab = (L_z / (2.0_dp * pi_dp)) * J_ab

    ! Scattering rates (1/s):
    !   Gamma = coupling * omega * (m0/hbar^2) * (N+1) * J_ij
    ! Units: [eV*AA] * [1/s] * [1/(eV*AA^2)] * [AA] = [1/s]
    if (delta_em > 0.0_dp) then
      rate_em = coupling_eVAA * omega_lo * m0_over_hbar2_AA &
        & * (N_lo + 1.0_dp) * J_em
    end if
    rate_ab = coupling_eVAA * omega_lo * m0_over_hbar2_AA &
      & * N_lo * J_ab

  end subroutine compute_pair_rates


  ! ------------------------------------------------------------------
  ! Form factor I_ij(q_z) = sum_n psi_i(z_n) * cos(q_z * z_n) * psi_j(z_n) * dz
  ! Uses the real (cosine) part of the exponential.
  ! ------------------------------------------------------------------
  function form_factor(phi_i, phi_j, z_grid, dz, fdstep, qz) result(I_ij)
    real(kind=dp), intent(in) :: phi_i(:), phi_j(:), z_grid(:)
    real(kind=dp), intent(in) :: dz, qz
    integer, intent(in)       :: fdstep
    real(kind=dp)             :: I_ij

    integer :: n

    I_ij = 0.0_dp
    do n = 1, fdstep
      I_ij = I_ij + phi_i(n) * cos(qz * z_grid(n)) * phi_j(n) * dz
    end do

  end function form_factor


  ! ------------------------------------------------------------------
  ! Simpson 1/3 rule weight for point iq of NQZ points with spacing h.
  ! Returns w*f(iq) so that sum w_i * f_i = integral.
  ! ------------------------------------------------------------------
  function simpson_weight(iq, nq, h) result(w)
    integer, intent(in)       :: iq, nq
    real(kind=dp), intent(in) :: h
    real(kind=dp)             :: w

    if (iq == 1 .or. iq == nq) then
      w = h / 3.0_dp
    else if (mod(iq, 2) == 0) then
      w = 4.0_dp * h / 3.0_dp
    else
      w = 2.0_dp * h / 3.0_dp
    end if

  end function simpson_weight


  ! ------------------------------------------------------------------
  ! Bose-Einstein distribution N = 1/(exp(E/(kB*T)) - 1).
  ! Returns 0 if T <= 0.  Uses careful handling for small arguments.
  ! ------------------------------------------------------------------
  function bose_einstein(E_eV, T_K) result(N_q)
    real(kind=dp), intent(in) :: E_eV    ! phonon energy (eV)
    real(kind=dp), intent(in) :: T_K      ! temperature (K)
    real(kind=dp)             :: N_q

    real(kind=dp) :: x

    if (T_K <= 0.0_dp) then
      N_q = 0.0_dp
      return
    end if

    x = E_eV / (kB_eV * T_K)

    if (x > 500.0_dp) then
      ! exp(x) is huge: N ~ exp(-x) ~ 0
      N_q = 0.0_dp
    else if (x < 1.0e-10_dp) then
      ! Classical limit: N ~ kB*T / E
      N_q = 1.0_dp / x
    else
      N_q = 1.0_dp / (exp(x) - 1.0_dp)
    end if

  end function bose_einstein


  ! ------------------------------------------------------------------
  ! Write scattering rates to output/scattering_rates.dat
  ! ------------------------------------------------------------------
  subroutine write_scattering_output(ncb, cb_state_idx, eigvals, rate_em, rate_ab)
    integer, intent(in)       :: ncb
    integer, intent(in)       :: cb_state_idx(:)
    real(kind=dp), intent(in) :: eigvals(:)
    real(kind=dp), intent(in), contiguous :: rate_em(:,:), rate_ab(:,:)

    integer :: iounit, ios, i, j
    real(kind=dp) :: dE, tau_em_ps, tau_ab_ps

    call ensure_output_dir()
    call get_unit(iounit)

    open(unit=iounit, file='output/scattering_rates.dat', status='replace', &
      & action='write', iostat=ios)

    write(iounit, '(A)') '# LO-phonon scattering rates (Froehlich, Ferreira & Bastard PRB 1989)'
    write(iounit, '(A)') '# i  j  E_ij(meV)  rate_emission(1/s)' // &
      & '  rate_absorption(1/s)  tau_emission(ps)  tau_absorption(ps)'

    do i = 1, ncb
      do j = i + 1, ncb
        dE = (eigvals(cb_state_idx(j)) - eigvals(cb_state_idx(i))) * 1000.0_dp  ! meV

        tau_em_ps = 0.0_dp
        tau_ab_ps = 0.0_dp
        if (rate_em(i, j) > 0.0_dp) &
          & tau_em_ps = 1.0_dp / rate_em(i, j) * 1.0e12_dp
        if (rate_ab(i, j) > 0.0_dp) &
          & tau_ab_ps = 1.0_dp / rate_ab(i, j) * 1.0e12_dp

        write(iounit, '(2(I3,1X),6(ES14.6,1X))') &
          & i, j, dE, rate_em(i, j), rate_ab(i, j), tau_em_ps, tau_ab_ps

      end do
    end do

    close(iounit)

  end subroutine write_scattering_output


end module scattering_solver
