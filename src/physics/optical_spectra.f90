module optical_spectra

  use definitions, only: IU, ZERO, c, dp, e, e0, hbar, hbar2O2m0, &
    kB_eV, m0, mu_B, optics_config, pi_dp
  use sparse_matrices, only: csr_matrix, csr_spmv
  use charge_density, only: fermi_dirac
  use spin_projection, only: spin_weights
  use utils, only: ensure_output_dir, get_unit
  use linalg, only: zdotc
  implicit none

  private

  public :: optics_engine, optics_init, optics_free
  public :: optics_accumulate, optics_accumulate_spontaneous
  public :: optics_apply_prefactor, optics_write_output
  public :: compute_gain_qw, gain_reset
  public :: compute_isbt_absorption, compute_intersubband_transitions

  ! Minimum transition energy threshold (eV)
  real(kind=dp), parameter :: DE_MIN = 1.0e-6_dp

  ! ------------------------------------------------------------------
  ! optics_engine: encapsulates all optical accumulation state.
  ! Public allocatable arrays for hot-path access (same pattern as
  ! csr_matrix).
  ! ------------------------------------------------------------------
  type optics_engine
    ! Stored configuration
    type(optics_config)  :: config
    ! Grid state
    integer              :: nE = 0
    ! Unconditional accumulation arrays
    real(kind=dp), allocatable :: E_grid(:)
    real(kind=dp), allocatable :: alpha_te(:)
    real(kind=dp), allocatable :: alpha_tm(:)
    real(kind=dp), allocatable :: alpha_isbt(:)
    ! Conditional: gain
    real(kind=dp), allocatable :: alpha_gain_te(:)
    real(kind=dp), allocatable :: alpha_gain_tm(:)
    ! Conditional: spontaneous emission
    real(kind=dp), allocatable :: spont_te(:)
    real(kind=dp), allocatable :: spont_tm(:)
    ! Conditional: spin-resolved
    real(kind=dp), allocatable :: alpha_te_up(:)
    real(kind=dp), allocatable :: alpha_te_dw(:)
    real(kind=dp), allocatable :: alpha_tm_up(:)
    real(kind=dp), allocatable :: alpha_tm_dw(:)
    ! Gain quasi-Fermi level state
    real(kind=dp) :: mu_e = 0.0_dp
    real(kind=dp) :: mu_h = 0.0_dp
    logical       :: gain_fermi_computed = .false.
    logical       :: was_freed = .false.
  contains
    procedure :: free => optics_free_tbp
    final :: optics_engine_finalize
  end type optics_engine

contains

  ! ------------------------------------------------------------------
  ! Private helper: compute z-dipole matrix element between two
  ! eigenstates of the 8-band k.p Hamiltonian.
  !
  ! z_ij = dz * sum_n [ sum_b conjg(eigvecs((n-1)*8+b, state_i))
  !               * z_grid(n) * eigvecs((n-1)*8+b, state_j)) ]
  !
  ! Returns the complex z_ij in units of AA.
  ! ------------------------------------------------------------------
  function z_dipole(eigvecs, z_grid, dz, fdstep, state_i, state_j) result(z_ij)
    complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)
    real(kind=dp), intent(in)    :: z_grid(:)
    real(kind=dp), intent(in)    :: dz
    integer, intent(in)          :: fdstep, state_i, state_j
    complex(kind=dp)             :: z_ij

    integer :: n, b, idx

    z_ij = ZERO
    do n = 1, fdstep
      do b = 1, 8
        idx = (b - 1) * fdstep + n
        z_ij = z_ij + conjg(eigvecs(idx, state_i)) &
          & * z_grid(n) * eigvecs(idx, state_j)
      end do
    end do
    z_ij = z_ij * dz

  end function z_dipole


  ! ------------------------------------------------------------------
  ! Pseudo-Voigt profile: weighted sum of Lorentzian and Gaussian.
  !
  ! V(E) = eta * L(E; E0, gamma_l) + (1-eta) * G(E; E0, gamma_g)
  !
  ! where eta = fwhm_L / (fwhm_L + fwhm_G) (Thompson mixing parameter).
  !
  ! Input gamma_l, gamma_g are HWHM (half-width at half-maximum).
  ! Returns the profile value (un-normalized in the sense that the
  ! peak value at E=E0 is ~1/(pi*gamma_l) for pure Lorentzian).
  ! ------------------------------------------------------------------
  function lineshape_voigt(E, E0, gamma_l, gamma_g) result(V)
    real(kind=dp), intent(in) :: E, E0, gamma_l, gamma_g
    real(kind=dp) :: V

    real(kind=dp) :: fwhm_l, fwhm_g, eta
    real(kind=dp) :: lorentz, gaussian
    real(kind=dp) :: x, sigma

    fwhm_l = 2.0_dp * gamma_l
    fwhm_g = 2.0_dp * gamma_g

    ! Mixing parameter (Thompson et al. approximation)
    if (fwhm_l + fwhm_g > 0.0_dp) then
      eta = fwhm_l / (fwhm_l + fwhm_g)
    else
      eta = 0.5_dp
    end if

    ! Lorentzian: L(E) = gamma_l / (pi * ((E-E0)^2 + gamma_l^2))
    if (gamma_l > 0.0_dp) then
      lorentz = gamma_l / (pi_dp * ((E - E0)**2 + gamma_l**2))
    else
      lorentz = 0.0_dp
    end if

    ! Gaussian: G(E) = (1/(sigma*sqrt(2*pi))) * exp(-(E-E0)^2/(2*sigma^2))
    ! sigma = fwhm_g / (2*sqrt(2*ln2))
    if (gamma_g > 0.0_dp) then
      sigma = fwhm_g / (2.0_dp * sqrt(2.0_dp * log(2.0_dp)))
      x = (E - E0) / sigma
      gaussian = exp(-0.5_dp * x**2) / (sigma * sqrt(2.0_dp * pi_dp))
    else
      gaussian = 0.0_dp
    end if

    V = eta * lorentz + (1.0_dp - eta) * gaussian

  end function lineshape_voigt

  ! ------------------------------------------------------------------
  ! Intersubband transitions: z-dipole matrix elements between CB
  ! subbands.  Outputs transition table to output/isbt_transitions.dat.
  !
  ! z_ij = sum_n [ sum_b conjg(eigvecs((n-1)*8+b, i)) * z_n
  !               * eigvecs((n-1)*8+b, j)) ] * dz
  !
  ! Oscillator strength: f_ij = E_ij * |z_ij|^2 / hbar2O2m0
  ! ------------------------------------------------------------------
  subroutine compute_intersubband_transitions(eigvals, eigvecs, z_grid, dz, &
    & numcb, numvb, fdstep, transitions_file)

    real(kind=dp), intent(in), contiguous    :: eigvals(:)           ! (nstates) eigenvalues ascending
    complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)         ! (8*fdstep, nstates)
    real(kind=dp), intent(in), contiguous    :: z_grid(:)            ! (fdstep) z-coords (AA)
    real(kind=dp), intent(in)    :: dz                   ! grid spacing (AA)
    integer, intent(in)          :: numcb, numvb, fdstep
    character(len=*), intent(in) :: transitions_file

    integer :: i, j, state_i, state_j
    integer :: iounit
    real(kind=dp) :: E_ij, z_ij_re, z_ij_im, z_ij_abs2, f_ij
    complex(kind=dp) :: z_ij

    call ensure_output_dir()
    call get_unit(iounit)
    open(unit=iounit, file=transitions_file, status='replace', action='write')
    write(iounit, '(a)') '# Intersubband transitions (ISBT)'
    write(iounit, '(a)') '# i  j  E_ij(eV)  Re(z_ij)(AA)  Im(z_ij)(AA)' &
      & // '  |z_ij|^2(AA^2)  f_ij(osc.str.)'
    write(iounit, '(a)') '# CB state indices are 1-based within CB manifold'

    ! Loop over all CB-CB pairs (i < j)
    do i = 1, numcb
      do j = i + 1, numcb
        ! Map to eigenvalue index: CB starts at numvb+1
        state_i = numvb + i
        state_j = numvb + j

        ! Transition energy
        E_ij = eigvals(state_j) - eigvals(state_i)
        if (E_ij < DE_MIN) cycle

        ! z-dipole matrix element via helper
        z_ij = z_dipole(eigvecs, z_grid, dz, fdstep, state_i, state_j)

        z_ij_re = real(z_ij, kind=dp)
        z_ij_im = aimag(z_ij)
        z_ij_abs2 = z_ij_re**2 + z_ij_im**2

        ! Oscillator strength: f_ij = E_ij * |z_ij|^2 / hbar2O2m0
        f_ij = E_ij * z_ij_abs2 / hbar2O2m0

        write(iounit, '(i4,1x,i4,1x,es14.6,1x,es14.6,1x,es14.6,1x,es14.6,1x,es14.6)') &
          & i, j, E_ij, z_ij_re, z_ij_im, z_ij_abs2, f_ij
      end do
    end do

    close(iounit)

  end subroutine compute_intersubband_transitions



  ! ------------------------------------------------------------------
  ! Find quasi-Fermi level by bisection for a given 2D carrier density.
  !
  ! For parabolic subbands in 2D:
  !   n_2D = sum_i (m_eff_i * kB*T) / (pi*hbar^2) * ln(1 + exp((mu-E_i)/(kBT)))
  !
  ! Since we work with the raw k.p eigenvalues (not known m_eff), we
  ! use a simplified model: approximate the 2D DOS per subband using
  ! the free-electron mass m0.  This gives a first estimate of mu;
  ! the gain spectrum itself is computed via k-space integration so
  ! the exact quasi-Fermi level only needs to be self-consistent
  ! within the k_sweep.
  !
  ! For the gain calculation, we use a simpler approach:
  !   n_2D = sum_i N_2D_i * ln(1 + exp((mu - E_i) / (kB*T)))
  ! where N_2D_i = m0 * kB * T / (pi * hbar^2) is the 2D DOS per subband.
  !
  ! NOTE: carrier_density is in cm^-2.  We convert to AA^-2 internally
  ! (1 cm^-2 = 1e-16 AA^-2).
  ! ------------------------------------------------------------------
  function find_quasi_fermi(eigvals, temperature, carrier_density_cm2, &
    & n_states, state_offset) result(mu)

    real(kind=dp), intent(in), contiguous :: eigvals(:)       ! all eigenvalues at k=0
    real(kind=dp), intent(in) :: temperature       ! K
    real(kind=dp), intent(in) :: carrier_density_cm2  ! cm^-2
    integer, intent(in)       :: n_states          ! number of subbands to sum
    integer, intent(in)       :: state_offset      ! index offset (numvb for CB)
    real(kind=dp)             :: mu                ! quasi-Fermi level (eV)

    real(kind=dp), parameter :: CM2_TO_AA2 = 1.0e-16_dp  ! cm^-2 -> AA^-2
    integer, parameter :: MAX_ITER = 200
    real(kind=dp), parameter :: BISECT_TOL = 1.0e-10_dp

    real(kind=dp) :: n_target, n_current, kBT
    real(kind=dp) :: N_2D_per_subband, x
    real(kind=dp) :: mu_lo, mu_hi, dmu
    integer :: iter, s

    kBT = kB_eV * temperature

    ! 2D DOS per subband (AA^-2): N_2D = m0 * kB*T / (pi * hbar^2)
    ! m0 [=] eV*s^2/AA^2, kB*T [=] eV, hbar^2 [=] (eV*s)^2
    ! N_2D [=] eV*s^2/AA^2 * eV / ((eV*s)^2) = 1/AA^2
    N_2D_per_subband = m0 * kBT / (pi_dp * hbar**2)

    ! Target carrier density in AA^-2
    n_target = carrier_density_cm2 * CM2_TO_AA2

    ! Bisection bounds: start wide around the subband energies
    mu_lo = minval(eigvals(state_offset+1:state_offset+n_states)) - 50.0_dp * kBT
    mu_hi = maxval(eigvals(state_offset+1:state_offset+n_states)) + 50.0_dp * kBT

    ! Bisection loop
    do iter = 1, MAX_ITER
      mu = 0.5_dp * (mu_lo + mu_hi)
      dmu = mu_hi - mu_lo

      ! Compute n_2D(mu) = sum over subbands
      n_current = 0.0_dp
      do s = 1, n_states
        x = (mu - eigvals(state_offset + s)) / kBT
        ! Clamp to avoid overflow in exp
        if (x > 500.0_dp) then
          n_current = n_current + N_2D_per_subband * x
        else
          n_current = n_current + N_2D_per_subband * log(1.0_dp + exp(x))
        end if
      end do

      if (n_current < n_target) then
        mu_lo = mu
      else
        mu_hi = mu
      end if

      if (dmu < BISECT_TOL * kBT) exit
    end do

  end function find_quasi_fermi


  ! ------------------------------------------------------------------
  ! Find quasi-Fermi level for holes in the VB manifold.
  !
  ! Hole density: p = sum_i (1 - f(E_i, mu_h)) = sum_i f(-E_i, -mu_h)
  ! We want p = carrier_density.
  !
  ! Using the identity 1 - f(E, mu) = f(-E, -mu), we negate the VB
  ! eigenvalues and search for -mu_h using the same 2D DOS formula:
  !   p = sum_i N_2D * ln(1 + exp((mu_h' - (-E_i)) / kBT))
  ! where mu_h' = -mu_h is the returned value, and -E_i are the
  ! negated VB energies.  The actual mu_h = -mu_h'.
  ! ------------------------------------------------------------------
  function find_quasi_fermi_holes(eigvals, temperature, carrier_density_cm2, &
    & numvb) result(mu_h)

    real(kind=dp), intent(in) :: eigvals(:)       ! all eigenvalues at k=0
    real(kind=dp), intent(in) :: temperature       ! K
    real(kind=dp), intent(in) :: carrier_density_cm2  ! cm^-2
    integer, intent(in)       :: numvb             ! number of VB states
    real(kind=dp)             :: mu_h              ! hole quasi-Fermi level (eV)

    real(kind=dp), parameter :: CM2_TO_AA2 = 1.0e-16_dp
    integer, parameter :: MAX_ITER = 200
    real(kind=dp), parameter :: BISECT_TOL = 1.0e-10_dp

    real(kind=dp) :: n_target, n_current, kBT
    real(kind=dp) :: N_2D_per_subband, x
    real(kind=dp) :: mu_lo, mu_hi, dmu, mu_prime
    real(kind=dp), allocatable :: neg_Evb(:)
    integer :: iter, s

    kBT = kB_eV * temperature
    N_2D_per_subband = m0 * kBT / (pi_dp * hbar**2)
    n_target = carrier_density_cm2 * CM2_TO_AA2

    ! Negate VB eigenvalues for hole calculation
    allocate(neg_Evb(numvb))
    do s = 1, numvb
      neg_Evb(s) = -eigvals(s)
    end do

    ! Bisection for mu' = -mu_h
    mu_lo = minval(neg_Evb(1:numvb)) - 50.0_dp * kBT
    mu_hi = maxval(neg_Evb(1:numvb)) + 50.0_dp * kBT

    do iter = 1, MAX_ITER
      mu_prime = 0.5_dp * (mu_lo + mu_hi)
      dmu = mu_hi - mu_lo

      n_current = 0.0_dp
      do s = 1, numvb
        x = (mu_prime - neg_Evb(s)) / kBT
        if (x > 500.0_dp) then
          n_current = n_current + N_2D_per_subband * x
        else
          n_current = n_current + N_2D_per_subband * log(1.0_dp + exp(x))
        end if
      end do

      if (n_current < n_target) then
        mu_lo = mu_prime
      else
        mu_hi = mu_prime
      end if

      if (dmu < BISECT_TOL * kBT) exit
    end do

    ! mu_h = -mu_prime
    mu_h = -mu_prime

    deallocate(neg_Evb)
  end function find_quasi_fermi_holes


  ! ==================================================================
  ! optics_engine lifecycle routines
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Initialize the optics_engine: store config, allocate and zero all
  ! arrays based on config flags.
  ! ------------------------------------------------------------------
  subroutine optics_init(oe, optcfg)
    type(optics_engine), intent(inout) :: oe
    type(optics_config), intent(in)    :: optcfg

    integer :: i
    real(kind=dp) :: dE

    ! Defensive cleanup in case of prior initialization
    call optics_free(oe)
    oe%was_freed = .false.

    ! Store config
    oe%config = optcfg

    ! Grid size
    oe%nE = optcfg%num_energy_points

    ! Unconditional arrays: E_grid, alpha_te, alpha_tm, alpha_isbt
    allocate(oe%E_grid(oe%nE))
    allocate(oe%alpha_te(oe%nE))
    allocate(oe%alpha_tm(oe%nE))
    allocate(oe%alpha_isbt(oe%nE))

    ! Build linear energy grid
    dE = (optcfg%E_max - optcfg%E_min) / max(oe%nE - 1, 1)
    do i = 1, oe%nE
      oe%E_grid(i) = optcfg%E_min + (i - 1) * dE
    end do

    ! Zero unconditional arrays
    oe%alpha_te = 0.0_dp
    oe%alpha_tm = 0.0_dp
    oe%alpha_isbt = 0.0_dp

    ! Conditional: gain
    if (optcfg%gain_enabled) then
      allocate(oe%alpha_gain_te(oe%nE))
      allocate(oe%alpha_gain_tm(oe%nE))
      oe%alpha_gain_te = 0.0_dp
      oe%alpha_gain_tm = 0.0_dp
    end if

    ! Conditional: spontaneous emission
    if (optcfg%spontaneous_enabled) then
      allocate(oe%spont_te(oe%nE))
      allocate(oe%spont_tm(oe%nE))
      oe%spont_te = 0.0_dp
      oe%spont_tm = 0.0_dp
    end if

    ! Conditional: spin-resolved
    if (optcfg%spin_resolved) then
      allocate(oe%alpha_te_up(oe%nE))
      allocate(oe%alpha_te_dw(oe%nE))
      allocate(oe%alpha_tm_up(oe%nE))
      allocate(oe%alpha_tm_dw(oe%nE))
      oe%alpha_te_up = 0.0_dp
      oe%alpha_te_dw = 0.0_dp
      oe%alpha_tm_up = 0.0_dp
      oe%alpha_tm_dw = 0.0_dp
    end if

    ! Reset gain state
    oe%gain_fermi_computed = .false.
    oe%mu_e = 0.0_dp
    oe%mu_h = 0.0_dp

  end subroutine optics_init


  ! ------------------------------------------------------------------
  ! Deallocate all arrays in the optics_engine and reset scalars.
  ! Safe to call multiple times (double-free safe).
  ! ------------------------------------------------------------------
  subroutine optics_free(oe)
    type(optics_engine), intent(inout) :: oe

    if (oe%was_freed) return
    oe%was_freed = .true.
    if (allocated(oe%E_grid))         deallocate(oe%E_grid)
    if (allocated(oe%alpha_te))       deallocate(oe%alpha_te)
    if (allocated(oe%alpha_tm))       deallocate(oe%alpha_tm)
    if (allocated(oe%alpha_isbt))     deallocate(oe%alpha_isbt)
    if (allocated(oe%alpha_gain_te))  deallocate(oe%alpha_gain_te)
    if (allocated(oe%alpha_gain_tm))  deallocate(oe%alpha_gain_tm)
    if (allocated(oe%spont_te))       deallocate(oe%spont_te)
    if (allocated(oe%spont_tm))       deallocate(oe%spont_tm)
    if (allocated(oe%alpha_te_up))    deallocate(oe%alpha_te_up)
    if (allocated(oe%alpha_te_dw))    deallocate(oe%alpha_te_dw)
    if (allocated(oe%alpha_tm_up))    deallocate(oe%alpha_tm_up)
    if (allocated(oe%alpha_tm_dw))    deallocate(oe%alpha_tm_dw)

    oe%nE = 0
    oe%gain_fermi_computed = .false.
    oe%mu_e = 0.0_dp
    oe%mu_h = 0.0_dp

  end subroutine optics_free


  ! ------------------------------------------------------------------
  ! Type-bound procedure wrapper for free (delegates to free subroutine)
  ! ------------------------------------------------------------------
  subroutine optics_free_tbp(oe)
    class(optics_engine), intent(inout) :: oe
    call optics_free(oe)
  end subroutine optics_free_tbp


  ! ------------------------------------------------------------------
  ! Finalizer: delegates to optics_free for automatic cleanup
  ! when an optics_engine variable goes out of scope.
  ! ------------------------------------------------------------------
  subroutine optics_engine_finalize(oe)
    type(optics_engine), intent(inout) :: oe
    call optics_free(oe)
  end subroutine optics_engine_finalize


  ! ==================================================================
  ! optics_engine accumulate/finalize routines
  !
  ! All state is read/written through the optics_engine type.
  ! ==================================================================

  ! ------------------------------------------------------------------
  ! Accumulate absorption contributions from one k_par point.
  ! Reads config from oe%config instead of optcfg argument.
  ! Reads/writes oe%E_grid, oe%alpha_te, oe%alpha_tm, etc.
  ! ------------------------------------------------------------------
  subroutine optics_accumulate(oe, eigvals, eigvecs, k_weight, &
    & vel, numcb, numvb, fermi_level)

    type(optics_engine), intent(inout) :: oe
    real(kind=dp), intent(in), contiguous :: eigvals(:)
    complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)
    real(kind=dp), intent(in) :: k_weight
    type(csr_matrix), intent(in) :: vel(3)
    integer, intent(in) :: numcb, numvb
    real(kind=dp), intent(in) :: fermi_level

    integer :: i, j, dir, ie, dim, Ngrid
    real(kind=dp) :: dE, f_c, f_v, occ_factor
    real(kind=dp) :: px, py, pz
    real(kind=dp) :: gamma_l, gamma_g
    real(kind=dp) :: w_up_i, w_dw_i, w_up_j, w_dw_j
    complex(kind=dp) :: Pele
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)

    if (oe%nE == 0) return

    dim = size(eigvecs, 1)
    Ngrid = dim / 8
    allocate(Ytmp(dim))

    ! Half-widths at half-maximum from FWHM
    gamma_l = oe%config%linewidth_lorentzian / 2.0_dp
    gamma_g = oe%config%linewidth_gaussian / 2.0_dp

    do i = 1, numvb
      f_v = fermi_dirac(eigvals(i), fermi_level, oe%config%temperature)

      do j = 1, numcb
        f_c = fermi_dirac(eigvals(numvb + j), fermi_level, oe%config%temperature)

        occ_factor = f_v - f_c
        if (occ_factor < 1.0e-30_dp) cycle

        dE = eigvals(numvb + j) - eigvals(i)
        if (dE < DE_MIN) cycle

        ! Velocity matrix elements
        px = 0.0_dp
        py = 0.0_dp
        pz = 0.0_dp

        do dir = 1, 3
          call csr_spmv(vel(dir), eigvecs(:,i), Ytmp, ONE, ZERO)
          Pele = zdotc(dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
          select case(dir)
          case(1)
            px = real(Pele * conjg(Pele), kind=dp)
          case(2)
            py = real(Pele * conjg(Pele), kind=dp)
          case(3)
            pz = real(Pele * conjg(Pele), kind=dp)
          end select
        end do

        ! Broaden and accumulate onto the energy grid
        if (oe%config%confinement == 'bulk') then
          ! Bulk: TE = TM = (px+py+pz)/2
          do ie = 1, oe%nE
            oe%alpha_te(ie) = oe%alpha_te(ie) + occ_factor * (px + py + pz) * 0.5_dp &
              & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            oe%alpha_tm(ie) = oe%alpha_tm(ie) + occ_factor * (px + py + pz) * 0.5_dp &
              & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          end do
        else
          ! QW/wire: TE = px+py, TM = pz
          do ie = 1, oe%nE
            oe%alpha_te(ie) = oe%alpha_te(ie) + occ_factor * (px + py) &
              & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            oe%alpha_tm(ie) = oe%alpha_tm(ie) + occ_factor * pz &
              & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          end do
        end if

        ! Spin-resolved accumulation
        if (oe%config%spin_resolved) then
          call spin_weights(eigvecs(:,numvb+j), Ngrid, w_up_j, w_dw_j)
          call spin_weights(eigvecs(:,i), Ngrid, w_up_i, w_dw_i)
          do ie = 1, oe%nE
            oe%alpha_te_up(ie) = oe%alpha_te_up(ie) + occ_factor * (px + py) &
              & * w_up_i * w_up_j * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            oe%alpha_tm_up(ie) = oe%alpha_tm_up(ie) + occ_factor * pz &
              & * w_up_i * w_up_j * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            oe%alpha_te_dw(ie) = oe%alpha_te_dw(ie) + occ_factor * (px + py) &
              & * w_dw_i * w_dw_j * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            oe%alpha_tm_dw(ie) = oe%alpha_tm_dw(ie) + occ_factor * pz &
              & * w_dw_i * w_dw_j * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          end do
        end if

      end do
    end do

    deallocate(Ytmp)

  end subroutine optics_accumulate


  ! ------------------------------------------------------------------
  ! Accumulate spontaneous emission contributions from one k_par point.
  ! ------------------------------------------------------------------
  subroutine optics_accumulate_spontaneous(oe, eigvals, eigvecs, k_weight, &
    & vel, numcb, numvb, fermi_level)

    type(optics_engine), intent(inout) :: oe
    real(kind=dp), intent(in), contiguous :: eigvals(:)
    complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)
    real(kind=dp), intent(in) :: k_weight
    type(csr_matrix), intent(in) :: vel(3)
    integer, intent(in) :: numcb, numvb
    real(kind=dp), intent(in) :: fermi_level

    integer :: i, j, dir, ie, dim
    real(kind=dp) :: dE, f_c, f_v, occ_factor
    real(kind=dp) :: px, py, pz
    real(kind=dp) :: gamma_l, gamma_g
    complex(kind=dp) :: Pele
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)

    if (oe%nE == 0) return
    if (.not. allocated(oe%spont_te)) return

    dim = size(eigvecs, 1)
    allocate(Ytmp(dim))

    gamma_l = oe%config%linewidth_lorentzian / 2.0_dp
    gamma_g = oe%config%linewidth_gaussian / 2.0_dp

    do i = 1, numvb
      f_v = fermi_dirac(eigvals(i), fermi_level, oe%config%temperature)

      do j = 1, numcb
        f_c = fermi_dirac(eigvals(numvb + j), fermi_level, oe%config%temperature)

        ! Spontaneous emission: f_c * (1 - f_v)
        occ_factor = f_c * (1.0_dp - f_v)
        if (occ_factor < 1.0e-30_dp) cycle

        dE = eigvals(numvb + j) - eigvals(i)
        if (dE < DE_MIN) cycle

        px = 0.0_dp
        py = 0.0_dp
        pz = 0.0_dp

        do dir = 1, 3
          call csr_spmv(vel(dir), eigvecs(:,i), Ytmp, ONE, ZERO)
          Pele = zdotc(dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
          select case(dir)
          case(1)
            px = real(Pele * conjg(Pele), kind=dp)
          case(2)
            py = real(Pele * conjg(Pele), kind=dp)
          case(3)
            pz = real(Pele * conjg(Pele), kind=dp)
          end select
        end do

        do ie = 1, oe%nE
          oe%spont_te(ie) = oe%spont_te(ie) + occ_factor * (px + py) &
            & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          oe%spont_tm(ie) = oe%spont_tm(ie) + occ_factor * pz &
            & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

    deallocate(Ytmp)

  end subroutine optics_accumulate_spontaneous


  ! ------------------------------------------------------------------
  ! Gain spectrum for a QW with population inversion (engine version).
  ! ------------------------------------------------------------------
  subroutine compute_gain_qw(oe, eigvals, eigvecs, k_weight, &
    & vel, numcb, numvb, carrier_density)

    type(optics_engine), intent(inout) :: oe
    real(kind=dp), intent(in), contiguous :: eigvals(:)
    complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)
    real(kind=dp), intent(in) :: k_weight
    type(csr_matrix), intent(in) :: vel(3)
    integer, intent(in) :: numcb, numvb
    real(kind=dp), intent(in) :: carrier_density

    integer :: i, j, dir, ie, dim
    real(kind=dp) :: dE, f_c, f_v, occ_factor
    real(kind=dp) :: px, py, pz
    real(kind=dp) :: gamma_l, gamma_g
    complex(kind=dp) :: Pele
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)

    if (oe%nE == 0) return
    if (.not. allocated(oe%alpha_gain_te)) return

    dim = size(eigvecs, 1)
    allocate(Ytmp(dim))

    ! Compute quasi-Fermi levels once on the first k-point call
    if (.not. oe%gain_fermi_computed) then
      oe%mu_e = find_quasi_fermi(eigvals, oe%config%temperature, &
        & carrier_density, numcb, numvb)
      oe%mu_h = find_quasi_fermi_holes(eigvals, oe%config%temperature, &
        & carrier_density, numvb)
      oe%gain_fermi_computed = .true.
      print '(a,es14.6,a)', '  Gain: mu_e = ', oe%mu_e, ' eV'
      print '(a,es14.6,a)', '  Gain: mu_h = ', oe%mu_h, ' eV'
      print '(a,es10.2,a)', '  Gain: carrier density = ', carrier_density, ' cm^-2'
    end if

    gamma_l = oe%config%linewidth_lorentzian / 2.0_dp
    gamma_g = oe%config%linewidth_gaussian / 2.0_dp

    do i = 1, numvb
      f_v = fermi_dirac(eigvals(i), oe%mu_h, oe%config%temperature)

      do j = 1, numcb
        f_c = fermi_dirac(eigvals(numvb + j), oe%mu_e, oe%config%temperature)

        occ_factor = f_v - f_c

        dE = eigvals(numvb + j) - eigvals(i)
        if (dE < DE_MIN) cycle

        if (abs(occ_factor) < 1.0e-30_dp) cycle

        px = 0.0_dp
        py = 0.0_dp
        pz = 0.0_dp

        do dir = 1, 3
          call csr_spmv(vel(dir), eigvecs(:,i), Ytmp, ONE, ZERO)
          Pele = zdotc(dim, eigvecs(1:dim,numvb+j), 1, Ytmp, 1)
          select case(dir)
          case(1)
            px = real(Pele * conjg(Pele), kind=dp)
          case(2)
            py = real(Pele * conjg(Pele), kind=dp)
          case(3)
            pz = real(Pele * conjg(Pele), kind=dp)
          end select
        end do

        do ie = 1, oe%nE
          oe%alpha_gain_te(ie) = oe%alpha_gain_te(ie) + occ_factor * (px + py) &
            & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          oe%alpha_gain_tm(ie) = oe%alpha_gain_tm(ie) + occ_factor * pz &
            & * lineshape_voigt(oe%E_grid(ie), dE, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

    deallocate(Ytmp)

  end subroutine compute_gain_qw


  ! ------------------------------------------------------------------
  ! ISBT absorption spectrum: TM-polarized (engine version).
  ! ------------------------------------------------------------------
  subroutine compute_isbt_absorption(oe, eigvals, eigvecs, &
    & vel, numcb, numvb, k_weight, fermi_level)

    type(optics_engine), intent(inout) :: oe
    real(kind=dp), intent(in), contiguous :: eigvals(:)
    complex(kind=dp), intent(in), contiguous :: eigvecs(:,:)
    type(csr_matrix), intent(in) :: vel(3)
    integer, intent(in)          :: numcb, numvb
    real(kind=dp), intent(in)    :: k_weight
    real(kind=dp), intent(in)    :: fermi_level

    integer :: i, j, ie, state_i, state_j, dim
    integer, parameter :: dir_isbt = 3
    real(kind=dp) :: E_ij, occ_factor, f_i, f_j
    real(kind=dp) :: gamma_l, gamma_g
    real(kind=dp) :: p_abs2
    complex(kind=dp) :: pele_ij
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)

    real(kind=dp) :: ef_cb

    if (oe%nE == 0) return

    dim = size(eigvecs, 1)
    allocate(Ytmp(dim))

    gamma_l = oe%config%linewidth_lorentzian / 2.0_dp
    gamma_g = oe%config%linewidth_gaussian / 2.0_dp

    if (oe%gain_fermi_computed) then
      ef_cb = oe%mu_e
    else
      ef_cb = fermi_level
    end if

    do i = 1, numcb
      state_i = numvb + i
      f_i = fermi_dirac(eigvals(state_i), ef_cb, oe%config%temperature)

      do j = i + 1, numcb
        state_j = numvb + j
        f_j = fermi_dirac(eigvals(state_j), ef_cb, oe%config%temperature)

        occ_factor = f_i - f_j
        if (occ_factor < 1.0e-30_dp) cycle

        E_ij = eigvals(state_j) - eigvals(state_i)
        if (E_ij < DE_MIN) cycle

        call csr_spmv(vel(dir_isbt), eigvecs(:,state_j), Ytmp, ONE, ZERO)
        pele_ij = zdotc(dim, eigvecs(1:dim,state_i), 1, Ytmp, 1)
        p_abs2 = real(pele_ij * conjg(pele_ij), kind=dp)

        do ie = 1, oe%nE
          oe%alpha_isbt(ie) = oe%alpha_isbt(ie) + occ_factor * p_abs2 &
            & * lineshape_voigt(oe%E_grid(ie), E_ij, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

    deallocate(Ytmp)

  end subroutine compute_isbt_absorption


  ! ------------------------------------------------------------------
  ! Reset the gain quasi-Fermi level state (engine version).
  ! ------------------------------------------------------------------
  subroutine gain_reset(oe)
    type(optics_engine), intent(inout) :: oe

    oe%gain_fermi_computed = .false.
    oe%mu_e = 0.0_dp
    oe%mu_h = 0.0_dp
  end subroutine gain_reset


  ! ------------------------------------------------------------------
  ! Apply the physical prefactor to all accumulation arrays.
  ! No file I/O.  Reads refractive_index from oe%config.
  ! ------------------------------------------------------------------
  subroutine optics_apply_prefactor(oe)
    type(optics_engine), intent(inout) :: oe

    integer :: ie
    real(kind=dp) :: prefactor_E

    ! eps0 in AA units
    real(kind=dp), parameter :: e0_AA = e0 / 10.0_dp
    real(kind=dp), parameter :: AA_TO_CM = 1.0e8_dp

    if (oe%nE == 0) return

    ! Apply prefactor to unconditional arrays
    do concurrent (ie = 1:oe%nE)
      if (oe%E_grid(ie) > 0.0_dp) then
        prefactor_E = (2.0_dp * pi_dp * e**2) &
          & / (oe%config%refractive_index * c * e0_AA * hbar**2 * oe%E_grid(ie))
      else
        prefactor_E = 0.0_dp
      end if
      oe%alpha_te(ie) = prefactor_E * oe%alpha_te(ie) * AA_TO_CM
      oe%alpha_tm(ie) = prefactor_E * oe%alpha_tm(ie) * AA_TO_CM
      oe%alpha_isbt(ie) = prefactor_E * oe%alpha_isbt(ie) * AA_TO_CM
    end do

    ! Apply prefactor to gain arrays
    if (oe%config%gain_enabled .and. allocated(oe%alpha_gain_te)) then
      do concurrent (ie = 1:oe%nE)
        if (oe%E_grid(ie) > 0.0_dp) then
          prefactor_E = (2.0_dp * pi_dp * e**2) &
            & / (oe%config%refractive_index * c * e0_AA * hbar**2 * oe%E_grid(ie))
        else
          prefactor_E = 0.0_dp
        end if
        oe%alpha_gain_te(ie) = prefactor_E * oe%alpha_gain_te(ie) * AA_TO_CM
        oe%alpha_gain_tm(ie) = prefactor_E * oe%alpha_gain_tm(ie) * AA_TO_CM
      end do
    end if

    ! Apply prefactor to spontaneous emission arrays
    if (oe%config%spontaneous_enabled .and. allocated(oe%spont_te)) then
      do concurrent (ie = 1:oe%nE)
        if (oe%E_grid(ie) > 0.0_dp) then
          prefactor_E = (2.0_dp * pi_dp * e**2) &
            & / (oe%config%refractive_index * c * e0_AA * hbar**2 * oe%E_grid(ie))
        else
          prefactor_E = 0.0_dp
        end if
        oe%spont_te(ie) = prefactor_E * oe%spont_te(ie) * AA_TO_CM
        oe%spont_tm(ie) = prefactor_E * oe%spont_tm(ie) * AA_TO_CM
      end do
    end if

    ! Apply prefactor to spin-resolved arrays
    if (oe%config%spin_resolved .and. allocated(oe%alpha_te_up)) then
      do concurrent (ie = 1:oe%nE)
        if (oe%E_grid(ie) > 0.0_dp) then
          prefactor_E = (2.0_dp * pi_dp * e**2) &
            & / (oe%config%refractive_index * c * e0_AA * hbar**2 * oe%E_grid(ie))
        else
          prefactor_E = 0.0_dp
        end if
        oe%alpha_te_up(ie) = prefactor_E * oe%alpha_te_up(ie) * AA_TO_CM
        oe%alpha_te_dw(ie) = prefactor_E * oe%alpha_te_dw(ie) * AA_TO_CM
        oe%alpha_tm_up(ie) = prefactor_E * oe%alpha_tm_up(ie) * AA_TO_CM
        oe%alpha_tm_dw(ie) = prefactor_E * oe%alpha_tm_dw(ie) * AA_TO_CM
      end do
    end if

  end subroutine optics_apply_prefactor


  ! ------------------------------------------------------------------
  ! Write all spectra files from oe% components.
  ! Uses oe%config for conditional output (gain, spontaneous, spin-resolved).
  ! ------------------------------------------------------------------
  subroutine optics_write_output(oe)
    type(optics_engine), intent(in) :: oe

    integer :: ie, iounit

    if (oe%nE == 0) return

    call ensure_output_dir()
    call get_unit(iounit)

    ! --- TE absorption ---
    open(unit=iounit, file='output/absorption_TE.dat', status='replace', &
      & action='write')
    write(iounit, '(a)') '# Interband TE absorption spectrum (alpha_TE vs E)'
    write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
    do ie = 1, oe%nE
      write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_te(ie)
    end do
    close(iounit)

    ! --- TM absorption ---
    open(unit=iounit, file='output/absorption_TM.dat', status='replace', &
      & action='write')
    write(iounit, '(a)') '# Interband TM absorption spectrum (alpha_TM vs E)'
    write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
    do ie = 1, oe%nE
      write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_tm(ie)
    end do
    close(iounit)

    ! --- ISBT absorption ---
    open(unit=iounit, file='output/absorption_ISBT.dat', status='replace', &
      & action='write')
    write(iounit, '(a)') '# ISBT absorption spectrum (TM-polarized, z-dipole only)'
    write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
    do ie = 1, oe%nE
      write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_isbt(ie)
    end do
    close(iounit)

    print '(a)', 'Optical spectra written to output/absorption_TE.dat, output/absorption_TM.dat, and output/absorption_ISBT.dat'

    ! --- Gain output ---
    if (oe%config%gain_enabled .and. allocated(oe%alpha_gain_te)) then
      open(unit=iounit, file='output/gain_TE.dat', status='replace', &
        & action='write')
      write(iounit, '(a)') '# Interband TE gain spectrum (gain_TE vs E)'
      write(iounit, '(a,es10.2,a)') '# Carrier density = ', &
        & oe%config%gain_carrier_density, ' cm^-2'
      write(iounit, '(a)') '# E(eV)  gain(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_gain_te(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/gain_TM.dat', status='replace', &
        & action='write')
      write(iounit, '(a)') '# Interband TM gain spectrum (gain_TM vs E)'
      write(iounit, '(a,es10.2,a)') '# Carrier density = ', &
        & oe%config%gain_carrier_density, ' cm^-2'
      write(iounit, '(a)') '# E(eV)  gain(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_gain_tm(ie)
      end do
      close(iounit)

      print '(a)', 'Gain spectra written to output/gain_TE.dat and output/gain_TM.dat'
    end if

    ! --- Spontaneous emission output ---
    if (oe%config%spontaneous_enabled .and. allocated(oe%spont_te)) then
      open(unit=iounit, file='output/spontaneous_TE.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spontaneous emission TE spectrum'
      write(iounit, '(a)') '# E(eV)  rate(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%spont_te(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/spontaneous_TM.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spontaneous emission TM spectrum'
      write(iounit, '(a)') '# E(eV)  rate(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%spont_tm(ie)
      end do
      close(iounit)

      print '(a)', 'Spontaneous emission written to output/spontaneous_TE.dat and output/spontaneous_TM.dat'
    end if

    ! --- Spin-resolved absorption output ---
    if (oe%config%spin_resolved .and. allocated(oe%alpha_te_up)) then
      open(unit=iounit, file='output/absorption_TE_up.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-up TE absorption spectrum (alpha_TE_up vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_te_up(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/absorption_TE_dw.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-down TE absorption spectrum (alpha_TE_dw vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_te_dw(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/absorption_TM_up.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-up TM absorption spectrum (alpha_TM_up vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_tm_up(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/absorption_TM_dw.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-down TM absorption spectrum (alpha_TM_dw vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, oe%nE
        write(iounit, '(es16.8, 2x, es16.8)') oe%E_grid(ie), oe%alpha_tm_dw(ie)
      end do
      close(iounit)

      print '(a)', 'Spin-resolved spectra written to output/absorption_TE_up.dat, output/absorption_TE_dw.dat, output/absorption_TM_up.dat, output/absorption_TM_dw.dat'
    end if

  end subroutine optics_write_output

end module optical_spectra
