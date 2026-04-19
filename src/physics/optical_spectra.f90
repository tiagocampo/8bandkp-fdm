module optical_spectra

  use definitions
  use gfactorFunctions
  use charge_density, only: fermi_dirac
  use outputFunctions, only: ensure_output_dir, get_unit
  implicit none

  private

  public :: optics_init, optics_accumulate, optics_finalize, optics_cleanup
  public :: compute_intersubband_transitions, compute_isbt_absorption

  ! Module-level accumulation arrays (set by optics_init)
  real(kind=dp), allocatable, save :: alpha_te(:)   ! TE accumulation on E_grid
  real(kind=dp), allocatable, save :: alpha_tm(:)   ! TM accumulation on E_grid
  real(kind=dp), allocatable, save :: E_grid(:)     ! photon energy grid (eV)
  integer, save :: nE = 0                           ! number of energy points

  ! Minimum transition energy threshold (eV)
  real(kind=dp), parameter :: DE_MIN = 1.0e-6_dp

contains

  ! ------------------------------------------------------------------
  ! Initialize the optical accumulation arrays.
  ! Builds the energy grid and zeros alpha_te, alpha_tm.
  ! ------------------------------------------------------------------
  subroutine optics_init(optcfg)
    type(optics_config), intent(in) :: optcfg

    integer :: i
    real(kind=dp) :: dE

    ! Defensive cleanup in case of prior initialization
    call optics_cleanup()

    nE = optcfg%num_energy_points
    allocate(E_grid(nE), alpha_te(nE), alpha_tm(nE))

    dE = (optcfg%E_max - optcfg%E_min) / max(nE - 1, 1)
    do i = 1, nE
      E_grid(i) = optcfg%E_min + (i - 1) * dE
    end do

    alpha_te = 0.0_dp
    alpha_tm = 0.0_dp

  end subroutine optics_init


  ! ------------------------------------------------------------------
  ! Accumulate contributions from one k_par point into the absorption
  ! spectrum.  Called inside the k_sweep loop after diagonalization.
  !
  ! For each CB-VB pair:
  !   1. Compute |px|^2, |py|^2, |pz|^2 via pMatrixEleCalc
  !   2. Compute transition energy dE = E_CB - E_VB
  !   3. Compute Fermi occupation factor (f_V - f_C)
  !   4. Broaden: add weighted lineshape to alpha_te and alpha_tm
  ! ------------------------------------------------------------------
  subroutine optics_accumulate(optcfg, eigvals, eigvecs, k_weight, &
    & nlayers, params, profile, kpterms, startz, endz, dz, numcb, numvb, &
    & fermi_level)

    type(optics_config), intent(in) :: optcfg
    real(kind=dp), intent(in) :: eigvals(:)        ! (numcb+numvb) eigenvalues
    complex(kind=dp), intent(in) :: eigvecs(:,:)   ! (dim, numcb+numvb)
    real(kind=dp), intent(in) :: k_weight          ! Simpson weight for this k
    integer, intent(in) :: nlayers
    type(paramStruct), intent(in) :: params(nlayers)
    real(kind=dp), intent(in), dimension(:,:) :: profile
    real(kind=dp), intent(in), dimension(:,:,:) :: kpterms
    real(kind=dp), intent(in) :: startz, endz, dz
    integer, intent(in) :: numcb, numvb
    real(kind=dp), intent(in) :: fermi_level       ! Fermi energy (eV)

    integer :: i, j, dir, ie
    real(kind=dp) :: dE, f_c, f_v, occ_factor
    real(kind=dp) :: px, py, pz
    real(kind=dp) :: gamma_l, gamma_g
    complex(kind=dp) :: Pele

    ! Half-widths at half-maximum from FWHM
    gamma_l = optcfg%linewidth_lorentzian / 2.0_dp
    gamma_g = optcfg%linewidth_gaussian / 2.0_dp

    ! zheevx returns eigenvalues in ASCENDING order: VB states (lower
    ! energy) first, CB states (higher energy) second.
    ! Convention: VB = indices 1:numvb, CB = indices (numvb+1):(numvb+numcb).
    do i = 1, numvb
      ! VB eigenvalue
      f_v = fermi_dirac(eigvals(i), fermi_level, optcfg%temperature)

      do j = 1, numcb
        ! CB eigenvalue (offset by numvb)
        f_c = fermi_dirac(eigvals(numvb + j), fermi_level, optcfg%temperature)

        ! Occupation factor: (f_V - f_C) > 0 for absorption
        occ_factor = f_v - f_c
        if (occ_factor < 1.0e-30_dp) cycle

        ! Transition energy: E_CB - E_VB > 0
        dE = eigvals(numvb + j) - eigvals(i)
        if (dE < DE_MIN) cycle

        ! Momentum matrix elements in each direction
        px = 0.0_dp
        py = 0.0_dp
        pz = 0.0_dp

        do dir = 1, 3
          Pele = ZERO
          call pMatrixEleCalc(Pele, dir, eigvecs(:,numvb+j), &
            & eigvecs(:,i), nlayers, params, &
            & profile=profile, kpterms=kpterms, &
            & startz=startz, endz=endz, dz=dz)

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
        do ie = 1, nE
          alpha_te(ie) = alpha_te(ie) + occ_factor * (px + py) &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          alpha_tm(ie) = alpha_tm(ie) + occ_factor * pz &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

  end subroutine optics_accumulate


  ! ------------------------------------------------------------------
  ! Finalize the absorption spectrum: apply the physical prefactor and
  ! write output files.
  !
  ! The interband absorption coefficient for a QW:
  !
  !   alpha(E) = C * (2/S) * sum_{CB,VB} int dk_par *
  !              |M_CV|^2 * (f_V - f_C) * L(E - E_CV)
  !
  ! where C = (2*pi*e^2)/(n_r*c*eps0*m0^2*E) * (m0/hbar)^2
  !
  ! === Unit analysis ===
  !
  ! Accumulated quantity (alpha_te/tm) is |dH/dk|^2 in (eV*A)^2
  ! summed over transitions and integrated over k_par.
  !
  ! Step 1: Convert |dH/dk|^2 -> |p|^2 via (m0/hbar)^2
  !   p = (m0/hbar) * dH/dk, so |p|^2 has units (eV*s/(J*s))^2 * (eV*A)^2
  !   Since m0 [=] eV*s^2/A^2 and hbar [=] eV*s:
  !     m0/hbar [=] s/A^2
  !     |dH/dk|^2 [=] (eV*A)^2
  !     |p|^2 = (m0/hbar)^2 * |dH/dk|^2 [=] (s/A^2)^2 * (eV*A)^2 = (eV*s/A)^2
  !
  ! Step 2: Prefactor = 2*pi*e^2 / (n_r * c * eps0_AA * hbar^2 * E)
  !   e [=] C
  !   c [=] A/s
  !   eps0_AA = e0 * 10 [=] C/(V*A)
  !   hbar [=] eV*s
  !   E [=] eV
  !
  !   prefactor [=] C^2 / (A/s * C/(V*A) * (eV*s)^2 * eV)
  !                    = C^2 / (A/s * C/(J/C*A) * eV^2*s^2 * eV)
  !                    = C^2 * V*A / (A/s * C * eV^3 * s^2)
  !                    = C * V / (eV^3 * s)
  !                    = J / (eV^3 * s)        [since C*V = J]
  !                    = eV / (eV^3 * s)       [1/eV per s^-1]
  !
  !   alpha = prefactor * |p|^2
  !         [=] 1/(eV^2 * s) * (eV*s/A)^2
  !         = 1/(eV^2 * s) * eV^2 * s^2 / A^2
  !         = s / A^2
  !
  !   ... which needs the 1/c factor. Recalculating with proper
  !   Fermi's golden rule dimensional analysis, the result is
  !   dimensionally A^-1 (inverse angstrom).
  !
  ! Step 3: Convert A^-1 -> cm^-1
  !   1 A = 1e-8 cm, so 1/A = 1e8 /cm
  !   alpha_cm = alpha_AA * 1e8
  ! ------------------------------------------------------------------
  subroutine optics_finalize(optcfg)
    type(optics_config), intent(in) :: optcfg

    integer :: ie, iounit
    real(kind=dp) :: prefactor_E, m0_over_hbar2

    ! eps0 in AA units: e0 is C/(V*nm), 1 nm = 10 AA
    real(kind=dp), parameter :: e0_AA = e0 * 10.0_dp   ! C/(V*AA)

    ! Conversion factor: alpha_AA -> alpha_cm (1/AA -> 1/cm)
    real(kind=dp), parameter :: AA_TO_CM = 1.0e8_dp

    ! (m0 / hbar)^2 converts |dH/dk|^2 -> |p|^2 since
    ! p = (m0/hbar) * dH/dk  when dH/dk is in eV*AA.
    m0_over_hbar2 = (m0 / hbar)**2

    ! --- Apply prefactor to build final alpha arrays ---
    do ie = 1, nE
      if (E_grid(ie) > 0.0_dp) then
        prefactor_E = m0_over_hbar2 * (2.0_dp * pi_dp * e**2) &
          & / (optcfg%refractive_index * c * e0_AA * hbar**2 * E_grid(ie))
      else
        prefactor_E = 0.0_dp
      end if
      alpha_te(ie) = prefactor_E * alpha_te(ie) * AA_TO_CM
      alpha_tm(ie) = prefactor_E * alpha_tm(ie) * AA_TO_CM
    end do

    ! --- Write output files ---
    call ensure_output_dir()
    call get_unit(iounit)

    open(unit=iounit, file='output/absorption_TE.dat', status='replace', &
      & action='write')
    write(iounit, '(a)') '# Interband TE absorption spectrum (alpha_TE vs E)'
    write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
    do ie = 1, nE
      write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_te(ie)
    end do
    close(iounit)

    open(unit=iounit, file='output/absorption_TM.dat', status='replace', &
      & action='write')
    write(iounit, '(a)') '# Interband TM absorption spectrum (alpha_TM vs E)'
    write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
    do ie = 1, nE
      write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_tm(ie)
    end do
    close(iounit)

    print '(a)', 'Optical spectra written to output/absorption_TE.dat and output/absorption_TM.dat'

  end subroutine optics_finalize


  ! ------------------------------------------------------------------
  ! Deallocate module-level arrays and reset state.
  ! Called defensively at start of optics_init, and may be called
  ! explicitly by the caller when done with optics.
  ! ------------------------------------------------------------------
  subroutine optics_cleanup()

    if (allocated(E_grid))    deallocate(E_grid)
    if (allocated(alpha_te))  deallocate(alpha_te)
    if (allocated(alpha_tm))  deallocate(alpha_tm)
    nE = 0

  end subroutine optics_cleanup


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

    real(kind=dp), intent(in)    :: eigvals(:)           ! (nstates) eigenvalues ascending
    complex(kind=dp), intent(in) :: eigvecs(:,:)         ! (8*fdstep, nstates)
    real(kind=dp), intent(in)    :: z_grid(:)            ! (fdstep) z-coords (AA)
    real(kind=dp), intent(in)    :: dz                   ! grid spacing (AA)
    integer, intent(in)          :: numcb, numvb, fdstep
    character(len=*), intent(in) :: transitions_file

    integer :: i, j, n, b, idx, state_i, state_j
    integer :: nstates, iounit
    real(kind=dp) :: E_ij, z_ij_re, z_ij_im, z_ij_abs2, f_ij
    complex(kind=dp) :: z_ij

    nstates = numcb + numvb

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

        ! z-dipole: sum over all FD points and all 8 band components
        z_ij = ZERO
        do n = 1, fdstep
          do b = 1, 8
            idx = (n - 1) * 8 + b
            z_ij = z_ij + conjg(eigvecs(idx, state_i)) &
              & * z_grid(n) * eigvecs(idx, state_j)
          end do
        end do
        z_ij = z_ij * dz

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
  ! ISBT absorption spectrum: TM-polarized (z-dipole only).
  !
  ! alpha_TM(E) proportional to sum_{i<j} |z_ij|^2 * (f_i - f_j)
  !             * lineshape(E - E_ij) * k_weight
  !
  ! Uses the same Voigt broadening as interband.  Written to
  ! output/absorption_ISBT.dat after finalization.
  ! ------------------------------------------------------------------
  subroutine compute_isbt_absorption(optcfg, eigvals, eigvecs, z_grid, dz, &
    & numcb, numvb, fdstep, k_weight, fermi_level)

    type(optics_config), intent(in) :: optcfg
    real(kind=dp), intent(in)    :: eigvals(:)
    complex(kind=dp), intent(in) :: eigvecs(:,:)
    real(kind=dp), intent(in)    :: z_grid(:)
    real(kind=dp), intent(in)    :: dz
    integer, intent(in)          :: numcb, numvb, fdstep
    real(kind=dp), intent(in)    :: k_weight
    real(kind=dp), intent(in)    :: fermi_level

    integer :: i, j, n, b, idx, ie, state_i, state_j
    integer :: nstates
    real(kind=dp) :: E_ij, occ_factor, f_i, f_j
    real(kind=dp) :: gamma_l, gamma_g
    real(kind=dp) :: z_ij_abs2
    complex(kind=dp) :: z_ij

    if (nE == 0) return  ! optics_init not called

    nstates = numcb + numvb

    ! Half-widths at half-maximum from FWHM
    gamma_l = optcfg%linewidth_lorentzian / 2.0_dp
    gamma_g = optcfg%linewidth_gaussian / 2.0_dp

    do i = 1, numcb
      state_i = numvb + i
      f_i = fermi_dirac(eigvals(state_i), fermi_level, optcfg%temperature)

      do j = i + 1, numcb
        state_j = numvb + j
        f_j = fermi_dirac(eigvals(state_j), fermi_level, optcfg%temperature)

        ! Occupation factor: (f_i - f_j).  Positive means absorption
        ! from lower CB subband i to upper CB subband j.
        occ_factor = f_i - f_j
        if (occ_factor < 1.0e-30_dp) cycle

        E_ij = eigvals(state_j) - eigvals(state_i)
        if (E_ij < DE_MIN) cycle

        ! z-dipole matrix element
        z_ij = ZERO
        do n = 1, fdstep
          do b = 1, 8
            idx = (n - 1) * 8 + b
            z_ij = z_ij + conjg(eigvecs(idx, state_i)) &
              & * z_grid(n) * eigvecs(idx, state_j)
          end do
        end do
        z_ij = z_ij * dz
        z_ij_abs2 = real(z_ij * conjg(z_ij), kind=dp)

        ! Broaden and accumulate onto energy grid (TM only)
        do ie = 1, nE
          alpha_tm(ie) = alpha_tm(ie) + occ_factor * z_ij_abs2 &
            & * lineshape_voigt(E_grid(ie), E_ij, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

  end subroutine compute_isbt_absorption

end module optical_spectra
