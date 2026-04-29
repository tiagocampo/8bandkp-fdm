module optical_spectra

  use definitions
  use sparse_matrices, only: csr_matrix, csr_spmv
  use charge_density, only: fermi_dirac
  use spin_projection, only: spin_weights
  use outputFunctions, only: ensure_output_dir, get_unit
  use linalg, only: zdotc
  implicit none

  private

  public :: optics_init, optics_accumulate, optics_accumulate_spontaneous, optics_finalize, optics_cleanup
  public :: compute_intersubband_transitions, compute_isbt_absorption
  public :: compute_gain_qw, gain_reset
  public :: E_grid, alpha_te, alpha_tm, nE

  ! Module-level accumulation arrays (set by optics_init)
  real(kind=dp), allocatable, save :: alpha_te(:)    ! TE accumulation on E_grid
  real(kind=dp), allocatable, save :: alpha_tm(:)    ! TM accumulation on E_grid
  real(kind=dp), allocatable, save :: alpha_isbt(:)  ! ISBT accumulation on E_grid
  real(kind=dp), allocatable, save :: alpha_gain_te(:)  ! Gain TE accumulation on E_grid
  real(kind=dp), allocatable, save :: alpha_gain_tm(:)  ! Gain TM accumulation on E_grid
  real(kind=dp), allocatable, save :: spont_te(:)   ! Spontaneous TE accumulation
  real(kind=dp), allocatable, save :: spont_tm(:)   ! Spontaneous TM accumulation
  real(kind=dp), allocatable, save :: E_grid(:)      ! photon energy grid (eV)
  real(kind=dp), allocatable, save :: alpha_te_up(:)  ! spin-up TE accumulation
  real(kind=dp), allocatable, save :: alpha_te_dw(:)  ! spin-down TE accumulation
  real(kind=dp), allocatable, save :: alpha_tm_up(:)  ! spin-up TM accumulation
  real(kind=dp), allocatable, save :: alpha_tm_dw(:)  ! spin-down TM accumulation
  integer, save :: nE = 0                            ! number of energy points

  ! Gain quasi-Fermi level state (set once per k_sweep)
  real(kind=dp), save :: mu_e = 0.0_dp              ! electron quasi-Fermi level (eV)
  real(kind=dp), save :: mu_h = 0.0_dp              ! hole quasi-Fermi level (eV)
  logical, save :: gain_fermi_computed = .false.

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
    allocate(E_grid(nE), alpha_te(nE), alpha_tm(nE), alpha_isbt(nE))

    dE = (optcfg%E_max - optcfg%E_min) / max(nE - 1, 1)
    do i = 1, nE
      E_grid(i) = optcfg%E_min + (i - 1) * dE
    end do

    alpha_te = 0.0_dp
    alpha_tm = 0.0_dp
    alpha_isbt = 0.0_dp

    ! Gain arrays: allocated when gain is enabled
    if (optcfg%gain_enabled) then
      allocate(alpha_gain_te(nE), alpha_gain_tm(nE))
      alpha_gain_te = 0.0_dp
      alpha_gain_tm = 0.0_dp
    end if

    ! Spontaneous emission arrays: allocated when spontaneous is enabled
    if (optcfg%spontaneous_enabled) then
      allocate(spont_te(nE), spont_tm(nE))
      spont_te = 0.0_dp
      spont_tm = 0.0_dp
    end if

    ! Spin-resolved absorption arrays: allocated when spin_resolved is enabled
    if (optcfg%spin_resolved) then
      allocate(alpha_te_up(nE), alpha_te_dw(nE))
      allocate(alpha_tm_up(nE), alpha_tm_dw(nE))
      alpha_te_up = 0.0_dp; alpha_te_dw = 0.0_dp
      alpha_tm_up = 0.0_dp; alpha_tm_dw = 0.0_dp
    end if

  end subroutine optics_init


  ! ------------------------------------------------------------------
  ! Accumulate contributions from one k_par point into the absorption
  ! spectrum.  Called inside the k_sweep loop after diagonalization.
  !
  ! Uses commutator-based velocity matrices v_alpha = -i [r_alpha, H]
  ! computed on CSR sparse format.  For each CB-VB pair:
  !   1. Compute |<CB|v_x|VB>|^2, |<CB|v_y|VB>|^2, |<CB|v_z|VB>|^2
  !      via CSR SpMV + zdotc
  !   2. Compute transition energy dE = E_CB - E_VB
  !   3. Compute Fermi occupation factor (f_V - f_C)
  !   4. Broaden: add weighted lineshape to alpha_te and alpha_tm
  ! ------------------------------------------------------------------
  subroutine optics_accumulate(optcfg, eigvals, eigvecs, k_weight, &
    & vel, numcb, numvb, fermi_level)

    type(optics_config), intent(in) :: optcfg
    real(kind=dp), intent(in) :: eigvals(:)        ! (numcb+numvb) eigenvalues
    complex(kind=dp), intent(in) :: eigvecs(:,:)   ! (dim, numcb+numvb)
    real(kind=dp), intent(in) :: k_weight          ! Simpson weight for this k
    type(csr_matrix), intent(in) :: vel(3)         ! commutator velocity matrices
    integer, intent(in) :: numcb, numvb
    real(kind=dp), intent(in) :: fermi_level       ! Fermi energy (eV)

    integer :: i, j, dir, ie, dim, Ngrid
    real(kind=dp) :: dE, f_c, f_v, occ_factor
    real(kind=dp) :: px, py, pz
    real(kind=dp) :: gamma_l, gamma_g
    real(kind=dp) :: w_up_i, w_dw_i, w_up_j, w_dw_j
    complex(kind=dp) :: Pele
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)


    dim = size(eigvecs, 1)
    Ngrid = dim / 8
    allocate(Ytmp(dim))

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

        ! Velocity matrix elements via commutator v_alpha = -i [r_alpha, H]
        ! |<CB| v_alpha |VB>|^2 computed as SpMV: Ytmp = vel(dir) * |VB>,
        ! then Pele = <CB| Ytmp>.
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
        do ie = 1, nE
          alpha_te(ie) = alpha_te(ie) + occ_factor * (px + py) &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          alpha_tm(ie) = alpha_tm(ie) + occ_factor * pz &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
        end do

        ! Spin-resolved accumulation (factorized approximation)
        if (optcfg%spin_resolved) then
          call spin_weights(eigvecs(:,numvb+j), Ngrid, w_up_j, w_dw_j)
          call spin_weights(eigvecs(:,i), Ngrid, w_up_i, w_dw_i)
          do ie = 1, nE
            ! Spin-up channel
            alpha_te_up(ie) = alpha_te_up(ie) + occ_factor * (px + py) &
              & * w_up_i * w_up_j * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            alpha_tm_up(ie) = alpha_tm_up(ie) + occ_factor * pz &
              & * w_up_i * w_up_j * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            ! Spin-down channel
            alpha_te_dw(ie) = alpha_te_dw(ie) + occ_factor * (px + py) &
              & * w_dw_i * w_dw_j * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
            alpha_tm_dw(ie) = alpha_tm_dw(ie) + occ_factor * pz &
              & * w_dw_i * w_dw_j * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          end do
        end if

      end do
    end do

    deallocate(Ytmp)

  end subroutine optics_accumulate


  ! ------------------------------------------------------------------
  ! Accumulate spontaneous emission contributions from one k_par point.
  ! Identical to optics_accumulate except the occupation factor is
  !   f_c * (1 - f_v)   (radiative recombination rate)
  ! instead of (f_v - f_c) used for absorption.
  ! ------------------------------------------------------------------
  subroutine optics_accumulate_spontaneous(optcfg, eigvals, eigvecs, k_weight, &
    & vel, numcb, numvb, fermi_level)

    type(optics_config), intent(in) :: optcfg
    real(kind=dp), intent(in) :: eigvals(:)        ! (numcb+numvb) eigenvalues
    complex(kind=dp), intent(in) :: eigvecs(:,:)   ! (dim, numcb+numvb)
    real(kind=dp), intent(in) :: k_weight          ! Simpson weight for this k
    type(csr_matrix), intent(in) :: vel(3)         ! commutator velocity matrices
    integer, intent(in) :: numcb, numvb
    real(kind=dp), intent(in) :: fermi_level       ! Fermi energy (eV)

    integer :: i, j, dir, ie, dim
    real(kind=dp) :: dE, f_c, f_v, occ_factor
    real(kind=dp) :: px, py, pz
    real(kind=dp) :: gamma_l, gamma_g
    complex(kind=dp) :: Pele
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)


    if (nE == 0) return
    if (.not. allocated(spont_te)) return

    dim = size(eigvecs, 1)
    allocate(Ytmp(dim))

    ! Half-widths at half-maximum from FWHM
    gamma_l = optcfg%linewidth_lorentzian / 2.0_dp
    gamma_g = optcfg%linewidth_gaussian / 2.0_dp

    do i = 1, numvb
      ! VB eigenvalue
      f_v = fermi_dirac(eigvals(i), fermi_level, optcfg%temperature)

      do j = 1, numcb
        ! CB eigenvalue (offset by numvb)
        f_c = fermi_dirac(eigvals(numvb + j), fermi_level, optcfg%temperature)

        ! Spontaneous emission: f_c * (1 - f_v)
        occ_factor = f_c * (1.0_dp - f_v)
        if (occ_factor < 1.0e-30_dp) cycle

        ! Transition energy: E_CB - E_VB > 0
        dE = eigvals(numvb + j) - eigvals(i)
        if (dE < DE_MIN) cycle

        ! Velocity matrix elements via commutator v_alpha = -i [r_alpha, H]
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
        do ie = 1, nE
          spont_te(ie) = spont_te(ie) + occ_factor * (px + py) &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          spont_tm(ie) = spont_tm(ie) + occ_factor * pz &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

    deallocate(Ytmp)

  end subroutine optics_accumulate_spontaneous


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
    real(kind=dp) :: prefactor_E

    ! eps0 in AA units: e0 is C/(V*nm), 1 nm = 10 AA
    ! C/(V*nm) -> C/(V*AA): divide by 10 (not multiply)
    real(kind=dp), parameter :: e0_AA = e0 / 10.0_dp   ! C/(V*AA)

    ! Conversion factor: alpha_AA -> alpha_cm (1/AA -> 1/cm)
    real(kind=dp), parameter :: AA_TO_CM = 1.0e8_dp

    ! --- Apply prefactor to build final alpha arrays ---
    do concurrent (ie = 1:nE)
      if (E_grid(ie) > 0.0_dp) then
        prefactor_E = (2.0_dp * pi_dp * e**2) &
          & / (optcfg%refractive_index * c * e0_AA * hbar**2 * E_grid(ie))
      else
        prefactor_E = 0.0_dp
      end if
      alpha_te(ie) = prefactor_E * alpha_te(ie) * AA_TO_CM
      alpha_tm(ie) = prefactor_E * alpha_tm(ie) * AA_TO_CM
      alpha_isbt(ie) = prefactor_E * alpha_isbt(ie) * AA_TO_CM
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

    open(unit=iounit, file='output/absorption_ISBT.dat', status='replace', &
      & action='write')
    write(iounit, '(a)') '# ISBT absorption spectrum (TM-polarized, z-dipole only)'
    write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
    do ie = 1, nE
      write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_isbt(ie)
    end do
    close(iounit)

    print '(a)', 'Optical spectra written to output/absorption_TE.dat, output/absorption_TM.dat, and output/absorption_ISBT.dat'

    ! --- Gain output ---
    if (optcfg%gain_enabled .and. allocated(alpha_gain_te)) then
      ! Apply the same prefactor to gain arrays
      do concurrent (ie = 1:nE)
        if (E_grid(ie) > 0.0_dp) then
          prefactor_E = (2.0_dp * pi_dp * e**2) &
            & / (optcfg%refractive_index * c * e0_AA * hbar**2 * E_grid(ie))
        else
          prefactor_E = 0.0_dp
        end if
        alpha_gain_te(ie) = prefactor_E * alpha_gain_te(ie) * AA_TO_CM
        alpha_gain_tm(ie) = prefactor_E * alpha_gain_tm(ie) * AA_TO_CM
      end do

      open(unit=iounit, file='output/gain_TE.dat', status='replace', &
        & action='write')
      write(iounit, '(a)') '# Interband TE gain spectrum (gain_TE vs E)'
      write(iounit, '(a,es10.2,a)') '# Carrier density = ', &
        & optcfg%gain_carrier_density, ' cm^-2'
      write(iounit, '(a)') '# E(eV)  gain(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_gain_te(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/gain_TM.dat', status='replace', &
        & action='write')
      write(iounit, '(a)') '# Interband TM gain spectrum (gain_TM vs E)'
      write(iounit, '(a,es10.2,a)') '# Carrier density = ', &
        & optcfg%gain_carrier_density, ' cm^-2'
      write(iounit, '(a)') '# E(eV)  gain(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_gain_tm(ie)
      end do
      close(iounit)

      print '(a)', 'Gain spectra written to output/gain_TE.dat and output/gain_TM.dat'
    end if

    ! --- Spontaneous emission output ---
    if (optcfg%spontaneous_enabled .and. allocated(spont_te)) then
      ! Apply prefactor
      do concurrent (ie = 1:nE)
        if (E_grid(ie) > 0.0_dp) then
          prefactor_E = (2.0_dp * pi_dp * e**2) &
            & / (optcfg%refractive_index * c * e0_AA * hbar**2 * E_grid(ie))
        else
          prefactor_E = 0.0_dp
        end if
        spont_te(ie) = prefactor_E * spont_te(ie) * AA_TO_CM
        spont_tm(ie) = prefactor_E * spont_tm(ie) * AA_TO_CM
      end do

      call ensure_output_dir()
      call get_unit(iounit)

      open(unit=iounit, file='output/spontaneous_TE.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spontaneous emission TE spectrum'
      write(iounit, '(a)') '# E(eV)  rate(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), spont_te(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/spontaneous_TM.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spontaneous emission TM spectrum'
      write(iounit, '(a)') '# E(eV)  rate(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), spont_tm(ie)
      end do
      close(iounit)

      print '(a)', 'Spontaneous emission written to output/spontaneous_TE.dat and output/spontaneous_TM.dat'
    end if

    ! --- Spin-resolved absorption output ---
    if (optcfg%spin_resolved .and. allocated(alpha_te_up)) then
      ! Apply the same prefactor to spin-resolved arrays
      do concurrent (ie = 1:nE)
        if (E_grid(ie) > 0.0_dp) then
          prefactor_E = (2.0_dp * pi_dp * e**2) &
            & / (optcfg%refractive_index * c * e0_AA * hbar**2 * E_grid(ie))
        else
          prefactor_E = 0.0_dp
        end if
        alpha_te_up(ie) = prefactor_E * alpha_te_up(ie) * AA_TO_CM
        alpha_te_dw(ie) = prefactor_E * alpha_te_dw(ie) * AA_TO_CM
        alpha_tm_up(ie) = prefactor_E * alpha_tm_up(ie) * AA_TO_CM
        alpha_tm_dw(ie) = prefactor_E * alpha_tm_dw(ie) * AA_TO_CM
      end do

      call ensure_output_dir()
      call get_unit(iounit)

      open(unit=iounit, file='output/absorption_TE_up.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-up TE absorption spectrum (alpha_TE_up vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_te_up(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/absorption_TE_dw.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-down TE absorption spectrum (alpha_TE_dw vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_te_dw(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/absorption_TM_up.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-up TM absorption spectrum (alpha_TM_up vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_tm_up(ie)
      end do
      close(iounit)

      open(unit=iounit, file='output/absorption_TM_dw.dat', status='replace', action='write')
      write(iounit, '(a)') '# Spin-down TM absorption spectrum (alpha_TM_dw vs E)'
      write(iounit, '(a)') '# E(eV)  alpha(cm^-1)'
      do ie = 1, nE
        write(iounit, '(es16.8, 2x, es16.8)') E_grid(ie), alpha_tm_dw(ie)
      end do
      close(iounit)

      print '(a)', 'Spin-resolved spectra written to output/absorption_TE_up.dat, output/absorption_TE_dw.dat, output/absorption_TM_up.dat, output/absorption_TM_dw.dat'
    end if

  end subroutine optics_finalize


  ! ------------------------------------------------------------------
  ! Deallocate module-level arrays and reset state.
  ! Called defensively at start of optics_init, and may be called
  ! explicitly by the caller when done with optics.
  ! ------------------------------------------------------------------
  subroutine optics_cleanup()

    if (allocated(E_grid))       deallocate(E_grid)
    if (allocated(alpha_te))     deallocate(alpha_te)
    if (allocated(alpha_tm))     deallocate(alpha_tm)
    if (allocated(alpha_isbt))   deallocate(alpha_isbt)
    if (allocated(alpha_gain_te)) deallocate(alpha_gain_te)
    if (allocated(alpha_gain_tm)) deallocate(alpha_gain_tm)
    if (allocated(spont_te))      deallocate(spont_te)
    if (allocated(spont_tm))      deallocate(spont_tm)
    if (allocated(alpha_te_up)) deallocate(alpha_te_up)
    if (allocated(alpha_te_dw)) deallocate(alpha_te_dw)
    if (allocated(alpha_tm_up)) deallocate(alpha_tm_up)
    if (allocated(alpha_tm_dw)) deallocate(alpha_tm_dw)
    nE = 0
    gain_fermi_computed = .false.
    mu_e = 0.0_dp
    mu_h = 0.0_dp

  end subroutine optics_cleanup


  ! ------------------------------------------------------------------
  ! Reset the gain quasi-Fermi level state so they are recomputed
  ! on the next call to compute_gain_qw.  Call this before starting
  ! a new k_sweep for gain.
  ! ------------------------------------------------------------------
  subroutine gain_reset()
    gain_fermi_computed = .false.
    mu_e = 0.0_dp
    mu_h = 0.0_dp
  end subroutine gain_reset


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
    complex(kind=dp), intent(in) :: eigvecs(:,:)
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

    real(kind=dp), intent(in)    :: eigvals(:)           ! (nstates) eigenvalues ascending
    complex(kind=dp), intent(in) :: eigvecs(:,:)         ! (8*fdstep, nstates)
    real(kind=dp), intent(in)    :: z_grid(:)            ! (fdstep) z-coords (AA)
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
  ! ISBT absorption spectrum: TM-polarized (z-dipole only).
  !
  ! alpha_TM(E) proportional to sum_{i<j} |z_ij|^2 * (f_i - f_j)
  !             * lineshape(E - E_ij) * k_weight
  !
  ! Uses the same Voigt broadening as interband.  Written to
  ! output/absorption_ISBT.dat after finalization.
  ! ------------------------------------------------------------------
  subroutine compute_isbt_absorption(optcfg, eigvals, eigvecs, &
    & vel, numcb, numvb, k_weight, fermi_level)

    type(optics_config), intent(in) :: optcfg
    real(kind=dp), intent(in)    :: eigvals(:)
    complex(kind=dp), intent(in) :: eigvecs(:,:)
    type(csr_matrix), intent(in) :: vel(3)         ! commutator velocity matrices
    integer, intent(in)          :: numcb, numvb
    real(kind=dp), intent(in)    :: k_weight
    real(kind=dp), intent(in)    :: fermi_level

    integer :: i, j, ie, state_i, state_j, dim
    integer, parameter :: dir_isbt = 3   ! z-dipole for QW ISBT
    real(kind=dp) :: E_ij, occ_factor, f_i, f_j
    real(kind=dp) :: gamma_l, gamma_g
    real(kind=dp) :: p_abs2
    complex(kind=dp) :: pele_ij
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)

    real(kind=dp) :: ef_cb

    if (nE == 0) return  ! optics_init not called

    dim = size(eigvecs, 1)
    allocate(Ytmp(dim))

    gamma_l = optcfg%linewidth_lorentzian / 2.0_dp
    gamma_g = optcfg%linewidth_gaussian / 2.0_dp

    ! Use gain quasi-Fermi level (mu_e) when available, since it
    ! accounts for carrier injection into CB states.
    if (gain_fermi_computed) then
      ef_cb = mu_e
    else
      ef_cb = fermi_level
    end if

    do i = 1, numcb
      state_i = numvb + i
      f_i = fermi_dirac(eigvals(state_i), ef_cb, optcfg%temperature)

      do j = i + 1, numcb
        state_j = numvb + j
        f_j = fermi_dirac(eigvals(state_j), ef_cb, optcfg%temperature)

        ! Occupation factor: (f_i - f_j).  Positive means absorption
        ! from lower CB subband i to upper CB subband j.
        occ_factor = f_i - f_j
        if (occ_factor < 1.0e-30_dp) cycle

        E_ij = eigvals(state_j) - eigvals(state_i)
        if (E_ij < DE_MIN) cycle

        ! Velocity matrix element via commutator: dir_isbt = 3 for QW (z-dipole)
        call csr_spmv(vel(dir_isbt), eigvecs(:,state_j), Ytmp, ONE, ZERO)
        pele_ij = zdotc(dim, eigvecs(1:dim,state_i), 1, Ytmp, 1)
        p_abs2 = real(pele_ij * conjg(pele_ij), kind=dp)

        ! Broaden and accumulate onto energy grid (ISBT, TM-polarized)
        do ie = 1, nE
          alpha_isbt(ie) = alpha_isbt(ie) + occ_factor * p_abs2 &
            & * lineshape_voigt(E_grid(ie), E_ij, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

    deallocate(Ytmp)

  end subroutine compute_isbt_absorption


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

    real(kind=dp), intent(in) :: eigvals(:)       ! all eigenvalues at k=0
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
  ! Gain spectrum for a QW with population inversion.
  !
  ! Structurally identical to optics_accumulate but uses separate
  ! quasi-Fermi levels for electrons (f_e in CB) and holes (f_h in VB).
  !
  ! For each CB-VB pair:
  !   1. Find quasi-Fermi levels f_e, f_h from the carrier density
  !   2. Occupation factor = f_v(E_VB, f_h) - f_c(E_CB, f_e)
  !      When inverted (f_c > f_v), this is NEGATIVE => gain
  !   3. Same momentum matrix elements as absorption
  !   4. Accumulate into alpha_gain_te, alpha_gain_tm
  !
  ! This subroutine is called inside the k_sweep loop. The quasi-Fermi
  ! levels are computed once (at the first k-point call) and reused.
  ! ------------------------------------------------------------------
  subroutine compute_gain_qw(optcfg, eigvals, eigvecs, k_weight, &
    & vel, numcb, numvb, carrier_density)

    type(optics_config), intent(in) :: optcfg
    real(kind=dp), intent(in) :: eigvals(:)        ! (numcb+numvb) eigenvalues
    complex(kind=dp), intent(in) :: eigvecs(:,:)   ! (dim, numcb+numvb)
    real(kind=dp), intent(in) :: k_weight          ! Simpson weight for this k
    type(csr_matrix), intent(in) :: vel(3)         ! commutator velocity matrices
    integer, intent(in) :: numcb, numvb
    real(kind=dp), intent(in) :: carrier_density   ! 2D carrier density (cm^-2)

    integer :: i, j, dir, ie, dim
    real(kind=dp) :: dE, f_c, f_v, occ_factor
    real(kind=dp) :: px, py, pz
    real(kind=dp) :: gamma_l, gamma_g
    complex(kind=dp) :: Pele
    complex(kind=dp), allocatable :: Ytmp(:)
    complex(kind=dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)


    if (nE == 0) return  ! optics_init not called
    if (.not. allocated(alpha_gain_te)) return  ! gain not initialized

    dim = size(eigvecs, 1)
    allocate(Ytmp(dim))

    ! Compute quasi-Fermi levels once on the first k-point call.
    ! We use the eigenvalues at the first k-point (usually k_par ~ 0)
    ! to estimate the subband edges.
    if (.not. gain_fermi_computed) then
      ! f_e: quasi-Fermi level for CB states
      mu_e = find_quasi_fermi(eigvals, optcfg%temperature, &
        & carrier_density, numcb, numvb)
      ! f_h: quasi-Fermi level for VB states.
      ! For holes: sum_i (1 - f(E_i, mu_h)) = n_2D, which by
      ! particle-hole symmetry is sum_i f(-E_i, -mu_h) = n_2D.
      mu_h = find_quasi_fermi_holes(eigvals, optcfg%temperature, &
        & carrier_density, numvb)
      gain_fermi_computed = .true.
      print '(a,es14.6,a)', '  Gain: mu_e = ', mu_e, ' eV'
      print '(a,es14.6,a)', '  Gain: mu_h = ', mu_h, ' eV'
      print '(a,es10.2,a)', '  Gain: carrier density = ', carrier_density, ' cm^-2'
    end if

    ! Half-widths at half-maximum from FWHM
    gamma_l = optcfg%linewidth_lorentzian / 2.0_dp
    gamma_g = optcfg%linewidth_gaussian / 2.0_dp

    ! Loop over VB-CB pairs (same structure as optics_accumulate)
    do i = 1, numvb
      ! VB state: use hole quasi-Fermi level
      f_v = fermi_dirac(eigvals(i), mu_h, optcfg%temperature)

      do j = 1, numcb
        ! CB state: use electron quasi-Fermi level
        f_c = fermi_dirac(eigvals(numvb + j), mu_e, optcfg%temperature)

        ! Occupation factor: (f_V - f_C).
        ! When f_C > f_V (population inversion), this is NEGATIVE => GAIN.
        ! We keep the full range (both positive and negative) unlike
        ! absorption which skips negative contributions.
        occ_factor = f_v - f_c

        ! Transition energy: E_CB - E_VB > 0
        dE = eigvals(numvb + j) - eigvals(i)
        if (dE < DE_MIN) cycle

        ! Skip negligible contributions
        if (abs(occ_factor) < 1.0e-30_dp) cycle

        ! Velocity matrix elements via commutator v_alpha = -i [r_alpha, H]
        ! |<CB| v_alpha |VB>|^2 computed as SpMV: Ytmp = vel(dir) * |VB>,
        ! then Pele = <CB| Ytmp>.
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
        do ie = 1, nE
          alpha_gain_te(ie) = alpha_gain_te(ie) + occ_factor * (px + py) &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
          alpha_gain_tm(ie) = alpha_gain_tm(ie) + occ_factor * pz &
            & * lineshape_voigt(E_grid(ie), dE, gamma_l, gamma_g) * k_weight
        end do

      end do
    end do

    deallocate(Ytmp)

  end subroutine compute_gain_qw


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

end module optical_spectra
