program kpfdm

  use definitions
  use parameters
  use hamiltonianConstructor
  use confinement_init, only: confinementInitialization_landau
  use OMP_lib
  use outputFunctions
  use input_parser
  use sc_loop
  use eigensolver, only: eigensolver_result, eigensolver_result_free
  use strain_solver, only: strain_result, strain_result_free
  use exciton_solver
  use scattering_solver
  use linalg, only: zheevx, mkl_set_num_threads_local, ilaenv, dlamch
  use spin_projection, only: compute_band_parts
  use simulation_setup_mod

  implicit none

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg
  real(kind=dp), allocatable, dimension(:,:) :: profile
  real(kind=dp), allocatable, dimension(:,:,:) :: kpterms

  ! wave vector
  type(wavevector), allocatable, dimension(:) :: smallk

  ! iteration consts
  integer :: i, k, ii, jj

  ! hamiltonian and LAPACK/BLAS
  integer :: info, NB, lwork, N, M, il, iuu
  real(kind=dp) :: abstol, vl, vu
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  complex(kind=dp), allocatable, dimension(:,:,:) :: eigv
  complex(kind=dp), allocatable, dimension(:,:) :: HT, HTmp
  integer, allocatable :: iwork(:), ifail(:)

  ! file handling
  integer(kind=4) :: iounit

  ! --- QW SC charge density output ---
  real(kind=dp), allocatable       :: sc_ne_qw(:), sc_nh_qw(:)
  real(kind=dp)                    :: sc_fermi_level
  real(kind=dp)                    :: sc_delta_phi
  logical                          :: sc_converged_flag
  integer                          :: sc_iterations

  ! --- Wire mode variables (local to wire block below) ---

  ! Shared setup: read input, initialize materials, confinement, external field
  call read_config(cfg)

  ! Semantic validation (no additional constraints for bandStructure)
  call validate_semantic(cfg, 'bandStructure')

  ! Bulk/QW initialization via simulation_setup (confinement, strain).
  ! SC is deferred until after wave vector setup — bandStructure needs the
  ! full k-point array for SC charge integration, while simulation_setup_init
  ! only uses a single k-point (designed for gfactor/optics executables).
  if (cfg%confinement == 'bulk' .or. cfg%confinement == 'qw') then
    block
      type(simulation_setup) :: setup
      call simulation_setup_init(cfg, setup, skip_sc=.true.)
      if (cfg%confinement == 'qw') then
        call move_alloc(setup%profile, profile)
        call move_alloc(setup%kpterms, kpterms)
      end if
      if (setup%has_strain) print *, '  QW strain calculation complete'
    end block
  end if

  ! Build wave vector array
  allocate(smallk(cfg%wave_vector%nsteps))
  smallk%kx = 0
  smallk%ky = 0
  smallk%kz = 0
  select case(cfg%wave_vector%mode)

    case ("kx")
      smallk%kx = [ ((i-1)*cfg%wave_vector%max/(cfg%wave_vector%nsteps-1), i=1,cfg%wave_vector%nsteps) ]

    case ("ky")
      smallk%ky = [ ((i-1)*cfg%wave_vector%max/(cfg%wave_vector%nsteps-1), i=1,cfg%wave_vector%nsteps) ]

    case ("kz")
      smallk%kz = [ ((i-1)*cfg%wave_vector%max/(cfg%wave_vector%nsteps-1), i=1,cfg%wave_vector%nsteps) ]

    case ("kxky")
      ! [110] direction: kx = ky = k, kz = 0
      do i = 1, cfg%wave_vector%nsteps
        smallk(i)%kx = (i-1) * cfg%wave_vector%max / (cfg%wave_vector%nsteps - 1)
        smallk(i)%ky = smallk(i)%kx
      end do

    case ("kxkz")
      ! [101] direction: kx = kz = k, ky = 0
      do i = 1, cfg%wave_vector%nsteps
        smallk(i)%kx = (i-1) * cfg%wave_vector%max / (cfg%wave_vector%nsteps - 1)
        smallk(i)%kz = smallk(i)%kx
      end do

    case ("kykz")
      ! [011] direction: ky = kz = k, kx = 0
      do i = 1, cfg%wave_vector%nsteps
        smallk(i)%ky = (i-1) * cfg%wave_vector%max / (cfg%wave_vector%nsteps - 1)
        smallk(i)%kz = smallk(i)%ky
      end do

    case ("k0")
      ! Already zeroed above

    case default
      stop "no such direction"

  end select

  ! ====================================================================
  ! Wire mode (confinement=2): sparse path via simulation_setup
  ! ====================================================================
  if (cfg%confinement == 'wire') then
    block
      type(simulation_setup) :: setup
      type(eigensolver_result) :: eigen_res
      real(kind=dp), allocatable :: eig_wire(:,:)
      real(kind=dp), allocatable :: prev_wire_eval(:)
      complex(kind=dp), allocatable :: prev_wire_evec(:,:)
      integer :: Ngrid, Ntot, nev_wire, max_nev_found
      ! Diagnostic outputs from simulation_setup_init
      type(strain_result), allocatable :: strain_diag
      real(kind=dp), allocatable :: sc_phi_diag(:,:), sc_ne_diag(:,:), sc_nh_diag(:,:)

      ! Delegate confinement, strain, SC, eigensolver to simulation_setup_init.
      ! Optional arguments capture strain/SC diagnostics for file output below.
      call simulation_setup_init(cfg, setup, strain_out=strain_diag, &
        sc_phi_out=sc_phi_diag, sc_ne_out=sc_ne_diag, sc_nh_out=sc_nh_diag)

      ! --- Write strain tensor (if strain was computed) ---
      if (allocated(strain_diag)) then
        call ensure_output_dir()
        call get_unit(iounit)
        open(unit=iounit, file='output/strain.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) eps_xx eps_yy eps_zz eps_xy eps_xz eps_yz'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            k = (jj - 1) * cfg%grid%nx + ii
            write(unit=iounit, fmt='(8(g14.6,1x))') &
              cfg%grid%x(ii), cfg%grid%z(jj), &
              strain_diag%eps_xx(k), strain_diag%eps_yy(k), &
              strain_diag%eps_zz(k), 0.0_dp, 0.0_dp, &
              strain_diag%eps_yz(k)
          end do
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  Wire strain tensor written to output/strain.dat'
        call strain_result_free(strain_diag)
      end if

      ! --- Write SC diagnostics (if SC was run) ---
      if (allocated(sc_phi_diag)) then
        call ensure_output_dir()
        call get_unit(iounit)
        open(unit=iounit, file='output/sc_potential_profile.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) EV EV_DeltaSO EC'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            k = (jj - 1) * cfg%grid%nx + ii
            write(unit=iounit, fmt='(5(g14.6,1x))') &
              cfg%grid%x(ii), cfg%grid%z(jj), &
              setup%profile_2d(k, 1), setup%profile_2d(k, 2), setup%profile_2d(k, 3)
          end do
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  SC band edge profile written to output/sc_potential_profile.dat'

        call get_unit(iounit)
        open(unit=iounit, file='output/sc_phi.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) phi'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            write(unit=iounit, fmt='(3(g14.6,1x))') &
              cfg%grid%x(ii), cfg%grid%z(jj), sc_phi_diag(ii, jj)
          end do
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  SC electrostatic potential written to output/sc_phi.dat'

        call get_unit(iounit)
        open(unit=iounit, file='output/sc_charge.dat', status='replace', action='write')
        write(iounit, '(A)') '# x(A) y(A) n_e n_h'
        do jj = 1, cfg%grid%ny
          do ii = 1, cfg%grid%nx
            write(unit=iounit, fmt='(4(g14.6,1x))') &
              cfg%grid%x(ii), cfg%grid%z(jj), sc_ne_diag(ii, jj), sc_nh_diag(ii, jj)
          end do
          write(iounit, '(A)') ''
        end do
        close(iounit)
        print *, '  SC charge density written to output/sc_charge.dat'

        deallocate(sc_phi_diag, sc_ne_diag, sc_nh_diag)
      end if

      Ngrid = setup%Ngrid
      Ntot = setup%Ntot
      nev_wire = setup%nev_wire

      print *, ''
      print *, '=== Wire mode (2D confinement) ==='
      print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', Ngrid
      print *, '  Matrix size: ', Ntot, 'x', Ntot
      print *, '  Requesting ', nev_wire, ' eigenvalues per k-point'
      print *, '  kz sweep: ', cfg%wave_vector%nsteps, ' points, kz_max=', cfg%wave_vector%max
      print *, ''

      if (setup%has_strain) print *, '  Wire strain calculation complete'
      if (setup%sc_was_run) then
        print *, '  Wire SC loop complete. profile_2d updated.'
      end if

      ! Allocate storage for eigenvalues across k-points
      allocate(eig_wire(Ntot, cfg%wave_vector%nsteps))
      eig_wire = 0.0_dp
      max_nev_found = 0

      ! --- Write 2D band edge profile (after strain + SC, before k-sweep) ---
      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/potential_profile.dat', status='replace', action='write')
      if (allocated(cfg%strain_blocks%delta_Ec)) then
        write(iounit, '(A)') '# x(A) y(A) EV_top_strained ESO_strained EC_strained'
      else
        write(iounit, '(A)') '# x(A) y(A) EV EV_DeltaSO EC'
      end if
      do jj = 1, cfg%grid%ny
        do ii = 1, cfg%grid%nx
          k = (jj - 1) * cfg%grid%nx + ii
          if (allocated(cfg%strain_blocks%delta_Ec)) then
            write(unit=iounit, fmt='(5(g14.6,1x))') &
              & cfg%grid%x(ii), cfg%grid%z(jj), &
              & max(setup%profile_2d(k, 1) + cfg%strain_blocks%delta_EHH(k), &
                    setup%profile_2d(k, 1) + cfg%strain_blocks%delta_ELH(k)), &
              & setup%profile_2d(k, 2) + cfg%strain_blocks%delta_ESO(k), &
              & setup%profile_2d(k, 3) + cfg%strain_blocks%delta_Ec(k)
          else
            write(unit=iounit, fmt='(5(g14.6,1x))') &
              & cfg%grid%x(ii), cfg%grid%z(jj), &
              & setup%profile_2d(k, 1), setup%profile_2d(k, 2), setup%profile_2d(k, 3)
          end if
        end do
        write(iounit, '(A)') ''
      end do
      close(iounit)
      print *, '  Wire band edge profile written to output/potential_profile.dat'

      ! ================================================================
      ! Wire kz sweep: serial k=1 (build CSR) + serial k=2..N
      ! ================================================================

      ! --- Serial k=1: build CSR, solve ---
      print *, 'k-point 1/', cfg%wave_vector%nsteps, ' kz=', smallk(1)%kz

      call setup_build_H(setup, cfg, smallk(1))

      if (cfg%feast%emin /= 0.0_dp .or. cfg%feast%emax /= 0.0_dp) then
        print *, '  Manual energy window: [', setup%eigen_cfg%emin, ',', setup%eigen_cfg%emax, ']'
      else
        print *, '  Auto energy window: [', setup%eigen_cfg%emin, ',', setup%eigen_cfg%emax, ']'
      end if

      call setup%eigen_solver%solve(setup%HT_csr_ptr, setup%eigen_cfg, eigen_res)

      if (.not. eigen_res%converged) then
        if (eigen_res%nev_found < nev_wire) then
          print *, '  ERROR: FEAST did not converge at k-point 1 and found only', &
            eigen_res%nev_found, 'eigenvalues (need', nev_wire, ')'
          stop 1
        end if
        print *, '  WARNING: FEAST subspace issue at k-point 1, but found enough eigenvalues'
      end if
      print *, '  Found ', eigen_res%nev_found, ' eigenvalues (requested ', &
        & nev_wire, ') iterations=', eigen_res%iterations

      if (eigen_res%nev_found > 0) then
        do i = 1, min(eigen_res%nev_found, Ntot)
          eig_wire(i, 1) = eigen_res%eigenvalues(i)
        end do
        max_nev_found = max(max_nev_found, eigen_res%nev_found)
      end if

      ! Write 2D eigenfunctions for k=1 (with band decomposition)
      if (eigen_res%nev_found > 0) then
        call writeEigenfunctions2d(cfg%grid, eigen_res%eigenvalues, &
          & eigen_res%eigenvectors, 1, eigen_res%nev_found, write_parts=.true.)
        print *, '  Wrote ', eigen_res%nev_found, ' 2D wavefunctions to output/eigenfunctions_k_00001_ev_*.dat'
      end if

      if (eigen_res%nev_found > 0) then
        allocate(prev_wire_eval(eigen_res%nev_found))
        allocate(prev_wire_evec(Ntot, eigen_res%nev_found))
        prev_wire_eval = eigen_res%eigenvalues
        prev_wire_evec = eigen_res%eigenvectors
      end if

      call eigensolver_result_free(eigen_res)

      ! --- Serial k=2..N with branch tracking against the previous k-point ---
      if (cfg%wave_vector%nsteps > 1) then
        info = mkl_set_num_threads_local(1)

        print '(A,I0,A)', ' Wire kz-sweep: k=2..', cfg%wave_vector%nsteps, &
          & ' (serial branch tracking)'

        do k = 2, cfg%wave_vector%nsteps
          print *, 'k-point ', k, '/', cfg%wave_vector%nsteps, ' kz=', smallk(k)%kz

          call setup_build_H(setup, cfg, smallk(k))
          call setup%eigen_solver%solve(setup%HT_csr_ptr, setup%eigen_cfg, eigen_res)

          if (.not. eigen_res%converged) then
            if (eigen_res%nev_found < nev_wire) then
              print *, '  ERROR: FEAST did not converge at k-point', k, &
                'and found only', eigen_res%nev_found, 'eigenvalues (need', nev_wire, ')'
              stop 1
            end if
            print *, '  WARNING: eigensolver subspace issue at k-point', k, &
              ' but found enough eigenvalues'
          end if

          if (eigen_res%nev_found > 0) then
            call reorder_wire_branches(prev_wire_eval, prev_wire_evec, &
              eigen_res%eigenvalues, eigen_res%eigenvectors)
            do i = 1, min(eigen_res%nev_found, Ntot)
              eig_wire(i, k) = eigen_res%eigenvalues(i)
            end do
            max_nev_found = max(max_nev_found, eigen_res%nev_found)

            prev_wire_eval = eigen_res%eigenvalues
            prev_wire_evec = eigen_res%eigenvectors

            if (k == cfg%wave_vector%nsteps/2 .or. k == cfg%wave_vector%nsteps) then
              call writeEigenfunctions2d(cfg%grid, eigen_res%eigenvalues, &
                & eigen_res%eigenvectors, k, eigen_res%nev_found, write_parts=.false.)
              print *, '  Wrote ', eigen_res%nev_found, &
                & ' 2D wavefunctions to output/eigenfunctions_k_', k, '_ev_*.dat'
            end if
          end if

          call eigensolver_result_free(eigen_res)
        end do
      end if

      ! Write eigenvalues to file
      if (max_nev_found > 0) then
        call writeEigenvalues(smallk, eig_wire(1:max_nev_found, :), cfg%wave_vector%nsteps, cfg)
      else
        call writeEigenvalues(smallk, eig_wire(1:nev_wire, :), cfg%wave_vector%nsteps, cfg)
      end if
      print *, ''
      print *, 'Wire band structure written to output/eigenvalues.dat'

      ! Clean up local wire allocations
      if (allocated(eig_wire)) deallocate(eig_wire)
      if (allocated(prev_wire_eval)) deallocate(prev_wire_eval)
      if (allocated(prev_wire_evec)) deallocate(prev_wire_evec)

      ! Free setup (profile_2d, kpterms_2d, eigen_solver, HT_csr_ptr, etc.)
      call simulation_setup_free(setup)

    end block

    ! Clean up common allocations and exit
    if (allocated(smallk))  deallocate(smallk)
    if (allocated(profile)) deallocate(profile)
    if (allocated(kpterms)) deallocate(kpterms)

    stop  ! wire mode complete

  end if

  ! ====================================================================
  ! Landau mode initialization (confinement=3)
  ! ====================================================================
  if (cfg%confinement == 'landau') then
    block
      character(len=255) :: lmat(1)
      type(paramStruct) :: lparams(1)

      lmat(1) = cfg%material_names(1)
      lparams(1) = cfg%params(1)

      N = cfg%landau%nx * 8
      allocate(kpterms(cfg%landau%nx, cfg%landau%nx, 10))
      kpterms = 0.0_dp
      call confinementInitialization_landau(cfg%grid%x, lmat, lparams, &
        & profile, kpterms, cfg%FDorder)
    end block

    print '(A,I0,A)', ' Landau mode: N=', N, ' (single material)'
  end if

  ! ====================================================================
  ! Bulk / QW mode (confinement=0 or 1): dense LAPACK path
  ! ====================================================================

  ! Set matrix dimensions and eigenvalue range
  N = cfg%ngrid * 8
  if (conf_direction(cfg%confinement) == 'n') then
    N = 8  ! For bulk, we only need 8x8
    if (cfg%evnum > 8) then
      print *, "Warning: requesting more bands than available in bulk (8). Using all 8 bands."
      cfg%evnum = 8
      cfg%bands%num_cb = 2
      cfg%bands%num_vb = 6
    end if
  end if

  ! For quantum well, check if requested bands are within limits
  if (conf_direction(cfg%confinement) == 'z') then
    if (cfg%bands%num_cb > NUM_CB_STATES*cfg%ngrid) then
      print *, "Warning: requesting more conduction bands than available. Limiting to ", NUM_CB_STATES*cfg%ngrid
      cfg%bands%num_cb = NUM_CB_STATES*cfg%ngrid
    end if
    if (cfg%bands%num_vb > NUM_VB_STATES*cfg%ngrid) then
      print *, "Warning: requesting more valence bands than available. Limiting to ", NUM_VB_STATES*cfg%ngrid
      cfg%bands%num_vb = NUM_VB_STATES*cfg%ngrid
    end if
    cfg%evnum = cfg%bands%num_cb + cfg%bands%num_vb
  end if

  vl = 0.0_dp
  vu = 0.0_dp
  if (cfg%confinement == 'bulk') then
    il = 1
    iuu = cfg%evnum
  else
    ! For quantum well, select the right range of states
    ! We want the highest numvb valence states and lowest numcb conduction states
    il = NUM_VB_STATES*cfg%ngrid - cfg%bands%num_vb + 1  ! Start from highest valence band
    iuu = NUM_VB_STATES*cfg%ngrid + cfg%bands%num_cb     ! Up to highest conduction band
    print *, "Computing states from index", il, "to", iuu
  end if

  NB = ILAENV(1, 'ZHETRD', 'UPLO', N, N, -1, -1)
  NB = MAX(NB,N)
  ABSTOL = DLAMCH('P')

  if (allocated(rwork)) deallocate(rwork)
  allocate(rwork(7*N))  ! For ZHEEVX
  if (allocated(iwork)) deallocate(iwork)
  allocate(iwork(5*N))
  if (allocated(ifail)) deallocate(ifail)
  allocate(ifail(N))
  if (allocated(eig)) deallocate(eig)
  allocate(eig(iuu-il+1,cfg%wave_vector%nsteps))  ! Only store the states we want
  if (allocated(eigv)) deallocate(eigv)
  if (conf_direction(cfg%confinement) == 'n') then
    allocate(eigv(8,cfg%evnum,cfg%wave_vector%nsteps))  ! 8x8 for bulk
  else
    allocate(eigv(N,iuu-il+1,cfg%wave_vector%nsteps))  ! Only store the states we want
  end if

  eig(:,:) = 0_dp
  eigv(:,:,:) = 0_dp

  allocate(HT(N,N))
  allocate(HTmp(8,8))  ! Temporary array for bulk diagonalization
  HT = 0.0_dp
  HTmp = 0.0_dp

  ! Initial workspace allocation
  if (allocated(work)) deallocate(work)
  allocate(work(N))  ! Will be resized after workspace query

  ! Print profile for QW and Landau modes
  if (conf_direction(cfg%confinement) == 'z') then
    call ensure_output_dir()
    call get_unit(iounit)
    open(unit=iounit, file='output/potential_profile.dat', status="replace", action="write")
    do i = 1, cfg%ngrid, 1
      write(iounit,*) cfg%z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)
  else if (conf_direction(cfg%confinement) == 'x') then
    ! Landau mode: print profile along x
    call ensure_output_dir()
    call get_unit(iounit)
    open(unit=iounit, file='output/potential_profile.dat', status="replace", action="write")
    do i = 1, cfg%landau%nx, 1
      write(iounit,*) cfg%grid%x(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)
  end if

  ! --- Self-consistent loop (QW only, after wave vector setup) ---
  if (conf_direction(cfg%confinement) == 'z' .and. cfg%sc%enabled == 1) then
    print *, ''
    print *, '=== Running self-consistent Schrodinger-Poisson loop ==='
    call self_consistent_loop(profile, cfg, kpterms, HT, eig, eigv, &
      & smallk, N, il, iuu, n_electron_out=sc_ne_qw, n_hole_out=sc_nh_qw, &
      & fermi_level_out=sc_fermi_level, converged_out=sc_converged_flag, &
      & iterations_out=sc_iterations, delta_phi_out=sc_delta_phi)

    ! Write updated profile after SC convergence
    call get_unit(iounit)
    open(unit=iounit, file='output/sc_potential_profile.dat', status="replace", action="write")
    do i = 1, cfg%ngrid, 1
      write(iounit,*) cfg%z(i), profile(i,1), profile(i,2), profile(i,3)
    end do
    close(iounit)
    print *, 'SC potential profile written to output/sc_potential_profile.dat'

    ! Write SC summary (converged flag, iterations, Fermi level)
    call ensure_output_dir()
    call get_unit(iounit)
    open(unit=iounit, file='output/sc_summary.dat', status="replace", action="write")
    write(iounit, '(A)') '# converged  iterations  |dPhi|  fermi_level(eV)'
    write(iounit, '(L1,1x,I6,1x,ES14.6,1x,ES14.6)') sc_converged_flag, sc_iterations, sc_delta_phi, sc_fermi_level
    close(iounit)
    print *, 'SC summary written to output/sc_summary.dat'

    ! Write charge density
    if (allocated(sc_ne_qw) .and. allocated(sc_nh_qw)) then
      call ensure_output_dir()
      call get_unit(iounit)
      open(unit=iounit, file='output/sc_charge.dat', status="replace", action="write")
      write(iounit, '(A)') '# z(A) n_e(cm^-3) n_h(cm^-3)'
      do i = 1, cfg%ngrid, 1
        write(iounit, '(3(g14.6,1x))') cfg%z(i), sc_ne_qw(i), sc_nh_qw(i)
      end do
      close(iounit)
      print *, 'SC charge density written to output/sc_charge.dat'
      deallocate(sc_ne_qw, sc_nh_qw)
    end if
  end if

  ! ====================================================================
  ! Workspace query (serial, k=1) — determines lwork for zheevx
  ! ====================================================================
  if (conf_direction(cfg%confinement) == 'n') then
    ! BULK: 8x8 workspace query
    call ZB8bandBulk(HT, smallk(1), cfg%params(1), cfg=cfg)
    HTmp = HT(1:8,1:8)
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:cfg%evnum,1), &
               HTmp, 8, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print *, 'Error: zheevx bulk workspace query failed, info =', info
      stop 1
    end if
    lwork = int(real(work(1)))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
  else if (conf_direction(cfg%confinement) == 'z') then
    ! QW: NxN workspace query
    call ZB8bandQW(HT, smallk(1), profile, kpterms, cfg=cfg)
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, eig(:,1), &
               HT, N, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print *, 'Error: zheevx QW workspace query failed, info =', info
      stop 1
    end if
    lwork = int(real(work(1)))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
  else if (conf_direction(cfg%confinement) == 'x') then
    ! LANDAU: NxN workspace query
    call ZB8bandLandau(HT, smallk(1), profile, kpterms, cfg%grid%x, cfg=cfg)
    if (allocated(work)) deallocate(work)
    allocate(work(1))
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, vl, vu, il, iuu, abstol, M, eig(:,1), &
               HT, N, work, lwork, rwork, iwork, ifail, info)
    if (info /= 0) then
      print *, 'Error: zheevx Landau workspace query failed, info =', info
      stop 1
    end if
    lwork = int(real(work(1)))
    if (allocated(work)) deallocate(work)
    allocate(work(lwork))
  end if

  ! ====================================================================
  ! k-vector sweep: sequential for bulk, OpenMP parallel for QW/Landau
  ! ====================================================================
  if (conf_direction(cfg%confinement) == 'n') then
    ! --- BULK (8x8, trivially fast — no parallelization needed) ---
    do k = 1, cfg%wave_vector%nsteps
      call ZB8bandBulk(HT, smallk(k), cfg%params(1), cfg=cfg)
      HTmp = HT(1:8,1:8)

      call zheevx('V', 'I', 'U', 8, HTmp, 8, vl, vu, il, iuu, abstol, M, eig(1:cfg%evnum,k), &
                 HTmp, 8, work, lwork, rwork, iwork, ifail, info)
      if (info /= 0) then
        print *, "Diagonalization error in bulk calculation, info = ", info
        if (info < 0) print *, "Parameter ", -info, " had illegal value"
        stop "error diag"
      end if

      eigv(:,:,k) = HTmp(:,1:cfg%evnum)
    end do

    ! Write bulk eigenfunctions at start, middle, end k-points
    do k = 1, cfg%wave_vector%nsteps
      if (k == 1 .or. k == int(cfg%wave_vector%nsteps/2) .or. k == cfg%wave_vector%nsteps) then
        call writeEigenfunctions(8, min(cfg%evnum,8), eigv(:,1:min(cfg%evnum,8),k), &
          & k, cfg%ngrid, cfg%z, .true., &
          & k_magnitude=sqrt(smallk(k)%kx**2 + smallk(k)%ky**2 + smallk(k)%kz**2))
      end if
    end do

  else if (conf_direction(cfg%confinement) == 'z') then
    ! --- QUANTUM WELL (NxN, OpenMP parallel) ---
    ! Disable MKL internal threading so each OpenMP thread calls
    ! zheevx in serial — avoids oversubscription with intel_thread MKL.
    ! mkl_set_num_threads_local returns previous setting (discarded).
    info = mkl_set_num_threads_local(1)

    print '(A,I0,A)', ' QW k-sweep: ', cfg%wave_vector%nsteps, ' k-points (OpenMP parallel)'

    block
      ! Thread-private temporaries for the parallel region
      complex(kind=dp), allocatable :: HT_loc(:,:), work_loc(:)
      real(kind=dp), allocatable    :: rwork_loc(:)
      integer, allocatable          :: iwork_loc(:), ifail_loc(:)
      integer :: info_loc, M_loc

      !$omp parallel private(k, HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc, info_loc, M_loc)
      ! Each thread allocates its own workspace
      allocate(HT_loc(N, N))
      allocate(work_loc(lwork))
      allocate(rwork_loc(7*N))
      allocate(iwork_loc(5*N))
      allocate(ifail_loc(N))
      HT_loc = (0.0_dp, 0.0_dp)

      !$omp do schedule(static)
      do k = 1, cfg%wave_vector%nsteps
        ! Build Hamiltonian for this k-point
        call ZB8bandQW(HT_loc, smallk(k), profile, kpterms, cfg=cfg)

        ! Diagonalize
        call zheevx('V', 'I', 'U', N, HT_loc, N, vl, vu, il, iuu, abstol, M_loc, &
                   eig(:,k), HT_loc, N, work_loc, lwork, rwork_loc, iwork_loc, &
                   ifail_loc, info_loc)
        if (info_loc /= 0) then
          !$omp critical
          print *, "ERROR: diagonalization failed at k=", k, "info=", info_loc
          if (info_loc < 0) print *, "  Parameter ", -info_loc, " had illegal value"
          !$omp end critical
          stop 1
        end if

        ! Store eigenvectors (HT_loc now holds them, zheevx overwrites input)
        eigv(:,:,k) = HT_loc(:, 1:iuu-il+1)
      end do
      !$omp end do

      deallocate(HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc)
      !$omp end parallel
    end block

    ! Write QW eigenfunctions at start, middle, end k-points (serial)
    do k = 1, cfg%wave_vector%nsteps
      if (k == 1 .or. k == int(cfg%wave_vector%nsteps/2) .or. k == cfg%wave_vector%nsteps) then
        call writeEigenfunctions(N, iuu-il+1, eigv(:,1:iuu-il+1,k), &
          & k, cfg%ngrid, cfg%z, .false.)
      end if
    end do

  else if (conf_direction(cfg%confinement) == 'x') then
    ! --- LANDAU (8NxN, OpenMP parallel) ---
    ! Disable MKL internal threading so each OpenMP thread calls
    ! zheevx in serial — avoids oversubscription.
    info = mkl_set_num_threads_local(1)

    print '(A,I0,A,I0,A)', ' Landau: N=', cfg%landau%nx, ', ', &
      & cfg%wave_vector%nsteps, ' k-points (OpenMP parallel)'

    block
      ! Thread-private temporaries for the parallel region
      complex(kind=dp), allocatable :: HT_loc(:,:), work_loc(:)
      real(kind=dp), allocatable    :: rwork_loc(:)
      integer, allocatable          :: iwork_loc(:), ifail_loc(:)
      integer :: info_loc, M_loc

      !$omp parallel private(k, HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc, info_loc, M_loc)
      ! Each thread allocates its own workspace
      allocate(HT_loc(N, N))
      allocate(work_loc(lwork))
      allocate(rwork_loc(7*N))
      allocate(iwork_loc(5*N))
      allocate(ifail_loc(N))
      HT_loc = (0.0_dp, 0.0_dp)

      !$omp do schedule(static)
      do k = 1, cfg%wave_vector%nsteps
        ! Build Landau Hamiltonian for this k-point
        call ZB8bandLandau(HT_loc, smallk(k), profile, kpterms, cfg%grid%x, cfg=cfg)

        ! Diagonalize
        call zheevx('V', 'I', 'U', N, HT_loc, N, vl, vu, il, iuu, abstol, M_loc, &
                   eig(:,k), HT_loc, N, work_loc, lwork, rwork_loc, iwork_loc, &
                   ifail_loc, info_loc)
        if (info_loc /= 0) then
          !$omp critical
          print *, "ERROR: Landau diagonalization failed at k=", k, "info=", info_loc
          if (info_loc < 0) print *, "  Parameter ", -info_loc, " had illegal value"
          !$omp end critical
          stop 1
        end if

        ! Store eigenvectors (HT_loc now holds them, zheevx overwrites input)
        eigv(:,:,k) = HT_loc(:, 1:iuu-il+1)
      end do
      !$omp end do

      deallocate(HT_loc, work_loc, rwork_loc, iwork_loc, ifail_loc)
    !$omp end parallel
    end block

    ! Write Landau eigenfunctions at start, middle, end k-points (serial)
    do k = 1, cfg%wave_vector%nsteps
      if (k == 1 .or. k == int(cfg%wave_vector%nsteps/2) .or. k == cfg%wave_vector%nsteps) then
        call writeEigenfunctions(N, iuu-il+1, eigv(:,1:iuu-il+1,k), &
          & k, cfg%landau%nx, cfg%grid%x, .false.)
      end if
    end do

    ! ==================================================================
    ! B-sweep fan diagram for Landau mode (landau_sweep='B')
    ! ==================================================================
    if (trim(cfg%landau%sweep) == 'B') then
      block
        real(kind=dp) :: B_min, B_max, B_step, B_val
        integer :: nB, iB, nL
        real(kind=dp), allocatable :: eig_B(:,:)
        complex(kind=dp), allocatable :: HT_B(:,:), work_B(:)
        real(kind=dp), allocatable :: rwork_B(:)
        integer, allocatable :: iwork_B(:), ifail_B(:)
        integer :: info_B, M_B
        character(len=64) :: fmt_str
        type(wavevector) :: wv0

        B_min  = cfg%bdg%B_sweep(1)
        B_max  = cfg%bdg%B_sweep(2)
        B_step = cfg%bdg%B_sweep(3)

        if (B_step <= 0.0_dp) then
          print *, 'ERROR: b_sweep step must be > 0, got:', B_step
          stop 1
        end if
        nB = int((B_max - B_min) / B_step) + 1
        nL = iuu - il + 1

        print '(A,I0,A,F6.2,A,F6.2,A,F6.2)', ' Landau B-sweep: ', nB, &
          & ' points, B=[', B_min, ',', B_max, '] T, step=', B_step

        allocate(eig_B(nL, nB))
        eig_B = 0.0_dp

        ! Fixed Gamma-point wavevector
        wv0%kx = 0.0_dp
        wv0%ky = 0.0_dp
        wv0%kz = 0.0_dp

        ! Workspace allocation
        allocate(HT_B(N, N))
        allocate(work_B(lwork))
        allocate(rwork_B(7*N))
        allocate(iwork_B(5*N))
        allocate(ifail_B(N))
        HT_B = (0.0_dp, 0.0_dp)

        do iB = 1, nB
          B_val = B_min + (iB - 1) * B_step

          ! Update B_vec for this sweep point
          cfg%bdg%B_vec(3) = B_val

          ! Build Landau Hamiltonian with updated B field
          call ZB8bandLandau(HT_B, wv0, profile, kpterms, cfg%grid%x, cfg=cfg)

          ! Diagonalize (eigenvalues only for speed)
          call zheevx('N', 'I', 'U', N, HT_B, N, vl, vu, il, iuu, abstol, M_B, &
            eig_B(:, iB), HT_B, N, work_B, lwork, rwork_B, iwork_B, &
            ifail_B, info_B)
          if (info_B /= 0) then
            print *, 'ERROR: B-sweep diagonalization failed at B=', B_val, &
              & ' info=', info_B
            stop 1
          end if
        end do

        ! Write fan diagram to output/landau_fan.dat
        call ensure_output_dir()
        call get_unit(iounit)
        open(unit=iounit, file='output/landau_fan.dat', status='replace', action='write')

        ! Header: # B[T] E_0 E_1 ... E_{nL-1}
        write(fmt_str, '(A,I0,A)') '(g14.6,', nL, '(1x,g14.6))'
        write(iounit, '(A)', advance='no') '# B[T]'
        do i = 1, nL
          write(iounit, '(A,I0)', advance='no') ' E_', i - 1
        end do
        write(iounit, '(A)') ''

        do iB = 1, nB
          B_val = B_min + (iB - 1) * B_step
          write(iounit, fmt_str) B_val, (eig_B(i, iB), i=1, nL)
        end do
        close(iounit)

        print '(A,I0,A)', ' Landau fan diagram written to output/landau_fan.dat (', nB, ' B-points)'

        deallocate(eig_B, HT_B, work_B, rwork_B, iwork_B, ifail_B)
      end block
    end if

  end if
  call writeEigenvalues(smallk, eig(:,:), cfg%wave_vector%nsteps)


  ! ====================================================================
  ! Exciton binding energy (standalone)
  ! ====================================================================
  if (cfg%exciton%enabled .and. conf_direction(cfg%confinement) == 'z') then
    block
      real(kind=dp) :: E_binding_sa, lambda_opt_sa
      call compute_exciton_binding(eig(:, 1), eigv(:, :, 1), &
        & cfg%z, cfg%dz, cfg%num_layers, cfg%params, &
        & cfg%bands%num_cb, cfg%bands%num_vb, cfg%ngrid, E_binding_sa, lambda_opt_sa, &
        & cfg%grid%material_id)
      print '(A,F8.3,A)', ' Exciton binding energy: ', E_binding_sa, ' meV'
      print '(A,F8.2,A)', ' Variational parameter: ', lambda_opt_sa, ' AA'
    end block
  end if

  ! ====================================================================
  ! LO-phonon scattering rates (QW only, k=0)
  ! ====================================================================
  if (cfg%scattering%enabled .and. conf_direction(cfg%confinement) == 'z') then
    call compute_phonon_scattering(cfg, eig(:, 1), eigv(:, :, 1), &
      & cfg%z, cfg%params, cfg%dz, cfg%bands%num_cb, cfg%bands%num_vb, cfg%ngrid)
    print '(A)', ' Scattering rates written to output/scattering_rates.dat'
  end if


  if (allocated(smallk)) deallocate(smallk)
  if (allocated(profile)) deallocate(profile)
  if (allocated(kpterms)) deallocate(kpterms)
  if (allocated(HT)) deallocate(HT)
  if (allocated(HTmp)) deallocate(HTmp)
  if (allocated(eig)) deallocate(eig)
  if (allocated(eigv)) deallocate(eigv)
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(ifail)) deallocate(ifail)


contains

  subroutine reorder_wire_branches(prev_eval, prev_evec, curr_eval, curr_evec)
    real(kind=dp), intent(in) :: prev_eval(:)
    complex(kind=dp), intent(in), contiguous :: prev_evec(:,:)
    real(kind=dp), intent(inout) :: curr_eval(:)
    complex(kind=dp), intent(inout), contiguous :: curr_evec(:,:)

    integer :: n_prev, n_curr, n_match, i, j, best_j, next_slot
    integer :: block_start, block_end
    real(kind=dp) :: best_score, score, overlap_mag, energy_delta, parts_similarity
    logical, allocatable :: used(:)
    integer, allocatable :: perm(:)
    real(kind=dp), allocatable :: eval_tmp(:)
    complex(kind=dp), allocatable :: evec_tmp(:,:)
    real(kind=dp), allocatable :: prev_parts(:,:), curr_parts(:,:)
    real(kind=dp), parameter :: degeneracy_tol = 1.0e-8_dp

    n_prev = size(prev_eval)
    n_curr = size(curr_eval)
    if (n_prev <= 0 .or. n_curr <= 0) return

    n_match = min(n_prev, n_curr)
    allocate(used(n_curr), perm(n_curr))
    allocate(prev_parts(8, n_prev), curr_parts(8, n_curr))
    used = .false.
    perm = 0

    do i = 1, n_prev
      call compute_band_parts(prev_evec(:, i), prev_parts(:, i))
    end do
    do j = 1, n_curr
      call compute_band_parts(curr_evec(:, j), curr_parts(:, j))
    end do

    do i = 1, n_match
      best_j = 0
      best_score = -huge(1.0_dp)
      do j = 1, n_curr
        if (used(j)) cycle
        overlap_mag = abs(sum(conjg(prev_evec(:, i)) * curr_evec(:, j)))
        energy_delta = abs(curr_eval(j) - prev_eval(i))
        parts_similarity = sum(prev_parts(:, i) * curr_parts(:, j))
        score = overlap_mag + 0.25_dp * parts_similarity - 1.0e-3_dp * energy_delta
        if (score > best_score) then
          best_score = score
          best_j = j
        end if
      end do
      if (best_j > 0) then
        perm(i) = best_j
        used(best_j) = .true.
      end if
    end do

    next_slot = n_match + 1
    do j = 1, n_curr
      if (.not. used(j)) then
        perm(next_slot) = j
        next_slot = next_slot + 1
      end if
    end do

    ! FEAST and dense LAPACK can choose different bases inside an exactly
    ! degenerate manifold at the previous k-point. Keep the matched current
    ! states in ascending-energy order inside those blocks so the first step
    ! away from k=0 remains deterministic across solvers.
    block_start = 1
    do while (block_start <= n_match)
      block_end = block_start
      do while (block_end < n_match)
        if (abs(prev_eval(block_end + 1) - prev_eval(block_start)) > degeneracy_tol) exit
        block_end = block_end + 1
      end do
      if (block_end > block_start) then
        call sort_perm_block_by_energy(perm(block_start:block_end), curr_eval)
      end if
      block_start = block_end + 1
    end do

    allocate(eval_tmp(n_curr))
    allocate(evec_tmp(size(curr_evec, 1), n_curr))
    do i = 1, n_curr
      eval_tmp(i) = curr_eval(perm(i))
      evec_tmp(:, i) = curr_evec(:, perm(i))
    end do

    curr_eval = eval_tmp
    curr_evec = evec_tmp

    deallocate(eval_tmp, evec_tmp, prev_parts, curr_parts, perm, used)
  end subroutine reorder_wire_branches

  subroutine sort_perm_block_by_energy(block_perm, curr_eval)
    integer, intent(inout) :: block_perm(:)
    real(kind=dp), intent(in) :: curr_eval(:)

    integer :: i, j, tmp_idx

    do i = 1, size(block_perm) - 1
      do j = i + 1, size(block_perm)
        if (curr_eval(block_perm(j)) < curr_eval(block_perm(i))) then
          tmp_idx = block_perm(i)
          block_perm(i) = block_perm(j)
          block_perm(j) = tmp_idx
        end if
      end do
    end do
  end subroutine sort_perm_block_by_energy

end program kpfdm
