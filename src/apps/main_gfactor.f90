program gfactor

  use definitions, only: NUM_CB_STATES, NUM_VB_STATES, conf_direction, &
    dp, optical_transition, simulation_config, validate_semantic, &
    wavevector
  use parameters
  use hamiltonianConstructor
  use outputFunctions
  use gfactorFunctions
  use utils, only: get_unit, ensure_output_dir
  use input_parser
  use sparse_matrices
  use eigensolver, only: eigensolver_result, eigensolver_result_free
  use linalg, only: mkl_set_num_threads_local
  use spin_projection, only: compute_band_parts
  use simulation_setup_mod

  implicit NONE

  ! Shared configuration from input_parser
  type(simulation_config) :: cfg

  ! wave vector
  type(wavevector), allocatable, dimension(:) :: smallk

  ! iteration consts
  integer :: i, j, k

  ! hamiltonian and LAPACK/BLAS
  integer :: info, lwork, N, lrwork, liwork
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  complex(kind=dp), allocatable :: work(:)
  integer, allocatable :: iwork(:)
  complex(kind=dp), allocatable, dimension(:,:) :: HT

  ! gfactor
  complex(kind=dp), allocatable, dimension(:,:) :: cb_state, vb_state
  real(kind=dp), allocatable, dimension(:) :: cb_value, vb_value
  integer :: whichBand, bandIdx
  complex(kind=dp), allocatable, dimension(:,:,:) :: tensor
  real(kind=dp) :: g_eff(3)

  ! file handling
  integer(kind=4) :: iounit

  ! --- Wire mode (confinement='wire') g-factor extraction variables ---
  type(eigensolver_result)         :: eigen_res
  integer                          :: gap_idx, cb_start, vb_start
  real(kind=dp), allocatable       :: gaps(:)
  real(kind=dp), allocatable       :: band_parts(:,:), cb_char(:)
  real(kind=dp), parameter         :: cb_threshold = 0.5_dp
  integer                          :: first_cb_idx
  logical                          :: used_character_gap

  ! Ensure MKL routines run single-threaded (sequential) regardless of
  ! MKL_THREADING setting.  Prevents thread oversubscription when
  ! MKL_THREADING=intel_thread.
  info = mkl_set_num_threads_local(1)

  ! Shared setup: read input, initialize materials
  call read_config(cfg)

  ! Semantic validation (gfactor requires k0 mode)
  call validate_semantic(cfg, 'gfactor')

  ! Initialize k for workspace query
  k = 1

  ! Allocate arrays for single k-point
  if (allocated(smallk)) deallocate(smallk)
  allocate(smallk(1))
  smallk(1)%kx = 0
  smallk(1)%ky = 0
  smallk(1)%kz = 0

  whichBand = cfg%which_band
  bandIdx = cfg%band_idx

  ! ====================================================================
  ! Branch: wire (confinement='wire') vs bulk/QW (confinement='bulk','qw')
  ! ====================================================================
  if (trim(cfg%confinement) == 'wire') then

    ! ----------------------------------------------------------------
    ! Wire mode: sparse Hamiltonian + FEAST eigensolver via setup
    ! ----------------------------------------------------------------

    ! For wire mode, numcb/numvb are used directly from input (not * fdStep)
    cfg%evnum = cfg%bands%num_cb + cfg%bands%num_vb

    block
      type(simulation_setup) :: setup

      ! Initialize: confinement, strain, SC, energy window, eigensolver
      call simulation_setup_init(cfg, setup)

      print *, ''
      print *, '=== Wire g-factor (2D confinement) ==='
      print *, '  Grid: nx=', cfg%grid%nx, ' ny=', cfg%grid%ny, ' Ngrid=', setup%Ngrid
      print *, '  Matrix size: ', setup%Ntot, 'x', setup%Ntot
      print *, '  numcb=', cfg%bands%num_cb, ' numvb=', cfg%bands%num_vb
      print *, ''

      ! Build wire Hamiltonian at kz=0
      call setup_build_H(setup, cfg, smallk(1))

      ! Build commutator-based velocity matrices
      call setup_build_velocity_matrices(setup, cfg)

      ! Solve eigenvalue problem
      call setup%eigen_solver%solve_sparse(setup%HT_csr_ptr, setup%eigen_cfg, eigen_res)

      if (eigen_res%nev_found == 0) then
        print *, 'Error: ', setup%eigen_solver%backend_name(), &
          ' found no eigenvalues. Check energy window.'
        print *, '  Window: [', setup%eigen_cfg%emin, ',', setup%eigen_cfg%emax, ']'
        error stop 'eigensolver found no eigenvalues'
      end if

      if (.not. eigen_res%converged) then
        print *, 'Warning: ', setup%eigen_solver%backend_name(), &
          ' did not converge for g-factor calculation'
      end if

      if (eigen_res%nev_found < cfg%bands%num_cb + cfg%bands%num_vb) then
        print *, 'Error: ', setup%eigen_solver%backend_name(), ' found', &
          eigen_res%nev_found, &
          ' eigenvalues but need', cfg%bands%num_cb + cfg%bands%num_vb
        print *, '  numcb=', cfg%bands%num_cb, ' numvb=', cfg%bands%num_vb
        error stop 'eigensolver eigenvalue count mismatch'
      end if

      N = setup%Ntot

      ! Write wire eigenfunctions at kz=0
      call writeEigenfunctions2d(cfg%grid, eigen_res%eigenvalues, &
        & eigen_res%eigenvectors, 1, eigen_res%nev_found, write_parts=.true.)

      ! Extract CB/VB states from the band edge, not from fixed sorted positions.
      ! FEAST may return many states in the search window, so the user-supplied
      ! numvb/numcb should be interpreted as counts around the actual gap.
      ! Prefer band character over a raw largest-gap split: wide wire windows can
      ! contain remote valence manifolds whose largest spectral gap is not the
      ! fundamental VB/CB edge.
      allocate(cb_value(cfg%bands%num_cb))
      allocate(vb_value(cfg%bands%num_vb))
      allocate(cb_state(N, cfg%bands%num_cb))
      allocate(vb_state(N, cfg%bands%num_vb))

      allocate(band_parts(8, eigen_res%nev_found))
      allocate(cb_char(eigen_res%nev_found))
      do i = 1, eigen_res%nev_found
        call compute_band_parts(eigen_res%eigenvectors(:, i), band_parts(:, i))
        cb_char(i) = band_parts(7, i) + band_parts(8, i)
      end do

      first_cb_idx = 0
      do i = 1, eigen_res%nev_found
        if (cb_char(i) >= cb_threshold) then
          first_cb_idx = i
          exit
        end if
      end do

      used_character_gap = .false.
      if (first_cb_idx > 0) then
        ! Character-based detection: CB edge identified at first_cb_idx.
        ! Use this as the band edge regardless of whether num_vb/num_cb fit
        ! exactly — the g-factor should use states nearest the edge.
        gap_idx = first_cb_idx - 1
        used_character_gap = .true.
      else
        allocate(gaps(eigen_res%nev_found - 1))
        gaps(:) = eigen_res%eigenvalues(2:eigen_res%nev_found) - &
          & eigen_res%eigenvalues(1:eigen_res%nev_found-1)
        gap_idx = 0
        ! Find the largest gap that can accommodate numvb states below and numcb above
        do while (gap_idx == 0)
          gap_idx = maxloc(gaps, dim=1)
          if (gap_idx - cfg%bands%num_vb + 1 >= 1 .and. &
            & gap_idx + cfg%bands%num_cb <= eigen_res%nev_found) exit
          gaps(gap_idx) = -1.0_dp  ! disqualify and retry
          gap_idx = 0
        end do
        deallocate(gaps)
      end if

      cb_start = gap_idx + 1
      vb_start = gap_idx - cfg%bands%num_vb + 1

      if (vb_start >= 1 .and. cb_start + cfg%bands%num_cb - 1 <= eigen_res%nev_found) then
        if (used_character_gap) then
          print *, '  Wire band edge selected by band character between eigenvalues ', &
            gap_idx, ' and ', gap_idx + 1
          print *, '  Edge characters: VB top CB-weight=', cb_char(gap_idx), &
            ' CB bottom CB-weight=', cb_char(cb_start)
        else
          print *, '  Wire band edge detected by largest spectral gap between eigenvalues ', &
            gap_idx, ' and ', gap_idx + 1
        end if
        print *, '  Selecting VB indices ', vb_start, ':', gap_idx, &
          & ' and CB indices ', cb_start, ':', cb_start + cfg%bands%num_cb - 1

        cb_value(:) = eigen_res%eigenvalues(cb_start:cb_start+cfg%bands%num_cb-1)
        vb_value(:) = eigen_res%eigenvalues(gap_idx:vb_start:-1)

        cb_state(:,:) = eigen_res%eigenvectors(:, cb_start:cb_start+cfg%bands%num_cb-1)
        vb_state(:,:) = eigen_res%eigenvectors(:, gap_idx:vb_start:-1)
      else
        print *, '  WARNING: automatic wire band-edge detection did not have enough'
        print *, '           states on both sides of the gap. Falling back to the'
        print *, '           legacy contiguous selection around numvb/numcb.'
        print *, '  gap_idx=', gap_idx, ' nev_found=', eigen_res%nev_found

        cb_value(:) = eigen_res%eigenvalues(cfg%bands%num_vb+1:cfg%bands%num_vb+cfg%bands%num_cb)
        vb_value(:) = eigen_res%eigenvalues(cfg%bands%num_vb:1:-1)

        cb_state(:,:) = eigen_res%eigenvectors(:, cfg%bands%num_vb+1:cfg%bands%num_vb+cfg%bands%num_cb)
        vb_state(:,:) = eigen_res%eigenvectors(:, cfg%bands%num_vb:1:-1)
      end if

      deallocate(band_parts, cb_char)

      ! Compute g-factor tensor
      allocate(tensor(2,2,3))
      tensor = 0

      call gfactorCalculation_wire(tensor, g_eff, whichBand, bandIdx, cfg%bands%num_cb, &
        & cfg%bands%num_vb, cb_state, vb_state, cb_value, vb_value, cfg, setup%vel)

      if (cfg%optics%enabled) then
        ! Compute optical transitions only for optics-enabled runs.
        block
          type(optical_transition), allocatable :: transitions(:)
          integer :: num_trans

          call compute_optical_matrix_wire(transitions, num_trans, &
            cb_state, vb_state, cb_value, vb_value, cfg%bands%num_cb, cfg%bands%num_vb, &
            setup%vel)

          ! Write to file
          call write_optical_transitions(transitions)

          deallocate(transitions)
        end block
      end if

      ! Free eigenresult and setup resources
      call eigensolver_result_free(eigen_res)
      call simulation_setup_free(setup)
    end block

  else

    ! ----------------------------------------------------------------
    ! Bulk / QW mode: dense Hamiltonian via simulation_setup
    ! ----------------------------------------------------------------

    print *, 'adjusting numcb and numvb to get all states'
    cfg%bands%num_cb = NUM_CB_STATES*cfg%grid%npoints()
    cfg%bands%num_vb = NUM_VB_STATES*cfg%grid%npoints()

    cfg%evnum = cfg%bands%num_cb + cfg%bands%num_vb

    N = cfg%grid%npoints()*8

    if (cfg%evnum /= N) error stop 'evnum not equal to total matrix size'

    ! Initialize simulation setup (confinement, profile, kpterms, workspace)
    block
      type(simulation_setup) :: setup

      call simulation_setup_init(cfg, setup)

      ! Print profile for QW mode
      if (conf_direction(cfg%confinement) == 'z') then
        call write_profile_1d(cfg%z, setup%profile)
      end if

      ! Solve at k=0
      allocate(eig(N,1))
      allocate(HT(N,N))
      HT = 0.0_dp
      eig = 0.0_dp

      call setup_solve_gamma_point(setup, cfg, smallk(1), eig(:,1), HT)

      ! Write eigenfunctions for multi-layer QW
      if (cfg%num_layers > 1) call writeEigenfunctions(N, N, HT, 1, cfg%grid%npoints(), cfg%z, cfg%num_layers==1)

      allocate(cb_state(N,cfg%bands%num_cb))
      allocate(vb_state(N,cfg%bands%num_vb))

      allocate(cb_value(cfg%bands%num_cb))
      allocate(vb_value(cfg%bands%num_vb))

      cb_state(:,:) = HT(:, cfg%bands%num_vb+1:cfg%bands%num_vb+cfg%bands%num_cb)
      vb_state(:,:) = HT(:, cfg%bands%num_vb:1:-1)

      cb_value(:) = eig(cfg%bands%num_vb+1:cfg%bands%num_vb+cfg%bands%num_cb,1)
      vb_value(:) = eig(cfg%bands%num_vb:1:-1,1)


      allocate(tensor(2,2,3))
      tensor = 0

      if (cfg%num_layers == 1) call gfactorCalculation(tensor, g_eff, whichBand, bandIdx, cfg%bands%num_cb, &
      & cfg%bands%num_vb, cb_state, vb_state, cb_value, vb_value, cfg%num_layers, cfg%params, &
      & cfg%z_min(1), cfg%z_max(1), dz=cfg%dz)
      if (cfg%num_layers > 1)  call gfactorCalculation(tensor, g_eff, whichBand, bandIdx, cfg%bands%num_cb, &
      & cfg%bands%num_vb, cb_state, vb_state, cb_value, vb_value, cfg%num_layers, cfg%params, &
      & cfg%z_min(1), cfg%z_max(1), profile=setup%profile, kpterms=setup%kpterms, dz=cfg%dz)

      ! QW optical transitions (must be before setup_free so profile/kpterms are alive)
      if (cfg%optics%enabled .and. conf_direction(cfg%confinement) == 'z' .and. cfg%num_layers > 1) then
        block
          type(optical_transition), allocatable :: transitions(:)
          integer :: num_trans

          call compute_optical_matrix_qw(transitions, num_trans, &
            cb_state, vb_state, cb_value, vb_value, cfg%bands%num_cb, cfg%bands%num_vb, &
            cfg%num_layers, cfg%params, setup%profile, setup%kpterms, &
            cfg%z_min(1), cfg%z_max(1), cfg%dz)

          call write_optical_transitions(transitions)
          deallocate(transitions)
        end block
      end if

      call simulation_setup_free(setup)
    end block

  end if  ! confinement /= 2

  print *, 'tensor'
  do i = 1, 2
    write(*,'(2(f9.4,1x,f9.4))') (tensor(i,j,1), j=1,2)
  end do
  print *, ' '
  do i = 1, 2
    write(*,'(2(f9.4,1x,f9.4))') (tensor(i,j,2), j=1,2)
  end do
  print *, ' '
  do i = 1, 2
    write(*,'(2(f9.4,1x,f9.4))') (tensor(i,j,3), j=1,2)
  end do

  print *, 'gx'
  print *, 0.0_dp
  print *, g_eff(1)

  print *, 'gy'
  print *, 0.0_dp
  print *, g_eff(2)

  print *, 'gz'
  print *, 0.0_dp
  print *, g_eff(3)

  call get_unit(iounit)
  open(unit=iounit, file='output/gfactor.dat', status="replace", action="write")
  write(iounit,*) g_eff(1), g_eff(2), g_eff(3)
  close(iounit)

  !----------------------------------------------------------------------------


  if (allocated(smallk)) deallocate(smallk)
  if (allocated(HT)) deallocate(HT)
  if (allocated(eig)) deallocate(eig)
  if (allocated(work)) deallocate(work)
  if (allocated(iwork)) deallocate(iwork)
  if (allocated(rwork)) deallocate(rwork)
  if (allocated(cb_state)) deallocate(cb_state)
  if (allocated(vb_state)) deallocate(vb_state)
  if (allocated(cb_value)) deallocate(cb_value)
  if (allocated(vb_value)) deallocate(vb_value)
  if (allocated(tensor)) deallocate(tensor)

contains

end program gfactor
