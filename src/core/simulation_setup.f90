module simulation_setup_mod

  use definitions
  use parameters
  use hamiltonianConstructor
  use confinement_init, only: confinementInitialization, confinementInitialization_2d
  use hamiltonian_wire, only: wire_workspace, wire_workspace_free, &
    & wire_coo_cache, wire_coo_cache_free, ZB8bandGeneralized, &
    & build_velocity_matrices
  use strain_solver, only: compute_strain, compute_bir_pikus_blocks, &
    & strain_result, strain_result_free
  use sc_loop, only: self_consistent_loop, self_consistent_loop_wire
  use eigensolver, only: eigensolver_base, make_eigensolver, eigensolver_config, &
    & eigensolver_result, eigensolver_result_free, auto_compute_energy_window
  use linalg, only: zheev, zheevd
  use utils
  use sparse_matrices
  use geometry, only: init_wire_from_config

  implicit none

  private
  public :: simulation_setup
  public :: simulation_setup_init, simulation_setup_free
  public :: setup_build_H, setup_solve_kpoint_serial
  public :: setup_build_velocity_matrices
  public :: thread_workspace, setup_alloc_sweep

  type :: simulation_setup
    character(len=8) :: confinement = 'none'
    logical :: has_strain = .false.
    logical :: sc_was_run = .false.
    real(kind=dp) :: fermi_level = 0.0_dp
    real(kind=dp), allocatable :: profile(:,:)
    real(kind=dp), allocatable :: kpterms(:,:,:)
    integer :: N = 0
    integer :: il = 0
    integer :: iuu = 0
    complex(kind=dp), allocatable :: HT(:,:)
    complex(kind=dp), allocatable :: work(:)
    real(kind=dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    integer :: lwork = 0
    type(csr_matrix) :: vel(3)
    logical :: vel_built = .false.
    real(kind=dp), allocatable    :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(csr_matrix), allocatable :: HT_csr_ptr
    type(wire_coo_cache), allocatable :: coo_cache_ptr
    type(wire_workspace), allocatable :: wire_ws_ptr
    type(eigensolver_config)      :: eigen_cfg
    class(eigensolver_base), allocatable :: eigen_solver
    integer :: Ngrid = 0, Ntot = 0, nev_wire = 0
  contains
    final :: simulation_setup_finalize
  end type simulation_setup

  type :: thread_workspace
    type(csr_matrix) :: HT_step
    type(eigensolver_config) :: cfg
    class(eigensolver_base), allocatable :: solver
  contains
    final :: thread_workspace_finalize
  end type thread_workspace

contains

  subroutine simulation_setup_init(cfg, setup, strain_out, sc_phi_out, sc_ne_out, sc_nh_out, skip_sc)
    type(simulation_config), intent(inout) :: cfg
    type(simulation_setup), intent(inout) :: setup
    type(strain_result), intent(out), optional, allocatable :: strain_out
    real(kind=dp), allocatable, intent(out), optional :: sc_phi_out(:,:), sc_ne_out(:,:), sc_nh_out(:,:)
    logical, intent(in), optional :: skip_sc

    integer :: info, lrwork, liwork, Ngrid_local, Ntot_local, nev
    real(kind=dp), allocatable :: eig_tmp(:,:)
    logical :: do_skip_sc

    do_skip_sc = .false.
    if (present(skip_sc)) do_skip_sc = skip_sc

    setup%confinement = cfg%confinement
    setup%has_strain = .false.
    setup%sc_was_run = .false.
    setup%fermi_level = 0.0_dp
    setup%vel_built = .false.

    select case(trim(cfg%confinement))
    case('bulk')
      if (cfg%sc%enabled == 1) then
        print *, 'Warning: SC with bulk mode requires confinement=qw.'
      end if
      setup%N = 8
      setup%il = 1
      setup%iuu = 8
      setup%lwork = 0

    case('qw')
      setup%N = cfg%grid%npoints() * 8
      setup%il = cfg%bands%num_vb + 1
      setup%iuu = setup%N
      allocate(setup%kpterms(grid_ngrid(cfg%grid), grid_ngrid(cfg%grid), 10))
      setup%kpterms = 0.0_dp
      call confinementInitialization(cfg, setup%profile, setup%kpterms)
      if (cfg%external_field%enabled .and. cfg%external_field%type == "EF") then
        if (abs(cfg%z(1)) < tolerance) then
          print *, 'Error: Electric field requires z(1) /= 0.'
          stop 1
        end if
        call externalFieldSetup_electricField(setup%profile, cfg%external_field%value, &
          cfg%totalSize, cfg%z)
      end if
      if (cfg%strain%enabled) then
        block
          type(strain_result), allocatable :: strain_local
          real(kind=dp) :: a0_ref
          allocate(strain_local)
          a0_ref = cfg%params(1)%a0
          call compute_strain(cfg%grid, cfg%params, cfg%grid%material_id, &
            cfg%strain, a0_ref, strain_local)
          call compute_bir_pikus_blocks(strain_local, cfg%params, &
            cfg%grid%material_id, cfg%grid, cfg%strain_blocks)
          if (present(strain_out)) then
            call move_alloc(strain_local, strain_out)
          else
            call strain_result_free(strain_local)
          end if
          setup%has_strain = .true.
        end block
      end if
      if (cfg%sc%enabled == 1 .and. .not. do_skip_sc) then
        block
          complex(kind=dp), allocatable :: eigv_sc(:,:,:)
          real(kind=dp), allocatable :: sc_ne_local(:), sc_nh_local(:)
          type(wavevector), allocatable :: smallk_arr(:)
          allocate(smallk_arr(1))
          smallk_arr(1)%kx = 0.0_dp
          smallk_arr(1)%ky = 0.0_dp
          smallk_arr(1)%kz = 0.0_dp
          allocate(eigv_sc(1, 1, 1))
          allocate(setup%HT(setup%N, setup%N))
          allocate(eig_tmp(setup%N, 1))
          call self_consistent_loop(setup%profile, cfg, setup%kpterms, &
            setup%HT, eig_tmp, eigv_sc, &
            smallk_arr, setup%N, 1, setup%N, &
            n_electron_out=sc_ne_local, n_hole_out=sc_nh_local, &
            fermi_level_out=setup%fermi_level)
          setup%sc_was_run = .true.
          if (allocated(eigv_sc)) deallocate(eigv_sc)
          if (allocated(sc_ne_local)) deallocate(sc_ne_local)
          if (allocated(sc_nh_local)) deallocate(sc_nh_local)
          if (allocated(eig_tmp)) deallocate(eig_tmp)
          if (allocated(smallk_arr)) deallocate(smallk_arr)
        end block
      end if
      if (.not. allocated(setup%HT)) then
        allocate(setup%HT(setup%N, setup%N))
      end if
      allocate(setup%work(1))
      allocate(setup%rwork(1))
      allocate(setup%iwork(1))
      setup%HT = 0.0_dp
      allocate(eig_tmp(setup%N, 1))
      eig_tmp = 0.0_dp
      if (cfg%num_layers == 1) then
        call zheev('V', 'L', setup%N, setup%HT, setup%N, eig_tmp(:,1), &
          setup%work, -1, setup%rwork, info)
        if (info /= 0) then
          print *, 'Error: zheev workspace query failed, info =', info
          stop 1
        end if
        setup%lwork = int(real(setup%work(1)))
        deallocate(setup%work)
        allocate(setup%work(setup%lwork))
        deallocate(setup%rwork)
        allocate(setup%rwork(max(1, 3*setup%N - 2)))
      else
        call zheevd('V', 'U', setup%N, setup%HT, setup%N, eig_tmp(:,1), &
          setup%work, -1, setup%rwork, -1, setup%iwork, -1, info)
        if (info /= 0) then
          print *, 'Error: zheevd workspace query failed, info =', info
          stop 1
        end if
        setup%lwork = int(real(setup%work(1)))
        lrwork = int(real(setup%rwork(1)))
        liwork = setup%iwork(1)
        deallocate(setup%work)
        allocate(setup%work(setup%lwork))
        deallocate(setup%rwork)
        allocate(setup%rwork(lrwork))
        deallocate(setup%iwork)
        allocate(setup%iwork(liwork))
      end if
      deallocate(eig_tmp)

    case('wire')
      if (.not. allocated(cfg%grid%x)) call init_wire_from_config(cfg)
      call confinementInitialization_2d(cfg%grid, cfg%params, cfg%wire%regions, &
        setup%profile_2d, setup%kpterms_2d, cfg%FDorder)
      if (cfg%strain%enabled) then
        block
          type(strain_result), allocatable :: strain_local
          real(kind=dp) :: a0_ref
          allocate(strain_local)
          a0_ref = cfg%params(1)%a0
          call compute_strain(cfg%grid, cfg%params, cfg%grid%material_id, &
            cfg%strain, a0_ref, strain_local)
          call compute_bir_pikus_blocks(strain_local, cfg%params, &
            cfg%grid%material_id, cfg%grid, cfg%strain_blocks)
          if (present(strain_out)) then
            call move_alloc(strain_local, strain_out)
          else
            call strain_result_free(strain_local)
          end if
          setup%has_strain = .true.
        end block
      end if
      Ngrid_local = grid_ngrid(cfg%grid)
      Ntot_local = 8 * Ngrid_local
      setup%Ngrid = Ngrid_local
      setup%Ntot = Ntot_local
      setup%N = Ntot_local
      setup%il = 1
      setup%iuu = Ntot_local
      nev = cfg%bands%num_cb + cfg%bands%num_vb
      if (nev > Ntot_local) nev = Ntot_local
      setup%nev_wire = nev
      if (cfg%feast%m0 < 0) then
        setup%eigen_cfg%method = 'DENSE'
      else
        setup%eigen_cfg%method = 'FEAST'
      end if
      setup%eigen_cfg%nev = nev
      setup%eigen_cfg%max_iter = 100
      setup%eigen_cfg%tol = 1.0e-10_dp
      setup%eigen_cfg%feast_m0 = cfg%feast%m0
      allocate(setup%HT_csr_ptr)
      allocate(setup%coo_cache_ptr)
      allocate(setup%wire_ws_ptr)
      ! Build preliminary Hamiltonian at kz=0 to estimate FEAST energy window
      call ZB8bandGeneralized(setup%HT_csr_ptr, 0.0_dp, setup%profile_2d, &
        setup%kpterms_2d, cfg, ws=setup%wire_ws_ptr)
      if (cfg%feast%emin /= 0.0_dp .or. cfg%feast%emax /= 0.0_dp) then
        setup%eigen_cfg%emin = cfg%feast%emin
        setup%eigen_cfg%emax = cfg%feast%emax
      else
        call auto_compute_energy_window(setup%HT_csr_ptr, &
          setup%eigen_cfg%emin, setup%eigen_cfg%emax)
        print *, '  Auto FEAST window: [', setup%eigen_cfg%emin, ',', setup%eigen_cfg%emax, ']'
      end if
      ! Free CSR data but preserve COO cache for fast rebuild by apps
      call csr_free(setup%HT_csr_ptr)
      setup%eigen_solver = make_eigensolver(setup%eigen_cfg)
      if (cfg%sc%enabled == 1 .and. .not. do_skip_sc) then
        block
          integer :: nk_sc, nev_sc
          real(kind=dp), allocatable :: eig_sc(:,:)
          complex(kind=dp), allocatable :: eigv_sc(:,:,:)
          real(kind=dp), allocatable :: sc_phi(:,:), sc_ne(:,:), sc_nh(:,:)
          nk_sc = cfg%sc%num_kpar
          if (mod(nk_sc, 2) == 0) nk_sc = nk_sc - 1
          nev_sc = nev
          allocate(eig_sc(nev_sc, nk_sc))
          allocate(eigv_sc(Ntot_local, nev_sc, nk_sc))
          eig_sc = 0.0_dp
          eigv_sc = cmplx(0.0_dp, 0.0_dp, kind=dp)
          call self_consistent_loop_wire(setup%profile_2d, cfg, &
            setup%kpterms_2d, cfg%grid, &
            setup%coo_cache_ptr, setup%eigen_cfg, eig_sc, eigv_sc, &
            phi_out=sc_phi, n_electron_out=sc_ne, n_hole_out=sc_nh)
          setup%sc_was_run = .true.
          deallocate(eig_sc, eigv_sc)
          if (present(sc_phi_out)) then
            call move_alloc(sc_phi, sc_phi_out)
          else
            if (allocated(sc_phi)) deallocate(sc_phi)
          end if
          if (present(sc_ne_out)) then
            call move_alloc(sc_ne, sc_ne_out)
          else
            if (allocated(sc_ne)) deallocate(sc_ne)
          end if
          if (present(sc_nh_out)) then
            call move_alloc(sc_nh, sc_nh_out)
          else
            if (allocated(sc_nh)) deallocate(sc_nh)
          end if
          call wire_coo_cache_free(setup%coo_cache_ptr)
          call wire_workspace_free(setup%wire_ws_ptr)
        end block
      end if

    case default
      print *, 'Error: simulation_setup_init unsupported confinement=', cfg%confinement
      stop 1
    end select
  end subroutine simulation_setup_init

  subroutine setup_build_H(setup, cfg, kvec, HT_out)
    type(simulation_setup), intent(inout) :: setup
    type(simulation_config), intent(in) :: cfg
    type(wavevector), intent(in) :: kvec
    complex(kind=dp), intent(inout), contiguous, optional :: HT_out(:,:)

    select case(setup%confinement)
    case('bulk')
      if (.not. present(HT_out)) then
        print *, 'Error: setup_build_H bulk requires HT_out'
        stop 1
      end if
      call ZB8bandBulk(HT_out, kvec, cfg%params(1), cfg=cfg)
    case('qw')
      if (.not. present(HT_out)) then
        print *, 'Error: setup_build_H QW requires HT_out'
        stop 1
      end if
      call ZB8bandQW(HT_out, kvec, setup%profile, setup%kpterms, cfg=cfg)
    case('wire')
      call ZB8bandGeneralized(setup%HT_csr_ptr, kvec%kz, &
        setup%profile_2d, setup%kpterms_2d, cfg, ws=setup%wire_ws_ptr)
    case default
      print *, 'Error: setup_build_H unsupported confinement=', setup%confinement
      stop 1
    end select
  end subroutine setup_build_H

  subroutine setup_solve_kpoint_serial(setup, cfg, kvec, evals, evecs)
    type(simulation_setup), intent(inout) :: setup
    type(simulation_config), intent(in) :: cfg
    type(wavevector), intent(in) :: kvec
    real(kind=dp), intent(out), contiguous :: evals(:)
    complex(kind=dp), intent(out), contiguous :: evecs(:,:)
    integer :: info
    real(kind=dp), allocatable :: rwork_local(:)

    select case(setup%confinement)
    case('bulk')
      if (.not. allocated(setup%HT)) allocate(setup%HT(8, 8))
      setup%HT = 0.0_dp
      call ZB8bandBulk(setup%HT, kvec, cfg%params(1), cfg=cfg)
      if (setup%lwork == 0) then
        allocate(setup%work(1))
        allocate(rwork_local(max(1, 3*8 - 2)))
        call zheev('V', 'L', 8, setup%HT, 8, evals, setup%work, -1, rwork_local, info)
        if (info /= 0) then
          print *, 'Error: zheev workspace query failed in setup_solve_kpoint_serial, info =', info
          stop 1
        end if
        setup%lwork = int(real(setup%work(1)))
        deallocate(setup%work)
        allocate(setup%work(setup%lwork))
        deallocate(rwork_local)
        allocate(setup%rwork(max(1, 3*8 - 2)))
      end if
      setup%HT = 0.0_dp
      call ZB8bandBulk(setup%HT, kvec, cfg%params(1), cfg=cfg)
      call zheev('V', 'L', 8, setup%HT, 8, evals, setup%work, setup%lwork, setup%rwork, info)
      if (info /= 0) then
        print *, 'Error: zheev failed in setup_solve_kpoint_serial, info =', info
        stop 1
      end if
      evecs = setup%HT
    case('qw')
      setup%HT = 0.0_dp
      call ZB8bandQW(setup%HT, kvec, setup%profile, setup%kpterms, cfg=cfg)
      if (cfg%num_layers == 1) then
        call zheev('V', 'L', setup%N, setup%HT, setup%N, evals, setup%work, setup%lwork, setup%rwork, info)
        if (info /= 0) then
          print *, 'Error: zheev failed in setup_solve_kpoint_serial, info =', info
          stop 1
        end if
      else
        call zheevd('V', 'U', setup%N, setup%HT, setup%N, evals, setup%work, setup%lwork, setup%rwork, size(setup%rwork), setup%iwork, size(setup%iwork), info)
        if (info /= 0) then
          print *, 'Error: zheevd failed in setup_solve_kpoint_serial, info =', info
          stop 1
        end if
      end if
      evecs = setup%HT
    case default
      print *, 'Error: setup_solve_kpoint_serial unsupported confinement=', setup%confinement
      stop 1
    end select
  end subroutine setup_solve_kpoint_serial

  subroutine setup_build_velocity_matrices(setup, cfg)
    type(simulation_setup), intent(inout) :: setup
    type(simulation_config), intent(in) :: cfg
    complex(kind=dp), allocatable :: vel_dense(:,:)
    type(csr_matrix) :: H_csr_tmp
    type(wavevector) :: kv_g

    if (setup%confinement == 'wire') then
      if (.not. allocated(setup%HT_csr_ptr%rowptr)) then
        print *, 'Error: call setup_build_H before setup_build_velocity_matrices for wire'
        stop 1
      end if
    end if

    select case(setup%confinement)
    case('bulk')
      allocate(vel_dense(8, 8))
      vel_dense = 0.0_dp; kv_g%kx = 1.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandBulk(vel_dense, kv_g, cfg%params(1), cfg=cfg, g='g')
      call dense_to_csr_8(vel_dense, setup%vel(1))
      vel_dense = 0.0_dp; kv_g%kx = 0.0_dp; kv_g%ky = 1.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandBulk(vel_dense, kv_g, cfg%params(1), cfg=cfg, g='g')
      call dense_to_csr_8(vel_dense, setup%vel(2))
      vel_dense = 0.0_dp; kv_g%kx = 0.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 1.0_dp
      call ZB8bandBulk(vel_dense, kv_g, cfg%params(1), cfg=cfg, g='g')
      call dense_to_csr_8(vel_dense, setup%vel(3))
      deallocate(vel_dense)
    case('qw')
      ! vel(3): commutator -i[r_z, H] via build_velocity_matrices (1D QW)
      ! vel(1), vel(2): overwritten with k-derivative dH/dkx, dH/dky
      allocate(vel_dense(setup%N, setup%N))
      vel_dense = 0.0_dp
      kv_g%kx = 0.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandQW(vel_dense, kv_g, setup%profile, setup%kpterms, cfg=cfg)
      call dense_to_csr_N(vel_dense, H_csr_tmp)
      call build_velocity_matrices(H_csr_tmp, cfg%grid, setup%vel)
      call csr_free(H_csr_tmp)
      ! Overwrite vel(1) with k-derivative (x-direction)
      vel_dense = 0.0_dp
      kv_g%kx = 1.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandQW(vel_dense, kv_g, setup%profile, setup%kpterms, cfg=cfg, g='g')
      call csr_free(setup%vel(1))
      call dense_to_csr_N(vel_dense, setup%vel(1))
      ! Overwrite vel(2) with k-derivative (y-direction)
      vel_dense = 0.0_dp
      kv_g%kx = 0.0_dp; kv_g%ky = 1.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandQW(vel_dense, kv_g, setup%profile, setup%kpterms, cfg=cfg, g='g')
      call csr_free(setup%vel(2))
      call dense_to_csr_N(vel_dense, setup%vel(2))
      deallocate(vel_dense)
    case('wire')
      call build_velocity_matrices(setup%HT_csr_ptr, cfg%grid, setup%vel(1), setup%vel(2))
      call ZB8bandGeneralized(setup%vel(3), 1.0_dp, setup%profile_2d, &
        setup%kpterms_2d, cfg, g='g3')
    case default
      print *, 'Error: setup_build_velocity_matrices unsupported confinement=', setup%confinement
      stop 1
    end select
    setup%vel_built = .true.
  end subroutine setup_build_velocity_matrices

  subroutine simulation_setup_free(setup)
    type(simulation_setup), intent(inout) :: setup
    integer :: i
    if (allocated(setup%profile)) deallocate(setup%profile)
    if (allocated(setup%kpterms)) deallocate(setup%kpterms)
    if (allocated(setup%HT)) deallocate(setup%HT)
    if (allocated(setup%work)) deallocate(setup%work)
    if (allocated(setup%rwork)) deallocate(setup%rwork)
    if (allocated(setup%iwork)) deallocate(setup%iwork)
    if (allocated(setup%profile_2d)) deallocate(setup%profile_2d)
    if (allocated(setup%kpterms_2d)) then
      do i = 1, size(setup%kpterms_2d)
        call csr_free(setup%kpterms_2d(i))
      end do
      deallocate(setup%kpterms_2d)
    end if
    if (allocated(setup%HT_csr_ptr)) then
      call csr_free(setup%HT_csr_ptr)
      deallocate(setup%HT_csr_ptr)
    end if
    if (allocated(setup%coo_cache_ptr)) then
      call wire_coo_cache_free(setup%coo_cache_ptr)
      deallocate(setup%coo_cache_ptr)
    end if
    if (allocated(setup%wire_ws_ptr)) then
      call wire_workspace_free(setup%wire_ws_ptr)
      deallocate(setup%wire_ws_ptr)
    end if
    if (allocated(setup%eigen_solver)) deallocate(setup%eigen_solver)
    if (setup%vel_built) then
      do i = 1, 3
        call csr_free(setup%vel(i))
      end do
    end if
    setup%N = 0
    setup%confinement = 'none'
    setup%lwork = 0
    setup%il = 0
    setup%iuu = 0
    setup%has_strain = .false.
    setup%sc_was_run = .false.
    setup%vel_built = .false.
    setup%fermi_level = 0.0_dp
    setup%Ngrid = 0
    setup%Ntot = 0
    setup%nev_wire = 0
  end subroutine simulation_setup_free



  subroutine simulation_setup_finalize(setup)
    type(simulation_setup), intent(inout) :: setup
    call simulation_setup_free(setup)
  end subroutine simulation_setup_finalize
  subroutine setup_alloc_sweep(setup, nthreads, tws)
    type(simulation_setup), intent(in) :: setup
    integer, intent(in) :: nthreads
    type(thread_workspace), allocatable, intent(out) :: tws(:)
    integer :: t
    if (setup%confinement /= 'wire') then
      print *, 'Error: setup_alloc_sweep requires confinement=2 (wire)'
      stop 1
    end if
    allocate(tws(nthreads))
    do t = 1, nthreads
      call csr_clone_structure(setup%HT_csr_ptr, tws(t)%HT_step)
      tws(t)%cfg = setup%eigen_cfg
      tws(t)%solver = make_eigensolver(tws(t)%cfg)
    end do
  end subroutine setup_alloc_sweep



  subroutine thread_workspace_finalize(tw)
    type(thread_workspace), intent(inout) :: tw
    call csr_free(tw%HT_step)
    if (allocated(tw%solver)) deallocate(tw%solver)
  end subroutine thread_workspace_finalize
  subroutine dense_to_csr_8(dns, csr_out)
    complex(kind=dp), intent(in), contiguous :: dns(:,:)
    type(csr_matrix), intent(out) :: csr_out
    integer :: nnz, i, j, idx
    nnz = 0
    do j = 1, 8
      do i = 1, 8
        if (abs(dns(i, j)) > 0.0_dp) nnz = nnz + 1
      end do
    end do
    csr_out%nrows = 8
    csr_out%ncols = 8
    csr_out%nnz = nnz
    allocate(csr_out%rowptr(9))
    allocate(csr_out%colind(nnz))
    allocate(csr_out%values(nnz))
    idx = 0
    csr_out%rowptr(1) = 1
    do i = 1, 8
      do j = 1, 8
        if (abs(dns(i, j)) > 0.0_dp) then
          idx = idx + 1
          csr_out%colind(idx) = j
          csr_out%values(idx) = dns(i, j)
        end if
      end do
      csr_out%rowptr(i + 1) = idx + 1
    end do
  end subroutine dense_to_csr_8

  subroutine dense_to_csr_N(dns, csr_out)
    complex(kind=dp), intent(in), contiguous :: dns(:,:)
    type(csr_matrix), intent(out) :: csr_out
    integer :: nnz, N, i, j, idx
    N = size(dns, 1)
    nnz = 0
    do j = 1, N
      do i = 1, N
        if (abs(dns(i, j)) > 0.0_dp) nnz = nnz + 1
      end do
    end do
    csr_out%nrows = N
    csr_out%ncols = N
    csr_out%nnz = nnz
    allocate(csr_out%rowptr(N + 1))
    allocate(csr_out%colind(nnz))
    allocate(csr_out%values(nnz))
    idx = 0
    csr_out%rowptr(1) = 1
    do i = 1, N
      do j = 1, N
        if (abs(dns(i, j)) > 0.0_dp) then
          idx = idx + 1
          csr_out%colind(idx) = j
          csr_out%values(idx) = dns(i, j)
        end if
      end do
      csr_out%rowptr(i + 1) = idx + 1
    end do
  end subroutine dense_to_csr_N

end module simulation_setup_mod
