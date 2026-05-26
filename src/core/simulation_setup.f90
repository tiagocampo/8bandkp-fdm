module simulation_setup_mod

  use definitions
  use parameters
  use hamiltonianConstructor
  use confinement_init, only: confinementInitialization
  use strain_solver
  use sc_loop, only: self_consistent_loop
  use linalg, only: zheev, zheevd
  use utils
  use sparse_matrices

  implicit none

  private
  public :: simulation_setup
  public :: simulation_setup_init, simulation_setup_free
  public :: setup_build_H, setup_solve_kpoint_serial
  public :: setup_build_velocity_matrices

  type :: simulation_setup
    integer :: confinement = -1   ! 0=bulk, 1=QW, 2=wire
    ! Common
    logical :: has_strain = .false.
    logical :: sc_was_run = .false.
    real(kind=dp) :: fermi_level = 0.0_dp
    ! Dense path (bulk/QW)
    real(kind=dp), allocatable :: profile(:,:)
    real(kind=dp), allocatable :: kpterms(:,:,:)
    integer :: N = 0
    integer :: il = 0
    integer :: iuu = 0
    complex(kind=dp), allocatable :: HT(:,:)
    ! LAPACK workspace
    complex(kind=dp), allocatable :: work(:)
    real(kind=dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    integer :: lwork = 0
    ! Velocity matrices (opt-in)
    type(csr_matrix) :: vel(3)
    logical :: vel_built = .false.
  contains
    final :: simulation_setup_finalize
  end type simulation_setup

contains

  ! ==================================================================
  ! Initialize simulation_setup from a parsed configuration.
  ! Dispatches on cfg%confinement: bulk (0), QW (1).
  ! Wire (2) is not yet supported — will be added in issue 03.
  ! ==================================================================
  subroutine simulation_setup_init(cfg, setup)
    type(simulation_config), intent(inout) :: cfg
    type(simulation_setup), intent(inout) :: setup

    integer :: info, lrwork, liwork
    real(kind=dp), allocatable :: eig_tmp(:,:)

    setup%confinement = cfg%confinement
    setup%has_strain = .false.
    setup%sc_was_run = .false.
    setup%fermi_level = 0.0_dp
    setup%vel_built = .false.

    select case(cfg%confinement)
    case(0)
      ! ---- Bulk mode: 8x8 matrix ----
      if (cfg%sc%enabled == 1) then
        print *, 'Warning: SC with bulk mode requires confinement=1, confDir=z.'
        print *, '  Skipping SC. Use confinement=1, confDir=z, numLayers=1 for delta-doped bulk.'
      end if
      setup%N = 8
      setup%il = 1
      setup%iuu = 8
      setup%lwork = 0  ! will be allocated on first solve

    case(1)
      ! ---- QW mode: 8*Nz dense matrix ----
      setup%N = cfg%fdStep * 8
      setup%il = cfg%numvb + 1
      setup%iuu = setup%N

      ! Confinement initialization (cfg is intent(inout) for the _cfg variant)
      allocate(setup%kpterms(grid_ngrid(cfg%grid), grid_ngrid(cfg%grid), 10))
      setup%kpterms = 0.0_dp
      call confinementInitialization(cfg, setup%profile, setup%kpterms)

      ! Electric field
      if (cfg%ExternalField == 1 .and. cfg%EFtype == "EF") then
        if (abs(cfg%z(1)) < tolerance) then
          print *, 'Error: Electric field requires z(1) /= 0.'
          print *, '  Adjust startPos/endPos so grid does not start at z=0.'
          stop 1
        end if
        call externalFieldSetup_electricField(setup%profile, cfg%Evalue, &
          cfg%totalSize, cfg%z)
      end if

      ! Strain
      if (cfg%strain%enabled) then
        block
          type(strain_result) :: strain_out
          real(kind=dp) :: a0_ref

          a0_ref = cfg%params(1)%a0
          call compute_strain(cfg%grid, cfg%params, cfg%grid%material_id, &
            cfg%strain, a0_ref, strain_out)
          call compute_bir_pikus_blocks(strain_out, cfg%params, &
            cfg%grid%material_id, cfg%grid, cfg%strain_blocks)
          call strain_result_free(strain_out)
          setup%has_strain = .true.
        end block
      end if

      ! Self-consistent loop
      if (cfg%sc%enabled == 1) then
        block
          complex(kind=dp), allocatable :: eigv_sc(:,:,:)
          real(kind=dp), allocatable :: sc_ne_out(:), sc_nh_out(:)
          type(wavevector), allocatable :: smallk_arr(:)

          allocate(smallk_arr(1))
          smallk_arr(1)%kx = 0.0_dp
          smallk_arr(1)%ky = 0.0_dp
          smallk_arr(1)%kz = 0.0_dp

          ! Allocate eigenvector dummy for SC
          allocate(eigv_sc(1, 1, 1))

          allocate(setup%HT(setup%N, setup%N))
          allocate(eig_tmp(setup%N, 1))

          print *, ''
          print *, '=== Running QW self-consistent SP loop ==='
          call self_consistent_loop(setup%profile, cfg, setup%kpterms, &
            setup%HT, eig_tmp, eigv_sc, &
            smallk_arr, setup%N, 1, setup%N, &
            n_electron_out=sc_ne_out, n_hole_out=sc_nh_out, &
            fermi_level_out=setup%fermi_level)
          print *, '  QW SC complete.'

          setup%sc_was_run = .true.
          if (allocated(eigv_sc)) deallocate(eigv_sc)
          if (allocated(sc_ne_out)) deallocate(sc_ne_out)
          if (allocated(sc_nh_out)) deallocate(sc_nh_out)
          if (allocated(eig_tmp)) deallocate(eig_tmp)
          if (allocated(smallk_arr)) deallocate(smallk_arr)
        end block
      end if

      ! Allocate HT workspace
      if (.not. allocated(setup%HT)) then
        allocate(setup%HT(setup%N, setup%N))
      end if

      ! LAPACK workspace query
      allocate(setup%work(1))
      allocate(setup%rwork(1))
      allocate(setup%iwork(1))

      setup%HT = 0.0_dp
      allocate(eig_tmp(setup%N, 1))
      eig_tmp = 0.0_dp

      if (cfg%numLayers == 1) then
        ! zheev workspace query
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
        ! zheevd workspace query
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

    case default
      ! Wire and other modes not yet supported
      print *, 'Error: simulation_setup_init does not yet support confinement=', &
        cfg%confinement
      stop 1
    end select

  end subroutine simulation_setup_init

  ! ==================================================================
  ! Build Hamiltonian for a given k-point.
  ! Bulk: ZB8bandBulk; QW: ZB8bandQW.
  ! ==================================================================
  subroutine setup_build_H(setup, cfg, kvec, HT_out)
    type(simulation_setup), intent(in) :: setup
    type(simulation_config), intent(in) :: cfg
    type(wavevector), intent(in) :: kvec
    complex(kind=dp), intent(inout), contiguous :: HT_out(:,:)

    select case(setup%confinement)
    case(0)
      ! Bulk
      call ZB8bandBulk(HT_out, kvec, cfg%params(1), cfg=cfg)
    case(1)
      ! QW
      call ZB8bandQW(HT_out, kvec, setup%profile, setup%kpterms, cfg=cfg)
    case default
      print *, 'Error: setup_build_H unsupported confinement=', setup%confinement
      stop 1
    end select

  end subroutine setup_build_H

  ! ==================================================================
  ! Build Hamiltonian and diagonalize via LAPACK for a single k-point.
  ! Returns eigenvalues (evals) and eigenvectors (evecs).
  ! ==================================================================
  subroutine setup_solve_kpoint_serial(setup, cfg, kvec, evals, evecs)
    type(simulation_setup), intent(inout) :: setup
    type(simulation_config), intent(in) :: cfg
    type(wavevector), intent(in) :: kvec
    real(kind=dp), intent(out), contiguous :: evals(:)
    complex(kind=dp), intent(out), contiguous :: evecs(:,:)

    integer :: info
    real(kind=dp), allocatable :: rwork_local(:)

    select case(setup%confinement)
    case(0)
      ! Bulk: always 8x8, use zheev
      if (.not. allocated(setup%HT)) then
        allocate(setup%HT(8, 8))
      end if
      setup%HT = 0.0_dp
      call ZB8bandBulk(setup%HT, kvec, cfg%params(1), cfg=cfg)

      ! Workspace query and solve for bulk
      if (setup%lwork == 0) then
        allocate(setup%work(1))
        allocate(rwork_local(max(1, 3*8 - 2)))
        call zheev('V', 'L', 8, setup%HT, 8, evals, setup%work, -1, &
          rwork_local, info)
        if (info /= 0) then
          print *, 'Error: zheev workspace query failed, info =', info
          stop 1
        end if
        setup%lwork = int(real(setup%work(1)))
        deallocate(setup%work)
        allocate(setup%work(setup%lwork))
        deallocate(rwork_local)
        allocate(setup%rwork(max(1, 3*8 - 2)))
      end if

      ! Re-build HT (workspace query destroys it)
      setup%HT = 0.0_dp
      call ZB8bandBulk(setup%HT, kvec, cfg%params(1), cfg=cfg)
      call zheev('V', 'L', 8, setup%HT, 8, evals, setup%work, setup%lwork, &
        setup%rwork, info)
      if (info /= 0) then
        print *, 'Error: zheev failed, info =', info
        stop 1
      end if
      evecs = setup%HT

    case(1)
      ! QW: use pre-allocated workspace
      setup%HT = 0.0_dp
      call ZB8bandQW(setup%HT, kvec, setup%profile, setup%kpterms, cfg=cfg)

      if (cfg%numLayers == 1) then
        call zheev('V', 'L', setup%N, setup%HT, setup%N, evals, &
          setup%work, setup%lwork, setup%rwork, info)
      else
        call zheevd('V', 'U', setup%N, setup%HT, setup%N, evals, &
          setup%work, setup%lwork, setup%rwork, size(setup%rwork), &
          setup%iwork, size(setup%iwork), info)
      end if
      if (info /= 0) then
        print *, 'Error: diagonalization failed, info =', info
        stop 1
      end if
      evecs = setup%HT

    case default
      print *, 'Error: setup_solve_kpoint_serial unsupported confinement=', &
        setup%confinement
      stop 1
    end select

  end subroutine setup_solve_kpoint_serial

  ! ==================================================================
  ! Build commutator-based velocity matrices for momentum matrix elements.
  ! Bulk: three ZB8bandBulk(g='g') calls, convert to CSR.
  ! QW: three ZB8bandQW(g='g') calls for x/y/z, convert to CSR.
  ! ==================================================================
  subroutine setup_build_velocity_matrices(setup, cfg)
    type(simulation_setup), intent(inout) :: setup
    type(simulation_config), intent(in) :: cfg

    complex(kind=dp), allocatable :: vel_dense(:,:)
    type(wavevector) :: kv_g

    select case(setup%confinement)
    case(0)
      ! Bulk: build three 8x8 velocity matrices from dH/dk
      allocate(vel_dense(8, 8))

      ! vx = dH/dkx
      vel_dense = 0.0_dp
      kv_g%kx = 1.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandBulk(vel_dense, kv_g, cfg%params(1), cfg=cfg, g='g')
      call dense_to_csr_8(vel_dense, setup%vel(1))

      ! vy = dH/dky
      vel_dense = 0.0_dp
      kv_g%kx = 0.0_dp; kv_g%ky = 1.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandBulk(vel_dense, kv_g, cfg%params(1), cfg=cfg, g='g')
      call dense_to_csr_8(vel_dense, setup%vel(2))

      ! vz = dH/dkz
      vel_dense = 0.0_dp
      kv_g%kx = 0.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 1.0_dp
      call ZB8bandBulk(vel_dense, kv_g, cfg%params(1), cfg=cfg, g='g')
      call dense_to_csr_8(vel_dense, setup%vel(3))

      deallocate(vel_dense)

    case(1)
      ! QW: x/y/z from ZB8bandQW(g='g')
      allocate(vel_dense(setup%N, setup%N))

      ! vx = dH/dkx
      vel_dense = 0.0_dp
      kv_g%kx = 1.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandQW(vel_dense, kv_g, &
        setup%profile, setup%kpterms, cfg=cfg, g='g')
      call dense_to_csr_N(vel_dense, setup%vel(1))

      ! vy = dH/dky
      vel_dense = 0.0_dp
      kv_g%kx = 0.0_dp; kv_g%ky = 1.0_dp; kv_g%kz = 0.0_dp
      call ZB8bandQW(vel_dense, kv_g, &
        setup%profile, setup%kpterms, cfg=cfg, g='g')
      call dense_to_csr_N(vel_dense, setup%vel(2))

      ! vz = dH/dkz
      vel_dense = 0.0_dp
      kv_g%kx = 0.0_dp; kv_g%ky = 0.0_dp; kv_g%kz = 1.0_dp
      call ZB8bandQW(vel_dense, kv_g, &
        setup%profile, setup%kpterms, cfg=cfg, g='g')
      call dense_to_csr_N(vel_dense, setup%vel(3))

      deallocate(vel_dense)

    case default
      print *, 'Error: setup_build_velocity_matrices unsupported confinement=', &
        setup%confinement
      stop 1
    end select

    setup%vel_built = .true.

  end subroutine setup_build_velocity_matrices

  ! ==================================================================
  ! Free all allocatable components of simulation_setup.
  ! The finalizer delegates here.
  ! ==================================================================
  subroutine simulation_setup_free(setup)
    type(simulation_setup), intent(inout) :: setup
    integer :: i

    if (allocated(setup%profile)) deallocate(setup%profile)
    if (allocated(setup%kpterms)) deallocate(setup%kpterms)
    if (allocated(setup%HT)) deallocate(setup%HT)
    if (allocated(setup%work)) deallocate(setup%work)
    if (allocated(setup%rwork)) deallocate(setup%rwork)
    if (allocated(setup%iwork)) deallocate(setup%iwork)

    ! Velocity matrices
    if (setup%vel_built) then
      do i = 1, 3
        call csr_free(setup%vel(i))
      end do
      setup%vel_built = .false.
    end if

    ! Reset scalars
    setup%N = 0
    setup%confinement = -1
    setup%lwork = 0
    setup%il = 0
    setup%iuu = 0
    setup%has_strain = .false.
    setup%sc_was_run = .false.
    setup%vel_built = .false.

  end subroutine simulation_setup_free

  ! ==================================================================
  ! Finalizer: delegates to simulation_setup_free.
  ! ==================================================================
  subroutine simulation_setup_finalize(setup)
    type(simulation_setup), intent(inout) :: setup
    call simulation_setup_free(setup)
  end subroutine simulation_setup_finalize

  ! ==================================================================
  ! Helper: convert 8x8 dense to CSR
  ! ==================================================================
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

  ! ==================================================================
  ! Helper: convert NxN dense to CSR
  ! ==================================================================
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
