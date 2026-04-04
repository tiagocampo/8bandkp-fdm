!===============================================================================
! GPU diagonalization test: hipSOLVER Zheevd via iso_c_binding
! Proof-of-concept for AMD GPU Hermitian eigenvalue decomposition
!
! Build:
!   gfortran -o gpu_diag_test tests/gpu_diag_test.f90 \
!     -I/opt/rocm/include \
!     -L/opt/rocm/lib -lhipsolver -lrocblas -lhipblas \
!     -Wl,-rpath,/opt/rocm/lib
!===============================================================================

program gpu_diag_test
  use, intrinsic :: iso_c_binding
  implicit none

  ! --- hipSOLVER C interface via iso_c_binding ---
  type, bind(c) :: hipsolverHandle
    integer(c_intptr_t) :: ptr = 0
  end type

  ! hipSOLVER enums (matching hipsolver-types.h)
  integer(c_int), parameter :: HIPSOLVER_EIG_MODE_VECTOR = 0  ! jobz='V'
  integer(c_int), parameter :: HIPSOLVER_FILL_MODE_UPPER  = 0 ! uplo='U'
  integer(c_int), parameter :: HIPSOLVER_STATUS_SUCCESS = 0

  ! hipStream_t is void*
  integer(c_intptr_t), parameter :: null_stream = 0

  interface
    function hipsolverCreate_c(handle) bind(c, name='hipsolverCreate')
      import :: c_int, c_ptr
      integer(c_int) :: hipsolverCreate_c
      type(c_ptr), intent(out) :: handle
    end function

    function hipsolverDestroy_c(handle) bind(c, name='hipsolverDestroy')
      import :: c_int, c_ptr
      integer(c_int) :: hipsolverDestroy_c
      type(c_ptr), value :: handle
    end function

    function hipsolverZheevd_bufferSize_c(handle, jobz, uplo, n, A, lda, &
                                          D, lwork) bind(c, name='hipsolverZheevd_bufferSize')
      import :: c_int, c_double, c_ptr
      integer(c_int) :: hipsolverZheevd_bufferSize_c
      type(c_ptr), value :: handle
      integer(c_int), value :: jobz, uplo, n, lda
      type(c_ptr), value :: A
      type(c_ptr), value :: D  ! double* on device
      type(c_ptr), value :: lwork  ! int* on host
    end function

    function hipsolverZheevd_c(handle, jobz, uplo, n, A, lda, &
                               D, work, lwork, devInfo) bind(c, name='hipsolverZheevd')
      import :: c_int, c_double, c_ptr
      integer(c_int) :: hipsolverZheevd_c
      type(c_ptr), value :: handle
      integer(c_int), value :: jobz, uplo, n, lda
      type(c_ptr), value :: A        ! hipDoubleComplex* on device
      type(c_ptr), value :: D        ! double* on device
      type(c_ptr), value :: work     ! hipDoubleComplex* on device
      integer(c_int), value :: lwork
      type(c_ptr), value :: devInfo  ! int* on device
    end function

    ! hipMalloc / hipFree
    function hipMalloc_c(ptr, size) bind(c, name='hipMalloc')
      import :: c_int, c_ptr, c_size_t
      integer(c_int) :: hipMalloc_c
      type(c_ptr), intent(out) :: ptr
      integer(c_size_t), value :: size
    end function

    function hipFree_c(ptr) bind(c, name='hipFree')
      import :: c_int, c_ptr
      integer(c_int) :: hipFree_c
      type(c_ptr), value :: ptr
    end function

    ! hipMemcpy (kind: hipMemcpyHostToDevice=1, hipMemcpyDeviceToHost=2)
    function hipMemcpy_c(dst, src, size, kind) bind(c, name='hipMemcpy')
      import :: c_int, c_ptr, c_size_t
      integer(c_int) :: hipMemcpy_c
      type(c_ptr), value :: dst, src
      integer(c_size_t), value :: size
      integer(c_int), value :: kind
    end function
  end interface

  ! --- Test matrix ---
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: N = 8  ! same as bulk k.p Hamiltonian size
  complex(dp), target :: H(N, N)
  real(dp), target :: eig(N)
  complex(dp), target :: H_host(N, N)  ! copy for GPU result
  integer :: i, j

  ! --- GPU variables ---
  type(c_ptr) :: handle = c_null_ptr
  type(c_ptr) :: d_A = c_null_ptr, d_D = c_null_ptr, d_work = c_null_ptr, d_info = c_null_ptr
  integer(c_int) :: lwork, status
  integer(c_int), target :: lwork_target
  integer, target :: host_info

  ! Build a simple Hermitian test matrix (tridiagonal)
  H = cmplx(0, 0, kind=dp)
  do i = 1, N
    H(i, i) = cmplx(real(i, dp), 0, kind=dp)
    if (i < N) then
      H(i, i+1) = cmplx(0.5_dp, 0.3_dp, kind=dp)
      H(i+1, i) = cmplx(0.5_dp, -0.3_dp, kind=dp)
    end if
  end do

  print *, '=== GPU Eigenvalue Test (hipSOLVER Zheevd) ==='
  print *, 'Matrix size: ', N, 'x', N
  print *, ''
  print *, 'Input Hermitian matrix:'
  do i = 1, N
    do j = 1, N
      if (abs(H(i,j)) > 1e-10_dp) then
        write(*, '(A,I2,A,I2,A,2F8.4)') '  H(', i, ',', j, ') = ', H(i,j)
      end if
    end do
  end do

  ! --- CPU reference: solve with LAPACK zheev ---
  call cpu_reference(H, eig, N)

  ! --- GPU path ---
  print *, ''
  print *, '--- GPU path ---'

  ! Create hipSOLVER handle
  status = hipsolverCreate_c(handle)
  if (status /= HIPSOLVER_STATUS_SUCCESS) then
    print *, 'ERROR: hipsolverCreate failed, status =', status
    print *, 'GPU may not be available or supported on this system.'
    stop 1
  end if
  print *, 'hipSOLVER handle created successfully.'

  ! Allocate device memory
  status = hipMalloc_c(d_A, int(N*N*16, c_size_t))  ! complex(dp) = 16 bytes
  if (status /= 0) then
    print *, 'ERROR: hipMalloc for A failed, status =', status
    stop 1
  end if
  status = hipMalloc_c(d_D, int(N*8, c_size_t))     ! real(dp) = 8 bytes
  status = hipMalloc_c(d_info, int(4, c_size_t))     ! int = 4 bytes

  ! Copy H to device
  status = hipMemcpy_c(d_A, c_loc(H), int(N*N*16, c_size_t), 1)
  print *, 'Matrix copied to GPU.'

  ! Query workspace size
  status = hipsolverZheevd_bufferSize_c(handle, &
    HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, &
    N, d_A, N, d_D, c_loc(lwork_target))
  lwork = lwork_target
  if (status /= HIPSOLVER_STATUS_SUCCESS) then
    print *, 'ERROR: bufferSize query failed, status =', status
    stop 1
  end if
  print *, 'Workspace size: ', lwork, ' bytes (', lwork/16, ' complex*16 elements)'

  ! Allocate workspace
  status = hipMalloc_c(d_work, int(lwork, c_size_t))

  ! Run eigenvalue decomposition
  print *, 'Running hipsolverZheevd...'
  status = hipsolverZheevd_c(handle, &
    HIPSOLVER_EIG_MODE_VECTOR, HIPSOLVER_FILL_MODE_UPPER, &
    N, d_A, N, d_D, d_work, lwork, d_info)
  if (status /= HIPSOLVER_STATUS_SUCCESS) then
    print *, 'ERROR: hipsolverZheevd failed, status =', status
    stop 1
  end if

  ! Copy results back
  status = hipMemcpy_c(c_loc(eig), d_D, int(N*8, c_size_t), 2)  ! device->host
  status = hipMemcpy_c(c_loc(H_host), d_A, int(N*N*16, c_size_t), 2)
  status = hipMemcpy_c(c_loc(host_info), d_info, int(4, c_size_t), 2)

  print *, 'devInfo:', host_info
  print *, ''
  print *, 'GPU eigenvalues:'
  do i = 1, N
    write(*, '(A,I2,A,F12.6)') '  eig(', i, ') = ', eig(i)
  end do

  ! Cleanup
  status = hipFree_c(d_A)
  status = hipFree_c(d_D)
  status = hipFree_c(d_work)
  status = hipFree_c(d_info)
  status = hipsolverDestroy_c(handle)
  print *, ''
  print *, 'GPU resources freed. Test complete!'

contains

  subroutine cpu_reference(H_in, eig_ref, n)
    integer, intent(in) :: n
    complex(dp), intent(inout) :: H_in(n, n)
    real(dp), intent(out) :: eig_ref(n)
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer :: info, lwork_local
    external :: zheev

    ! Workspace query
    allocate(work(1), rwork(3*n-2))
    call zheev('N', 'U', n, H_in, n, eig_ref, work, -1, rwork, info)
    lwork_local = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork_local))

    ! Solve
    call zheev('N', 'U', n, H_in, n, eig_ref, work, lwork_local, rwork, info)

    print *, ''
    print *, '--- CPU reference (LAPACK zheev) ---'
    print *, 'CPU eigenvalues:'
    do i = 1, n
      write(*, '(A,I2,A,F12.6)') '  eig(', i, ') = ', eig_ref(i)
    end do

    deallocate(work, rwork)
  end subroutine

end program
