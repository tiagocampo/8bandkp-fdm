!=====================================================================
! Coarray QW k-point sweep prototype
!
! Standalone program demonstrating the coarray parallelism pattern for
! distributing k-points across images in a QW band structure sweep.
!
! Compile (single-image, validates correctness):
!   gfortran -fcoarray=single -o test_coarray_qw_sweep \
!           tests/integration/test_coarray_qw_sweep.f90
!
! Compile (multi-image via OpenCoarrays):
!   caf -o test_coarray_qw_sweep tests/integration/test_coarray_qw_sweep.f90
!   cafrun -n 4 ./test_coarray_qw_sweep
!
! This prototype demonstrates the distribution pattern only.  Full
! integration (calling ZB8bandQW + zheevx) requires exposing smallk
! via simulation_config, which is out of scope for Phase 4 exploration.
!=====================================================================
program test_coarray_qw_sweep
  use, intrinsic :: iso_fortran_env, only: real64, output_unit
  implicit none

  integer, parameter :: dp = real64

  ! --- Simulation parameters (hardcoded; normally from input.cfg) ---
  integer, parameter :: nk_total    = 100     ! total k-points in sweep
  integer, parameter :: n_bands     = 8       ! bands per k-point (8-band k.p)
  integer, parameter :: n_fd        = 40      ! FD grid points (FDstep)
  integer, parameter :: ham_size    = n_bands * n_fd  ! Hamiltonian dimension

  real(kind=dp), parameter :: k_max = 0.05_dp  ! max k-vector (1/A)
  character(len=4), parameter :: sweep_dir = 'kz  '   ! sweep direction

  ! --- Coarray distribution ---
  ! Each image owns a contiguous block of k-points.
  ! n_local is the number of k-points this image handles.
  ! k_lo/k_hi are the 1-based indices into the global k-array.
  integer :: n_images, img, n_local, k_lo, k_hi, remainder, i, k

  ! --- k-vector array (local portion) ---
  real(kind=dp), allocatable :: k_values(:)

  ! --- Eigenvalue result: each image computes its block independently ---
  ! coarray dimension [*] means one copy per image.
  ! Image p owns eig_local(:, k_lo:k_hi).
  real(kind=dp), allocatable :: eig_local(:,:)[:]   ! (n_bands, nk_total)

  ! --- Collected results on image 1 ---
  real(kind=dp), allocatable :: eig_global(:,:)      ! (n_bands, nk_total)

  ! ------------------------------------------------------------------
  ! Distribution: block-decompose k-points across images
  ! ------------------------------------------------------------------
  img      = this_image()
  n_images = num_images()

  ! Block size: nk_total / n_images, with remainder distributed to
  ! the first 'remainder' images (standard block-cyclic distribution).
  n_local  = nk_total / n_images
  remainder = mod(nk_total, n_images)

  if (img <= remainder) then
    n_local = n_local + 1
    k_lo   = (img - 1) * n_local + 1
  else
    k_lo   = remainder * (n_local + 1) + (img - remainder - 1) * n_local + 1
  end if
  k_hi = k_lo + n_local - 1

  ! ------------------------------------------------------------------
  ! Allocate coarray eigenvalue array (same size on every image so
  ! image 1 can pull any slice via coarray subscript)
  ! ------------------------------------------------------------------
  allocate(eig_local(n_bands, nk_total)[*])
  eig_local = 0.0_dp

  ! Build local k-vector values
  allocate(k_values(k_lo:k_hi))
  do k = k_lo, k_hi
    ! Linear spacing from 0 to k_max, matching main.f90 convention
    k_values(k) = real(k - 1, dp) * k_max / real(nk_total - 1, dp)
  end do

  ! ------------------------------------------------------------------
  ! Simulate k-point sweep: each image "computes" eigenvalues for its
  ! block.  In the real implementation this would be:
  !   call ZB8bandQW(H_local, smallk(k), profile, kpterms, cfg=cfg)
  !   call zheevx(...)
  ! Here we fill with a deterministic pattern so we can verify collection.
  ! ------------------------------------------------------------------
  do k = k_lo, k_hi
    do i = 1, n_bands
      ! Deterministic eigenvalue: band-dependent offset + k-dependent slope
      ! E_n(k) = n * 0.1 + k * 10.0  (arbitrary units, pattern only)
      eig_local(i, k) = real(i, dp) * 0.1_dp + k_values(k) * 10.0_dp
    end do
  end do

  ! ------------------------------------------------------------------
  ! Synchronize: ensure all images have computed their blocks
  ! ------------------------------------------------------------------
  sync all

  ! ------------------------------------------------------------------
  ! Image 1: collect results from all images
  ! ------------------------------------------------------------------
  if (img == 1) then
    write(output_unit, '(A)')        '=== Coarray QW k-point sweep prototype ==='
    write(output_unit, '(A,I0)')     'Total k-points:    ', nk_total
    write(output_unit, '(A,I0)')     'Hamiltonian size:  ', ham_size
    write(output_unit, '(A,I0)')     'Sweep direction:   ' // trim(sweep_dir)
    write(output_unit, '(A,I0)')     'Number of images:  ', n_images
    write(output_unit, '(A)')        ''

    ! Report per-image distribution
    do i = 1, n_images
      block
        integer :: nl, klo, khi, rem
        nl  = nk_total / n_images
        rem = mod(nk_total, n_images)
        if (i <= rem) then
          nl  = nl + 1
          klo = (i - 1) * nl + 1
        else
          klo = rem * (nl + 1) + (i - rem - 1) * nl + 1
        end if
        khi = klo + nl - 1
        write(output_unit, '(A,I0,A,I0,A,I0)') &
          '  Image ', i, ': k-points ', klo, ' -- ', khi
      end block
    end do
    write(output_unit, '(A)') ''

    ! Collect: copy from each image's coarray into local buffer.
    ! coarray subscript eig_local(:,:)[p] reads image p's data.
    allocate(eig_global(n_bands, nk_total))
    do i = 1, n_images
      ! Pull the block owned by image i.
      ! We recompute k_lo/k_hi for image i (same formula as above).
      block
        integer :: nl, klo, khi, rem
        nl  = nk_total / n_images
        rem = mod(nk_total, n_images)
        if (i <= rem) then
          nl  = nl + 1
          klo = (i - 1) * nl + 1
        else
          klo = rem * (nl + 1) + (i - rem - 1) * nl + 1
        end if
        khi = klo + nl - 1
        eig_global(:, klo:khi) = eig_local(:, klo:khi)[i]
      end block
    end do

    ! Verify: check first and last k-points
    write(output_unit, '(A)') 'Eigenvalues at k=0 (band 1-8):'
    write(output_unit, '(8F8.3)') eig_global(:, 1)
    write(output_unit, '(A)') 'Eigenvalues at k=k_max (band 1-8):'
    write(output_unit, '(8F8.3)') eig_global(:, nk_total)

    ! Consistency check: eig_global(:,1) should equal eig_local(:,1)
    ! since image 1 owns k=1 when remainder >= 0 (always true).
    if (all(abs(eig_global(:, 1) - eig_local(:, 1)) < 1.0e-12_dp)) then
      write(output_unit, '(A)') ''
      write(output_unit, '(A)') 'Collection consistency: OK'
    else
      write(output_unit, '(A)') 'Collection consistency: FAILED'
      error stop 1
    end if

    ! Expected values check
    ! E_1(k=0) = 1*0.1 + 0*10 = 0.1
    ! E_8(k=0) = 8*0.1 + 0*10 = 0.8
    ! E_1(k_max) = 0.1 + 0.05*10 = 0.6
    ! E_8(k_max) = 0.8 + 0.05*10 = 1.3
    if (abs(eig_global(1, 1) - 0.1_dp) < 1.0e-12_dp .and. &
        abs(eig_global(8, 1) - 0.8_dp) < 1.0e-12_dp .and. &
        abs(eig_global(1, nk_total) - 0.6_dp) < 1.0e-12_dp .and. &
        abs(eig_global(8, nk_total) - 1.3_dp) < 1.0e-12_dp) then
      write(output_unit, '(A)') 'Eigenvalue values:    OK'
    else
      write(output_unit, '(A)') 'Eigenvalue values:    FAILED'
      error stop 1
    end if

    write(output_unit, '(A)') ''
    write(output_unit, '(A)') 'PASS: coarray QW k-point sweep prototype'
  end if

  ! Cleanup
  deallocate(eig_local)
  if (allocated(k_values)) deallocate(k_values)
  if (allocated(eig_global)) deallocate(eig_global)

end program test_coarray_qw_sweep
