program fd_sign_compare

  !> Compare order-2 tridiagonal path vs higher-order FD matrix path.
  !> Verifies kpterms construction produces consistent results.
  !>
  !> For FDorder=2, both paths must produce identical kpterms.
  !> Any mismatch indicates a sign/factor bug.

  use definitions
  use parameters
  use hamiltonianConstructor
  use finitedifferences

  implicit none

  integer, parameter :: NL = 41
  real(kind=dp) :: delta, totalSize
  real(kind=dp) :: z(NL)
  real(kind=dp), allocatable :: profile(:,:), kpterms_o2(:,:,:), kpterms_ho(:,:,:)

  ! startPos/endPos are real(dp) in the interface
  real(kind=dp) :: startPos(2), endPos(2)
  integer :: intStartPos(2), intEndPos(2)
  character(len=255) :: material(2)
  type(paramStruct) :: params(2)

  integer :: i, term, ii, jj
  real(kind=dp) :: maxdiff, val_o2, val_ho
  logical :: all_pass

  ! Setup grid: z from -10 to +10 nm
  totalSize = 20.0_dp
  delta = totalSize / real(NL - 1, kind=dp)
  z = [( -10.0_dp + (i-1)*delta, i=1,NL )]

  ! 2-layer QW: AlGaAs barriers + GaAs well
  startPos(1) = -10.0_dp
  endPos(1) = 10.0_dp
  startPos(2) = -5.0_dp
  endPos(2) = 5.0_dp

  intStartPos(1) = 1
  intEndPos(1) = NL
  intStartPos(2) = int(abs(startPos(1) - startPos(2)) / delta) + 1
  intEndPos(2) = intStartPos(2) + int(abs(endPos(2) - startPos(2)) / delta)

  material(1) = "Al30Ga70As"
  material(2) = "GaAs"

  call paramDatabase(material, 2, params)

  ! =================== ORDER-2 PATH ===================
  allocate(kpterms_o2(NL, NL, 10))
  kpterms_o2 = 0.0_dp
  allocate(profile(NL, 3))
  profile = 0.0_dp
  call confinementInitialization(z, intStartPos, intEndPos, material, 2, params, &
    & 'z', profile, kpterms_o2, FDorder=2)

  ! =================== HIGHER-ORDER PATH (FDorder=4) ===================
  allocate(kpterms_ho(NL, NL, 10))
  kpterms_ho = 0.0_dp
  deallocate(profile)
  allocate(profile(NL, 3))
  profile = 0.0_dp
  call confinementInitialization(z, intStartPos, intEndPos, material, 2, params, &
    & 'z', profile, kpterms_ho, FDorder=4)

  ! =================== COMPARISON ===================
  print *, '=== kpterms element-by-element comparison ==='
  print *, 'Grid: N=', NL, ' delta=', delta
  print *, ''

  ! Dump tridiagonal structure for term 5 (A*kz^2 = 2nd derivative)
  print *, '--- Term 5 (A*kz^2) interior [row 10, col 8..12] ---'
  print *, 'Order-2: ', (kpterms_o2(10,jj,5), jj=8,12)
  print *, 'Order-4: ', (kpterms_ho(10,jj,5), jj=8,12)
  print *, ''

  ! Dump tridiagonal structure for term 6 (P*kz = 1st derivative)
  print *, '--- Term 6 (P*kz) interior [row 10, col 8..12] ---'
  print *, 'Order-2: ', (kpterms_o2(10,jj,6), jj=8,12)
  print *, 'Order-4: ', (kpterms_ho(10,jj,6), jj=8,12)
  print *, ''

  ! Element-by-element comparison for interior points [5:NL-4]
  ! Skip boundary regions where stencils differ
  print *, '=== Interior point comparison (rows 5:37) ==='
  all_pass = .true.
  do term = 1, 10
    maxdiff = 0.0_dp
    do ii = 5, NL-4
      do jj = max(1, ii-4), min(NL, ii+4)
        val_o2 = kpterms_o2(jj, ii, term)
        val_ho = kpterms_ho(jj, ii, term)
        if (abs(val_o2) > 1.0e-15 .or. abs(val_ho) > 1.0e-15) then
          maxdiff = max(maxdiff, abs(val_o2 - val_ho))
        end if
      end do
    end do

    if (maxdiff > 1.0e-10) then
      print *, '  term ', term, ': DIFFER  max_diff=', maxdiff
      all_pass = .false.
      ! Show first few differing elements
      do ii = 10, 15
        do jj = ii-2, ii+2
          val_o2 = kpterms_o2(jj, ii, term)
          val_ho = kpterms_ho(jj, ii, term)
          if (abs(val_o2 - val_ho) > 1.0e-12 .and. abs(val_o2) > 1.0e-15) then
            if (abs(val_o2) > 1.0e-15) then
              print *, '    [', jj, ',', ii, '] o2=', val_o2, ' o4=', val_ho, &
                & ' ratio=', val_ho/val_o2
            else
              print *, '    [', jj, ',', ii, '] o2=', val_o2, ' o4=', val_ho
            end if
          end if
        end do
      end do
    else
      print *, '  term ', term, ': MATCH  max_diff=', maxdiff
    end if
  end do

  print *, ''
  if (all_pass) then
    print *, 'ALL INTERIOR POINTS MATCH (order-2 vs order-4 agree)'
  else
    print *, 'DIFFERENCES FOUND - check sign/factor convention'
    print *, ''
    print *, 'NOTE: order-4 uses different stencil, so values will differ'
    print *, '      from order-2. Check that the ratio is consistent and'
    print *, '      that signs match for corresponding matrix elements.'
  end if

  ! Also dump diagonal values for visual comparison
  print *, ''
  print *, '=== Diagonal values comparison [points 10:15] ==='
  do term = 5, 9
    print *, 'Term ', term, ':'
    print *, '  O2 diag: ', (kpterms_o2(ii,ii,term), ii=10,15)
    print *, '  O4 diag: ', (kpterms_ho(ii,ii,term), ii=10,15)
  end do

  deallocate(kpterms_o2, kpterms_ho, profile)

end program fd_sign_compare
