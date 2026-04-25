program opticalProperties
  use definitions
  use input_parser
  implicit none
  type(simulation_config) :: cfg
  real(kind=dp), allocatable, dimension(:,:) :: profile
  real(kind=dp), allocatable, dimension(:,:,:) :: kpterms

  call read_and_setup(cfg, profile, kpterms)
  print '(a)', 'opticalProperties: input parsed successfully'
  print '(a,L)', '  optics enabled = ', cfg%optics%enabled

end program opticalProperties
