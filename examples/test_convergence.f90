! Example file demonstrating convergence testing and step size optimization
! This example shows how to use the new gfactorCalculationNumericalConverged subroutine

program test_convergence

  use gfactorFunctions
  use definitions
  implicit none

  ! Configuration structures
  type(NumericalConfig) :: num_config
  type(ConvergenceConfig) :: conv_config
  type(ConvergenceResults) :: conv_results

  ! Variables for g-factor calculation
  complex(kind=dp), allocatable :: tensor(:,:,:)
  complex(kind=dp), allocatable :: cb_state(:,:), vb_state(:,:)
  real(kind=dp), allocatable :: cb_value(:), vb_value(:)
  type(paramStruct) :: params(1)

  ! Test parameters
  integer :: nlayers = 1
  real(kind=dp) :: startz = 0.0_dp, endz = 10.0_dp, dz = 0.1_dp
  integer :: whichBand = 0, bandIdx = 1, numcb = 2, numvb = 6

  print *, 'Testing numerical g-factor convergence framework'

  ! Set up numerical configuration
  num_config%method = 'numerical'
  num_config%tolerance = 1e-12_dp
  num_config%perturbation_step = 1e-3_dp  ! Start with relatively large step
  num_config%use_sparse_solver = .false.
  num_config%validate_with_analytical = .true.
  num_config%verbose_output = .true.

  ! Set up convergence configuration
  conv_config%energy_tolerance = 1e-6_dp          ! Energy difference < 1e-6 eV
  conv_config%gfactor_tolerance = 1e-3_dp         ! g-factor change < 0.1%
  conv_config%min_step_size = 1e-12_dp           ! Minimum perturbation step
  conv_config%step_reduction_factor = 0.1_dp     ! Reduce step by factor of 10
  conv_config%max_iterations = 10                ! Maximum convergence iterations
  conv_config%enable_adaptive_stepping = .true.  ! Enable adaptive step sizing
  conv_config%verbose_convergence = .true.       ! Detailed convergence output

  ! Allocate arrays (example sizes)
  allocate(tensor(2,2,3))
  allocate(cb_state(8,2))
  allocate(vb_state(8,6))
  allocate(cb_value(2))
  allocate(vb_value(6))

  ! Initialize with dummy data (in real usage, these would come from actual calculations)
  cb_state = (0.0_dp, 0.0_dp)
  vb_state = (0.0_dp, 0.0_dp)
  cb_value = 0.0_dp
  vb_value = 0.0_dp

  ! Set up basic material parameters (example for InAs)
  params(1)%m0 = 0.023_dp
  params(1)%Ep = 21.5_dp
  params(1)%Eg = 0.417_dp
  params(1)%Es0 = 19.7_dp
  params(1)%Delta_so = 0.39_dp
  params(1)%gamma1 = 20.0_dp
  params(1)%gamma2 = 8.5_dp
  params(1)%gamma3 = 9.2_dp
  params(1)%C11 = 8.329e6_dp
  params(1)%C12 = 4.526e6_dp
  params(1)%C44 = 3.959e6_dp
  params(1)%a0 = 6.0583_dp
  params(1)%ac = -5.08_dp
  params(1)%av = 1.00_dp
  params(1)%b = -1.8_dp
  params(1)%d = -3.6_dp
  params(1)%kappa = 1.0_dp
  params(1)%g_e = -14.6_dp

  print *, 'Configuration set up successfully'
  print *, 'Initial perturbation step:', num_config%perturbation_step
  print *, 'Energy tolerance:', conv_config%energy_tolerance, 'eV'
  print *, 'G-factor tolerance:', conv_config%gfactor_tolerance
  print *, 'Maximum iterations:', conv_config%max_iterations

  ! Test convergence framework
  print *, 'Running converged g-factor calculation...'

  call gfactorCalculationNumericalConverged(tensor, whichBand, bandIdx, numcb, numvb, &
    & cb_state, vb_state, cb_value, vb_value, nlayers, params, startz, endz, &
    & dz, num_config, conv_config, conv_results)

  ! Output results
  print *, ''
  print *, '=== Convergence Results ==='
  print *, 'Converged:', conv_results%converged
  print *, 'Iterations used:', conv_results%iterations_used
  print *, 'Final step size:', conv_results%final_step_size
  print *, 'Final energy difference:', conv_results%final_energy_diff, 'eV'
  print *, 'Final g-factor difference:', conv_results%final_gfactor_diff
  print *, 'Convergence status:', trim(conv_results%convergence_status)
  print *, 'Numerical stability achieved:', conv_results%numerical_stability_achieved

  if (conv_results%converged) then
    print *, ''
    print *, 'Final g-factor tensor components:'
    print *, 'g_xx:', real(tensor(1,1,1))
    print *, 'g_yy:', real(tensor(2,2,2))
    print *, 'g_zz:', real(tensor(1,1,3))
    print *, 'Initial g-factor:', conv_results%initial_gfactor
    print *, 'Final g-factor:', conv_results%final_gfactor
  else
    print *, ''
    print *, 'Convergence failed. See status message above.'
  end if

  ! Clean up
  deallocate(tensor)
  deallocate(cb_state)
  deallocate(vb_state)
  deallocate(cb_value)
  deallocate(vb_value)

  print *, ''
  print *, 'Convergence test completed.'

end program test_convergence