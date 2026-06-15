module strain_types

  ! ==============================================================================
  ! Strain result type (Issue #06, ADR 0005).
  !
  ! Extracted from strain_solver.f90 as the natural shared dependency between
  ! the Bir-Pikus / strain-table concern (strain_solver.f90) and the Navier-
  ! Cauchy PDE concern (strain_pde.f90). Both modules read/write the strain
  ! tensor through this type; keeping it in a leaf module avoids a circular
  ! `use` between the two strain concerns and lets the PDE solver be split
  ! out without touching defs.f90.
  !
  ! Nothing here is physics — just the strain-tensor container + its
  ! finalizer. The Bir-Pikus SSOT (compute_bp_scalar) and the strain table
  ! SSOT (get_strain_table) remain in strain_solver.f90.
  ! ==============================================================================

  use definitions, only: dp
  implicit none

  private

  public :: strain_result
  public :: strain_result_free

  ! ------------------------------------------------------------------
  ! Strain tensor at each grid point.
  ! For QW: arrays have length grid%ny.
  ! For wire: arrays have length grid%nx * grid%ny.
  ! ------------------------------------------------------------------
  type :: strain_result
    real(kind=dp), allocatable :: eps_xx(:)   ! (Ngrid) strain along wire axis / growth
    real(kind=dp), allocatable :: eps_yy(:)   ! (Ngrid) in-plane strain
    real(kind=dp), allocatable :: eps_zz(:)   ! (Ngrid) in-plane / growth strain
    real(kind=dp), allocatable :: eps_yz(:)   ! (Ngrid) shear strain
  contains
    final :: strain_result_finalize
  end type strain_result

contains

  ! ==================================================================
  ! Finalizer: automatically called when a strain_result goes out of
  ! scope. Delegates to strain_result_free so existing manual frees
  ! remain valid (double-free safe via the implicit allocated() checks).
  ! ==================================================================
  subroutine strain_result_finalize(sr)
    type(strain_result), intent(inout) :: sr
    call strain_result_free(sr)
  end subroutine strain_result_finalize

  ! ==================================================================
  ! Deallocate strain_result arrays (manual free; idempotent).
  ! ==================================================================
  subroutine strain_result_free(s)
    type(strain_result), intent(inout) :: s

    if (allocated(s%eps_xx)) deallocate(s%eps_xx)
    if (allocated(s%eps_yy)) deallocate(s%eps_yy)
    if (allocated(s%eps_zz)) deallocate(s%eps_zz)
    if (allocated(s%eps_yz)) deallocate(s%eps_yz)
  end subroutine strain_result_free

end module strain_types
