module wire_setup_mod

  ! ==============================================================================
  ! wire_setup: strain-aware wire initialization + cleanup type.
  !
  ! Owns the 2D-wire confinement state (profile_2d, sparse k.p terms kpterms_2d,
  ! wire workspace, COO cache) and runs the SAME strain setup as the canonical
  ! QW/wire path in simulation_setup_init (case 'wire'): confinementInitialization_2d
  ! followed, when cfg%strain%enabled, by compute_strain + compute_bir_pikus_blocks.
  !
  ! Why this exists (Issue #04, ADRs 0001 + 0005):
  !   The topology wire subroutines (run_bdg_wire, run_qshe_wire, eval_wire_bdg_gap_app
  !   in main_topology.f90) and the spectral-function path (compute_spectral_function_wire
  !   in green_functions.f90) each copy-pasted confinementInitialization_2d but
  !   SKIPPED compute_strain + compute_bir_pikus_blocks. The strain insertion in
  !   ZB8bandGeneralized is gated on `allocated(cfg%strain_blocks%delta_Ec)`, so
  !   the strain field was silently dropped on every strained wire topology/BdG/
  !   spectral run. Routing those paths through this type fixes the omission by
  !   construction: init always runs the strain step when cfg%strain%enabled.
  !
  ! Design:
  !   - Concrete type (ADR 0001: no polymorphic builder hierarchy).
  !   - Single init entry point:
  !       wire_setup_init             -- strain-applying (the canonical path).
  !   - Idempotent free via was_freed flag (CLAUDE.md pattern).
  !   - finalizer delegates to wire_setup_free.
  !
  ! NOTE: BHZ wire (build_bhz_wire_hamiltonian) is a separate 4-band model that
  !   does NOT use profile_2d/kpterms_2d; 8-band k.p strain does not apply to it.
  !   Do not route BHZ through this type.
  ! ==============================================================================

  use definitions, only: dp, grid_ngrid, simulation_config
  use sparse_matrices, only: csr_matrix, csr_free
  use geometry, only: init_wire_from_config
  use confinement_init, only: confinementInitialization_2d
  use hamiltonian_wire, only: wire_workspace, wire_workspace_free, &
    & wire_coo_cache, wire_coo_cache_free
  use strain_solver, only: compute_strain, compute_bir_pikus_blocks, &
    & strain_result, strain_result_free

  implicit none
  private

  public :: wire_setup
  public :: wire_setup_init, wire_setup_free

  type :: wire_setup
    ! --- 2D-wire confinement state (mirrors simulation_setup wire fields) ---
    real(kind=dp), allocatable    :: profile_2d(:,:)
    type(csr_matrix), allocatable :: kpterms_2d(:)
    type(wire_coo_cache), allocatable :: coo_cache
    type(wire_workspace), allocatable :: ws
    logical :: has_strain = .false.
    logical :: was_freed = .false.
  contains
    final :: wire_setup_finalize
  end type wire_setup

contains

  ! ==============================================================================
  ! wire_setup_init: the strain-applying canonical wire initialization.
  !
  ! Mirrors simulation_setup_init case('wire') (lines ~296-315): runs
  ! confinementInitialization_2d, then (if cfg%strain%enabled) compute_strain
  ! + compute_bir_pikus_blocks which populates cfg%strain_blocks. That
  ! population is what unblocks the strain insertion in ZB8bandGeneralized.
  !
  ! Mutates cfg: populates cfg%strain_blocks when strain is enabled, and
  ! ensures cfg%grid%x is allocated (init_wire_from_config). This matches the
  ! canonical path exactly.
  ! ==============================================================================
  subroutine wire_setup_init(self, cfg)
    type(wire_setup), intent(inout) :: self
    type(simulation_config), intent(inout) :: cfg

    type(strain_result), allocatable :: strain_local
    real(kind=dp) :: a0_ref

    if (self%was_freed) then
      print *, 'Error: wire_setup_init called on a freed wire_setup'
      error stop 'wire_setup_init: instance already freed'
    end if

    ! --- Ensure the wire grid is built (mirrors canonical case('wire')) ---
    if (.not. allocated(cfg%grid%x)) call init_wire_from_config(cfg)

    ! --- 2D confinement init (kpterms FD operator matrices for the wire) ---
    call confinementInitialization_2d(cfg%grid, cfg%params, cfg%wire%regions, &
      & self%profile_2d, self%kpterms_2d, cfg%FDorder)

    ! --- Strain: the step the copy-pasted paths silently skipped (#04) ---
    ! Identical to simulation_setup_init case('wire') strain block.
    if (cfg%strain%enabled) then
      allocate(strain_local)
      a0_ref = cfg%params(1)%a0
      call compute_strain(cfg%grid, cfg%params, cfg%grid%material_id, &
        cfg%strain, a0_ref, strain_local)
      call compute_bir_pikus_blocks(strain_local, cfg%params, &
        cfg%grid%material_id, cfg%grid, cfg%strain_blocks)
      call strain_result_free(strain_local)
      deallocate(strain_local)
      self%has_strain = .true.
    else
      self%has_strain = .false.
    end if

    ! --- Allocate workspace handles (callers populate/use as needed) ---
    allocate(self%coo_cache)
    allocate(self%ws)
  end subroutine wire_setup_init

  ! ==============================================================================
  ! wire_setup_free: idempotent cleanup of all wire state.
  ! ==============================================================================
  subroutine wire_setup_free(self)
    type(wire_setup), intent(inout) :: self
    integer :: i

    if (self%was_freed) return
    self%was_freed = .true.

    if (allocated(self%profile_2d)) deallocate(self%profile_2d)
    if (allocated(self%kpterms_2d)) then
      do i = 1, size(self%kpterms_2d)
        call csr_free(self%kpterms_2d(i))
      end do
      deallocate(self%kpterms_2d)
    end if
    if (allocated(self%coo_cache)) then
      call wire_coo_cache_free(self%coo_cache)
      deallocate(self%coo_cache)
    end if
    if (allocated(self%ws)) then
      call wire_workspace_free(self%ws)
      deallocate(self%ws)
    end if
    self%has_strain = .false.
  end subroutine wire_setup_free

  subroutine wire_setup_finalize(self)
    type(wire_setup), intent(inout) :: self
    call wire_setup_free(self)
  end subroutine wire_setup_finalize

end module wire_setup_mod
