module bdg_observables

  ! ==============================================================================
  ! BdG per-point observables — the foundational seam extracted from app glue.
  !
  ! All three per-point call sites in main_topology (run_bdg_wire, run_bdg_qw,
  ! eval_wire_bdg_gap) used to inline:
  !   - minigap: 2 * minval(abs(eigvals_bdg))
  !   - near-zero threshold: 0.001 * delta_0
  !   - invariant flag: count(|E| < threshold) >= 2
  !
  ! This module folds those three steps into one pure-function call so the
  ! build-and-solve stays in main_topology (per ADR 0003) while the per-point
  ! physics decision lives in one place. Downstream slices (Pfaffian wrapper,
  ! polarization, LDOS) consume the same contract.
  !
  ! Pure only — no I/O, no allocations, no state.
  ! ==============================================================================

  use definitions, only: dp
  use sparse_matrices, only: csr_matrix
  use pfaffian, only: complex_pfaffian, kitaev_majorana_number
  use topological_analysis, only: wire_pfaffian_witness_sweep

  implicit none

  private

  public :: bdg_eval_params_t
  public :: bdg_eval_result_t
  public :: eval_bdg_point
  public :: q_zero_tol
  public :: bdg_eval_params_with_delta
  public :: eval_bdg_pfaffian_witness_csr

  ! Module-level defaults (SSOT for the near-zero literals). Extracted from
  ! inline 0.001_dp / 1.0e-10_dp at the 5 call sites in main_topology so the
  ! magic numbers live in one place.
  real(kind=dp), parameter, public :: bdg_default_near_zero_frac = 0.001_dp
  real(kind=dp), parameter, public :: bdg_default_min_threshold = 1.0e-10_dp
  ! Magic-number SSOT for the slim Pfaffian floor (replaces the literal at
  ! topological_analysis.f90:1663, 1681, 1766). Promoted per ticket 02 of
  ! `.scratch/bdg-evaluator-pfaffian/`.
  real(kind=dp), parameter, public :: bdg_default_pfaffian_floor = 1.0e-12_dp

  ! Parameters for a single BdG evaluation.
  type :: bdg_eval_params_t
    real(kind=dp) :: delta_0        ! SC gap magnitude (eV) — scale for near-zero band
    real(kind=dp) :: near_zero_frac ! default 0.001; |E| < near_zero_frac*delta_0 counts as near-zero
    real(kind=dp) :: pfaffian_floor ! default 1.0e-12; |Pf| below this counts as zero (slim witness)
  end type

  ! Result of a single BdG evaluation.
  type :: bdg_eval_result_t
    real(kind=dp) :: minigap         ! 2 * minval(|E|)
    integer       :: near_zero_count ! count of |E| < near_zero_frac * delta_0
    integer       :: invariant_flag  ! 1 if near_zero_count >= 2, else 0
  end type

contains

  ! ==============================================================================
  ! Per-point BdG evaluator. Pure: identical input → identical output.
  !
  ! Returns the SC minigap, the count of eigenvalues inside the near-zero band,
  ! and a heuristic invariant flag (1 if at least one ±E pair sits in the band).
  !
  ! Empty spectrum is handled defensively (zero values), not as a fatal error;
  ! it is not on the hot path.
  ! ==============================================================================
  pure function eval_bdg_point(eigenvalues, params) result(r)
    real(kind=dp), intent(in), contiguous :: eigenvalues(:)
    type(bdg_eval_params_t), intent(in) :: params
    type(bdg_eval_result_t) :: r

    real(kind=dp) :: near_zero_threshold
    integer :: i

    r%minigap = 0.0_dp
    r%near_zero_count = 0
    r%invariant_flag = 0

    if (size(eigenvalues) < 1) return

    r%minigap = 2.0_dp * minval(abs(eigenvalues))

    near_zero_threshold = params%near_zero_frac * params%delta_0
    do i = 1, size(eigenvalues)
      if (abs(eigenvalues(i)) < near_zero_threshold) then
        r%near_zero_count = r%near_zero_count + 1
      end if
    end do

    if (r%near_zero_count >= 2) r%invariant_flag = 1
  end function eval_bdg_point

  ! ==============================================================================
  ! Factory: build a bdg_eval_params_t from a single delta_0 using the module
  ! defaults for near_zero_frac. Collapses the 5 call sites in main_topology
  ! from `bdg_eval_params_t(cfg%bdg%delta_0, 0.001_dp)` to a one-liner.
  !
  ! EXACTLY equivalent to bdg_eval_params_t(delta_0, bdg_default_near_zero_frac)
  ! — same defaults as the previous direct construction.
  ! ==============================================================================
  pure function bdg_eval_params_with_delta(delta_0) result(p)
    real(kind=dp), intent(in) :: delta_0
    type(bdg_eval_params_t) :: p
    p%delta_0        = delta_0
    p%near_zero_frac = bdg_default_near_zero_frac
    p%pfaffian_floor = bdg_default_pfaffian_floor
  end function bdg_eval_params_with_delta

  ! ==============================================================================
  ! QW near-zero tolerance helper. Returns the BdG near-zero threshold with a
  ! numerical floor (1e-10 eV) so callers that need an absolute precision
  ! floor (e.g., QW Majorana profile extraction at very small delta_0) get
  ! a consistent value via the seam rather than re-deriving the literal.
  !
  ! Equivalent to: max(1.0e-10_dp, params%near_zero_frac * abs(params%delta_0)).
  ! Abs protects against negative delta_0 (shouldn't occur but defensive).
  ! ==============================================================================
  pure function q_zero_tol(params) result(t)
    type(bdg_eval_params_t), intent(in) :: params
    real(kind=dp) :: t
    t = max(bdg_default_min_threshold, params%near_zero_frac * abs(params%delta_0))
  end function q_zero_tol

  ! ==============================================================================
  ! CSR BdG slim projected Pfaffian witness — seam sibling (wire-rung invariant).
  !
  ! Returns the S2-projected Pfaffian sign (s2_sign ∈ {-1, 0, +1}) of the
  ! BdG matrix H_bdg_csr. S2 = bands 7-8 per the k.p block table SSOT.
  !
  ! Per design decision 2026-07-13 (option b): this seam sibling is a thin
  ! wrapper that delegates the CSR-aware S2 extraction to
  ! wire_pfaffian_witness_sweep in topological_analysis.f90. That subroutine
  ! already operates on CSR input (per ticket 04 — `main_topology.f90:1371`
  ! migration pattern) and is the existing CSR-aware dense-path witness;
  ! reusing it keeps the seam thin (single subroutine import) without
  ! re-implementing the S2 row-extraction for CSR.
  !
  ! The params%pfaffian_floor field is the SSOT for future U13 work where
  ! the floor will be wired through; for now wire_pfaffian_witness_sweep
  ! uses its own internal floor (`bdg_default_pfaffian_floor`).
  !
  ! DEP EDGE: this is the one symbol `bdg_observables` imports from `topological_analysis`
  ! (Layer 3). Necessary because the helper encapsulates the CSR-aware S2
  ! indexing logic in one place; documented in src/physics/AGENTS.md Task 3.1.
  ! ==============================================================================
  function eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) result(s2_sign)
    type(csr_matrix), intent(in)        :: H_bdg_csr
    integer,          intent(in)        :: Nbdg
    type(bdg_eval_params_t), intent(in) :: params
    integer                              :: s2_sign

    ! Delegate to the CSR-aware dense-path witness. params is unused for
    ! now (params%pfaffian_floor plumbing is the U13 follow-up — the SSOT
    ! exists, the helper doesn't accept it yet).
    call wire_pfaffian_witness_sweep(H_bdg_csr, Nbdg, s2_sign)
  end function eval_bdg_pfaffian_witness_csr

end module bdg_observables
