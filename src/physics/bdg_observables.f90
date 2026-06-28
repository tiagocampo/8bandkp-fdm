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

  implicit none

  private

  public :: bdg_eval_params_t
  public :: bdg_eval_result_t
  public :: eval_bdg_point

  ! Parameters for a single BdG evaluation.
  type :: bdg_eval_params_t
    real(kind=dp) :: delta_0        ! SC gap magnitude (eV) — scale for near-zero band
    real(kind=dp) :: near_zero_frac ! default 0.001; |E| < near_zero_frac*delta_0 counts as near-zero
    real(kind=dp) :: invariant_tol  ! reserved for future use (current logic is exact); default 1e-10
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

end module bdg_observables