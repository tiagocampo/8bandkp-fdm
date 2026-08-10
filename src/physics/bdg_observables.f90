module bdg_observables

  ! ==============================================================================
  ! BdG per-point observables — the foundational seam extracted from app glue.
  !
  ! Three seams, three consumers:
  !   - eval_bdg_point (per-point minigap/near-zero-count/invariant_flag):
  !     four call sites in main_topology (run_bdg_wire ×2, run_bdg_qw,
  !     eval_wire_bdg_gap) used to inline
  !       minigap: 2 * minval(abs(eigvals_bdg))
  !       near-zero threshold: 0.001 * delta_0
  !       invariant flag: count(|E| < threshold) >= 2
  !   - eval_bdg_pfaffian_witness_csr (slim projected Pfaffian S2, wire rung):
  !     one call site in main_topology (eval_wire_bdg_gap). PR #42 retired
  !     the dense wire_pfaffian_witness (S1+S2) and folded production into
  !     the seam sibling; U13 is the destination for full Bloch-Pfaffian.
  !   - eval_bdg_pfaffian_witness_product_csr (strict S1×S2 product, U13 T2):
  !     new sibling taking T1's Bloch stack and returning the L3-mapped
  !     strict z2 ∈ {-1, 0, +1} with optional disagreement_reason out-arg.
  !     Zero call sites in main_topology today (T3's dispatch wires it).
  !
  ! This module folds the per-point decision into one pure-function call so
  ! the build-and-solve stays in main_topology (per ADR 0003) while the
  ! per-point physics lives in one place. Downstream slices (Kitaev wrapper,
  ! polarization, LDOS) consume the same contract.
  !
  ! No I/O, no global state. The T2 product seam allocates a CSR scratch
  ! buffer for the dense→CSR conversion of slice 1 (private to the seam);
  ! the seam's other branches are pure.
  ! ==============================================================================

  use definitions, only: dp
  use sparse_matrices, only: csr_matrix, csr_build_from_coo, csr_free, dense_to_csr
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
  public :: eval_bdg_pfaffian_witness_product_csr
  public :: eval_bdg_kitaev_majorana
  public :: bdg_pfaffian_params_t
  public :: bdg_pfaffian_params_with_floor

  ! Module-level defaults (SSOT for the near-zero literals). Extracted from
  ! inline 0.001_dp / 1.0e-10_dp at the 5 call sites in main_topology so the
  ! magic numbers live in one place.
  real(kind=dp), parameter, public :: bdg_default_near_zero_frac = 0.001_dp
  real(kind=dp), parameter, public :: bdg_default_min_threshold = 1.0e-10_dp
  ! Magic-number SSOT for the slim Pfaffian floor (replaces the literal at
  ! topological_analysis.f90:1663, 1681, 1766). Promoted per ticket 02 of
  ! `.scratch/archive/bdg-evaluator-pfaffian/`.
  real(kind=dp), parameter, public :: bdg_default_pfaffian_floor = 1.0e-12_dp

  ! Parameters for a single BdG evaluation (per-point minigap/near-zero-count).
  ! No Pfaffian floor here — that belongs to the Pfaffian witness, not the
  ! per-point evaluator (SRP ticket 01 of .scratch/bdg-u2-actual-ship/).
  type :: bdg_eval_params_t
    real(kind=dp) :: delta_0        ! SC gap magnitude (eV) — scale for near-zero band
    real(kind=dp) :: near_zero_frac ! default 0.001; |E| < near_zero_frac*delta_0 counts as near-zero
  end type

  ! Parameters for the slim projected Pfaffian witness seam sibling. Single
  ! field — the |Pf| floor below which a Pfaffian counts as zero. Defaulted to
  ! the SSOT so the bare structure constructor `bdg_pfaffian_params_t()`
  ! returns the safe default; the factory `bdg_pfaffian_params_with_floor` is
  ! the validation site (rejects zero/negative floors).
  type :: bdg_pfaffian_params_t
    real(kind=dp) :: pfaffian_floor = bdg_default_pfaffian_floor
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
  ! Factory: build a bdg_pfaffian_params_t validated against the SSOT.
  ! Optional pfaffian_floor overrides the SSOT; zero/negative floors are a
  ! fatal config error (error stop) — a non-positive floor would silently
  ! classify every Pfaffian as non-zero. The bare structure constructor
  ! `bdg_pfaffian_params_t()` is safe (returns the SSOT) but does NOT validate;
  ! this factory is the validation site.
  ! ==============================================================================
  function bdg_pfaffian_params_with_floor(pfaffian_floor) result(p)
    real(kind=dp), intent(in), optional :: pfaffian_floor
    type(bdg_pfaffian_params_t) :: p

    if (present(pfaffian_floor)) then
      if (pfaffian_floor <= 0.0_dp) then
        error stop 'bdg_pfaffian_params_t: pfaffian_floor must be > 0'
      end if
      p%pfaffian_floor = pfaffian_floor
    else
      p%pfaffian_floor = bdg_default_pfaffian_floor
    end if
  end function bdg_pfaffian_params_with_floor

  ! ==============================================================================
  ! CSR BdG slim projected Pfaffian witness — seam sibling (wire-rung invariant).
  !
  ! Returns the S2-projected Pfaffian sign (s2_sign ∈ {-1, 0, +1}) of the
  ! BdG matrix H_bdg_csr. S2 = bands 7-8 per the k.p block table SSOT.
  !
  ! User Story 1 seam contract: the wire rung's invariant_flag is the slim
  ! projected Pfaffian witness, S2 = bands 7-8 (per hamiltonian_blocks.f90
  ! SSOT). Callers pass the Pfaffian floor via bdg_pfaffian_params_t; this is
  ! the route the declared-but-not-consumed SSOT at :45 was missing in PR #42.
  !
  ! Per design decision 2026-07-13 (option b): this seam sibling is a thin
  ! wrapper that delegates the CSR-aware S2 extraction to
  ! wire_pfaffian_witness_sweep in topological_analysis.f90. That subroutine
  ! already operates on CSR input (per ticket 04 — `main_topology.f90:1371`
  ! migration pattern) and is the existing CSR-aware dense-path witness;
  ! reusing it keeps the seam thin (single subroutine import) without
  ! re-implementing the S2 row-extraction for CSR.
  !
  ! Non-pure by intent: the call chain
  !   eval_bdg_pfaffian_witness_csr → wire_pfaffian_witness_sweep → complex_pfaffian
  ! ends in src/math/pfaffian.f90, whose helpers are not pure. Making this
  ! sibling pure would require pure-marking that whole module — a separate PR
  ! out of this map's scope (ticket 01 sub-decision 3; the L3-import exception
  ! is the topological_analysis import below, documented in AGENTS.md Task 3.1
  ! and forwarded to Codacy triage ticket 04).
  !
  ! U13 forward reference: the slim S2 witness is the interim stand-in for the
  ! full Bloch-Pfaffian sweep (S1+S2 strict sign agreement). S1 needs the
  ! periodic/Bloch BdG construction deferred to U13 — Issue 05 of the parent
  ! plan. Until then the seam accepts s2 ∈ {-1, 0, +1} with s2 /= 0 on
  ! non-diagonal synthetic fixtures as the GREEN contract (User Story 5).
  ! ==============================================================================
  function eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params, best_pf_abs) result(s2_sign)
    type(csr_matrix), intent(in)            :: H_bdg_csr
    integer,          intent(in)            :: Nbdg
    type(bdg_pfaffian_params_t), intent(in) :: params
    ! Optional `best_pf_abs` out-arg threads the per-(B, mu) max-site |Pf|
    ! magnitude through the seam for callers accumulating a per-B `min-|Pf|`
    ! proxy (U10 T1a/T1b); absent -> magnitude discarded, existing API contract
    ! preserved.
    real(kind=dp), intent(out), optional    :: best_pf_abs
    integer                                  :: s2_sign
    real(kind=dp)                           :: best_pf_abs_local

    ! Delegate to the CSR-aware dense-path witness, threading the SSOT
    ! pfaffian_floor through (PR #42's declared-but-not-consumed gap, closed
    ! by ticket 01 of .scratch/bdg-u2-actual-ship/).
    if (present(best_pf_abs)) then
      call wire_pfaffian_witness_sweep(H_bdg_csr, Nbdg, params%pfaffian_floor, s2_sign, best_pf_abs)
    else
      call wire_pfaffian_witness_sweep(H_bdg_csr, Nbdg, params%pfaffian_floor, s2_sign, best_pf_abs_local)
    end if
  end function eval_bdg_pfaffian_witness_csr

  ! ==============================================================================
  ! QW+Kitaev-rung Majorana number seam sibling.
  !
  ! Wraps `pfaffian.f90:kitaev_majorana_number` with a seam-shape signature.
  ! Returns majorana_number ∈ {-1, 0, +1} (Kitaev 2001, Eq. 26 convention:
  ! M = -1 topological, M = +1 trivial). The wrapped helper uses the
  ! Lutchyn-Oreg sign-of-det (ADR 0008 §1) formula via polar decomposition.
  !
  ! The advanced `omega_struct` argument is intentionally NOT exposed on the
  ! seam — consumers that need a custom particle-hole structure (non-canonical
  ! BdG) can fall back to the helper directly. The canonical Kitaev form
  ! (`omega = σ_y ⊗ I_N`) is the only one in scope for U2/U3 consumers.
  ! ==============================================================================
  function eval_bdg_kitaev_majorana(H_k_array, k_par_values) result(majorana_number)
    complex(kind=dp), intent(in) :: H_k_array(:,:,:)
    real(kind=dp),    intent(in) :: k_par_values(:)
    integer                       :: majorana_number

    majorana_number = kitaev_majorana_number(H_k_array, k_par_values)
  end function eval_bdg_kitaev_majorana

  ! ==============================================================================
  ! U13 T2 — strict S1×S2 product seam (wire rung, Bloch-periodic).
  !
  ! Wires the k-product Majorana number S1 (from `kitaev_majorana_number` over
  ! T1's Bloch stack) with the slim projected Pfaffian S2 (from
  ! `eval_bdg_pfaffian_witness_csr` at slice 1, the canonical reference kz
  ! `k_par_values(1)`), and applies the L3 disagreement mapping from T6's chart
  ! (`.scratch/bdg-u13-bloch-pfaffian/issues/06-grill-loose-fog.md`):
  !
  !   case i   (S1 = S2 = ±1)              → z2 = S1, reason = 0 (none)
  !   case ii  (S1 = S2 = 0)               → z2 = 0,  reason = 0 (none)
  !   case iii (S1 ≠ S2 in sign)           → z2 = 0,  reason = 1 (s1_s2_sign_split)
  !   case iv  (S1 = ±1, S2 = 0)           → z2 = 0,  reason = 2 (s2_fest_saturation)
  !   case v   (S1 = 0,  S2 = ±1)          → z2 = 0,  reason = 3 (s1_closure_s2_misses)
  !
  ! The disagreement_reason integer is the writer-facing API for T3's witness
  ! file (per the L3 inheritance note in T6); absent out-arg discards it. The
  ! best_pf_abs out-arg mirrors the slim seam's contract — per-(B, μ) max-site
  ! |Pf| from the band-7/8 subblock projection, accumulated by the gate for the
  ! degeneracy-detection branch U10 T6 added.
  !
  ! Slice 1 is the canonical S2 reference kz (k_par_values(1), typically Γ).
  ! The wire's free-z is the Bloch-periodic axis per T6 L2; slice 1 is the
  ! reference point of the lattice the gate compares the strict product
  ! against. This is the slim-seam "one point" approximation lifted to the
  ! strict k-product invariant.
  !
  ! The dense→CSR conversion of slice 1 uses the public `dense_to_csr`
  ! helper in `sparse_matrices` (the inverse of `csr_to_dense_work`),
  ! shared with U13 T3's `eval_wire_bdg_gap_bloch` dispatch path.
  ! ==============================================================================
  function eval_bdg_pfaffian_witness_product_csr(H_k_array, k_par_values, params, &
                                                   best_pf_abs, disagreement_reason) result(z2)
    complex(kind=dp), intent(in), contiguous :: H_k_array(:,:,:)
    real(kind=dp),    intent(in)             :: k_par_values(:)
    type(bdg_pfaffian_params_t), intent(in)  :: params
    real(kind=dp), intent(out), optional     :: best_pf_abs
    integer, intent(out), optional           :: disagreement_reason
    integer :: z2

    integer :: s1, s2, n_full
    type(csr_matrix) :: H_csr_slice
    real(kind=dp) :: best_pf_abs_local

    ! --- S1: k-product Majorana number from the dense Bloch stack ---
    s1 = kitaev_majorana_number(H_k_array, k_par_values)

    ! --- S2: slim projected Pfaffian from slice 1 (canonical reference kz) ---
    ! Dense→CSR conversion of the first slice. The slice is square (n_full×n_full)
    ! and dense (16N×16N for an N-site wire); only the nonzero entries are emitted.
    n_full = size(H_k_array, 1)
    if (n_full /= size(H_k_array, 2)) then
      ! Shape contract violated — return z2 = 0 with sign-split reason (the
      ! most conservative disagreement branch). The T2 unit test never hits this;
      ! it's a defensive guard for malformed input.
      z2 = 0
      if (present(disagreement_reason)) disagreement_reason = 1
      if (present(best_pf_abs)) best_pf_abs = 0.0_dp
      return
    end if

    call dense_to_csr(H_csr_slice, H_k_array(:, :, 1))

    s2 = eval_bdg_pfaffian_witness_csr(H_csr_slice, n_full, params, best_pf_abs_local)

    call csr_free(H_csr_slice)

    if (present(best_pf_abs)) best_pf_abs = best_pf_abs_local

    ! --- L3 disagreement mapping ---
    z2 = map_s1_s2_to_z2(s1, s2, disagreement_reason)

  end function eval_bdg_pfaffian_witness_product_csr

  ! ==============================================================================
  ! L3 disagreement mapping — single source of truth for the S1×S2 → z2 logic.
  !
  ! Extracted so the seam's body stays short (under the 50-line budget) and the
  ! mapping itself can be unit-tested in isolation if needed.
  !
  ! Argument conventions match `kitaev_majorana_number` ({-1, 0, +1}) and
  ! `eval_bdg_pfaffian_witness_csr` ({-1, 0, +1}). Per T6 L3, disagreement
  ! emits z2 = 0 with a disagreement_reason code; strict agreement yields
  ! z2 = S1 (case i) or z2 = 0 (case ii). The five L3 cases are enumerated
  ! in the function header at `eval_bdg_pfaffian_witness_product_csr` above.
  ! ==============================================================================
  function map_s1_s2_to_z2(s1, s2, disagreement_reason) result(z2)
    integer, intent(in) :: s1, s2
    integer, intent(out), optional :: disagreement_reason
    integer :: z2

    if (s1 == s2) then
      ! case i (both ±1) or case ii (both 0): strict agreement, no warn.
      z2 = s1
      if (present(disagreement_reason)) disagreement_reason = 0
    else if (s1 == 0) then
      ! case v: S1 detects closure, S2 has a clean sign at slice 1. z2 = 0
      ! (closure is indeterminate). disagreement_reason = s1_closure_s2_misses.
      z2 = 0
      if (present(disagreement_reason)) disagreement_reason = 3
    else if (s2 == 0) then
      ! case iv: S1 strict-clean, S2 at floor (FEST saturation or Pf ≈ 0).
      ! z2 = 0 (don't trust S2's silent zero as a clean answer).
      ! disagreement_reason = s2_fest_saturation.
      z2 = 0
      if (present(disagreement_reason)) disagreement_reason = 2
    else
      ! case iii: S1 and S2 disagree on sign. Physically impossible in
      ! canonical PHS; signals sampling/build bug. z2 = 0 with reason =
      ! s1_s2_sign_split (the only L3 branch reserved for a possible follow-up
      ! `error stop` guard per T6's case-iii reservation).
      z2 = 0
      if (present(disagreement_reason)) disagreement_reason = 1
    end if
  end function map_s1_s2_to_z2

end module bdg_observables
