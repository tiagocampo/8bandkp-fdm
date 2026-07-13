export const meta = {
  name: 'eigensolver-hardening-implement',
  description: 'Harden eigensolver standardization: 6 surgical fixes, zero physics change',
  phases: [
    { title: 'Sequenced', detail: 'slices 01->02 (shared defs.f90/main.f90/simulation_setup.f90)' },
    { title: 'Parallel', detail: 'slices 03/04/05/06 (disjoint modules, edit-only)' },
  ],
}

// ---------------------------------------------------------------------------
// Shared context baked into every implementer prompt.
// ---------------------------------------------------------------------------
const PREAMBLE = `
You are implementing ONE slice of the "Eigensolver Standardization — Hardening (Review Fixes Round 2)" refactor.
Repo: /data/8bandkp-fdm  (cwd is the repo root). Branch: refactor/eigensolver-standardization (already checked out — do NOT switch branches).

INVARIABLE RULE — ZERO physics-value change. The golden regression suite MUST stay bit-identical. Your change is a robustness/consistency fix only. Do NOT alter any computed eigenvalue, eigenvector, g-factor, optical spectrum, or SC potential on any supported configuration.

GLOBAL CONSTRAINTS (from CLAUDE.md — follow exactly):
- Fortran F2018 enforced (-std=f2018, NO GNU extensions). Prefer generic intrinsics (sqrt, not dsqrt/dble). Prefer do/do concurrent over forall. Use error stop (with a descriptive message) for fatal errors, never bare 'stop 1'.
- Modules use 'private' default with explicit 'public ::' exports. Add 'public :: resolve_solver_defaults' (slice 01) etc. when you add a public symbol.
- Scalar pure helpers should be 'elemental pure' where applicable (F2008 requires both keywords).
- 'contiguous' attribute on hot-path assumed-shape array args (not on optional/allocatable).
- declaration-before-use in array dimension expressions.
- Follow file/function length guidelines: ~300 lines/file, ~50 lines/function.

READ-BEFORE-EDIT (mandatory):
- BEFORE editing anything under src/physics/, read src/physics/AGENTS.md (basis ordering, kpterms sign convention, SSOT tables, anti-patterns).
- BEFORE editing anything under tests/integration/, read tests/integration/AGENTS.md (test patterns, COVERAGE annotations, tolerance tiers, anti-patterns: never create configs at runtime in shell tests — use tests/regression/configs/; put COVERAGE annotations in verify_*.py not .sh wrappers).

defs.f90 HARD RULE (approval-gated per CLAUDE.md): you may add a PURE function/subroutine and/or a validate() check ONLY. NO derived-type changes. defs.f90 is the ROOT module — it must NOT 'use' eigensolver or any downstream module (the dependency DAG flows defs <- ... <- eigensolver). So any helper you add there operates on primitive args (character strings, integers) and returns primitives.

OUT OF SCOPE — do NOT touch: material parameters (parameters.f90), basis ordering, Bir-Pikus signs, the Zeeman table, FD stencil coefficients (finitedifferences.f90), the Poisson solver (poisson.f90), SC-loop CONVERGENCE logic (the DIIS/mixing/iteration logic — slice 05 only clamps the RESULT WRITE, not convergence), the wire path's own value-only fast-path logic (hamiltonian_wire.f90).

CRITICAL — DO NOT RUN cmake OR ctest. The build/test is centralized AFTER all six slices land. Reason: Fortran builds share one build/ directory; concurrent 'cmake --build' calls clobber each other. Other slices are editing concurrently right now. Instead: author carefully, then RE-READ your complete diff (git diff), and reason explicitly about why every supported configuration's output stays bit-identical. You MAY use grep/Read to confirm symbols, call sites, and that no golden config is affected.

FILE OWNERSHIP — edit ONLY the files named in YOUR slice below. Other slices own other files. Do not edit files outside your list.

METHODOLOGY — test-driven-development discipline:
- For slices adding a NEW test (01, 02, 06): write the test first expressing the desired behavior, then the implementation that makes it pass.
- For surgical slices with no new test (03, 04, 05): the existing golden regression + equivalence suite is your GREEN gate. Reason step-by-step about why the surgical edit cannot change any computed value.

When done, return your structured result (the schema forces it). Be precise and honest in 'bit_identical_reasoning' and 'acceptance_criteria'.
`

const SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['slice', 'files_changed', 'summary', 'acceptance_criteria', 'bit_identical_reasoning', 'risks_or_notes'],
  properties: {
    slice: { type: 'string' },
    files_changed: { type: 'array', items: { type: 'string' } },
    summary: { type: 'string', description: 'What you changed and why, 3-8 sentences.' },
    acceptance_criteria: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['criterion', 'met', 'evidence'],
        properties: {
          criterion: { type: 'string' },
          met: { type: 'boolean' },
          evidence: { type: 'string', description: 'file:line or test name proving it' },
        },
      },
    },
    bit_identical_reasoning: { type: 'string', description: 'Step-by-step why no supported config output changes. Name configs you checked.' },
    risks_or_notes: { type: 'string' },
  },
}

// ---------------------------------------------------------------------------
// SLICE 01 — Method-aware AUTO resolution (deep module). RISKIEST. Sequenced.
// ---------------------------------------------------------------------------
const SLICE01 = PREAMBLE + `
=== YOUR SLICE: 01 — Method-aware AUTO resolution (deep module) ===
Parent issue: .scratch/eigensolver-hardening/issues/01-method-aware-auto-resolution.md  (READ IT). Covers PRD user stories 5-10. This is the riskiest slice.

GOAL: Extract a single PURE resolver that turns (confinement, method, mode) into a resolved (method, mode), make it the single source of truth for solver dispatch, and Landau stops forcing INDEX.

THE CONTRACT (codified in root CONTEXT.md under "AUTO"):
- Method resolves FIRST: if method=='AUTO' -> confinement default (bulk->DENSE, qw->DENSE, landau->DENSE, wire->FEAST); else -> the user's method (trimmed).
- Mode then resolves METHOD-AWARE: if mode!='AUTO' -> honor it (passthrough, trimmed/uppercased); if mode=='AUTO' -> a default COMPATIBLE with the resolved method:
    * resolved method FEAST  -> 'ENERGY'
    * resolved method DENSE  -> confinement native mode: bulk->'FULL', qw->'INDEX', landau->'INDEX', wire->'ENERGY'
- INVARIANT: AUTO never yields an invalid combination. FEAST+INDEX is unreachable through AUTO. (Explicit FEAST+INDEX set by the user is still rejected — by existing validate() check I15 in defs.f90 ~line 815 and eigensolver_config_validate — but AUTO must never produce it.)

DESIGN DECISION (you MUST follow this so the unit test and dispatch agree):
- Place 'resolve_solver_defaults' in src/core/defs.f90 as a PUBLIC PURE SUBROUTINE (defs.f90 cannot 'use' eigensolver, so it cannot return EIGEN_MODE_* integer constants — it returns resolved UPPERCASE character strings). Signature:
    pure subroutine resolve_solver_defaults(confinement, method, mode, out_method, out_mode)
      character(len=*), intent(in)  :: confinement, method, mode
      character(len=*), intent(out) :: out_method, out_mode   ! 'DENSE'/'FEAST' and 'FULL'/'INDEX'/'ENERGY'
  (Rationale: returning a pair of strings via a pure subroutine matches the "(method,mode)" pair semantics and avoids adding any derived type to defs.f90, which is approval-gated. Normalize inputs with trim()/uppercasing consistent with existing validate() string handling.)
- Each dispatch block then calls it and maps out_mode string -> EIGEN_MODE_* constant. To avoid re-duplicating that 4-way map in four places, add ONE small private helper in src/core/simulation_setup.f90, e.g. a pure function eigen_mode_from_string(s) -> integer, used by all blocks (this helper MAY live in simulation_setup.f90 which already 'use's the eigensolver module and its EIGEN_MODE_* constants).

CURRENT CODE — the four duplicated dispatch blocks you are consolidating (src/core/simulation_setup.f90):
- BULK  ~lines 100-118: method AUTO->DENSE; mode select: AUTO/FULL->FULL, INDEX->INDEX, ENERGY->ENERGY; then nev/il/iu=8/1/8.  (Bulk AUTO mode resolves to FULL.)
- QW    ~lines 178-196: method AUTO->DENSE; mode select: AUTO/INDEX->INDEX, FULL->FULL, ENERGY->ENERGY; nev=min(iuu-il+1,N). (QW AUTO -> INDEX.)
- WIRE  ~lines 255-269: method AUTO->FEAST; mode select: AUTO/ENERGY->ENERGY, FULL->FULL, INDEX->INDEX. (Wire AUTO -> ENERGY.)
- LANDAU~lines 341-352: method AUTO->DENSE; mode FORCED to EIGEN_MODE_INDEX unconditionally (line ~347) — THIS IS THE BUG: it ignores the user's [solver] mode. Your resolver fixes it: landau AUTO->DENSE, AUTO mode -> INDEX (DENSE native for landau), but an explicit mode (e.g. ENERGY) is now honored.
Replace each block's method/mode resolution with a call to resolve_solver_defaults(cfg%confinement, cfg%solver%method, cfg%solver%mode, m, mm) then eigen_mode_from_string(mm). KEEP all the confinement-specific nev/il/iu bookkeeping EXACTLY as-is. KEEP the subsequent eigensolver_config_validate + make_eigensolver calls.

ALSO route the 5th resolution site through it for consistency (PRD user stories 9-10: single source of truth): src/apps/main.f90 ~lines 668-700 (the band-structure sweep config: bulk FULL / QW INDEX). Replace its method+mode selects with the same resolver call + eigen_mode_from_string. IMPORTANT: main.f90 does NOT 'use definitions' for this necessarily — check its 'use' statements; resolve_solver_defaults is in defs.f90's module 'definitions', so add the 'use' if missing. Verify the resolved bulk=DENSE/FULL and QW=DENSE/INDEX match today's behavior EXACTLY (they must — bulk sweep uses FULL, QW sweep uses INDEX).

BIT-IDENTICAL ANALYSIS YOU MUST DO AND DOCUMENT:
The resolution change only affects FEAST+AUTO combinations:
- wire FEAST+AUTO: old -> ENERGY, new -> ENERGY  (UNCHANGED).
- qw  FEAST+AUTO: old -> INDEX => FEAST+INDEX => ABORT today (not a working config); new -> ENERGY (the FIX, user story 5).
- landau FEAST+AUTO: old -> forced INDEX => ABORT today; new -> ENERGY (the FIX, user story 6).
- bulk FEAST+AUTO: old -> FULL (all 8); new -> ENERGY.
VERIFY no golden/regression config uses bulk+FEAST or otherwise relies on the bulk-FEAST-FULL behavior: grep tests/regression/configs/ for method = "FEAST" and confirm every hit is confinement = "wire". If all FEAST configs are wire, bit-identical is preserved (wire unchanged; qw/landau/bulk FEAST+AUTO were either aborting or absent from the golden suite). Document this explicitly.

UNIT TEST — land it WITH the extraction (mandatory). Add a table-driven test covering every confinement x method x mode:
- confinements: bulk, qw, wire, landau.
- methods: AUTO, DENSE, FEAST.
- modes: AUTO, FULL, INDEX, ENERGY.
For each combo, assert the resolved (out_method, out_mode) matches the contract above, AND assert the INVARIANT: whenever mode=='AUTO', the resolved pair is NEVER ('FEAST','INDEX'). (Explicit FEAST+INDEX via mode='INDEX' method='FEAST' returns ('FEAST','INDEX') — note in the test that this is the caller/validate's responsibility to reject, not the resolver's; the resolver faithfully honors explicit input.)
Place the test by EXTENDING tests/unit/test_eigensolver.pf (it already 'use definitions' and is already registered in tests/CMakeLists.txt ~line 70 — this avoids CMake edits and CMake-ownership conflicts). Add @test subroutines there. If you strongly prefer a separate file test_resolve_solver_defaults.pf, you MUST also register it in tests/CMakeLists.txt (you are the only slice touching tests/CMakeLists.txt, so that is allowed) — but extending test_eigensolver.pf is preferred.

ACCEPTANCE CRITERIA (from the issue — each must map to your diff):
1. Pure resolve_solver_defaults resolves (confinement,method,mode)->(method,mode); AUTO is method-aware as specified.
2. All four dispatch blocks (bulk/qw/wire/landau in simulation_setup.f90) call it; duplicated per-block method/mode resolution removed. (Plus main.f90:668-700 routed through it.)
3. Landau honors an explicit mode (no longer forces INDEX).
4. (Verified later by central ctest) QW/Landau with method=FEAST and no mode runs successfully (FEAST+ENERGY).
5. Unit test table-drives every confinement x method x mode and asserts the no-invalid-combo invariant.
6. (Verified later) Golden regression green; existing FEAST+INDEX rejection tests still pass.

OWNED FILES (edit ONLY these): src/core/defs.f90, src/core/simulation_setup.f90, src/apps/main.f90, tests/unit/test_eigensolver.pf  (+ tests/CMakeLists.txt ONLY IF you create a new .pf file).
`

// ---------------------------------------------------------------------------
// SLICE 02 — Bulk full-diagonalization completion. Sequenced AFTER 01.
// ---------------------------------------------------------------------------
const SLICE02 = PREAMBLE + `
=== YOUR SLICE: 02 — Bulk full-diagonalization completion ===
Parent issue: .scratch/eigensolver-hardening/issues/02-bulk-full-diagonalization.md  (READ IT). Covers PRD user stories 1-4. Blocked by slice 01 (both edit main.f90 / defs.f90) — slice 01 has ALREADY landed; build on its result.

GOAL: Finish the bulk full-diagonalization migration the refactor left half-done. Bulk now fully diagonalizes (returns all 8 eigenvalues) but storage is still sized by the *requested band window* (cfg%evnum), so a bulk config asking for <8 bands writes out of bounds.

TWO PARTS:

PART A — Size bulk storage to 8 (src/apps/main.f90):
- ~line 577-579: 'if bulk then il=1; iuu=cfg%evnum'. For bulk, set iuu=8 (bulk is always 8 bands). Keep the QW branch (il/iuu from num_vb/num_cb) UNTOUCHED.
- ~line 590: allocate(eig(iuu-il+1, nsteps)) — for bulk this is now 8.
- ~line 592-596: bulk branch allocates eigv(8, cfg%evnum, nsteps). Change the bulk eigv allocation to size 8 in the band dimension: eigv(8, 8, nsteps) (i.e. use 8, not cfg%evnum, for the band index extent). Keep the QW branch (eigv(N, iuu-il+1, nsteps)) UNTOUCHED.
- The bulk write at ~line 711-712: eig(1:nev_found,k)=result_bs%eigenvalues(1:nev_found); eigv(:,1:nev_found,k)=result_bs%eigenvectors. With FULL mode nev_found=8 and storage now 8, this is safe.
- The bulk writeEigenfunctions at ~line 719 uses min(cfg%evnum,8) — leave as-is (harmless; with evnum=8 it's 8). Do not change QW write paths.

PART B — Tighten validation (src/core/defs.f90), per ADR 0002 inside validate():
- The existing check "V1: bulk evnum must not exceed 8" (~defs.f90 lines 666-678) currently only rejects evnum > 8. TIGHTEN it to require evnum == 8 for bulk (reject evnum /= 8 with a contextual error stop). Keep the message style consistent with neighboring checks (use the character buffer write pattern). This converts a silent out-of-bounds write into an up-front error stop.
- IMPORTANT ORDERING: place/keep this check so it does NOT preempt the earlier num_cb>=1 / num_vb>=1 / evnum==num_cb+num_vb checks (those must still fire first for the existing rejection-test cases V2/V3 which use num_cb=0 or num_vb=0). Verify the existing rejection cases still hit their intended messages.

PART C — Extend the rejection shell script (tests/integration/test_validate_rejects_bad_configs.sh):
- Add a new case (e.g. V16) after V15: a bulk config with num_cb + num_vb /= 8 (e.g. num_cb=2, num_vb=5 => evnum=7) that must now be rejected. Follow the EXACT existing pattern (cat > input.toml heredoc; run_test "V16_bulk_band_count_ne8" "evnum" — match a substring of your new error message). Verify NONE of the existing bulk cases (which use num_cb=2,num_vb=6 => evnum=8) newly fail.
- Note: this script is already registered in tests/CMakeLists.txt (~line 1103) — NO CMake edit needed.

BIT-IDENTICAL ANALYSIS: Issue acceptance states existing bulk configs use band count 8. VERIFY: grep tests/regression/configs/ for confinement = "bulk" and confirm every bulk config has num_cb+num_vb == 8 (evnum 8). If so, Part A (sizing to 8) changes nothing for them (storage was already effectively 8), and Part B (reject !=8) never fires on supported configs. Document the configs you checked.

ACCEPTANCE CRITERIA:
1. Bulk eigenvalue/eigenvector arrays sized to 8.
2. Bulk write path stores all 8 eigenpairs safely regardless of requested window.
3. validate() rejects a bulk config whose band count /= 8 with a contextual error stop.
4. Bad-config rejection script extended with a case asserting the new error fires.
5. (Central ctest) Golden regression green; existing bulk configs (evnum 8) unchanged.

OWNED FILES (edit ONLY these): src/apps/main.f90, src/core/defs.f90, tests/integration/test_validate_rejects_bad_configs.sh.
(Slice 01 already finished editing main.f90's resolution block ~668-700 and defs.f90's resolver; your main.f90 edit is the bulk storage sizing ~577-596, and your defs.f90 edit is the V1 check ~666-678 — distinct regions. Do not touch slice 01's resolver or resolution blocks.)
`

// ---------------------------------------------------------------------------
// SLICE 03 — QW fast path -> single rebuild path. Parallel.
// ---------------------------------------------------------------------------
const SLICE03 = PREAMBLE + `
=== YOUR SLICE: 03 — QW fast path -> single rebuild path ===
Parent issue: .scratch/eigensolver-hardening/issues/03-qw-fast-path-single-rebuild.md  (READ IT). Covers PRD user stories 11-13. Isolated to ONE module.

GOAL: In the QW CSR k-sweep fast path, DELETE the dead value-only branch; keep the SINGLE rebuild path that reuses the cached COO->CSR sort map (the real O(NNZ) vs O(NNZ log NNZ) win). The wire path's own value-only logic is UNTOUCHED (separate file, proven in production).

TARGET FILE: src/physics/hamiltonian_qw.f90  (READ src/physics/AGENTS.md FIRST).

CURRENT CODE — the fast-path rebuild dispatch (~lines 374-388):
    if (HT_csr%nnz > 0 .and. ws%coo_cache%initialized) then
      call csr_set_values_from_coo(HT_csr, coo_idx, ws%coo_cache%coo_to_csr(1:coo_idx), ws%coo_vals(1:coo_idx))   ! <-- VALUE-ONLY branch: DELETE
    else
      call csr_build_from_coo_cached(HT_csr, Ntot, Ntot, coo_idx, ws%coo_rows(1:coo_idx), ws%coo_cols(1:coo_idx), ws%coo_vals(1:coo_idx), ws%coo_cache%coo_to_csr)   ! <-- KEEP (rebuild via cached sort map)
      ws%coo_cache%coo_nnz_in = coo_idx
      ws%coo_cache%initialized = .true.
    end if
Replace the if/else with the rebuild path UNCONDITIONALLY (always csr_build_from_coo_cached + the two cache bookkeeping lines). Update the preceding comment block (~lines 319-330 and ~374-378) to describe the single rebuild path and remove references to the value-only branch / the nnz>0 distinction.

WHY THIS IS SAFE (the dead-branch argument — verify and document):
- The active QW FEAST k-sweep (src/apps/main.f90 ~lines 742-758) calls csr_free(HT_csr_loc) every iteration, so HT_csr%nnz==0 on every entry to ZB8bandQW_csr => it ALWAYS takes the else (rebuild) branch today. The value-only (nnz>0) branch is never executed in any current call path.
- The value-only branch is a latent footgun: a future caller reusing a persistent HT_csr across a Gamma->off-Gamma transition would hit csr_set_values_from_coo with a structure-mismatched sort map (Gamma-point structure is sparser than the sentinel-recorded map) => corruption or stop 1. Removing it makes the fast path safe-by-construction.
- The cached SORT MAP (ws%coo_cache%coo_to_csr) is STILL reused by csr_build_from_coo_cached — that is the actual performance win; you are only dropping the marginal "skip rowptr/colind rebuild", which no current caller avoids anyway.

DO NOT TOUCH:
- The SLOW path (~lines 485-525, the finalize_coo_to_csr dispatch and the sentinel-k workspace initialization). That is a different dispatch and is correct as-is.
- src/physics/hamiltonian_wire.f90 (the wire path's value-only logic) — out of scope.
- The helper subroutines update_kp_term_values, insert_*, csr_build_from_coo_cached, csr_set_values_from_coo themselves (do not delete csr_set_values_from_coo if other code references it — grep first; if it becomes unused after your edit, leaving it is fine, or note it).

NO NEW TEST (per PRD: covered by the existing QW FEAST fast-vs-dense equivalence test, verification_qw_sparse_solver). Your GREEN gate is that equivalence test staying green — reason about why deleting an unreachable branch cannot change any computed CSR matrix (the rebuild path it keeps is byte-identical to what runs today).

ACCEPTANCE CRITERIA:
1. Value-only branch removed from the QW CSR fast path; it always rebuilds via the cached sort map.
2. (Central ctest) Existing fast-vs-dense equivalence test (verification_qw_sparse_solver) stays green — proves no behavior change.
3. (Central ctest) Golden regression green.

OWNED FILES (edit ONLY this): src/physics/hamiltonian_qw.f90.
`

// ---------------------------------------------------------------------------
// SLICE 04 — FEAST subspace-size persistence. Parallel.
// ---------------------------------------------------------------------------
const SLICE04 = PREAMBLE + `
=== YOUR SLICE: 04 — FEAST subspace-size persistence ===
Parent issue: .scratch/eigensolver-hardening/issues/04-feast-m0-persistence.md  (READ IT). Covers PRD user stories 14-15. Isolated to the eigensolver module.

GOAL: Stop solve_feast from throwing away its hard-won subspace size M0 on every k-point. Today M0 is recomputed from config%m0 on every entry (~line 225-228); after an info=3 retry grows M0 and succeeds, the next call resets to the original => the workspace cache never matches and the retry repeats from scratch on every point.

TARGET FILE: src/math/eigensolver.f90  (subroutine solve_feast, ~lines 188-371; type feast_workspace ~lines 101-111).

FIX:
- Persist the WINNING M0 on the feast_workspace so the next call seeds from it. The feast_workspace type already has an 'M0' field (~line 105) used for the pattern-matching cache (fw%M0 in feast_workspace_matches_pattern ~line 157). Add a SEPARATE field for the persisted seed (e.g. 'last_successful_m0') — this is a derived-type change to eigensolver.f90's feast_workspace, which is ALLOWED (the defs.f90 derived-type restriction does NOT apply to eigensolver.f90). Initialize it to 0. Keep the existing finalizer/feast_workspace_free correct (reset the new field there too ~line 458-469).
- On ENTRY to solve_feast, after computing the base M0 from config (~line 225-228: M0=config%m0; if <=0 then 2*nev; max(nev+1); min(N)), SEED from the persisted value: if present(fw) and fw%last_successful_m0 > 0, set M0 = max(M0_base, fw%last_successful_m0), then re-apply the caps (max(nev+1), min(N)). So a later call starts at least as large as the size that last worked, never smaller than the user's setting.
- On SUCCESS (a converged solve, info==0 or info==2, after the retry loop exits ~line 334/340): store fw%last_successful_m0 = M0 (the M0 that actually succeeded). Only when present(fw).
- KEEP the existing info=3 retry loop (~lines 236-338) as the BACKSTOP: if a later k-point needs an EVEN larger subspace, the retry still doubles M0 (so no missing eigenvalues — user story 15). Do not weaken it.

SOUNDNESS (verify & document): the energy window [emin,emax] is fixed across a k-sweep, so the eigenvalue count (hence needed subspace) is ~constant; seeding from the last winning M0 avoids re-paying the retry tax. The retry loop remains for the rare larger-subspace point.

BIT-IDENTICAL ANALYSIS: This is a PERFORMANCE/caching change. The CONVERGED eigenvalues FEAST returns for a given (matrix, window) are independent of M0 as long as M0 is large enough (M0 only bounds the subspace; FEAST converges the same eigenpairs). Seeding a larger M0 does not change WHICH eigenvalues are found, only how fast. Therefore no computed eigenvalue/eigenvector changes on any config that already converged. The only behavioral effect: configs that hit info=3 today and abort may now converge (strictly an improvement) — but verify no golden config is in that aborting state today (if it were, it wouldn't be a passing golden test). Document this.

ACCEPTANCE CRITERIA:
1. Successful M0 stored on the feast_workspace and used as the starting size on the next call (capped at N, >= nev+1).
2. Retry loop still recovers if a later k-point requires a larger subspace (no missing eigenvalues).
3. (Central ctest) Existing QW FEAST sweep test (verification_qw_sparse_solver) stays green.

OWNED FILES (edit ONLY this): src/math/eigensolver.f90.
`

// ---------------------------------------------------------------------------
// SLICE 05 — SC-loop result clamp. Parallel.
// ---------------------------------------------------------------------------
const SLICE05 = PREAMBLE + `
=== YOUR SLICE: 05 — Self-consistent loop result clamp ===
Parent issue: .scratch/eigensolver-hardening/issues/05-sc-loop-result-clamp.md  (READ IT). Covers PRD user story 16. Isolated to the SC-loop module.

GOAL: Clamp the SC-loop's eigensolve RESULT WRITE to min(found, requested) and zero-fill the tail, matching every other dispatch site (band-structure sweep main.f90 ~748/797, Landau, optics all clamp). Today it is the lone holdout assuming the solver returns exactly the requested band count.

TARGET FILE: src/physics/sc_loop.f90  (READ src/physics/AGENTS.md FIRST). The QW SC loop self_consistent_loop, the result write ~lines 219-223:
    M_out = sc_result%nev_found
    eig_kpar(:, k_idx) = sc_result%eigenvalues(1:M_out)
    eigv_kpar(:, :, k_idx) = sc_result%eigenvectors(:, 1:num_subbands)
The bug: M_out = nev_found (FOUND), but eigv_kpar is sized (N, num_subbands, nk) and the eigenvector write indexes 1:num_subbands (REQUESTED). If nev_found < num_subbands, the eigenvalue write into eig_kpar(:,k_idx) (sized num_subbands) leaves a garbage tail, AND if nev_found > num_subbands the eigenvector read 1:num_subbands is fine but the eigenvalue write eig_kpar(:,k_idx)=eigenvalues(1:M_out) could overflow if M_out>num_subbands (eig_kpar sized num_subbands — CHECK the dimension of eig_kpar; it is (num_subbands, nk) per the allocation ~line ~176 area — read it).

FIX (mirror the pattern at sc_loop.f90 ~738-742 for the wire path, which ALREADY clamps correctly with num_subbands_actual = min(nev_found, nev_sc) and zero-fills the tail):
- M_out = min(sc_result%nev_found, num_subbands)
- eig_kpar(1:M_out, k_idx) = sc_result%eigenvalues(1:M_out)
- eigv_kpar(:, 1:M_out, k_idx) = sc_result%eigenvectors(:, 1:M_out)
- zero-fill the tail: if M_out < num_subbands: eig_kpar(M_out+1:num_subbands, k_idx) = 0.0_dp ; eigv_kpar(:, M_out+1:num_subbands, k_idx) = cmplx(0,0,dp)
Confirm the exact dimension/shape of eig_kpar and eigv_kpar by reading their allocations before writing the slice, so the indices are correct.

DO NOT TOUCH the SC-loop CONVERGENCE logic (DIIS/mixing/iteration count/Fermi bisection) — out of scope (approval-gated). Only the result WRITE clamps.

BIT-IDENTICAL ANALYSIS: On every valid config, the solver returns EXACTLY num_subbands (nev_found == num_subbands), so M_out = num_subbands and the clamp + zero-fill tail are no-ops (tail is empty). Verify: sc_cfg%nev is set to num_subbands (read where sc_cfg is built, ~lines before 200), and solve_dense in DENSE+INDEX mode returns exactly the requested nev. Therefore no computed SC potential / charge density / Fermi level changes. Document this. (The clamp is purely defensive for a future config where the solver might return fewer — user story 16.)

ACCEPTANCE CRITERIA:
1. SC-loop result write clamped to min(found, requested); tail zero-filled.
2. (Central ctest) SC regression tests stay green (clamp is a no-op on valid configs where count matches).

OWNED FILES (edit ONLY this): src/physics/sc_loop.f90.
`

// ---------------------------------------------------------------------------
// SLICE 06 — Truncation-warning integration test. Parallel.
// ---------------------------------------------------------------------------
const SLICE06 = PREAMBLE + `
=== YOUR SLICE: 06 — Truncation-warning integration test ===
Parent issue: .scratch/eigensolver-hardening/issues/06-truncation-warning-test.md  (READ IT). Covers PRD user stories 17-18. TEST-ONLY (the truncation warning code already exists).

GOAL: Close the testing gap from round 1. The QW FEAST truncation warning (printed when FEAST returns more eigenvalues than the retained band window) is implemented but untested. Add an integration check that runs a wide-window + narrow-band-range QW FEAST config and greps stdout for the warning.

THE WARNING TEXT (src/apps/main.f90 ~lines 749-753):
    print '(A,I0,A,I0,A)', '  Warning: FEAST returned ', result_bs%nev_found, ' eigenvalues at k-point ', k, '; only the lowest will be kept (widen bands or narrow energy window).'
So grep for the stable substring: "only the lowest will be kept".

TARGET FILE: tests/integration/verify_qw_sparse_solver.py  (READ tests/integration/AGENTS.md FIRST). This script already runs QW FEAST configs and asserts fast-vs-dense eigenvalue equivalence. It uses star_helpers.run_exe(BUILD_DIR, "bandStructure", config, workdir, timeout=...) and parses output. (It is already registered in tests/CMakeLists.txt ~line 903 — NO CMake edit needed.)

ADD a new check function (follow the existing verify_qw_* helper pattern in the file):
- Create a QW FEAST config TOML (per AGENTS.md anti-pattern, configs live in tests/regression/configs/ — create a NEW dedicated config file there, e.g. tests/regression/configs/qw_feast_truncation_warning.toml; do NOT create configs at runtime in a shell test). The config must trigger truncation: an ENERGY window [emin,emax] WIDER than the requested band range so FEAST finds MORE eigenvalues than (iuu-il+1). I.e. num_cb + num_vb (requested) small, but [solver] emin/emax spanning many subbands. Use confinement="qw", [solver] method="FEAST", explicit emin/emax, modest fd_step (e.g. 60-100) and nsteps (e.g. 5) to keep it fast. Model it on an existing QW config (look at tests/regression/configs/qw_*.toml for the required sections: [wave_vector], [bands], [[material]] with z_min/z_max, [solver]).
- Run bandStructure via run_exe, capture stdout (run_exe returns (rc, output_dir); you may need to capture stdout separately — check how the script currently captures program stdout; if run_exe doesn't surface stdout, run the exe directly via subprocess capturing stdout text, consistent with the script's existing style).
- PASS if "only the lowest will be kept" appears in stdout; FAIL otherwise. Print a clear PASS/FAIL line. Add a '# COVERAGE:' annotation (observable=truncation_warning geometry=qw material=GaAs) per AGENTS.md.
- Wire the new check into the script's main flow / exit code so ctest picks up failures (the script should exit nonzero if the check fails, matching how its other checks aggregate).

THE TEST MUST FAIL IF THE WARNING IS REMOVED (acceptance criterion 2): because it greps for the exact substring, deleting the print statement from main.f90 makes the grep fail. (You don't need to prove this by editing main.f90 — just ensure the grep is strict on the substring.)

NO source changes — test-only. The truncation code at main.f90:749-753 already exists.

ACCEPTANCE CRITERIA:
1. New check runs the wide-window + narrow-band-range config and greps stdout for "only the lowest will be kept".
2. Check fails if the warning text is removed (strict substring grep).
3. Check follows the existing verify_qw_* helper pattern in the script; new config in tests/regression/configs/.

OWNED FILES (edit ONLY these): tests/integration/verify_qw_sparse_solver.py, and the new config tests/regression/configs/qw_feast_truncation_warning.toml.
`

// ---------------------------------------------------------------------------
// ORCHESTRATION: sequenced 01->02 chain, concurrent with parallel 03/04/05/06.
// All edit-only (no builds). Disjoint file ownership => safe concurrency;
// 01->02 serialized because they share defs.f90 / main.f90 / simulation_setup.f90.
// ---------------------------------------------------------------------------
phase('Sequenced')
const r1 = await agent(SLICE01, { label: 'slice-01-resolver', phase: 'Sequenced', schema: SCHEMA })

phase('Parallel')
const five = await parallel([
  // The sequenced chain: 02 must wait for 01 (shared files). Runs concurrently
  // with 03/04/05/06 because its files are disjoint from theirs.
  async () => {
    const r2 = await agent(SLICE02, { label: 'slice-02-bulk', phase: 'Parallel', schema: SCHEMA })
    return { kind: 'chain', r2 }
  },
  () => agent(SLICE03, { label: 'slice-03-fastpath', phase: 'Parallel', schema: SCHEMA }).then(r => ({ kind: 'solo', r })),
  () => agent(SLICE04, { label: 'slice-04-m0', phase: 'Parallel', schema: SCHEMA }).then(r => ({ kind: 'solo', r })),
  () => agent(SLICE05, { label: 'slice-05-sc-clamp', phase: 'Parallel', schema: SCHEMA }).then(r => ({ kind: 'solo', r })),
  () => agent(SLICE06, { label: 'slice-06-trunc-test', phase: 'Parallel', schema: SCHEMA }).then(r => ({ kind: 'solo', r })),
])

// five[0] is the chain result (r2); five[1..4] are 03..06.
const r2  = five[0].r2
const r3  = five[1].r
const r4  = five[2].r
const r5  = five[3].r
const r6  = five[4].r

return {
  slices: { '01': r1, '02': r2, '03': r3, '04': r4, '05': r5, '06': r6 },
  all_files_changed: [r1, r2, r3, r4, r5, r6].filter(Boolean).flatMap(r => r.files_changed || []),
}
