export const meta = {
  name: 'eigensolver-hardening-review',
  description: '3 parallel read-only reviewers: spec compliance, code quality, invariants',
  phases: [
    { title: 'Review', detail: 'spec / quality / invariants in parallel' },
    { title: 'Verify', detail: 'adversarially re-check each high-severity finding' },
  ],
}

const CONTEXT = `
Review the "Eigensolver Standardization — Hardening (Review Fixes Round 2)" changes in repo /data/8bandkp-fdm (cwd = repo root).

WHAT CHANGED (six slices, all on branch refactor/eigensolver-standardization):
- 01: pure subroutine resolve_solver_defaults(confinement,method,mode,out_method,out_mode) added to src/core/defs.f90; the four confinement dispatch blocks in src/core/simulation_setup.f90 (bulk/QW/wire/Landau) + the band-structure sweep block in src/apps/main.f90 route through it; Landau no longer forces INDEX; table-driven unit test in tests/unit/test_eigensolver.pf.
- 02: bulk eigenpair storage sized to 8 (src/apps/main.f90); validate() tightened so bulk band count must == 8 (src/core/defs.f90); new rejection case in tests/integration/test_validate_rejects_bad_configs.sh.
- 03: QW CSR fast-path value-only branch deleted; single rebuild path kept (src/physics/hamiltonian_qw.f90).
- 04: FEAST persists the winning M0 on the workspace and seeds the next call (src/math/eigensolver.f90).
- 05: SC-loop result write clamped to min(found,requested) with zero-filled tail (src/physics/sc_loop.f90).
- 06: truncation-warning integration check (tests/integration/verify_qw_sparse_solver.py + new config tests/regression/configs/qw_feast_truncation_warning.toml).

INVARIABLE RULE that MUST hold: ZERO physics-value change — the golden regression suite is bit-identical. FEAST+INDEX must be unreachable via AUTO.

Use these read-only commands to gather evidence: 'git diff --stat', 'git diff -- <file>', 'git log --oneline -10'. Read CLAUDE.md, src/physics/AGENTS.md, tests/integration/AGENTS.md, .scratch/eigensolver-hardening/PRD.md and issues/01-06 as needed. Do NOT edit any files.
`

const FINDING_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['dimension', 'summary', 'findings', 'verdict'],
  properties: {
    dimension: { type: 'string' },
    summary: { type: 'string' },
    findings: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['severity', 'title', 'location', 'detail', 'actionable_fix'],
        properties: {
          severity: { type: 'string', enum: ['critical', 'high', 'medium', 'low', 'nit'] },
          title: { type: 'string' },
          location: { type: 'string', description: 'file:line or test name' },
          detail: { type: 'string', description: 'what is wrong / missing, with evidence' },
          actionable_fix: { type: 'string', description: 'concrete suggested change, or "none — informational"' },
        },
      },
    },
    verdict: { type: 'string', enum: ['pass', 'pass-with-notes', 'fail'], description: 'overall for this dimension' },
  },
}

const VERDICT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['title', 'confirmed_real', 'reasoning'],
  properties: {
    title: { type: 'string' },
    confirmed_real: { type: 'boolean', description: 'true if the finding is a real defect that must be fixed before merge' },
    reasoning: { type: 'string', description: 'evidence-based reasoning, cite file:line' },
  },
}

const REVIEWERS = [
  {
    key: 'spec',
    prompt: CONTEXT + `
YOUR DIMENSION: SPEC COMPLIANCE.
Build a checklist from EVERY user story in .scratch/eigensolver-hardening/PRD.md (there are 21, grouped) and EVERY acceptance-criterion checkbox in issues/01-06. For each, determine whether the diff actually delivers it, citing the file:line / test that proves it.
Pay special attention to:
- User stories 5-8 (FEAST+AUTO works on QW/Landau; Landau honors explicit mode) — does resolve_solver_defaults make these true? Is the no-invalid-combo invariant unit-tested?
- User story 1-4 (bulk always 8 bands; narrow range discarded safely; invalid count rejected up front; window is a pure display filter for bulk).
- User story 14-15 (M0 persists; retry still recovers larger subspace).
- User story 17-18 (truncation warning tested; test fails if warning removed).
- User story 19-20 (bit-identical golden).
Flag any acceptance criterion with NO corresponding diff as a HIGH finding. Flag any PRD user story not delivered as HIGH/CRITICAL.`,
  },
  {
    key: 'quality',
    prompt: CONTEXT + `
YOUR DIMENSION: CODE QUALITY & CONVENTIONS.
Check the diff against CLAUDE.md "Code Conventions" and the two AGENTS.md files (src/physics/AGENTS.md, tests/integration/AGENTS.md). Specifically:
- F2018 compliance: no GNU extensions, generic intrinsics (sqrt not dsqrt), error stop with messages (no bare 'stop 1'), private-default + explicit public exports for any new public symbol, elemental pure for scalar helpers, contiguous on hot-path assumed-shape args, no goto.
- defs.f90 hard rule: only a pure function/subroutine + validate() check added; NO derived-type changes; defs.f90 did not gain a 'use' of a downstream module.
- New test follows pFUnit patterns (@test, @assertEqual) and the rejection-script / verify_qw_*.py patterns; COVERAGE annotations present where required; no runtime config creation in shell tests (configs live in tests/regression/configs/).
- File/function length guidelines (~300/file, ~50/function) not egregiously blown.
- No dead code left behind (e.g., if slice 03 made csr_set_values_from_coo unused, note it; an unused helper is a low-severity note, not a blocker, unless it indicates an incomplete edit).
- Comments updated to match the new single-path / persisted-M0 / clamped logic (no stale comments describing removed branches).
Flag convention violations by severity. A genuine F2018 violation or a missed error stop is HIGH; style nits are low/nit.`,
  },
  {
    key: 'invariants',
    prompt: CONTEXT + `
YOUR DIMENSION: INVARIANTS & CORRECTNESS (the highest-stakes dimension).
Independently verify, from the source diff (not from agent claims):
1. BIT-IDENTICAL GOLDEN: For each slice, trace whether ANY supported configuration's computed output could change.
   - Slice 01: confirm the only resolution changes are FEAST+AUTO combos, that wire FEAST+AUTO is unchanged (ENERGY->ENERGY), and that no bulk/qw/landau config in tests/regression/configs/ uses method="FEAST" (so qw/landau FEAST+AUTO was aborting, bulk FEAST absent). If a golden config COULD change, that is CRITICAL.
   - Slice 02: confirm bulk storage sizing to 8 + the !=8 rejection cannot change any existing bulk config (all must have evnum 8). Check tests/regression/configs/ bulk configs.
   - Slice 03: confirm the deleted value-only branch was truly unreachable (the QW sweep frees HT_csr each iteration => nnz==0 => always the rebuild branch). If any current caller took the value-only branch, deleting it changes behavior => CRITICAL.
   - Slice 04: confirm persisted M0 changes only convergence SPEED, not the converged eigenvalues (M0 bounds subspace size; FEAST returns the same eigenpairs given a large-enough M0). Confirm the retry backstop is intact.
   - Slice 05: confirm min(found,requested) + zero-fill is a no-op when found==requested (the valid case).
2. FEAST+INDEX UNREACHABLE VIA AUTO: trace resolve_solver_defaults for every confinement with mode AUTO — none may yield (FEAST,INDEX). Confirm explicit FEAST+INDEX is still rejected by validate() check I15.
3. OUT-OF-BOUNDS: slice 02's sizing and slice 05's clamp must actually be correct (right array dimensions, no off-by-one in the zero-fill tail).
4. SCOPE DISCIPLINE: confirm nothing out-of-scope was touched (material params, basis ordering, Bir-Pikus signs, Zeeman table, FD stencils, Poisson, SC convergence logic, wire value-only path). 'git diff --stat' should show ONLY the owned files.
Any risk to bit-identical golden or to the AUTO invariant is CRITICAL/HIGH and must block.`,
  },
]

phase('Review')
const reviews = await parallel(
  REVIEWERS.map(r => () => agent(r.prompt, { label: 'review:' + r.key, phase: 'Review', schema: FINDING_SCHEMA }))
)

// Adversarially verify every high/critical finding before surfacing it —
// kills plausible-but-wrong reviewer claims.
phase('Verify')
const toVerify = reviews.filter(Boolean).flatMap(r => (r.findings || [])
  .filter(f => f.severity === 'critical' || f.severity === 'high')
  .map(f => ({ ...f, dimension: r.dimension })))

const verified = await parallel(toVerify.map(f => () =>
  agent(
    CONTEXT + `\nA reviewer (dimension: ${f.dimension}) raised this ${f.severity} finding:\nTITLE: ${f.title}\nLOCATION: ${f.location}\nDETAIL: ${f.detail}\nPROPOSED FIX: ${f.actionable_fix}\n\nAdversarially verify whether this is a REAL defect that must block merge. Read the actual code at the location and the surrounding context. Default to confirmed_real=false unless you can cite concrete evidence (file:line) that the defect exists and affects behavior or violates a stated invariant. If it is real, say so plainly; if the reviewer is wrong, explain why.`,
    { label: 'verify:' + (f.title || '').slice(0, 30), phase: 'Verify', schema: VERDICT_SCHEMA }
  ).then(v => ({ ...v, original: f })).catch(() => null)
))

const confirmed = verified.filter(Boolean).filter(v => v.confirmed_real)

return {
  dimensions: reviews.filter(Boolean).map(r => ({ dimension: r.dimension, verdict: r.verdict, finding_count: (r.findings || []).length })),
  confirmed_blockers: confirmed.map(v => ({ title: v.title, location: v.original.location, detail: v.original.detail, reasoning: v.reasoning, proposed_fix: v.original.actionable_fix })),
  all_findings: reviews.filter(Boolean),
}
