# T5 — TDD-red strict S1×S2 agreement test

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch;
> the red half of T4's green. Materialized from the map's T5 entry.

Type: task
Status: resolved
Owner: this-session (2026-08-09 — claimed after T2 closure at a864f98)
Blocked by: 02 (closed 2026-08-09 — T2 at 449c7fb)

## Answer (2026-08-09)

Extended `tests/unit/test_bdg_pfaffian_witness_product_csr.pf` with 3 new
`@test` subroutines covering the remaining L3 branches (cases ii, iii, v).
The 5-branch L3 mapping in `bdg_observables.f90:361-389` is now
end-to-end pinned by 6 `@test` subroutines (3 from T2 + 3 from T5):

  - case (i) — `test_product_seam_strict_agreement_sign_split` (T2) — z2=-1, reason=0
  - case (i) — `test_product_seam_strict_topological_sign`     (T2) — z2=-1, reason=0
  - case (ii) — `test_product_seam_strict_agreement_both_zero`  (T5) — z2=0,  reason=0
  - case (iii) — `test_product_seam_strict_sign_split`          (T5) — z2=0,  reason=1
  - case (iv) — `test_product_seam_s2_fest_saturation_branch`   (T2) — z2=0,  reason=2
  - case (v) — `test_product_seam_s1_closure_s2_misses`         (T5) — z2=0,  reason=3

### Design refinement vs T2's forward-pointer

T2's ticket body suggested 4×4 (n_odd=1) Kitaev fixtures with controlled
det for cases (iii) and (v). Resolved in T5 to 16×16 stacks exploiting
two cleaner mechanisms:

  - **n_k=1 stacks** trigger `kitaev_majorana_number`'s `n_k < 2` guard
    (`src/math/pfaffian.f90:132`) → S1=0 immediately. Used for cases
    (ii) and (v).
  - **Band-7/8 coupling sign** controls S2 sign on a 16×16 slice-1 CSR
    (empirically: +0.3 → s2=-1, -0.3 → s2=+1). Used for case (iii).

No 4×4 detour needed; all 5 branches pinned on the existing 16×16 stack
shape T2 already validated.

### TDD red/green cycle (executed)

RED: stashed T5 test changes + reverted `src/physics/bdg_observables.f90`
to pre-T2 baseline `1a14489`. Build failed:
```
Error: Symbol 'eval_bdg_pfaffian_witness_product_csr' referenced at (1)
not found in module 'bdg_observables'
```
GREEN: restored T2 seam at HEAD + restored T5 tests. `ctest -V` reports
`(6 tests) OK`. Full cycle: red on pre-T2, green at HEAD.

### Fixture math notes (empirical, not derived)

T5 fixtures do not analytically derive S1; they trigger the polar
decomposition + slim Pfaffian empirically and assert the seam's
disagreement_reason mapping. The seam's `map_s1_s2_to_z2` is the SSOT
for the S1×S2 → z2 mapping, and the tests pin it end-to-end.

### Out of scope (deferred to T4)

- B-sweep full-fixture strict-agreement assertion (the gate-flip ticket).
- `tests/integration/test_wire_bdg_topological_phase.sh` extension —
  T5 is unit-only per the seam/gate split the map §State records.
- Updating the file header's `RED on main` → `GREEN at HEAD` line
  (T2 still says "RED on main"; factual at the time T2 wrote it, but
  no longer accurate).

T5 unblocks T4 (the destination row of the chart).

## Question

## Question

The **red** test, written first under TDD before any of T1–T4 machinery exists,
that pins the destination: it asserts **strict S1×S2 agreement at μ≈0.6601
across B** and it **must fail on `main`** (red), then turn green once T1–T4 land.
It straddles T2 (the seam) — claim it right after T2's wiring so the red probe
exists before T3's dispatch and T4's gate flip.

Contract a session should nail down (drawing on the map's §State and the U10
memory chain `project_bdg_u10_t3_execution` / `project_bdg_u10_t4_execution`,
which established a TDD-red tin `test_wire_bdg_topological_phase.sh` for the native
`{-1,+1,0}` convention):

- **What it asserts.** For the canonical wire-BdG fixture
  (`tests/.../wire_inas_gaas_bdg_topological_phase.toml`, μ≈0.6601 eV over a B
  grid), once the strict seam (T2) is wired through T1's stack, S1 (the
  k-product Majorana number) and S2 (the slim projected Pfaffian) **agree** —
  equal signs in the open/closed region, and the product `s1 × s2` reproduces the
  gap-closure pattern `bcrit_2d` already encodes (the gap-min strategy at
  `lecture_13_topological.py:178-196`). Native `{-1, +1, 0}` per the U10 T3 wire
  path convention.
- **Why it fails on `main` (red).** On `main`, the only Pfaffian seam is
  the slim single-point S2 (`eval_bdg_pfaffian_witness_csr`); S1 is never
  computed because there is no Bloch stack. The assertion "S1 and S2 agree"
  has no S1 to read → the test errors/fails until T1+T2 supply it. That's the
  intended red.
- **Where it lives.** A pFUnit `@test` in `tests/unit/test_wire_pfaffian_witness.pf`
  (or a sibling `.pf`) for the unit-level seam probe, AND a regression/integration
  hook in `tests/integration/test_wire_bdg_topological_phase.sh` for the
  full-fixture strict-agreement gate — mirror the U10 T3 red-tin pattern. Mind
  the pFUnit gotchas in project memory: `@assertEqual`/`@assertTrue` must be
  **single-line** (no `&` continuation on `@`), and a `@test` hitting an
  `error stop` kills the whole ctest process — so avoid authoring a
  failing-input `@test` for the rejection branch; put rejection assertions in
  the `.sh` exit-code path.
- **COVERAGE annotation.** Include a `# COVERAGE: observable=majorana_number
  geometry=wire material=InAs tier=verification` line per `tests/integration/AGENTS.md`
  so the coverage matrix tracks the new strict observable.

Output: the red test, committed on `feat/bdg-u13-bloch-pfaffian` **before** T3
and T4. It turns green when T2 (its dependency) + T1 (T2's dep) are in; T4
flips the label the test doesn't directly check. T5 is thus blocked by T2 only
(the seam it probes), not by T4 (the label T4 flips). Read
`tests/integration/AGENTS.md` and `tests/unit/` conventions before authoring.
