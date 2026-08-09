# T5 — TDD-red strict S1×S2 agreement test

> Child of [.scratch/bdg-u13-bloch-pfaffian/map.md](../map.md). Build branch;
> the red half of T4's green. Materialized from the map's T5 entry.

Type: task
Status: pending
Blocked by: 02

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
