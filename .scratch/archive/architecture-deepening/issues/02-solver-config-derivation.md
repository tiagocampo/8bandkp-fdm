# 02 — Solver-config derivation

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 1, U-C). ADRs 0004, 0005.

## What to build

The eigensolver config — method/mode from the method-aware `AUTO` resolver, then
`nev`/`il`/`iu`, subspace size, validation, and construction — is hand-assembled
in each confinement dispatch block and re-shorthand across the apps. Replace that
with **one derivation entry point** that all four confinement cases and the apps
cross. It **calls** the existing `AUTO` resolver (ADR 0004); it does not replace
it — two layers, not two functions racing to own `AUTO`.

Also rename the single-k serial solver used only by the g-factor Γ computation so
its name reflects that role (it is not a general k-sweep entry point — its only
caller is the g-factor app), and update that caller.

## Acceptance criteria

- [ ] One derivation entry point is the sole path to a built, validated
      eigensolver config for all four geometries and the apps.
- [ ] The derivation calls the existing method-aware `AUTO` resolver; `AUTO`
      semantics and its table-driven unit test are unchanged.
- [ ] New unit test drives `AUTO` resolution per geometry × method × mode
      directly (including the no-invalid-combination invariant) without running
      a full executable — the case that is currently untestable through init.
- [ ] The g-factor Γ solver is renamed; its single caller updated; behavior
      unchanged.
- [ ] Golden regression outputs unchanged.

## Blocked by

None — can start immediately.
