# 04 — Wire setup type + strain-omission fix

## Parent / Source

PRD: `.scratch/architecture-deepening/PRD.md` (Thread 2, U-C5 / U-C1-remainder).
ADRs 0001, 0005.

## What to build

The 2D-wire init (profile, sparse k.p terms, wire workspace, COO cache,
eigensolver) and cleanup (per-element free, workspace/cache free) is copy-pasted
across the topology wire subroutines and the spectral-function path. Those copies
also **skip the strain and electric-field setup the QW path gets** — so strained
wire topology / BdG calculations silently ignore strain. This is a **latent
physics defect**, not a style issue.

Encapsulate wire init + cleanup in **one type** whose init runs the **same strain
and electric-field setup as the QW path** (fixing the bug by construction), with
a variant that accepts pre-computed profile/terms for sink-style callers that
receive that data as arguments. Route the topology wire subroutines and the
spectral-function path through it.

## Acceptance criteria

- [ ] One wire setup type owns init + idempotent free for all wire callers; the
      copy-pasted boilerplate is gone.
- [ ] Sink-style callers accept pre-computed profile/terms instead of
      re-initializing.
- [ ] **Regression test: a strained wire topology/BdG config whose spectrum
      shifts once strain is applied** — proves the strain-omission fix.
- [ ] A documented solution entry is added once the strain fix lands.
- [ ] Golden regressions unchanged where strain was already applied.

## Blocked by

- #02 (solver-config derivation) — the wire setup builds its eigensolver config
  through the derivation seam rather than re-assembling it.
