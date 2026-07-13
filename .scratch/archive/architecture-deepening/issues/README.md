# Architecture Deepening — Issue Index

Parent/source: `../PRD.md`. Branch: `refactor/architecture-deepening`.
All issues are **AFK-ready** (decisions locked in ADRs 0001–0005 + the
2026-06-13/14 grill session; no in-flight design questions). Created in
dependency order so blockers reference lower numbers.

| # | Title | Type | Blocked by | User stories |
|---|---|---|---|---|
| 01 | [Eigensolver dispatch cleanup](01-eigensolver-dispatch-cleanup.md) | AFK | — | 20, 21, 22 |
| 02 | [Solver-config derivation](02-solver-config-derivation.md) | AFK | — | 9, 10, 28 |
| 03 | [Energy-window authority (#10 gate)](03-energy-window-authority.md) | AFK | #02 | 2, 3, 11, 12, 13 |
| 04 | [Wire setup type + strain-omission fix](04-wire-setup-strain-fix.md) | AFK | #02 | 14, 15, 16 |
| 05 | [Confinement-assembly unit tests (C7 gate)](05-confinement-tests.md) | AFK | — | 17 |
| 06 | [Block-formula descriptor + strain split](06-block-formula-strain-split.md) | AFK | #05 | 18, 19, 30 |
| 07 | [Output writers + Simpson absorption](07-output-writers-simpson.md) | AFK | — | 23, 24 |
| 08 | [FEAST-enable QW g-factor](08-feast-qw-gfactor.md) | AFK | #03 | 8 |
| 09 | [FEAST-enable QW optics](09-feast-qw-optics.md) | AFK | #02, #03 | 1, 4, 6 |
| 10 | [FEAST-enable Landau](10-feast-landau.md) | AFK | #03 | 5 |

**Critical paths:** #02 → #03 → {#08, #09, #10} (the FEAST-parity flip-ons) and
#05 → #06 (the physics-critical block-table track). Starters with no blockers:
**#01, #02, #05, #07**.
