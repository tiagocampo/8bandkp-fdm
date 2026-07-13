# Issue 3: Consolidate app-semantic validation into validate_semantic()

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21
**GitHub issue**: #22

## What to build

Add 10 app-specific checks to the existing `validate_semantic(cfg, app_name)` subroutine in `defs.f90`. These checks require knowledge of which executable is running. All use `error stop` with contextual messages.

Checks to add (S1–S10 from ADR 0002):

| ID | App | Check | Error message pattern |
|----|-----|-------|-----------------------|
| S1 | gfactor | `bandIdx` in range [1, num_cb] for wire | `bandIdx (=X) out of range [1, Y] for wire gfactor` |
| S2 | topologicalAnalysis | confinement must be QW or wire for QSHE Z2 | `QSHE Z2 requires QW or wire confinement, got 'X'` |
| S3 | topologicalAnalysis | BdG mode requires `[bdg]` section enabled | `topology mode 'bdg' requires [bdg] section` |
| S4 | topologicalAnalysis | BdG confinement must be QW or wire | `BdG requires QW or wire confinement, got 'X'` |
| S5 | topologicalAnalysis | `spectral_eta > 0` | `spectral_eta must be positive, got X` |
| S6 | topologicalAnalysis | `spectral_nk >= 1` | `spectral_nk must be >= 1, got X` |
| S7 | topologicalAnalysis | `spectral_nE >= 1` | `spectral_nE must be >= 1, got X` |
| S8 | topologicalAnalysis | sweep_model matches confinement | `sweep_model 'X' requires confinement='Y', got 'Z'` |
| S9 | topologicalAnalysis | `conductance_method` is valid enum | `conductance_method 'X' not recognized (expected: landauer, kutzner)` |
| S10 | topologicalAnalysis | topology mode is valid enum | `topology mode 'X' not recognized (expected: qhe, qshe, bdg)` |

These checks are added to the existing `select case(app_name)` dispatch in `validate_semantic`. The existing checks (gfactor↔k0, optics↔enabled, topology↔enabled+mode) remain unchanged.

Key implementation notes:
- S1: Only fires when `app_name == 'gfactor'` and `cfg%confinement == 'wire'`. Check `cfg%band_idx >= 1 .and. cfg%band_idx <= cfg%bands%num_cb`.
- S2: Only fires for QSHE Z2 mode. Check confinement is QW or wire.
- S3–S4: Only fire when `cfg%topo%mode == 'bdg'`.
- S5–S7: Only fire when spectral mode is active (check `cfg%topo%mode` or spectral params present).
- S8: Check each sweep_model value against its required confinement.
- S9: Enumerate valid conductance methods from existing code.
- S10: Enumerate valid topology modes from existing code.

## Acceptance criteria

- [ ] All 10 checks added to `validate_semantic` in the appropriate `select case` branches
- [ ] Each check uses `error stop` with contextual message
- [ ] Existing 5 checks in `validate_semantic()` remain unchanged
- [ ] Invalid app-specific configs fail with clear messages
- [ ] All 34 unit tests pass

## Blocked by

None — can start immediately.
