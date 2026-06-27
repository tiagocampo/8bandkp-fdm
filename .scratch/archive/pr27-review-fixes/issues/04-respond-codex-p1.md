# Issue 04: Respond to Codex P1 comments on PR #27

**Priority:** Important (communication)
**No code changes** — GitHub PR comments only

## Task

Post inline comments on the 2 Codex P1 findings explaining why they are false positives:

- [ ] **Comment on `src/io/input_parser.f90:703` (b_field copy):**
  Explain that lines 699-700 (`cfg%bdg%B_vec = cfg%b_field%components` and `cfg%bdg%g_factor = cfg%b_field%g_factor`) execute BEFORE the early return at line 703 (`if (.not. associated(bdg_tbl)) return`). The data flow is correct for Landau configs that have `[b_field]` but no `[bdg]`.

- [ ] **Comment on `src/io/input_parser.f90:871` (strainSubstrate):**
  Explain that `parse_strain` writes to `cfg%params(1)%strainSubstrate` at line 870. The subsequent `paramDatabase()` call replaces parameter structs but does NOT touch `strainSubstrate` — the field survives. Fragile coupling (noted as suggestion), but not a bug.

## Verification

- [ ] Comments posted on PR #27 with 👍 reactions
