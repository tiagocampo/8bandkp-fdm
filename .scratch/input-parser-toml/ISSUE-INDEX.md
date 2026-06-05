# Issue Tracker: TOML Input Parser Refactor

Parent PRD: `.scratch/input-parser-toml/PRD.md`

## Dependency Graph

```
001 (toml-f build) ─────────────────────────────────────────────────────┐
                                                                         │
002 (bulk tracer bullet) ◄──────────────────────────────────────────────┘
 │
 ├── 003 (QW geometry)
 ├── 004 (wire geometry)
 ├── 005 (Landau geometry)
 ├── 006 (external field + B-field + strain + FEAST + g-factor)
 ├── 007 (SC + doping)
 ├── 008 (topology + BdG)
 ├── 009 (optics + exciton + scattering) ──────────────────────┐
 └── 010 (config converter script)                              │
      │                                                         │
      011 (test infrastructure)                                 012 (docs overhaul) ◄─┘
      │                                                         │
      └─────────────── 013 (full verification + cleanup) ◄──────┘
```

## Issues

| # | Title | Type | Blocked by | User Stories |
|---|-------|------|-----------|-------------|
| 001 | toml-f build integration | AFK | None | 28 |
| 002 | Bulk mode tracer bullet — type restructure + parser + consumers | AFK | 001 | 1-4, 6, 10, 16-18, 20-28 |
| 003 | QW geometry parser | AFK | 002 | 5, 23 |
| 004 | Wire geometry parser | AFK | 002 | 7 |
| 005 | Landau geometry parser | AFK | 002 | 8 |
| 006 | External field + B-field + strain + FEAST + g-factor parsers | AFK | 002 | 9, 10, 14, 16 |
| 007 | SC + doping parser | AFK | 002 | 11, 13 |
| 008 | Topology + BdG parser | AFK | 002 | 12 |
| 009 | Optics + exciton + scattering parsers | AFK | 002 | 15 |
| 010 | Config converter script | AFK | 002 | 27 |
| 011 | Test infrastructure update | AFK | 010 | 26, 27 |
| 012 | Documentation overhaul | AFK | 009 | 29-33 |
| 013 | Full verification + cleanup | AFK | 011, 012 | 27 |

## Parallelism after Issue 002

Issues 003-010 can all be worked in parallel after Issue 002 lands. The critical path is:

**001 → 002 → 003-010 (parallel) → 011 + 012 → 013**
