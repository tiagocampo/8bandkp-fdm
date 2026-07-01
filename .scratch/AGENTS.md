# .scratch/ — Active Planning Artifacts

## Convention

- **Active** planning artifacts (PRDs, issue briefs, workflow scripts) for an in-flight
  PR/issue live in **`.scratch/<name>/`** at the repo root.
- **Completed** PRs' planning artifacts are moved to **`.scratch/archive/<name>/`** to
  preserve the record of what was decided without cluttering the active namespace.

## Why `git mv` (not `mv` + `git add`)

Use `git mv <src> .scratch/archive/<dst>/` to preserve file history. A plain
`mv` followed by `git add` / `git rm` rewrites history as a delete + add and
loses rename detection — `git log --follow` no longer works on those files.

## Precedent

This convention is already established in two prior cleanup commits:

- `6de7c67 chore(scratch): archive pr27-review-fixes to .scratch/archive/`
- `fa84912 chore(scratch): archive 14 completed-PRD scratch directories to .scratch/archive/`

Both moved multiple subtrees with `git mv` so the archive entries keep a
contiguous history back to their original creation.

## Tracking status (intentionally tracked)

`.scratch/` and `.scratch/archive/` are **tracked** in git. This is intentional
and is the established pattern (see the two precedent commits above plus the
116 files currently tracked under `.scratch/`).

This is **different** from `tmp/`, which is gitignored. Use `tmp/` for ephemeral
scratch you don't want in the repo; use `.scratch/` for planning artifacts that
need to travel with the PR/issue they describe.

## Recovery: stale tracked build directories

If a build directory such as `build/`, `build-nolto/`, or `build-clang/` ever
ends up tracked in git (because the directory was created **before** the
`build*/` `.gitignore` rule was added), untrack it without deleting it from disk:

```bash
git rm -r --cached build-nolto/   # or whichever stale build dir
# on-disk files are preserved; .gitignore blocks future re-tracking
rm -rf build-nolto/               # free the disk space only after the git rm
```

The CMake-generated `.o`/`.mod`/test-driver files have no business being in
source control — untrack them and let `cmake` regenerate as needed.
