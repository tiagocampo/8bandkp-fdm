# Issue 001: toml-f build integration

**Type:** AFK
**Blocked by:** None — can start immediately
**User stories:** US 28

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add the `toml-f` library as a git submodule and integrate it into the CMake build system. After this issue, the project should compile and link with `toml-f` available as a dependency, but no code changes to the parser or types yet.

Specifically:
- Add `toml-f` v0.4.2 as a git submodule in `subprojects/toml-f/`
- Update the root `CMakeLists.txt` to `add_subdirectory(subprojects/toml-f)`
- Link `toml-f-lib` to the `input_parser` module (or the executables that use it)
- Verify the build compiles with `-std=f2018` — toml-f is pure Fortran 2008 which is a subset
- Verify existing tests still pass (no behavior changes yet)

## Acceptance criteria

- [ ] `git submodule status` shows `subprojects/toml-f` checked out
- [ ] `cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl` succeeds
- [ ] `cmake --build build` succeeds (all four executables built)
- [ ] `ctest --test-dir build -j4` passes with same results as before
- [ ] No changes to any source files except `CMakeLists.txt`

## Blocked by

None — can start immediately.
