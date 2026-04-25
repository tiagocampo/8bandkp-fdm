# Chapter 12: Extending the Code

A practical guide for developers who want to add materials, physics, output, or tests to the 8-band k.p finite-difference code. This chapter assumes you have read the earlier chapters and are comfortable navigating the source tree. The goal is to give you enough map to find your way without a companion.

## Module Dependency Graph

Every module lives under `src/` in one of four subdirectories. The dependency graph below shows which module `use`s which. Arrows point from a module toward the modules it depends on. Read it bottom-up: if you change a lower module, everything above it may need recompilation.

```
src/apps/
  main.f90 (program kpfdm)
  main_gfactor.f90 (program gfactor)
    |-> input_parser       (reads input.cfg, calls paramDatabase, init_grid)
    |-> hamiltonianConstructor  (builds H for each k-point)
    |-> outputFunctions    (writes eigenvalues, eigenfunctions)
    |-> sc_loop            (self-consistent SP, only main.f90)
    |-> strain_solver      (elastic + Pikus-Bir, both mains)
    |-> eigensolver        (FEAST wrapper, wire mode only)
    |-> sparse_matrices    (CSR type, Kronecker, COO assembly)
    |-> linalg             (zheevx wrapper, ILAENV, DLAMCH)

src/physics/
  hamiltonianConstructor.f90
    |-> finitedifferences  (FD stencil matrices)
    |-> sparse_matrices    (CSR operations for wire mode)
    |-> strain_solver      (Bir-Pikus strain blocks)
    |-> utils              (dense-to-sparse conversion)
  gfactor_functions.f90
    |-> hamiltonianConstructor
    |-> sparse_matrices    (CSR SpMV for QW/wire perturbation)
  poisson.f90
    |-> definitions        (types, constants only)
  charge_density.f90
    |-> utils
    |-> definitions
  sc_loop.f90
    |-> poisson
    |-> charge_density
    |-> hamiltonianConstructor
  strain_solver.f90
    |-> definitions        (types, constants only)

src/math/
  finitedifferences.f90
    |-> definitions
  sparse_matrices.f90       (CSR type, Kronecker, cut-cell ops)
    |-> definitions
  linalg.f90                (zheevx wrapper, ILAENV, DLAMCH)
    |-> definitions
  eigensolver.f90           (FEAST contour eigensolver)
    |-> definitions
    |-> sparse_matrices
  geometry.f90              (wire grid, cut-cells, immersed boundary)
    |-> definitions

src/io/
  outputFunctions.f90
    |-> utils
  input_parser.f90
    |-> definitions, parameters, geometry
    |-> hamiltonianConstructor (for confinementInitialization)
    |-> outputFunctions

src/core/
  defs.f90                   (kinds, constants, derived types -- no deps)
  parameters.f90             (material database)
    |-> definitions
  utils.f90                  (dense-to-sparse, Simpson integration)
    |-> definitions
```

The key takeaway: `defs.f90` is the root. Every other module depends on it, directly or transitively. The `simulation_config` derived type defined there is the single data structure that threads through the entire program. If you add a new field to the simulation, start in `defs.f90`.

## Data Flow

Understanding the data flow is the most important thing for extending the code. Everything passes through a single pipeline:

```
input.cfg
  |
  v
input_parser.f90 :: read_and_setup(cfg, profile, kpterms)
  |   reads label-value pairs
  |   calls paramDatabase() to fill cfg%params(:)
  |   calls init_grid_from_config() to build cfg%grid
  |   calls confinementInitialization() for QW mode
  |         or init_wire_from_config() for wire mode
  v
simulation_config cfg   <-- lives in defs.f90, holds everything
  |
  +-- cfg%params(:)        paramStruct array, one per layer/region
  +-- cfg%grid             spatial_grid (coordinates, material_id, cut-cells)
  +-- cfg%sc               sc_config (self-consistent parameters)
  +-- cfg%strain           strain_config
  +-- profile(:,:)         band edges at each grid point (EV, EV-DeltaSO, EC)
  +-- kpterms(:,:,:)       FD operators weighted by material params (QW mode)
  or
  +-- profile_2d(:,:)      band edges on 2D grid (wire mode)
  +-- kpterms_2d(:)        CSR matrices for 2D FD operators (wire mode)
  |
  v
main.f90: k-vector sweep
  for each k-point:
    hamiltonianConstructor :: ZB8bandBulk(HT, wv, params)
         or ZB8bandQW(HT, wv, profile, kpterms)
         or ZB8bandGeneralized(HT_csr, kz, profile_2d, kpterms_2d, cfg)
    |
    v
  solver: zheevx (dense, bulk/QW) or FEAST (sparse, wire)
    |
    v
  eigenvalues(:,:), eigenvectors(:,:,:)
    |
    v
  outputFunctions :: writeEigenvalues(), writeEigenfunctions()
    writes to output/eigenvalues.dat, output/eigenfunctions_k_*.dat
```

The self-consistent loop wraps the eigenvalue solve. In QW mode, `sc_loop.f90` calls `poisson.f90` and `charge_density.f90` between iterations, updating the `profile` array before the next Hamiltonian construction. The Hamiltonian assembly itself never changes -- only the band-edge profile shifts.

## Adding a New Material

All materials live in `src/core/parameters.f90` inside the `paramDatabase` subroutine. It is a large `select case` block, one case per material name string. The `paramStruct` type is defined in `defs.f90` and has these fields:

| Field | Units | Description |
|-------|-------|-------------|
| `meff` | m0 | Effective mass (informational, not used directly) |
| `gamma1`, `gamma2`, `gamma3` | dimensionless | Luttinger parameters |
| `P` | eV^(1/2) A^(-1) | Momentum matrix element (from EP: P = sqrt(EP*const)) |
| `A` | dimensionless | Conduction band kinetic coefficient (A = 1/meff) |
| `EP` | eV | Kane energy |
| `Eg` | eV | Band gap |
| `deltaSO` | eV | Spin-orbit splitting |
| `EV` | eV | Valence band edge (absolute, reference-dependent) |
| `EC` | eV | Conduction band edge (EC = EV + Eg) |
| `eps0` | dimensionless | Static dielectric constant |
| `C11`, `C12`, `C44` | GPa | Elastic constants |
| `a0` | Angstrom | Lattice constant |
| `ac` | eV | CB hydrostatic deformation potential |
| `av` | eV | VB hydrostatic deformation potential (positive sign convention) |
| `b_dp` | eV | Tetragonal shear deformation potential |
| `d_dp` | eV | Rhombohedral shear deformation potential |

Step-by-step to add a new material, say "InP":

1. Open `src/core/parameters.f90` and find the `select case(material(i))` block.

2. Add a new `case ("InP")` block before the `case default` clause. It must set every field listed above. For example:

```fortran
case ("InP")
  params(i)%meff    = 0.077_dp
  params(i)%EP      = 20.7_dp
  params(i)%Eg      = 1.424_dp
  params(i)%deltaSO = 0.11_dp
  params(i)%gamma1  = 5.33_dp
  params(i)%gamma2  = 1.57_dp
  params(i)%gamma3  = 2.11_dp
  params(i)%EV      = -0.94_dp
  params(i)%EC      = 0.484_dp
  params(i)%eps0    = 12.61_dp
  params(i)%C11   = 1011.0_dp
  params(i)%C12   = 561.0_dp
  params(i)%C44   = 456.0_dp
  params(i)%a0    = 5.8697_dp
  params(i)%ac    = -6.0_dp
  params(i)%av    = 1.7_dp
  params(i)%b_dp  = -2.0_dp
  params(i)%d_dp  = -4.5_dp
```

3. Compute `A` and `P` from `EP`:
   - `P = sqrt(EP * const)` where `const` is defined in `defs.f90` (3.809982... eV A^2).
   - `A = 1/meff` (the default; for W-variant materials, A follows a more complex formula involving EP, Eg, deltaSO, and the gamma parameters).
   These are used inside `confinementInitialization`, which reads `params(i)%P` and `params(i)%A` from the `paramStruct`. If you set them explicitly in the case block, make sure they are consistent with EP.

4. Set `EV` and `EC` with the correct band offsets relative to the other materials in your heterostructure. The code uses the convention `EC = EV + Eg`. For W-variant materials (suffix `W`), EV uses Winkler's InSb reference. For non-W materials, Vurgaftman's reference is used.

5. Verify your parameters against Vurgaftman (2001) or Winkler (2003). The project policy requires a published reference for every parameter set.

6. Test with a bulk calculation. Create an `input.cfg` with `confinement 0`, `materialN InP`, and check that the band gap, effective masses, and Luttinger parameters reproduce published values.

For alloy materials (e.g., AlGaAs with specific Al fraction), add a dedicated case with Vegard-interpolated parameters. The existing `Al20Ga80As` and `Al15Ga85As` entries show the pattern.

## Adding New Physics: Where to Hook In

The Hamiltonian construction is the primary hook point for new physics. There are three levels where you can intervene:

### Level 1: Before the Hamiltonian (modify the band-edge profile)

If your physics changes the on-site energies but not the kinetic terms, modify the `profile` array before the k-sweep. The existing examples are:

- **Electric field**: `externalFieldSetup_electricField` in `hamiltonianConstructor.f90` adds a linear potential `-E * z` to the profile.
- **Strain**: `compute_bir_pikus_blocks` in `strain_solver.f90` computes Bir-Pikus shifts from the strain tensor.
- **Self-consistent SP**: The SC loop updates `profile` iteratively via the Poisson potential.

The profile is a 2D array with shape `(Ngrid, 3)` where column 1 is EV, column 2 is EV minus DeltaSO (SO band edge), and column 3 is EC. In QW mode it has shape `(fdStep, 3)`. Modify it any time before the Hamiltonian is built for a given k-point.

### Level 2: Inside the Hamiltonian (modify k.p terms)

If you need to change the kinetic terms or add new couplings, you must edit `hamiltonianConstructor.f90`. The relevant entry points are:

- **Bulk**: `ZB8bandBulk` -- 8x8 matrix, all terms are explicit scalars. Add new off-diagonal elements directly.
- **QW**: `ZB8bandQW` -- 8N x 8N block matrix. The `kpterms(:,:,:)` array holds the precomputed position-dependent FD operators. Index 1-4 are diagonal (gamma1, gamma2, gamma3, P), 5 is A times d^2/dz^2, 6 is P times d/dz, 7 is (gamma1-2gamma2) times d^2/dz^2, 8 is (gamma1+2gamma2) times d^2/dz^2, 9 is gamma3 times d/dz, 10 is A (diagonal).
- **Wire**: `ZB8bandGeneralized` -- 8*Ngrid x 8*Ngrid sparse CSR. The `kpterms_2d(:)` array holds 17 CSR matrices (indices 1-17). Indices 1-4 diagonal, 5 A*Laplacian, 6 P*gradient, 7-8 Q/T kinetic terms, 9 gamma3*gradient (legacy), 10 A diagonal, 11 cross-derivative, 12-13 P*x/y gradients for g-factor, 14-15 gamma3*x/y gradients for S/SC terms and g-factor, 16 gamma2*(D2x-D2y) anisotropic Laplacian for R term, 17 placeholder (unused).

To add a new k.p term (say, a strain-induced coupling), you would:

1. Add a new entry to `kpterms` (QW) or `kpterms_2d` (wire) during `confinementInitialization`.
2. Compute the spatial profile of your coupling using the material parameters and FD stencils already in place.
3. Insert the new block into the 8x8 block topology in the assembly routines.

### Level 3: New post-processing (after eigenvalue solve)

If you need to compute derived quantities from the eigenvalues and eigenvectors, add a new module under `src/physics/`. The existing `gfactor_functions.f90` is the template: it reads eigenvectors, constructs perturbation Hamiltonians (by calling `ZB8bandQW` with `g='g'`), and computes matrix elements via Lowdin partitioning.

The pattern is:

1. Solve the main Hamiltonian to get eigenvalues and eigenvectors.
2. Build perturbation Hamiltonians using the same assembly routines but with different k-parameters or the `g` flag.
3. Compute observables from the perturbation matrix elements.

## Adding New Output

The output pipeline is straightforward. All output goes through `src/io/outputFunctions.f90`. The public routines are:

- `writeEigenvalues(smallk, eig, nsteps)` or `writeEigenvalues(smallk, eig, nsteps, cfg)` -- writes the band structure to `output/eigenvalues.dat`.
- `writeEigenfunctions(N, evnum, A, k, fdstep, z, is_bulk)` -- writes wavefunctions for a single k-point.
- `writeEigenfunctions2d(grid, evals, evecs, k, nev, write_parts)` -- writes 2D wavefunctions for wire mode.

To add a new output file:

1. Create a new subroutine in `outputFunctions.f90` (or a new module if it involves substantial physics).
2. Call `ensure_output_dir()` at the top to create `output/` if needed.
3. Use `call get_unit(iounit)` to get a free file unit number.
4. Write your data in a space-separated column format, with a `#` comment header line.
5. Call your routine from the appropriate place in `main.f90` or `main_gfactor.f90`, after the eigenvalue solve.

For 2D data intended for gnuplot `splot`, insert blank lines between y-rows (see the wire band edge profile output in `main.f90` for the pattern).

For consistency, use the format specifier `(g14.6)` or `(Ng14.6)` for numerical columns. This gives 14-character fields with 6 significant digits, which is sufficient for most visualization purposes.

## Adding New Tests

The test infrastructure has two layers: pFUnit unit tests and shell-based regression tests.

### Unit Tests (pFUnit)

Unit tests live in `tests/unit/` as `.pf` files. pFUnit uses a Python-like decorator syntax. Here is a minimal template for testing a new module:

```fortran
module test_myfeature
  use funit
  use definitions
  use my_new_module
  implicit none

contains

  @test
  subroutine test_basic_property()
    real(kind=dp) :: result
    result = my_function(1.0_dp, 2.0_dp)
    @assertEqual(3.0_dp, result, tolerance=1.0e-12_dp)
  end subroutine test_basic_property

  @test
  subroutine test_edge_case()
    ! Test boundary behavior
    @assertTrue(my_check(0.0_dp))
  end subroutine test_edge_case

end module test_myfeature
```

To wire it into the build system, add to `tests/CMakeLists.txt`:

```cmake
add_pfunit_ctest(test_myfeature
    TEST_SOURCES unit/test_myfeature.pf
    LINK_LIBRARIES 8bandkp_common
    LABELS "unit"
)
```

The `8bandkp_common` library is the static library that bundles all `src/` modules. Your test module can `use` any public entity from it.

### Regression Tests

Regression tests are shell scripts in `tests/integration/` that run the full executable and compare output against reference data. The template:

```bash
#!/bin/bash
# Integration test: <description>
set -euo pipefail

EXE="$1"           # Path to bandStructure or gfactorCalculation
CONFIG="$2"        # Path to .cfg file in tests/regression/configs/
REF_DIR="$3"       # Path to reference data in tests/regression/data/
COMPARE="$4"       # Path to compare_output.py

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

# Setup: copy config as input.cfg
/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"

# Run
cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: executable returned exit code $RC"
    cat test_output.log
    exit 1
fi

# Compare output
python3 "$COMPARE" "$REF_DIR/eigenvalues.dat" "$WORKDIR/output/eigenvalues.dat" --tolerance 1e-8
```

To add a new regression test:

1. Create a config file in `tests/regression/configs/` (e.g., `bulk_inp_kx.cfg`).
2. Run the executable with this config to generate reference output.
3. Copy the reference output to `tests/regression/data/<test_name>/`.
4. Create the shell script in `tests/integration/`.
5. Register it in `tests/CMakeLists.txt` with `add_test()`.

Run the full test suite with `ctest --test-dir build`. Use `-L unit` or `-L regression` to run only one category. Use `-V` for verbose output when debugging failures.

## Known Limitations

The following capabilities are deliberately excluded from the current codebase and would require substantial new physics modules:

**Non-resonant photoluminescence (NRPL).** nextnano Tutorial 5.9.13 computes NRPL spectra by combining the k.p eigenstates with carrier transport from a drift-diffusion solver. The drift-diffusion equation couples the electron and hole continuity equations with Poisson's equation and the current density ($J_n = q n \mu_n E + q D_n \nabla n$), requiring spatially resolved carrier mobilities, recombination rates, and generation profiles. This goes beyond the equilibrium Schrodinger-Poisson solver described in Chapter 7, which finds the ground-state charge distribution but does not simulate carrier dynamics or steady-state transport under optical excitation. Implementing NRPL support would require adding a drift-diffusion module to `src/physics/`, coupling it to the existing Poisson solver and optical absorption modules, and is deferred to future work.

## Feature Roadmap

The table below lists the major features needed to reproduce the full set of results from the target publications. Each feature is assigned a priority (1 = highest), the papers it unlocks (P4 = Paper 4, etc.), and an estimated effort.

| Priority | Feature | Unlocks | Effort |
|----------|---------|---------|--------|
| 1 | Wurtzite 8-band k.p Hamiltonian | P4, P6, P7 | Major |
| 2 | WZ material parameters | P4, P6, P7 | Medium |
| 3 | Circular wire geometry | P4, P8 | Medium |
| 4 | Interband absorption Im[epsilon] | P6 | Medium |
| 5 | DLP calculation | P6 | Small |
| 6 | Spin relaxation (DP, EY) | P8 | Major |
| 7 | SC self-energy (BdG) | P5 | Major |
| 8 | Unified WZ+ZB Hamiltonian | P6 | Medium |
| 9 | Ab initio fitting pipeline | P7 | Out of scope |

**Wurtzite 8-band k.p Hamiltonian** (Priority 1): The current code only implements the zincblende (ZB) 8-band Hamiltonian. Wurtzite requires a different basis and different coupling terms: the crystal field splitting replaces the cubic symmetry terms, and the Luttinger-like parameters have a different structure (gamma1_parallel, gamma1_perp, etc.). The implementation would follow the same pattern as `ZB8bandQW` but with a different 8x8 block topology. The `confinementInitialization` routine needs new k.p term indices for the additional anisotropy terms. A new subroutine `WZ8bandQW` in `hamiltonianConstructor.f90` is the natural home.

**WZ material parameters** (Priority 2): A wurtzite parameter database paralleling the existing `paramDatabase`. Sources: Vurgaftman (2001) for nitride parameters, Chuang's textbook for the wurtzite Hamiltonian structure. Add new `case` entries with names like "GaN_WZ", "AlN_WZ", "InN_WZ". The `paramStruct` type needs additional fields for the wurtzite-specific parameters (delta_cr, gamma1_par, gamma1_perp, gamma2_par, gamma2_perp, gamma3, etc.).

**Circular wire geometry** (Priority 3): The cut-cell immersed boundary framework is already in place (`geometry.f90`). The `wire_geometry` type supports `shape = 'circle'`. What remains is validating the cut-cell accuracy for cylindrical structures and adding the corresponding input parsing. The finite-difference stencil handles curved boundaries via the fractional cell volumes and face fractions.

**Interband absorption** (Priority 4): Compute the imaginary part of the dielectric function from the k.p eigenstates. Requires momentum matrix elements between conduction and valence band states, which can be extracted from the `dH/dk` perturbation already computed for g-factors. The optical_transition type in `defs.f90` is already scaffolded for this.

**DLP calculation** (Priority 5): Deformation potential theory for piezoelectric coupling. A post-processing step using strain and wavefunctions already available from the strain solver and eigenvalue solve. Small effort because it reuses existing data.

**Spin relaxation** (Priority 6): D'yakonov-Perel' and Elliott-Yafet mechanisms. Requires the k-dependent g-factor across the Brillouin zone (already partially available) plus spin-flip matrix elements. Major effort due to the complexity of the spin-orbit coupling terms and the need for time-dependent perturbation theory.

**SC self-energy (BdG)** (Priority 7): Bogoliubov-de Gennes extension for proximity-induced superconductivity in nanowires. Would add a particle-hole doubling of the Hamiltonian and a self-consistent gap equation. Major effort, as it requires a new diagonalization scheme (non-Hermitian) and coupling to a superconductor Green's function.

**Unified WZ+ZB Hamiltonian** (Priority 8): A single parameterized Hamiltonian that reduces to ZB or WZ depending on input. Medium effort once both individual Hamiltonians are validated.

**Ab initio fitting pipeline** (Priority 9): Tooling to fit k.p parameters from DFT band structures. Marked as out of scope for the core code -- this would be a separate utility.

## Build and Test Checklist

Before submitting any change, run through this checklist:

1. `cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl && cmake --build build` -- clean build.
2. `ctest --test-dir build` -- all tests pass (currently 38 tests: 15 unit + 23 regression).
3. Check for stale `.mod` files in the project root: `rm -f *.mod` if you see type mismatch errors.
4. Verify that `input.cfg` is not committed with personal test configs (use `tests/regression/configs/` instead).
5. If you changed `defs.f90` derived types or `hamiltonianConstructor.f90` Hamiltonian construction, flag for review per project policy.
