# 8bandkp-fdm

Fortran 2018 solver for the **8-band zinc-blende k.p Hamiltonian** via finite differences.
Computes electronic band structures for bulk semiconductors, quantum wells, and quantum wires. Calculates Landau g-factors via second-order Lowdin partitioning with commutator-based velocity operators. Computes optical absorption, gain, spontaneous emission, and intersubband transitions using commutator-based velocity matrices. Includes a self-consistent Schrodinger-Poisson solver with DIIS acceleration.

Based on Chuang & Chang (1997) with Foreman renormalization and Winkler g-factor formalism. GPL v3.0, authored by Tiago de Campos.

If you use this code, please cite: [10.1021/acsaelm.0c00269](https://doi.org/10.1021/acsaelm.0c00269)

## Build

**Prerequisites:** Fortran compiler (gfortran, supports F2018), BLAS/LAPACK, FFTW3, CMake >= 3.15. For best performance use Intel MKL. Optional: [fortran-stdlib](https://github.com/fortran-lang/stdlib) (detected automatically by CMake).

Ubuntu:
```bash
sudo apt install gfortran gcc g++ liblapack-dev libblas-dev libfftw3-dev
```

```bash
# Configure (with MKL)
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl

# Build
cmake --build build

# Or use the Make wrapper
make all       # Configure + build all executables
make run       # Build and run bandStructure
make clean     # Remove build/ directory
```

Executables: `build/src/bandStructure`, `build/src/gfactorCalculation`, `build/src/opticalProperties`, and `build/src/topologicalAnalysis`.

**Experimental fpm manifest** -- `fpm.toml` is provided as a project manifest but does not encode MKL/FFTW3 link flags. Use CMake for production builds.

## Quick Start

Write one of the configs below to `input.toml`, then run `./build/src/bandStructure`, `./build/src/gfactorCalculation`, `./build/src/opticalProperties`, or `./build/src/topologicalAnalysis`. See `tests/regression/configs/` for more examples and [docs/reference/input-reference.md](docs/reference/input-reference.md) for the full format.

**Bulk** (GaAs, k along x):
```toml
confinement = "bulk"
FDorder = 2
fd_step = 101

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 11

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "GaAs"
```

**Quantum well** (GaAs/AlGaAs):
```toml
confinement = "qw"
FDorder = 2
fd_step = 201

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 21

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "Al30Ga70As"
z_min = -200
z_max = 200

[[material]]
name = "GaAs"
z_min = -50
z_max = 50
```

**g-factor** (bulk GaAs conduction band):
```toml
confinement = "bulk"
FDorder = 2
fd_step = 1
which_band = 0
band_idx = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "GaAs"
```

**Self-consistent SP** (GaAs/AlAs QW, n-doped well):
```toml
confinement = "qw"
FDorder = 2
fd_step = 101
which_band = 0
band_idx = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "AlAs"
z_min = -150
z_max = 150

[[material]]
name = "GaAs"
z_min = -50
z_max = 50

[sc]
max_iterations = 50
tolerance = 1.0e-6
mixing_alpha = 0.3
diis_history = 7
temperature = 300.0
fermi_mode = "charge_neutrality"
num_kpar = 41
kpar_max = 0.2
bc_type = "DD"

[[doping]]
ND = 0.0
NA = 0.0

[[doping]]
ND = 1.0e18
NA = 0.0
```

**Quantum wire** (GaAs rectangular cross section):
```toml
confinement = "wire"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "kz"
max = 0.1
nsteps = 5

[bands]
num_cb = 8
num_vb = 16

[wire]
nx = 21
ny = 21
dx = 3.0
dy = 3.0

[wire.geometry]
shape = "rectangle"
width = 63.0
height = 63.0

[[region]]
material = "GaAs"
inner = 0.0
outer = 100.0

[feast]
emin = -1.5
emax = 2.0
m0 = 128
```

## Documentation

**Lectures** (`docs/lecture/`) -- pedagogical progression from k.p theory to quantum wires:

[00 Quickstart](docs/lecture/00-quickstart.md) | [01 Bulk](docs/lecture/01-bulk-band-structure.md) | [02 QW](docs/lecture/02-quantum-well.md) | [03 Wavefunctions](docs/lecture/03-wavefunctions.md) | [04 Strain](docs/lecture/04-strain.md) | [05 g-Factor](docs/lecture/05-gfactor.md) | [06 Optical](docs/lecture/06-optical-properties.md) | [07 SC-SP](docs/lecture/07-self-consistent-sp.md) | [08 Wire](docs/lecture/08-quantum-wire.md) | [09 Numerics](docs/lecture/09-numerical-methods.md) | [10 QCSE](docs/lecture/10-qcse.md) | [11 Convergence](docs/lecture/11-convergence.md) | [12 Extending](docs/lecture/12-extending-the-code.md) | [13 Topological](docs/lecture/13-topological-superconductivity.md) | [14 Excitons](docs/lecture/14-excitons-scattering.md)

**Reference** (`docs/reference/`): [Input parameters](docs/reference/input-reference.md) | [Output files](docs/reference/output-reference.md) | [Benchmarks](docs/reference/benchmarks.md)

Config files are TOML format (`input.toml`). See the [input reference](docs/reference/input-reference.md) for the complete schema.

## Testing

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-<ver>/cmake
cmake --build build
make test
# or: ctest --test-dir build              # all
#     ctest --test-dir build -L unit      # pFUnit unit tests
#     ctest --test-dir build -L regression  # regression vs golden output
```

## Architecture

```
defs -> parameters -> utils -> finitedifferences -> hamiltonian_blocks -> hamiltonianConstructor -> gfactor_functions
                                  finitedifferences -> hamiltonian_blocks -> hamiltonian_wire
                                  finitedifferences -> poisson
                                  finitedifferences -> charge_density -> sc_loop
                   input_parser    outputFunctions
```

Executables: `bandStructure` (band structure vs k), `gfactorCalculation` (g-factor at k=0), `opticalProperties` (optical spectra for bulk/QW/wire), `topologicalAnalysis` (topological invariants).

## License

GNU General Public License v3.0
