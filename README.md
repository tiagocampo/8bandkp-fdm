# 8bandkp-fdm

Fortran 90 solver for the **8-band zinc-blende k.p Hamiltonian** via finite differences.

Computes electronic band structures for bulk semiconductors, quantum wells, and quantum wires. Calculates Landau g-factors via second-order Lowdin partitioning. Includes a self-consistent Schrodinger-Poisson solver with DIIS acceleration.

Based on Chuang & Chang (1997) with Foreman renormalization and Winkler g-factor formalism. GPL v3.0, authored by Tiago de Campos.

If you use this code, please cite: [10.1021/acsaelm.0c00269](https://doi.org/10.1021/acsaelm.0c00269)

## Build

**Prerequisites:** Fortran compiler (gfortran), BLAS/LAPACK, FFTW3, CMake >= 3.15. For best performance use Intel MKL.

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
make all       # Configure + build both executables
make run       # Build and run bandStructure
make clean     # Remove build/ directory
```

Executables: `build/src/bandStructure` and `build/src/gfactorCalculation`.

## Quick Start

Write one of the configs below to `input.cfg`, then run `./build/src/bandStructure` or `./build/src/gfactorCalculation`. See `tests/regression/configs/` for more examples and [docs/reference/input-reference.md](docs/reference/input-reference.md) for the full format.

**Bulk** (GaAs, k along x):
```
waveVector: kx            waveVectorMax: 0.1      waveVectorStep: 11
confinement: 0            FDstep: 101             FDorder: 2
numLayers: 1              material1: GaAs
numcb: 2                  numvb: 6
ExternalField: 0 EF       EFParams: 0.0005
```

**Quantum well** (AlSbW/GaSbW/InAsW):
```
waveVector: kx            waveVectorMax: 0.1      waveVectorStep: 11
confinement: 1            FDstep: 101             FDorder: 2
numLayers: 3              numcb: 32               numvb: 32
material1: AlSbW -250  250 0
material2: GaSbW -135  135 0.2414
material3: InAsW  -35   35 -0.0914
ExternalField: 0 EF       EFParams: 0.0005
```

**g-factor** (bulk GaAs conduction band):
```
waveVector: k0            waveVectorMax: 0.1      waveVectorStep: 0
confinement: 0            FDstep: 1               FDorder: 2
numLayers: 1              material1: GaAs
numcb: 2                  numvb: 6                whichBand: 0   bandIdx: 1
ExternalField: 0 EF       EFParams: 0.0005
```

**Self-consistent SP** (GaAs/AlAs QW, n-doped well):
```
waveVector: k0            waveVectorMax: 0.0      waveVectorStep: 1
confinement: 1            FDstep: 101             FDorder: 2
numLayers: 3              numcb: 4                numvb: 8
material1: AlAs -150 150 0   material2: GaAs -50 50 0   material3: AlAs -150 150 0
ExternalField: 0 EF       EFParams: 0.0005
SC: 1   max_iter: 50   tolerance: 1.0e-6   mixing_alpha: 0.3   diis_history: 7
temperature: 300.0   fermi_mode: 1   fermi_level: 1.5   num_kpar: 21   kpar_max: 0.1
bc_type: DD   doping1: 0.0 0.0   doping2: 1.0e18 0.0   doping3: 0.0 0.0
```

**Quantum wire** (GaAs rectangular cross section):
```
waveVector: kz            waveVectorMax: 0.1      waveVectorStep: 5
confinement: 2            wire_nx: 11             wire_ny: 11
wire_dx: 2.0              wire_dy: 2.0            wire_shape: rectangle
wire_width: 22.0          wire_height: 22.0
numRegions: 1             region: GaAs 0.0 100.0  numcb: 4   numvb: 4
```

## Documentation

**Lectures** (`docs/lecture/`) -- pedagogical progression from k.p theory to quantum wires:

[00 Quickstart](docs/lecture/00-quickstart.md) | [01 Bulk](docs/lecture/01-bulk-band-structure.md) | [02 QW](docs/lecture/02-quantum-well.md) | [03 Wavefunctions](docs/lecture/03-wavefunctions.md) | [04 Strain](docs/lecture/04-strain.md) | [05 g-Factor](docs/lecture/05-gfactor.md) | [06 Optical](docs/lecture/06-optical-properties.md) | [07 SC-SP](docs/lecture/07-self-consistent-sp.md) | [08 Wire](docs/lecture/08-quantum-wire.md) | [09 Numerics](docs/lecture/09-numerical-methods.md) | [10 QCSE](docs/lecture/10-qcse.md) | [11 Convergence](docs/lecture/11-convergence.md) | [12 Extending](docs/lecture/12-extending-the-code.md)

**Reference** (`docs/reference/`): [Input parameters](docs/reference/input-reference.md) | [Output files](docs/reference/output-reference.md) | [Benchmarks](docs/reference/benchmarks.md)

`input.cfg` guidance:
- Follow the canonical block order used in `tests/regression/configs/`.
- Optional block entry labels are name-based (`optics:`, `exciton:`, `scattering:`, `feast_emin:`, `strain:`), but parameters inside each block still use the documented order.

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
defs -> parameters -> utils -> finitedifferences -> hamiltonianConstructor -> gfactor_functions
                                  finitedifferences -> poisson
                                  finitedifferences -> charge_density -> sc_loop
                   input_parser    outputFunctions
```

Executables: `bandStructure` (band structure vs k), `gfactorCalculation` (g-factor at k=0).

## License

GNU General Public License v3.0
