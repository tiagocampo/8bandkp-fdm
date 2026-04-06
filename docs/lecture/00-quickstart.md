# Chapter 00: Quickstart Guide

This chapter gets you from zero to your first band-structure calculation in
under five minutes. We build the code, run a bulk GaAs calculation, and walk
through every line of output so you know exactly what the numbers mean.

If you hit a snag, jump straight to [Section 6: Common Issues](#6-common-issues).

---

## 1. Prerequisites

| Package | Minimum version | How to check |
|---|---|---|
| GNU Fortran (`gfortran`) | 9.x | `gfortran --version` |
| Intel MKL | 2019 or later | `ls $MKLROOT/lib/cmake/mkl` |
| CMake | 3.15 | `cmake --version` |
| Ninja (optional) | any | `ninja --version` |
| FFTW3 | 3.x | `pkg-config --cflags fftw3` |

### Environment setup

MKL must be discoverable via the `MKLROOT` environment variable. If you use
Intel's `setvars.sh` or oneAPI installer, `MKLROOT` is set automatically.
Otherwise, export it yourself:

```bash
export MKLROOT=/opt/intel/oneapi/mkl/latest    # adjust path to your install
```

FFTW3 headers are needed at compile time only. On Debian/Ubuntu:

```bash
sudo apt install libfftw3-dev
```

On Arch/Manjaro:

```bash
sudo pacman -S fftw
```

---

## 2. Clone and Build

Three commands:

```bash
git clone <repository-url> 8bandkp-fdm
cd 8bandkp-fdm
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl
cmake --build build
```

If you prefer Makefiles over Ninja, drop the `-G Ninja` flag.

When the build succeeds you will find two executables:

```
build/src/bandStructure        # band-structure sweeps
build/src/gfactorCalculation   # Landau g-factor at Gamma
```

You can safely ignore the `gfactorCalculation` executable for now; we will
only use `bandStructure` in this guide.

---

## 3. First Run: Bulk GaAs

The program reads its input from a file called `input.cfg` in the project root.
Create it with the following contents:

```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 11
confinement:  0
FDstep: 101
FDorder: 2
numLayers:  1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0005
```

**What this input means:**

| Parameter | Value | Meaning |
|---|---|---|
| `waveVector` | `kx` | Sweep along the [100] direction |
| `waveVectorMax` | `0.1` | Maximum k in 1/Angstrom |
| `waveVectorStep` | `11` | Number of k-points (including k=0) |
| `confinement` | `0` | Bulk mode: 8x8 Hamiltonian |
| `FDstep` | `101` | Ignored in bulk; code forces it to 1 |
| `FDorder` | `2` | Finite-difference order (irrelevant for bulk) |
| `numLayers` | `1` | Single material layer |
| `material1` | `GaAs` | Gallium arsenide from the built-in database |
| `numcb` | `2` | Request 2 conduction-band eigenvalues |
| `numvb` | `6` | Request 6 valence-band eigenvalues |
| `ExternalField` | `0 EF` | No external electric field |
| `EFParams` | `0.0005` | Field strength (unused when field is off) |

Now run the program:

```bash
./build/src/bandStructure
```

### 3.1 Standard output

The program prints a summary of parsed parameters followed by the material
database entry for GaAs:

```
 waveVector:kx
 waveVectorMax:  0.10000000000000001
 waveVectorStep:          11
 confinement:           0
 FDstep:         101
 FDorder:           2
 numLayers:           1
 Warning: bulk mode requires fdStep=1. Forcing fdStep=1.
 material1:GaAs
 numcb:           2
 numvb:           6
 ExternalField:           0 EF
 EFParams:   5.0000000000000001E-004

 Material: GaAs
 Parameters
 EP :   28.800000000000001
 P  :   10.475088634541221
 A  :   14.925373134328357
 gamma1 :   6.9800000000000004
 gamma2 :   2.0600000000000001
 gamma3 :   2.9300000000000002
```

The "Warning" about `fdStep` is harmless: bulk mode solves an 8x8 matrix at
each k-point and does not use finite differences, so the code overrides the
value to 1.

The material parameters block shows the Kane energy $E_P = 28.8$ eV, the
interband momentum matrix element $P = \sqrt{E_P \cdot \hbar^2/(2m_0)}
\approx 10.48$ eV-Angstrom, and the Luttinger parameters $\gamma_1$, $\gamma_2$,
$\gamma_3$. These are the Vurgaftman 2001 values for GaAs.

### 3.2 Output file: `output/eigenvalues.dat`

The program writes eigenvalues to `output/eigenvalues.dat`. The first few lines
look like this:

```
 #k, values
    0.00000      -0.341000      -0.341000       -0.00000       -0.00000       -0.00000       -0.00000        1.51900        1.51900
   0.100000E-01  -0.343693      -0.343693      -0.585596E-02  -0.585596E-02  -0.286000E-03  -0.286000E-03    1.52723        1.52723
   0.200000E-01  -0.352123      -0.352123      -0.226020E-01  -0.226020E-01  -0.114400E-02  -0.114400E-02    1.55146        1.55146
```

**Column format:**

| Column | Content |
|---|---|
| 1 | Wave vector $k$ in 1/Angstrom |
| 2--9 | Eight eigenvalues in eV, sorted from lowest to highest |

The file has 12 rows: a header line (`#k, values`) plus 11 k-points (the
`waveVectorStep` value).

---

## 4. Interpreting the Output

### 4.1 Eigenvalues at the Gamma point (k = 0)

Focusing on the first data row (k = 0):

| Band index | Eigenvalue (eV) | Identity |
|---|---|---|
| 1 | -0.341 | Split-off hole (SO), spin down |
| 2 | -0.341 | Split-off hole (SO), spin up |
| 3 | 0.000 | Heavy hole (HH), spin up |
| 4 | 0.000 | Light hole (LH), spin up |
| 5 | 0.000 | Light hole (LH), spin down |
| 6 | 0.000 | Heavy hole (HH), spin down |
| 7 | +1.519 | Conduction band (CB), spin down |
| 8 | +1.519 | Conduction band (CB), spin up |

Three features to notice:

1. **Band gap.** The valence-band top (HH/LH at 0.000 eV) and the conduction-band
   bottom (CB at +1.519 eV) are separated by exactly the GaAs band gap
   $E_g = 1.519$ eV. This confirms the code is using Vurgaftman's published
   value.

2. **Spin-orbit splitting.** The split-off bands sit at $-0.341$ eV, which is
   $\Delta_{\text{SO}} = 0.341$ eV below the valence-band edge. Again, this is
   the accepted GaAs value.

3. **HH/LH degeneracy at Gamma.** The four valence bands (HH and LH, both spins)
   are all exactly degenerate at k = 0. This is a consequence of the
   $\Gamma_8$ irreducible representation of the zincblende point group. The
   degeneracy lifts as soon as k departs from zero.

### 4.2 The full dispersion

The following figure shows all eight bands across the full k-range:

![Bulk GaAs 8-band E(k) dispersion](../figures/bulk_gaas_bands.png)

*Figure 1: Bulk GaAs 8-band E(k) dispersion along [100], computed with the
input above. The conduction band (upper curve) has a small effective mass
($m^* \approx 0.067\,m_0$), so it curves gently. The heavy-hole band is
nearly flat (large effective mass). The light-hole band has a much smaller
effective mass and bends away more quickly.*

---

## 5. Next Steps

Now that you have a working build and understand the output format, you can
explore further:

| If you want to... | Read... |
|---|---|
| Understand the 8-band Hamiltonian and why the bands look like this | [Chapter 01: Bulk Band Structure](01-bulk-band-structure.md) |
| Simulate a quantum well with confinement and subbands | [Chapter 02: Quantum Well](02-quantum-well.md) |
| Visualize wavefunctions and probability densities | [Chapter 03: Wavefunctions](03-wavefunctions.md) |
| Add a self-consistent Schrodinger-Poisson loop | [Chapter 07: Self-Consistent SP](07-self-consistent-sp.md) |
| Compute Landau g-factors | [Chapter 05: g-Factor](05-gfactor.md) |
| Understand the finite-difference machinery | [Chapter 09: Numerical Methods](09-numerical-methods.md) |

The remaining input files in `tests/regression/configs/` are ready-to-run
examples for quantum wells, g-factors, self-consistent calculations, and more.
Copy any of them to `input.cfg` and re-run `./build/src/bandStructure`.

---

## 6. Common Issues

### Stale `.mod` files

If you see type-mismatch errors like `Symbol 'some_type' at (1) has no IMPLICIT
type`, old `.mod` files in the project root may be shadowing the fresh ones in
`build/`. Fix:

```bash
rm -f *.mod
cmake --build build
```

### `cp -i` shell alias

If your shell aliases `cp` to `cp -i` (interactive), copying `input.cfg` via a
script will hang waiting for confirmation. Use your editor or write the file
directly instead of `cp`.

### MKL not found

CMake error: `Could not find MKL`. Make sure `MKLROOT` points to a valid MKL
installation and pass the cmake hint:

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl
```

If you do not have Intel MKL installed, see Intel's oneAPI documentation for
your distribution.

### OMP library not found

CMake error: `OMP_LIBRARY has an invalid value "OMP_LIBRARY-NOTFOUND"`. This
occurs when MKL is configured for `intel_thread` threading but Intel's OpenMP
runtime (`libiomp5`) is not installed. If you use gfortran without the Intel
compiler suite, switch MKL to sequential threading:

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl -DMKL_THREADING=sequential
```

### FFTW3 missing

CMake error: `Could not find FFTW3`. Install the development headers:

```bash
# Debian / Ubuntu
sudo apt install libfftw3-dev

# Arch / Manjaro
sudo pacman -S fftw

# Fedora
sudo dnf install fftw-devel
```

### Build directory is stale

If you switch compilers or change the MKL path, a clean reconfigure is safest:

```bash
rm -rf build
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl
cmake --build build
```

### Wrong number of eigenvalues

The code returns `numcb + numvb` eigenvalues. For bulk (8-band), the maximum
is 8. Requesting `numcb=2, numvb=6` gives all eight bands. Requesting fewer
(e.g., `numcb=1, numvb=3`) returns only the four lowest eigenvalues. If you
need all bands, keep the sum at 8.

### No `output/` directory

The program creates `output/` automatically on first run. If for some reason it
does not (permissions issue), create it manually:

```bash
mkdir -p output
```
