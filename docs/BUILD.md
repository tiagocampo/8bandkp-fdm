# Build Documentation: 8bandkp-fdm-ai

**Date**: 2025-01-27  
**Status**: Complete

## Prerequisites

### Required Software
- **Fortran Compiler**: gfortran 15.2.1+ (or ifort/ifx)
- **Intel MKL**: 2025.0+ (for sparse matrix operations)
- **FFTW3**: 3.3+ (for Fourier transforms)
- **Make**: Build system
- **gnuplot**: 5.0+ (for plotting, optional)

### Installation Examples

**Ubuntu/Debian**:
```bash
sudo apt install gfortran gcc g++ gnuplot make
# Install Intel OneAPI Base Toolkit for MKL
```

**CentOS/RHEL**:
```bash
sudo yum install gcc-gfortran gnuplot make
# Install Intel OneAPI Base Toolkit for MKL
```

**Intel OneAPI Installation**:
```bash
# Download and install Intel OneAPI Base Toolkit
# Source environment: source /opt/intel/oneapi/setvars.sh
```

## Build Process

### Step 1: Verify Dependencies
```bash
# Check Fortran compiler
gfortran --version

# Check MKL availability
ls /opt/intel/oneapi/mkl/2025.0/lib/intel64/

# Check FFTW3
pkg-config --exists fftw3 && echo "FFTW3 found"
```

### Step 2: Configure Build
The Makefile is pre-configured for Intel MKL. If you need to use standard BLAS/LAPACK instead:

1. Comment out the MKL LDFLAGS section
2. Uncomment the standard BLAS/LAPACK section
3. Update library paths as needed

### Step 3: Compile
```bash
# Clean previous builds
make clean_all

# Build both executables
make all

# Verify executables created
ls -la bandStructure gfactorCalculation
```

### Step 4: Test Build
```bash
# Test band structure calculation (current behavior uses input.cfg by default)
cp examples/bulk_InAs60Sb40.example input.cfg
./bandStructure

# OR (after CLI filename support is added)
# ./bandStructure examples/bulk_InAs60Sb40.example --out outputs/$(date +%Y%m%d-%H%M%S)

# Test g-factor calculation
cp examples/gfactor_quantum_well.example input.cfg
./gfactorCalculation
# OR: ./gfactorCalculation examples/gfactor_quantum_well.example --out outputs/$(date +%Y%m%d-%H%M%S)
```

## Makefile Configuration

### Current MKL Configuration
```makefile
LDFLAGS=-I/usr/include -L/usr/lib64 -L/opt/intel/oneapi/mkl/2025.0/lib/intel64 \
        -I/opt/intel/oneapi/mkl/2025.0/include \
        -I/opt/intel/oneapi/mkl/2025.0/include/fftw \
        -L/opt/intel/oneapi/compiler/2025.0/lib \
        -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread \
        -lmkl_core -liomp5 -lpthread -lm -ldl -lfftw3
```

### Alternative BLAS/LAPACK Configuration
```makefile
LDFLAGS=-I/usr/include -L/usr/lib64 -lm -lfftw3 -llapack -lblas
```

## Build Targets

- `make all`: Build both executables
- `make bandStructure`: Build band structure calculator only
- `make gfactor`: Build g-factor calculator only
- `make clean`: Remove object files and data files
- `make clean_all`: Remove all build artifacts including executables
- `make run`: Build and run bandStructure

## Troubleshooting

### Common Issues

**1. MKL Library Not Found**
```
/usr/bin/ld: cannot find -lmkl_gf_lp64
```
**Solution**: Verify Intel MKL installation and update paths in Makefile

**2. OpenMP Library Not Found**
```
/usr/bin/ld: cannot find -liomp5
```
**Solution**: Add OpenMP library path: `-L/opt/intel/oneapi/compiler/2025.0/lib`

**3. FFTW3 Not Found**
```
/usr/bin/ld: cannot find -lfftw3
```
**Solution**: Install FFTW3 development package or update library path

**4. Undefined MKL Sparse Functions**
```
undefined reference to `MKL_SPARSE_Z_MV'
```
**Solution**: Ensure MKL sparse matrix libraries are linked correctly

### Build Verification

After successful build, verify:
1. Both executables are created and executable
2. No linking errors in build output
3. Test runs complete without segmentation faults
4. Output files are generated under `outputs/<run-id>/` (or repository root for legacy runs) and include `run.meta`

## Performance Notes

- **Build Time**: ~30 seconds on modern systems
- **Executable Size**: bandStructure (~343KB), gfactorCalculation (~497KB)
- **Memory Usage**: Depends on problem size and discretization
- **Optimization**: Uses `-O2 -flto -funroll-loops` for performance

## Dependencies

### Required Libraries
- **Intel MKL**: Sparse matrix operations, BLAS, LAPACK
- **FFTW3**: Fast Fourier transforms
- **OpenMP**: Parallel processing support
- **Standard C Library**: Basic system functions

### Optional Libraries
- **gnuplot**: For result visualization
- **Intel TBB**: Alternative threading library

## Next Steps

1. Run example calculations to verify functionality
2. Create custom input files for specific materials
3. Use plotting scripts for result visualization
4. Refer to VALIDATION.md for result verification
