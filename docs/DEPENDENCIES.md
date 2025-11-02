# Dependency Requirements: 8bandkp-fdm-ai

**Date**: 2025-01-27  
**Status**: Complete

## System Requirements

### Operating System
- **Linux/Unix**: Primary target platform
- **Architecture**: x86_64 (64-bit)
- **Kernel**: Linux 6.12+ (tested on Manjaro)

### Compiler Requirements
- **Fortran Compiler**: gfortran 15.2.1+ (recommended)
- **Alternative**: Intel Fortran (ifort/ifx) 2021.0+
- **C Compiler**: gcc 12.0+ (for linking)
- **Make**: GNU Make 4.0+

## Core Dependencies

### Intel MKL (Required)
**Purpose**: Sparse matrix operations, BLAS, LAPACK  
**Version**: 2025.0+  
**Installation**: Intel OneAPI Base Toolkit  
**Location**: `/opt/intel/oneapi/mkl/2025.0/`

**Required Libraries**:
- `libmkl_gf_lp64.so` - 64-bit integer interface
- `libmkl_intel_thread.so` - Intel threading
- `libmkl_core.so` - Core MKL functions
- `libiomp5.so` - Intel OpenMP runtime

**Sparse Matrix Functions Used**:
- `MKL_SPARSE_Z_CREATE_COO` - Create coordinate format matrix
- `MKL_SPARSE_CONVERT_CSR` - Convert to CSR format
- `MKL_SPARSE_Z_MV` - Sparse matrix-vector multiplication

### FFTW3 (Required)
**Purpose**: Fast Fourier transforms  
**Version**: 3.3+  
**Installation**: System package manager or source  
**Library**: `libfftw3.so`

### OpenMP (Required)
**Purpose**: Parallel processing support  
**Version**: 5.0+  
**Implementation**: Intel OpenMP (liomp5) or GNU OpenMP (gomp)

## Optional Dependencies

### gnuplot (Optional)
**Purpose**: Result visualization and plotting  
**Version**: 5.0+  
**Installation**: System package manager  
**Usage**: For creating publication-quality plots

### Intel TBB (Optional)
**Purpose**: Alternative threading library  
**Version**: 2021.0+  
**Installation**: Intel OneAPI Base Toolkit  
**Usage**: Alternative to OpenMP for threading

## Installation Instructions

### Ubuntu/Debian
```bash
# Install basic tools
sudo apt update
sudo apt install gfortran gcc g++ make gnuplot

# Install FFTW3
sudo apt install libfftw3-dev

# Install Intel OneAPI Base Toolkit
# Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html
wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/...
sudo sh ./l_BaseKit_p_2025.0.1_offline.sh
```

### CentOS/RHEL
```bash
# Install basic tools
sudo yum install gcc-gfortran gcc-c++ make gnuplot

# Install FFTW3
sudo yum install fftw3-devel

# Install Intel OneAPI Base Toolkit
# Download and install from Intel website
```

### Arch/Manjaro
```bash
# Install basic tools
sudo pacman -S gcc-fortran gcc make gnuplot

# Install FFTW3
sudo pacman -S fftw

# Install Intel OneAPI Base Toolkit
# Use AUR package or download from Intel
yay -S intel-oneapi-basekit
```

## Environment Setup

### Intel OneAPI Environment
```bash
# Source Intel OneAPI environment
source /opt/intel/oneapi/setvars.sh

# Verify MKL is available
echo $MKLROOT
ls $MKLROOT/lib/intel64/
```

### Library Path Configuration
```bash
# Add MKL to library path
export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2025.0/lib/intel64:$LD_LIBRARY_PATH

# Add OpenMP to library path
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2025.0/lib:$LD_LIBRARY_PATH
```

## Verification Commands

### Check Compiler
```bash
gfortran --version
# Expected: GNU Fortran (GCC) 15.2.1+
```

### Check MKL
```bash
ls /opt/intel/oneapi/mkl/2025.0/lib/intel64/
# Expected: libmkl_*.so files present
```

### Check FFTW3
```bash
pkg-config --exists fftw3 && echo "FFTW3 found"
# Expected: FFTW3 found
```

### Check OpenMP
```bash
ls /opt/intel/oneapi/compiler/2025.0/lib/libiomp5.so
# Expected: File exists
```

## Troubleshooting

### Common Issues

**1. MKL Not Found**
- Verify Intel OneAPI installation
- Check environment variables
- Update library paths in Makefile

**2. OpenMP Issues**
- Ensure Intel OpenMP is installed
- Check library path configuration
- Verify threading support

**3. FFTW3 Missing**
- Install development package
- Check pkg-config configuration
- Verify library linking

**4. Compiler Issues**
- Update to supported Fortran compiler
- Check compiler flags compatibility
- Verify C++ standard library

### Build Verification
```bash
# Test build process
make clean_all
make all

# Verify executables
ls -la bandStructure gfactorCalculation

# Test basic functionality with explicit input and outputs directory
./bandStructure examples/bulk_InAs60Sb40.example --out outputs/$(date +%Y%m%d-%H%M%S)/bulk
```

## Performance Considerations

### Memory Requirements
- **Minimum**: 4GB RAM
- **Recommended**: 8GB+ RAM
- **Large Problems**: 16GB+ RAM

### CPU Requirements
- **Minimum**: 2 cores
- **Recommended**: 4+ cores
- **Parallel Processing**: OpenMP support required

### Storage Requirements
- **Source Code**: ~50MB
- **Build Artifacts**: ~100MB
- **Example Data**: ~10MB
- **Output Files**: Variable (depends on problem size)

## Version Compatibility

### Tested Combinations
- **gfortran 15.2.1** + **Intel MKL 2025.0** + **FFTW3 3.3.10**
- **gfortran 12.0** + **Intel MKL 2024.0** + **FFTW3 3.3.8**
- **ifort 2021.0** + **Intel MKL 2021.0** + **FFTW3 3.3.7**

### Known Issues
- **gfortran < 10.0**: May have compatibility issues with MKL
- **MKL < 2020.0**: Missing some sparse matrix functions
- **FFTW3 < 3.3**: May have performance issues

## Next Steps

1. Verify all dependencies are installed correctly
2. Run build verification tests
3. Test with example input files
4. Refer to BUILD.md for compilation instructions
5. Use VALIDATION.md for result verification
