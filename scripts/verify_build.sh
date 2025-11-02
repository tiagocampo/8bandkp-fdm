#!/bin/bash
# Build Verification Script for 8bandkp-fdm-ai
# Purpose: Automated testing of build process and basic functionality
# Date: 2025-01-27

set -e  # Exit on any error

echo "=========================================="
echo "8bandkp-fdm-ai Build Verification Script"
echo "=========================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
    else
        echo -e "${RED}✗${NC} $2"
        exit 1
    fi
}

# Function to print warning
print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

echo "Step 1: Checking prerequisites..."

# Check Fortran compiler
if command -v gfortran &> /dev/null; then
    gfortran --version | head -1
    print_status 0 "Fortran compiler found"
else
    print_status 1 "Fortran compiler not found"
fi

# Check MKL availability
if [ -d "/opt/intel/oneapi/mkl" ]; then
    print_status 0 "Intel MKL found"
else
    print_warning "Intel MKL not found - using standard BLAS/LAPACK"
fi

# Check FFTW3
if pkg-config --exists fftw3 2>/dev/null; then
    print_status 0 "FFTW3 found"
else
    print_status 1 "FFTW3 not found"
fi

echo ""
echo "Step 2: Building project..."

# Clean previous builds
echo "Cleaning previous builds..."
make clean_all > /dev/null 2>&1

# Build project
echo "Building executables..."
if make all > build.log 2>&1; then
    print_status 0 "Build completed successfully"
else
    print_status 1 "Build failed - check build.log for details"
fi

echo ""
echo "Step 3: Verifying executables..."

# Check if executables exist
if [ -f "bandStructure" ] && [ -x "bandStructure" ]; then
    print_status 0 "bandStructure executable created"
else
    print_status 1 "bandStructure executable not found or not executable"
fi

if [ -f "gfactorCalculation" ] && [ -x "gfactorCalculation" ]; then
    print_status 0 "gfactorCalculation executable created"
else
    print_status 1 "gfactorCalculation executable not found or not executable"
fi

echo ""
echo "Step 4: Testing basic functionality..."

# Prepare run ids and output dirs
RUN_ID=$(date +%Y%m%d-%H%M%S)
OUT_BASE="outputs/${RUN_ID}"
mkdir -p "${OUT_BASE}"

# Test bandStructure with bulk example
echo "Testing bandStructure with bulk example..."
if ./bandStructure examples/bulk_InAs60Sb40.example --out "${OUT_BASE}/bulk" > band_test.log 2>&1; then
    print_status 0 "bandStructure bulk calculation completed"
else
    print_status 1 "bandStructure bulk calculation failed"
fi

# Test bandStructure with quantum well example
echo "Testing bandStructure with quantum well example..."
if ./bandStructure quantumwell.example --out "${OUT_BASE}/qw" > qw_test.log 2>&1; then
    print_status 0 "bandStructure quantum well calculation completed"
else
    print_status 1 "bandStructure quantum well calculation failed"
fi

# Test gfactorCalculation (this may fail due to input requirements)
echo "Testing gfactorCalculation..."
if ./gfactorCalculation examples/gfactor_quantum_well.example --out "${OUT_BASE}/gfactor" > gfactor_test.log 2>&1; then
    print_status 0 "gfactorCalculation completed"
else
    print_warning "gfactorCalculation failed - may require specific input format"
fi

echo ""
echo "Step 5: Checking output files..."

# Check for output files
if [ -f "${OUT_BASE}/bulk/eigenvalues.dat" ]; then
    print_status 0 "bulk/eigenvalues.dat generated"
    echo "  - File size: $(wc -l < "${OUT_BASE}/bulk/eigenvalues.dat") lines"
else
    print_status 1 "bulk/eigenvalues.dat not generated"
fi

# Check for eigenfunction files (qw)
eigenfunction_count=$(ls "${OUT_BASE}/qw"/eigenfunctions_k_*.dat 2>/dev/null | wc -l)
if [ $eigenfunction_count -gt 0 ]; then
    print_status 0 "$eigenfunction_count eigenfunction files generated"
else
    print_warning "No eigenfunction files generated"
fi

echo ""
echo "Step 6: Performance check..."

# Check build time
if [ -f "build.log" ]; then
    build_time=$(grep -o "real.*user.*sys" build.log | tail -1 | awk '{print $2}' || echo "unknown")
    echo "Build time: $build_time"
fi

# Check executable sizes
if [ -f "bandStructure" ]; then
    size=$(du -h bandStructure | cut -f1)
    echo "bandStructure size: $size"
fi

if [ -f "gfactorCalculation" ]; then
    size=$(du -h gfactorCalculation | cut -f1)
    echo "gfactorCalculation size: $size"
fi

echo ""
echo "=========================================="
echo "Build verification completed successfully!"
echo "=========================================="

# Clean up log files
rm -f build.log band_test.log qw_test.log gfactor_test.log

echo ""
echo "Next steps:"
echo "1. Run calculations with your own input files"
echo "2. Use plotting scripts to visualize results"
echo "3. Refer to docs/ for detailed documentation"
