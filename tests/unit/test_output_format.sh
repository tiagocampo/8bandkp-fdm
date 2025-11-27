#!/bin/bash
# Unit Test: Output File Format Validation
# Tests that output files have correct format and no invalid values

set -euo pipefail

TEST_OUTPUT_DIR="${1:-outputs/test-unit}"
mkdir -p "${TEST_OUTPUT_DIR}"

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

print_status() {
    if [ "$1" -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
    else
        echo -e "${RED}✗${NC} $2"
        return 1
    fi
}

echo "Testing output file formats..."

# Run a simple bulk calculation to generate outputs
echo "Generating test outputs..."
./bandStructure examples/bulk_GaAs.example --out . > "${TEST_OUTPUT_DIR}/band_run.log" 2>&1

# Test 1: eigenvalues.dat exists and has content
echo "Test 1: eigenvalues.dat format"
if [ -f "eigenvalues.dat" ] && [ -s "eigenvalues.dat" ]; then
    print_status 0 "eigenvalues.dat exists"
else
    print_status 1 "eigenvalues.dat missing or empty"
fi

# Test 2: Check for NaN or Inf values in eigenvalues.dat
if [ -f "eigenvalues.dat" ]; then
    if grep -qi "nan\|inf" eigenvalues.dat; then
        print_status 1 "NaN or Inf values found in eigenvalues.dat"
    else
        print_status 0 "No NaN or Inf values in eigenvalues.dat"
    fi
fi

# Test 3: Check eigenvalues.dat has proper columns
if [ -f "eigenvalues.dat" ]; then
    header_cols=$(head -1 eigenvalues.dat | wc -w)
    data_cols=$(head -2 eigenvalues.dat | tail -1 | wc -w)
    if [ "$header_cols" -eq "$data_cols" ]; then
        print_status 0 "eigenvalues.dat has consistent column count ($data_cols)"
    else
        print_status 1 "Inconsistent column count in eigenvalues.dat"
    fi
fi

# Test 4: Energy values in reasonable range
if [ -f "eigenvalues.dat" ]; then
    min_energy=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' eigenvalues.dat | sort -n | head -1)
    max_energy=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' eigenvalues.dat | sort -n | tail -1)
    
    # Check if energies are in -20 to +20 eV range (very generous)
    if (( $(echo "$min_energy > -20" | bc -l) )) && (( $(echo "$max_energy < 20" | bc -l) )); then
        print_status 0 "Energy values in reasonable range ($min_energy to $max_energy eV)"
    else
        print_status 1 "Energy values outside reasonable range"
    fi
fi

# Test 5: Run g-factor calculation and check gfactor.dat format
echo "Test 5: gfactor.dat format"
./gfactorCalculation examples/gfactor_bulk_GaAs_numerical.example > "${TEST_OUTPUT_DIR}/gfactor_run.log" 2>&1

if [ -f "gfactor.dat" ]; then
    # Should have exactly 3 values: gx gy gz
    value_count=$(wc -w < gfactor.dat)
    if [ "$value_count" -eq 3 ]; then
        print_status 0 "gfactor.dat has correct format (3 values)"
    else
        print_status 1 "gfactor.dat has wrong number of values ($value_count, expected 3)"
    fi
    
    # Check for NaN or Inf
    if grep -qi "nan\|inf" gfactor.dat; then
        print_status 1 "NaN or Inf values found in gfactor.dat"
    else
        print_status 0 "No NaN or Inf values in gfactor.dat"
    fi
else
    print_status 1 "gfactor.dat not created"
fi

# Clean up
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat gfactor.dat fort.* *.log 2>/dev/null || true

echo ""
echo "Output format tests completed"
exit 0
