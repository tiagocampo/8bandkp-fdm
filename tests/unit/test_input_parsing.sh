#!/bin/bash
# Unit Test: Input File Parsing
# Tests that executables correctly parse various input file formats

set -euo pipefail

# Get test output directory from argument
TEST_OUTPUT_DIR="${1:-outputs/test-unit}"
mkdir -p "${TEST_OUTPUT_DIR}"

# Colors
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

echo "Testing input file parsing..."

# Test 1: Valid bulk example
echo "Test 1: Valid bulk input (GaAs)"
if ./bandStructure examples/bulk_GaAs.example > "${TEST_OUTPUT_DIR}/bulk_parse.log" 2>&1; then
    print_status 0 "Bulk input parsed successfully"
else
    print_status 1 "Failed to parse bulk input"
fi

# Test 2: Valid QW example
echo "Test 2: Valid QW input (GaAs/AlGaAs)"
if ./bandStructure examples/qw_gaas_algaas.example > "${TEST_OUTPUT_DIR}/qw_parse.log" 2>&1; then
    print_status 0 "QW input parsed successfully"
else
    print_status 1 "Failed to parse QW input"
fi

# Test 3: Valid g-factor numerical input
echo "Test 3: Valid g-factor numerical input"
if ./gfactorCalculation examples/gfactor_bulk_GaAs_numerical.example > "${TEST_OUTPUT_DIR}/gfactor_num_parse.log" 2>&1; then
    print_status 0 "G-factor numerical input parsed successfully"
else
    print_status 1 "Failed to parse g-factor numerical input"
fi

# Test 4: Valid g-factor analytical input
echo "Test 4: Valid g-factor analytical input"
if ./gfactorCalculation examples/gfactor_analytical_bulk_GaAs.example > "${TEST_OUTPUT_DIR}/gfactor_ana_parse.log" 2>&1; then
    print_status 0 "G-factor analytical input parsed successfully"
else
    print_status 1 "Failed to parse g-factor analytical input"
fi

# Test 5: Invalid input should fail gracefully
echo "Test 5: Invalid input handling"
echo "invalid_parameter: bad_value" > "${TEST_OUTPUT_DIR}/invalid.example"
if ./bandStructure "${TEST_OUTPUT_DIR}/invalid.example" > "${TEST_OUTPUT_DIR}/invalid_parse.log" 2>&1; then
    print_status 1 "Should have failed on invalid input"
else
    print_status 0 "Correctly rejected invalid input"
fi

# Clean up temporary files
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat gfactor.dat fort.* *.log 2>/dev/null || true

echo ""
echo "Input parsing tests completed"
exit 0
