#!/bin/bash
# Integration Test: Bulk G-Factor Numerical Method
# Tests numerical g-factor calculation for bulk materials

set -euo pipefail

TEST_OUTPUT_DIR="${1:-outputs/test-integration}"
mkdir -p "${TEST_OUTPUT_DIR}/gfactor_num"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

TESTS_PASSED=0
TESTS_FAILED=0

print_status() {
    if [ "$1" -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗${NC} $2"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Validate g-factor values
validate_gfactor() {
    local gfactor_file="$1"
    local expected_min="$2"
    local expected_max="$3"
    local material_name="$4"
    
    if [ ! -f "$gfactor_file" ]; then
        echo "  ERROR: gfactor.dat not found"
        return 1
    fi
    
    # Read gx, gy, gz
    read gx gy gz < "$gfactor_file"
    
    echo "  g-factors: gx=$gx, gy=$gy, gz=$gz"
    
    # Check for NaN
    if [[ "$gx" == *"nan"* ]] || [[ "$gy" == *"nan"* ]] || [[ "$gz" == *"nan"* ]]; then
        echo "  ERROR: NaN values detected"
        return 1
    fi
    
    # Check sign consistency (for cubic bulk: gx = gy = gz)
    local gx_abs=$(echo "$gx" | awk '{print ($1<0)?-$1:$1}')
    local gy_abs=$(echo "$gy" | awk '{print ($1<0)?-$1:$1}')
    local gz_abs=$(echo "$gz" | awk '{print ($1<0)?-$1:$1}')
    
    local diff_xy=$(echo "$gx_abs - $gy_abs" | bc -l | awk '{print ($1<0)?-$1:$1}')
    local diff_xz=$(echo "$gx_abs - $gz_abs" | bc -l | awk '{print ($1<0)?-$1:$1}')
    
    # For cubic bulk, should be isotropic (within 1%)
    if (( $(echo "$diff_xy < 0.01 * $gx_abs && $diff_xz < 0.01 * $gx_abs" | bc -l) )); then
        echo "  Isotropy check: PASS (cubic symmetry)"
    else
        echo "  WARNING: Anisotropy detected (unexpected for bulk)"
    fi
    
    # Check if in expected range
    local g_avg=$(echo "($gx_abs + $gy_abs + $gz_abs) / 3" | bc -l)
    
    if (( $(echo "$g_avg >= $expected_min && $g_avg <= $expected_max" | bc -l) )); then
        echo "  Magnitude check: PASS"
        return 0
    else
        echo "  WARNING: |g| = $g_avg outside expected range [$expected_min, $expected_max]"
        return 1
    fi
}

echo "Testing bulk g-factor numerical calculations..."
echo ""

# Test 1: GaAs CB (g ~ -0.31 to -0.32)
echo "Test 1: Bulk GaAs CB (numerical)"
if ./gfactorCalculation examples/gfactor_bulk_GaAs_numerical.example > "${TEST_OUTPUT_DIR}/gfactor_num/gaas.log" 2>&1; then
    if validate_gfactor "gfactor.dat" "0.25" "0.40" "GaAs"; then
        print_status 0 "GaAs CB g-factor valid"
    else
        print_status 1 "GaAs CB g-factor out of range"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_num/gaas_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "GaAs CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 2: InAs CB (g ~ -14.6)
echo "Test 2: Bulk InAs CB (numerical)"
if ./gfactorCalculation examples/gfactor_bulk_InAs_numerical.example > "${TEST_OUTPUT_DIR}/gfactor_num/inas.log" 2>&1; then
    if validate_gfactor "gfactor.dat" "14.0" "15.5" "InAs"; then
        print_status 0 "InAs CB g-factor valid"
    else
        print_status 1 "InAs CB g-factor out of range"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_num/inas_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "InAs CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 3: GaSb CB (g ~ -8.7)
echo "Test 3: Bulk GaSb CB (numerical)"
if ./gfactorCalculation examples/gfactor_bulk_GaSb_numerical.example > "${TEST_OUTPUT_DIR}/gfactor_num/gasb.log" 2>&1; then
    if validate_gfactor "gfactor.dat" "8.0" "9.5" "GaSb"; then
        print_status 0 "GaSb CB g-factor valid"
    else
        print_status 1 "GaSb CB g-factor out of range"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_num/gasb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "GaSb CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 4: InSb CB (g ~ -49 to -51)
echo "Test 4: Bulk InSb CB (numerical)"
if ./gfactorCalculation examples/gfactor_bulk_InSb_numerical.example > "${TEST_OUTPUT_DIR}/gfactor_num/insb.log" 2>&1; then
    if validate_gfactor "gfactor.dat" "45.0" "55.0" "InSb"; then
        print_status 0 "InSb CB g-factor valid"
    else
        print_status 1 "InSb CB g-factor out of range"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_num/insb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "InSb CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Final cleanup
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* *.log 2>/dev/null || true

echo ""
echo "Bulk g-factor numerical tests: $TESTS_PASSED passed, $TESTS_FAILED failed"

if [ $TESTS_FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
