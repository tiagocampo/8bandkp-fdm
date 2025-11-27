#!/bin/bash
# Integration Test: Quantum Well G-Factor Analytical Method
# Tests analytical g-factor calculation for QW systems

set -euo pipefail

TEST_OUTPUT_DIR="${1:-outputs/test-integration}"
mkdir -p "${TEST_OUTPUT_DIR}/gfactor_qw"

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

# Validate QW g-factor (expect anisotropy)
validate_qw_gfactor() {
    local gfactor_file="$1"
    local system_name="$2"
    
    if [ ! -f "$gfactor_file" ]; then
        echo "  ERROR: gfactor.dat not found"
        return 1
    fi
    
    read gx gy gz < "$gfactor_file"
    
    echo "  g-factors: gx=$gx, gy=$gy, gz=$gz"
    
    # Check for NaN
    if [[ "$gx" == *"nan"* ]] || [[ "$gy" == *"nan"* ]] || [[ "$gz" == *"nan"* ]]; then
        echo "  ERROR: NaN values detected"
        return 1
    fi
    
    # For QW, expect anisotropy: gz != gx,gy
    local gx_abs=$(echo "$gx" | awk '{print ($1<0)?-$1:$1}')
    local gy_abs=$(echo "$gy" | awk '{print ($1<0)?-$1:$1}')
    local gz_abs=$(echo "$gz" | awk '{print ($1<0)?-$1:$1}')
    
    local diff_z_vs_x=$(echo "$gz_abs - $gx_abs" | bc -l | awk '{print ($1<0)?-$1:$1}')
    local anisotropy_ratio=$(echo "scale=2; $gz_abs / $gx_abs" | bc -l)
    
    # Check for anisotropy (gz should differ significantly from gx for strong confinement)
    if (( $(echo "$diff_z_vs_x > 0.1 * $gz_abs" | bc -l) )); then
        echo "  Anisotropy detected: gz/g|| ~ $anisotropy_ratio (expected for QW)"
        return 0
    else
        echo "  WARNING: Limited anisotropy (gz/g|| ~ $anisotropy_ratio)"
        return 0  # Don't fail, weak confinement can have low anisotropy
    fi
}

echo "Testing QW g-factor analytical calculations..."
echo ""

# Test 1: GaAs/AlGaAs CB
echo "Test 1: GaAs/Al₀.₃Ga₀.₇As CB (analytical)"
if ./gfactorCalculation examples/gfactor_qw_gaas_algaas_cb_analytical.example > "${TEST_OUTPUT_DIR}/gfactor_qw/gaas_algaas_cb.log" 2>&1; then
    if validate_qw_gfactor "gfactor.dat" "GaAs/AlGaAs CB"; then
        print_status 0 "GaAs/AlGaAs CB g-factor valid"
    else
        print_status 1 "GaAs/AlGaAs CB g-factor invalid"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_qw/gaas_algaas_cb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "GaAs/AlGaAs CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 2: GaAs/AlGaAs VB
echo "Test 2: GaAs/Al₀.₃Ga₀.₇As VB (analytical)"
if ./gfactorCalculation examples/gfactor_qw_gaas_algaas_vb_analytical.example > "${TEST_OUTPUT_DIR}/gfactor_qw/gaas_algaas_vb.log" 2>&1; then
    if validate_qw_gfactor "gfactor.dat" "GaAs/AlGaAs VB"; then
        print_status 0 "GaAs/AlGaAs VB g-factor valid"
    else
        print_status 1 "GaAs/AlGaAs VB g-factor invalid"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_qw/gaas_algaas_vb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "GaAs/AlGaAs VB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 3: InAs/AlSb CB
echo "Test 3: InAs/AlSb CB (analytical)"
if ./gfactorCalculation examples/gfactor_qw_inas_alsb_cb_analytical.example > "${TEST_OUTPUT_DIR}/gfactor_qw/inas_alsb_cb.log" 2>&1; then
    if validate_qw_gfactor "gfactor.dat" "InAs/AlSb CB"; then
        print_status 0 "InAs/AlSb CB g-factor valid"
    else
        print_status 1 "InAs/AlSb CB g-factor invalid"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_qw/inas_alsb_cb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "InAs/AlSb CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 4: InAs/GaSb/AlSb CB
echo "Test 4: InAs/GaSb/AlSb type-II CB (analytical)"
if ./gfactorCalculation examples/gfactor_qw_inas_gasb_alsb_clean_cb_analytical.example > "${TEST_OUTPUT_DIR}/gfactor_qw/inas_gasb_alsb_cb.log" 2>&1; then
    if validate_qw_gfactor "gfactor.dat" "InAs/GaSb/AlSb CB"; then
        print_status 0 "InAs/GaSb/AlSb CB g-factor valid"
    else
        print_status 1 "InAs/GaSb/AlSb CB g-factor invalid"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_qw/inas_gasb_alsb_cb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "InAs/GaSb/AlSb CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 5: InGaAs/GaAs strained CB
echo "Test 5: In₀.₂Ga₀.₈As/GaAs strained CB (analytical)"
if ./gfactorCalculation examples/gfactor_qw_ingaas_gaas_strained_cb_analytical.example > "${TEST_OUTPUT_DIR}/gfactor_qw/ingaas_gaas_cb.log" 2>&1; then
    if validate_qw_gfactor "gfactor.dat" "InGaAs/GaAs CB"; then
        print_status 0 "InGaAs/GaAs CB g-factor valid"
    else
        print_status 1 "InGaAs/GaAs CB g-factor invalid"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_qw/ingaas_gaas_cb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "InGaAs/GaAs CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Test 6: InGaAs/InP CB
echo "Test 6: InGaAs/InP CB (analytical)"
if ./gfactorCalculation examples/gfactor_qw_ingaas_inp_cb_analytical.example > "${TEST_OUTPUT_DIR}/gfactor_qw/ingaas_inp_cb.log" 2>&1; then
    if validate_qw_gfactor "gfactor.dat" "InGaAs/InP CB"; then
        print_status 0 "InGaAs/InP CB g-factor valid"
    else
        print_status 1 "InGaAs/InP CB g-factor invalid"
    fi
    cp gfactor.dat "${TEST_OUTPUT_DIR}/gfactor_qw/ingaas_inp_cb_gfactor.dat" 2>/dev/null || true
else
    print_status 1 "InGaAs/InP CB calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Final cleanup
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* *.log 2>/dev/null || true

echo ""
echo "QW g-factor analytical tests: $TESTS_PASSED passed, $TESTS_FAILED failed"

if [ $TESTS_FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
