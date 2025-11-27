#!/bin/bash
# Integration Test: Quantum Well Band Structure Calculations
# Tests QW band structure for multiple systems

set -euo pipefail

TEST_OUTPUT_DIR="${1:-outputs/test-integration}"
mkdir -p "${TEST_OUTPUT_DIR}/qw"

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

# Validate QW behavior
validate_qw_structure() {
    local eigenvalues_file="$1"
    local system_name="$2"
    
    if [ ! -f "$eigenvalues_file" ]; then
        return 1
    fi
    
    # Check for NaN or Inf
    if grep -qi "nan\|inf" "$eigenvalues_file"; then
        echo "  ERROR: NaN or Inf values found"
        return 1
    fi
    
    # Check for quantized levels (multiple unique energy values)
    local unique_energies=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' "$eigenvalues_file" | sort -n | uniq | wc -l)
    
    if [ "$unique_energies" -lt 5 ]; then
        echo "  WARNING: Limited energy level diversity ($unique_energies levels)"
    else
        echo "  Energy level diversity: $unique_energies levels"
    fi
    
    # Check energy range is reasonable
    local min_energy=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' "$eigenvalues_file" | sort -n | head -1)
    local max_energy=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' "$eigenvalues_file" | sort -n | tail -1)
    
    echo "  Energy range: $min_energy to $max_energy eV"
    
    if (( $(echo "$min_energy < -10 || $max_energy > 10" | bc -l) )); then
        echo "  WARNING: Energy range seems very large"
    fi
    
    return 0
}

echo "Testing quantum well band structure calculations..."
echo ""

# Test 1: GaAs/AlGaAs QW
echo "Test 1: GaAs/Al₀.₃Ga₀.₇As QW"
if ./bandStructure examples/qw_gaas_algaas.example > "${TEST_OUTPUT_DIR}/qw/gaas_algaas.log" 2>&1; then
    if validate_qw_structure "eigenvalues.dat" "GaAs/AlGaAs"; then
        print_status 0 "GaAs/AlGaAs QW structure valid"
    else
        print_status 1 "GaAs/AlGaAs QW validation failed"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/qw/gaas_algaas_eigenvalues.dat" 2>/dev/null || true
    mv parts.dat "${TEST_OUTPUT_DIR}/qw/gaas_algaas_parts.dat" 2>/dev/null || true
else
    print_status 1 "GaAs/AlGaAs QW calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 2: InAs/AlSb QW
echo "Test 2: InAs/AlSb QW"
if ./bandStructure examples/qw_inas_alsb.example > "${TEST_OUTPUT_DIR}/qw/inas_alsb.log" 2>&1; then
    if validate_qw_structure "eigenvalues.dat" "InAs/AlSb"; then
        print_status 0 "InAs/AlSb QW structure valid"
    else
        print_status 1 "InAs/AlSb QW validation failed"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/qw/inas_alsb_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "InAs/AlSb QW calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 3: InAs/GaSb/AlSb type-II QW
echo "Test 3: InAs/GaSb/AlSb type-II QW"
if ./bandStructure examples/qw_inas_gasb_alsb_clean.example > "${TEST_OUTPUT_DIR}/qw/inas_gasb_alsb.log" 2>&1; then
    if validate_qw_structure "eigenvalues.dat" "InAs/GaSb/AlSb"; then
        print_status 0 "InAs/GaSb/AlSb QW structure valid"
    else
        print_status 1 "InAs/GaSb/AlSb QW validation failed"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/qw/inas_gasb_alsb_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "InAs/GaSb/AlSb QW calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 4: Strained InGaAs/GaAs QW
echo "Test 4: Strained In₀.₂Ga₀.₈As/GaAs QW"
if ./bandStructure examples/qw_ingaas_gaas_strained.example > "${TEST_OUTPUT_DIR}/qw/ingaas_gaas.log" 2>&1; then
    if validate_qw_structure "eigenvalues.dat" "InGaAs/GaAs"; then
        print_status 0 "InGaAs/GaAs QW structure valid"
    else
        print_status 1 "InGaAs/GaAs QW validation failed"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/qw/ingaas_gaas_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "InGaAs/GaAs QW calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 5: InGaAs/InP QW
echo "Test 5: InGaAs/InP QW"
if ./bandStructure examples/qw_ingaas_inp.example > "${TEST_OUTPUT_DIR}/qw/ingaas_inp.log" 2>&1; then
    if validate_qw_structure "eigenvalues.dat" "InGaAs/InP"; then
        print_status 0 "InGaAs/InP QW structure valid"
    else
        print_status 1 "InGaAs/InP QW validation failed"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/qw/ingaas_inp_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "InGaAs/InP QW calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Final cleanup
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* *.log 2>/dev/null || true

echo ""
echo "QW band structure tests: $TESTS_PASSED passed, $TESTS_FAILED failed"

if [ $TESTS_FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
