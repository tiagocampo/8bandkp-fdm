#!/bin/bash
# Integration Test: Bulk Band Structure Calculations
# Tests bulk material band structure for multiple materials

set -euo pipefail

TEST_OUTPUT_DIR="${1:-outputs/test-integration}"
mkdir -p "${TEST_OUTPUT_DIR}/bulk"

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

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Validate band gap is in expected range
validate_band_gap() {
    local eigenvalues_file="$1"
    local expected_min="$2"
    local expected_max="$3"
    local material_name="$4"
    
    if [ ! -f "$eigenvalues_file" ]; then
        return 1
    fi
    
    # Find valence band max and conduction band min at k=0
    local k0_line=$(head -2 "$eigenvalues_file" | tail -1)
    local valence_max=$(echo "$k0_line" | awk '{for(i=2;i<=NF;i++) if($i<=0) print $i}' | sort -n | tail -1)
    local conduction_min=$(echo "$k0_line" | awk '{for(i=2;i<=NF;i++) if($i>0) print $i}' | sort -n | head -1)
    
    if [ -z "$valence_max" ] || [ -z "$conduction_min" ]; then
        print_warning "Could not find band gap for $material_name"
        return 1
    fi
    
    local band_gap=$(echo "$conduction_min - $valence_max" | bc -l)
    
    echo "  Band gap: $band_gap eV (expected: $expected_min - $expected_max eV)"
    
    # Check if in range
    if (( $(echo "$band_gap >= $expected_min && $band_gap <= $expected_max" | bc -l) )); then
        return 0
    else
        return 1
    fi
}

echo "Testing bulk band structure calculations..."
echo ""

# Test 1: GaAs (Eg ~ 1.42 eV)
echo "Test 1: Bulk GaAs"
if ./bandStructure examples/bulk_GaAs.example --out . > "${TEST_OUTPUT_DIR}/bulk/gaas.log" 2>&1; then
    if validate_band_gap "eigenvalues.dat" "1.3" "1.6" "GaAs"; then
        print_status 0 "GaAs band structure valid"
    else
        print_status 1 "GaAs band gap out of range"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/bulk/gaas_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "GaAs calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 2: InAs (Eg ~ 0.35 eV)
echo "Test 2: Bulk InAs"
if ./bandStructure examples/validation_bulk_InAs.example --out . > "${TEST_OUTPUT_DIR}/bulk/inas.log" 2>&1; then
    if validate_band_gap "eigenvalues.dat" "0.25" "0.45" "InAs"; then
        print_status 0 "InAs band structure valid"
    else
        print_status 1 "InAs band gap out of range"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/bulk/inas_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "InAs calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 3: GaSb (Eg ~ 0.72 eV)
echo "Test 3: Bulk GaSb"
if ./bandStructure examples/validation_bulk_GaSb.example --out . > "${TEST_OUTPUT_DIR}/bulk/gasb.log" 2>&1; then
    if validate_band_gap "eigenvalues.dat" "0.6" "0.85" "GaSb"; then
        print_status 0 "GaSb band structure valid"
    else
        print_status 1 "GaSb band gap out of range"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/bulk/gasb_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "GaSb calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 4: InSb (Eg ~ 0.17 eV)
echo "Test 4: Bulk InSb"
if ./bandStructure examples/validation_bulk_InSb.example --out . > "${TEST_OUTPUT_DIR}/bulk/insb.log" 2>&1; then
    if validate_band_gap "eigenvalues.dat" "0.10" "0.25" "InSb"; then
        print_status 0 "InSb band structure valid"
    else
        print_status 1 "InSb band gap out of range"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/bulk/insb_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "InSb calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Test 5: InAs60Sb40 (Eg ~ 0.11 eV)
echo "Test 5: Bulk InAs60Sb40"
if ./bandStructure examples/validation_bulk_InAs60Sb40.example --out . > "${TEST_OUTPUT_DIR}/bulk/inassb.log" 2>&1; then
    if validate_band_gap "eigenvalues.dat" "0.05" "0.20" "InAs60Sb40"; then
        print_status 0 "InAs60Sb40 band structure valid"
    else
        print_status 1 "InAs60Sb40 band gap out of range"
    fi
    mv eigenvalues.dat "${TEST_OUTPUT_DIR}/bulk/inassb_eigenvalues.dat" 2>/dev/null || true
else
    print_status 1 "InAs60Sb40 calculation failed"
fi
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* 2>/dev/null || true

# Final cleanup
rm -f eigenvalues.dat parts.dat eigenfunctions_k_*.dat fort.* *.log 2>/dev/null || true

echo ""
echo "Bulk band structure tests: $TESTS_PASSED passed, $TESTS_FAILED failed"

if [ $TESTS_FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
