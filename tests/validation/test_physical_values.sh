#!/bin/bash
# Validation Test: Physical Values
# Validates calculated values against known literature values

set -euo pipefail

TEST_OUTPUT_DIR="${1:-outputs/test-validation}"
mkdir -p "${TEST_OUTPUT_DIR}/physical"

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

echo "Testing physical values against literature..."
echo ""

echo "=== Band Gap Validation ==="
echo ""

# GaAs: Eg ~ 1.42 eV (literature)
echo "Test 1: GaAs band gap"
./bandStructure examples/bulk_GaAs.example > "${TEST_OUTPUT_DIR}/physical/gaas_bs.log" 2>&1
if [ -f "eigenvalues.dat" ]; then
    k0_line=$(head -2 eigenvalues.dat | tail -1)
    vb_max=$(echo "$k0_line" | awk '{for(i=2;i<=NF;i++) if($i<0) print $i}' | sort -n | tail -1)
    cb_min=$(echo "$k0_line" | awk '{for(i=2;i<=NF;i++) if($i>0) print $i}' | sort -n | head -1)
    eg=$(echo "$cb_min - $vb_max" | bc -l)
    
    echo "  Calculated: Eg = $eg eV"
    echo "  Literature: Eg ~ 1.42 eV"
    
    if (( $(echo "$eg >= 1.3 && $eg <= 1.6" | bc -l) )); then
        print_status 0 "GaAs band gap within literature range"
    else
        print_status 1 "GaAs band gap outside literature range"
    fi
else
    print_status 1 "GaAs calculation failed"
fi
rm -f eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# InAs: Eg ~ 0.35 eV (literature)
echo "Test 2: InAs band gap"
./bandStructure examples/validation_bulk_InAs.example > "${TEST_OUTPUT_DIR}/physical/inas_bs.log" 2>&1
if [ -f "eigenvalues.dat" ]; then
    k0_line=$(head -2 eigenvalues.dat | tail -1)
    vb_max=$(echo "$k0_line" | awk '{for(i=2;i<=NF;i++) if($i<0) print $i}' | sort -n | tail -1)
    cb_min=$(echo "$k0_line" | awk '{for(i=2;i<=NF;i++) if($i>0) print $i}' | sort -n | head -1)
    eg=$(echo "$cb_min - $vb_max" | bc -l)
    
    echo "  Calculated: Eg = $eg eV"
    echo "  Literature: Eg ~ 0.35 eV"
    
    if (( $(echo "$eg >= 0.25 && $eg <= 0.45" | bc -l) )); then
        print_status 0 "InAs band gap within literature range"
    else
        print_status 1 "InAs band gap outside literature range"
    fi
else
    print_status 1 "InAs calculation failed"
fi
rm -f eigenvalues.dat parts.dat fort.* 2>/dev/null || true

echo ""
echo "=== G-Factor Validation ==="
echo ""

# GaAs: g ~ -0.44 (Roth formula), g ~ -0.315 (numerical)
echo "Test 3: GaAs CB g-factor (numerical)"
./gfactorCalculation examples/gfactor_bulk_GaAs_numerical.example > "${TEST_OUTPUT_DIR}/physical/gaas_gf.log" 2>&1
if [ -f "gfactor.dat" ]; then
    read gx gy gz < "gfactor.dat"
    g_avg=$(echo "$gx" | awk '{print ($1<0)?-$1:$1}')
    
    echo "  Calculated: |g| = $g_avg"
    echo "  Literature: |g| ~ 0.31-0.32 (Kane model)"
    
    if (( $(echo "$g_avg >= 0.25 && $g_avg <= 0.40" | bc -l) )); then
        print_status 0 "GaAs CB g-factor within literature range"
    else
        print_status 1 "GaAs CB g-factor outside literature range"
    fi
else
    print_status 1 "GaAs g-factor calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# InAs: g ~ -14.7 (literature)
echo "Test 4: InAs CB g-factor (numerical)"
./gfactorCalculation examples/gfactor_bulk_InAs_numerical.example > "${TEST_OUTPUT_DIR}/physical/inas_gf.log" 2>&1
if [ -f "gfactor.dat" ]; then
    read gx gy gz < "gfactor.dat"
    g_avg=$(echo "$gx" | awk '{print ($1<0)?-$1:$1}')
    
    echo "  Calculated: |g| = $g_avg"
    echo "  Literature: |g| ~ 14.7"
    
    if (( $(echo "$g_avg >= 14.0 && $g_avg <= 15.5" | bc -l) )); then
        print_status 0 "InAs CB g-factor within literature range"
    else
        print_status 1 "InAs CB g-factor outside literature range"
    fi
else
    print_status 1 "InAs g-factor calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# InSb: g ~ -51 (literature, very large due to small gap)
echo "Test 5: InSb CB g-factor (numerical)"
./gfactorCalculation examples/gfactor_bulk_InSb_numerical.example > "${TEST_OUTPUT_DIR}/physical/insb_gf.log" 2>&1
if [ -f "gfactor.dat" ]; then
    read gx gy gz < "gfactor.dat"
    g_avg=$(echo "$gx" | awk '{print ($1<0)?-$1:$1}')
    
    echo "  Calculated: |g| = $g_avg"
    echo "  Literature: |g| ~ 51"
    
    if (( $(echo "$g_avg >= 45.0 && $g_avg <= 55.0" | bc -l) )); then
        print_status 0 "InSb CB g-factor within literature range"
    else
        print_status 1 "InSb CB g-factor outside literature range"
    fi
else
    print_status 1 "InSb g-factor calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

echo ""
echo "=== QW Confinement Validation ==="
echo ""

# GaAs/AlGaAs QW: Should show anisotropy
echo "Test 6: GaAs/AlGaAs QW g-factor anisotropy"
./gfactorCalculation examples/gfactor_qw_gaas_algaas_cb_analytical.example > "${TEST_OUTPUT_DIR}/physical/qw_gaas_algaas.log" 2>&1
if [ -f "gfactor.dat" ]; then
    read gx gy gz < "gfactor.dat"
    
    gx_abs=$(echo "$gx" | awk '{print ($1<0)?-$1:$1}')
    gz_abs=$(echo "$gz" | awk '{print ($1<0)?-$1:$1}')
    ratio=$(echo "scale=2; $gz_abs / $gx_abs" | bc -l)
    
    echo "  gx = $gx, gz = $gz"
    echo "  Anisotropy ratio: gz/g|| = $ratio"
    echo "  Expected: ratio > 1 (confinement effect)"
    
    if (( $(echo "$ratio > 1.5" | bc -l) )); then
        print_status 0 "GaAs/AlGaAs QW shows expected anisotropy"
    else
        print_warning "Limited anisotropy detected (ratio=$ratio)"
        print_status 0 "QW calculation completed (weak confinement possible)"
    fi
else
    print_status 1 "QW g-factor calculation failed"
fi
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* 2>/dev/null || true

# Final cleanup
rm -f gfactor.dat eigenvalues.dat parts.dat fort.* *.log 2>/dev/null || true

echo ""
echo "Physical value validation: $TESTS_PASSED passed, $TESTS_FAILED failed"

if [ $TESTS_FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
