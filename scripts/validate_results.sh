#!/bin/bash
# Result Validation Script for 8bandkp-fdm-ai
# Purpose: Automated validation of calculation results against physical criteria
# Date: 2025-01-27

set -e  # Exit on any error

echo "=========================================="
echo "8bandkp-fdm-ai Result Validation Script"
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
    fi
}

# Function to print warning
print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Sanitize example file: remove comments and inline notes
sanitize_example() {
    local src=$1
    local dest=$2
    awk '
        BEGIN{FS=""}
        {
            line=$0
            # Remove full-line comments
            if (line ~ /^\s*#/) next
            # Strip inline comments starting with #
            sub(/#.*/, "", line)
            # Trim whitespace
            sub(/^\s+/, "", line); sub(/\s+$/, "", line)
            if (length(line)>0) print line
        }
    ' "$src" > "$dest"
}

# Function to validate eigenvalues file
validate_eigenvalues() {
    local file=$1
    local test_name=$2
    
    echo "Validating $test_name..."
    
    if [ ! -f "$file" ]; then
        print_status 1 "eigenvalues.dat not found"
        return 1
    fi
    
    # Check file format
    local lines=$(wc -l < "$file")
    if [ $lines -lt 2 ]; then
        print_status 1 "eigenvalues.dat has insufficient data"
        return 1
    fi
    
    # Check for reasonable energy values
    local min_energy=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' "$file" | sort -n | head -1)
    local max_energy=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' "$file" | sort -n | tail -1)
    
    # Check for unphysical values
    if (( $(echo "$min_energy < -10" | bc -l) )); then
        print_warning "Very low energy values detected: $min_energy eV"
    fi
    
    if (( $(echo "$max_energy > 10" | bc -l) )); then
        print_warning "Very high energy values detected: $max_energy eV"
    fi
    
    # Check for NaN or Inf values
    if grep -q "nan\|inf\|NaN\|Inf" "$file"; then
        print_status 1 "NaN or Inf values found in eigenvalues.dat"
        return 1
    fi
    
    print_status 0 "eigenvalues.dat format and values look reasonable"
    echo "  - Energy range: $min_energy to $max_energy eV"
    echo "  - Data points: $lines lines"
    
    return 0
}

# Function to validate quantum well behavior
validate_quantum_well() {
    local file=$1
    
    echo "Validating quantum well behavior..."
    
    # Check for quantized levels (should have discrete energy levels)
    local unique_energies=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' "$file" | sort -n | uniq | wc -l)
    local total_energies=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="") print $i}' "$file" | wc -l)
    
    if [ $unique_energies -lt $((total_energies / 2)) ]; then
        print_warning "Limited energy level diversity - may indicate insufficient discretization"
    else
        print_status 0 "Good energy level diversity: $unique_energies unique levels"
    fi
    
    # Check for band gap
    local valence_max=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="" && $i<0) print $i}' "$file" | sort -n | tail -1)
    local conduction_min=$(awk 'NR>1 {for(i=2;i<=NF;i++) if($i!="" && $i>0) print $i}' "$file" | sort -n | head -1)
    
    if [ -n "$valence_max" ] && [ -n "$conduction_min" ]; then
        local band_gap=$(echo "$conduction_min - $valence_max" | bc -l)
        echo "  - Band gap: $band_gap eV"
        
        if (( $(echo "$band_gap > 0.05 && $band_gap < 2.0" | bc -l) )); then
            print_status 0 "Band gap in reasonable range"
        else
            print_warning "Band gap may be outside expected range"
        fi
    fi
}

# Function to validate bulk behavior
validate_bulk() {
    local file=$1
    
    echo "Validating bulk semiconductor behavior..."
    
    # Check for proper band structure at k=0
    local k0_line=$(head -2 "$file" | tail -1)
    local k0_energies=($(echo "$k0_line" | awk '{for(i=2;i<=NF;i++) print $i}'))
    
    # Find valence band maximum and conduction band minimum
    local valence_max=""
    local conduction_min=""
    
    for energy in "${k0_energies[@]}"; do
        if (( $(echo "$energy < 0" | bc -l) )); then
            if [ -z "$valence_max" ] || (( $(echo "$energy > $valence_max" | bc -l) )); then
                valence_max=$energy
            fi
        elif (( $(echo "$energy > 0" | bc -l) )); then
            if [ -z "$conduction_min" ] || (( $(echo "$energy < $conduction_min" | bc -l) )); then
                conduction_min=$energy
            fi
        fi
    done
    
    if [ -n "$valence_max" ] && [ -n "$conduction_min" ]; then
        local band_gap=$(echo "$conduction_min - $valence_max" | bc -l)
        echo "  - Band gap at k=0: $band_gap eV"
        
        if (( $(echo "$band_gap > 0.05 && $band_gap < 2.0" | bc -l) )); then
            print_status 0 "Bulk band gap in reasonable range"
        else
            print_warning "Bulk band gap may be outside expected range"
        fi
    fi
}

echo "Step 1: Running validation test cases..."

# Outputs base
RUN_ID=$(date +%Y%m%d-%H%M%S)
OUT_BASE="outputs/${RUN_ID}"
mkdir -p "${OUT_BASE}"

# Test bulk calculation
echo "Running bulk InAs60Sb40 validation..."
mkdir -p "${OUT_BASE}/bulk"
sanitize_example examples/validation_bulk_InAs60Sb40.example "${OUT_BASE}/bulk/input_clean.cfg"
if ( ./bandStructure "${OUT_BASE}/bulk/input_clean.cfg" ) > "${OUT_BASE}/bulk/bulk_validation.log" 2>&1; then
    print_status 0 "Bulk calculation completed"
    BULK_RUN_DIR=$(grep -l "input_file=${OUT_BASE}/bulk/input_clean.cfg" outputs/*/run.meta 2>/dev/null | xargs -r dirname | tail -n 1)
    validate_eigenvalues "${BULK_RUN_DIR}/eigenvalues.dat" "Bulk InAs60Sb40"
    validate_bulk "${BULK_RUN_DIR}/eigenvalues.dat"
else
    print_status 1 "Bulk calculation failed"
fi

# Test quantum well calculation
echo ""
echo "Running quantum well GaSb/InAs/AlSb validation..."
mkdir -p "${OUT_BASE}/qw"
sanitize_example examples/validation_quantum_well_GaSb_InAs_AlSb.example "${OUT_BASE}/qw/input_clean.cfg"
if ( ./bandStructure "${OUT_BASE}/qw/input_clean.cfg" ) > "${OUT_BASE}/qw/qw_validation.log" 2>&1; then
    print_status 0 "Quantum well calculation completed"
    QW_RUN_DIR=$(grep -l "input_file=${OUT_BASE}/qw/input_clean.cfg" outputs/*/run.meta 2>/dev/null | xargs -r dirname | tail -n 1)
    validate_eigenvalues "${QW_RUN_DIR}/eigenvalues.dat" "Quantum Well GaSb/InAs/AlSb"
    validate_quantum_well "${QW_RUN_DIR}/eigenvalues.dat"
else
    print_status 1 "Quantum well calculation failed"
fi

echo ""
echo "Step 2: Checking output file completeness..."

# Check for required output files
if [ -n "$BULK_RUN_DIR" ] && [ -f "$BULK_RUN_DIR/eigenvalues.dat" ]; then
    print_status 0 "bulk eigenvalues.dat generated"
else
    print_status 1 "bulk eigenvalues.dat missing"
fi

if [ -n "$QW_RUN_DIR" ] && [ -f "$QW_RUN_DIR/parts.dat" ]; then
    print_status 0 "qw parts.dat generated"
else
    print_warning "parts.dat missing"
fi

# Check for eigenfunction files
eigenfunction_count=$(ls "$QW_RUN_DIR"/eigenfunctions_k_*.dat 2>/dev/null | wc -l)
if [ $eigenfunction_count -gt 0 ]; then
    print_status 0 "$eigenfunction_count eigenfunction files generated"
else
    print_warning "No eigenfunction files generated"
fi

echo ""
echo "Step 3: Performance validation..."

# Check calculation time
if [ -f "bulk_validation.log" ]; then
    echo "Bulk calculation completed successfully"
fi

if [ -f "qw_validation.log" ]; then
    echo "Quantum well calculation completed successfully"
fi

echo ""
echo "=========================================="
echo "Result validation completed!"
echo "=========================================="

# Clean up log files
rm -f bulk_validation.log qw_validation.log

echo ""
echo "Validation summary:"
echo "- Check eigenvalues.dat for physical behavior"
echo "- Verify band gaps are in reasonable range"
echo "- Confirm quantum well effects (if applicable)"
echo "- Review any warnings or errors above"
