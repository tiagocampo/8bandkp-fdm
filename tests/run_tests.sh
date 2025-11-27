#!/bin/bash
# Main Test Runner for 8bandkp-fdm-ai
# Purpose: Discover and run all test scripts, generate summary report
# Usage: ./run_tests.sh [category] [options]
#   category: all (default), unit, integration, validation, quick
#   options: --verbose, --parallel, --stop-on-fail

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_SKIPPED=0

# Configuration
VERBOSE=false
PARALLEL=false
STOP_ON_FAIL=false
TEST_CATEGORY="all"
TEST_OUTPUT_DIR=""

# Print functions
print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_status() {
    if [ "$1" -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
    else
        echo -e "${RED}✗${NC} $2"
    fi
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

print_info() {
    echo -e "${BLUE}ℹ${NC} $1"
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            unit|integration|validation|all|quick)
                TEST_CATEGORY="$1"
                shift
                ;;
            --verbose|-v)
                VERBOSE=true
                shift
                ;;
            --parallel|-p)
                PARALLEL=true
                shift
                ;;
            --stop-on-fail|-s)
                STOP_ON_FAIL=true
                shift
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

show_help() {
    cat <<EOF
Usage: ./run_tests.sh [category] [options]

Categories:
  all           Run all tests (default)
  unit          Run only unit tests
  integration   Run only integration tests
  validation    Run only validation tests
  quick         Run quick smoke tests only

Options:
  --verbose, -v      Show detailed test output
  --parallel, -p     Run tests in parallel (experimental)
  --stop-on-fail, -s Stop on first test failure
  --help, -h         Show this help message

Examples:
  ./run_tests.sh                    # Run all tests
  ./run_tests.sh unit --verbose     # Run unit tests with verbose output
  ./run_tests.sh integration        # Run only integration tests
EOF
}

# Setup test environment
setup_test_env() {
    # Create timestamped test output directory
    local timestamp=$(date +%Y%m%d-%H%M%S)
    TEST_OUTPUT_DIR="outputs/test-${timestamp}"
    mkdir -p "${TEST_OUTPUT_DIR}"
    
    # Create subdirectories
    mkdir -p "${TEST_OUTPUT_DIR}/unit"
    mkdir -p "${TEST_OUTPUT_DIR}/integration"
    mkdir -p "${TEST_OUTPUT_DIR}/validation"
    
    print_info "Test output directory: ${TEST_OUTPUT_DIR}"
}

# Run a single test script
run_test() {
    local test_script="$1"
    local test_name=$(basename "${test_script}" .sh)
    local test_category=$(basename "$(dirname "${test_script}")")
    
    TESTS_RUN=$((TESTS_RUN + 1))
    
    echo ""
    echo "Running: ${test_name} (${test_category})"
    echo "---"
    
    local test_log="${TEST_OUTPUT_DIR}/${test_category}/${test_name}.log"
    local test_start=$(date +%s)
    
    # Run the test
    local exit_code=0
    if [ "$VERBOSE" = true ]; then
        bash "${test_script}" "${TEST_OUTPUT_DIR}/${test_category}" 2>&1 | tee "${test_log}"
        exit_code=${PIPESTATUS[0]}
    else
        bash "${test_script}" "${TEST_OUTPUT_DIR}/${test_category}" > "${test_log}" 2>&1
        exit_code=$?
    fi
    
    local test_end=$(date +%s)
    local test_duration=$((test_end - test_start))
    
    # Record result
    if [ $exit_code -eq 0 ]; then
        TESTS_PASSED=$((TESTS_PASSED + 1))
        print_status 0 "${test_name} (${test_duration}s)"
        echo "PASS" > "${TEST_OUTPUT_DIR}/${test_category}/${test_name}.result"
    else
        TESTS_FAILED=$((TESTS_FAILED + 1))
        print_status 1 "${test_name} (${test_duration}s)"
        echo "FAIL: Exit code ${exit_code}" > "${TEST_OUTPUT_DIR}/${test_category}/${test_name}.result"
        
        # Show last few lines of failed test in non-verbose mode
        if [ "$VERBOSE" = false ]; then
            echo "Last 10 lines of output:"
            tail -10 "${test_log}"
        fi
        
        if [ "$STOP_ON_FAIL" = true ]; then
            print_warning "Stopping due to test failure (--stop-on-fail)"
            exit 1
        fi
    fi
}

# Discover test scripts
discover_tests() {
    local category="$1"
    local test_dir="tests/${category}"
    
    if [ ! -d "${test_dir}" ]; then
        return
    fi
    
    find "${test_dir}" -name "test_*.sh" -type f | sort
}

# Run tests in a category
run_test_category() {
    local category="$1"
    
    print_header "Running ${category} tests"
    
    local test_scripts=$(discover_tests "${category}")
    
    if [ -z "${test_scripts}" ]; then
        print_warning "No tests found in ${category}"
        return
    fi
    
    for test_script in ${test_scripts}; do
        run_test "${test_script}"
    done
}

# Generate test summary
generate_summary() {
    local summary_file="${TEST_OUTPUT_DIR}/summary.txt"
    
    cat > "${summary_file}" <<EOF
8bandkp-fdm-ai Test Suite Results
==================================

Run Time: $(date)
Category: ${TEST_CATEGORY}

Test Results:
  Total:   ${TESTS_RUN}
  Passed:  ${TESTS_PASSED}
  Failed:  ${TESTS_FAILED}
  Skipped: ${TESTS_SKIPPED}

EOF
    
    if [ ${TESTS_FAILED} -gt 0 ]; then
        echo "Failed Tests:" >> "${summary_file}"
        find "${TEST_OUTPUT_DIR}" -name "*.result" -type f | while read result_file; do
            if grep -q "FAIL" "${result_file}"; then
                local test_name=$(basename "${result_file}" .result)
                echo "  - ${test_name}" >> "${summary_file}"
            fi
        done
        echo "" >> "${summary_file}"
    fi
    
    echo "Detailed logs available in: ${TEST_OUTPUT_DIR}" >> "${summary_file}"
    
    # Display summary
    echo ""
    print_header "Test Summary"
    cat "${summary_file}"
}

# Main execution
main() {
    parse_args "$@"
    
    print_header "8bandkp-fdm-ai Test Suite"
    echo "Category: ${TEST_CATEGORY}"
    echo ""
    
    # Check if executables exist
    if [ ! -f "bandStructure" ] || [ ! -f "gfactorCalculation" ]; then
        print_warning "Executables not found. Building..."
        make all > "${TEST_OUTPUT_DIR}/build.log" 2>&1 || {
            echo "Build failed! Check ${TEST_OUTPUT_DIR}/build.log"
            exit 1
        }
        print_status 0 "Build completed"
    fi
    
    setup_test_env
    
    local start_time=$(date +%s)
    
    # Run tests based on category
    case "${TEST_CATEGORY}" in
        all)
            run_test_category "unit"
            run_test_category "integration"
            run_test_category "validation"
            ;;
        quick)
            # Run a subset of fast tests
            print_info "Quick mode: running subset of tests"
            run_test_category "unit"
            ;;
        unit|integration|validation)
            run_test_category "${TEST_CATEGORY}"
            ;;
        *)
            echo "Invalid category: ${TEST_CATEGORY}"
            exit 1
            ;;
    esac
    
    local end_time=$(date +%s)
    local total_duration=$((end_time - start_time))
    
    generate_summary
    
    echo ""
    echo "Total execution time: ${total_duration}s"
    echo "Results saved to: ${TEST_OUTPUT_DIR}"
    
    # Exit with appropriate code
    if [ ${TESTS_FAILED} -gt 0 ]; then
        exit 1
    else
        exit 0
    fi
}

# Run main function
main "$@"
