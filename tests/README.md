# Test Suite Documentation

## Overview

This directory contains a comprehensive automated test suite for the 8bandkp-fdm-ai project. The tests validate band structure calculations and g-factor calculations for both bulk materials and quantum well systems.

## Test Organization

```
tests/
├── run_tests.sh           # Main test runner
├── unit/                  # Unit tests
│   ├── test_input_parsing.sh
│   └── test_output_format.sh
├── integration/           # Integration tests
│   ├── test_bulk_bandstructure.sh
│   ├── test_qw_bandstructure.sh
│   ├── test_bulk_gfactor_numerical.sh
│   ├── test_bulk_gfactor_analytical.sh
│   └── test_qw_gfactor_analytical.sh
└── validation/            # Validation tests
    ├── test_gfactor_consistency.sh
    └── test_physical_values.sh
```

## Running Tests

### Quick Start

```bash
# Run all tests
make test

# Run specific test categories
make test-unit           # Unit tests only
make test-integration    # Integration tests only
make test-validation     # Validation tests only
make test-quick          # Quick smoke tests
```

### Advanced Usage

```bash
# Run tests with verbose output
./tests/run_tests.sh all --verbose

# Stop on first failure
./tests/run_tests.sh all --stop-on-fail

# Run specific test category
./tests/run_tests.sh integration

# Get help
./tests/run_tests.sh --help
```

### Run Individual Tests

```bash
# Run a single test script directly
bash tests/unit/test_input_parsing.sh
bash tests/integration/test_bulk_bandstructure.sh
bash tests/validation/test_physical_values.sh
```

## Test Coverage

### Unit Tests (2 tests)

**test_input_parsing.sh**
- Validates input file parsing for bulk and QW examples
- Tests both bandStructure and gfactorCalculation executables
- Verifies graceful handling of invalid inputs

**test_output_format.sh**
- Validates output file formats (eigenvalues.dat, gfactor.dat)
- Checks for NaN and Inf values
- Verifies column consistency and data integrity

### Integration Tests (5 tests)

**test_bulk_bandstructure.sh**
- Tests 5 bulk materials: GaAs, InAs, GaSb, InSb, InAs60Sb40
- Validates band gap values against literature
- Checks valence/conduction band positions

**test_qw_bandstructure.sh**
- Tests 5 QW systems:
  - GaAs/Al₀.₃Ga₀.₇As
  - InAs/AlSb
  - InAs/GaSb/AlSb (type-II)
  - In₀.₂Ga₀.₈As/GaAs (strained)
  - InGaAs/InP
- Validates quantized energy levels
- Checks for confinement effects

**test_bulk_gfactor_numerical.sh**
- Tests numerical g-factor for 4 materials
- Validates against expected ranges
- Checks isotropy for cubic bulk

**test_bulk_gfactor_analytical.sh**
- Tests analytical g-factor for CB and VB
- Validates 5 materials
- Checks for reasonable values

**test_qw_gfactor_analytical.sh**
- Tests analytical g-factor for 6 QW systems (CB and VB)
- Validates anisotropy (g_z ≠ g_x, g_y)
- Checks confinement-induced effects

### Validation Tests (2 tests)

**test_gfactor_consistency.sh**
- Compares analytical vs numerical methods for bulk
- Tests 4 materials: GaAs, InAs, GaSb, InSb
- Validates agreement within tolerance

**test_physical_values.sh**
- Compares calculated values against literature
- Validates band gaps for bulk materials
- Validates g-factors against experimental values
- Checks QW anisotropy

## Expected Results

### Band Gaps (eV)
- GaAs: 1.30 - 1.60
- InAs: 0.25 - 0.45
- GaSb: 0.60 - 0.85
- InSb: 0.10 - 0.25
- InAs60Sb40: 0.05 - 0.20

### G-Factors (bulk CB, absolute values)
- GaAs: 0.25 - 0.40
- InAs: 14.0 - 15.5
- GaSb: 8.0 - 9.5
- InSb: 45.0 - 55.0

### QW G-Factors
- Anisotropy: |g_z| > |g_x|, |g_y| for strong confinement
- Example (GaAs/AlGaAs): g_z ~ -140, g_|| ~ -45

## Test Output

All test outputs are saved in timestamped directories:
```
outputs/test-YYYYMMDD-HHMMSS/
├── summary.txt          # Test results summary
├── unit/                # Unit test logs and results
├── integration/         # Integration test outputs
│   ├── bulk/           # Bulk bandstructure results
│   ├── qw/             # QW bandstructure results
│   ├── gfactor_num/    # Numerical g-factor results
│   ├── gfactor_ana/    # Analytical g-factor results
│   └── gfactor_qw/     # QW g-factor results
└── validation/          # Validation test outputs
```

## Adding New Tests

### Test Script Template

```bash
#!/bin/bash
# Test: Description
# Purpose: What this test validates

set -euo pipefail

TEST_OUTPUT_DIR="${1:-outputs/test-category}"
mkdir -p "${TEST_OUTPUT_DIR}"

# Test implementation here
# ...

# Exit with appropriate code
if [ $TESTS_FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
```

### Naming Conventions

- Test files: `test_*.sh`
- Test names: Descriptive and lowercase with underscores
- Output files: Save to `${TEST_OUTPUT_DIR}` subdirectories

### Best Practices

1. **Always clean up**: Remove temporary files (eigenvalues.dat, etc.)
2. **Use proper exit codes**: 0 for success, 1 for failure
3. **Provide informative output**: Echo what's being tested
4. **Validate thoroughly**: Check for NaN, Inf, and range errors
5. **Save artifacts**: Copy important outputs to TEST_OUTPUT_DIR

## Continuous Integration

The test suite is designed for easy CI/CD integration:

```bash
# In CI pipeline
make clean_all
make all
make test  # Exit code 0 = all pass, 1 = failures
```

## Troubleshooting

### Tests Fail Due to Missing Executables
```bash
# Solution: Build first
make clean_all
make all
make test
```

### Tests Leave Temporary Files
```bash
# Solution: Tests should clean up, but manual cleanup:
rm -f eigenvalues.dat gfactor.dat parts.dat eigenfunctions_k_*.dat fort.*
```

### Test Output Directory Full
```bash
# Solution: Clean old test runs
rm -rf outputs/test-*
```

### Individual Test Fails
```bash
# Run test with verbose output
bash tests/path/to/test_name.sh outputs/debug

# Check logs
cat outputs/debug/*.log
```

## Test Maintenance

- **Update expected values**: When material parameters change
- **Add new materials**: Follow existing test patterns
- **Update tolerances**: If physics changes significantly
- **Review periodically**: Ensure tests stay relevant

## Performance

- **Full test suite**: ~3-5 minutes
- **Unit tests only**: ~30 seconds
- **Quick tests**: ~30 seconds
- **Integration tests**: ~2-3 minutes
- **Validation tests**: ~1-2 minutes

## Contact

For questions about the test suite:
- Review this documentation
- Check test script comments
- See main project README.md
