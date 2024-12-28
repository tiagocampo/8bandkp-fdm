# Tests

This directory contains the test suite for the 8bandkp-fdm project.

## Structure

* `unit/`: Unit tests for individual components
  - Tests for individual modules and functions
  - Verifies basic functionality
  - Checks edge cases and error handling

* `integration/`: Integration tests
  - Tests interaction between modules
  - Verifies combined functionality
  - Checks data flow between components

* `regression/`: Regression tests
  - Tests against known good results
  - Verifies physical correctness
  - Ensures consistency across changes

## Running Tests

Tests can be enabled during CMake configuration:
```bash
cmake -DBUILD_TESTING=ON ..
make
ctest
```

## Test Data

* Reference data for regression tests is stored in `regression/data/`
* Test configurations are stored in `unit/configs/` and `integration/configs/`
* Results are compared against known good outputs

## Adding Tests

1. Create a new test file in the appropriate directory
2. Add test to CMakeLists.txt in the test directory
3. Add any necessary test data
4. Document test purpose and requirements
5. Verify test coverage 