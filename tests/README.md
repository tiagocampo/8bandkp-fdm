# Tests

Test suite for the 8bandkp-fdm project using pFUnit (unit tests) and shell/Python scripts (regression tests).

## Structure

```
tests/
  unit/              pFUnit .pf test files
    test_defs.pf            Tests for definitions module (kronij, constants)
    test_finitedifferences.pf  Tests for FD stencils, Toeplitz, buildFD2ndDerivMatrix
    test_utils.pf            Tests for Simpson integration, COO insertion/finalization
    test_parameters.pf       Tests for material parameter database
    test_hamiltonian.pf      Tests for bulk Hamiltonian (Hermiticity, Gamma point)
  integration/       Shell scripts for regression tests
    test_bulk_bandstructure.sh   Run bulk config, compare eigenvalues
    test_qw_bandstructure.sh     Run QW config, compare eigenvalues
    test_gfactor.sh              Run gfactorCalculation, compare g-values
  regression/
    configs/           Input config files for regression tests
    data/              Reference (golden) output data
    compare_output.py  Python script for numerical comparison (tolerance-aware)
```

## Running Tests

### Configure with testing enabled

```bash
# pFUnit must be installed first
cmake -G Ninja -B build \
    -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON \
    -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-<ver>/cmake
cmake --build build
```

### Run all tests

```bash
ctest --test-dir build              # all tests
ctest --test-dir build -L unit      # pFUnit unit tests only
ctest --test-dir build -L regression  # regression/golden-output tests only
ctest --test-dir build -V           # verbose output
```

### Or via Make wrapper

```bash
make test
```

## Test Categories

| Label | Description | Framework |
|---|---|---|
| `unit` | Individual module tests (FD stencils, Simpson, parameters, Hamiltonian) | pFUnit |
| `regression` | Full-executable tests comparing against reference data | Shell + Python |

## Adding New Tests

### Unit test (pFUnit)

1. Create `tests/unit/test_<module>.pf` with `@test` subroutines
2. Add `add_pfunit_ctest(test_<module> ...)` to `tests/CMakeLists.txt`
3. Build and run: `cmake --build build && ctest --test-dir build -L unit`

### Regression test

1. Create a config file in `tests/regression/configs/`
2. Run the executable, save reference output to `tests/regression/data/<name>/`
3. Add a shell script in `tests/integration/` or an `add_test()` block in `tests/CMakeLists.txt`

## Prerequisites

- **pFUnit** >= 4.9: Build from source at https://github.com/Goddard-Fortran-Ecosystem/pFUnit
  ```bash
  git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit
  cd pFUnit && mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local/pfunit
  make install
  ```
- **Python 3**: For `compare_output.py` regression comparison
