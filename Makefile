# Makefile wrapper for CMake build system
# Preserves familiar make targets: all, run, clean, test

BUILD_DIR := build

all: $(BUILD_DIR)/build.ninja
	cmake --build $(BUILD_DIR)

$(BUILD_DIR)/build.ninja:
	CMPLR_ROOT=/opt/intel/oneapi/compiler/2025.0 cmake -G Ninja -B $(BUILD_DIR) -DMKL_DIR=$(MKLROOT)/lib/cmake/mkl

run: all
	./$(BUILD_DIR)/src/bandStructure

test:
	CMPLR_ROOT=/opt/intel/oneapi/compiler/2025.0 cmake -G Ninja -B $(BUILD_DIR) -DMKL_DIR=$(MKLROOT)/lib/cmake/mkl -DBUILD_TESTING=ON
	cmake --build $(BUILD_DIR)
	cd $(BUILD_DIR) && ctest --output-on-failure

gfactor: all
	./$(BUILD_DIR)/src/gfactorCalculation

clean:
	rm -rf $(BUILD_DIR)

clean_all: clean
	rm -f *.txt fort.* *dat

test-fd-sign:
	@echo "Build and run manually: gfortran tests/fd_sign_compare.f90 -o build/fd_sign_compare && ./build/fd_sign_compare"

lecture-00: all
	python3 scripts/lecture_00_quickstart.py

lecture-01: all
	python3 scripts/lecture_01_bulk.py

lecture-02: all
	python3 scripts/lecture_02_qw.py

lecture-03: all
	python3 scripts/lecture_03_wavefunctions.py

lecture-04: all
	python3 scripts/lecture_04_strain.py

lecture-05: all
	python3 scripts/lecture_05_gfactor.py

lecture-06: all
	python3 scripts/lecture_06_optical.py

lecture-07: all
	python3 scripts/lecture_07_scsp.py

lecture-08: all
	python3 scripts/lecture_08_wire.py

lecture-09: all
	python3 scripts/lecture_09_numerical.py

lecture-10: all
	python3 scripts/lecture_10_qcse.py

lecture-11: all
	python3 scripts/lecture_11_convergence.py

lecture-12: all
	python3 scripts/lecture_12_extending.py

lecture-13: all
	python3 scripts/lecture_13_topological.py

lectures: lecture-00 lecture-01 lecture-02 lecture-03 lecture-04 lecture-05 lecture-06 lecture-07 lecture-08 lecture-09 lecture-10 lecture-11 lecture-12 lecture-13

.PHONY: all run test gfactor clean clean_all test-fd-sign lecture-00 lecture-01 lecture-02 lecture-03 lecture-04 lecture-05 lecture-06 lecture-07 lecture-08 lecture-09 lecture-10 lecture-11 lecture-12 lecture-13 lectures
