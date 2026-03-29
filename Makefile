# Makefile wrapper for CMake build system
# Preserves familiar make targets: all, run, clean, test

BUILD_DIR := build

all: $(BUILD_DIR)/build.ninja
	cmake --build $(BUILD_DIR)

$(BUILD_DIR)/build.ninja:
	cmake -G Ninja -B $(BUILD_DIR)

run: all
	./$(BUILD_DIR)/bandStructure

test: all
	cmake -B $(BUILD_DIR) -DBUILD_TESTING=ON
	cmake --build $(BUILD_DIR)
	cd $(BUILD_DIR) && ctest --output-on-failure

gfactor: all
	./$(BUILD_DIR)/gfactorCalculation

clean:
	rm -rf $(BUILD_DIR)

clean_all: clean
	rm -f bandStructure gfactorCalculation
	rm -f *.txt fort.* *dat

.PHONY: all run test gfactor clean clean_all
