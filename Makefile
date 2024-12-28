FC    := gfortran

EXENAME		:= bandStructure
EXENAME_gfactor := gfactorCalculation
HOMEDIR		:= .
BUILD_DIR   := build
SRC_D		:= $(HOMEDIR)/src
vpath  %.f90 $(SRC_D)
vpath  %.mod $(BUILD_DIR)
vpath  %.o   $(BUILD_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

DEBUG_FLAG := -g -fbacktrace -fbounds-check -fcheck=all -fopenmp #-ffpe-trap=zero,overflow,underflow,denormal #-Wall #-Wextra -pedantic -Wsurprising

FPP := -x f95-cpp-input

FFLAGS=-ffree-line-length-none -mtune=native -march=native\
          -fno-second-underscore -ffree-form -fimplicit-none \
          -m64 -O2 -flto -funroll-loops

# Using Intel MKL (sequential version) for sparse matrix operations
LDFLAGS=-I/usr/include -L/usr/lib64 -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include \
      -I$(MKLROOT)/include/fftw -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential \
      -lmkl_core -lpthread -lm -ldl -lfftw3

all :  exe gfactor
.PHONY: all

gfactor: gfactorCalculation
.PHONY: gfactor

clean_obs:
	rm -rf $(BUILD_DIR)

clean_exe: clean
	rm -f $(EXENAME)*
	rm -f $(EXENAME_gfactor)*

clean_files:
	rm -f *.txt
	rm -f fort.*
	rm -f *dat

clean: clean_obs clean_files

clean_all: clean clean_exe

run: exe
	./$(EXENAME)

OBJS :=\
defs.o\
mkl_spblas.o\
mkl_sparse_handle.o\
parameters.o\
utils.o\
finitedifferences.o\
outputFunctions.o\
hamiltonianConstructor.o

OBJS_gfactor :=\
main_gfactor.o\
gfactor_functions.o

# Special rule for mkl_spblas.o
$(BUILD_DIR)/mkl_spblas.o: mkl_spblas.f90 | $(BUILD_DIR)
	$(FC) -w -fcray-pointer -cpp -J$(BUILD_DIR) -c $< -o $@

# Pattern rule for other object files
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -J$(BUILD_DIR) -c $< -o $@

gfactorCalculation: $(addprefix $(BUILD_DIR)/,$(OBJS)) $(addprefix $(BUILD_DIR)/,$(OBJS_gfactor)) $(BUILD_DIR)/main_gfactor.o
	$(FC) $(FFLAGS) -o $(EXENAME_gfactor) $^ $(LDFLAGS) $(CPPFLAGS) $(DEBUG_FLAG)

exe: $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/main.o
	$(FC) $(FFLAGS) -o $(EXENAME) $^ $(LDFLAGS) $(CPPFLAGS) $(DEBUG_FLAG)

# Dependencies
$(BUILD_DIR)/defs.o: defs.f90
$(BUILD_DIR)/mkl_sparse_handle.o: mkl_sparse_handle.f90 $(BUILD_DIR)/mkl_spblas.o
$(BUILD_DIR)/parameters.o: parameters.f90 $(BUILD_DIR)/defs.o
$(BUILD_DIR)/utils.o: utils.f90 $(BUILD_DIR)/defs.o $(BUILD_DIR)/mkl_spblas.o
$(BUILD_DIR)/finitedifferences.o: finitedifferences.f90 $(BUILD_DIR)/defs.o $(BUILD_DIR)/utils.o
$(BUILD_DIR)/outputFunctions.o: outputFunctions.f90 $(BUILD_DIR)/defs.o
$(BUILD_DIR)/hamiltonianConstructor.o: hamiltonianConstructor.f90 $(BUILD_DIR)/defs.o $(BUILD_DIR)/finitedifferences.o $(BUILD_DIR)/mkl_spblas.o $(BUILD_DIR)/utils.o
$(BUILD_DIR)/main.o: main.f90 $(addprefix $(BUILD_DIR)/,$(OBJS))
$(BUILD_DIR)/gfactor_functions.o: gfactor_functions.f90 $(BUILD_DIR)/defs.o $(BUILD_DIR)/hamiltonianConstructor.o $(BUILD_DIR)/mkl_spblas.o
$(BUILD_DIR)/main_gfactor.o: main_gfactor.f90 $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/gfactor_functions.o

build/utils.o: src/utils.f90
	$(FC) -x f95-cpp-input -ffree-line-length-none -mtune=native -march=native -fno-second-underscore -ffree-form -fimplicit-none -m64 -O2 -flto -funroll-loops -g -fbacktrace -fbounds-check -fcheck=all -fopenmp -Jbuild -c ./src/utils.f90 -o build/utils.o
