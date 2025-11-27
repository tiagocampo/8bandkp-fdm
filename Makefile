FC    := gfortran

EXENAME		:= bandStructure
EXENAME_gfactor := gfactorCalculation
HOMEDIR		:= .
SRC_D		:= $(HOMEDIR)/src
SRC_D		+= $(HOMEDIR)/src
vpath  %.f90 $(SRC_D)
vpath  %.mod $(HOMEDIR)
vpath  %.o   $(HOMEDIR)


DEBUG_FLAG := -g -fbacktrace -fbounds-check -fcheck=all -fopenmp #-ffpe-trap=zero,overflow,underflow,denormal #-Wall #-Wextra -pedantic -Wsurprising

FPP := -x f95-cpp-input

# FFLAGS=-ffree-line-length-none -mtune=native -march=skylake\
#           -fno-second-underscore -ffree-form -fimplicit-none \
#           -O2 -m64 -ffast-math -fno-unsafe-math-optimizations \
# 	        -frounding-math -fsignaling-nans

FFLAGS=-ffree-line-length-none -mtune=native -march=native\
          -fno-second-underscore -ffree-form -fimplicit-none \
          -m64 -O2 -flto -funroll-loops

# if you have installed mkl use this ldflgas
LDFLAGS=-I/usr/include -L/usr/lib64 -L/opt/intel/oneapi/mkl/2025.0/lib/intel64 -I/opt/intel/oneapi/mkl/2025.0/include \
      -I/opt/intel/oneapi/mkl/2025.0/include/fftw -L/opt/intel/oneapi/compiler/2025.0/lib -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread \
      -lmkl_core -liomp5 -lpthread -lm -ldl -lfftw3 #-lfftw3f -lfftw3l #-larpack

# if you dont have mkl, but have normal lapack and blass use this
#LDFLAGS=-I/usr/include -L/usr/lib64 -lm -lfftw3 -llapack -lblas

# #
# LDFLAGS=-I/usr/include -L/usr/lib64 -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include \
#         -I$(MKLROOT)/include/fftw -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential \
# 				-lmkl_core -lpthread -lm -ldl -lfftw3 #-lfftw3f #-lfftw3l -lfftw3q -lquadmath

all :  exe gfactor
.PHONY: all

gfactor: gfactorCalculation
.PHONY: gfactor

clean_obs:
	rm -f $(HOMEDIR)/*.o
	rm -f $(HOMEDIR)/*.mod

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
parameters.o\
hamiltonianConstructor.o\
finitedifferences.o\
outputFunctions.o\
mkl_spblas.o\
utils.o

OBJS_gfactor :=\
main_gfactor.o\
gfactor_functions.o\
perturbation_errors.o

mkl_spblas.o: mkl_spblas.f90 defs.o
	$(FC) -w -fcray-pointer -cpp -c $< -o $@

gfactorCalculation: $(OBJS) $(OBJS_gfactor) main_gfactor.o
	$(FC) $(FFLAGS) -o $(EXENAME_gfactor) $^ $(LDFLAGS) $(CPPFLAGS) $(DEBUG_FLAG)

gfactor_functions.o: gfactor_functions.f90 defs.o hamiltonianConstructor.o mkl_spblas.o perturbation_errors.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

main_gfactor.o: main_gfactor.f90 $(OBJS) gfactor_functions.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

# Perturbation error handling dependency
perturbation_errors.o: perturbation_errors.f90 defs.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

finitedifferences.o: finitedifferences.f90 defs.o utils.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

outputFunctions.o: outputFunctions.f90 defs.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

defs.o: defs.f90
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

utils.o: utils.f90 defs.o mkl_spblas.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

parameters.o: parameters.f90 defs.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

hamiltonianConstructor.o: hamiltonianConstructor.f90 defs.o finitedifferences.o mkl_spblas.o utils.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

main.o: main.f90 $(OBJS)
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

exe: $(OBJS) main.o
	$(FC) $(FFLAGS) -o $(EXENAME) $^ $(LDFLAGS) $(CPPFLAGS) $(DEBUG_FLAG)

# Test targets
test: all
	@echo "Running all tests..."
	@bash tests/run_tests.sh all

test-unit: all
	@echo "Running unit tests..."
	@bash tests/run_tests.sh unit

test-integration: all
	@echo "Running integration tests..."
	@bash tests/run_tests.sh integration

test-validation: all
	@echo "Running validation tests..."
	@bash tests/run_tests.sh validation

test-quick: all
	@echo "Running quick tests..."
	@bash tests/run_tests.sh quick

.PHONY: test test-unit test-integration test-validation test-quick

