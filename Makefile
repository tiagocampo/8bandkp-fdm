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

LDFLAGS=-I/usr/include -L/usr/lib64 -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include \
      -I$(MKLROOT)/include/fftw -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread \
      -lmkl_core -liomp5 -lpthread -lm -ldl -lfftw3 #-lfftw3f -lfftw3l
# #
# LDFLAGS=-I/usr/include -L/usr/lib64 -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include \
#         -I$(MKLROOT)/include/fftw -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential \
# 				-lmkl_core -lpthread -lm -ldl -lfftw3 #-lfftw3f #-lfftw3l -lfftw3q -lquadmath

all :  exe
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
gfactor_functions.o

mkl_spblas.o: mkl_spblas.f90
	$(FC) -w -fcray-pointer -cpp -c $< -o $@

gfactorCalculation: $(OBJS) $(OBJS_gfactor) main_gfactor.o
	$(FC) $(FFLAGS) -o $(EXENAME_gfactor) $^ $(LDFLAGS) $(CPPFLAGS) $(DEBUG_FLAG)

gfactor_functions.o: gfactor_functions.f90 defs.o hamiltonianConstructor.o mkl_spblas.o
	$(FC) $(FPP) $(FFLAGS) $(DEBUG_FLAG) -c $< -o $@

main_gfactor.o: main_gfactor.f90 $(OBJS) gfactor_functions.o
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
