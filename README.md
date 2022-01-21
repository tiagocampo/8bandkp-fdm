# 8bandkp-fdm
Fortran implementation of 8band ZB kp-method using finite difference method, working for bulk and quantum-wells

## Disclaim:
 * This software is for my personal usage only. It comes with no warranty at all.
Author: Tiago de Campos

## Description:
 * Solves k.p Hamiltonian by finite differences method;
 * Its based on Chuang, S. L. and Chang, C. S., Semiconductor Science and Technology, 12, 252, 1997.
 * For a description of the derivatives see https://www.wias-berlin.de/people/john/LEHRE/NUM_PDE_FUB/num_pde_fub_2.pdf
 * Use of renormalization of the Kane interband momentum matrix element, P, suggested by Foreman B A 1997 Physical Review B 56 R12748.
 * g-factor calculation using second order Lowdin partitioning. Winkler, R. Spin-orbit coupling effects in two-dimensional electron and hole systems; Physics and Astronomy online Library 191; Springer, 2003. Tadjine, A; Niquet, Y.-M.; Delerue, C. Universal behavior of electron g-factors in semiconductor nanostructures. Physical Review B 2017, 95, 235437.
 * Addition of electric field, but no self-consistency calculation.

## Use:

### General info

 * makefile provided
 * example input.cfg provided
 * some materials have built-in parameters, check parameters.f90 to see availability. You can add more if so pleased.
 * Material parameters can be found at Vurgaftman, I. Meyer, J. R. and Ram-Mohan, L. R., Journal of Applied Physics, 11, 5815, 2001.

### Compilation steps

Steps to compile and run the program

 * modify Makefile to reflect either intel mkl library or just un-optimized lapack and blas. You should choose the correct LDFLAGS.
 * make all

This will generate two executable files: bandStructure and gfactorCalculation.

To compute the band structure just run ./bandStructure

To compute the gfactors just run ./gfactorCalculation

### Input file structure

 * waveVector: reciprocal space direction kx or ky or kz. For confined system should be kx or ky since kz is the confined direction.
 * waveVectorMax: percentage of Brillouin Zone to compute. Should not be too far due to k.p being an approximation to the Brillouin zone center.
 * waveVectorStep: number of points between k=0 and k=kmax
 * confinement: 0 for bulk 1 for quantum well
 * FDStep: number of discretization points.
 * numLayers: number of layers that your quantum well is formed. Centered at 0.
 * material
	* material1: initial and final position of material. Usually material 1 is the host and the rest of materials make the well inside. Band offset. See input.cfg for example
	* material2: inital and final position and band offset
	* etc
 * numcb: number of conduction band. Should be 2 for bulk and not exceed 2 x fdstep for confined system
 * numvb: number of valence band. Should be 6 for bulk and not exceed 6 x fdstep for confined system.
 * ExternalField: 0 or 1 and EF. Sets filed on aor off and set type. EF means electric field and it is the only one implemented
 * EFParams: strenght of the field.

### Pre-requisits

 * fortran compiler, blas, lapack and fftw3
 * On ubuntu: sudo apt install gfortran gcc g++ liblapack-dev libblas-dev libfftw3-dev 

The above will install un-optimized lapack and blas. For the intel version, search for intel mkl to see how to install it.


# Citation

 I used this code to compute results of 10.1021/acsaelm.0c00269 and part of 10.1088/1361-648X/ab38a1.

 Please, add me to you citation list!

## If you use this software and wish to contribute to my citation count, here is a list of articles you can cite.

 * 10.1021/acsaelm.0c00269 
 * 10.1063/1.5096970 
 * 10.1088/1361-648X/ab38a1
 * 10.1103/PhysRevB.97.245402
 * 10.1088/0268-1242/31/10/105002
 * 10.1103/PhysRevB.93.235204
 * 10.1063/1.4901209

# LICENSE
GNU General Public License v3.0
just please cite 10.1021/acsaelm.0c00269 or one of the other listed articles if it suits more your work.



