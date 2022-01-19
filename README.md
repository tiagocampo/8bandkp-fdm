# 8bandkp-fdm
Fortran implementation of 8band ZB kp-method using finite difference method, working for bulk and quantum-wells

## Disclaim:
 * This software is for my personal usage only. It comes with no warranty at all.
Author: Tiago de Campos

### Description:
 * Solves k.p Hamiltonian by finite differences method;
 * Its based on Chuang, S. L. and Chang, C. S., Semiconductor Science and Technology, 12, 252, 1997.
 * For a description of the derivatives see https://www.wias-berlin.de/people/john/LEHRE/NUM_PDE_FUB/num_pde_fub_2.pdf
 * Use of renormalization of the Kane interband momentum matrix element, P, suggested by Foreman B A 1997 Physical Review B 56 R12748.
 * g-factor calculation using second order Lowdin partitioning. Winkler, R. Spin-orbit coupling effects in two-dimensional electron and hole systems; Physics and Astronomy online Library 191; Springer, 2003. Tadjine, A; Niquet, Y.-M.; Delerue, C. Universal behavior of electron g-factors in semiconductor nanostructures. Physical Review B 2017, 95, 235437.
 * Addition of electric field, but no self-consistency calculation.

### Use:
 * makefile provided
 * example input.cfg provided
 * some materials have built-in parameters, check parameters.f90 to see availability. You can add more if so pleased.
 * Material parameters can be found at Vurgaftman, I. Meyer, J. R. and Ram-Mohan, L. R., Journal of Applied Physics, 11, 5815, 2001.


# Citation

 I used this code to compute results of 10.1021/acsaelm.0c00269    

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



