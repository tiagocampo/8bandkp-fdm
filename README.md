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

 I used this code to compute results of

    MAYER, WILLIAM ; SCHIELA, WILLIAM F. ; YUAN, JOSEPH ; HATEFIPOUR, MEHDI ; SARNEY, WENDY L. ; SVENSSON, STEFAN P. ; LEFF, ASHER C. ; CAMPOS, TIAGO ; WICKRAMASINGHE, KAUSHINI S. ; DARTIAILH, MATTHIEU C. ; 'UTI', IGOR ; SHABANI, JAVAD . Superconducting Proximity Effect in InAsSb Surface Quantum Wells with In Situ Al Contacts . ACS Applied Electronic Materials, v. 2, p. 2351-2356, 2020.

Please, add me to you citation list!

## If you use this software and wish to contribute to my citation count, here is a list of articles you can cite.

    MAYER, WILLIAM ; SCHIELA, WILLIAM F. ; YUAN, JOSEPH ; HATEFIPOUR, MEHDI ; SARNEY, WENDY L. ; SVENSSON, STEFAN P. ; LEFF, ASHER C. ; CAMPOS, TIAGO ; WICKRAMASINGHE, KAUSHINI S. ; DARTIAILH, MATTHIEU C. ; 'UTI', IGOR ; SHABANI, JAVAD . Superconducting Proximity Effect in InAsSb Surface Quantum Wells with In Situ Al Contacts	. ACS Applied Electronic Materials, v. 2, p. 2351-2356, 2020.
    DIRNBERGER, FLORIAN ; KAMMERMEIER, MICHAEL ; KÖNIG, JAN ; FORSCH, MORITZ ; FARIA JUNIOR, PAULO E. ; CAMPOS, TIAGO ; FABIAN, JAROSLAV ; SCHLIEMANN, JOHN ; SCHÜLLER, CHRISTIAN ; KORN, TOBIAS ; WENK, PAUL ; BOUGEARD, DOMINIQUE . Ultralong spin lifetimes in one-dimensional semiconductor nanowires . APPLIED PHYSICS LETTERS, v. 114, p. 202101, 2019.
    CAMPOS, T; TOLOZA SANDOVAL, M A ; DIAGO-CISNEROS, L ; SIPAHI, G M . Electrical tuning of helical edge states in topological multilayers . JOURNAL OF PHYSICS-CONDENSED MATTER, v. 31, p. 495501, 2019.
    CAMPOS, T.; FARIA JUNIOR, PAULO E. ; GMITRA, MARTIN ; SIPAHI, GUILHERME M. ; FABIAN, JAROSLAV . Spin-orbit coupling effects in zinc-blende InSb and wurtzite InAs nanowires: Realistic calculations with multiband k · p method. PHYSICAL REVIEW B, v. 97, p. 245402, 2018.
    BASTOS, C. M. O. ; SABINO, F. ; FARIA JUNIOR, PAULO E. ; CAMPOS, T. ; SILVA, J. L. F. ; SIPAHI, G. M. . Stability and accuracy control of k · p parameters. Semiconductor Science and Technology (Print), v. 31, p. 105002, 2016.
    FARIA JUNIOR, PAULO E. ; CAMPOS, TIAGO ; BASTOS, CARLOS M. O. ; GMITRA, MARTIN ; FABIAN, JAROSLAV ; SIPAHI, GUILHERME M. . Realistic multiband approach from and spin-orbit coupling effects of InAs and InP in wurtzite phase. Physical Review B, v. 93, p. 235204, 2016.
    FARIA JUNIOR, P. E. ; CAMPOS, T. ; SIPAHI, G. M. . Interband polarized absorption in InP polytypic superlattices. JOURNAL OF APPLIED PHYSICS, v. 116, p. 193501, 2014.


