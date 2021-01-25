## atomic-displacement
# Efficacious symmetry-adapted atomic displacement method for lattice dynamical studies

For theory see Computer Physics Communications, volume 259 (2021) 107635 (or  https://arxiv.org/abs/2007.06916 ).

To build the binary `fm-forces', modify Makefile accordingly. 
The only external library is LAPACK/BLAS.
We have compiled the codes successfully under gfortran (GNU Fortran 9.3.0, 7.5.0, 6.3.0, 4.4.7) 
and ifort (11.1).

Once the binary `fm-forces' is built, go to example1-Si8 with a 2x2x2 supercell.

Issue "../fm-forces ." (omit the double quotes. Note ".", the current directory, is
actually an argument to fm-forces). fm-forces will read the input files ineq.vasp, 
p.vasp, supercell.vasp, fm-forces.par, etc and it will
attempt to read fm-forces.dat. Under normal circumstances fm-forces
will obviously fail the first time since it is looking for fm-forces.dat but
DFT force calculations have not been carried out. However, even though fm-forces
halts in the middle of execution, all the necessary POSCAR files with finite displacements
will be saved in the directory fmdisp-poscar-dir. Run VASP on all of them and collect forces
under TOTAL-FORCES in OUTCAR to produce fm-forces.dat.
Now, a second run of "../fm-forces ." will produce partial-forces.dat for further processing. 
Note that if you were to run "../fm-forces ." in example1-Si8, fm-forces will
not crash but it will reach the end of a complete execution since a good fm-forces.dat 
has been prepared in the same directory.

A second example2-Grapene is with a 4x4x1 supercell for graphene.

A third example3-hcp-Mg is with a 3x3x2 supercell for Mg (the lattice parameters are 
arbitrarily chosen to test the code only).

The stdout gives a sense of the code and contains many useful information including the 
displacement patterns.

It is hoped that a complete writeup about the input files will be available when the
whole phonon package called QPHONON is distributed.
