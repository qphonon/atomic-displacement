all: fm-forces

#F90=ifort -fPIC -ftrapuv -warn all -CB 
F90=gfortran -fPIC -ffree-line-length-0

#lapack=-L /install-linux-ifort/lapack-3.8.0 -llapack -lrefblas
#lapack=/usr/lib/x86_64-linux-gnu/lapack/liblapack.a /usr/lib/x86_64-linux-gnu/blas/libblas.a 
lapack=-L /install-linux-gfortran/lapack-3.8.0 -llapack -lrefblas

MYPROG=.

types_and_constants.o: $(MYPROG)/types_and_constants.f90
	$(F90) -c types_and_constants.f90

derived_constants.o: types_and_constants.o $(MYPROG)/derived_constants.f90
	$(F90) -c derived_constants.f90
	
sgsym.o: types_and_constants.o $(MYPROG)/sgsym.f90
	$(F90) -c sgsym.f90

commod.o: types_and_constants.o derived_constants.o sgsym.o $(MYPROG)/commod.f90
	$(F90) -c commod.f90

fm-forces: types_and_constants.o derived_constants.o sgsym.o commod.o fm-forces.f90
	$(F90) -o $@ types_and_constants.o sgsym.o commod.o fm-forces.f90 $(lapack)

clean:
	- rm fm-forces *.mod types_and_constants.o derived_constants.o sgsym.o commod.o
