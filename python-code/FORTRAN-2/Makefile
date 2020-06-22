FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
SRC=precision_m.f90 navier_stokes_module.f90 navier_stokes_solver2.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

NS_solve: $(OBJ)
	@rm -fr FORTRAN-data/FORTRAN-32-data/*.dat
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	@rm -f *.mod *.o NS_solve
