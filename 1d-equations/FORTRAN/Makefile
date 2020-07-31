FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
SRC=precision_m.f90 1d-navier-stokes-library.f90 main_NS_1d.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
NS_solve: $(OBJ)
	@rm -fr data/*.dat
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	@rm -f *.mod *.o NS_solve
