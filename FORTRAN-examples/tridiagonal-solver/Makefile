FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
SRC=precision_m.f90 library.f90 laasonen_example.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
btcs_solve: $(OBJ)
	@rm -fr data/*.dat
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	@rm -f *.mod *.o btcs_solve
