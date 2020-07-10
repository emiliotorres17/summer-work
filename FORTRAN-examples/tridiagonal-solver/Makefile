FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
SRC=precision_m.f90 tri_solver.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
tri_solve: $(OBJ)
	@rm -fr data/*.dat
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	@rm -f *.mod *.o tri_solve
