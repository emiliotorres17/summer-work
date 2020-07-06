FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
SRC=precision_m.f90 SOR.f90
OBJ=${SRC:.f90=.o}


%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
SOR_solve: $(OBJ)
	@rm -fr data/*.dat
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	@rm -f *.mod *.o NS_solve
