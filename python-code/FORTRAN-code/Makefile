FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
SRC=precision_m.f90 navier_stokes.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
	
NS_solve: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -rf *.mod *.o *.out NS-solve
