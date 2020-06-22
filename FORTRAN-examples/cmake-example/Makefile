FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
VAR3=planet
SRC=geometry.f90 math.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

maths: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	@rm -f *.mod *.o maths
