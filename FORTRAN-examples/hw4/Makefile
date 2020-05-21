FFTW_DIR = /opt/fftw
F90 = gfortran


FFLAGS = -c -fbacktrace -cpp -DLIT_SKIP_HOOKS -Wextra -W \
	      -Wno-unused-function -fbounds-check -fopenmp -ffpe-trap=invalid,zero,overflow,denormal \
	      -fdump-core -Wunderflow
FFTW_INC = -I$(FFTW_DIR)/include/
FFTW_LIB = -L$(FFTW_DIR)/lib/ -lfftw3 -lm
# $(F90) $(FFLAGS) fftw_test.f90

objects = fftw_test.o

fftwmake: $(objects)
		$(F90) -o main $(objects) $(FFTW_LIB)
fftw_test.o: fftw_test.f90
		$(F90) $(FFTW_INC) $(FFLAGS) fftw_test.f90
clean:
		rm $(objects) main
