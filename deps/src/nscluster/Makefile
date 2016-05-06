# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>

ERR = $(shell which ifort>/dev/null; echo $$?)
LIB=nscluster
ifeq "$(ERR)" "0"
	FC=ifort
	FFLAGS=-fPIC -O3
else
	FC=gfortran
	FFLAGS=-fPIC -O3 -ffast-math
endif
SUFFIX=so


OBJ=$(patsubst %.f,%.o,$(shell ls *.f))
all: lib$(LIB).$(SUFFIX)
$(OBJ):
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
lib$(LIB).$(SUFFIX): $(OBJ)
	$(FC) -o lib$(LIB).$(SUFFIX) -shared $?
clean:
	rm *.o *.$(SUFFIX)
