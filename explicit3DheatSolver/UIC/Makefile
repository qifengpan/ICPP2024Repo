GFORTRAN    = gfortran -Wall -O3 -g -march=native -mprefer-vector-width=512 -funroll-loops
NFORT       = nfort -Wall -O3 -floop-unroll-complete-nest=3 -ftrace
CC          = gcc -Wall -O3
MPI_FORTRAN = mpif90 -Wall -O3 -fallow-argument-mismatch -g
MPI_CC      = mpicc -Wall -O3
LD = -lm

.SUFFIXES : .o .f90

all: explicitSeq explicitPar

explicitSeq : explicitSeq.o explUtilSeq.o
	$(GFORTRAN) -o $@ explicitSeq.o explUtilSeq.o 

explicitSeq.o : explicitSeq.f90 
	$(GFORTRAN) -c  $(*F).f90

explUtilSeq.o : explUtilSeq.f90
	$(GFORTRAN) -c $(*F).f90	

explicitPar : explicitPar.o explUtilPar.o updateBound.o readParam.o
	$(MPI_FORTRAN) -o $@ explicitPar.o explUtilPar.o updateBound.o readParam.o 

.f90.o :
	$(MPI_FORTRAN) -c -g  $(*F).f90

clean : 
	/bin/rm -f *.o explicitSeq explicitPar
