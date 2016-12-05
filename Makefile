
CC=gcc
FC=gfortran-4.8
NVCC=/usr/local/cuda-5.5/bin/nvcc
#FFLAGS= -cpp -Wall -pedantic -O3
#FFLAGS= -cpp -march=native -mtune=native -O3 
FFLAGS= -cpp -O3 
CFLAGS= -O3 -std=c99
NVFLAGS= -O3 -arch=sm_20

OPENMP= -fopenmp
INCS= -I. 
NVINCS= -I /usr/local/cuda/include
OPENMPLIBS= 
LIBS= -lm 
NVLIBS= -L/usr/local/cuda/lib64 -lcuda -lcudart 

F_OBJECTS=  ctable.o srchnd.o krige.o 
C_OBJECTS=  srchndvectorized.o 
CU_OBJECTS= 

GSLIB= gslib/gslib.a 

all: seq-for seq-lev-for par-lev-for

seq-for: 
	cd gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd ..
	$(FC) -c $(FFLAGS) sisim.for
	$(FC) $(FFLAGS) $(INCS) sisim.o -o sisimFortranSeq.exe $(GSLIB) $(LIBS)

par-lev-for:
	cd gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd ..
	$(FC) -c $(FFLAGS) $(OPENMP) ctable.for
	$(FC) -c $(FFLAGS) $(OPENMP) srchnd.for
	$(FC) -c $(FFLAGS) $(OPENMP) krige.for
	$(FC) -DFORTRAN -c $(FFLAGS) $(OPENMP) sisimLevels.for
	gcc-4.8 -c -O3 -msse4.1 srchndvectorized.c
	$(FC) -DFORTRAN -DLEVELS $(FFLAGS) $(OPENMP) $(INCS) sisimLevels.o $(C_OBJECTS) $(F_OBJECTS) -o sisimFortranLevelsPar.exe $(GSLIB) $(LIBS) $(OPENMPLIBS)

seq-lev-for:
	cd gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd ..
	$(FC) -c $(FFLAGS) ctable.for
	$(FC) -c $(FFLAGS) srchnd.for
	$(FC) -c $(FFLAGS) krige.for
	$(FC) -DFORTRAN -c $(FFLAGS) sisimLevels.for
	gcc-4.8 -c -O3 -msse4.1 srchndvectorized.c
	$(FC) -DFORTRAN -DLEVELS $(FFLAGS) $(INCS) sisimLevels.o $(C_OBJECTS) $(F_OBJECTS) -o sisimFortranLevelsSeq.exe $(GSLIB) $(LIBS)

clean:
	rm *.exe *.o *.mod gslib/*.o gslib/*.a gslib/*.mod
