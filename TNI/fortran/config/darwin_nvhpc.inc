# vim: set filetype=make:
# Compilers
#FC := nvfortran
#FC := mpif90  -ffree-line-length-none 
FC := mpif90
MPIRUN := mpirun

# Optimization flags
#OPT    = -cpp -g -O2 -fopenmp
#OPT    = -cpp -g -traceback -O2 -fopenmp
#OPT    = -cpp -O1
OPT    = -cpp

# Multiple-precision calculations enable MPFR library below
#MPFLIB  = -lmpfr
#MPFLIB  = -I/vast/home/iboureima/apps/mpfr/4.2.0/include/lmpfr.h

# select BLAS and LAPACK libraries in your system
BLASLIB = -llapack -lblas
#BLASLIB = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core

# CUDA LIB FLAGS
#CULIBS = -stdpar=gpu,multicore
#CULIBS = -stdpar=multicore

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f


