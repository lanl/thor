# vim: set filetype=make:

# Compilers
FC := mpifort
MPIRUN := mpirun

# Optimization flags
#OPT := -cpp -g -O2 -fopenmp
OPT := -cpp -g -O2 -fopenmp -ffpe-trap=overflow,underflow,invalid -Wall -Wextra -Wno-uninitialized -Wno-unused-variable -Wno-unused-parameter -fbacktrace -fimplicit-none -fdefault-double-8 -ffree-line-length-none

# Multiple-precision calculations enable MPFR library below
#MPFLIB  = -lmpfr
#MPFLIB = -I$(HOME)/local/openmpi/3.1.1/include/mpif.h

# select BLAS and LAPACK libraries in your system
#BLASLIB = -I$(HOME)/local/include -L$(HOME)/local/lib -llapack -lrefblas
#BLASLIB = -I$(HOME)/local/include -L$(HOME)/local/lib -llapack -lblas
BLASLIB = -L$(HOME)/local/share/openblas/arm -lopenblas

# CUDA LIB FLAGS
#CULIBS = -stdpar=gpu,multicore
CULIBS =

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f

