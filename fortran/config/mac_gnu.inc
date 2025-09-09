# vim: set filetype=make:

# Compilers
FC := mpifort
MPIRUN := mpirun

# Optimization flags
#OPT := -cpp -g -O2 -fopenmp
OPT := -cpp -g -O2 -fopenmp -fallow-argument-mismatch -ffpe-trap=overflow,underflow,invalid -Wall -Wextra -Wno-uninitialized -Wno-unused-variable -Wno-unused-parameter -fbacktrace -fimplicit-none -fdefault-double-8 -ffree-line-length-none -fexternal-blas

# Multiple-precision calculations enable MPFR library below
#MPFLIB  = -lmpfr
#MPFLIB = -I$(HOME)/local/openmpi/3.1.1/include/mpif.h

# select BLAS and LAPACK libraries in your system
# Define environment variable OPENBLAS_DIR in your rc file (.zshrc or .bash_profile)
BLASLIB = -I$(OPENBLAS_DIR)/include -L$(OPENBLAS_DIR)/lib -llapack -lblas

# CUDA LIB FLAGS
#CULIBS = -stdpar=gpu,multicore
CULIBS =

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f

