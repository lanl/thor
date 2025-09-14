# vim: set filetype=make:

# Compilers
FC := mpif90
MPIRUN := mpirun

# Compiler options
OPT := -cpp -g -O2 -fopenmp -ffpe-trap=overflow -Wall -Wextra -Wno-uninitialized -Wno-unused-variable -Wno-unused-parameter -fbacktrace -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -ffree-line-length-none -fallow-argument-mismatch

# select BLAS and LAPACK libraries in your system
BLASLIB = -L$(HOME)/local/share/openblas-gnu13-2 -lopenblas

# CUDA LIB FLAGS
CULIBS =

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f

