# vim: set filetype=make:

# Compilers
FC := mpifort
MPIRUN := mpirun

# Optimization flags
#OPT := -cpp -g -O2 -fopenmp
OPT := -cpp -g -O2 -fopenmp -ffpe-trap=overflow -Wall -Wextra -Wno-uninitialized -Wno-unused-variable -Wno-unused-parameter -Wno-argument-mismatch -fbacktrace -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -ffree-line-length-none -fexternal-blas

# Multiple-precision calculations enable MPFR library below
#MPFLIB  = -lmpfr
#MPFLIB = -I$(HOME)/local/openmpi/3.1.1/include/mpif.h

# select BLAS and LAPACK libraries in your system
#BLASLIB = -I$(HOME)/local/include -L$(HOME)/local/lib -llapack -lrefblas
#BLASLIB = -I$(HOME)/local/include -L$(HOME)/local/lib -llapack -lblas
BLASLIB = -I$(HOME)/local/share/openblas/include -L$(HOME)/local/share/openblas/lib -lopenblas

# HDF5 library
HDF5_DIR = $(HOME)/local/share/hdf5/1.14.4
HDF5_INC = -I$(HDF5_DIR)/include
HDF5_LIB = -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -Wl,-rpath -Wl,$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran

# CUDA LIB FLAGS
#CULIBS = -stdpar=gpu,multicore
CULIBS =

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f

