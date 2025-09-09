# vim: set filetype=make:

# Compilers
FC := mpif90
MPIRUN := mpirun

# Compiler options
#OPT := -cpp -g -O2 -fopenmp
#OPT := -cpp -g -O2 -fopenmp -ffree-line-length-none -fexternal-blas -fallow-argument-mismatch
OPT := -cpp -g -O2 -fopenmp -ffpe-trap=overflow -Wall -Wextra -Wno-uninitialized -Wno-unused-variable -Wno-unused-parameter -fbacktrace -fimplicit-none -fdefault-integer-8 -fdefault-double-8 -ffree-line-length-none -fexternal-blas -fallow-argument-mismatch

# select BLAS and LAPACK libraries in your system
#BLASLIB = -llapack -lblas
#BLASLIB = -L/projects/opt/centos8/x86_64/oneAPI/2023.1.0.46401/mkl/2023.1.0/lib/intel64/ -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
#BLASLIB = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core # works with intel-mpi/2021.9.0 intel-mkl/2021.3.0
#BLASLIB = -L$(HOME)/local/share/openblas/haswell/lib64 -lopenblas
BLASLIB = -L/vast/home/korobkin/local/share/openblas-gnu13-2 -lopenblas
#BLASLIB = -L$(HOME)/local/share/openblas/milan/lib64 -lopenblas
#BLASLIB = -L$(HOME)/local/share/openblas/skylake_platinum -lopenblas
#BLASLIB = -L$(HOME)/local/share/openblas/redstone/lib64 -lopenblas

# CUDA LIB FLAGS
CULIBS =

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f

