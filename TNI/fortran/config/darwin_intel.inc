# vim: set filetype=make:

# Compilers
FC := mpifort
MPIRUN := mpirun

# Compiler options
OPT := -cpp -g -O2 -fopenmp -traceback

# select BLAS and LAPACK libraries in your system
#BLASLIB = -llapack -lblas
BLASLIB = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core # works with intel-mpi/2021.9.0 intel-mkl/2021.3.0
#BLASLIB = -L$(HOME)/local/share/openblas/haswell/lib64 -lopenblas
#BLASLIB = -L$(HOME)/local/share/openblas/milan/lib64 -lopenblas
#BLASLIB = -L$(HOME)/local/share/openblas/skylake_platinum -lopenblas
#BLASLIB = -L$(HOME)/local/share/openblas/redstone/lib64 -lopenblas

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f

