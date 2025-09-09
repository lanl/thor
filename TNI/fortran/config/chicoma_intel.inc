# vim: set filetype=make:

# Compilers
FC := mpifort
MPIRUN := srun

# Compiler options
#OPT := -cpp -g -O2 -fopenmp -fdefault-integer-8 -ffree-line-length-none -fexternal-blas
#OPT := -cpp -g -O2 -ffree-line-length-none -fexternal-blas
OPT := -cpp -g -O2 -qopenmp -traceback

# select BLAS and LAPACK libraries in your system
#BLASLIB = -llapack -lblas
#BLASLIB = -L/projects/opt/centos8/x86_64/oneAPI/2023.1.0.46401/mkl/2023.1.0/lib/intel64/ -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
#BLASLIB = -L$(HOME)/local/share/openblas/haswell/lib64 -lopenblas
BLASLIB = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core # works with intel-mpi/2021.9.0 intel-mkl/2021.3.0

# CUDA LIB FLAGS
CULIBS =

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f

