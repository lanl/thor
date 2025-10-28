# vim: set filetype=make:
# Compilers
#FC := nvfortran
#FC := mpif90  -ffree-line-length-none 
FC := mpif90
NVCC :=nvcc

# Optimization flags
OPT    = -cpp

# select BLAS and LAPACK libraries in your system
BLASLIB = -llapack -lblas
#BLASLIB = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core

# CUDA LIB FLAGS
MPI_PATH=/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/comm_libs/11.8/hpcx/hpcx-2.14/ompi
UCC_PATH=/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/comm_libs/11.8/hpcx/hpcx-2.14/ucc
CUDA_PATH=/vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/cuda

CUDA_CAL_PATH       = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/math_libs/11.8/targets/x86_64-linux
CUBLASMP_PATH       = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/math_libs/11.8/include
CUSOLVERMP_PATH     = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/math_libs/11.8/include

CUDA_LIB_PATH       = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/cuda/lib64
CUSOLVER_LIB_PATH   = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/math_libs/11.8/lib64
CUBLASMP_LIB_PATH   = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/math_libs/11.8/lib64
CUSOLVERMP_LIB_PATH = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/math_libs/11.8/lib64
CUDAMP_LIB_PATH     = /vast/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/math_libs/11.8/lib64

CUMP_INCLUDES     = -I$(MPI_PATH)/include  -I$(UCC_PATH)/include -I$(CUDA_PATH)/include  -I$(CUDA_CAL_PATH)/include
CUMP_LIBS         = -L$(CUDAMP_LIB_PATH) -lcusolverMp -lcusolver -lcublasmp -lcublas  -L$(CUDA_CAL_PATH)/lib -lcal -L$(CUDA_LIB_PATH) -lcudart -L$(MPI_PATH)/lib -lmpi -L$(UCC_PATH)/lib/ -lucc


#NCCL_LIBS = -L/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/23.3/comm_libs/nccl/lib -lnccl
#NCCL_LIBS = -L/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.5/comm_libs/nccl/lib -lnccl
NCCL_LIBS = -L/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/comm_libs/nccl/lib -lnccl
#NCCL_LIBS = -Wl,-rpath=/projects/darwin-nv/rhel8/x86_64/packages/nvhpc//Linux_x86_64/24.7/comm_libs/nccl/lib/libnccl.so.2

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(NCCL_LIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f


