# vim: set filetype=make:
# Compilers
#FC := nvfortran
#FC := mpif90  -ffree-line-length-none 
FC := mpif90
NVCC :=nvcc
# NVCC :=nvcc -gencode arch=compute_80,code=sm_80

CUFC := mpif90
CUCC := nvcc
CUFC_FLAGS= -cuda -cpp
#CUFC_FLAGS= -cuda -cpp -gpu=cc80  ## adding -gpu=cc80 screws up kernel compilation

# Optimization flags
OPT    = -cpp

# select BLAS and LAPACK libraries in your system
BLASLIB = -llapack -lblas
#BLASLIB = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core

# CUDA LIB FLAGS
NVHPC_ROOT=/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/25.1
MPI_PATH=$(NVHPC_ROOT)/comm_libs/12.6/hpcx/hpcx-2.21/ompi
UCC_PATH=$(NVHPC_ROOT)/comm_libs/12.6/hpcx/hpcx-2.21/ucc
CUDA_PATH=$(NVHPC_ROOT)/cuda

CUDA_CAL_PATH       = $(NVHPC_ROOT)/math_libs/12.6/targets/x86_64-linux
CUBLASMP_PATH       = $(NVHPC_ROOT)/math_libs/12.6/include
CUSOLVERMP_PATH     = $(HOME)/src/cusolvermp0.7.1/usr/include

CUDA_LIB_PATH       = $(NVHPC_ROOT)/cuda/lib64
CUSOLVER_LIB_PATH   = $(NVHPC_ROOT)/math_libs/12.6/lib64
CUBLASMP_LIB_PATH   = $(NVHPC_ROOT)/math_libs/12.6/lib64
CUSOLVERMP_LIB_PATH = $(HOME)/src/cusolvermp0.7.1/usr/lib64/libcusolvermp/12
CUDAMP_LIB_PATH     = $(NVHPC_ROOT)/math_libs/12.6/lib64

CUMP_INCLUDES     = -I$(MPI_PATH)/include  -I$(UCC_PATH)/include -I$(CUSOLVERMP_PATH) -I$(CUDA_PATH)/include  -I$(CUDA_CAL_PATH)/include
CUMP_LIBS         = -Wl,-rpath,$(CUSOLVERMP_LIB_PATH) -L$(CUDAMP_LIB_PATH) $(CUSOLVERMP_LIB_PATH)/libcusolverMp.so -lcusolver -lcublasmp -lcublas -L$(CUDA_CAL_PATH)/lib -lcal -L$(CUDA_LIB_PATH) -lcudart -L$(MPI_PATH)/lib -lmpi -L$(UCC_PATH)/lib/ -lucc


#NCCL_LIBS = -L/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/23.3/comm_libs/nccl/lib -lnccl
#NCCL_LIBS = -L/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.5/comm_libs/nccl/lib -lnccl
#NCCL_LIBS = -L/projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/comm_libs/nccl/lib -lnccl
NCCL_LIBS = -L$(NVHPC_ROOT)/comm_libs/nccl/lib -lnccl
#NCCL_LIBS = -Wl,-rpath=/projects/darwin-nv/rhel8/x86_64/packages/nvhpc//Linux_x86_64/24.7/comm_libs/nccl/lib/libnccl.so.2

# Linker
AR := ar rcs
LD := $(FC) $(CULIBS) $(NCCL_LIBS) $(BLASLIB) $(MPFLIB) $(OPT)

# Cleanup
RM := rm -f


