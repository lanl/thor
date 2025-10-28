# Project name
NAME := thor

# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# System-specific include
ifeq (,$(wildcard ./make.inc))
    $(error "make.inc is missing; run '. load.sh <option>'")
endif
include make.inc

# CUDA LIB FLAGS
CULIBS = -stdpar=gpu,multicore

#CUFLAGS = -O2 -Minfo -Mcudalib=cusolver,cublas -Kieee -acc -Mcuda $(NCCL_LIBS)
CUFLAGS =  -O2 -Minfo -cudalib=cusolver,cublas -Kieee -acc -cuda -cpp $(NCCL_LIBS)
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)
MPIRUN := /projects/darwin-nv/rhel8/x86_64/packages/nvhpc/Linux_x86_64/24.7/comm_libs/mpi/bin/mpirun

# List of all source files
SRCS := Tensor_Network_Integrators/fortran/src/utils/rnd.f90\
    Tensor_Network_Integrators/fortran/src/utils/string.f90\
    Tensor_Network_Integrators/fortran/src/utils/timef.f90\
    Tensor_Network_Integrators/fortran/src/utils/matrix_util.f90\
    Tensor_Network_Integrators/fortran/src/utils/cell_arrays_module.f90\
    Tensor_Network_Integrators/fortran/src/utils/nan.f90\
    Tensor_Network_Integrators/fortran/src/utils/linalg_utils.f90\
    Tensor_Network_Integrators/fortran/src/utils/mat.f90\
    src/communicators/comm_types.f90\
    src/communicators/mpi_allreduce.f90\
    src/communicators/mpi_bcast.f90\
    src/communicators/mpi_recv.f90\
    src/communicators/mpi_send.f90\
    src/communicators/communicators.f90\
    src/communicators/distributed_comm.f90\

CUMP_C_SRCS := cuda_interface/cuda_cusolvermp_aux.c\
               cuda_interface/cuda_cal_aux.c\
               cuda_interface/cuda_cublasmp_aux.c\
               cuda_interface/cuda_cublas_aux.c\

CUMP_F_SRCS := cuda_interface/cuda_comm_aux.f90\
               cuda_interface/cublas_aux.f90\
               cuda_interface/cublasmp_aux.f90\
               cuda_interface/cusolvermp_aux.f90\
               cuda_interface/mat_ops_kernels.f90\
               src/structs/dray.f90\
               src/thor.f90\
               src/thorio.f90\

TEST_SRCS :=  \
    tests/test_cuda_comm_aux.f90\
    tests/test_cusolvermp_aux.f90\
    tests/test_thor_super_comm.f90\
    tests/test_mat_ops_kernels.f90\
    tests/test_dray.f90\
    tests/test_cublas_aux.f90\
    tests/test_cublasmp_aux.f90\
    tests/test_dthor_constructor.f90\

# Create lists of the build artefacts in this project
OBJS := $(SRCS:.f90=.o)
#CUMP_C_OBJS := $(addsuffix .o, $(CUMP_C_SRCS))
#CUMP_F_OBJS := $(addsuffix .o, $(CUMP_F_SRCS))
CUMP_C_OBJS := $(CUMP_C_SRCS:.c=.o)
CUMP_F_OBJS := $(CUMP_F_SRCS:.f90=.o)
TEST_OBJS := $(TEST_SRCS:.f90=.o)
LIB := $(patsubst %, lib%.a, $(NAME))
LIBWCUDAMP := $(patsubst %, lib%_cuda_mp.a, $(NAME))
TEST_EXE := $(TEST_SRCS:.f90=.exe)

# Declare all public targets
.PHONY: all clean cuda purge cudamp
all: cudamp $(LIB) $(TEST_EXE) # $(EXAMPLE_EXE)

cudamp: purge $(LIBWCUDAMP) $(CUMP_TEST_EXE)

$(TEST_EXE) : $(LIB)

test : $(TEST_EXE)
	@bash -c 'passed=0; total=0;'\
	'echo "Running a test suite:";'\
	'for t in $(TEST_EXE); do '\
	'  $(MPIRUN) -np 1 ./$$t >/dev/null 2>/dev/null;'\
	'  j=$$[$$?>0];'\
	'  if [ $$j -eq 0 ]; then '\
	'    printf "  %-50s : passed\n" $$t;'\
	'    passed=$$[ $$passed + 1 ];'\
	'  else '\
	'    printf "  %-50s : FAILED\n" $$t;'\
	'  fi;'\
	'  total=$$[ $$total + 1 ];'\
	'done;'\
	'echo "- - -";'\
	'echo Summary: $$passed/$$total tests passed'

# Create the static library from the object files
$(LIB): $(OBJS)
	$(AR) $@ $^

$(LIBWCUDAMP) : $(OBJS) $(CUMP_C_OBJS) $(CUMP_F_OBJS)
	$(AR) $@ $^

# Define all module interdependencies
cuda_comm_aux.mod := cuda_interface/cuda_comm_aux.o
cublas_aux.mod := cuda_interface/cublas_aux.o
cublasmp_aux.mod := cuda_interface/cublasmp_aux.o
cusolvermp_aux.mod := cuda_interface/cusolvermp_aux.o
matrix_ops.mod := cuda_interface/mat_ops_kernels.o

comm_types.mod := src/communicators/comm_types.o
comm_allreduce.mod := src/communicators/mpi_allreduce.o
comm_bcast.mod := src/communicators/mpi_bcast.o
comm_recv.mod := src/communicators/mpi_recv.o
comm_send.mod := src/communicators/mpi_send.o
communicators.mod := src/communicators/communicators.o
distributed_comms.mod := src/communicators/distributed_comm.o
distributed_arrays.mod := src/structs/dray.o

string_lib.mod := Tensor_Network_Integrators/fortran/src/utils/string.o
rnd_lib.mod := Tensor_Network_Integrators/fortran/src/utils/rnd.o
nan_lib.mod := Tensor_Network_Integrators/fortran/src/utils/nan.o
linalg_module.mod := Tensor_Network_Integrators/fortran/src/utils/linalg_utils.o
cell_arrays_module.mod := Tensor_Network_Integrators/fortran/src/utils/cell_arrays_module.o
mat_lib.mod := Tensor_Network_Integrators/fortran/src/utils/mat.o
matrix_util.mod := Tensor_Network_Integrators/fortran/src/utils/matrix_util.o
time_lib.mod := Tensor_Network_Integrators/fortran/src/utils/timef.o
thor.mod := src/thor.o
dthorio_lib.mod := src/thorio.o

$(linalg_module.mod)      : $(string_lib.mod)
$(mat_lib.mod)            : $(nan_lib.mod) $(linalg_module.mod)
$(matrix_util.mod)        : $(mat_lib.mod)
$(comm_allreduce.mod) $(comm_bcast.mod) $(comm_send.mod) $(comm_recv.mod) : $(comm_types.mod)
$(communicators.mod)      : $(comm_types.mod) $(comm_allreduce.mod) $(comm_bcast.mod) $(comm_recv.mod) $(comm_send.mod)
$(distributed_comms.mod)  : $(string_lib.mod) $(comm_types.mod) $(communicators.mod) $(cuda_comm_aux.mod) $(cublasmp_aux.mod) $(cusolvermp_aux.mod)
$(distributed_arrays.mod) : $(distributed_comms.mod)  $(cusolvermp_aux.mod) $(cublasmp_aux.mod) $(cublas_aux.mod) $(matrix_ops.mod) $(matrix_util.mod)
$(thor.mod)               : $(cell_arrays_module.mod) $(distributed_comms.mod) $(distributed_arrays.mod)
$(dthorio_lib.mod)        : $(thor.mod) $(distributed_comms.mod) $(distributed_arrays.mod)
$(cublasmp_aux.mod) $(cusolvermp_aux.mod) : $(cuda_comm_aux.mod)
#
tests/test_thor_super_comm.f90.o : $(distributed_comms.mod) $(cuda_comm_aux.mod) $(cusolvermp_aux.mod) $(cublasmp_aux.mod)
tests/test_dray.f90.o : $(distributed_comms.mod) $(matrix_ops.mod)


$(CUMP_C_OBJS): %.o: %.c
	$(CUCC) $(CUMP_INCLUDES) -c -o $@ $<

$(CUMP_F_OBJS): %.o: %.f90
	$(CUFC) $(CUFC_FLAGS) $(CUMP_INCLUDES) -c -o $@ $<

$(TEST_EXE): %.exe: %.f90  $(LIBWCUDAMP)
	$(CUFC) $(CUMP_INCLUDES) $(CUMP_LIBS) $(BLASLIB) -o $@ $^ $(CUFC_FLAGS) $(CUFLAGS)


# Suffix rules
.SUFFIXES: .f90 .o .exe

# %.f90 -> %.exe for the test executables

.f90.exe:
	$(FC) $(OPT) -o $@ $^ $(LIB) $(MPFLIB) $(BLASLIB)

# Create object files from Fortran source
.f90.o:
	$(FC) $(OPT) -c -o $@ $<

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) core.* $(filter %.o, $(OBJS) $(CUMP_C_OBJS) $(CUMP_F_OBJS) $(TEST_OBJS) $(BENCH_OBJ) $(CUMP_TEST_OBJS)) $(filter %.exe, $(TEST_EXE) $(CUMP_TEST_EXE)) $(LIB) $(LIBWCUDAMP) $(wildcard *.mod) $(wildcard *.o)

purge:
	find . -name "core.*" -exec rm -rf \{\} \;
