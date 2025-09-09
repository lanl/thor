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
#CULIBS =
CUFC = mpif90
CUFLAGS = -O2 -Minfo -Mcudalib=cusolver,cublas -Kieee -acc -Mcuda
#FC := pgf90
LD := $(FC) $(CULIBS) $(BLASLIB) $(MPFLIB) $(OPT)


# List of all source files
SRCS := src/utils/rnd.f90\
	src/utils/default.f90\
	src/utils/string.f90\
	src/utils/timef.f90\
	src/utils/matlab_struct_module.f90\
	src/utils/say.f90\
	src/utils/nan.f90\
	src/utils/linalg_utils.f90\
	src/utils/mat.f90\
	src/utils/matrix_util.f90\
	src/utils/ort.f90\
	src/utils/lr.f90\
	src/utils/trans.f90\
	src/utils/ttlocsolve.f90\
	src/thor_pointers.f90\
    src/thor.f90\
	src/ttop.f90\
	src/thorio.f90\
	src/qtt.f90\
	src/ttnodeop.f90\
	src/cross/dmrgg.f90\
	src/solvers/ttamen.f90\

CUSRCS := src/utils/cuda_solvers.f90\

TEST_SRCS := \
	tests/test_thor_tensor_zeros.f90\
	tests/test_thor_tensor_ones.f90\
    tests/test_thor_tensor_api.f90\
	tests/test_thor_tensor_rand.f90\
	tests/test_thor_tensor_ops.f90\
	tests/test_thor_tensor_inv.f90\
	tests/test_thor_matrix_constructor.f90\
	tests/test_thor_matrix_ops.f90\
	tests/test_mat_lib.f90\
	tests/test_rnd_lib.f90\
	tests/test_amen_ops.f90\
	tests/test_amen_mm.f90\
	tests/test_amen_mv.f90\
	tests/test_thor_io_lib.f90\
	tests/test_qtt_tensor_ops.f90\

BENCH_SRCS := \
	bench/dtt_ones.f90 \
	bench/dtt_tensor.f90 \
	bench/dtt_tensor_rounding.f90 \
	bench/dtt_amen.f90 \
	bench/dtt_matmul.f90 \
	bench/dmrgg_inv.f90 \
	bench/low_rank.f90  \
	bench/rndlib.f90  \
	bench/matmuls.f90  \

CUTEST_SRCS := tests/test_cuda_solvers.f90\

# Create lists of the build artefacts in this project
OBJS := $(SRCS:.f90=.o)
CUOBJS := $(addsuffix .o, $(CUSRCS))
TEST_OBJS := $(TEST_SRCS:.f90=.o)
BENCH_OBJS := $(BENCH_SRCS:.f90=.o)
CUTEST_OBJS := $(addsuffix .o, $(CUTEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
LIBWCUDA := $(patsubst %, lib%_cuda.a, $(NAME))
TEST_EXE := $(TEST_SRCS:.f90=.exe)
BENCH_EXE := $(BENCH_SRCS:.f90=.exe)
CUTEST_EXE := $(patsubst %.f90, %.exe, $(CUTEST_SRCS))

# Declare all public targets
.PHONY: all clean cuda
all: $(LIB) $(TEST_EXE) # $(EXAMPLE_EXE) 

bench: $(BENCH_EXE)

cuda: $(CUTEST_EXE) $(LIBWCUDA)

$(TEST_EXE) $(BENCH_EXE) : $(LIB)

test : $(TEST_EXE)
	@bash -c 'passed=0; total=0;'\
	'echo "Running a test suite:";'\
	'for t in $(TEST_EXE); do '\
	'  $(MPIRUN) -n 1 ./$$t >/dev/null 2>/dev/null;'\
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

$(LIBWCUDA) : $(OBJS) $(CUOBJS)
	$(AR) $@ $^

# Define all module interdependencies
rnd_lib.mod := src/utils/rnd.o
string_lib.mod := src/utils/string.o
say_lib.mod := src/utils/say.o
nan_lib.mod := src/utils/nan.o
linalg_module.mod := src/utils/linalg_utils.o
matlab_struct_module.mod := src/utils/matlab_struct_module.o
mat_lib.mod := src/utils/mat.o
matrix_util.mod := src/utils/matrix_util.o
ort_lib.mod := src/utils/ort.o
lr_lib.mod := src/utils/lr.o
time_lib.mod := src/utils/timef.o
default_lib.mod := src/utils/default.o
trans_lib.mod := src/utils/trans.o
ttlocsolve_lib.mod := src/utils/ttlocsolve.o
cuda_solvers.mod := src/utils/cuda_solvers.f90.o
thor_pointers.mod := src/thor_pointers.o
thor_lib.mod := src/thor.o
ttnodeop_lib.mod := src/ttnodeop.o
ttop.mod := src/ttop.o
thor_io.mod := src/thorio.o
dmrgg_lib.mod := src/cross/dmrgg.o
qtt_lib.mod := src/qtt.o
ttamen_lib.mod := src/solvers/ttamen.o

src/utils/linalg_utils.o: $(string_lib.mod)
src/utils/mat.o: $(nan_lib.mod) $(linalg_module.mod)
src/utils/ort.o: $(nan_lib.mod)
src/utils/lr.o: $(ort_lib.mod)
src/utils/ttlocsolve.o: $(trans_lib.mod)
src/thor.o: $(rnd_lib.mod) $(say_lib.mod) $(lr_lib.mod) $(mat_lib.mod) $(thor_pointers.mod) $(matlab_struct_module.mod)
src/ttnodeop.o: $(thor_pointers.mod) $(thor_lib.mod) $(default_lib.mod) $(trans_lib.mod) $(mat_lib.mod) $(ttlocsolve_lib.mod)
src/ttop.o: $(qtt_lib.mod)
src/thorio.o: $(thor_lib.mod) $(string_lib.mod)
src/qtt.o: $(thor_lib.mod)
src/cross/dmrgg.o: $(rnd_lib.mod) $(default_lib.mod) $(time_lib.mod) $(lr_lib.mod) $(thor_pointers.mod) $(thor_lib.mod)
src/solvers/ttamen.o : $(thor_pointers.mod) $(thor_lib.mod) $(ttlocsolve_lib.mod) $(ttnodeop_lib.mod) $(default_lib.mod)
tests/test_cuda_solvers.f90.o : $(cuda_solvers.mod)

# Create cuda object files from Fortran source
$(CUOBJS) $(CUTEST_OBJS): %.o: %
	$(CUFC) $(CUFLAGS) -c -o $@ $<

$(CUTEST_EXE): %.exe: %.f90.o $(LIBWCUDA)
	$(CUFC) $(CUFLAGS) -o $@ $^

# Suffix rules
.SUFFIXES: .f90 .o .exe

# %.f90 -> %.exe for the test executables

.f90.exe:
	$(FC) $(OPT) -o $@ $^ $(LIB) $(MPFLIB) $(BLASLIB)

# Create object files from Fortran source
.f90.o:
	$(FC) $(OPT) -c -o $@ $<

src/solvers/ttamen.f90.o : $(thor_pointers.mod) $(thor.mod) $(ttlocsolve_lib.mod) $(ttnodeop_lib.mod) $(default_lib.mod)

# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(CUOBJS) $(TEST_OBJS) $(BENCH_OBJ) $(CUTEST_OBJS)) $(filter %.exe, $(TEST_EXE) $(CUTEST_EXE)) $(LIB) $(LIBWCUDA) $(wildcard *.mod) $(wildcard *.o)

