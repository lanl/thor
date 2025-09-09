# Compiler
FC = mpifort

# Optimization & safety/debug flags (from original)
OPT = -cpp -g -O2 -fopenmp -fallow-argument-mismatch \
      -ffpe-trap=overflow,underflow,invalid -Wall -Wextra \
      -Wno-uninitialized -Wno-unused-variable -Wno-unused-parameter \
      -fbacktrace -fimplicit-none -fdefault-double-8 \
      -ffree-line-length-none -fexternal-blas

# BLAS/LAPACK (you must define OPENBLAS_DIR in your shell environment)
BLASLIB = -I$(OPENBLAS_DIR)/include -L$(OPENBLAS_DIR)/lib -llapack -lblas

# Source files
SRC = main.f90 

# Object files
OBJ = $(SRC:.f90=.o)

# Output executable
EXE = 01-hconvergence.exe

# Directories with .mod files
MOD_DIRS = ../../fortran/utils/stdlib/build/gfortran_BDF09C972676FF1B ../../fortran
MOD_FLAGS = $(foreach dir,$(MOD_DIRS),-I$(dir))

# Static libraries to link (compiled previously)
LIBS = ../../fortran/utils/stdlib/build/gfortran_37DB36B31CAE42A6/stdlib/libstdlib.a ../../fortran/libthor.a

# Rules
all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(OPT) $(MOD_FLAGS) $(OBJ) $(LIBS) $(BLASLIB) -o $(EXE)

%.o: %.f90
	$(FC) $(OPT) $(MOD_FLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE)
