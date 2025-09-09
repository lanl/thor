# Compiler
FC = mpifort

# Optimization & safety/debug flags
OPT = -cpp -g -O2 -fopenmp -Wno-argument-mismatch \
      -ffpe-trap=overflow -Wall -Wextra \
      -Wno-uninitialized -Wno-unused-variable -Wno-unused-parameter \
      -fbacktrace -fimplicit-none -fdefault-double-8 \
      -ffree-line-length-none -fexternal-blas

# BLAS/LAPACK
BLASLIB = -I$(OPENBLAS_DIR)/include -L$(OPENBLAS_DIR)/lib -llapack -lblas

# Directories with .mod files
MOD_DIRS = ../../../fortran
MOD_FLAGS = $(foreach dir,$(MOD_DIRS),-I$(dir))

# Static libraries to link
LIBS = ../../../fortran/libthor.a

# Source files
SRC = TT_3D_Burger.f90

# Object files
OBJ = $(SRC:.f90=.o)

# Executables
EXE = TT_3D_Burger.exe

# Default rule builds both
all: $(EXE) 

# Build rules
$(EXE): $(OBJ) Makefile
	$(FC) $(OPT) $(MOD_FLAGS) $(OBJ) $(LIBS) $(BLASLIB) -o $(EXE)


%.o: %.f90
	$(FC) $(OPT) $(MOD_FLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJ) $(EXE)
