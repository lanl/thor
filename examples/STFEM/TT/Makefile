# Compiler
FC = mpifort

# Optimization & safety/debug flags
OPT = -cpp -g -O2 -fopenmp -fallow-argument-mismatch \
      -ffpe-trap=overflow,underflow,invalid -Wall -Wextra \
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
SRC_CDR = TT_3D_CDR.f90
SRC_POISSON = TT_3D_Poisson.f90

# Object files
OBJ_CDR = $(SRC_CDR:.f90=.o)
OBJ_POISSON = $(SRC_POISSON:.f90=.o)

# Executables
EXE_CDR = TT_3D_CDR.exe
EXE_POISSON = TT_3D_Poisson.exe

# Default rule builds both
all: $(EXE_CDR) $(EXE_POISSON)

# Build rules
$(EXE_CDR): $(OBJ_CDR)
	$(FC) $(OPT) $(MOD_FLAGS) $(OBJ_CDR) $(LIBS) $(BLASLIB) -o $(EXE_CDR)

$(EXE_POISSON): $(OBJ_POISSON)
	$(FC) $(OPT) $(MOD_FLAGS) $(OBJ_POISSON) $(LIBS) $(BLASLIB) -o $(EXE_POISSON)

%.o: %.f90
	$(FC) $(OPT) $(MOD_FLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJ_CDR) $(OBJ_POISSON) $(EXE_CDR) $(EXE_POISSON)
