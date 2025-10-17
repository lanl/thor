CC = gcc
MPI_INCLUDE = .

# Adjust accordingly 
HDF5_DIR=${HOME}/HOMEBREW

INC = -I. -I${HDF5_DIR}/include -I${MPI_INCLUDE}
LIBPATH = -L. -L${HDF5_DIR}/lib
FLAGS = -g -fPIC

# Toggle H5DUMP_DEBUG using a Makefile variable
# ifeq ($(H5DEBUG),1)
# FLAGS += -DH5DUMP_DEBUG
# $(info 🔧 Compiling with H5DUMP_DEBUG enabled)
# endif

HDRS = minip.h globals.h mesh.h geom.h vof3d.h vof3d.h \
	xdmf.h h5aux.h util.h mat.h h5aux.h xdmf.h eos.h debug.h

# if ('uname' = "Darwin") set CC = gcc
all : main.x output libvof3d.so #record_tests.x

#main.x: main.o control.o init.o geom.o util.o eos.o mesh.o update.o \
	bdry.o advection.o vof3d.o vof2d.o mat.o bdry.o xdmf.o h5aux.o io.o libmio_serial.a \
	${HDRS} Makefile
#	${CC} ${FLAGS} -o main.x main.o control.o init.o eos.o geom.o \
				util.o update.o bdry.o advection.o \
				vof3d.o vof2d.o mesh.o mat.o xdmf.o h5aux.o io.o \
				$(INC) $(LIBPATH) \
				-lmio_serial -lhdf5 -lm
main.x: main.o control.o init.o geom.o advection.o util.o mesh.o update.o \
	eos.o bdry.o vof3d.o vof2d.o mat.o eos.o xdmf.o h5aux.o io.o io_asci.o \
	${HDRS} Makefile
	${CC} ${FLAGS} -o main.x main.o control.o init.o geom.o util.o update.o \
	   eos.o bdry.o advection.o vof3d.o vof2d.o mesh.o mat.o xdmf.o h5aux.o io.o \
	 $(INC) $(LIBPATH) -lhdf5 -lm

record_tests.x: testdata tests3d.o tests.o control.o init.o geom.o util.o eos.o mesh.o update.o \
	bdry.o advection.o vof3d.o vof2d.o mat.o bdry.o xdmf.o io.o libmio_serial.a \
	${HDRS} Makefile
	${CC} ${FLAGS} -o record_tests.x tests3d.o tests.o control.o init.o eos.o geom.o \
				util.o update.o bdry.o advection.o \
				vof3d.o vof2d.o mesh.o mat.o xdmf.o h5aux.o io.o \
				$(INC) $(LIBPATH) \
				-lmio_serial -lhdf5 -lm

libvof3d.so: vof3d_lib.o util.o 
	${CC} -shared -o $@ $^


output:
	mkdir -p output

testdata:
	mkdir -p testdata

libmio_serial.a : mio.o bio.o
	ar rcs $@ $^

.PHONY: all clean

clean :
	rm -rf output cs_output.h5 cs_output.xdmf matlab/output.* *.o libmio_serial.a main.x record_tests.x file_000* file_initial libvof3d.so

.SUFFIXES: .c .o
.c.o:
	$(CC) $(FLAGS) $(INC) -c -o $(notdir $*.o) $(INC) $<

