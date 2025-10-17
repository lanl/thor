#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "geom.h"
#include "h5aux.h"
#include "xdmf.h"
#include "minip.h"
#include "bdry.h"
#include "eos.h"
#include "update.h"
#include "vof2d.h"
#include "util.h"
#include "mesh.h"
#include "mat.h"
#include "tests.h"

void random_field_3dmesh(const int *dims, double ****field, 
                         const double xmin, const double xmax) {

    double **field2d, *field1d;
    *field = NULL;
    field2d = NULL;
    field1d = NULL;

    ASSIGN_3D_FORM(double, (*field), field2d, field1d, dims)      
    
    for (int k = 0; k < dims[2]; k++)
    for (int j = 0; j < dims[1]; j++)
    for (int i = 0; i < dims[0]; i++)
        (*field)[k][j][i] = random_double(xmin, xmax);

}


void random_field_3dmesh_int(const int *dims, int ****field, 
                             const int imin, const int imax) {
    int **field2d, *field1d;
    *field = NULL;
    field2d = NULL;
    field1d = NULL;

    ASSIGN_3D_FORM(int, (*field), field2d, field1d, dims)      
    
    for (int k = 0; k < dims[2]; k++)
    for (int j = 0; j < dims[1]; j++)
    for (int i = 0; i < dims[0]; i++)
        (*field)[k][j][i] = random_int(imin, imax);
}


// Function to allocate a 3D array (e.g. for velocities)
void allocate_velocity_3dmesh(const int *dims, double *****vel3d) {

    double ***v3d, **v2d, *v1d;
    (*vel3d) = NULL;
    v3d = NULL;
    v2d = NULL; 
    v1d = NULL;
    const int sizes[4] = {dims[0], dims[1], dims[2], 3};
     sizes[0], 
     sizes[1], 
     sizes[2], 
     sizes[3]);

    ASSIGN_4D_FORM(double, sizes, (*vel3d), v3d, v2d, v1d)

}


