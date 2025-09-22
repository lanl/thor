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
#include "control.h"
#include "mat.h"
#include "tests.h"


//#define VERBOSE

double random_int(const int imin, const int imax) {
    return imin + rand() % (imax - imin + 1);
}


double random_double(const double xmin, const double xmax) {
    return xmin + (xmax - xmin)*((double)rand() / (double)RAND_MAX);
}


// Function to create a random mesh
void create_random_mesh(int dim, int *ncell, int *nbdry, int *ncell_ext, int *nnode, int *nnode_ext,
                        double *xl, double *xr, double *dx, const int verbose) {
    *nbdry = random_int(0, 4);
    for (int i = 0; i < dim; i++) {
        ncell[i] = random_int(1, 10);  // Random number of cells between 1 and 10
        ncell_ext[i] = ncell[i] + 2*(*nbdry);
        nnode[i] = ncell[i] + 1;
        nnode_ext[i] = ncell_ext[i] + 1;
    }

    // Generate random lower and upper bounds (xl and xr) and cell sizes (dx)
    for (int i = 0; i < dim; i++) {
        xl[i] = random_double(-10.0, 10.0);
        xr[i] = random_double(xl[i], xl[i] + 10.0);
        dx[i] = (xr[i] - xl[i])/(double)(ncell[i]);
    }

    if (verbose) {
        // Print the generated mesh properties
        switch (dim) {
        case 2:
            printf("Generated 2D Mesh:\n");
            printf(" - number of cells (ncell): %d x %d\n", ncell[0], ncell[1]);
            printf(" - number of boundary layers (nbdry): %d\n", *nbdry);
            printf(" - lower bounds (xl): %f, %f\n", xl[0], xl[1]);
            printf(" - upper bounds (xr): %f, %f\n", xr[0], xr[1]);
            printf(" - cell sizes (dx): %f, %f\n", dx[0], dx[1]);
            break;

        case 3:
            printf("Generated 3D Mesh:\n <<< TODO >>> \n");
            printf(" - number of cells (ncell): %d x %d x %d\n", ncell[0], ncell[1], ncell[2]);
            printf(" - number of boundary layers (nbdry): %d\n", *nbdry);
            printf(" - lower bounds (xl): %f, %f, %f\n", xl[0], xl[1], xl[2]);
            printf(" - upper bounds (xr): %f, %f, %f\n", xr[0], xr[1], xr[2]);
            printf(" - cell sizes (dx): %f, %f, %f\n", dx[0], dx[1], dx[2]);
            break;

        default:
            printf("Generated %d Mesh.\n", dim);
        }
    }

}

void random_field_2dmesh(const int *dims, double ***field, const double xmin, const double xmax) {

    int i, j;

    *field = (double **) malloc(dims[1] * sizeof(double *));
    for (j = 0; j < dims[1]; j++) {
        if (j == 0)
            (*field)[0] = (double *) malloc(dims[0]*dims[1] * sizeof(double));
        else
            (*field)[j]  =  (*field)[j-1] + dims[0];
        for (i = 0; i < dims[0]; i++) {
            (*field)[j][i] = random_double(xmin, xmax);
        }
    }
}


void random_field_2dmesh_int(const int *dims, int ***field, const int imin, const int imax) {
    int i, j;

    *field = (int **) malloc(dims[1] * sizeof(int *));
    for (j = 0; j < dims[1]; j++) {
        if (j == 0)
            (*field)[0] = (int *) malloc(dims[0]*dims[1] * sizeof(int));
        else
            (*field)[j]  =  (*field)[j-1] + dims[0];
        for (i = 0; i < dims[0]; i++) {
            (*field)[j][i] = random_int(imin, imax);
        }
    }
}


// Function to allocate a 3D array (e.g. for velocities)
void allocate_velocity_2dmesh(const int *dims, double ****vel2d) {

    const int dim = 2;
    int i, j, offset, lsize_node;
    double **v2d, *v1d;

    lsize_node = dims[0]*dims[1];
    (*vel2d) = (double ***) malloc(dims[1] * sizeof(double **));
    v2d     = (double **) malloc(lsize_node * sizeof(double *));
    v1d     = (double  *) malloc(lsize_node * dim * sizeof(double));
    offset = 0;

    for (j = 0; j < dims[1]; j++) {
        v2d[0] = v1d + offset;
        for (i = 1; i < dims[0]; i++) {
            v2d[i] = v2d[i-1] + dim;
        }
        (*vel2d)[j] = v2d;
        offset += (dims[0] * dim);
        v2d    +=  dims[0];
    }

}


void copy_2dvel_components(const int *dims, double * const * const vx,
        double * const * const vy, double ***vel2d) {
    int i, j;
    for (j = 0; j < dims[1]; ++j) {
        for (i = 0; i < dims[0]; ++i) {
            vel2d[j][i][0] = vx[j][i];
            vel2d[j][i][1] = vy[j][i];
        }
    }
}

void copy_components_from_2dvel(const int *dims,
        double *** const vel2d, double **vx, double **vy) {
    int i, j;
    for (j = 0; j < dims[1]; ++j) {
        for (i = 0; i < dims[0]; ++i) {
            vx[j][i] = vel2d[j][i][0];
            vy[j][i] = vel2d[j][i][1];
        }
    }
}


void generate_random_meshes(const int dim, const int num_meshes,
        const char *meshes_fname_h5, const int verbose) {

    int imesh, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    double xl[3], xr[3], dx[3];
    double **rho_for_2dcell, **ei_for_2dcell, **pres_for_2dcell, **cs_for_2dcell;
    double **vel_x, **vel_y, **vav_x, **vav_y, ***vel_for_2dnode, ***velav_for_2dnode;
    int **nmat_for_2dcell, **matid_for_2dcell, i, j;

    delete_file_if_exists(meshes_fname_h5);
    hid_t file_meshes_id = h5_new_file(meshes_fname_h5);

    char group_name[50];
    hid_t mesh_gid, attr_id, test_gid;

    for (imesh = 0; imesh < num_meshes; ++imesh) {
        // Create a random 2D mesh
        create_random_mesh(dim, ncell, &nbdry, ncell_ext, nnode, nnode_ext, xl, xr, dx, verbose);

        // Create a new group for each mesh
        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // Attribute: dim (single scalar)
        h5_write_attr_int(mesh_gid, "dim", dim);
        h5_write_attr_int(mesh_gid, "nbdry", nbdry);
        h5_write_attr_1d(mesh_gid, "xl", xl, dim);
        h5_write_attr_1d(mesh_gid, "xr", xr, dim);
        h5_write_attr_1d(mesh_gid, "dx", dx, dim);
        h5_write_attr_1d_int(mesh_gid, "ncell", ncell, dim);

        if (dim == 2) {
            // Mesh fields
            random_field_2dmesh (ncell_ext, &rho_for_2dcell, 0., 10.);
            random_field_2dmesh (ncell_ext, &ei_for_2dcell, 0., 10.);
            random_field_2dmesh (ncell_ext, &pres_for_2dcell, 0., 10.);
            random_field_2dmesh_int (ncell_ext, &nmat_for_2dcell, 1, 5);
            random_field_2dmesh_int (ncell_ext, &matid_for_2dcell, 0, 4);
            random_field_2dmesh(ncell_ext, &cs_for_2dcell, 0., 10.);


            h5_write_2d(mesh_gid, "rho", rho_for_2dcell, ncell_ext);
            h5_write_2d(mesh_gid, "ei",   ei_for_2dcell, ncell_ext);
            h5_write_2d(mesh_gid, "pres", pres_for_2dcell, ncell_ext);
            h5_write_2d_int(mesh_gid, "nmat", nmat_for_2dcell, ncell_ext);
            h5_write_2d_int(mesh_gid, "matid", matid_for_2dcell, ncell_ext);
            h5_write_2d(mesh_gid, "cs", cs_for_2dcell, ncell_ext);


            // Velocities
            random_field_2dmesh (nnode_ext, &vel_x, -10., 10.);
            random_field_2dmesh (nnode_ext, &vel_y, -10., 10.);
            allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
            copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);
            h5_write_2d(mesh_gid, "vel_x", vel_x, nnode_ext);
            h5_write_2d(mesh_gid, "vel_y", vel_y, nnode_ext);

            random_field_2dmesh (nnode_ext, &vav_x, -10., 10.);
            random_field_2dmesh (nnode_ext, &vav_y, -10., 10.);
            allocate_velocity_2dmesh (nnode_ext, &velav_for_2dnode);
            copy_2dvel_components (nnode_ext, vav_x, vav_y, velav_for_2dnode);
            h5_write_2d(mesh_gid, "vav_x", vav_x, nnode_ext);
            h5_write_2d(mesh_gid, "vav_y", vav_y, nnode_ext);

            // for (int j=0; j<ncell_ext[1]; ++j) {
            //     for (int i=0; i<ncell_ext[0]; ++i) {
            //         printf ("%2d ", nmat_for_2dcell[j][i]);
            //     }
            //     printf("\n");
            // }

            free(rho_for_2dcell[0]); free(rho_for_2dcell);
            free(ei_for_2dcell[0]); free(ei_for_2dcell);
            free(pres_for_2dcell[0]); free(pres_for_2dcell);
            free(nmat_for_2dcell[0]); free(nmat_for_2dcell);
            free(matid_for_2dcell[0]); free(matid_for_2dcell);
            free(vel_x[0]); free(vel_x);
            free(vav_x[0]); free(vav_x);
            free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
            free(velav_for_2dnode[0][0]); free(velav_for_2dnode[0]); free(velav_for_2dnode);
        }
        else {
            printf ("ERROR in generate_random_meshes: dim value not implemented\n");
            exit(-1);
        }

        // Cleanup: close the group, free memory
        H5Gclose(mesh_gid);
    }

    H5Fclose(file_meshes_id);
    printf("written: %s\n", meshes_fname_h5);

} // generate_random_meshes


void record_test__rec_rec(int N, const char *filename, const int verbose) {


    // Allocate memory for inputs and outputs
    int *ifinquiry = (int *)malloc(N * sizeof(int));
    int *ifmixed = (int *)malloc(N * sizeof(int));
    int *geop = (int *)malloc(N * sizeof(int));
    double *vcell = (double *)malloc(N * sizeof(double));
    int *dim = (int *)malloc(N * sizeof(int));
    double *xxl = (double *)malloc(N * 2 * sizeof(double));  // 2D array for lower bound
    double *xxr = (double *)malloc(N * 2 * sizeof(double));  // 2D array for upper bound
    double *xl = (double *)malloc(N * 2 * sizeof(double));   // 2D array for cell coordinates
    double *dx = (double *)malloc(N * 2 * sizeof(double));   // 2D array for cell size
    double *vol = (double *)malloc(N * sizeof(double));      // Output: cell volume

    // Randomly generate input arguments and call rec_rec
    for (int i = 0; i < N; i++) {
        ifinquiry[i] = random_int(0, 1);    // Random 0 or 1
        geop[i] =  random_int(0, 2);        // Random 0, 1, or 2
        vcell[i] = random_double(0, 10);    // Random double [0, 10]
        dim[i] = 2;                         // Fixed 2D for this example
        xxl[2*i]     = random_double(-10.0, 10.0);
        xxl[2*i + 1] = -1.0 + random_double(-0.1, 0.1);
        xxr[2*i]     = xxl[2*i] + random_double(0.1, 5.0);
        xxr[2*i + 1] = +1.0 + random_double(-0.1, 0.1);
        xl[2*i]      = random_double(-10.0, 10.0);
        xl[2*i + 1]  = -1.0 + random_double(-0.2, 0.2);
        dx[2*i]      = random_double(3.0, 8.0);
        dx[2*i + 1]  = random_double(0.8, 1.5);

        // Call the target function rec_rec
        rec_rec(ifinquiry[i], &ifmixed[i], geop[i], vcell[i], dim[i],
                &xxl[2*i], &xxr[2*i], &xl[2 * i], &dx[2 * i], &vol[i]);
        if (verbose) {
            if (i==0) printf("dx[0] = {%f, %f}\n", dx[0], dx[1]);
        }
    }

    // HDF5 file and dataset handles
    hid_t file_id, dataspace_N, dataspace_2N, dataset_id;
    hsize_t dims[2];
    int idims[2];
    dims[0] = N; dims[1] = 2;
    idims[0] = N; idims[1] = 2;

    // Create HDF5 file
    delete_file_if_exists(filename);
    file_id = h5_new_file(filename);

    // Write each input and output dataset to the file
    dataspace_N  = H5Screate_simple(1, dims, NULL);
    dataspace_2N = H5Screate_simple(2, dims, NULL);

    h5_write_1d_int(file_id, "ifinquiry", ifinquiry, N);
    h5_write_1d_int(file_id, "ifmixed", ifmixed, N);
    h5_write_1d_int(file_id, "geop", geop, N);
    h5_write_1d(file_id,     "vcell", vcell, N);
    h5_write_1d_int(file_id, "dim", dim, N);

    h5_write_2d(file_id, "xxl", &xxl, idims);
    h5_write_2d(file_id, "xxr", &xxr, idims);
    h5_write_2d(file_id, "xl",  &xl, idims);
    h5_write_2d(file_id, "dx",  &dx, idims);
    h5_write_1d(file_id, "vol", vol, N);

    H5Sclose(dataspace_N);
    H5Sclose(dataspace_2N);
    H5Fclose(file_id);

    // Report success
    printf ("written: %s\n", filename);

    // Free allocated memory
    free(ifinquiry);
    free(ifmixed);
    free(geop);
    free(vcell);
    free(dim);
    free(xxl);
    free(xxr);
    free(xl);
    free(dx);
    free(vol);
} // record_test__rec_rec


void record_test__bdry_cell_2d(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    double **rho_for_2dcell, **ei_for_2dcell, **pres_for_2dcell;
    int **nmat_for_2dcell, **matid_for_2dcell, i, j;
    Bdry_Type btype[2] = {bdry_transmitted, bdry_transmitted};

    char group_name[50];
    hid_t mesh_gid, attr_id, test_gid;
    int num_meshes, imesh, d;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);
        ei_for_2dcell   = h5_read_2d(mesh_gid, "ei",  ncell_ext);
        pres_for_2dcell = h5_read_2d(mesh_gid, "pres",ncell_ext);
        nmat_for_2dcell = h5_read_2d_int(mesh_gid, "nmat",ncell_ext);
        matid_for_2dcell = h5_read_2d_int(mesh_gid, "matid",ncell_ext);

        // Call the function under test
        bdry_cell_2d(ncell, nbdry, btype, btype, nmat_for_2dcell,
                     matid_for_2dcell, rho_for_2dcell, ei_for_2dcell, pres_for_2dcell);

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_2d(test_gid, "rho", rho_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "ei",   ei_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "pres", pres_for_2dcell, ncell_ext);
        h5_write_2d_int(test_gid, "nmat", nmat_for_2dcell, ncell_ext);
        h5_write_2d_int(test_gid, "matid", matid_for_2dcell, ncell_ext);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);
        free(rho_for_2dcell[0]); free(rho_for_2dcell);
        free(ei_for_2dcell[0]); free(ei_for_2dcell);
        free(pres_for_2dcell[0]); free(pres_for_2dcell);
        free(nmat_for_2dcell[0]); free(nmat_for_2dcell);
        free(matid_for_2dcell[0]); free(matid_for_2dcell);
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__bdry_cell_2d


void record_test__bdry_node_2d(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    double **vel_x, **vel_y, **vav_x, **vav_y, ***vel_for_2dnode, ***velav_for_2dnode;
    Bdry_Type btype[2] = {bdry_transmitted, bdry_transmitted};

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        vel_x = h5_read_2d(mesh_gid, "vel_x", nnode_ext);
        vel_y = h5_read_2d(mesh_gid, "vel_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
        copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);

        vav_x = h5_read_2d(mesh_gid, "vav_x", nnode_ext);
        vav_y = h5_read_2d(mesh_gid, "vav_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &velav_for_2dnode);
        copy_2dvel_components (nnode_ext, vav_x, vav_y, velav_for_2dnode);

        // Call the function under test
        bdry_node_2d(ncell, nbdry, btype, btype, vel_for_2dnode, velav_for_2dnode);

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        copy_components_from_2dvel (nnode_ext, vel_for_2dnode, vel_x, vel_y);
        copy_components_from_2dvel (nnode_ext, velav_for_2dnode, vav_x, vav_y);
        h5_write_2d(test_gid, "vel_x", vel_x, nnode_ext);
        h5_write_2d(test_gid, "vel_y", vel_y, nnode_ext);
        h5_write_2d(test_gid, "vav_x", vav_x, nnode_ext);
        h5_write_2d(test_gid, "vav_y", vav_y, nnode_ext);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);
        free(vel_x[0]); free(vel_x);
        free(vav_x[0]); free(vav_x);
        free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
        free(velav_for_2dnode[0][0]); free(velav_for_2dnode[0]); free(velav_for_2dnode);
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__bdry_node_2d


void record_test__sound_speed_solid(int N, const char *filename, const int verbose) {

    // Allocate memory for inputs and outputs
    double *rho = (double *)malloc(N * sizeof(double));
    double *cs  = (double *)malloc(N * sizeof(double));

    // Randomly generate input arguments and call rec_rec
    for (int i = 0; i < N; i++) {
        rho[i] = exp(random_double(-8.0, 12.0));  // uniformly random in log space

        // Call the target function rec_rec
        sound_speed_solid(rho[i], &(cs[i]));
        if (verbose) {
            printf("rho[%3d], cs[%3d] = {%14.7e, %14.7e}\n", i, i, rho[i], cs[i]);
        }
    }
    // Create HDF5 file
    delete_file_if_exists(filename);
    hid_t file_id = h5_new_file(filename);

    // HDF5 file and dataset handles
    h5_write_1d(file_id, "rho", rho, N);
    h5_write_1d(file_id, "cs", cs, N);

    // Report success
    printf ("written: %s\n", filename);

    // Free allocated memory
    free(rho);
    free(cs);
} // record_test__sound_speed_solid

void record_test__courant_from_cs(const char *meshes_h5, const char *test_h5, const int verbose){
    int dim, ncell[3], nbdry, ncell_ext[3]; 
    double dx[3], **cs_for_2dcell, dt, courant_cs;

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5); 
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
        }

        dt  = random_double(1.0e-9,1.0); // random double between 1e-9 and 1
        cs_for_2dcell  = h5_read_2d(mesh_gid, "cs", ncell_ext);

        // Call function under test
        courant_from_cs(dim, ncell, nbdry, dx, cs_for_2dcell, NULL, dt, &courant_cs);

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_1d(test_gid, "dt", &dt, 1);
        h5_write_1d(test_gid, "courant", &courant_cs, 1);

        H5Gclose(test_gid);

        H5Gclose(mesh_gid);

        free(cs_for_2dcell[0]); free(cs_for_2dcell);

    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);
}// record_test__courant_from_cs

void record_test__compute_divu(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    double **vel_x, **vel_y, ***vel_for_2dnode, ***ptr_divu_for_2dcell;
    double dx[3];

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        vel_x = h5_read_2d(mesh_gid, "vel_x", nnode_ext);
        vel_y = h5_read_2d(mesh_gid, "vel_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
        copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);

        // Call the function under test
        compute_divu(dim, ncell, nbdry, dx, vel_for_2dnode, NULL);

        // Get divu_out
        ptr_divu_for_2dcell = get_ptr_divu();

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        copy_components_from_2dvel (nnode_ext, vel_for_2dnode, vel_x, vel_y);
        h5_write_2d(test_gid, "vel_x", vel_x, nnode_ext);
        h5_write_2d(test_gid, "vel_y", vel_y, nnode_ext);
        h5_write_2d(test_gid, "divu", *ptr_divu_for_2dcell, ncell_ext);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);
        free(vel_x[0]); free(vel_x); free(vel_y[0]); free(vel_y);
        free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
        free((*ptr_divu_for_2dcell)[0]); free(*ptr_divu_for_2dcell);
        *ptr_divu_for_2dcell = NULL;

    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__compute_divu

void record_test__compute_qvis(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    double **vel_x, **vel_y, ***vel_for_2dnode, **rho_for_2dcell, **cs_for_2dcell;
    double ***ptr_divu_for_2dcell, ***ptr_qvis_for_2dcell;
    double dx[3];

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        cs_for_2dcell  = h5_read_2d(mesh_gid, "cs", ncell_ext);
        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);
        vel_x = h5_read_2d(mesh_gid, "vel_x", nnode_ext);
        vel_y = h5_read_2d(mesh_gid, "vel_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
        copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);

        // Call the function under test (compute_qvis depends on output from compute_divu)
        compute_divu(dim, ncell, nbdry, dx, vel_for_2dnode, NULL);
        compute_qvis(dim, ncell, nbdry, dx[0], rho_for_2dcell, cs_for_2dcell, vel_for_2dnode, NULL, NULL, NULL);

        // Get outputs for divu and qvis
        ptr_divu_for_2dcell = get_ptr_divu();
        ptr_qvis_for_2dcell = get_ptr_qvis();

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        copy_components_from_2dvel (nnode_ext, vel_for_2dnode, vel_x, vel_y);
        h5_write_2d(test_gid, "vel_x", vel_x, nnode_ext);
        h5_write_2d(test_gid, "vel_y", vel_y, nnode_ext);
        h5_write_2d(test_gid, "divu", *ptr_divu_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "qvis", *ptr_qvis_for_2dcell, ncell_ext);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);
        free(vel_x[0]); free(vel_x); free(vel_y[0]); free(vel_y);
        free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
        free(cs_for_2dcell[0]); free(cs_for_2dcell);
        free(rho_for_2dcell[0]); free(rho_for_2dcell);
        free((*ptr_divu_for_2dcell)[0]); free(*ptr_divu_for_2dcell);
        *ptr_divu_for_2dcell = NULL;
        free((*ptr_qvis_for_2dcell)[0]); free(*ptr_qvis_for_2dcell);
        *ptr_qvis_for_2dcell = NULL;
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__compute_qvis

void record_test__compute_force(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3], direction;
    double **vel_x, **vel_y, ***vel_for_2dnode, **rho_for_2dcell, **cs_for_2dcell, **pres_for_2dcell;
    double ***ptr_divu_for_2dcell, ***ptr_qvis_for_2dcell, ***ptr_force_for_2dnode;
    double dx[3];

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        cs_for_2dcell  = h5_read_2d(mesh_gid, "cs", ncell_ext);
        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);
        pres_for_2dcell  = h5_read_2d(mesh_gid, "pres", ncell_ext);
        vel_x = h5_read_2d(mesh_gid, "vel_x", nnode_ext);
        vel_y = h5_read_2d(mesh_gid, "vel_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
        copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);
        direction = random_int(0,1); // set direction to be either 0 or 1

        // Call the function under test (compute_force depends on both compute_divu and compute_qvis)
        compute_divu(dim, ncell, nbdry, dx, vel_for_2dnode, NULL);
        compute_qvis(dim, ncell, nbdry, dx[0], rho_for_2dcell, cs_for_2dcell, vel_for_2dnode, NULL, NULL, NULL);
        compute_force(dim, ncell, nbdry, dx, direction, pres_for_2dcell, NULL, rho_for_2dcell, NULL);

        // Get outputs for divu, qvis, force
        ptr_divu_for_2dcell = get_ptr_divu();
        ptr_qvis_for_2dcell = get_ptr_qvis();
        ptr_force_for_2dnode = get_ptr_force();

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        copy_components_from_2dvel (nnode_ext, vel_for_2dnode, vel_x, vel_y);
        h5_write_2d(test_gid, "vel_x", vel_x, nnode_ext);
        h5_write_2d(test_gid, "vel_y", vel_y, nnode_ext);
        h5_write_2d(test_gid, "divu", *ptr_divu_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "qvis", *ptr_qvis_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "force", *ptr_force_for_2dnode, nnode_ext);
        h5_write_1d_int(test_gid, "direction", &direction, 1);

        H5Gclose(test_gid);

        H5Gclose(mesh_gid);
        free(vel_x[0]); free(vel_x); free(vel_y[0]); free(vel_y);
        free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
        free(cs_for_2dcell[0]); free(cs_for_2dcell);
        free(rho_for_2dcell[0]); free(rho_for_2dcell);
        free(pres_for_2dcell[0]); free(pres_for_2dcell);
        free((*ptr_divu_for_2dcell)[0]); free(*ptr_divu_for_2dcell);
        *ptr_divu_for_2dcell = NULL;
        free((*ptr_qvis_for_2dcell)[0]); free(*ptr_qvis_for_2dcell);
        *ptr_qvis_for_2dcell = NULL;
        free((*ptr_force_for_2dnode)[0]); free(*ptr_force_for_2dnode);
        *ptr_force_for_2dnode = NULL;
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__compute_force

void record_test__update_vel_comp(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3], direction[1];
    double **vel_x, **vel_y, **vav_x, **vav_y, ***vel_for_2dnode, ***velav_for_2dnode;
    double**rho_for_2dcell, **cs_for_2dcell, **pres_for_2dcell;
    double ***ptr_divu_for_2dcell, ***ptr_qvis_for_2dcell, ***ptr_force_for_2dnode;
    double dx[3], dt[1];

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        vel_x = h5_read_2d(mesh_gid, "vel_x", nnode_ext);
        vel_y = h5_read_2d(mesh_gid, "vel_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
        copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);

        vav_x = h5_read_2d(mesh_gid, "vav_x", nnode_ext);
        vav_y = h5_read_2d(mesh_gid, "vav_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &velav_for_2dnode);
        copy_2dvel_components (nnode_ext, vav_x, vav_y, velav_for_2dnode);

        direction[0] = random_int(0,1); // set direction to be either 0 or 1
        dt[0]  = random_double(1.0e-9,1.0); // random double between 1e-9 and 1

        cs_for_2dcell  = h5_read_2d(mesh_gid, "cs", ncell_ext);
        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);
        pres_for_2dcell  = h5_read_2d(mesh_gid, "pres", ncell_ext);

        // Call the function under test (depends on compute_force, which depends on compute_divu and compute_qvis)
        compute_divu(dim, ncell, nbdry, dx, vel_for_2dnode, NULL);
        compute_qvis(dim, ncell, nbdry, dx[0], rho_for_2dcell, cs_for_2dcell, vel_for_2dnode, NULL, NULL, NULL);
        compute_force(dim, ncell, nbdry, dx, direction[0], pres_for_2dcell, NULL, rho_for_2dcell, NULL);
        update_vel_comp(dim, ncell, nbdry, dt[0], direction[0], vel_for_2dnode, velav_for_2dnode, NULL, NULL);

        // Get outputs for divu, qvis, force
        ptr_divu_for_2dcell = get_ptr_divu();
        ptr_qvis_for_2dcell = get_ptr_qvis();
        ptr_force_for_2dnode = get_ptr_force();


        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        copy_components_from_2dvel (nnode_ext, vel_for_2dnode, vel_x, vel_y);
        copy_components_from_2dvel (nnode_ext, velav_for_2dnode, vav_x, vav_y);
        h5_write_2d(test_gid, "vel_x", vel_x, nnode_ext);
        h5_write_2d(test_gid, "vel_y", vel_y, nnode_ext);
        h5_write_2d(test_gid, "vav_x", vav_x, nnode_ext);
        h5_write_2d(test_gid, "vav_y", vav_y, nnode_ext);
        h5_write_2d(test_gid, "divu", *ptr_divu_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "qvis", *ptr_qvis_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "force", *ptr_force_for_2dnode, nnode_ext);
        h5_write_1d_int(test_gid, "direction", direction, 1);
        h5_write_1d(test_gid, "dt", dt, 1);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);
        free(vel_x[0]); free(vel_x); free(vel_y[0]); free(vel_y);
        free(vav_x[0]); free(vav_x); free(vav_y[0]); free(vav_y);
        free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
        free(velav_for_2dnode[0][0]); free(velav_for_2dnode[0]); free(velav_for_2dnode);
        free(cs_for_2dcell[0]); free(cs_for_2dcell);
        free(rho_for_2dcell[0]); free(rho_for_2dcell);
        free(pres_for_2dcell[0]); free(pres_for_2dcell);
        free((*ptr_divu_for_2dcell)[0]); free(*ptr_divu_for_2dcell);
        *ptr_divu_for_2dcell = NULL;
        free((*ptr_qvis_for_2dcell)[0]); free(*ptr_qvis_for_2dcell);
        *ptr_qvis_for_2dcell = NULL;
        free((*ptr_force_for_2dnode)[0]); free(*ptr_force_for_2dnode);
        *ptr_force_for_2dnode = NULL;
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__update_vel_comp

void record_test__update_density(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    int ijksize[2], mxsize[2];
    int nmixcell, **ijk_in_mixcell, *nmat_in_mixcell , **nmat_for_2dcell;
    double **vel_x, **vel_y, ***vel_for_2dnode;
    double **rho_for_2dcell;
    double ***ptr_divu_for_2dcell, **vf_in_mixcell, **rho_in_mixcell, cournt;
    double dx[3], dt;

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;
    int nm, mx, m;
    int lsize;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        vel_x = h5_read_2d(mesh_gid, "vel_x", nnode_ext);
        vel_y = h5_read_2d(mesh_gid, "vel_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
        copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);
        nmat_for_2dcell = h5_read_2d_int(mesh_gid, "nmat", ncell_ext);
        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);

        dt  = random_double(1.0e-9,1.0); // random double between 1e-9 and 1

        // Count nmixcell based on how many cells in nmat_for_2dcell have more than 1 material
        nmixcell = 0;
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {  // Ignore boundary cells
                if (nmat_for_2dcell[j][i] > 1) {
                    nmixcell++;
                }
            }
        }

        if (nmixcell > 0) {

            // 1. allocate and fill nmat_in_mixcell and ijk_in_mixcell
            nmat_in_mixcell = (int *)malloc(nmixcell * sizeof(int));
            ijk_in_mixcell  = (int **)malloc(nmixcell * sizeof(int*));
            ijk_in_mixcell[0] = (int *)malloc(dim * nmixcell * sizeof(int));
            for (mx = 1; mx < nmixcell; mx++)
                ijk_in_mixcell[mx] = ijk_in_mixcell[mx-1] + dim;

            mx = 0;
            for (j = 0; j < ncell_ext[1]; j++) {
                for (i = 0; i < ncell_ext[0]; i++) {
                    if (nmat_for_2dcell[j][i] > 1) {
                        nmat_in_mixcell[mx] = nmat_for_2dcell[j][i];
                        ijk_in_mixcell[mx][0] = i;
                        ijk_in_mixcell[mx][1] = j;
                        ++mx;
                    }
                }
            }

            // 2. count the size of arrays for vf_in_mixcell and rho_in_mixcell
            lsize = 0;
            for (mx = 0; mx < nmixcell; mx++)
                lsize += nmat_in_mixcell[mx];

            // 3. allocate vf_in_mixcell and rho_in_mixcell
            vf_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            rho_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            vf_in_mixcell[0] = (double *)malloc(lsize * sizeof(double));
            rho_in_mixcell[0] = (double *)malloc(lsize * sizeof(double));

            int offset = nmat_in_mixcell[0];
            for (mx = 1; mx < nmixcell; mx++) {
                vf_in_mixcell[mx]  = vf_in_mixcell[0] + offset;
                rho_in_mixcell[mx] = rho_in_mixcell[0] + offset;
                offset += nmat_in_mixcell[mx];
            }

            // 4. fill out the vf and rho arrays in mixcell
            for (mx = 0; mx < nmixcell; mx++) {
                // Initialize vf_in_mixcell (random doubles between 0 and 1, sum = 1.0 for each cell)
                nm = nmat_in_mixcell[mx];
                double sum_vf = 0.;
                for (m = 0; m < nm; m++) {
                    rho_in_mixcell[mx][m] = random_double(0.0, 10.0);
                    vf_in_mixcell[mx][m] = random_double(0.0, 1.0);
                    sum_vf += vf_in_mixcell[mx][m];
                }
                // Normalize so that the sum equals 1.0
                for (m = 0; m < nm; m++) {
                    vf_in_mixcell[mx][m] /= sum_vf;
                }
            }

        }
        else {
            // If no mixed cells, set pointers to NULL
            ijk_in_mixcell = NULL;
            nmat_in_mixcell = NULL;
            vf_in_mixcell = NULL;
            rho_in_mixcell = NULL;
        }

        // Store rho_in_mixcell input before updating density 
        test_gid = h5_open_group(testfile_id, group_name);
        if (nmixcell > 0) {
            h5_write_1d(test_gid, "rho_mixcell_old", rho_in_mixcell[0], lsize);
        }

        // Call the function under test (depends compute_divu)
        compute_divu(dim, ncell, nbdry, dx, vel_for_2dnode, NULL);
        update_density(dim, ncell, nbdry, dt, nmat_for_2dcell, rho_for_2dcell, NULL, NULL,
                       nmixcell, ijk_in_mixcell, nmat_in_mixcell, vf_in_mixcell, rho_in_mixcell, &cournt);
        // Get outputs for divu
        ptr_divu_for_2dcell = get_ptr_divu();

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        copy_components_from_2dvel (nnode_ext, vel_for_2dnode, vel_x, vel_y);
        h5_write_2d(test_gid, "vel_x", vel_x, nnode_ext);
        h5_write_2d(test_gid, "vel_y", vel_y, nnode_ext);
        h5_write_2d(test_gid, "divu", *ptr_divu_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "rho", rho_for_2dcell, ncell_ext);
        h5_write_1d(test_gid, "dt", &dt, 1);
        if (nmixcell > 0) {
            h5_write_1d_int(test_gid, "ijk_mixcell", ijk_in_mixcell[0], dim*nmixcell);
            h5_write_1d_int(test_gid, "nmat_mixcell", nmat_in_mixcell, nmixcell);
            h5_write_1d(test_gid, "vf_mixcell", vf_in_mixcell[0], lsize);
            h5_write_1d(test_gid, "rho_mixcell_new", rho_in_mixcell[0], lsize);
        }
        h5_write_1d(test_gid, "courant", &cournt, 1);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);

        free(vel_x[0]); free(vel_x); free(vel_y[0]); free(vel_y);
        free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
        free(rho_for_2dcell[0]); free(rho_for_2dcell);
        free((*ptr_divu_for_2dcell)[0]); free(*ptr_divu_for_2dcell);
        *ptr_divu_for_2dcell = NULL;
        if (nmixcell > 0) {
            // Free the outer arrays
            free(ijk_in_mixcell[0]); free(ijk_in_mixcell);
            free(vf_in_mixcell[0]);  free(vf_in_mixcell);
            free(rho_in_mixcell[0]); free(rho_in_mixcell);
            free(nmat_in_mixcell);
        }
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__update_density
void record_test__update_energy(const char *meshes_h5, const char *test_h5, const int verbose){

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    int nmixcell, **ijk_in_mixcell, *nmat_in_mixcell , **nmat_for_2dcell;
    double **vel_x, **vel_y, ***vel_for_2dnode;
    double **rho_for_2dcell_old, **rho_for_2dcell;
    double ***ptr_divu_for_2dcell, **vf_in_mixcell, **rho_in_mixcell_old, **rho_in_mixcell;
    double **pres_for_2dcell, **es_for_2dcell_old, **ei_for_2dcell;
    double **pres_in_mixcell, **es_in_mixcell_old, **ei_in_mixcell;

    double dx[3], dt, cournt;

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;
    int nm, mx, m, l;
    long long lsize; 
    int nm_tot;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        vel_x = h5_read_2d(mesh_gid, "vel_x", nnode_ext);
        vel_y = h5_read_2d(mesh_gid, "vel_y", nnode_ext);
        allocate_velocity_2dmesh (nnode_ext, &vel_for_2dnode);
        copy_2dvel_components (nnode_ext, vel_x, vel_y, vel_for_2dnode);
        nmat_for_2dcell = h5_read_2d_int(mesh_gid, "nmat", ncell_ext);
        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);
        ei_for_2dcell   = h5_read_2d(mesh_gid, "ei",  ncell_ext);
        pres_for_2dcell = h5_read_2d(mesh_gid, "pres",ncell_ext);
        dt  = random_double(1.0e-9,1.0); // random double between 1e-9 and 1

        // Count nmixcell based on how many cells in nmat_for_2dcell have more than 1 material
        nmixcell = 0;
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {  // Ignore boundary cells
                if (nmat_for_2dcell[j][i] > 1) {
                    nmixcell++;
                }
            }
        }

        if (nmixcell > 0) {

            // 1. allocate and fill nmat_in_mixcell and ijk_in_mixcell
            nmat_in_mixcell = (int *)malloc(nmixcell * sizeof(int));
            ijk_in_mixcell  = (int **)malloc(nmixcell * sizeof(int*));
            ijk_in_mixcell[0] = (int *)malloc(dim * nmixcell * sizeof(int));
            for (mx = 1; mx < nmixcell; mx++)
                ijk_in_mixcell[mx] = ijk_in_mixcell[mx-1] + dim;

            mx = 0;
            for (j = 0; j < ncell_ext[1]; j++) {
                for (i = 0; i < ncell_ext[0]; i++) {
                    if (nmat_for_2dcell[j][i] > 1) {
                        nmat_in_mixcell[mx] = nmat_for_2dcell[j][i];
                        ijk_in_mixcell[mx][0] = i;
                        ijk_in_mixcell[mx][1] = j;
                        ++mx;
                    }
                }
            }

            // 2. count the size of arrays for vf_in_mixcell and rho_in_mixcell
            nm_tot = 0;
            for (mx = 0; mx < nmixcell; mx++)
                nm_tot += nmat_in_mixcell[mx];

            // 3. allocate vf_in_mixcell, rho_in_mixcell, pres_in_mixcell
            vf_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            rho_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            pres_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            ei_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));

            vf_in_mixcell[0] = (double *)malloc(nm_tot * sizeof(double));
            rho_in_mixcell[0] = (double *)malloc(nm_tot * sizeof(double));
            pres_in_mixcell[0] = (double *)malloc(nm_tot * sizeof(double));
            ei_in_mixcell[0] = (double *)malloc(nm_tot * sizeof(double));

            int offset = nmat_in_mixcell[0];
            for (mx = 1; mx < nmixcell; mx++) {
                vf_in_mixcell[mx]  = vf_in_mixcell[0] + offset;
                rho_in_mixcell[mx] = rho_in_mixcell[0] + offset;
                pres_in_mixcell[mx] = pres_in_mixcell[0] + offset;
                ei_in_mixcell[mx] = ei_in_mixcell[0] + offset;
                offset += nmat_in_mixcell[mx];
            }

            // 4. fill out the vf, rho, pres arrays in mixcell
            for (mx = 0; mx < nmixcell; mx++) {
                // Initialize vf_in_mixcell (random doubles between 0 and 1, sum = 1.0 for each cell)
                nm = nmat_in_mixcell[mx];
                double sum_vf = 0.;
                for (m = 0; m < nm; m++) {
                    rho_in_mixcell[mx][m] = random_double(0.0, 10.0);
                    pres_in_mixcell[mx][m] = random_double(0.0, 10.0);
                    vf_in_mixcell[mx][m] = random_double(0.0, 1.0);
                    ei_in_mixcell[mx][m] = random_double(0.0, 10.0);
                    sum_vf += vf_in_mixcell[mx][m];
                }
                // Normalize so that the sum equals 1.0
                for (m = 0; m < nm; m++) {
                    vf_in_mixcell[mx][m] /= sum_vf;
                }
            }

            // 5. Allocate memory for rho and es in mixcell old variables 
            rho_in_mixcell_old = (double **) malloc((nmixcell) * sizeof(double *));
            es_in_mixcell_old  = (double **) malloc((nmixcell) * sizeof(double *));
            rho_in_mixcell_old[0] = (double *) malloc((nm_tot) * sizeof(double));
            es_in_mixcell_old[0]  = (double *) malloc((nm_tot) * sizeof(double));    
            for (mx = 1; mx < nmixcell; mx++) {
                rho_in_mixcell_old[mx] = rho_in_mixcell_old[mx-1] + nmat_in_mixcell[mx-1];
                es_in_mixcell_old[mx]  = es_in_mixcell_old[mx-1]  + nmat_in_mixcell[mx-1];
            }

            // Copy variables, compute es_in_mixcell
            memcpy(rho_in_mixcell_old[0], rho_in_mixcell[0], (size_t)(nm_tot * sizeof(double)));
            for (m = 0; m < nm_tot; m++) {
                es_in_mixcell_old[0][m] = ei_in_mixcell[0][m]/(rho_in_mixcell[0][m] + tiny);
            }
        }
        else {
            // If no mixed cells, set pointers to NULL
            ijk_in_mixcell = NULL;
            nmat_in_mixcell = NULL;
            vf_in_mixcell = NULL;
            rho_in_mixcell = NULL;
            pres_in_mixcell = NULL;
            rho_in_mixcell_old = NULL; 
            es_in_mixcell_old = NULL;  
            ei_in_mixcell = NULL;     
        }

        // Set rho_for_2dcell_old and es_for_2dcell_old
        lsize = 1.0;
        for (i = 0; i < dim; i++) {
            lsize *= ncell_ext[i];
        }
        rho_for_2dcell_old    = (double **) malloc(ncell_ext[1] * sizeof(double *));
        rho_for_2dcell_old[0] = (double  *) malloc(lsize * sizeof(double));
        es_for_2dcell_old     = (double **) malloc(ncell_ext[1] * sizeof(double *));
        es_for_2dcell_old[0]  = (double  *) malloc(lsize * sizeof(double));
        for (j = 1; j < ncell_ext[1]; j++) {
            rho_for_2dcell_old[j]  = rho_for_2dcell_old[j-1] + ncell_ext[0];
            es_for_2dcell_old[j]   = es_for_2dcell_old[j-1]  + ncell_ext[0];
        }
        if (dim == 2) {
            for (l = 0; l < lsize; l++) {
                rho_for_2dcell_old[0][l] = rho_for_2dcell[0][l];
                es_for_2dcell_old[0][l]  = ei_for_2dcell[0][l]/(rho_for_2dcell[0][l] + tiny);
            }
        }
       
        // Store ei_in_mixcell input before updating density 
        test_gid = h5_open_group(testfile_id, group_name);
        if (nmixcell > 0) {
            h5_write_1d(test_gid, "ei_mixcell_old", ei_in_mixcell[0], nm_tot);
        }

        // Call the function under test (depends compute_divu and update_density)
        compute_divu(dim, ncell, nbdry, dx, vel_for_2dnode, NULL);
        update_density(dim, ncell, nbdry, dt, nmat_for_2dcell, rho_for_2dcell, NULL, NULL,
                nmixcell, ijk_in_mixcell, nmat_in_mixcell, vf_in_mixcell, rho_in_mixcell, &cournt);

        update_energy(dim, ncell, nbdry, dt, nmat_for_2dcell, rho_for_2dcell_old, rho_for_2dcell,
		    pres_for_2dcell, es_for_2dcell_old, ei_for_2dcell,
            NULL, NULL,  NULL,NULL, NULL, NULL,
		    nmixcell, ijk_in_mixcell, nmat_in_mixcell, vf_in_mixcell,
		    rho_in_mixcell_old, rho_in_mixcell,
		    pres_in_mixcell, es_in_mixcell_old, ei_in_mixcell);

        // Get outputs for divu
        ptr_divu_for_2dcell = get_ptr_divu();

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        copy_components_from_2dvel (nnode_ext, vel_for_2dnode, vel_x, vel_y);
        h5_write_2d(test_gid, "vel_x", vel_x, nnode_ext);
        h5_write_2d(test_gid, "vel_y", vel_y, nnode_ext);
        h5_write_2d(test_gid, "divu", *ptr_divu_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "rho", rho_for_2dcell, ncell_ext);
        h5_write_2d(test_gid, "es", es_for_2dcell_old, ncell_ext);
        h5_write_2d(test_gid, "ei", ei_for_2dcell, ncell_ext);

        h5_write_1d(test_gid, "dt", &dt, 1);
        if (nmixcell > 0) {
            h5_write_1d_int(test_gid, "ijk_mixcell", ijk_in_mixcell[0], dim*nmixcell);
            h5_write_1d_int(test_gid, "nmat_mixcell", nmat_in_mixcell, nmixcell);
            h5_write_1d(test_gid, "vf_mixcell", vf_in_mixcell[0], nm_tot);
            h5_write_1d(test_gid, "rho_mixcell_old", rho_in_mixcell_old[0], nm_tot);
            h5_write_1d(test_gid, "rho_mixcell_new", rho_in_mixcell[0], nm_tot);
            h5_write_1d(test_gid, "pres_mixcell", pres_in_mixcell[0], nm_tot);
            h5_write_1d(test_gid, "ei_mixcell", ei_in_mixcell[0], nm_tot);
            h5_write_1d(test_gid, "es_mixcell_old", es_in_mixcell_old[0], nm_tot);
        }
        h5_write_1d(test_gid, "courant", &cournt, 1);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);

        free(vel_x[0]); free(vel_x); free(vel_y[0]); free(vel_y);
        free(vel_for_2dnode[0][0]); free(vel_for_2dnode[0]); free(vel_for_2dnode);
        free(rho_for_2dcell[0]); free(rho_for_2dcell);
        free(rho_for_2dcell_old[0]); free(rho_for_2dcell_old);
        free(es_for_2dcell_old[0]); free(es_for_2dcell_old);
        free(pres_for_2dcell[0]); free(pres_for_2dcell);
        free(ei_for_2dcell[0]); free(ei_for_2dcell);
        free(nmat_for_2dcell[0]); free(nmat_for_2dcell);
        free((*ptr_divu_for_2dcell)[0]); free(*ptr_divu_for_2dcell);
        *ptr_divu_for_2dcell = NULL;
        if (nmixcell > 0) {
            // Free the outer arrays
            free(ijk_in_mixcell[0]); free(ijk_in_mixcell);
            free(vf_in_mixcell[0]);  free(vf_in_mixcell);
            free(rho_in_mixcell[0]); free(rho_in_mixcell);
            free(rho_in_mixcell_old[0]); free(rho_in_mixcell_old);
            free(nmat_in_mixcell);
            free(pres_in_mixcell[0]); free(pres_in_mixcell);
            free(ei_in_mixcell[0]); free(ei_in_mixcell);
            free(es_in_mixcell_old[0]); free(es_in_mixcell_old);
        }
    }
    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);
} // record_test__update_energy


void record_test__mesh_sound_speed(const char *meshes_h5, const char *test_h5, const int verbose){
    int dim, ncell[3], nbdry, ncell_ext[3];
    int  nmat, *matids, *is_solid, **nmat_for_2dcell, **matid_for_2dcell;
    double **cs_for_2dcell, *gamma_ea_mat;
    double **rho_for_2dcell, **pres_for_2dcell;


    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5); 
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
        }

        cs_for_2dcell  = h5_read_2d(mesh_gid, "cs", ncell_ext);
        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);
        pres_for_2dcell  = h5_read_2d(mesh_gid, "pres", ncell_ext);
        nmat_for_2dcell  = h5_read_2d_int(mesh_gid, "nmat", ncell_ext);
        matid_for_2dcell = h5_read_2d_int(mesh_gid, "matid", ncell_ext);

        // Find the maximum number of materials allocated 
        nmat = 1; 
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {
                nmat = MAX(nmat, nmat_for_2dcell[j][i]);
            }
        }
        matids = (int *)malloc(nmat * sizeof(int));
        is_solid = (int *)malloc(nmat * sizeof(int));
        gamma_ea_mat = (double *)malloc(nmat * sizeof(double));

        for (int m=0; m < nmat; m++) {
            gamma_ea_mat[m] = random_double(1.1,1.9);
            matids[m] = m;
            is_solid[m] = random_int(0,1);
        }

         // Call function under test
        mesh_put_mesh_data(dim, ncell, nbdry, nmat, matids,
                         nmat_for_2dcell, matid_for_2dcell, rho_for_2dcell,
                         pres_for_2dcell, NULL,
                         NULL, NULL,
                         NULL, NULL, NULL,
                         NULL, NULL,
                         NULL, NULL);
        mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat, cs_for_2dcell, NULL);

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_2d(test_gid, "cs", cs_for_2dcell, ncell_ext);
        h5_write_1d_int(test_gid, "nmat", &nmat, 1);
        h5_write_1d_int(test_gid, "matids", matids, nmat); 
        h5_write_1d_int(test_gid, "issolid", is_solid, nmat); 
        h5_write_1d(test_gid, "gamma", gamma_ea_mat, nmat);

        H5Gclose(test_gid);

        H5Gclose(mesh_gid);

        free(cs_for_2dcell[0]); free(cs_for_2dcell);
        free(rho_for_2dcell[0]); free(rho_for_2dcell);
        free(pres_for_2dcell[0]); free(pres_for_2dcell);
        free(nmat_for_2dcell[0]); free(nmat_for_2dcell);
        free(matid_for_2dcell[0]); free(matid_for_2dcell);
        free(is_solid); 
        free(matids); 
        free(gamma_ea_mat);

        mesh_put_mesh_data(0, NULL, 0, 0, NULL,
                         NULL, NULL, NULL, NULL, NULL,
                         NULL, NULL,
                         NULL, NULL, NULL,
                         NULL, NULL,
                         NULL, NULL);
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} //record_test__mesh_sound_speed

void record_test__mat_sound_speed(const char *meshes_h5, const char *test_h5, const int verbose){
    int dim, ncell[3], nbdry, ncell_ext[3];
    int  nmat, *matids, *is_solid, **nmat_for_2dcell;
    double **cs_for_2dcell, *gamma_ea_mat;
    double **rho_in_mixcell, **pres_in_mixcell, **vf_in_mixcell; 
    int **ijk_in_mixcell, *nmat_in_mixcell, nmixcell, nm_tot, **matids_in_mixcell;

    char group_name[50];
    hid_t mesh_gid, test_gid;
    int num_meshes, imesh, d, i, j;
    int mx, nm, m;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5); 
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        if (dim != 2){
            printf("%s: dimension not implemented\n",group_name);
            H5Gclose(mesh_gid);
            continue;
        }
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
        }

        cs_for_2dcell  = h5_read_2d(mesh_gid, "cs", ncell_ext);
        nmat_for_2dcell  = h5_read_2d_int(mesh_gid, "nmat", ncell_ext);

        // Find the maximum number of materials allocated 
        nmat = 1; 
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {
                nmat = MAX(nmat, nmat_for_2dcell[j][i]);
            }
        }

        matids = (int *)malloc(nmat * sizeof(int));
        is_solid = (int *)malloc(nmat * sizeof(int));
        gamma_ea_mat = (double *)malloc(nmat * sizeof(double));

        for (int m=0; m < nmat; m++) {
            gamma_ea_mat[m] = 1.0;
            matids[m] = m;
            is_solid[m] = random_int(0,1);
        }

        // Count nmixcell based on how many cells in nmat_for_2dcell have more than 1 material
        nmixcell = 0;
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {  // Ignore boundary cells
                if (nmat_for_2dcell[j][i] > 1) {
                    nmixcell++;
                }
            }
        }

        if (nmixcell > 0) {

            // 1. allocate and fill nmat_in_mixcell and ijk_in_mixcell
            nmat_in_mixcell = (int *)malloc(nmixcell * sizeof(int));
            ijk_in_mixcell  = (int **)malloc(nmixcell * sizeof(int*));
            ijk_in_mixcell[0] = (int *)malloc(dim * nmixcell * sizeof(int));
            for (mx = 1; mx < nmixcell; mx++)
                ijk_in_mixcell[mx] = ijk_in_mixcell[mx-1] + dim;

            mx = 0;
            for (j = 0; j < ncell_ext[1]; j++) {
                for (i = 0; i < ncell_ext[0]; i++) {
                    if (nmat_for_2dcell[j][i] > 1) {
                        nmat_in_mixcell[mx] = nmat_for_2dcell[j][i];
                        ijk_in_mixcell[mx][0] = i;
                        ijk_in_mixcell[mx][1] = j;
                        ++mx;
                    }
                }
            }

            // 2. count the size of arrays for vf_in_mixcell and rho_in_mixcell
            nm_tot = 0;
            for (mx = 0; mx < nmixcell; mx++)
                nm_tot += nmat_in_mixcell[mx];

            // 3. allocate vf_in_mixcell, rho_in_mixcell, pres_in_mixcell, matids_in_mixcell
            vf_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            rho_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            pres_in_mixcell = (double **)malloc(nmixcell * sizeof(double *));
            matids_in_mixcell = (int **)malloc(nmixcell * sizeof(int *));

            vf_in_mixcell[0] = (double *)malloc(nm_tot * sizeof(double));
            rho_in_mixcell[0] = (double *)malloc(nm_tot * sizeof(double));
            pres_in_mixcell[0] = (double *)malloc(nm_tot * sizeof(double));
            matids_in_mixcell[0] = (int *)malloc(nm_tot * sizeof(int));

            int offset = nmat_in_mixcell[0];
            for (mx = 1; mx < nmixcell; mx++) {
                vf_in_mixcell[mx]  = vf_in_mixcell[0] + offset;
                rho_in_mixcell[mx] = rho_in_mixcell[0] + offset;
                pres_in_mixcell[mx] = pres_in_mixcell[0] + offset;
                matids_in_mixcell[mx] = matids_in_mixcell[0] + offset;
                offset += nmat_in_mixcell[mx];
            }

            // 4. fill out the vf, rho, pres arrays in mixcell
            for (mx = 0; mx < nmixcell; mx++) {
                // Initialize vf_in_mixcell (random doubles between 0 and 1, sum = 1.0 for each cell)
                nm = nmat_in_mixcell[mx];
                double sum_vf = 0.;
                for (m = 0; m < nm; m++) {
                    rho_in_mixcell[mx][m] = random_double(0.0, 10.0);
                    pres_in_mixcell[mx][m] = random_double(0.0, 10.0);
                    vf_in_mixcell[mx][m] = random_double(0.0, 1.0);
                    matids_in_mixcell[mx][m] = random_int(0,nmat-1);
                    sum_vf += vf_in_mixcell[mx][m];
                }
                // Normalize so that the sum equals 1.0
                for (m = 0; m < nm; m++) {
                    vf_in_mixcell[mx][m] /= sum_vf;
                }
            }
    
        }
        else {
            // If no mixed cells, set pointers to NULL
            ijk_in_mixcell = NULL;
            nmat_in_mixcell = NULL;
            vf_in_mixcell = NULL;
            rho_in_mixcell = NULL;
            pres_in_mixcell = NULL;

        }

        // Call function under test
        // Set nmat_prob, matids_prob, nmixcell, nmat_in_mixcell, matids_in_mixcell
        // pres_in_mixcell, rho_in_mixcell, vf_in_mixcell, ijk_in_mixcell
        mat_put_mix_data(nmat, matids,
                         nmixcell, 0, 0,
                         nmat_in_mixcell, ijk_in_mixcell, matids_in_mixcell, 
                         vf_in_mixcell, rho_in_mixcell,
                         pres_in_mixcell, NULL);
        mat_sound_speed(dim, ncell, nbdry, nmat, matids, is_solid, gamma_ea_mat, cs_for_2dcell, NULL);

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_2d(test_gid, "cs", cs_for_2dcell, ncell_ext);
        h5_write_1d_int(test_gid, "nmat", &nmat, 1);
        h5_write_1d_int(test_gid, "matids", matids, nmat); 
        h5_write_1d_int(test_gid, "issolid", is_solid, nmat); 
        h5_write_1d(test_gid, "gamma", gamma_ea_mat, nmat);
        if (nmixcell > 0) {
            h5_write_1d_int(test_gid, "ijk_mixcell", ijk_in_mixcell[0], dim*nmixcell);
            h5_write_1d_int(test_gid, "nmat_mixcell", nmat_in_mixcell, nmixcell);
            h5_write_1d_int(test_gid, "matids_mixcell", matids_in_mixcell[0], nm_tot);
            h5_write_1d(test_gid, "vf_mixcell", vf_in_mixcell[0], nm_tot);
            h5_write_1d(test_gid, "rho_mixcell", rho_in_mixcell[0], nm_tot);
            h5_write_1d(test_gid, "pres_mixcell", pres_in_mixcell[0], nm_tot);
        }

        H5Gclose(test_gid);

        H5Gclose(mesh_gid);

        free(cs_for_2dcell[0]); free(cs_for_2dcell);
        free(is_solid); 
        free(matids); 
        free(gamma_ea_mat);
        if (nmixcell > 0) {
            free(ijk_in_mixcell[0]); free(ijk_in_mixcell);
            free(nmat_in_mixcell);
            free(matids_in_mixcell[0]); free(matids_in_mixcell);
            free(vf_in_mixcell[0]); free(vf_in_mixcell);
            free(rho_in_mixcell[0]); free(rho_in_mixcell); 
            free(pres_in_mixcell[0]); free(pres_in_mixcell);
        }

        mat_put_mix_data(0, NULL,
                         0, 0, 0,
                         NULL, NULL, NULL, 
                         NULL, NULL,
                         NULL, NULL);

    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} //record_test__mat_sound_speed


// Function to generate random polygon with the correct node ordering
//  - outputs node coordinates in the "coords" array;
//  - assumes that coords has already been allocatex.
//  - TODO: make sure to generate only convex polygons
//
void generate_random_convex_2dpoly (int nnode, double xrange, double yrange, double *coords) {
    double angle[100], xch, ax, ay, bx, by;
    int i, j, ix, iy, jx, jy;
    bool convex = false;

    while (!convex) {
    for (i = 0; i < nnode; i++) {
        ix = 2*i;
        iy = 2*i + 1;
        coords[ix] = random_double(-xrange/2, xrange/2);
        coords[iy] = random_double(-yrange/2, yrange/2);
        angle[i] = atan2(coords[iy], coords[ix]);
    }

    // insertion sort
    for (i = 0; i < nnode; i++) {
        ix = 2*i;
        iy = 2*i + 1;
        for (j = i + 1; j < nnode; j++) {
            if (angle[j] < angle[i]) {
                jx = 2*j;
                jy = 2*j + 1;
                xch = angle[i]; angle[i] = angle[j]; angle[j] = xch;
                xch = coords[ix]; coords[ix] = coords[jx]; coords[jx] = xch;
                xch = coords[iy]; coords[iy] = coords[jy]; coords[jy] = xch;
            }
        }
    }

    // check convexity
    convex = true;
    for (i = 0; i < nnode; i++) {
        ix = 2*i;
        iy = 2*i + 1;

        jx = (ix + 2) % (2*nnode);
        jy = (iy + 2) % (2*nnode);
        ax = coords[jx] - coords[ix];
        ay = coords[jy] - coords[iy];

        ix = (jx + 2) % (2*nnode);
        iy = (jy + 2) % (2*nnode);
        bx = coords[ix] - coords[jx];
        by = coords[iy] - coords[jy];

        convex = convex && (ax*by - ay*bx > 0);
    }} // convex

    // printf("# nnode = %d, convex: %d\n", nnode, convex);
    // for (i = 0; i < nnode; i++) {
    //     ix = 2*i;
    //     iy = 2*i + 1;
    //     printf ("%7.4f  %7.4f\n", coords[ix], coords[iy]);
    // }
    // printf("\n\n");

} // generate_random_convex_2dpoly


// Function to record test cases for find_interface2d into an HDF5 file
void record_test__find_interface2d(int num_tests, const char *test_h5, const int verbose) {

    // HDF5 group handles
    const int nnode_max = 6;
    const double range = 5.0;
    double coords[2*nnode_max];
    int nodelist[nnode_max], node_loc[nnode_max];

    // outputs
    int nnode_new, nnode_interface, nnode_lower, nnode_upper;
    double coords_new[4*nnode_max];
    int nodes_interface[2*nnode_max];

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    for (int t = 0; t < num_tests; t++) {
        // Randomly generate input parameters
        int nnode = random_int(3, 6);  // Number of nodes in the polygon
        generate_random_convex_2dpoly(nnode, range, range, coords);

        // every once in a while for nnode=4, randomly decide to pick a unit square
        if (nnode == 4 && random_int(0, 3) == 0) {
            coords[0] = 0.; coords[1] = 0.;
            coords[2] = 1.; coords[3] = 0.;
            coords[4] = 1.; coords[5] = 1.;
            coords[6] = 0.; coords[7] = 1.;
        }

        for (int i = 0; i < nnode; i++) {
            nodelist[i] = i;  // Simply assign node indices
        }

        double norm[2], norm_sc, distance;

        // To test the variousness of node_loc, we randomly select a line
        // that can go through one or two of the nodes, or not

        double pt_a[2], pt_b[2]; // two points that the line will go through
        int n_a, n_b;            // potential node indices
        n_a = n_b = -1;

        // first point of the line
        if (random_int(0, 2) == 0) {
            n_a = random_int(0, nnode - 1);
            pt_a[0] = coords[2*n_a];
            pt_a[1] = coords[2*n_a + 1];
        }
        else {
            pt_a[0] = random_double(-range/4., range/4.);
            pt_a[1] = random_double(-range/4., range/4.);
        }

        // second point of the line
        if (random_int(0, 2) == 0) {
            int i;
            for (i = 0; i < 100; i++) {
                n_b = random_int(0, nnode - 1);
                if (n_b != n_a) break;
            }
            pt_b[0] = coords[2*n_b];
            pt_b[1] = coords[2*n_b + 1];
        }
        else {
            pt_b[0] = random_double(-range/4., range/4.);
            pt_b[1] = random_double(-range/4., range/4.);
        }


        // determine the normal vector
        norm[0] = pt_b[1] - pt_a[1];
        norm[1] = pt_a[0] - pt_b[0];
        norm_sc = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
        norm[0] /= norm_sc;
        norm[1] /= norm_sc;

        // determine the distance
        // equation of the line: (\vec{r} - \vec{r}_a)\cdot\vec{n} = 0
        // effectively states that the distance from point \vec[r} to the line
        // must be zero;
        // if \vec{r} == 0, then the distance is \vec{r}_a\cdot\vec{n}
        distance = norm[0]*pt_a[0] + norm[1]*pt_a[1];

        if (verbose) {
            printf("next test:\n");
            printf(" - n = {%12.5f, %12.5f}\n", norm[0], norm[1]);
            printf(" - d = %12.5f\n", distance);
        }
        for (int i = 0; i < nnode; i++) {
            if (i == n_b || i == n_a)
                node_loc[i] = 0;
            else {
                double dist = coords[2*i]*norm[0] + coords[2*i+1]*norm[1];
                node_loc[i] = (dist > 0) ? 1 : -1;
            }
            if (verbose) {
                printf(" - p[%d] = {%12.5f, %12.5f}, loc = %d\n",
                       i, coords[2*i], coords[2*i+1], node_loc[i]);
            }
        }

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs to the HDF5 file
        h5_write_1d(test_gid, "coords", coords, 2*nnode);       // Coordinates of the nodes
        h5_write_1d_int(test_gid, "nodelist", nodelist, nnode); // List of node indices
        h5_write_1d(test_gid, "norm", norm, 2);                 // Normal vector
        h5_write_1d(test_gid, "distance", &distance, 1);        // Distance of the interface
        h5_write_1d_int(test_gid, "node_loc", node_loc, nnode); // Classification of nodes

        // Outputs
        int *nodelist_lower = NULL, *nodelist_upper = NULL;

        // Call the actual function find_interface2d (assuming it exists)
        find_interface2d(nnode, coords, nodelist, norm, distance, node_loc,
                         &nnode_new, coords_new, &nnode_interface, nodes_interface,
                         &nnode_lower, &nodelist_lower, &nnode_upper, &nodelist_upper);

        // Write outputs to the HDF5 file
        h5_write_1d(test_gid, "coords_new", coords_new, 2*nnode_new);
        h5_write_1d_int(test_gid, "nodes_interface", nodes_interface, nnode_interface);
        h5_write_1d_int(test_gid, "nodelist_lower", nodelist_lower, nnode_lower);
        h5_write_1d_int(test_gid, "nodelist_upper", nodelist_upper, nnode_upper);

        if (verbose) {
            printf (" - nodes_upper: {");
            for (int i = 0; i < nnode_upper; ++i)
                printf("%3d", nodelist_upper[i]);
            printf ("}\n");
            printf (" - nodes_lower: {");
            for (int i = 0; i < nnode_lower; ++i)
                printf("%3d", nodelist_lower[i]);
            printf ("}\n");
            printf (" - nodes_interface: {");
            for (int i = 0; i < nnode_interface; ++i)
                printf("%3d", nodes_interface[i]);
            printf ("}\n");
            printf(" - new nodes:");
            if (nnode_new > 0) {
                printf("\n");
                for (int i = 0; i < nnode_new; ++i)
                    printf(" - n[%d] = {%12.5f, %12.5f}\n",
                           i + nnode, coords_new[2*i], coords_new[2*i+1]);
            }
            else {
                printf(" none.\n");
            }
        }

        // Free dynamically allocated memory
        free(nodelist_lower);
        free(nodelist_upper);

        // Close the group
        H5Gclose(test_gid);
    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);
}


// Function to record test cases for cal_poly_area
void record_test__cal_poly_area(int num_tests, const char *test_h5, const int verbose) {

    // HDF5 group handles
    hid_t group_id;
    const int nnode_max = 6;
    const double range = 5.0;
    double coords[2*nnode_max];
    int nodelist[nnode_max], node_loc[nnode_max];

    // outputs
    double area;

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    for (int t = 0; t < num_tests; t++) {
        // Randomly generate input parameters
        int nnode = random_int(3, 6);  // Number of nodes in the polygon
        generate_random_convex_2dpoly(nnode, range, range/1.1, coords);

        for (int i = 0; i < nnode; i++) {
            nodelist[i] = i;  // Simply assign node indices
        }

        // Call the function to record test results
        cal_poly_area(nnode, coords, nnode, nodelist, &area);

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs to the HDF5 file
        h5_write_1d(test_gid, "coords", coords, 2*nnode);       // Coordinates of the nodes
        h5_write_1d_int(test_gid, "nodelist", nodelist, nnode); // List of node indices

        // Write outputs to the HDF5 file
        h5_write_1d(test_gid, "area", &area, 1);                // Area of the polygon

        // Close the group
        H5Gclose(test_gid);
    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);
}


// Function to record test cases for rz_area
void record_test__rz_area(int num_tests, const char *test_h5, const int verbose) {

    // HDF5 group handles
    hid_t group_id;
    const int nn_max = 6;
    const double range = 5.0;
    double rz[2*nn_max];

    // outputs
    double vol;

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    for (int t = 0; t < num_tests; t++) {
        // Randomly generate input parameters
        int nn = random_int(3, 6);  // Number of nodes in the polygon
        generate_random_convex_2dpoly(nn, range, range, rz);

        // Call the function to record test results
        rz_area(nn, rz, &vol);

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs to the HDF5 file
        h5_write_1d(test_gid, "rz", rz, 2*nn);       // Coordinates of the nodes

        // Write outputs to the HDF5 file
        h5_write_1d(test_gid, "vol", &vol, 1);                // Area of the polygon

        // Close the group
        H5Gclose(test_gid);
    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);
}


// Function to record test cases for bounds_2d into an HDF5 file
void record_test__bounds_2d(int num_tests, const char *test_h5, const int verbose) {

    // HDF5 group handles
    const int nnode_max = 6;
    const double range = 5.0;
    double coords[2*nnode_max];
    int nodelist[nnode_max], node_loc[nnode_max];

    // outputs
    double coords_new[4*nnode_max];
    int nodes_interface[2*nnode_max];
    int node_order_for_ds[nnode_max];
    double ds_ea_node[nnode_max];

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    for (int t = 0; t < num_tests; t++) {

        // Randomly generate input polygon
        int geop = 1;
        double vf_to_match = random_double(0., 1.);
        double volume = random_double(1e-3, 10);
        int nnode = random_int(3, 6);  // Number of nodes in the polygon
        generate_random_convex_2dpoly(nnode, range, range, coords);

        // every once in a while for nnode=4, randomly decide to pick a unit square
        if (nnode == 4 && random_int(0, 3) == 0) {
            coords[0] = 0.; coords[1] = 0.;
            coords[2] = 1.; coords[3] = 0.;
            coords[4] = 1.; coords[5] = 1.;
            coords[6] = 0.; coords[7] = 1.;
        }

        // sometimes randomly decide to make vf_to_match 0 or 1
        if (random_int(0, 3) == 0)
            vf_to_match = (double)random_int(0, 1);

        for (int i = 0; i < nnode; i++) {
            nodelist[i] = i;  // Simply assign node indices
        }

        double norm[2], norm_sc;

        // determine the normal vector
        norm[0] = random_double(-1., 1.);
        norm[1] = random_double(-1., 1.);
        norm_sc = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
        norm[0] /= norm_sc;
        norm[1] /= norm_sc;

        order_nodes_along_norm(2, norm, nnode, coords, node_order_for_ds, ds_ea_node);

        if (verbose) {
            printf("next test:\n");
            printf(" - vf_to_match = %14.5e\n", vf_to_match);
            printf(" - volume = %14.5e\n", volume);
            printf(" - n = {%12.5f, %12.5f}\n", norm[0], norm[1]);
            for (int i = 0; i < nnode; i++) {
                printf(" - p[%d] = {%12.5f, %12.5f}, ds = %12.5f\n",
                       i, coords[2*i], coords[2*i+1], ds_ea_node[i]);
            }
            printf(" - node_order_for_ds: {");
            for (int i = 0; i < nnode; i++) {
                printf("%3d ", node_order_for_ds[i]);
            }
            printf("}\n");
        }
        // Outputs
        int *nodelist_lower = NULL, *nodelist_upper = NULL;
        int nnode_new, nnode_interface, nnode_lower, nnode_upper;
        double ds_lower, ds_upper, vf_lower, vf_upper;
        int vol_matched;

        // Call the function
        bounds_2d(geop, vf_to_match, volume, norm,
                       nnode, coords, nodelist,
                       node_order_for_ds, ds_ea_node,
                       &ds_lower, &ds_upper,
                       &vf_lower, &vf_upper,
                       &nnode_new, coords_new,
                       &vol_matched,
                       &nnode_interface, nodes_interface,
                       &nnode_lower, &nodelist_lower,
                       &nnode_upper, &nodelist_upper);

        if (verbose) {
            printf ("outputs:\n");
            printf (" - ds(lower, upper) = {%12.5f, %12.5f}\n", ds_lower, ds_upper);
            printf (" - vf(lower, upper) = {%12.5f, %12.5f}\n", vf_lower, vf_upper);
            printf (" - vol_matched = %d\n", vol_matched);
            printf (" - nodelist_upper: {");
            for (int i = 0; i < nnode_upper; ++i)
                printf("%3d", nodelist_upper[i]);
            printf ("}\n");
            printf (" - nodelist_lower: {");
            for (int i = 0; i < nnode_lower; ++i)
                printf("%3d", nodelist_lower[i]);
            printf ("}\n");
            printf (" - nodes_interface: {");
            for (int i = 0; i < nnode_interface; ++i)
                printf("%3d", nodes_interface[i]);
            printf ("}\n");
            printf(" - new nodes:");
            if (nnode_new > 0) {
                printf("\n");
                for (int i = 0; i < nnode_new; ++i)
                    printf(" - n[%d] = {%12.5f, %12.5f}\n",
                           i + nnode, coords_new[2*i], coords_new[2*i+1]);
            }
            else {
                printf(" none.\n");
            }
        }

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs to the HDF5 file
        h5_write_1d_int(test_gid, "geop", &geop, 1);
        h5_write_1d(test_gid, "vf_to_match", &vf_to_match, 1);
        h5_write_1d(test_gid, "volume", &volume, 1);
        h5_write_1d(test_gid, "norm", norm, 2);
        h5_write_1d_int(test_gid, "nnode", &nnode, 1);
        h5_write_1d(test_gid, "coords", coords, 2*nnode);
        h5_write_1d_int(test_gid, "nodelist", nodelist, nnode);
        h5_write_1d_int(test_gid, "node_order_for_ds", node_order_for_ds, nnode);
        h5_write_1d(test_gid, "ds_ea_node", ds_ea_node, nnode);

        // Write outputs to the HDF5 file
        h5_write_1d(test_gid, "ds_lower", &ds_lower, 1);
        h5_write_1d(test_gid, "ds_upper", &ds_upper, 1);
        h5_write_1d(test_gid, "vf_lower", &vf_lower, 1);
        h5_write_1d(test_gid, "vf_upper", &vf_upper, 1);
        h5_write_1d(test_gid, "coords_new", coords_new, 2*nnode_new);
        h5_write_1d_int(test_gid, "vol_matched", &vol_matched, 1);
        h5_write_1d_int(test_gid, "nodes_interface", nodes_interface, nnode_interface);
        h5_write_1d_int(test_gid, "nodelist_lower", nodelist_lower, nnode_lower);
        h5_write_1d_int(test_gid, "nodelist_upper", nodelist_upper, nnode_upper);

        // Free dynamically allocated memory
        free(nodelist_lower);
        free(nodelist_upper);

        // Close the group
        H5Gclose(test_gid);
    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);
} // record_test__bounds_2d


// Function to record test cases for order_nodes_along_norm
void record_test__order_nodes_along_norm(int num_tests, const char *test_h5, const int verbose) {

    // HDF5 group handles
    hid_t group_id;
    const int nnode_max = 6;
    const double range = 5.0;
    double coords[2*nnode_max];

    // outputs
    int node_order[nnode_max];
    double ds_ea_node[nnode_max];

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    for (int t = 0; t < num_tests; t++) {
        // Randomly generate input parameters
        int dim = 2;
        double norm[2];
        norm[0] = random_double(-1.0, 1.0);  // Random normal vector
        norm[1] = random_double(-1.0, 1.0);

        int nnode = random_int(3, 6);  // Number of nodes in the polygon
        generate_random_convex_2dpoly(nnode, range, range/1.5, coords);

        // Call the function to record test results
        order_nodes_along_norm(dim, norm, nnode, coords, node_order, ds_ea_node);

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs to the HDF5 file
        h5_write_1d_int(test_gid, "dim", &dim, 1);
        h5_write_1d(test_gid, "norm", norm, 2);
        h5_write_1d_int(test_gid, "nnode", &nnode, 1);
        h5_write_1d(test_gid, "coords", coords, 2*nnode);

        // Write outputs to the HDF5 file
        h5_write_1d_int(test_gid, "node_order", node_order, nnode);
        h5_write_1d(test_gid, "ds_ea_node", ds_ea_node, nnode);

        // Close the group
        H5Gclose(test_gid);
    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);
}


// Function to record test cases for cal_distance2d into an HDF5 file
void record_test__cal_distance2d(int num_tests, const char *test_h5, const int verbose) {

    // HDF5 group handles
    const int nnode_max = 8;
    const double range = 5.0;
    double coords[2*nnode_max];
    int nodelist[nnode_max], node_loc[nnode_max];

    // outputs
    double coords_new[2*nnode_max];
    int nodes_interface[2*nnode_max];

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    for (int t = 0; t < num_tests; t++) {
        // Randomly generate input parameters
        int geop = 1;
        double vf_to_match = random_double(0., 1.);
        double volume = random_double(1e-3, 10);

        double norm[2], norm_sc;

        // determine the normal vector
        norm[0] = random_double(-1., 1.);
        norm[1] = random_double(-1., 1.);
        norm_sc = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
        norm[0] /= norm_sc;
        norm[1] /= norm_sc;


        int nnode = random_int(3, 4);  // Number of nodes in the polygon
        generate_random_convex_2dpoly(nnode, range, range*0.7, coords);

        // every once in a while for nnode=4, randomly decide to pick a unit square
        if (nnode == 4 && random_int(0, 3) == 0) {
            coords[0] = 0.; coords[1] = 0.;
            coords[2] = 1.; coords[3] = 0.;
            coords[4] = 1.; coords[5] = 1.;
            coords[6] = 0.; coords[7] = 1.;
        }
        for (int i = 0; i < nnode; i++) {
            nodelist[i] = i;
        }

        if (verbose) {
            printf("next test:\n");
            printf(" - vf_to_match = %14.5e\n", vf_to_match);
            printf(" - volume = %14.5e\n", volume);
            printf(" - n = {%12.5f, %12.5f}\n", norm[0], norm[1]);
            for (int i = 0; i < nnode; i++) {
                printf(" - p[%d] = {%12.5f, %12.5f}\n",
                       i, coords[2*i], coords[2*i+1]);
            }
        }

        // Outputs
        int nnode_new, nnode_interface, nnode_lower, nnode_upper;
        double distance;
        int *nodelist_lower = NULL, *nodelist_upper = NULL;

        cal_distance2d(geop, vf_to_match, volume, norm,
                       nnode, coords, nodelist,
                       &nnode_new, coords_new,
                       &distance,
                       &nnode_interface, nodes_interface,
                       &nnode_lower, &nodelist_lower,
                       &nnode_upper, &nodelist_upper);

        if (verbose) {
            printf ("outputs:\n");
            printf(" - new nodes:");
            if (nnode_new > 0) {
                printf("\n");
                for (int i = 0; i < nnode_new; ++i)
                    printf(" - n[%d] = {%12.5f, %12.5f}\n",
                           i + nnode, coords_new[2*i], coords_new[2*i+1]);
            }
            else {
                printf(" none.\n");
            }
            printf (" - distance = %12.5f\n", distance);
            printf (" - nodes_interface: {");
            for (int i = 0; i < nnode_interface; ++i)
                printf("%3d", nodes_interface[i]);
            printf ("}\n");
            printf (" - nodelist_upper: {");
            for (int i = 0; i < nnode_upper; ++i)
                printf("%3d", nodelist_upper[i]);
            printf ("}\n");
            printf (" - nodelist_lower: {");
            for (int i = 0; i < nnode_lower; ++i)
                printf("%3d", nodelist_lower[i]);
            printf ("}\n");
        }

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs to the HDF5 file
        h5_write_1d_int(test_gid, "geop", &geop, 1);
        h5_write_1d(test_gid, "vf_to_match", &vf_to_match, 1);
        h5_write_1d(test_gid, "volume", &volume, 1);
        h5_write_1d(test_gid, "norm", norm, 2);
        h5_write_1d_int(test_gid, "nnode", &nnode, 1);
        h5_write_1d(test_gid, "coords", coords, 2*nnode);
        h5_write_1d_int(test_gid, "nodelist", nodelist, nnode);

        // Write outputs to the HDF5 file
        h5_write_1d(test_gid, "coords_new", coords_new, 2*nnode_new);
        h5_write_1d(test_gid, "distance", &distance, 1);
        h5_write_1d_int(test_gid, "nodes_interface", nodes_interface, nnode_interface);
        h5_write_1d_int(test_gid, "nodelist_lower", nodelist_lower, nnode_lower);
        h5_write_1d_int(test_gid, "nodelist_upper", nodelist_upper, nnode_upper);

        // Free dynamically allocated memory
        free(nodelist_lower);
        free(nodelist_upper);

        // Close the group
        H5Gclose(test_gid);
    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);
} // record_test__cal_distance2d


// Function to record test cases for reconstruct2d_nmat_pagosa into an HDF5 file
void record_test__reconstruct2d_nmat_pagosa(int num_tests, const char *test_h5, const int verbose) {

    // HDF5 group handles
    const int nnode_max = 8;
    const double range = 5.0;
    double coords[2*nnode_max];
    int nodelist[nnode_max], node_loc[nnode_max];

    // outputs
    double coords_new[2*nnode_max];
    int nodes_interface[2*nnode_max];

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    for (int t = 0; t < num_tests; t++) {
        // Randomly generate input parameters
        int geop = 1;
        double xl[2], dx[2];

        xl[0] = 0.;//random_double(-10., 10.);
        xl[1] = 0.;//random_double(-10., 10.);
        dx[0] = random_double(0., 1.);
        dx[1] = dx[0];

        // generate random reference point for a straight line
        double xr[2];
        xr[0] = random_double(0., dx[0]);
        xr[1] = random_double(0., dx[1]);
        int dice_roll = random_int(0, 12000);
        switch (dice_roll) {
        case 0:
            xr[0] = 0.; xr[1] = 0.;
            break;
        case 1:
            xr[0] = 0.; xr[1] = dx[1];
            break;
        case 2:
            xr[0] = dx[0]; xr[1] = 0.;
            break;
        case 3:
            xr[0] = dx[0]; xr[1] = dx[1];
            break;
        }

        // determine the normal vector
        double norm[2], norm_sc;
        norm[0] = random_double(-1., 1.);
        norm[1] = random_double(-1., 1.);

        // randomly pick one components to be zero
        dice_roll = random_int(0, 5000);
        switch (dice_roll) {
        case 0:
            norm[0] = 0.;
            break;
        case 1:
            norm[1] = 0.;
            break;
        }

        norm_sc = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
        norm[0] /= norm_sc;
        norm[1] /= norm_sc;


        // corner points
        double cn[4][4];
        double cn0 = xr[0]*norm[0] + xr[1]*norm[1];
        for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            cn[i][j] = (i - 1)*norm[0]*dx[0] + (j - 1)*norm[1]*dx[1] - cn0;

        int nmmax = 2;
        int ***matid_2dsmesh = (int ***) malloc(3 * sizeof(int **));
        int  **ibuffer2d     = (int  **) malloc(9 * sizeof(int *));
        int   *ibuffer1d     = (int   *) malloc(9 * nmmax * sizeof(int));
        int sizes[] = {nmmax, 3, 3};
        ASSIGN_3D_FORM(int, matid_2dsmesh, ibuffer2d, ibuffer1d, sizes);

        int **nmat_2dsmesh = (int **) malloc(3 * sizeof(int *));
        nmat_2dsmesh[0] = (int *) malloc(9 * sizeof(int));
        for (int i = 1; i < 3; i++) {
            nmat_2dsmesh[i] = nmat_2dsmesh[i-1] + 3;
        }

        double ***vf_2dsmesh = (double ***) malloc(3 * sizeof(double **));
        double  **dbuffer2d  = (double  **) malloc(9 * sizeof(double *));
        double   *dbuffer1d  = (double   *) malloc(9 * nmmax * sizeof(double));
        ASSIGN_3D_FORM(double, vf_2dsmesh, dbuffer2d, dbuffer1d, sizes);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                // nmat_2dsmesh = ...
                // matid_2dsmesh = ...
                // vf_2dsmesh = ...
                double corners[4], aux;
                int corners_order[4] = {0, 1, 2, 3};
                corners[0] = cn[i][j];
                corners[1] = cn[i+1][j];
                corners[2] = cn[i][j+1];
                corners[3] = cn[i+1][j+1];
                vf_2dsmesh[i][j][0] = 0.;

                if (corners[0] > 0 && corners[1] > 0 && corners[2] > 0 && corners[3] > 0) {
                    nmat_2dsmesh[i][j] = 1;
                    vf_2dsmesh[i][j][0] = 1.;
                    matid_2dsmesh[i][j][0] = 1.;
                }
                else if (corners[0] < 0 && corners[1] < 0 && corners[2] < 0 && corners[3] < 0) {
                    nmat_2dsmesh[i][j] = 1;
                    vf_2dsmesh[i][j][0] = 1.;
                    matid_2dsmesh[i][j][0] = 0.;
                }
                else {
                    nmat_2dsmesh[i][j] = 2;
                    if (corners[0] > corners[1]) {
                        aux = corners[0]; corners[0] = corners[1]; corners[1] = aux;
                        aux = corners[2]; corners[2] = corners[3]; corners[3] = aux;
                    }
                    if (corners[0] > corners[2]) {
                        aux = corners[0]; corners[0] = corners[2]; corners[2] = aux;
                        aux = corners[1]; corners[1] = corners[3]; corners[3] = aux;
                    }
                    if (corners[0] > corners[3]) {
                        aux = corners[0]; corners[0] = corners[3]; corners[3] = aux;
                        aux = corners[2]; corners[2] = corners[1]; corners[1] = aux;
                    }
                    if (corners[1] > corners[2]) {
                        aux = corners[1]; corners[1] = corners[2]; corners[2] = aux;
                    }
                    double h = 1./(corners[2] - corners[0]);
                    double dx = min(0, corners[1]) - corners[0];
                    double tol = 1e-12;
                    vf_2dsmesh[i][j][0] = 0.5*dx*dx/(corners[1] - corners[0] + tol)*h;
                    if (corners[1] < tol) {
                        dx = min(0., corners[2]) - corners[1];
                        vf_2dsmesh[i][j][0] += dx*h;
                    }
                    if (corners[2] < tol) {
                        dx = corners[3];
                        vf_2dsmesh[i][j][0] = 1. - .5*dx*dx/(corners[3] - corners[2] + tol)*h;
                    }
                    vf_2dsmesh[i][j][1] = 1. - vf_2dsmesh[i][j][0];

                    matid_2dsmesh[i][j][0] = 0;
                    matid_2dsmesh[i][j][1] = 1;
                }
            }
        }

        if (verbose) {
            printf("next test (%d):\n", t);
            printf(" - xr = {%9.5f, %9.5f}, norm = {%9.5f, %9.5f}\n", xr[0], xr[1], norm[0], norm[1]);
            printf(" - xl = {%7.4f, %7.4f}\n", xl[0], xl[1]);
            printf(" - dx = {%7.4f, %7.4f}\n", dx[0], dx[1]);
            printf(" - nmat_2dsmesh:\n");
            for (int i = 0; i < 3; i++) {
                printf ("               %2d %2d %2d\n",
                    nmat_2dsmesh[i][0], nmat_2dsmesh[i][1], nmat_2dsmesh[i][2]);
            }
            printf(" - vf_2dsmesh:\n");
            for (int i = 0; i < 3; i++) {
                printf ("               %7.4f %7.4f %7.4f\n",
                    vf_2dsmesh[i][0][0], vf_2dsmesh[i][1][0], vf_2dsmesh[i][2][0]);
            }
        }

        // Outputs
        int nnode_tot = 0;
        double *coords_tot = NULL;
        int *nnode_for_minterface, **nodes_for_minterface;
        int *nnode_for_mpoly, **nodes_for_mpoly;

        int nm = nmmax;
        nnode_for_minterface    = (int  *) malloc((nm - 1) * sizeof(int));
        nodes_for_minterface    = (int **) malloc((nm - 1) * sizeof(int *));
        nodes_for_minterface[0] = (int  *) malloc((nm - 1) * 2 * sizeof(int));
        for (int m = 1; m < nm-1; m++) {
            nodes_for_minterface[m] = nodes_for_minterface[m-1] + 2;
        }
        nnode_for_mpoly = (int *) malloc(nm * sizeof(int));
        nodes_for_mpoly = NULL;

        reconstruct2d_nmat_pagosa(geop, xl, dx, nmat_2dsmesh, matid_2dsmesh, vf_2dsmesh,
                   &nnode_tot, &coords_tot,
                   nnode_for_minterface, nodes_for_minterface,
                   nnode_for_mpoly, &nodes_for_mpoly);

        if (verbose) {
            printf("outputs:\n");
            printf(" - nnode_tot = %d\n", nnode_tot);
            printf(" - coords_tot:");
            if (nnode_tot > 0) {
                for (int i = 0; i < nnode_tot; ++i)
                    printf("\n   %d -> {%12.5f, %12.5f}",
                           i, coords_tot[2*i], coords_tot[2*i+1]);
                printf("\n");
            }
            else {
                printf(" none.\n");
            }
            printf(" - nnode_for_minterface = %d\n", nnode_for_minterface[0]);
            printf(" - nodes_for_minterface: {");
            for (int m = 0; m < nm-1; m++)
                printf("%3d %3d; ", nodes_for_minterface[m][0], nodes_for_minterface[m][1]);
            printf ("}\n");
            printf(" - nnode_for_mpoly = {");
            for (int m = 0; m < nm; m++)
                printf("%3d", nnode_for_mpoly[m]);
            printf("}\n");

            for (int m = 0; m < nm; m++) {
                printf(" - nodes_for_mpoly %d: {", m);
                for (int k = 0; k < nnode_for_mpoly[m]; k++)
                    printf("%3d ", nodes_for_mpoly[m][k]);
                printf ("}\n");
            }

        }

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs to the HDF5 file
        h5_write_1d_int(test_gid, "geop", &geop, 1);
        h5_write_1d(test_gid, "xl", xl, 2);
        h5_write_1d(test_gid, "dx", dx, 2);
        h5_write_2d_int(test_gid, "nmat_2dsmesh", nmat_2dsmesh, &(sizes[1]));
        h5_write_1d(test_gid, "vf_2dsmesh", &(vf_2dsmesh[0][0][0]), 3*3*nm);
        h5_write_1d_int(test_gid, "matid_2dsmesh", &(matid_2dsmesh[0][0][0]), 3*3*nm);

        h5_write_1d_int(test_gid, "nnode_tot", &nnode_tot, 1);
        h5_write_1d(test_gid, "coords_tot", coords_tot, 2*nnode_tot);
        h5_write_1d_int(test_gid, "nnode_for_minterface", nnode_for_minterface, nm-1);
        h5_write_1d_int(test_gid, "nodes_for_minterface", nodes_for_minterface[0], (nm-1)*2);
        h5_write_1d_int(test_gid, "nnode_for_mpoly", nnode_for_mpoly, nm);
        int nmtot = 0;
        for (int m = 0; m < nm; m++)
            nmtot += nnode_for_mpoly[m];
        h5_write_1d_int(test_gid, "nodes_for_mpoly", nodes_for_mpoly[0], nmtot);

        // Close the group
        H5Gclose(test_gid);

        // Free dynamically allocated memory
        free(vf_2dsmesh[0][0]); free(vf_2dsmesh[0]); free(vf_2dsmesh);
        free(matid_2dsmesh[0][0]); free(matid_2dsmesh[0]); free(matid_2dsmesh);
        free(nmat_2dsmesh[0]); free(nmat_2dsmesh);

        free(coords_tot);
        free(nnode_for_minterface);
        free(nodes_for_minterface[0]); free(nodes_for_minterface);
        free(nnode_for_mpoly);

    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);

} // record_test__reconstruct2d_nmat_pagosa


void record_test__cal_cell_zgrad2d(int num_tests, const char *test_h5, const int verbose) {

    // Allocate memory for inputs and outputs
    double dx[2], var[3][3], var9[9], grad[2];

    // Create HDF5 file
    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);
    char group_name[50];

    // Randomly generate input arguments and call rec_rec
    for (int t = 0; t < num_tests; t++) {

        dx[0] = random_double(1e-9, 1.);
        dx[1] = random_double(1e-9, 1.);
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            var[i][j] = random_double(-10., 10.);
            var9[i + 3*j] = var[i][j];
        }

        // Create a group for this test case
        sprintf(group_name, "/test_%d", t);
        hid_t test_gid = h5_open_group(testfile_id, group_name);

        // Write inputs
        h5_write_1d(test_gid, "dx", dx, 2);
        h5_write_1d(test_gid, "var", var9, 9);

        cal_cell_zgrad2d(dx, var, grad);
        if (verbose) {
            printf("dx[0] = {%f, %f}\n", dx[0], dx[1]);
        }

        // Write outputs
        h5_write_1d(test_gid, "grad", grad, 2);

        // Close the group
        H5Gclose(test_gid);
    }

    // Close HDF5 file
    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);

} // record_test__cal_cell_zgrad2d


void record_test__bdry_cell_1var_2d(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    int **nmat_for_2dcell, i, j;
    Bdry_Type btype[2] = {bdry_transmitted, bdry_transmitted};

    char group_name[50];
    hid_t mesh_gid, attr_id, test_gid;
    int num_meshes, imesh, d;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        nmat_for_2dcell = h5_read_2d_int(mesh_gid, "nmat",ncell_ext);

        // Call the function under test
        bdry_cell_1var_2d(ncell, nbdry, btype, btype, nmat_for_2dcell);

        // Create a new group for each mesh
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_2d_int(test_gid, "nmat", nmat_for_2dcell, ncell_ext);
        H5Gclose(test_gid);

        H5Gclose(mesh_gid);
        free(nmat_for_2dcell[0]); free(nmat_for_2dcell);
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__bdry_cell_1var_2d


void record_test__bdry_cell_ragged_2d(const char *meshes_h5, const char *test_h5, const int verbose) {

    int dim, ncell[3], nbdry, ncell_ext[3], nnode[3], nnode_ext[3];
    double dx[3], dvol;
    int **nmat_for_cell, i, j;
    Bdry_Type btype[2] = {bdry_transmitted, bdry_transmitted};

    char group_name[50];
    hid_t mesh_gid, attr_id, test_gid;
    int num_meshes, imesh, d, mxi, idx;
    int nmat_max, lsize_mat, lsize_cell, nm_cell;

    int nm_new, offset;
    int *matid_list, ***matid_for_cell, **matid_for_2dcell;
    double *vol_list, *mass_list, *ener_list, vol_cell;
    double ***vol_for_cell, ***mass_for_cell, ***ener_for_cell;
    double **rho_for_2dcell, **ei_for_2dcell, **pres_for_2dcell;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    hid_t file_meshes_id = h5_open_existing_rdonly(meshes_h5);
    num_meshes = h5_count_groups(file_meshes_id);
    for (imesh = 0; imesh < num_meshes; ++imesh) {

        snprintf(group_name, sizeof(group_name), "mesh_%d", imesh);
        mesh_gid = h5_open_group(file_meshes_id, group_name);

        // read mesh attributes
        dim = h5_read_attr_int(mesh_gid, "dim");
        nbdry = h5_read_attr_int(mesh_gid, "nbdry");
        h5_read_attr_1d_int(mesh_gid, "ncell", ncell, dim);
        h5_read_attr_1d(mesh_gid, "dx", dx, dim);

        for (d = 0; d < dim; ++d) {
            ncell_ext[d] = ncell[d] + 2*nbdry;
            nnode[d] = ncell[d] + 1;
            nnode_ext[d] = nnode[d] + 2*nbdry;
        }

        nmat_for_cell = h5_read_2d_int(mesh_gid, "nmat", ncell_ext);
        matid_for_2dcell = h5_read_2d_int(mesh_gid, "matid", ncell_ext);
        rho_for_2dcell  = h5_read_2d(mesh_gid, "rho", ncell_ext);
        ei_for_2dcell   = h5_read_2d(mesh_gid, "ei",  ncell_ext);

        // Apply boundary conditions for nmat so that boundary cells are correctly allocated
        bdry_cell_1var_2d(ncell, nbdry, btype, btype, nmat_for_cell);

        // determine the maximum number of materials and length of the ragged array
        lsize_mat = 0;
        nmat_max = 0;
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {
                lsize_mat += nmat_for_cell[j][i];
                if (nmat_for_cell[j][i] > nmat_max) nmat_max = nmat_for_cell[j][i];
            }
        }
        matid_for_cell  = (int    ***) malloc(ncell_ext[1] * sizeof(int **));
        vol_for_cell    = (double ***) malloc(ncell_ext[1] * sizeof(double **));
        mass_for_cell   = (double ***) malloc(ncell_ext[1] * sizeof(double **));
        ener_for_cell   = (double ***) malloc(ncell_ext[1] * sizeof(double **));

        lsize_cell = ncell_ext[1]*ncell_ext[0];
        matid_for_cell[0] = (int    **) malloc(lsize_cell * sizeof(int *));
        vol_for_cell[0]   = (double **) malloc(lsize_cell * sizeof(double *));
        mass_for_cell[0]  = (double **) malloc(lsize_cell * sizeof(double *));
        ener_for_cell[0]  = (double **) malloc(lsize_cell * sizeof(double *));
        for (j = 1; j < ncell_ext[1]; j++) {
            matid_for_cell[j] = matid_for_cell[j-1] + ncell_ext[0];
            vol_for_cell[  j] = vol_for_cell[  j-1] + ncell_ext[0];
            mass_for_cell[ j] = mass_for_cell[ j-1] + ncell_ext[0];
            ener_for_cell[ j] = ener_for_cell[ j-1] + ncell_ext[0];
        }
        matid_for_cell[0][0] = (int    *) malloc(lsize_mat * sizeof(int));
        vol_for_cell[0][0]   = (double *) malloc(lsize_mat * sizeof(double));
        mass_for_cell[0][0]  = (double *) malloc(lsize_mat * sizeof(double));
        ener_for_cell[0][0]  = (double *) malloc(lsize_mat * sizeof(double));

        offset = 0;
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {
                if (i > 0 || j > 0) {
                   matid_for_cell[j][i] = matid_for_cell[0][0] + offset;
                     vol_for_cell[j][i] =   vol_for_cell[0][0] + offset;
                    mass_for_cell[j][i] =  mass_for_cell[0][0] + offset;
                    ener_for_cell[j][i] =  ener_for_cell[0][0] + offset;
                }
                offset += nmat_for_cell[j][i];
            }
        }

        vol_cell = dx[0]*dx[1];
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {
                matid_for_cell[j][i][0] = matid_for_2dcell[j][i];
                vol_for_cell[j][i][0]   = vol_cell;
                mass_for_cell[j][i][0]  = rho_for_2dcell[j][i] * vol_cell;
                ener_for_cell[j][i][0]  = ei_for_2dcell[j][i]  * vol_cell;
                nm_cell = nmat_for_cell[j][i];

                for (idx = 1; idx < nm_cell; idx++) {
                    matid_for_cell[j][i][idx] = random_int(0, nmat_max);
                    dvol                      = vol_cell * random_double(0, 1);
                    vol_for_cell[ j][i][idx]  = dvol;
                    mass_for_cell[j][i][idx]  = dvol * rho_for_2dcell[j][i];
                    ener_for_cell[j][i][idx]  = dvol * ei_for_2dcell[j][i];
                }
            }
        }

        // Record the input data
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_2d_int(test_gid, "nmat_before", nmat_for_cell, ncell_ext);
        h5_write_1d_int(test_gid, "matid_before", matid_for_cell[0][0], lsize_mat);
        h5_write_1d(test_gid, "vol_before", vol_for_cell[0][0], lsize_mat);
        h5_write_1d(test_gid, "mass_before", mass_for_cell[0][0], lsize_mat);
        h5_write_1d(test_gid, "ener_before", ener_for_cell[0][0], lsize_mat);

        // Call the function under test
        bdry_cell_ragged_2d(ncell, nbdry, btype, btype,
                         nmat_for_cell, matid_for_cell, vol_for_cell, mass_for_cell, ener_for_cell);

        // Record function output
        h5_write_2d_int(test_gid, "nmat_after", nmat_for_cell, ncell_ext);
        h5_write_1d_int(test_gid, "matid_after", matid_for_cell[0][0], lsize_mat);
        h5_write_1d(test_gid, "vol_after", vol_for_cell[0][0], lsize_mat);
        h5_write_1d(test_gid, "mass_after", mass_for_cell[0][0], lsize_mat);
        h5_write_1d(test_gid, "ener_after", ener_for_cell[0][0], lsize_mat);


        H5Gclose(test_gid);
        H5Gclose(mesh_gid);
        free(nmat_for_cell[0]);    free(nmat_for_cell);
        free(rho_for_2dcell[0]);   free(rho_for_2dcell);
        free(ei_for_2dcell[0]);    free(ei_for_2dcell);
        free(matid_for_2dcell[0]); free(matid_for_2dcell);

        free(vol_for_cell[0][0]);
        free(mass_for_cell[0][0]);
        free(ener_for_cell[0][0]);
        free(matid_for_cell[0][0]);

        free(vol_for_cell[0]);
        free(mass_for_cell[0]);
        free(ener_for_cell[0]);
        free(matid_for_cell[0]);

        free(vol_for_cell);
        free(mass_for_cell);
        free(ener_for_cell);
        free(matid_for_cell);
    }

    H5Fclose(testfile_id);
    H5Fclose(file_meshes_id);
    printf("written: %s\n", test_h5);

} // record_test__bdry_cell_ragged_2d


void generate_two_material_mesh(int *dim, int *nbdry, int *ncell, int *nnode,
        int *ncell_ext, int *nnode_ext, double *dx, double *xl,
        int ***nmat_for_2dcell, int ***matid_for_2dcell, int ***mixcell_for_2dcell,
        double ***vf_for_2dcell, int *nmixcell, int **nmat_in_mixcell, int ***matids_in_mixcell,
        double ***vf_in_mixcell, int ***ijk_in_mixcell, double ****vfgrad_in_mixcell) {

    int d, i, j, mx, offset;
    double **cn, cn0;

    // read mesh attributes
    *dim = 2;
    *nbdry = 2; // don't bother with other possibilities for now
    ncell[0] = random_int(1,10);
    ncell[1] = random_int(1,10);
    dx[0] = random_double(0., 1.);
    dx[1] = dx[0]; // set dx[0] == dx[1] for simplicity
    xl[0] = 0.; xl[1] = 0.;

    for (d = 0; d < *dim; ++d) {
        ncell_ext[d] = ncell[d] + 2*(*nbdry);
        nnode[d] = ncell[d] + 1;
        nnode_ext[d] = nnode[d] + 2*(*nbdry);
    }

    // allocate nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell,
    // and vf_for_2dcell
    *nmat_for_2dcell = (int **) malloc(ncell_ext[1]*sizeof(int *));
    *matid_for_2dcell = (int **) malloc(ncell_ext[1]*sizeof(int *));
    *mixcell_for_2dcell = (int **) malloc(ncell_ext[1]*sizeof(int *));
    *vf_for_2dcell = (double **) malloc(ncell_ext[1]*sizeof(double *));

    int lsize_mat = ncell_ext[1]*ncell_ext[0];
    (*nmat_for_2dcell)[0] = (int *) malloc(lsize_mat*sizeof(int));
    (*matid_for_2dcell)[0] = (int *) malloc(lsize_mat*sizeof(int));
    (*mixcell_for_2dcell)[0] = (int *) malloc(lsize_mat*sizeof(int));
    (*vf_for_2dcell)[0] = (double *) malloc(lsize_mat*sizeof(double));

    for (j = 1; j < ncell_ext[1]; j++) {
        (*nmat_for_2dcell)[j] = (*nmat_for_2dcell)[j-1] + ncell_ext[0];
        (*matid_for_2dcell)[j] = (*matid_for_2dcell)[j-1] + ncell_ext[0];
        (*mixcell_for_2dcell)[j] = (*mixcell_for_2dcell)[j-1] + ncell_ext[0];
        (*vf_for_2dcell)[j] = (*vf_for_2dcell)[j-1] + ncell_ext[0];
    }

    // use a random straight line and to divide two materials
    double xr[2]; // reference point for the line to go through
    xr[0] = random_double(0., dx[0]*ncell[0]);
    xr[1] = random_double(0., dx[1]*ncell[1]);

    // random normal vector
    double norm[2], norm_sc;
    norm[0] = random_double(-1., 1.);
    norm[1] = random_double(-1., 1.);

    // randomly pick one components to be zero
    int dice_roll = random_int(0, 5000);
    switch (dice_roll) {
    case 0:
        norm[0] = 0.;
        break;
    case 1:
        norm[1] = 0.;
        break;
    }

    norm_sc = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0] /= norm_sc;
    norm[1] /= norm_sc;

    cn0 = xr[0]*norm[0] + xr[1]*norm[1];
    cn = (double**) malloc(nnode_ext[1]*sizeof(double*));
    for (j = 0; j < nnode_ext[1]; j++){
        if (j == 0)
            cn[0] = (double *) malloc(nnode_ext[1]*nnode_ext[0]*sizeof(double));
        else
            cn[j] = cn[j-1] + nnode_ext[0];

        for (i = 0; i < nnode_ext[0]; i++)
            cn[j][i] = (i - *nbdry)*norm[0]*dx[0] + (j - *nbdry)*norm[1]*dx[1] - cn0;
    }

    (*nmixcell) = 0; // also count the mixcells
    for (j = 0; j < ncell_ext[1]; j++) {
        for (i = 0; i < ncell_ext[0]; i++) {
            // nmat_2dsmesh = ...
            // matid_2dsmesh = ...
            // vf_2dsmesh = ...
            double corners[4], aux;
            int corners_order[4] = {0, 1, 2, 3};
            corners[0] = cn[j][i];
            corners[1] = cn[j+1][i];
            corners[2] = cn[j][i+1];
            corners[3] = cn[j+1][i+1];
            (*vf_for_2dcell)[j][i] = 0.;

            if (corners[0] > 0 && corners[1] > 0 && corners[2] > 0 && corners[3] > 0) {
                (*nmat_for_2dcell)[j][i] = 1;
                (*matid_for_2dcell)[j][i] = 1;
                (*vf_for_2dcell)[j][i] = 1.;
            }
            else if (corners[0] < 0 && corners[1] < 0 && corners[2] < 0 && corners[3] < 0) {
                (*nmat_for_2dcell)[j][i] = 1;
                (*matid_for_2dcell)[j][i] = 1;
                (*vf_for_2dcell)[j][i] = 0.;
            }
            else {
                (*nmat_for_2dcell)[j][i] = 2;
                (*matid_for_2dcell)[j][i] = 0;
                if (corners[0] > corners[1]) {
                    aux = corners[0]; corners[0] = corners[1]; corners[1] = aux;
                    aux = corners[2]; corners[2] = corners[3]; corners[3] = aux;
                }
                if (corners[0] > corners[2]) {
                    aux = corners[0]; corners[0] = corners[2]; corners[2] = aux;
                    aux = corners[1]; corners[1] = corners[3]; corners[3] = aux;
                }
                if (corners[0] > corners[3]) {
                    aux = corners[0]; corners[0] = corners[3]; corners[3] = aux;
                    aux = corners[2]; corners[2] = corners[1]; corners[1] = aux;
                }
                if (corners[1] > corners[2]) {
                    aux = corners[1]; corners[1] = corners[2]; corners[2] = aux;
                }
                double h = 1./(corners[2] - corners[0]);
                double dx = min(0, corners[1]) - corners[0];
                double tol = 1e-12;
                (*vf_for_2dcell)[j][i] = 0.5*dx*dx/(corners[1] - corners[0] + tol)*h;
                if (corners[1] < tol) {
                    dx = min(0., corners[2]) - corners[1];
                    (*vf_for_2dcell)[j][i] += dx*h;
                }
                if (corners[2] < tol) {
                    dx = corners[3];
                    (*vf_for_2dcell)[j][i] = 1. - .5*dx*dx/(corners[3] - corners[2] + tol)*h;
                }

                ++(*nmixcell); // don't forget to count the mixcells
            }
        }
    }

    if ((*nmixcell) == 0) {
        printf ("ERROR in record_test__cal_mixcell_zgrad: nmixcell should never be zero\n");
        exit(0);
    }
    // allocate the mixcell structures
    *nmat_in_mixcell = (int*) malloc((*nmixcell)*sizeof(int));
    *matids_in_mixcell = (int**) malloc((*nmixcell)*sizeof(int*));
    *vf_in_mixcell = (double**) malloc((*nmixcell)*sizeof(double*));
    *ijk_in_mixcell = (int**)malloc((*nmixcell)*sizeof(int*));
    *vfgrad_in_mixcell = (double***) malloc((*nmixcell)*sizeof(double**));

    // there are only 2 materials in this test, so lengths of ragged arrays are known
    for (mx = 0; mx < (*nmixcell); mx++) {
        if (mx == 0) {
            (*matids_in_mixcell)[0] = (int*) malloc(2*(*nmixcell)*sizeof(int));
            (*vf_in_mixcell)[0] = (double*) malloc(2*(*nmixcell)*sizeof(double));
            (*ijk_in_mixcell)[0] = (int*)malloc(2*(*nmixcell)*sizeof(int));
            (*vfgrad_in_mixcell)[0] = (double**) malloc(2*(*nmixcell)*sizeof(double*));

            (*vfgrad_in_mixcell)[0][0] = (double*) malloc(2*2*(*nmixcell)*sizeof(double));
            (*vfgrad_in_mixcell)[0][1] = (*vfgrad_in_mixcell)[0][0] + 2;
            offset = 4;
        }
        else {
            (*matids_in_mixcell)[mx] = (*matids_in_mixcell)[mx-1] + 2;
            (*vf_in_mixcell)[mx] = (*vf_in_mixcell)[mx-1] + 2;
            (*ijk_in_mixcell)[mx] = (*ijk_in_mixcell)[mx-1] + 2;
            (*vfgrad_in_mixcell)[mx] = (*vfgrad_in_mixcell)[mx-1] + 2;

            (*vfgrad_in_mixcell)[mx][0] = (*vfgrad_in_mixcell)[0][0] + offset;
            offset += 2;
            (*vfgrad_in_mixcell)[mx][1] = (*vfgrad_in_mixcell)[0][0] + offset;
            offset += 2;
        }
    }
    memset((*nmat_in_mixcell), 0, (*nmixcell)*sizeof(int));
    memset((*matids_in_mixcell)[0], 0, 2*(*nmixcell)*sizeof(int));
    memset((*ijk_in_mixcell)[0], 0, 2*(*nmixcell)*sizeof(int));
    memset((*vf_in_mixcell)[0], 0., 2*(*nmixcell)*sizeof(double));
    memset((*vfgrad_in_mixcell)[0][0], 0., 2*2*(*nmixcell)*sizeof(double));

    // assign the values in mixcell
    mx = 0;
    for (j = 0; j < ncell_ext[1]; j++) {
        for (i = 0; i < ncell_ext[0]; i++) {
            double vf = (*vf_for_2dcell)[j][i];
            if (0. < vf && vf < 1.) {
                (*nmat_in_mixcell)[mx] = 2;

                (*matids_in_mixcell)[mx][0] = 0;
                (*matids_in_mixcell)[mx][1] = 1;

                (*ijk_in_mixcell)[mx][0] = i; // WARNING: these two can be flipped!
                (*ijk_in_mixcell)[mx][1] = j;

                (*vf_in_mixcell)[mx][0] = vf;
                (*vf_in_mixcell)[mx][1] = 1. - vf;

                (*mixcell_for_2dcell)[j][i] = mx;
                ++mx;
            }
            else {
                (*mixcell_for_2dcell)[j][i] = -1;
            }
        }
    }
    free(cn[0]);    free(cn);

} // generate_two_material_mesh


void record_test__cal_mixcell_zgrad2d(int num_tests, const char *test_h5, const int verbose) {

    int dim, ncell[2], nbdry, ncell_ext[2], nnode[2], nnode_ext[2];
    double dx[2], xl[2];

    char group_name[50];
    hid_t attr_id, test_gid;

    int **nmat_for_2dcell, **matid_for_2dcell, **mixcell_for_2dcell;
    double **vf_for_2dcell;

    int nmixcell, *nmat_in_mixcell, **ijk_in_mixcell, **matids_in_mixcell;
    double **vf_in_mixcell, ***vfgrad_in_mixcell;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    for (int t = 0; t < num_tests; t++) {

        generate_two_material_mesh(&dim, &nbdry, ncell, nnode, ncell_ext, nnode_ext, dx, xl,
            &nmat_for_2dcell, &matid_for_2dcell, &mixcell_for_2dcell, &vf_for_2dcell,
            &nmixcell, &nmat_in_mixcell, &matids_in_mixcell, &vf_in_mixcell, &ijk_in_mixcell,
            &vfgrad_in_mixcell);

        // Record the input data
        sprintf(group_name, "/test_%d", t);
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_1d_int(test_gid, "ncell", ncell, 2);
        h5_write_1d_int(test_gid, "nbdry", &nbdry, 1);
        h5_write_1d(test_gid, "xl", xl, 2);
        h5_write_1d(test_gid, "dx", dx, 2);
        h5_write_2d_int(test_gid, "nmat_for_2dcell", nmat_for_2dcell, ncell_ext);
        h5_write_2d_int(test_gid, "matid_for_2dcell", matid_for_2dcell, ncell_ext);
        h5_write_2d_int(test_gid, "mixcell_for_2dcell", mixcell_for_2dcell, ncell_ext);
        //h5_write_2d(test_gid, "vf_for_2dcell", vf_for_2dcell, ncell_ext);

        h5_write_1d_int(test_gid, "nmixcell", &nmixcell, 1);
        h5_write_1d_int(test_gid, "nmat_in_mixcell", nmat_in_mixcell, nmixcell);

        int sizes[2] = {nmixcell, 2};
        h5_write_2d_int(test_gid, "ijk_in_mixcell", ijk_in_mixcell, sizes);
        h5_write_2d(test_gid, "vf_in_mixcell", vf_in_mixcell, sizes);
        h5_write_2d_int(test_gid, "matids_in_mixcell", matids_in_mixcell, sizes);

        // Call the function under test
        cal_mixcell_zgrad2d(ncell, nbdry, xl, dx,
                            nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell,
                            nmixcell, nmat_in_mixcell, matids_in_mixcell, vf_in_mixcell,
                            ijk_in_mixcell, vfgrad_in_mixcell);

        // Record function output
        h5_write_1d(test_gid, "vfgrad_in_mixcell", vfgrad_in_mixcell[0][0], 4*nmixcell);

        H5Gclose(test_gid);
        free(nmat_for_2dcell[0]);    free(nmat_for_2dcell);
        free(matid_for_2dcell[0]); free(matid_for_2dcell);
        free(mixcell_for_2dcell[0]);    free(mixcell_for_2dcell);
        free(vf_for_2dcell[0]);    free(vf_for_2dcell);

        if (nmixcell > 0) {
            free(nmat_in_mixcell);
            free(ijk_in_mixcell[0]); free(ijk_in_mixcell);
            free(vf_in_mixcell[0]); free(vf_in_mixcell);
            free(matids_in_mixcell[0]); free(matids_in_mixcell);
            free(vfgrad_in_mixcell[0][0]); free(vfgrad_in_mixcell[0]); free(vfgrad_in_mixcell);
        }

    }

    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);

} // record_test__cal_mixcell_zgrad2d


void record_test__get_mpoly(int num_tests, const char *test_h5, const int verbose) {

    int dim, ncell[2], nbdry, ncell_ext[2], nnode[2], nnode_ext[2];
    double dx[2], xl[2];

    char group_name[50];
    hid_t attr_id, test_gid;

    int **nmat_for_2dcell, **matid_for_2dcell, **mixcell_for_2dcell;
    double **vf_for_2dcell;

    int nmixcell, *nmat_in_mixcell, **ijk_in_mixcell, **matids_in_mixcell;
    double **vf_in_mixcell, ***vfgrad_in_mixcell;

    delete_file_if_exists(test_h5);
    hid_t testfile_id = h5_new_file(test_h5);

    for (int t = 0; t < num_tests; t++) {

        generate_two_material_mesh(&dim, &nbdry, ncell, nnode, ncell_ext, nnode_ext, dx, xl,
            &nmat_for_2dcell, &matid_for_2dcell, &mixcell_for_2dcell, &vf_for_2dcell,
            &nmixcell, &nmat_in_mixcell, &matids_in_mixcell, &vf_in_mixcell, &ijk_in_mixcell,
            &vfgrad_in_mixcell);

        // Record the input data
        sprintf(group_name, "/test_%d", t);
        test_gid = h5_open_group(testfile_id, group_name);
        h5_write_1d_int(test_gid, "ncell", ncell, 2);
        h5_write_1d_int(test_gid, "nbdry", &nbdry, 1);
        h5_write_1d(test_gid, "xl", xl, 2);
        h5_write_1d(test_gid, "dx", dx, 2);
        h5_write_2d_int(test_gid, "nmat_for_2dcell", nmat_for_2dcell, ncell_ext);
        h5_write_2d_int(test_gid, "matid_for_2dcell", matid_for_2dcell, ncell_ext);
        h5_write_2d_int(test_gid, "mixcell_for_2dcell", mixcell_for_2dcell, ncell_ext);
        //h5_write_2d(test_gid, "vf_for_2dcell", vf_for_2dcell, ncell_ext);

        h5_write_1d_int(test_gid, "nmixcell", &nmixcell, 1);
        h5_write_1d_int(test_gid, "nmat_in_mixcell", nmat_in_mixcell, nmixcell);

        int sizes[2] = {nmixcell, 2};
        h5_write_2d_int(test_gid, "ijk_in_mixcell", ijk_in_mixcell, sizes);
        h5_write_2d(test_gid, "vf_in_mixcell", vf_in_mixcell, sizes);
        h5_write_2d_int(test_gid, "matids_in_mixcell", matids_in_mixcell, sizes);

        // Call the function under test
        mix_data_pass_mat(nmixcell, nmixcell, nmixcell, nmat_in_mixcell, ijk_in_mixcell,
                          matids_in_mixcell, vf_in_mixcell, NULL, NULL, NULL);
        get_mpoly(dim, ncell, nbdry, xl, dx, nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell,
                  NULL, NULL, NULL);

        // Retrieve the data from mpoly structures
        double **coords_for_mixcell;
        int _nmixcell, **_ijk_for_mixcell, *_nmat_for_mixcell,
            **_matids_for_mixcell, *nnode_for_mixcell,
            **nnode_for_mpoly_in_mixcell, ***nodes_for_mpoly_in_mixcell,
            **nface_for_mpoly_for_mixcell, ***nnode_for_face_ea_mpoly_for_mixcell,
            ***nodelist_for_face_ea_mpoly_for_mixcell,
            *mix_to_mpoly_mix, *mix_mpoly_to_mix, nmixcell_mpoly,
            **nnode_for_minterface_in_mixcell, ***nodes_for_minterface_in_mixcell;
        mat_pass_mix(&_nmixcell, &_ijk_for_mixcell, &_nmat_for_mixcell,
                  &_matids_for_mixcell, &nnode_for_mixcell,
                  &coords_for_mixcell, &nnode_for_mpoly_in_mixcell, &nodes_for_mpoly_in_mixcell,
                  &nface_for_mpoly_for_mixcell, &nnode_for_face_ea_mpoly_for_mixcell,
                  &nodelist_for_face_ea_mpoly_for_mixcell);

        mat_pass_mix_extra(&nmixcell_mpoly, &mix_to_mpoly_mix, &mix_mpoly_to_mix,
            &nnode_for_minterface_in_mixcell, &nodes_for_minterface_in_mixcell);

        // Record function output
        h5_write_1d_int(test_gid, "nnode_for_mixcell", nnode_for_mixcell, nmixcell_mpoly);
        h5_write_1d_int(test_gid, "nmixcell_mpoly", &nmixcell_mpoly, 1);
        h5_write_1d_int(test_gid, "mix_to_mpoly_mix", mix_to_mpoly_mix, nmixcell);
        h5_write_1d_int(test_gid, "mix_mpoly_to_mix", mix_mpoly_to_mix, nmixcell_mpoly);

        int nm, mx, mxp, nmat_total, nmix_total, nnod_total, nint_total, k, nn;
        nmat_total = 0;
        nmix_total = 0;
        nnod_total = 0;
        nint_total = 0;
        for (mxp = 0; mxp < nmixcell_mpoly; mxp++) {
            mx = mix_mpoly_to_mix[mxp];
            nmat_total += nmat_in_mixcell[mx];
            nnod_total += nnode_for_mixcell[mxp];
            nint_total += nmat_in_mixcell[mx] - 1;
            for (nm = 0; nm < nmat_in_mixcell[mx]; nm++) {
                nmix_total += nnode_for_mpoly_in_mixcell[mxp][nm];
            }
        }
        if (verbose) {
            printf ("next test:\n");
            printf (" - total number of nodes in all mixcells: %d\n", nnod_total);
            printf (" - total number of nodes in all mixcell material polygons: %d\n", nmix_total);
            printf (" - total nodes in each mixcell: {");
            for (mxp = 0; mxp < nmixcell_mpoly; mxp++)
                printf ("%d ", nnode_for_mixcell[mxp]);
            printf ("}\n");

            printf (" - total polygons (= total materials in mixcell with mpoly) = %d\n", nmat_total);
            for (mxp = 0; mxp < nmixcell_mpoly; mxp++) {
                mx = mix_mpoly_to_mix[mxp];
                printf ("   [");
                for (nm = 0; nm < nmat_in_mixcell[mx]; nm++) {
                    printf("%d ", nnode_for_mpoly_in_mixcell[mxp][nm]);
                }
                printf ("] ");
            }
            printf("\n");

            printf (" - total interfaces = %d\n", nint_total);
            for (mxp = 0; mxp < nmixcell_mpoly; mxp++) {
                mx = mix_mpoly_to_mix[mxp];
                printf ("   [");
                for (nm = 0; nm < nmat_in_mixcell[mx]-1; nm++) {
                    printf("%d ", nnode_for_minterface_in_mixcell[mxp][nm]);
                }
                printf ("] ");
            }
            printf("\n");
        }

        int *nnode_for_mpoly_in_mixcell_flat = (int*)malloc(nmat_total*sizeof(int));
        k = 0;
        for (mxp = 0; mxp < nmixcell_mpoly; mxp++) {
            mx = mix_mpoly_to_mix[mxp];
            for (nm = 0; nm < nmat_in_mixcell[mx]; nm++) {
                nnode_for_mpoly_in_mixcell_flat[k++] = nnode_for_mpoly_in_mixcell[mxp][nm];
            }
        }
        h5_write_1d_int(test_gid, "nnode_for_mpoly_in_mixcell",
                                   nnode_for_mpoly_in_mixcell_flat, nmat_total);
        free(nnode_for_mpoly_in_mixcell_flat);
        h5_write_1d_int(test_gid, "nnode_for_mixcell_total", &nnod_total, 1);
        h5_write_1d_int(test_gid, "nnode_for_mix_mpoly_total", &nmix_total, 1);


        int *nodes_for_mpoly_in_mixcell_flat = (int*)malloc(nmix_total*sizeof(int));
        k = 0;
        for (mxp = 0; mxp < nmixcell_mpoly; mxp++) {
            mx = mix_mpoly_to_mix[mxp];
            for (nm = 0; nm < nmat_in_mixcell[mx]; nm++) {
                for (nn = 0; nn < nnode_for_mpoly_in_mixcell[mxp][nm]; nn++)
                    nodes_for_mpoly_in_mixcell_flat[k++] = nodes_for_mpoly_in_mixcell[mxp][nm][nn];
            }
        }
        h5_write_1d_int(test_gid, "nodes_for_mpoly_in_mixcell",
                                   nodes_for_mpoly_in_mixcell_flat, nmix_total);
        free(nodes_for_mpoly_in_mixcell_flat);


        int *nodes_for_minterface_in_mixcell_flat = (int*)malloc(2*nint_total*sizeof(int));
        k = 0;
        for (mxp = 0; mxp < nmixcell_mpoly; mxp++) {
            mx = mix_mpoly_to_mix[mxp];
            for (nm = 0; nm < nmat_in_mixcell[mx]-1; nm++) {
                for (nn = 0; nn < 2; nn++)
                    nodes_for_minterface_in_mixcell_flat[k++] =
                        nodes_for_minterface_in_mixcell[mxp][nm][nn];
            }
        }
        h5_write_1d_int(test_gid, "nodes_for_minterface_in_mixcell",
                                   nodes_for_minterface_in_mixcell_flat, 2*nint_total);
        free(nodes_for_minterface_in_mixcell_flat);


        double *coords_for_mixcell_flat = (double*)malloc(dim*nnod_total*sizeof(double));
        k = 0;
        for (mxp = 0; mxp < nmixcell_mpoly; mxp++) {
            mx = mix_mpoly_to_mix[mxp];
            for (nn = 0; nn < nnode_for_mixcell[mxp]; nn++) {
                coords_for_mixcell_flat[k++] = coords_for_mixcell[mxp][2*nn];
                coords_for_mixcell_flat[k++] = coords_for_mixcell[mxp][2*nn+1];
            }
        }
        h5_write_1d(test_gid, "coords_for_mixcell", coords_for_mixcell_flat, dim*nnod_total);
        free(coords_for_mixcell_flat);

        //h5_write_1d(test_gid, "vfgrad_in_mixcell", vfgrad_in_mixcell[0][0], 4*nmixcell);

        H5Gclose(test_gid);
        free(nmat_for_2dcell[0]);    free(nmat_for_2dcell);
        free(matid_for_2dcell[0]);   free(matid_for_2dcell);
        free(mixcell_for_2dcell[0]); free(mixcell_for_2dcell);
        free(vf_for_2dcell[0]);      free(vf_for_2dcell);

        free(nmat_in_mixcell);
        free(ijk_in_mixcell[0]); free(ijk_in_mixcell);
        free(vf_in_mixcell[0]); free(vf_in_mixcell);
        free(matids_in_mixcell[0]); free(matids_in_mixcell);
        free(vfgrad_in_mixcell[0][0]); free(vfgrad_in_mixcell[0]); free(vfgrad_in_mixcell);

        clean_mpoly();

    }

    H5Fclose(testfile_id);
    printf("written: %s\n", test_h5);

} // record_test__get_mpoly


int main (int argc, char *argv[]) {

    const int verbose = 0;
    srand(time(NULL));  // Seed random number generator

    /*** testing some tests3d.c functions:
    */
    const int dims[3] = {3, 5, 2};
    double ****vv;
    allocate_velocity_3dmesh(dims, &vv);

    for (int m = 0; m < 3; m++) {
        for (int k = 0; k < dims[2]; k++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int i = 0; i < dims[0]; i++) {
                    printf ("%4.1f ", vv[m][k][j][i]);
                }
                printf("\n");
            }
            printf("----------\n");
        }
        printf("==========\n");
    }
    exit(0);
    /***/

    //--- simple unit tests ----------------------------------------------------------
    record_test__rec_rec(100,           "testdata/test_rec_rec.h5",      verbose);
    record_test__sound_speed_solid(243, "testdata/sound_speed_solid.h5", verbose);
    record_test__find_interface2d(120,  "testdata/find_interface2d.h5",  verbose);
    record_test__cal_poly_area(37,      "testdata/cal_poly_area.h5",     verbose);
    record_test__rz_area(123,           "testdata/rz_area.h5",           verbose);
    record_test__bounds_2d(150,         "testdata/bounds_2d.h5",         verbose);
    record_test__order_nodes_along_norm(100, "testdata/order_nodes_along_norm.h5", verbose);
    record_test__cal_distance2d(100,    "testdata/cal_distance2d.h5",    verbose);
    record_test__cal_cell_zgrad2d(100,  "testdata/cal_cell_zgrad2d.h5",  verbose);
    record_test__reconstruct2d_nmat_pagosa(100,"testdata/reconstruct2d_nmat_pagosa.h5", verbose);

    //--- tests on two-material semirandom meshes  -----------------------------------
    record_test__cal_mixcell_zgrad2d(10,"testdata/cal_mixcell_zgrad2d.h5", verbose);
    record_test__get_mpoly(10,          "testdata/get_mpoly.h5",         verbose);

    //--- tests on random meshes  ----------------------------------------------------
    const char meshes_h5[] = "testdata/test_meshes.h5";
    generate_random_meshes(2, 10, meshes_h5, verbose);
    record_test__bdry_cell_2d(meshes_h5, "testdata/bdry_cell_2d.h5", verbose);
    record_test__bdry_node_2d(meshes_h5, "testdata/bdry_node_2d.h5", verbose);
    record_test__bdry_cell_1var_2d(meshes_h5, "testdata/bdry_cell_1var_2d.h5", verbose);
    record_test__bdry_cell_ragged_2d(meshes_h5, "testdata/bdry_cell_ragged_2d.h5", verbose);
    record_test__compute_divu(meshes_h5, "testdata/compute_divu.h5", verbose);
    record_test__compute_qvis(meshes_h5, "testdata/compute_qvis.h5", verbose);
    record_test__compute_force(meshes_h5, "testdata/compute_force.h5", verbose);
    record_test__update_vel_comp(meshes_h5, "testdata/update_vel_comp.h5", verbose);
    record_test__update_density(meshes_h5, "testdata/update_density.h5", verbose);
    record_test__update_energy(meshes_h5, "testdata/update_energy.h5", verbose);
    record_test__courant_from_cs(meshes_h5, "testdata/courant_from_cs.h5", verbose);
    record_test__mesh_sound_speed(meshes_h5, "testdata/mesh_sound_speed.h5", verbose);
    record_test__mat_sound_speed(meshes_h5, "testdata/mat_sound_speed.h5", verbose);

    return 0;
}
