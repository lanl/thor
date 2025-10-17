
#ifdef MPI
#include "mpi.h"
#endif

#include <sys/stat.h> 
#include "minip.h"
#include "globals.h" 

#include "init.h"
#include "control.h" 

int main(int argc, char **argv) {

    struct stat stats;
    char   *filename, *probname;
    int    iret, dim, nbdry, nmat, ncell[3], err;
    int    *solid_num_ea_mat, *if_fixed_state_ea_mat, *matid_ea_mat;
    double *rho_fixed_state_ea_mat, *ei_fixed_state_ea_mat, *pres_fixed_state_ea_mat;
    double xl_prob[3], xr_prob[3];
    double *gamma_ea_mat;
    double tmax = 10000.0;
    int ncycle_max = 1;
    double dt_initial = 0.001;
    double courant = 0.25;

    int mesh_scaling = 0;

    int    ncycle_viz = 1;
    double dt_viz     = 1.0e+06;

    Bdry_Type btype_lower[3], btype_upper[3];

    if (argc < 2) {
        printf("usage: xtest inputFileName\n");
        return 0;
    }
    filename = argv[1];

//! Check to get mesh_scaling
    if (argc >= 3)
    {
        mesh_scaling = atoi(argv[2]);
        printf("mesh_scaling = %d\n", mesh_scaling);
    }

    printf("Checking file: '%s'\n", filename);
    
    iret = stat(filename, &stats);
    // assert(!iret);
    
    solid_num_ea_mat = NULL;
    if_fixed_state_ea_mat = NULL;
    matid_ea_mat = NULL;
    gamma_ea_mat  = NULL;
    probname = NULL;

    rho_fixed_state_ea_mat = NULL;
    ei_fixed_state_ea_mat  = NULL;
    pres_fixed_state_ea_mat = NULL;

    err = init(filename,
               &dim, ncell, &nbdry, btype_lower, btype_upper, xl_prob, xr_prob,
               &nmat, &matid_ea_mat,  &solid_num_ea_mat, &gamma_ea_mat,
               &if_fixed_state_ea_mat, &rho_fixed_state_ea_mat,
               &ei_fixed_state_ea_mat, &pres_fixed_state_ea_mat,
               &tmax, &ncycle_max,
               &dt_viz, &ncycle_viz,
               &courant, &dt_initial,  &probname, mesh_scaling);

    if (err) return 0;

    int ncycle_final = 0;
    
    control(probname, dim, ncell, nbdry, btype_lower, btype_upper, xl_prob, xr_prob,
	    nmat, matid_ea_mat, solid_num_ea_mat, gamma_ea_mat,
            if_fixed_state_ea_mat, 
            rho_fixed_state_ea_mat, ei_fixed_state_ea_mat, pres_fixed_state_ea_mat,
            dt_initial, tmax, ncycle_max, courant, ncycle_viz, dt_viz,
            &ncycle_final, mesh_scaling);  
    
    return 0;
} 
