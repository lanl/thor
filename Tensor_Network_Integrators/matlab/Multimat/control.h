#ifndef _CONTROL_
#define _CONTROL_

#ifdef __cplusplus
extern "C" {
#endif

#include "globals.h"
#include "minip.h" 

void control(const char *probname,
             int dim, int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
             double *xl_prob, double *xr_prob,
             int nmat, int *matids, int *solid_num_ea_mat, double *gamma_ea_mat,
             int *if_fixed_state_ea_mat,
             double *rho_fixed_state_ea_mat, double *ei_fixed_state_ea_mat, double *pres_fixed_state_ea_mat,
             double dt_initial, double tmax, int ncycle_max, double courant,
             int ncycle_viz_freq, double dt_viz_freq, int *ncycle_final, int mesh_scaling); 

#ifdef __cplusplus
}
#endif
#endif
