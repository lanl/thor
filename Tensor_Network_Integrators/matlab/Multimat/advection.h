#ifndef _ADVECTION_
#define _ADVECTION_

#ifdef __cplusplus
extern "C" {
#endif

#include "globals.h"
#include "minip.h" 

void advection(int fileid, int dim, int *ncell, int nbdry, double *xl_prob, double *dx_prob,
              int nmat,  int *solid_num_ea_mat, double *gamma_ea_mat,
              int *if_fixed_state_ea_mat,
              double *rho_fixed_state_ea_mat, double *ei_fixed_state_ea_mat, double *pres_fixed_state_ea_mat,
              int ncycle, int dir, double dt,
              Bdry_Type *btype_lower, Bdry_Type *btype_upper,
              double ***vf_2dmat, double ***rho_2dmat, double ***ei_2dmat, double ***pres_2dmat,
              double **rho_2dcell, double **ei_2dcell, double **pres_2dcell,
              double ***vel_2dnode, double ***vav_for_2dnode,
//
              double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat, double ****pres_3dmat,
              double ***rho_3dcell, double ***ei_3dcell, double ***pres_3dcell,
              double ****vel_3dnode, double ****vav_for_3dnode,
              double *courant_adv);
// #define H5DUMP_DEBUG
// #define ASCII_OUTPUT
#ifdef __cplusplus
}
#endif
#endif
