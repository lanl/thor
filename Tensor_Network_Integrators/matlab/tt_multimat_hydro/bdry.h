#ifndef _BDRY_
#define _BDRY_

#ifdef __cplusplus
extern "C" {
#endif

#include "globals.h"
#include "minip.h" 

void bdry_cell_2d(int nmat, int *ncell, int nbdry,
                    Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                    double ***vf_2dmat, double ***rho_2dmat,
                    double ***ei_2dmat, double ***pres_2dmat,
                    double **rho_2dcell, double **ei_2dcell, double **pres_2dcell);

void bdry_cell_3d(int nmat, int *ncell, int nbdry,
                Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                double ****vf_3dmat, double ****rho_3dmat,
                double ****ei_3dmat, double ****pres_3dmat, 
                double ***rho_3dcell, double ***ei_3dcell, double ***pres_3dcell);

void bdry_node_2d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                  double ***vel_for_2dnode);

void bdry_node_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                  double ****vel_for_3dnode);

void bdry_cell_vel_2d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                      double ***vel_for_2dcell);

void bdry_cell_vel_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                      double ****vel_for_3dcell);

#ifdef __cplusplus
}
#endif
#endif
