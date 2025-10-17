#ifndef _VOF2_
#define _VOF2_

#ifdef __cplusplus
extern "C" {
#endif

void reconstruct2d_nmat_pagosa(int geop,
                        double *xl, double *dx,
                        int **nmat_mesh, int ***matid_mesh, double ***vf_mesh,
                        int *nnode_final, double *coords_final,
                        int *nnode_for_interface, int **nodes_for_interface,
                        int *nnode_for_mpoly, int **nodelist_for_mpoly);

void find_interface2d(int nnode, double *coords, int *nodelist,
                      double *norm, double distance,
                      int *node_loc,
                      int *nnode_new, double *coords_new,
                      int *nnode_interface, int *nodes_interface,
                      int *nnode_lower, int **nodelist_lower,
                      int *nnode_upper, int **nodelist_upper);

void remap2d_scaled(int ifdump,
                    int nnode, double *coords, double *inward_norm_ea_face,
                    int nnode_m, double *coords_m,
                    int *nnodex, double **coordsx);

void bounds_2d(int geop,
               double vf_to_match, double volume, double *norm,
               int nnode, double *coords, int *nodelist,
               int *node_order_for_ds, double *ds_ea_node,
               double *ds_lower, double *ds_upper,
               double *vf_lower, double *vf_upper,
               int *nnode_new,   double *coords_new,
               int *vol_matched,
               int *nnode_interface, int *nodes_interface,
               int *nnode_lower, int **nodelist_lower,
               int *nnode_upper, int **nodelist_upper);

void cal_distance2d(int geop,
                    double vf_to_match, double volume, double *norm,
                    int nnode, double *coords, int *nodelist,
                    int *nnode_new, double *coords_new,
                    double *distance,
                    int *nnode_interface, int *nodes_interface,
                    int *nnode_lower, int **nodelist_lower,
                    int *nnode_upper, int **nodelist_upper);

#ifdef __cplusplus
}
#endif
#endif

