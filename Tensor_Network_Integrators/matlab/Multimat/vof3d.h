/* ======================= copyright begin ========================  */
/* Copyright (C) 2010 Los Alamos National Laboratory.                */
/* All rights Reserved.                                              */
/* Export Controlled Information                                     */
/* ======================== copyright end =========================  */

#ifndef _VOF3D_
#define _VOF3D_

#ifdef __cplusplus
extern "C" {
#endif

void reconstruct3d_nmat_pagosa(double *xl, double *dx,
                        int ***nmat_mesh, int ****matid_mesh, double ****vf_mesh,
                        int *nnode_final, double **coords_final,
                        int *nnode_for_minterface, int ***nodes_for_minterface,
                        int *nface_for_mat, int ***nnode_for_face_ea_mat,
                        int ***nodelist_for_face_ea_mat);


void polyhedron_plane(int nface, int nnode, double *coords,
                      int *nnode_ea_face, int *nodelist_for_face,
                      double *norm_plane, double ds_plane,

          int *nface_lower, int *nnode_lower, double **coords_lower,
          int **nnode_ea_face_lower, int **nodelist_for_face_lower,
          int *nface_upper, int *nnode_upper, double **coords_upper,
          int **nnode_ea_face_upper, int **nodelist_for_face_upper);

void cal_vol(int nface, int nnode, double *coords,
             int *nnode_for_face, int *nodelist_for_face,
             double *vol);

void get_edges(int nface, int nnode, int *nnode_for_face,
                int *nodelist_for_face, int *nedge,
                int **edgelist_for_face, int **nodelist_for_edge);

///////////////////////////////////////////////////////////////////////////////////////
void write_poly_3d(int fileid, char *meshname, int nface, int nnode, double *coords,
                   int *nnode_for_face, int *nodelist_for_face);

void write_poly_2d(int fileid,  char *meshname, int nnode_in_coord, double *coords,
                   int nnode, int *nodelist);

#ifdef __cplusplus
}
#endif
#endif

