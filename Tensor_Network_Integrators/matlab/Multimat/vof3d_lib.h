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

void testf(double x, double *xp, double **xpp, int n, int *np, int **npp);

void matlab_interface3d(double dx_scale,
                 int nface, int nnode, double *coords,
                 int *nnode_for_face, int *nodelist_for_face,
                 double cell_vol, double vf_to_match, double *normal,
                 int *out_nnode_new,
                 double *out_coords_new,
                 int *out_nnode_for_interface,
                 int *out_nodelist_for_interface,
                 int *out_faceindexlow_of_interface,
                 int *out_faceindexhgh_of_interface,
                 int *out_nface_lower,
                 int *out_nnode_for_face_lower,
                 int *out_nodelist_for_face_lower,
                 int *out_nface_upper,
                 int *out_nnode_for_face_upper,
                 int *out_nodelist_for_face_upper);

void matlab_polyhedron_plane(int nface, int nnode, double *coords,
                      int *nnode_ea_face, int *nodelist_for_face,
                      double *norm_plane, double ds_plane,
          int *out_nface_lower, int *out_nnode_lower, double *out_coords_lower,
          int *out_nnode_ea_face_lower, int *out_nodelist_for_face_lower,
          int *out_nface_upper, int *out_nnode_upper, double *out_coords_upper,
          int *out_nnode_ea_face_upper, int *out_nodelist_for_face_upper);
#ifdef __cplusplus
}
#endif
#endif

