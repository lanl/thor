/* ======================= copyright begin ========================  */
/* Copyright (C) 2010 Los Alamos National Security, LLC.             */
/* All rights Reserved.  See Copyright Notice File.                  */
/* Export Controlled Information                                     */
/*                                                                   */
/* Author:  William W. Dai, dai@lanl.gov                             */
/* ======================== copyright end =========================  */
 
#include "globals.h"
#include "geom.h"
#include "util.h" 
#include "vof3d.h" 
#include "geom.h" 
// #include "mio.h"

enum amr_Loc {
 
    // not bdry
 
     not_located = -1,
     inner       = 0,
 
     // faces
 
     xlow = 1,
     xhgh = 2,
     ylow = 3,
     yhgh = 4,
     zlow = 5,
     zhgh = 6,
 
     // 3D coarners
 
     corner_x0y0z0 = 7,
     corner_x1y0z0 = 8,
     corner_x0y1z0 = 9,
     corner_x1y1z0 = 10,
 
     corner_x0y0z1 = 11,
     corner_x1y0z1 = 12,
     corner_x0y1z1 = 13,
     corner_x1y1z1 = 14,
 
     // edges
 
     zedge_x0y0   = 15,
     zedge_x1y0   = 16,
     zedge_x0y1   = 17,
     zedge_x1y1   = 18,
 
     yedge_z0x0   = 19,
     yedge_z1x0   = 20,
     yedge_z0x1   = 21,
     yedge_z1x1   = 22,
 
     xedge_y0z0   = 23,
     xedge_y1z0   = 24,
     xedge_y0z1   = 25,
     xedge_y1z1   = 26
};
typedef enum amr_Loc amr_Loc;
 
enum geom_obj_type {
     polygon = 1,
     sphere  = 2
};
typedef enum geom_obj_type geom_obj_type;

 
static double tiny   = 1.0e-30;
static double  small = 1.0e-12;
static double  tol   = 1.0e-10;
 
static int if_decompose_concave = 0;
 
static int if_write_debug_viz = 0;
 
// reserved for OSO file
 
static int  nfile_oso = 0;
static int  reg_oso[100];
static int  nxy_oso[100];
static double *xys_oso[100];
static FILE *fp_oso[100];
 
static int    npart_oso[100];
static int    *nnode_part_oso[100];
static double *factor_part_oso[100];
static double **coords_part_oso[100];
 
static int     *nnode_grp_oso[100];
static int     *nedge_grp_oso[100];
static int     *nface_grp_oso[100];
static double  *ctr_grp_oso[100];
static double  **coords_grp_oso[100];
 
static int     **nodelist_for_edge_oso[100];
static int     **edgelist_for_face_oso[100];
static int     **nodelist_for_face_oso[100];
static int     **nnode_for_face_oso[100];
static double  **normal_of_face_oso[100];
 
 
// static int fid_for_viz_oso = -1;
 
 
 
// reserved for obligue ellipse
 
static int    aligned_sav    = 1;
static int    major_axis_sav = -1;
static double sinphi_sav     = 0.0;
static double cosphi_sav     = 1.0;
static double sintheta_sav   = 0.0;
static double costheta_sav   = 1.0;
 
//  reserved for conic
static double sinphi_conic    = 0.0;
static double cosphi_conic    = 1.0;
static double sintheta_conic  = 0.0;
static double costheta_conic  = 1.0;
 
////////////////////////////////////////////////////////////////////////////////////
void hex_plane(double *point, double *normal, double *vertices, double *coord_res, int *nnode);
////////////////////////////////////////////////////////////////////////////////////
 
void check_align(int dim, double *pt1, double *pt2, double *rads);
 
void get_direction(int dim, double *pt1, double *pt2);
 
void check_plane(int dim, double *xyz1, double *xyz2, int nn0, int *ifcyl,
                 int *nnode, int *nedge, int *nface,
                 int *nodelist_for_edge, int *edgelist_for_face,
                 int *nedge_ea_face, int *nodelist_for_face,
                 double *coords, double *ctr, double *normal_of_face);
 
void gellipse_rec(int ifinquiry, int *mixed, int geop, double vcell, int dim,
                  double *ctr, double *r, double *xl, double *dx, double *vol);
 
void gellipse_rec_old(int ifinquiry, int *mixed, int geop, double vcell, int dim,
                  double *ctr, double *r, double *xl, double *dx, double *vol);
 
void gconic_rec_rotate(int ifinquiry, int *mixed,
                int dim, double r0, double r1, double *ctr0, double *ctr1,
                double *xl, double *dx, double *vol);
 
void gconic_rec_cyl(int ifinquiry, int *mixed, int ifcyl,
                int dim, double r0, double r1, double *cl, double *cr,
                double *xl, double *dx, double *vol);
 
void conic_rec(int dim, double r0, double r1, double *cl, double *cr,
             double *xl, double *dx, double *vol);
 
void conic_rec0(int dim, double r, double *cl, double *cr,
             double *xl, double *dx, double *vol);
 
void gpoly3d_rec(int nreg, int ifinquiry, int *mixed, int ifcyl, int dim,
                 double *coords, int nnode,
                 int nedge, int nface, int *nodelist_for_edge,
                 int *edgelist_for_face, int *nedge_ea_face,
                 int *nodelist_for_face, int *num_nodes_for_face,
                 double *ctr, double *norm_of_face,
                 double *xl, double *dx, double *vol);
 
void poly3dg_rec(int ifinquiry, int *mixed, int ifcyl, int dim,
                 double *coords, int nnode,
                 int nedge, int nface, int *nodelist_for_edge,
                 int *edgelist_for_face, int *nedge_ea_face,
                 int *nodelist_for_face, int *num_nodes_for_face,
                 double *ctr, double *norm_of_face,
                 double *xl, double *dx, double *vol);
 
void poly3d_rec(int ifinquiry, int *mixed, int ifcyl, int dim,
                 int nnode, double *xys1,
                 double *xys2,
                 double z0, double z1, double *xl, double *dx, double *vol);
 
void poly_rec3d(int dim, int nnode, double *xys1,
                double *xys2,
                double z0, double z1,
                double *xl, double *dx, double *vol);
 
void poly_rec3d0(int dim, int nnode, double *xys,
                 double *x2min, double *x2max, double z0, double z1,
                 double *xl, double *dx, double *vol);
 
void poly_rec2d_overlap(int dim, int nnode, double *pts, double *xl, double *dx,
                        int *overlapped);
 
void poly_rec2d(int geop,
                int dim, int nnode, double *pts, double *xl, double *dx,
                double *vol);
 
void sph_poly3d(int ifinquiry, int dim, double *ctr, double r,
           double *ptrs, int nnode, int nedge, int nface,
           int *nedge_for_face, int *nodelist_for_edge, int *edgelist_for_face,
           double vcell, int *mixed, double *vol);
 
void ellipse_poly3d(int ifinquiry, int dim, double *ctr, double *r,
           double *ptrs, int nnode, int nedge, int nface,
           int *nedge_for_face, int *nodelist_for_edge, int *edgelist_for_face,
           double vcell, int *mixed, double *vol);
 
void conic_poly3d(int ifinquiry, int dim,
           double x0, double r0, double x1, double r1,
           double *pts, int nnode, int nedge, int nface,
           int *nedge_for_face, int *nodelist_for_edge, int *edgelist_for_face,
           double vcell, int *mixed, double *vol);
 
void plane_poly3d(int dim, double xint, double *ptrs,
                  int nnode, int nedge, int nface,
                  int *nedge_for_face, int *nodelist_for_edge,
                  int *edgelist_for_face,
                  int *nint, double *ptrs_int);
 
void sph_poly2d(int ifinquiry, int geop, int dim, double *ctr, double r,
                double *ptrs, int nnode, int *mixed, double *vol, double *ctroid);
 
void sph_rec(int geop, int dim, double *ctr, double r, double *xl, double *dx, double vcell, double *vol);
 
void sph_rec2d_ctr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
 
void sph_rec2d0_ll(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dx_ll(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dy_ll(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dxy_ll(int geop, int dim, double r,double *ctr, double *x0, double *dx, double *vol);

void sph_rec2d0_ul(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dx_ul(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dy_ul(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dxy_ul(int geop, int dim, double r,double *ctr, double *x0, double *dx, double *vol);
 
void sph_rec2d0_lr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dx_lr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dy_lr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dxy_lr(int geop, int dim, double r,double *ctr, double *x0, double *dx, double *vol);

void sph_rec2d0_ur(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dx_ur(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dy_ur(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol);
void sph_rec2dxy_ur(int geop, int dim, double r,double *ctr, double *x0, double *dx, double *vol);

double integral_l(double r, double *ctr, double *x0, double *dx, double y0, double y1); 
double integral_r(double r, double *ctr, double *x0, double *dx, double y1, double y2); 

void sph_rec3d(int dim, double r, double *x0, double *dx, double vcell, double *vol);
void sph_rec3dx(int dim, double r, double *x0, double *dx, double vcell, double *vol);
 
void sph_rec3d_b(int dim, double r, double *xl, double *dx, double vcell, double *vol);
 
void bintegral(double rsphere, double *d3, double *val, int *np);
 
void sintegral(double rsphere, double *d3, double x1, double *val, int *np);
 
void ssphere_int(double x0, double x1, double rsphere, double d, double *vol);
 
void sphere_int(double x0, double x1, double rsphere, double d, double *vol);
 
void dereax(double rsphere, double d2, double fixedx, double y, double *f);
 
void sph_rec3d0_old(int dim, double r, double *x0, double *dx, double *vol);
 
// void plane_rec3d(int dim, double *p, double *norm, double *xl, double *dx,
//                 double *vol);
 
void cube_plane_3(double *p, double *norm, double *vol);
 
void cube_plane_ppp(double *p, double *norm, double *vol);
 
void cube_plane_n2n1n0(double *p, double *norm, double *vol);
 
void cube_plane_n2n0n1(double *p, double *norm, double *vol);
 
void cube_plane_2(double * p, double *norm, double *vol);
 
void square_line_pp(double *p, double *norm, double *area);
 
void square_line_n1n0(double *p, double *norm, double *area);
 
void square_line_n0n1(double *p, double *norm, double *area);
 
void aintegral(double rsphere, double x0, double x1, double p0,
               double *val);
 
void sintegral_minus(double rsphere, double x0, double x1,
                    double p0, double m0,
                    double vcell, double *val, int *np);
 
void sintegral_plus(double rsphere, double x0, double x1,
                    double z0, double y0, double f0,
                    double vcell, double *val, int *np);
 
void ssph_rec3d_a1(int dim, double rsphere, double *xl, double *dx,
                   double xlower, double xupper, double vcell, double *vol);
 
void ssph_rec3d_a2(int dim, double rsphere, double *xl, double *dx,
                   double xlower, double xupper, double vcell, double *vol);
 
void ssph_rec3d_a3(int dim, double rsphere, double *xl, double *dx, double *vol);
 
void ssph_rec3d_general(int dim, double rsphere, double *xl, double *dx,
                        double volc, double *vol);
 
void ssph_rec3d_s(int dim, double rsphere, double *xl, double *dx,
                        double volc, double *vol);
void poly_rec3d_split(int nreg, int ifinquiry, int dim,
                double *xl, double *xr,
                int nface, int nedge,  int nnode,  double *coords1d,
                int *num_node_for_face, int *nodelist_for_face,
                int *edgelist_for_face, int *nodelist_for_edge,
                int *nnode_out, double *coords1d_out,
                int *nface_out, int *nnode_for_face_out,
                int *nodelist_for_face_out);
 
void poly_rec3d_unsplit(int nz, int ifinquiry, int dim, double *xl, double *xr,
                int nface, int nedge, int nnode,  double *coords1d,
                int *num_node_for_face, int *nodelist_for_face,
                int *edgelist_for_face, int *nodelist_for_edge,
                int *npoly, amr_Loc *poly_pos,
                int *nnode_out, double **coords1d_out,
                int *nface_for_zone_out, int **nnode_for_face_out,
                int **nodelist_for_face_out);
 
void poly_poly3d(int ifinquiry, int dim,
                int nface1, int nedge1,  int nnode1,  double *coords1d1,
                int *num_node_for_face1, int *nodelist_for_face1,
                int *edgelist_for_face1, int *nodelist_for_edge1,
                int nface2, int nedge2,  int nnode2,  double *coords1d2,
                int *num_node_for_face2, int *nodelist_for_face2,
                int *edgelist_for_face2, int *nodelist_for_edge2,
                int *nnode_out, double *coords1d_out,
                int *nface_out, int *nnode_for_face_out,
                int *nodelist_for_face_out);
 
void poly3d_cut_by_xplane(int dim, double xplane,
                 int nface, int nedge, int nnode,  double *coords1d,
                 int *num_edge_for_face, int *edgelist_for_face,
                 int *nodelist_for_edge, int *nodelist_for_face,
                 int *nnode_out, double *coords1d_out,
                 int *nface_for_zone_out, int **nnode_for_face_out,
                 int **nodelist_for_face_out,
                 int *nnode_intface, int *nodelist_intface);
 
void mychull_sort(double *ptr, int *map, int npoint, int *failed);
 
void clean_oso();
 
void getnline(char *filename, int nreg, int nskip_header, int *npair);
 
void get_xys(char *regname, int npes, int mype,
             int nreg, double *xys, int *npair, int nskip_header, double scale);
 
void oso_rec(int mype,
             int nreg, int ifinquiry, int *ifmixed, int geop, double vcell,
             int dim, double *xl, double *dx, double *vol);
 
void check_intersected_sphere(int dim, double *xl, double *xr,
                              double *ctr, double r, int *intersected);
 
void write_hull(char *fname, double *r1, int n1, double *r2, int n2, double *ri, int ni);
 
void grevolve_obj(geom_obj_type obj,
                  int ifinquiry, int *mixed,
                  double *pt_axis, int axis_rotate,
                  double angle_theta, double angle_phi,
                  double angle_begin, double angle_end,
            double *ctr, double radius,                   // for obj = 2, sphere
            int if_helix, double h0,   double dhdphi,     // for helix
            int nnode, double *xys, double *xys_scratch,  // for obj = 1, poly
            double *x2min, double *x2max,                 // for obj = 1, poly
            double *xl, double *dx, double *vol);
 
void axis_decompose(double *coords, int nnode, int nedge, int nface,
               int *nnode_for_face, int *nodelist_for_face,
               int *nodelist_for_edge, int *edgelist_for_face,
               double dx_scale, int *npoly_out,
               int *nnode_out, int *nedge_out, int *nface_out,
               double ***coords_out, int ***nnode_for_face_out,
               int ***nodelist_for_face_out,
               int ***nodelist_for_edge_out, int ***edgelist_for_face_out);
 
void grevolve_obj_z(geom_obj_type  obj,
              int ifinquiry, int *mixed,
              double angle_begin, double angle_end, double *pt_axis,
              double *ctr, double radius,
              int if_helix, double h0, double dhdphi,
              int nnode, double *xys, double *xys_scratch,
              double *x2min, double *x2max,
              double *xl, double *dx, double *vol);
 
void grevolve_obj_p(geom_obj_type  obj,
                    int ifinquiry, int *mixed,
                    double angle_begin, double angle_end,
                    double *ctr, double radius,
                    int if_helix, double h0, double dhdphi,
                    int nnode, double *coords2d, double *x2min, double *x2max,
                    double *coords, int nnode8, int nedge8, int nface8,
                    int *nnode_for_face, int *nodelist_for_face,
                    int *nodelist_for_edge, int *edgelist_for_face,
                    double scale, double vcell, double *vol);
 
void revolve_obj_p(geom_obj_type  obj, double angle_begin, double angle_end,
                   double *ctr, double radius,
                   int if_helix, double h0, double dhdphi,
                   int nn, double *coords2dp, double *x2min, double *x2max,
                   double *coords8, double *myphis,
                   int *inside_v, int *vofphi,
                   int nnode8, int nedge, int nface,
                   int *nodelist_for_edge, int *edgelist_for_face,
                   int *nodelist_for_face, int *nnode_for_face,
                   double scale,  double vcell, double *vol);
 
void revolve_poly_poly(geom_obj_type obj,
                      double angle_begin, double angle_end,
                      int ifstart, int ifend,
                      double *ctr, double radius,
                      int nn2d, double *coords2d,
                      double *x2min, double *x2max,
                      double *coords8, int nnode8, int nedge, int nface,
                      int    *inside_v,
                      int *nodelist_for_edge, int *edgelist_for_face,
                      int *nodelist_for_face, int *nnode_for_face,
                      double *vol, int *is_new_nnodes_inside);
 
void check_inside(geom_obj_type obj,
                  int nface, int *nnode_for_face, int *nodelist_for_face,
                  double *coords3d, int nnode_tot, int nnode_checkd, int *inside_v,
                  double *ctr, double radius,
                  double *coords2d, int nn2d,
                  double *x2min, double *x2max,
                  int *is_new_nnodes_inside);
 
void revolve_obj_slice(int obj,
                      double angle_begin, double angle_end,
                      double *ctr0, double radius,
                      int if_helix, double h0, double dhdphi,
                      int nn, double *coords2dp, double *x2min, double *x2max,
                      double *coords0, int nnode, int nedge, int nface,
                      int *nodelist_for_edge, int *edgelist_for_face,
                      int *nodelist_for_face, int *nnode_for_face,
                      double scale, double vcell, double *vol);
 
void polygon_fr_poly3d_cut_by_xplane(double xplane, double *coords,
                                     int nnode, int nedge, int *nodelist_for_edge,
                                     int *nnode_intface, double *coords2d_intface);
 
void decompose_2d(double *xl, double *xr, double *ctr, int *nrec,
                  double *xl4, double *xr4);
 
void decompose_3d(double *xl, double *xr, double *ctr, int *nrec,
                  double *xl8, double *xr8);
 
void psph_rec(int geop, int dim, double *ctr, double radius,
              double ang_start, double ang_end,
              double *xl, double *dx, double *vol);
 
void cal_vf_cyl_cord(double radius, double *ctr, double *p0, double *p1,
                     double *vol);
 
void get_3d_face_description(int nreg, int *i2dto3d, double bottom, double top);
 
void gpoly3d_rec_oso(int nreg, int ifinquiry, int *ifmixed, double vcell,
                         int ifcyl, double *xl, double *dx, double *vol);
 
 
void for_concave( int if_decompose, int mype, int nreg,  double *coords, int nnode);
void for_concave2(int if_decompose, int mype, int nreg, double *coords0, int dim, int nnode);
 
void concave_decompose0(double *coords0, int nnode0,
                       int *npart, int *nnode_part,
                       double *factors, double *coords_part[2]);
 
void concave_decompose(int mype,
                       int nreg, double *coords, int nnode, int *npart_max,
                       int *npart, int **nnode_part,
                       double **factor_part, double ***coords_part);
 


void is_point_inside_polygon(double *p, double *coords, int nnode, int *isinside);
 
 
void sph_rec_int(int dim, int geop, double *ctr, double radius,
                 double vcell, double *xl, double *dx, double *vol)
{
     int  i, k, r, nrec, nn, mixed;
     double  r2, ds2mn, ds2mx, myvol;
     double  xr[3], ds2[8], xl8[24], xr8[24], xc[3], mydx[3];
     double  c4[4][2], c8[8][3], *cnode;
     double  *c, *myxl, *myxr;
 
     r2 = radius * radius;
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     if (dim == 2) {
         c4[0][0] =  xl[0] - ctr[0];
         c4[0][1] =  xl[1] - ctr[1];
 
         c4[1][0] =  xr[0] - ctr[0];
         c4[1][1] =  xl[1] - ctr[1];
 
         c4[2][0] =  xr[0] - ctr[0];
         c4[2][1] =  xr[1] - ctr[1];
 
         c4[3][0] =  xl[0] - ctr[0];
         c4[3][1] =  xr[1] - ctr[1];
 
         cnode = c4[0];
         nn    = 4;
     }
     else if (dim == 3) {
         c8[0][0] =  xl[0] - ctr[0];
         c8[0][1] =  xl[1] - ctr[1];
         c8[0][2] =  xl[2] - ctr[2];
 
         c8[1][0] =  xr[0] - ctr[0];
         c8[1][1] =  xl[1] - ctr[1];
         c8[1][2] =  xl[2] - ctr[2];
 
         c8[2][0] =  xr[0] - ctr[0];
         c8[2][1] =  xr[1] - ctr[1];
         c8[2][2] =  xl[2] - ctr[2];
 
         c8[3][0] =  xl[0] - ctr[0];
         c8[3][1] =  xr[1] - ctr[1];
         c8[3][2] =  xl[2] - ctr[2];

         c8[4][0] =  xl[0] - ctr[0];
         c8[4][1] =  xl[1] - ctr[1];
         c8[4][2] =  xr[2] - ctr[2];

         c8[5][0] =  xr[0] - ctr[0];
         c8[5][1] =  xl[1] - ctr[1];
         c8[5][2] =  xr[2] - ctr[2];

         c8[6][0] =  xr[0] - ctr[0];
         c8[6][1] =  xr[1] - ctr[1];
         c8[6][2] =  xr[2] - ctr[2];

         c8[7][0] =  xl[0] - ctr[0];
         c8[7][1] =  xr[1] - ctr[1];
         c8[7][2] =  xr[2] - ctr[2];

         cnode = c8[0];
         nn    = 8;
     }
     for (i = 0; i < nn; i++) {
         ds2[i] = 0.0;
         c      = cnode + (dim * i);
         for (k = 0; k < dim; k++) {
             ds2[i] += (c[k] * c[k]);
         }
     }
     ds2mn = ds2[0];
     ds2mx = ds2[0];
     for (i = 1; i < nn; i++) {
         if (ds2[i] < ds2mn) {
             ds2mn = ds2[i];
         }
         else if (ds2[i] > ds2mx) {
             ds2mx = ds2[i];
         }
     }
     if (ds2mx <= r2) {
         *vol = vcell;
         return;
     }
     else if (ds2mn >= r2) {
         *vol = 0.0;
         return;
     }
     if (dim == 2) {
         decompose_2d(xl, xr, ctr, &nrec, xl8, xr8);
     }
     else if (dim == 3) {
         decompose_3d(xl, xr, ctr, &nrec, xl8, xr8);
     }
     *vol = 0.0;
     for (r = 0; r < nrec; r++) {
         myxl = xl8 + (dim * r);
         myxr = xr8 + (dim * r);
         for (k = 0; k < dim; k++) {
             xc[k] = 0.5 * (myxl[k] + myxr[k]);
             mydx[k] = myxr[k] - myxl[k];
         }
         if (dim == 2) {
             c4[0][0] =  myxl[0] - ctr[0];
             c4[0][1] =  myxl[1] - ctr[1];
 
             c4[1][0] =  myxr[0] - ctr[0];
             c4[1][1] =  myxl[1] - ctr[1];
 
             c4[2][0] =  myxr[0] - ctr[0];
             c4[2][1] =  myxr[1] - ctr[1];
 
             c4[3][0] =  myxl[0] - ctr[0];
             c4[3][1] =  myxr[1] - ctr[1];
 
             cnode = c4[0];
             nn = 4;
         }
         else if (dim == 3) {
             c8[0][0] =  myxl[0] - ctr[0];
             c8[0][1] =  myxl[1] - ctr[1];
             c8[0][2] =  myxl[2] - ctr[2];
 
             c8[1][0] =  myxr[0] - ctr[0];
             c8[1][1] =  myxl[1] - ctr[1];
             c8[1][2] =  myxl[2] - ctr[2];
 
             c8[2][0] =  myxr[0] - ctr[0];
             c8[2][1] =  myxr[1] - ctr[1];
             c8[2][2] =  myxl[2] - ctr[2];
 
             c8[3][0] =  myxl[0] - ctr[0];
             c8[3][1] =  myxr[1] - ctr[1];
             c8[3][2] =  myxl[2] - ctr[2];

             c8[4][0] =  myxl[0] - ctr[0];
             c8[4][1] =  myxl[1] - ctr[1];
             c8[4][2] =  myxr[2] - ctr[2];

             c8[5][0] =  myxr[0] - ctr[0];
             c8[5][1] =  myxl[1] - ctr[1];
             c8[5][2] =  myxr[2] - ctr[2];

             c8[6][0] =  myxr[0] - ctr[0];
             c8[6][1] =  myxr[1] - ctr[1];
             c8[6][2] =  myxr[2] - ctr[2];

             c8[7][0] =  myxl[0] - ctr[0];
             c8[7][1] =  myxr[1] - ctr[1];
             c8[7][2] =  myxr[2] - ctr[2];
 
             cnode = c8[0];
             nn = 8;
         }
         for (i = 0; i < nn; i++) {
             ds2[i] = 0.0;
             c      = cnode + (dim * i);
             for (k = 0; k < dim; k++) {
                 ds2[i] += (c[k] * c[k]);
             }
         }
         ds2mn = ds2[0];
         ds2mx = ds2[0];
         for (i = 1; i < nn; i++) {
             if (ds2[i] < ds2mn) {
                 ds2mn = ds2[i];
             }
             else if (ds2[i] > ds2mx) {
                 ds2mx = ds2[i];
             }
         }
         if (ds2mx <= r2) {
             myvol = mydx[0];
             for (i = 1; i < dim; i++) {
                 myvol *= mydx[i];
             }
             if ((geop == 2) && (dim == 2)) {
                 myvol *= fabs(myxl[0] + myxr[0]);
             }
             else if (!((geop == 1)||(dim == 3))) {
                 printf("ERROR: geop and dim in ..\n");
                 exit(1);
             }
         }
         else if (ds2mn >= r2) {
             myvol = 0.0;
         }
         else {
             gsph_rec(0, &mixed, geop, dim, ctr, radius, myxl, mydx, &myvol);
         }
         *vol += (myvol);
     }
     return;
 }
 
void decompose_3d(double *xl, double *xr, double *ctr, int *nrec,
                  double *xl8, double *xr8)
{
     double *myxl, *myxr;
     double small;
 
     small = 1.0e-12 *(xr[0] - xl[0]);
 
     if ((xl[0] + small < ctr[0]) && (xr[0] - small > ctr[0]) &&
         (xl[1] + small < ctr[1]) && (xr[1] - small > ctr[1]) &&
         (xl[2] + small < ctr[2]) && (xr[2] - small > ctr[2]) ) {
        *nrec = 8;
 
         myxl = xl8;
         myxr = xr8;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = ctr[0];
         myxr[1] = ctr[1];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 3;
         myxr = xr8 + 3;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
         myxr[0] = xr[0];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 6;
         myxr = xr8 + 6;
         myxl[0] = ctr[0];
         myxl[1] = ctr[1];
         myxl[2] =  xl[2];
         myxr[0] =  xr[0];
         myxr[1] =  xr[1];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 9;
         myxr = xr8 + 9;
         myxl[0] =  xl[0];
         myxl[1] = ctr[1];
         myxl[2] =  xl[2];
         myxr[0] = ctr[0];
         myxr[1] =  xr[1];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 12;
         myxr = xr8 + 12;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = ctr[2];
         myxr[0] = ctr[0];
         myxr[1] = ctr[1];
         myxr[2] =  xr[2];
 
         myxl = xl8 + 3;
         myxr = xr8 + 3;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxl[2] = ctr[2];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
         myxr[0] = xr[0];
         myxr[2] =  xr[2];
 
         myxl = xl8 + 6;
         myxr = xr8 + 6;
         myxl[0] = ctr[0];
         myxl[1] = ctr[1];
         myxl[2] = ctr[2];
         myxr[0] =  xr[0];
         myxr[1] =  xr[1];
         myxr[2] =  xr[2];
 
         myxl = xl8 + 9;
         myxr = xr8 + 9;
         myxl[0] =  xl[0];
         myxl[1] = ctr[1];
         myxl[2] = ctr[2];
         myxr[0] = ctr[0];
         myxr[1] =  xr[1];
         myxr[2] =  xr[2];
     }
     else if ((xl[0] + small < ctr[0]) && (xr[0] - small > ctr[0]) &&
              (xl[2] + small < ctr[2]) && (xr[2] - small > ctr[2]) ) {
 
         *nrec = 4;
 
         myxl = xl8;
         myxr = xr8;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = ctr[0];
         myxr[1] =  xr[1];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 3;
         myxr = xr8 + 3;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 6;
         myxr = xr8 + 6;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = ctr[2];
         myxr[0] = ctr[0];
         myxr[1] =  xr[1];
         myxr[2] =  xr[2];
 
         myxl = xl8 + 9;
         myxr = xr8 + 9;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxl[2] = ctr[2];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
         myxr[2] = xr[2];
 
     }
     else if ((xl[1] + small < ctr[1]) && (xr[1] - small > ctr[1]) &&
              (xl[2] + small < ctr[2]) && (xr[2] - small > ctr[2])) {
         *nrec = 4;
 
         myxl = xl8;
         myxr = xr8;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 3;
         myxr = xr8 + 3;
         myxl[0] = xl[0];
         myxl[1] = ctr[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
         myxr[2] = ctr[2];
 
         myxl = xl8 + 6;
         myxr = xr8 + 6;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = ctr[2];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
         myxr[2] = xr[2];
 
         myxl = xl8 + 9;
         myxr = xr8 + 9;
         myxl[0] = xl[0];
         myxl[1] = ctr[1];
         myxl[2] = ctr[2];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
         myxr[2] = xr[2];
     }
     else if ((xl[0] + small < ctr[0]) && (xr[0] - small > ctr[0]) &&
              (xl[1] + small < ctr[1]) && (xr[1] - small > ctr[1]))   {
 
        *nrec = 4;
 
         myxl = xl8;
         myxr = xr8;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = ctr[0];
         myxr[1] = ctr[1];
         myxr[2] = xr[2];
 
         myxl = xl8 + 3;
         myxr = xr8 + 3;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
         myxr[0] = xr[0];
         myxr[2] = xr[2];
 
         myxl = xl8 + 6;
         myxr = xr8 + 6;
         myxl[0] = ctr[0];
         myxl[1] = ctr[1];
         myxl[2] =  xl[2];
         myxr[0] =  xr[0];
         myxr[1] =  xr[1];
         myxr[2] =  xr[2];
 
         myxl = xl8 + 9;
         myxr = xr8 + 9;
         myxl[0] =  xl[0];
         myxl[1] = ctr[1];
         myxl[2] =  xl[2];
         myxr[0] = ctr[0];
         myxr[1] =  xr[1];
         myxr[2] =  xr[2];
     }
     else if ((xl[0] + small < ctr[0]) && (xr[0] - small > ctr[0])) {
 
         *nrec = 2;
 
         myxl = xl8;
         myxr = xr8;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = ctr[0];
         myxr[1] =  xr[1];
         myxr[2] =  xr[2];
 
         myxl = xl8 + 3;
         myxr = xr8 + 3;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
         myxr[2] = xr[2];
     }
     else if ((xl[1] + small < ctr[1]) && (xr[1] - small > ctr[1])) {
         *nrec = 2;
 
         myxl = xl8;
         myxr = xr8;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
         myxr[2] = xr[2];
 
         myxl = xl8 + 3;
         myxr = xr8 + 3;
         myxl[0] = xl[0];
         myxl[1] = ctr[1];
         myxl[2] = xl[2];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
         myxr[2] = xr[2];
     }
     else {
         *nrec = 1;
 
         xl8[0] = xl[0];
         xl8[1] = xl[1];
         xl8[2] = xl[2];
         xr8[0] = xr[0];
         xr8[1] = xr[1];
         xr8[2] = xr[2];
     }
     return;
 }
 
 
void decompose_2d(double *xl, double *xr, double *ctr, int *nrec,
                  double *xl4, double *xr4)
{
     double *myxl, *myxr;
     double small;
 
     small = 1.0e-10 *(xr[0] - xl[0]);
 
     if ((xl[0] + small < ctr[0]) && (xr[0] - small > ctr[0]) &&
         (xl[1] + small < ctr[1]) && (xr[1] - small > ctr[1])) {
        *nrec = 4;
 
         myxl = xl4;
         myxr = xr4;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxr[0] = ctr[0];
         myxr[1] = ctr[1];
 
         myxl = xl4 + 2;
         myxr = xr4 + 2;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
 
         myxl = xl4 + 4;
         myxr = xr4 + 4;
         myxl[0] = ctr[0];
         myxl[1] = ctr[1];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
 
         myxl = xl4 + 6;
         myxr = xr4 + 6;
         myxl[0] = xl[0];
         myxl[1] = ctr[1];
         myxr[0] = ctr[0];
         myxr[1] = xr[1];
     }
     else if ((xl[0] + small < ctr[0]) && (xr[0] - small > ctr[0])) {
 
         *nrec = 2;
 
         myxl = xl4;
         myxr = xr4;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxr[0] = ctr[0];
         myxr[1] = xr[1];
 
         myxl = xl4 + 2;
         myxr = xr4 + 2;
         myxl[0] = ctr[0];
         myxl[1] = xl[1];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
     }
     else if ((xl[1] + small < ctr[1]) && (xr[1] - small > ctr[1])) {
         *nrec = 2;
 
         myxl = xl4;
         myxr = xr4;
         myxl[0] = xl[0];
         myxl[1] = xl[1];
         myxr[0] = xr[0];
         myxr[1] = ctr[1];
 
         myxl = xl4 + 2;
         myxr = xr4 + 2;
         myxl[0] = xl[0];
         myxl[1] = ctr[1];
         myxr[0] = xr[0];
         myxr[1] = xr[1];
     }
     else {
         *nrec = 1;
 
         xl4[0] = xl[0];
         xl4[1] = xl[1];
         xr4[0] = xr[0];
         xr4[1] = xr[1];
     }
     return;
 }
 
 
void clean_oso()
{
     int i, g, ngrp;
     FILE *fp;
 
     for (i = 0; i < nfile_oso; i++) {
         fp = fp_oso[i];
         if (fp) fclose(fp);
         fp_oso[i] = NULL;
         if (xys_oso[i]) {
             free(xys_oso[i]);
             xys_oso[i] = NULL;
         }
         nxy_oso[i] = 0;
 
 
     }
 
 }
void getnline(char *filename, int nreg, int nskip_header, int *npair)
{
     int i, i0, c;
     FILE *fp;
 
     i0 = -1;
     for (i = 0; i < nfile_oso; i++) {
         if (reg_oso[i] == nreg) {
             i0 = i;
             break;
         }
     }
     if (i0 >= 0) {
         *npair = nxy_oso[i0];
     }
     else {
         *npair = 0;
         fp = fopen(filename, "r");
         if (!fp) { 
             printf("file %s not found\n", filename);
             assert(fp);
         }
         reg_oso[nfile_oso] = nreg;
         fp_oso[nfile_oso]  = fp;
         xys_oso[nfile_oso] = NULL;
         nxy_oso[nfile_oso] = 0;
 
         c = getc(fp);
         while (c != EOF) {
             if (c == '\n') {
                 (*npair)++;
             }
             c = getc(fp);
         }
         (*npair) -= (nskip_header + 1);  // 1 for ending symbol
         nxy_oso[nfile_oso] = *npair;
         xys_oso[nfile_oso] = NULL;
 
         nfile_oso++;
     }
     return;
}
 
void get_xys(char *regname, int npes, int mype,
             int nreg, double *xys, int *npair, int nskip_header, double scale)
{
     char name[128];
     int ifile, i, i0, i1, i2, k, istart, szdim2;
     int nl, nn, c, nxy, alreadyread, nnode;
     int *selected;
     double small, scaleinv, cross, *p, *p0, *p1, *p2, dp0[2], dp1[2];
     FILE *fp;
 
     small = 1.0e-12;
     szdim2 = 2 * sizeof(double);
 
     fp = NULL;
     ifile = -1;
     for (i = 0; i < nfile_oso; i++) {
         if (reg_oso[i] == nreg) {
             fp = fp_oso[i];
             ifile = i;
             break;
         }
     }
     assert(ifile >= 0);
     alreadyread = 0;
//   assert(*npair == nxy_oso[ifile]);
 
     nn = *npair;
     if (xys_oso[ifile]) {
         alreadyread = 1;
     }
     if (alreadyread) {
         nxy = nxy_oso[ifile];
         memcpy(xys, xys_oso[ifile], (size_t)((nxy + nxy + 2) * sizeof(double)));
     }
     else {
         assert(fp);
         rewind(fp);
 
         nl = 0;
         c = getc(fp);
         while ((c != EOF) && (nl < nskip_header)) {
             if (c == '\n') {
                 nl++;
             }
             if (nl < nskip_header) {
                 c = getc(fp);
             }
         }
         p  = xys;
         for (i = 0; i < nn; i++) {
             fscanf(fp, "%le%le", p, p + 1);
             p += 2;
         }
         p = xys + (nn + nn - 2);  // the last data in the file
         if ((p[0] == xys[0]) && (p[1] == xys[1])) {
             nn--;
         }
         else {
             p = xys + (nn + nn);
             p[0] = xys[0];
             p[1] = xys[1];
         }
         selected = (int *) malloc(nn * sizeof(int));
         for (i = 0; i < nn; i++) {
             selected[i] = 1;
         }
         istart = -1;
         scaleinv = 1.0/scale;
         for (i = 0; i < nn; i++) {
             i1 = (i  + 1) % nn;
             i2 = (i1 + 1) % nn;
             p0 = xys + (i  + i);
             p1 = xys + (i1 + i1);
             p2 = xys + (i2 + i2);
             for (k = 0; k < 2; k++) {
                 dp0[k] = scaleinv *(p1[k] - p0[k]);
                 dp1[k] = scaleinv *(p2[k] - p1[k]);
             }
             cross = dp0[0] * dp1[1] - dp0[1] * dp1[0];
             if (fabs(cross) > small) {
                 istart = i;
                 break;
             }
         }
         assert(istart >= 0);
         i = 0;
         i1 = (istart + 1) % nn;
         while ((i < nn) && (i1 != istart)) {
             i0 = (i + istart) % nn;
             i1 = (i0 + 1) % nn;
             while ((selected[i1] == 0) && (i1 != i0)) {
                i1 = (i1 + 1) % nn;
             }
             i2 = (i1 + 1) % nn;
             p0 = xys + (i0 + i0);
             p1 = xys + (i1 + i1);
             p2 = xys + (i2 + i2);
             for (k = 0; k < 2; k++) {
                 dp0[k] = scaleinv *(p1[k] - p0[k]);
                 dp1[k] = scaleinv *(p2[k] - p1[k]);
             }
             cross = dp0[0] * dp1[1] - dp0[1] * dp1[0];
 
             if (fabs(cross) < small) {
                 selected[i1] = 0;
             }
             else {
                 selected[i1] = 1;
                 i++;
             }
         }
         *npair = 0;
         for (i = 0; i < nn; i++) {
             if (selected[i]) (*npair)++;
         }
         xys_oso[ifile] = (double *) malloc(  (*npair + *npair + 2) * sizeof(double));
         if (nn == *npair) {
             memcpy(xys_oso[ifile], xys, (size_t)((*npair + *npair + 2) * sizeof(double)));
         }
         else {
             *npair = 0;
             p1 = xys_oso[ifile];
             for (i = 0; i < nn; i++) {
                 if (selected[i]) {
                     p0 = xys + (i + i);
                     memcpy(p1, p0, (size_t)szdim2);
                     p1 += 2;
                     (*npair)++;
                 }
             }
             memcpy(p1, xys_oso[ifile], (size_t)szdim2);
         }
         memcpy(xys, xys_oso[ifile], (size_t)((*npair + *npair + 2) * sizeof(double)));
 
         nxy_oso[ifile]  = *npair;
 
         free(selected);
     }
     return;
}
 
void for_concave2(int if_decompose, int mype, int nreg, double *coords0, int dim, int nnode)
{
     int szdim, i;
     double *coords, *mycoords, *p0, *p1;
 
     coords = NULL;
     szdim  = 2 * sizeof(double);
 
     if (dim == 3) {
         coords   = (double *) malloc((nnode + nnode + 2) * sizeof(double));
         mycoords = coords;
 
         p0 = coords;
         p1 = coords0;
         for (i = 0; i < nnode; i++) {
             memcpy(p0, p1, (size_t)szdim);
             p0 += 2;
             p1 += dim;
         }
         p0 = coords + (nnode + nnode);
         memcpy(p0, coords, (size_t)szdim);
     }
     else if (dim == 2) {
         mycoords = coords0;
     }
     for_concave(if_decompose, mype, nreg, mycoords, nnode);
 
     if (coords) free(coords);
 
     return;
 }
 
 
void for_concave(int if_decompose, int mype, int nreg,  double *coords, int nnode)
{
     char name[128];
     int npart_max, ifile, i, i1, npart, nedge, nn, p, fileid, meshid;
     int *nnode_part, *nodelist_for_edge;
     double *factor_part, *mycoords;
     double **coords_part;
 
     if_decompose_concave = if_decompose;
 
     npart_max    = 64;
     nnode_part   = (int     *) malloc(npart_max * sizeof(int));
     factor_part  = (double  *) malloc(npart_max * sizeof(double));
     coords_part  = (double **) malloc(npart_max * sizeof(double *));
 
     if (if_decompose_concave) {
         concave_decompose(mype,
                       nreg, coords, nnode, &npart_max, &npart, &nnode_part,
                       &factor_part, &coords_part);
     }
     else {
         npart = 1;
         nnode_part[0] = nnode;
         factor_part[0] = 1.0;
         coords_part[0] = (double *) malloc((nnode + nnode + 2) * sizeof(double));
         memcpy(coords_part[0], coords, (size_t)((nnode + nnode) * sizeof(double)));
     }
     ifile  = -1;
     for (i = 0; i < nfile_oso; i++) {
         if (reg_oso[i] == nreg) {
             ifile = i;
             break;
         }
     }
     if (ifile >= 0) {
         if (xys_oso[ifile]) {
             free(xys_oso[ifile]);
             xys_oso[ifile] = NULL;
         }
     }
     else {
         reg_oso[nfile_oso] = nreg;
         nxy_oso[nfile_oso] = nnode;
         xys_oso[nfile_oso] = NULL;
         fp_oso[nfile_oso]  = NULL;
         ifile = nfile_oso;
         nfile_oso++;
     }
     npart_oso[ifile] = npart;
     factor_part_oso[ifile] = factor_part;
     nnode_part_oso[ifile]  = nnode_part;
     coords_part_oso[ifile] = coords_part;
 
     nodelist_for_edge = (int *) malloc((nnode + nnode) * sizeof(int));
     sprintf(name, "region_%d", nreg);
//     mio_open_file(name,mio_file_create, &fileid);
 
     nedge = 0;
     for (i = 0; i < nnode; i++) {
         i1 = (i + 1) % nnode;
         nodelist_for_edge[i+i  ] = i;
         nodelist_for_edge[i+i+1] = i1;
     }
     if (mype == 0) {
         nedge = nnode;
     }
//     write_edge(fileid, 2, name, nedge, nnode, coords, nodelist_for_edge, &meshid);
 
     for (p = 0; p < npart; p++) {
         nn = nnode_part_oso[ifile][p];
         mycoords = coords_part_oso[ifile][p];
 
 
         for (i = 0; i < nn; i++) {
             i1 = (i + 1) % nn;
             nodelist_for_edge[i+i  ] = i;
             nodelist_for_edge[i+i+1] = i1;
         }
         if (mype == 0) {
             nedge = nn;
         }
         sprintf(name, "region_%d_part_%d", nreg, p);
//         write_edge(fileid, 2, name, nedge, nn, mycoords, nodelist_for_edge, &meshid);
     }
//     mio_close_file(fileid);
 
     return;
 }
 
 
void concave_decompose(int mype,
                       int nreg, double *coords, int nnode, int *npart_max,
                       int *npart, int **nnode_part,
                       double **factor_part, double ***coords_part)
{
    int np, p, i, i0, i1, imx, k, dim, szdim;
    int nn, n, new_npart_max, iter, mynnode, nedge, fileid, meshid;
    int nnode2[2], *new_nnode_part, *nodelist_for_edge;
    double factor2[2], *coord, *coords2[2], *new_factor_part;
    double *rm, *r0, *rp, dsm, dsp, cosa_mx, cosa, dm[2], dp[2];
    double **new_coords_part, *mycoords;
    char name[128];
 
    dim   = 2;
    szdim = dim * sizeof(double);
 
    *npart = 0;
 
    np   = 0;
    iter = 0;
 
    coords2[0] = NULL;
    coords2[1] = NULL;
 
    coord = coords;
    nn    = nnode;
 
    sprintf(name, "region_%d_debug", nreg);
//    mio_open_file(name,mio_file_create, &fileid);
    nodelist_for_edge = (int *) malloc((nnode + nnode + 1) * sizeof(int));
 
    while (np != 1) {
          concave_decompose0(coord, nn, &np, nnode2, factor2, coords2);
 
          if (*npart + np > *npart_max) {
              new_npart_max   = *npart_max + *npart_max;
              new_nnode_part  = (int     *) malloc(new_npart_max * sizeof(int));
              new_factor_part = (double  *) malloc(new_npart_max * sizeof(double));
              new_coords_part = (double **) malloc(new_npart_max * sizeof(double *));
              for (i = 0; i < *npart; i++) {
                  n = (*nnode_part)[i];
                  new_nnode_part[i]  = n;
                  new_coords_part[i] = (*coords_part)[i];
                  new_factor_part[i] = (*factor_part)[i];
              }
              free(*coords_part);
              free(*nnode_part);
              free(*factor_part);
              *coords_part = new_coords_part;
              *nnode_part  = new_nnode_part;
              *factor_part = new_factor_part;
              *npart_max   = new_npart_max;
          }
          if (np > 1) {
              if ((factor2[0] < 0.0) || (factor2[1] < 0.0)) {
                  free(coords2[0]);
                  free(coords2[1]);
 
                  cosa_mx = -2.0;
                  for (i = 0; i < nn; i++) {
                      i0 = (i - 1 + nn) % nn;
                      i1 = (i + 1)      % nn;
                      rm = coord + (i0 + i0);
                      r0 = coord + (i  + i);
                      rp = coord + (i1 + i1);
                      for (k = 0; k < dim; k++) {
                          dm[k] = rm[k] - r0[k];
                          dp[k] = rp[k] - r0[k];
                      }
                      dsm = 1.0/sqrt(dm[0]*dm[0] + dm[1]*dm[1]);
                      dsp = 1.0/sqrt(dp[0]*dp[0] + dp[1]*dp[1]);
                      for (k = 0; k < dim; k++) {
                          dm[k] *= dsm;
                          dp[k] *= dsp;
                      }
                      cosa = dm[0] * dp[0] + dm[1] * dp[1];
                      if (cosa > cosa_mx) {
                          imx = i;
                      }
                  }
                  if (imx != 1) {
                      mycoords = (double *) malloc((nn + nn + 2) * sizeof(double));
                      for (i = 0; i < nn; i++) {
                          i0 = (i + imx - 1 + nn) % nn;
                          memcpy(mycoords + (i + i), coord + (i0 + i0), (size_t)szdim);
                      }
                      memcpy(mycoords + (nn + nn), mycoords, (size_t)szdim);
                  }
                  free(coord);
                  coord = mycoords;
                  continue;
              }
              (*nnode_part)[ *npart] = nnode2[1];
              (*coords_part)[*npart] = coords2[1];
              (*factor_part)[*npart] = factor2[1];
 
//////
//            for debug
 
              for (i = 0; i < nn; i++) {
                  i1 = (i + 1) % nn;
                  nodelist_for_edge[i+i  ] = i;
                  nodelist_for_edge[i+i+1] = i1;
              }
              sprintf(name, "contour_split_%d_original", *npart);
              if (mype == 0) {
                  nedge = nn;
              }
              else {
                  nedge = 0;
              }
//              write_edge(fileid, 2, name, nedge, nedge, coord, nodelist_for_edge, &meshid);
 
              for (p = 0; p < np; p++) {
                  mynnode = nnode2[p];
                  mycoords = coords2[p];
 
                  for (i = 0; i < mynnode; i++) {
                      i1 = (i + 1) % mynnode;
                      nodelist_for_edge[i+i  ] = i;
                      nodelist_for_edge[i+i+1] = i1;
                  }
                  if (mype == 0) {
                      nedge = mynnode;
                  }
                  else {
                      nedge = 0;
                  }
                  sprintf(name, "contour_split_%d_part_%d", *npart, p);
//                  write_edge(fileid, 2, name, nedge, nedge, mycoords, nodelist_for_edge, &meshid);
              }
////////
 
              coord      = coords2[0];
              nn         = nnode2[0];
              coords2[0] = NULL;
              coords2[1] = NULL;
              (*npart)++;
          }
     }
     sprintf(name, "contour_split_%d_original", *npart);
     for (i = 0; i < nnode2[0]; i++) {
         i1 = (i + 1) % nn;
         nodelist_for_edge[i+i  ] = i;
         nodelist_for_edge[i+i+1] = i1;
     }
     if (mype == 0) {
         nedge = nnode2[0];
     }
     else {
         nedge = 0;
     }
//     write_edge(fileid, 2, name, nedge, nedge, coords2[0], nodelist_for_edge, &meshid);
     free(nodelist_for_edge);
//     mio_close_file(fileid);
 
     (*nnode_part)[ *npart] = nnode2[0];
     (*coords_part)[*npart] = coords2[0];
     (*factor_part)[*npart] = factor2[0];
     (*npart)++;
 
     return;
 }
 
void concave_decompose0(double *coords0, int nnode0,
                       int *npart, int *nnode_part,
                       double *factors, double *coords_part[2])
{
//   coords is (nnode + 1) long, and circular.
 
     int itest, p, nn2[2], isinside, isconcave2[2];
     int dim, szdim, nnode, nn, nn0, n1, n2, i1_save;
     int just_turned, islasti0, i02[2];
     int ncave, i, j, k, i0, i1, i2, istart, nalloc;
     int istart_old, i0_old, i_concave, i_convex;
     int *i0s, *mark, *choose_i1;
     double small, dx, factor, myfactor, cross0, cross_last, cross;;
     double *coords, *crosses, *mycoords, *mycrosses, *p0, *p1, *p2;
     double dp0[2], dp1[2], xmid[2];
     double *coords2[2];
 
     small = 1.0e-12;
 
     dim    = 2;
     szdim  = dim * sizeof(double);
 
     nalloc = 64;
     i0s = (int *) malloc(nalloc * sizeof(int));
 
//   remove the middle point in any line segment.
 
     coords = (double *) malloc((nnode0 + 1) * szdim);
     mark   = (int    *) malloc( nnode0 * sizeof(int));
     for (i = 0; i < nnode0; i++) {
         mark[i] = 0;
     }
 
     i  = 0;
     i1 = 1;
     while (i1 < nnode0) {
         i2 = (i1 + 1) % nnode0;
         p0 = coords0 + (i  + i);
         p1 = coords0 + (i1 + i1);
         p2 = coords0 + (i2 + i2);
         dx = 0.0;
         for (k = 0; k < 2; k++) {
             dp0[k] = p1[k] - p0[k];
             dp1[k] = p2[k] - p1[k];
             dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
         }
         factor = 1.0/MAX(dx, small);
         for (k = 0; k < 2; k++) {
             dp0[k] *= factor;
             dp1[k] *= factor;
         }
         cross0 = dp0[0] * dp1[1] - dp0[1] * dp1[0];
         if (fabs(cross0) < small) {
             mark[i1] = 1;
         }
         else {
             i = i1;
         }
         i1++;
     }
     nnode = 0;
     for (i = 0; i < nnode0; i++) {
         if (mark[i]) continue;
         memcpy(coords + (nnode * dim), coords0 + (i * dim), (size_t)szdim);
         nnode++;
     }
     memcpy(coords + (nnode * dim), coords, (size_t)szdim);
     free(mark);
 
     for (i = 0; i < nnode; i++) {
         i1 = (i  + 1) % nnode;
         i2 = (i1 + 1) % nnode;
         p0 = coords + (i  + i);
         p1 = coords + (i1 + i1);
         p2 = coords + (i2 + i2);
         dx = 0.0;
         for (k = 0; k < 2; k++) {
             dp0[k] = p1[k] - p0[k];
             dp1[k] = p2[k] - p1[k];
             dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
         }
         factor = 1.0/MAX(dx, small);
         for (k = 0; k < 2; k++) {
             dp0[k] *= factor;
             dp1[k] *= factor;
         }
         cross0 = dp0[0] * dp1[1] - dp0[1] * dp1[0];
//       if (fabs(cross0) > small) {
         if (cross0 != 0.0) {
             istart = i;
             break;
         }
     }
     crosses = (double *) malloc(nnode * sizeof(double));
     for (i = 0; i < nnode; i++) {
         i0 = (i + istart) % nnode;
         i1 = (i0 + 1) % nnode;
         i2 = (i1 + 1) % nnode;
         p0 = coords + (i0 + i0);
         p1 = coords + (i1 + i1);
         p2 = coords + (i2 + i2);
         for (k = 0; k < 2; k++) {
             dp0[k] = p1[k] - p0[k];
             dp1[k] = p2[k] - p1[k];
             dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
         }
         factor = 1.0/MAX(dx, small);
         for (k = 0; k < 2; k++) {
             dp0[k] *= factor;
             dp1[k] *= factor;
         }
         cross = dp0[0] * dp1[1] - dp0[1] * dp1[0];
//       if (fabs(cross) < small) {
//           crosses[i0] = 0.0;
//       }
         crosses[i0] = cross;
     }
     choose_i1 = (int *) malloc(nnode * sizeof(int));
     for (i = 0; i < nnode; i++) {
         choose_i1[i] = -1;
     }
     just_turned = 0;
     cross_last  = cross0;
     for (i = 0; i <= nnode+1; i++) {
         i0 = (i  + istart) % nnode;
         i1 = (i0 + 1)      % nnode;
         cross = crosses[i0];
//       if (fabs(cross) < small) continue;
 
         if (cross_last * cross < 0.0) {
             if (just_turned) {
                 choose_i1[i0] = 1;
                 just_turned = 0;
             }
             else {
                 just_turned = 1;
             }
             cross_last = cross;
         }
         else {
             just_turned = 0;
         }
     }
     cross_last = cross0;
 
     i1_save = -1;
     ncave = 0;
     islasti0 = 0;
     for (i = 0; i <= nnode; i++) {
         i0 = (i  + istart) % nnode;
         i1 = (i0 + 1) % nnode;
         cross = crosses[i0];
 
//       if (fabs(cross) < small) continue;
 
         if (cross * cross_last < 0.0) {
             if (choose_i1[i1] >= 0) {
                 if (islasti0) {
                     if (ncave > 0) {
                         if (i0s[ncave-1] != i1) {
                             i0s[ncave] = i1;
                             ncave++;
                         }
                     }
                     else {
                         i0s[ncave] = i1;
                         ncave++;
                     }
                     islasti0   = 0;
                 }
                 else {
                     i1_save = i1;
                 }
             }
             else {
                 if (islasti0) {
                     if (ncave > 0) {
                         if (i0s[ncave-1] != i1) {
                             i0s[ncave] = i1;
                             ncave++;
                             islasti0 = 0;
                         }
                     }
                     else {
                         i0s[ncave] = i1;
                         ncave++;
                         islasti0 = 0;
                     }
                 }
                 else {
                     if (ncave > 0) {
                         if (i0s[(ncave)-1] != i0) {
                             i0s[ncave] = i0;
                             ncave++;
                             islasti0 = 1;
                         }
                     }
                     else {
                         i0s[ncave] = i0;
                         ncave++;
                         islasti0 = 1;
                     }
                 }
             }
             cross_last = cross;
             if (ncave == 1) break;
//           if (ncave == 2) break;
         }
     }
 
 
 
/*****
 
     ncave       = 0;
     just_turned = 0;
     cross_last  = cross0;
     for (i = 0; i < nnode; i++) {
         i0 = (i  + istart) % nnode;
         i1 = (i0 + 1)      % nnode;
         cross = crosses[i0];
//       if (fabs(cross) < small) continue;
 
         if (cross_last * cross < 0.0) {
             if (!just_turned) {
                 i0s[ncave] = i1;
                 ncave++;
                 just_turned = 1;
             }
             else {
                 just_turned = 0;
             }
             cross_last = cross;
         }
         else {
             just_turned = 0;
         }
     }
 
*****/
 
 
     if ((ncave == 0) && (i1_save >= 0)) {
         i0s[ncave] = i1;
         ncave++;
     }
//     if (ncave == 1) {
//         itest  = (nnode -1 + istart) % nnode;
//         if ((itest == i0s[ncave-1]+1) || (itest + 1 == i0s[ncave-1])) {
//             itest = (nnode + istart) % (nnode + 1);
//         }
//         i0s[ncave] =  itest;
//         ncave++;
//     }
     *npart = 0;
     if (ncave == 0) {
         *npart = 1;
         nnode_part[0]  = nnode;
         factors[0]     = 1.0;
         coords_part[0] = (double *) malloc((nnode + 1) * szdim);
         memcpy(coords_part[0], coords, (size_t)((nnode + 1) * szdim));
     }
     else if (ncave == 1) {
         if ((i0s[0] + 1 == istart) || (istart + 1 == i0s[0])) {
             i0s[0] = (i0s[0] + 1) % nnode;
         }
         *npart = 2;
 
//       check the edge of (istart, i0s[0]) is outside of the polygon
 
         p0 = coords + (istart * dim);
         p1 = coords + (i0s[0] * dim);
         xmid[0] = 0.5 *(p0[0] + p1[0]);
         xmid[1] = 0.5 *(p0[1] + p1[1]);
         is_point_inside_polygon(xmid, coords, nnode, &isinside);
         if (isinside) {
 
//           check whether (istart, i0s[0]) is a convex polygon
 
             if (i0s[0] > istart) {
                 nn2[1] = i0s[0] - istart + 1;
             }
             else {
                 nn2[1] = (nnode - istart) + i0s[0] + 1;
             }
             mycoords  = (double *) malloc((nn2[1] + 1) * szdim);
             mycrosses = (double *) malloc( nn2[1] * sizeof(double));
 
             p0 = mycoords;
             for (i = 0; i < nn2[1]; i++) {
                 i0 = (nn2[1] - 1 - i + istart) % nnode;    // reverse order
//               i0 = (i + istart) % nnode;   // original order
                 memcpy(p0 + (i * dim), coords + (i0 * dim), (size_t)szdim);
             }
             for (i = 0; i < nn2[1]; i++) {
                 i1 = (i  + 1) % nn2[1];
                 i2 = (i1 + 1) % nn2[1];
                 p0 = coords + (i  + i );
                 p1 = coords + (i1 + i1);
                 p2 = coords + (i2 + i2);
                 for (k = 0; k < 2; k++) {
                     dp0[k] = p1[k] - p0[k];
                     dp1[k] = p2[k] - p1[k];
                     dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
                 }
                 factor = 1.0/MAX(dx, small);
                 for (k = 0; k < 2; k++) {
                     dp0[k] *= factor;
                     dp1[k] *= factor;
                 }
                 cross = dp0[0] * dp1[1] - dp0[1] * dp1[0];
                 mycrosses[i] = cross;
             }
             i0 = -1;
             for (i = 0; i < nn2[1]; i++) {
                 if (mycrosses[i] * mycrosses[0] < 0.0) {
                     i0 = i;
                     break;
                 }
             }
             free(mycoords);
             free(mycrosses);
             if (i0 < 0) {
                 myfactor = 1.0;
             }
             else {
                 i = nn2[1] - i0;
                 istart = (istart + i) % nnode;
                 if (i0s[0] > istart) {
                     nn2[1] = i0s[0] - istart + 1;
                 }
                 else {
                     nn2[1] = (nnode - istart) + i0s[0] + 1;
                 }
                 if (nn2[1] < 3) {
                     myfactor = -1.0;   // failed
                 }
                 else {
                      myfactor = 1.0;
                 }
             }
         }
         else {
             istart_old = istart;
             i0_old     = i0s[0];
             istart = (istart + nnode - 1) % nnode;   // istart - 1
             i0s[0] = (i0s[0] + nnode - 1) % nnode;   // i0 - 1
             p0 = coords + (istart * dim);
             p1 = coords + (i0s[0] * dim);
             xmid[0] = 0.5 *(p0[0] + p1[0]);
             xmid[1] = 0.5 *(p0[1] + p1[1]);
             is_point_inside_polygon(xmid, coords, nnode, &isinside);
             if (isinside) {
                 myfactor = 1.0;
             }
             else {
                 istart = (istart_old + 1) % nnode;   // istart + 1
                 i0s[0] = (i0_old     + 1) % nnode;
                 p0 = coords + (istart * dim);
                 p1 = coords + (i0s[0] * dim);
                 xmid[0] = 0.5 *(p0[0] + p1[0]);
                 xmid[1] = 0.5 *(p0[1] + p1[1]);
                 is_point_inside_polygon(xmid, coords, nnode, &isinside);
                 if (isinside) {
                     myfactor = 1.0;
                 }
                 else {
                     myfactor = -1.0;   // nested this will cause trouble
                 }
             }
         }
         if (i0s[0] > istart) {
             nn2[1] = i0s[0] - istart + 1;
         }
         else {
             nn2[1] = (nnode - istart) + i0s[0] + 1;
         }
//         nn2[1] = i0s[0] - istart + 1;
         nn2[0] = nnode  - nn2[1] + 2;
 
         coords2[0] = (double *) malloc((nn2[0] + 1) * szdim);
         coords2[1] = (double *) malloc((nn2[1] + 1) * szdim);
 
         p0 = coords2[1];
         for (i = 0; i < nn2[1]; i++) {
             i0 = (i + istart) % nnode;
             memcpy(p0 + (i * dim), coords + (i0 * dim), (size_t)szdim);
         }
         p0 = coords2[0];
         memcpy(p0, coords + (istart * dim), (size_t)szdim);
         p0 += 2;
         for (i = 0; i < nn2[0]-1; i++) {
//           i0 = (i + istart + i0s[0]) % nnode;
             i0 = (i + i0s[0]) % nnode;
             memcpy(p0 + (i * dim), coords + (i0 * dim), (size_t)szdim);
         }
 
/*****
         for (p = 0; p < 2; p++) {
             for (i = 0; i < nn2[p]; i++) {
                 i1 = (i  + 1) % nn2[p];
                 i2 = (i1 + 1) % nn2[p];
                 p0 = coords2[p] + (i  + i);
                 p1 = coords2[p] + (i1 + i1);
                 p2 = coords2[p] + (i2 + i2);
                 dx = 0.0;
                 for (k = 0; k < 2; k++) {
                     dp0[k] = p1[k] - p0[k];
                     dp1[k] = p2[k] - p1[k];
                     dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
                 }
                 factor = 1.0/MAX(dx, small);
                 for (k = 0; k < 2; k++) {
                     dp0[k] *= factor;
                     dp1[k] *= factor;
                 }
                 crosses[i] = dp0[0] * dp1[1] - dp0[1] * dp1[0];
             }
             isconcave2[p] = 0;
             for (i = 0; i < nn2[p]; i++) {
                 i1 = (i + 1) % nn2[p];
                 if (crosses[i] * crosses[i1] < 0.0) {
                     isconcave2[p] = 1;
                     break;
                 }
             }
         }
         if (isconcave2[0] && isconcave2[1]) {
             printf("ERROR: both concave in concave_decompose0\n");
             return;
         }
         if (isconcave2[0]) {
             i_concave = 0;
             i_convex  = 1;
         }
         else {
             i_concave = 1;
             i_convex  = 0;
         }
*****/
         i_concave      = 0;
         i_convex       = 1;
 
         coords_part[0] = coords2[i_concave];  // tp be further worked on.
         nnode_part[0]  = nn2[i_concave];
         factors[0]     = 1.0;
         coords_part[1] = coords2[i_convex];
         nnode_part[1]  = nn2[i_convex];
         factors[1]     = myfactor;
 
         memcpy(coords_part[0] + (nnode_part[0] * dim), coords_part[0], (size_t)szdim);
         memcpy(coords_part[1] + (nnode_part[1] * dim), coords_part[1], (size_t)szdim);
     }
     free(coords);
     free(i0s);
     free(crosses);
     free(choose_i1);
 
     return;
  }
 
void concave_decompose0_minus(double *coords0, int nnode0,
                       int *npart, int *nnode_part,
                       double *factors, double *coords_part[2])
{
//   coords is (nnode + 1) long, and circular.
 
     int dim, szdim, nnode, nn, nn0, n1, n2, done, i1_save;
     int just_turned, islasti0, itest, i02[2];
     int ncave, i, j, k, i0, i1, i2, istart, iadd, nalloc;
     int *i0s, *mark, *choose_i1;
     double small, dx, factor, cross0, cross_last, cross;;
     double *coords, *crosses, *p0, *p1, *p2, dp0[2], dp1[2];
 
     small = 1.0e-12;
 
     dim    = 2;
     szdim  = dim * sizeof(double);
 
     nalloc = 64;
     i0s = (int *) malloc(nalloc * sizeof(int));
 
//   remove the middle point in any line segment.
 
     coords = (double *) malloc((nnode0 + 1) * szdim);
     mark   = (int    *) malloc( nnode0 * sizeof(int));
     for (i = 0; i < nnode0; i++) {
         mark[i] = 0;
     }
 
     i  = 0;
     i1 = 1;
     while (i1 < nnode0) {
         i2 = (i1 + 1) % nnode0;
         p0 = coords0 + (i  + i);
         p1 = coords0 + (i1 + i1);
         p2 = coords0 + (i2 + i2);
         dx = 0.0;
         for (k = 0; k < 2; k++) {
             dp0[k] = p1[k] - p0[k];
             dp1[k] = p2[k] - p1[k];
             dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
         }
         factor = 1.0/MAX(dx, small);
         for (k = 0; k < 2; k++) {
             dp0[k] *= factor;
             dp1[k] *= factor;
         }
         cross0 = dp0[0] * dp1[1] - dp0[1] * dp1[0];
         if (fabs(cross0) < small) {
             mark[i1] = 1;
         }
         else {
             i = i1;
         }
         i1++;
     }
     nnode = 0;
     for (i = 0; i < nnode0; i++) {
         if (mark[i]) continue;
         memcpy(coords + (nnode * dim), coords0 + (i * dim), (size_t)szdim);
         nnode++;
     }
     memcpy(coords + (nnode * dim), coords, (size_t)szdim);
     free(mark);
 
     for (i = 0; i < nnode; i++) {
         i1 = (i  + 1) % nnode;
         i2 = (i1 + 1) % nnode;
         p0 = coords + (i  + i);
         p1 = coords + (i1 + i1);
         p2 = coords + (i2 + i2);
         dx = 0.0;
         for (k = 0; k < 2; k++) {
             dp0[k] = p1[k] - p0[k];
             dp1[k] = p2[k] - p1[k];
             dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
         }
         factor = 1.0/MAX(dx, small);
         for (k = 0; k < 2; k++) {
             dp0[k] *= factor;
             dp1[k] *= factor;
         }
         cross0 = dp0[0] * dp1[1] - dp0[1] * dp1[0];
//       if (fabs(cross0) > small) {
         if (cross0 != 0.0) {
             istart = i;
             break;
         }
     }
     crosses = (double *) malloc(nnode * sizeof(double));
     for (i = 0; i < nnode; i++) {
         i0 = (i + istart) % nnode;
         i1 = (i0 + 1) % nnode;
         i2 = (i1 + 1) % nnode;
         p0 = coords + (i0 + i0);
         p1 = coords + (i1 + i1);
         p2 = coords + (i2 + i2);
         for (k = 0; k < 2; k++) {
             dp0[k] = p1[k] - p0[k];
             dp1[k] = p2[k] - p1[k];
             dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
         }
         factor = 1.0/MAX(dx, small);
         for (k = 0; k < 2; k++) {
             dp0[k] *= factor;
             dp1[k] *= factor;
         }
         cross = dp0[0] * dp1[1] - dp0[1] * dp1[0];
//       if (fabs(cross) < small) {
//           crosses[i0] = 0.0;
//       }
         crosses[i0] = cross;
     }
     choose_i1 = (int *) malloc(nnode * sizeof(int));
     for (i = 0; i < nnode; i++) {
         choose_i1[i] = -1;
     }
     just_turned = 0;
     cross_last  = cross0;
     for (i = 0; i <= nnode+1; i++) {
         i0 = (i  + istart) % nnode;
         i1 = (i0 + 1)      % nnode;
         cross = crosses[i0];
//       if (fabs(cross) < small) continue;
 
         if (cross_last * cross < 0.0) {
             if (just_turned) {
                 choose_i1[i0] = 1;
                 just_turned = 0;
             }
             else {
                 just_turned = 1;
             }
             cross_last = cross;
         }
         else {
             just_turned = 0;
         }
     }
     cross_last = cross0;
 
     i1_save = -1;
     ncave = 0;
     islasti0 = 0;
     for (i = 0; i <= nnode; i++) {
         i0 = (i  + istart) % nnode;
         i1 = (i0 + 1) % nnode;
         cross = crosses[i0];
 
//       if (fabs(cross) < small) continue;
 
         if (cross * cross_last < 0.0) {
             if (choose_i1[i1] >= 0) {
                 if (islasti0) {
                     if (ncave > 0) {
                         if (i0s[ncave-1] != i1) {
                             i0s[ncave] = i1;
                             ncave++;
                         }
                     }
                     else {
                         i0s[ncave] = i1;
                         ncave++;
                     }
                     islasti0   = 0;
                 }
                 else {
                     i1_save = i1;
                 }
             }
             else {
                 if (islasti0) {
                     if (ncave > 0) {
                         if (i0s[ncave-1] != i1) {
                             i0s[ncave] = i1;
                             ncave++;
                             islasti0 = 0;
                         }
                     }
                     else {
                         i0s[ncave] = i1;
                         ncave++;
                         islasti0 = 0;
                     }
                 }
                 else {
                     if (ncave > 0) {
                         if (i0s[(ncave)-1] != i0) {
                             i0s[ncave] = i0;
                             ncave++;
                             islasti0 = 1;
                         }
                     }
                     else {
                         i0s[ncave] = i0;
                         ncave++;
                         islasti0 = 1;
                     }
                 }
             }
             cross_last = cross;
             if (ncave == 2) break;
         }
     }
     if ((ncave == 0) && (i1_save >= 0)) {
         i0s[ncave] = i1;
         ncave++;
     }
     if (ncave == 1) {
//       itest  = (nnode -1 + istart) % nnode;
         itest = (nnode + istart) % (nnode + 1);
         i0s[ncave] =  itest;
         ncave++;
     }
     *npart = 0;
     if (ncave == 0) {
         *npart = 1;
         nnode_part[0]  = nnode;
         factors[0]     = 1.0;
         coords_part[0] = (double *) malloc((nnode + 1) * szdim);
         memcpy(coords_part[0], coords, (size_t)((nnode + 1) * szdim));
     }
     else if (ncave == 2) {
         *npart = 2;
         nn  = i0s[1] - i0s[0] + 1;
         nn0 = nnode - nn + 2;
         coords_part[0] = (double *) malloc((nn0 + 1) * szdim);
         coords_part[1] = (double *) malloc((nn  + 1) * szdim);
         nnode_part[0] = nn0;
         nnode_part[1] = nn;
 
         i0 = i0s[0];
         i1 = i0s[1];
         if (i0 > istart) {
             n1 = i0 - istart + 1;
         }
         else {
             n1  = MIN(i1, (i0 + nnode - istart + 1));
         }
         p0 = coords_part[0];
         for (i = 0; i < n1; i++) {
             i0 = (i + istart) % nnode;
             memcpy(p0 + (i * dim), coords + (i0 * dim), (size_t)szdim);
         }
         p0 = coords_part[0] + (dim * n1);
         n2 = nnode - n1 - (nn - 2);
         assert(n1 + n2 == nn0);
 
         for (i = 0; i < n2; i++) {
             i0 = (i + i1) % nnode;
             memcpy(p0 + (i * dim), coords + (i0 * dim), (size_t)szdim);
         }
         i0 = i0s[0];
         p1 = coords_part[1];
         for (i = 0; i < nn; i++) {
             i1 = (i + i0) % nnode;
             memcpy(p1 + (i * dim), coords + (i1 * dim), (size_t)szdim);
         }
         done = 0;
         for (i = 0; i < nn0; i++) {
             i1 = (i  + 1) % nn0;
             i2 = (i1 + 1) % nn0;
             p0 = coords_part[0] + (i  + i);
             p1 = coords_part[0] + (i1 + i1);
             p2 = coords_part[0] + (i2 + i2);
             dx = 0.0;
             for (k = 0; k < 2; k++) {
                 dp0[k] = p1[k] - p0[k];
                 dp1[k] = p2[k] - p1[k];
                 dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
             }
             factor = 1.0/MAX(dx, small);
             for (k = 0; k < 2; k++) {
                 dp0[k] *= factor;
                 dp1[k] *= factor;
             }
             cross0 = dp0[0] * dp1[1] - dp0[1] * dp1[0];
             if (cross0 != 0.0) {
                 done = 1;
                 break;
             }
         }
         assert(done);
         factors[0] = 1.0;
 
         done = 0;
         for (i = 0; i < nn; i++) {
             i1 = (i  + 1) % nn;
             i2 = (i1 + 1) % nn;
             p0 = coords_part[1] + (i  + i);
             p1 = coords_part[1] + (i1 + i1);
             p2 = coords_part[1] + (i2 + i2);
             dx = 0.0;
             for (k = 0; k < 2; k++) {
                 dp0[k] = p1[k] - p0[k];
                 dp1[k] = p2[k] - p1[k];
                 dx = MAX(MAX(dx, fabs(dp0[k])), fabs(dp0[0]));
             }
             factor = 1.0/MAX(dx, small);
             for (k = 0; k < 2; k++) {
                 dp0[k] *= factor;
                 dp1[k] *= factor;
             }
             cross = dp0[0] * dp1[1] - dp0[1] * dp1[0];
             if (cross * cross0 > 0.0) {
                 factors[1] = 1.0;
                 done = 1;
                 break;
             }
             else if (cross * cross0 < 0.0) {
                 factors[1] = -1.0;
                 done = 1;
                 break;
             }
         }
         assert(done);
         memcpy(coords_part[0] + (nn0 * dim), coords_part[0], (size_t)szdim);
         memcpy(coords_part[1] + (nn  * dim), coords_part[1], (size_t)szdim);
     }
     free(coords);
     free(i0s);
     free(crosses);
     free(choose_i1);
 
     return;
  }
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
void check_align(int dim, double *pt1, double *pt2, double *rads)
{
     int n0, i;
     double  s2, r2, tmp, rxy, rxyz, phi, theta;
     double  dr[3];
 
     aligned_sav  = 1;
     sinphi_sav   = 0.0;
     cosphi_sav   = 1.0;
     sintheta_sav = 0.0;
     costheta_sav = 1.0;
 
     n0 = 0;
     for (i = 0; i < dim; i++) {
         if (pt2[i] != pt1[i]) n0++;
         dr[i] = pt2[i] - pt1[i];
     }
     if (n0 <= 1) {
         aligned_sav = 1;
     }
     else {
         aligned_sav = 0;
     }
     if (!aligned_sav) {
         if (dim == 2) {
             s2 = 0.0;
             for (i = 0; i < 2; i++)  {
                 tmp = dr[i] * dr[i];
                 s2 += tmp;
             }
             rxy        = sqrt(s2);
             sinphi_sav = dr[1]/rxy;
             cosphi_sav = dr[0]/rxy;
         }
         else if (dim == 3) {
             s2 = 0.0;
             r2 = dr[2] * dr[2];
             for (i = 0; i < 2; i++)  {
                 tmp = dr[i] * dr[i];
                 s2 += tmp;
                 r2 += tmp;
             }
             rxy          = sqrt(s2);
             rxyz         = sqrt(r2);
             sintheta_sav = rxy/rxyz;
             theta        = asin(sintheta_sav);
             if (dr[2] < 0.0) {
                 theta = PI - theta;
             }
             costheta_sav = cos(theta);
             sinphi_sav   = dr[1]/rxy;
             phi          = asin(sinphi_sav);
             if ((dr[0] < 0.0) && (dr[1] < 0.0)) {
                 phi = PI + fabs(phi);
             }
             else if (dr[0] < 0.0) {
                 phi = PI - fabs(phi);
             }
             cosphi_sav = cos(phi);
 
             if (rads[0] > rads[2]) {    // 0-axis is the major axis.
                 major_axis_sav = 0;
                 tmp          =  sintheta_sav;
                 sintheta_sav =  -costheta_sav;
                 costheta_sav =  tmp;
             }
             else {
                 major_axis_sav = 2;
             }
             assert(rads[1] <= rads[major_axis_sav]);
         }
     }
     return;
 }
 
void get_direction(int dim, double *pt1, double *pt2)
{
     int     i;
     double  s2, r2, tmp, rxy, rxyz, phi, theta;
     double  dr[3];
 
     for (i = 0; i < dim; i++) {
         dr[i] = pt2[i] - pt1[i];
     }
     s2 = 0.0;
     r2 = dr[2] * dr[2];
     for (i = 0; i < 2; i++)  {
         tmp = dr[i] * dr[i];
         s2 += tmp;
         r2 += tmp;
     }
     rxy            = sqrt(s2);
     rxyz           = sqrt(r2);
     sintheta_conic = rxy/rxyz;
     theta          = asin(sintheta_conic);
     if (dr[2] < 0.0) {
         theta = PI - theta;
     }
     costheta_conic = cos(theta);
     sinphi_conic   = dr[1]/rxy;
     phi            = asin(sinphi_conic);
     if ((dr[0] < 0.0) && (dr[1] < 0.0)) {
         phi = PI + fabs(phi);
     }
     else if (dr[0] < 0.0) {
         phi = PI - fabs(phi);
     }
     cosphi_conic   = cos(phi);
 
     tmp            =  sintheta_conic;
     sintheta_conic =  -costheta_conic;
     costheta_conic =  tmp;
 
     return;
 }
 
 
void get_3d_face_description(int nreg, int *i2dto3d, double bottom, double top)
{
//   This function construct 3D faces from a 2D-palne, bootom, and top.
//   If necessary, it will decompose a concave polygons to a group of
//   convex polygons and construct 3D faces for the group.
 
     int    dim, szdim, ifile, nn_ea_face;
     int    ngrp, g, i, i1, i2, i3, j, k, k1, sz;
     int    nn0, nn, n, n1, n2, offset;
     int    nnode, nedge, nface;
     int    *nodelist_for_edge, *edgelist_for_face, *nodelist_for_face;
     int    *nnode_for_face, *facelist_for_zone;
     int    *edgelist, *nodelist;
     int    *nn_grp, nng;
     double *coords2d;
     double *xyz1, *xyz2, *coords, **xys_grp, **normal_grp, *xys;
     double *normal_of_face, norm[3], dp0[3], dp1[3], ctr[3];
     double tmp, vol;
     double *p0, *p1, *p2;
 
 
//   check xyz1 and xyz2 are two planes
 
     dim   = 3;
     szdim = dim * sizeof(double);
 
     nn_ea_face = 4;
 
     ifile  = -1;
     for (i = 0; i < nfile_oso; i++) {
         if (reg_oso[i] == nreg) {
             ifile = i;
             break;
         }
     }
     assert(ifile >= 0);
     ngrp = npart_oso[ifile];
 
     nnode_grp_oso[ifile]  = (int     *) malloc(ngrp * sizeof(int));
     nedge_grp_oso[ifile]  = (int     *) malloc(ngrp * sizeof(int));
     nface_grp_oso[ifile]  = (int     *) malloc(ngrp * sizeof(int));
     ctr_grp_oso[ifile]    = (double  *) malloc(ngrp * sizeof(double) * 3);
     coords_grp_oso[ifile] = (double **) malloc(ngrp * sizeof(double *));
 
     nodelist_for_edge_oso[ifile] = (int    **) malloc(ngrp * sizeof(int *));
     edgelist_for_face_oso[ifile] = (int    **) malloc(ngrp * sizeof(int *));
     nodelist_for_face_oso[ifile] = (int    **) malloc(ngrp * sizeof(int *));
     nnode_for_face_oso[ifile]    = (int    **) malloc(ngrp * sizeof(int *));
     normal_of_face_oso[ifile]    = (double **) malloc(ngrp * sizeof(double *));
 
 
     if (ngrp == 1) {
         coords2d = coords_part_oso[ifile][0];
         nn0      = nnode_part_oso[ifile][0];
         xyz1     = (double *) malloc((nn0 + 1) * szdim);
         xyz2     = (double *) malloc((nn0 + 1) * szdim);
 
         for (i = 0; i < nn0; i++) {
             i2 = i + i;
             i3 = i + i2;
             for (k = 0; k < 2; k++) {
                 k1 = i2dto3d[k] - 1;
                 xyz1[i3 + k1] = coords2d[i2 + k];
                 xyz2[i3 + k1] = coords2d[i2 + k];
             }
             k1 = i2dto3d[2] - 1;
             xyz1[i3 + k1] = bottom;
             xyz2[i3 + k1] = top;
         }
         nnode = nn0 + nn0;
         nedge = 3 * nn0;
         nface = nn0 + 2;
 
         nodelist_for_edge = (int *) malloc((nedge + nedge) * sizeof(int));
 
    //   bottom edges
         nodelist = nodelist_for_edge;
         for (i = 0; i < nn0; i++) {
             i1 = (i + 1) % nn0;
             nodelist[0] = i;
             nodelist[1] = i1;
             nodelist += 2;
         }
    //   top edge
         for (i = 0; i < nn0; i++) {
             i1 = (i + 1) % nn0;
             nodelist[0] = nn0 + i;
             nodelist[1] = nn0 + i1;
             nodelist += 2;
         }
    //   edge at side
         for (i = 0; i < nn0; i++) {
             nodelist[0] = i;
             nodelist[1] = i + nn0;
             nodelist += 2;
         }
    //   faces at side
 
         sz = (nface - 2) * nn_ea_face + 2 * nn0;
         edgelist_for_face = (int *) malloc(sz * sizeof(int));
         nodelist_for_face = (int *) malloc(sz * sizeof(int));
         nnode_for_face    = (int *) malloc(nface * sizeof(int));
 
         edgelist = edgelist_for_face;
         nodelist = nodelist_for_face;
         for (i = 0; i < nn0; i++) {
             i1 = (i + 1) % nn0;
             edgelist[0] = i;                   // bottom edge
             edgelist[1] = nn0 + nn0 + i1;  // side edge
             edgelist[2] = nn0 + i;           // top edge
             edgelist[3] = nn0 + nn0 + i;   // side edge
             edgelist += 4;
 
             nodelist[0] = i;
             nodelist[1] = i1;
             nodelist[2] = nn0 + i1;
             nodelist[3] = nn0 + i;
             nnode_for_face[i] = 4;
             nodelist += 4;
         }
    //   bottom face
         for (i = 0; i < nn0; i++) {
             i1 = nn0 + i;
             edgelist[i]  = i;
             edgelist[i1] = i1;
 
             nodelist[i]  = i;
             nodelist[i1] = i1;
         }
         nnode_for_face[nn0]   = nn0;
         nnode_for_face[nn0+1] = nn0;
 
         coords = (double *) malloc((nnode + 1) * szdim);
         sz = nn0 * szdim;
         memcpy(coords, xyz1, (size_t)sz);
         memcpy(coords + (dim * nn0), xyz2, (size_t)sz);
 
         facelist_for_zone = (int *) malloc(nface * sizeof(int));
         for (i = 0; i < nface; i++) {
             facelist_for_zone[i] = i;
         }
         cal_ctr0(nface, nnode, facelist_for_zone, nnode_for_face,
                  nodelist_for_face, coords, ctr, &vol);
 
         free(facelist_for_zone);
 
    //   calculate normals of each face
 
         normal_of_face = (double *) malloc(nface * szdim);
 
         offset = 0;
         for (k = 0; k < nface; k++) {
             nn = nnode_for_face[k];
             nodelist = nodelist_for_face + offset;
             i  = 0;
             i1 = (i + 1) % nn;
             i2 = (i1+ 1) % nn;
             n = nodelist[i];
             n1 = nodelist[i1];
             n2 = nodelist[i2];
             p0 = coords + (n  * dim);
             p1 = coords + (n1 * dim);
             p2 = coords + (n2 * dim);
             for (i = 0; i < dim; i++) {
                 dp0[i] = p1[i] - p0[i];
                 dp1[i] = p2[i] - p1[i];
             }
             norm[0] = dp0[1] * dp1[2] - dp0[2] * dp1[1];
             norm[1] = dp0[2] * dp1[0] - dp0[0] * dp1[2];
             norm[2] = dp0[0] * dp1[1] - dp0[1] * dp1[0];
             tmp = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
             tmp = 1.0/sqrt(tmp);
             for (j = 0; j < dim; j++) {
                 norm[j] *= tmp;
             }
             memcpy(normal_of_face + (k * dim), norm, (size_t)szdim);
 
             offset += nn;
         }
         nnode_grp_oso[ifile][0] = nnode;
         nedge_grp_oso[ifile][0] = nedge;
         nface_grp_oso[ifile][0] = nface;
         memcpy(ctr_grp_oso[ifile], ctr, (size_t)szdim);
 
         nodelist_for_edge_oso[ifile][0] = nodelist_for_edge;
         edgelist_for_face_oso[ifile][0] = edgelist_for_face;
         nodelist_for_face_oso[ifile][0] = nodelist_for_face;
         nnode_for_face_oso[ifile][0]    = nnode_for_face;
         coords_grp_oso[ifile][0]        = coords;
 
         normal_of_face_oso[ifile][0]    = normal_of_face;
 
         free(xyz1);
         free(xyz2);
     }    // endof ngrp = 1
     else if (ngrp > 1) {
 
         xys_grp  = coords_part_oso[ifile];
         nn_grp   = nnode_part_oso[ifile];
         normal_grp  = normal_of_face_oso[ifile];
 
         for (g = 0; g < ngrp; g++) {
 
             xys = xys_grp[g];
             nng  = nn_grp[g];
 
             xyz1 = (double *) malloc((nng + 1) * szdim);
             xyz2 = (double *) malloc((nng + 1) * szdim);
 
             for (i = 0; i < nng; i++) {
                 i2 = i + i;
                 i3 = i + i2;
                 for (k = 0; k < 2; k++) {
                     k1 = i2dto3d[k] - 1;
                     xyz1[i3 + k1] = xys[i2 + k];
                     xyz2[i3 + k1] = xys[i2 + k];
                 }
                 k1 = i2dto3d[2] - 1;
                 xyz1[i3 + k1] = bottom;
                 xyz2[i3 + k1] = top;
             }
             nnode = nng + nng;
             nedge = 3 * nng;
             nface = nng + 2;
 
             nodelist_for_edge = (int *) malloc((nedge + nedge) * sizeof(int));
 
        //   bottom edges
             nodelist = nodelist_for_edge;
             for (i = 0; i < nng; i++) {
                 i1 = (i + 1) % nng;
                 nodelist[0] = i;
                 nodelist[1] = i1;
                 nodelist += 2;
             }
        //   top edge
             for (i = 0; i < nng; i++) {
                 i1 = (i + 1) % nng;
                 nodelist[0] = nng + i;
                 nodelist[1] = nng + i1;
                 nodelist += 2;
             }
        //   edge at side
             for (i = 0; i < nng; i++) {
                 nodelist[0] = i;
                 nodelist[1] = i + nng;
                 nodelist += 2;
             }
        //   faces at side
 
             sz = (nface - 2) * nn_ea_face + 2 * nng;
             edgelist_for_face = (int *) malloc(sz * sizeof(int));
             nodelist_for_face = (int *) malloc(sz * sizeof(int));
             nnode_for_face = (int *) malloc(nface * sizeof(int));
 
             edgelist = edgelist_for_face;
             nodelist = nodelist_for_face;
             for (i = 0; i < nng; i++) {
                 i1 = (i + 1) % nng;
                 edgelist[0] = i;                   // bottom edge
                 edgelist[1] = nng + nng + i1;  // side edge
                 edgelist[2] = nng + i;           // top edge
                 edgelist[3] = nng + nng + i;   // side edge
                 edgelist += 4;
 
                 nodelist[0] = i;
                 nodelist[1] = i1;
                 nodelist[2] = nng + i1;
                 nodelist[3] = nng + i;
                 nnode_for_face[i] = 4;
                 nodelist += 4;
             }
        //   bottom face
             for (i = 0; i < nng; i++) {
                 i1 = nng + i;
                 edgelist[i]  = i;
                 edgelist[i1] = i1;
 
                 nodelist[i]  = i;
                 nodelist[i1] = i1;
             }
             nnode_for_face[nng]   = nng;
             nnode_for_face[nng+1] = nng;
 
             coords = (double *) malloc((nnode + 1) * szdim);
             sz = nng * szdim;
             memcpy(coords, xyz1, (size_t)sz);
             memcpy(coords + (dim * nng), xyz2, (size_t)sz);
 
             facelist_for_zone = (int *) malloc(nface * sizeof(int));
             for (i = 0; i < nface; i++) {
                 facelist_for_zone[i] = i;
             }
             cal_ctr0(nface, nnode, facelist_for_zone, nnode_for_face,
                      nodelist_for_face, coords, ctr, &vol);
 
             free(facelist_for_zone);
 
        //   calculate normals of each face
 
             normal_of_face = (double *) malloc(nface * szdim);
 
             offset = 0;
             for (k = 0; k < nface; k++) {
                 nn = nnode_for_face[k];
                 nodelist = nodelist_for_face + offset;
                 i  = 0;
                 i1 = (i + 1) % nn;
                 i2 = (i1+ 1) % nn;
                 n = nodelist[i];
                 n1 = nodelist[i1];
                 n2 = nodelist[i2];
                 p0 = coords + (n  * dim);
                 p1 = coords + (n1 * dim);
                 p2 = coords + (n2 * dim);
                 for (i = 0; i < dim; i++) {
                     dp0[i] = p1[i] - p0[i];
                     dp1[i] = p2[i] - p1[i];
                 }
                 norm[0] = dp0[1] * dp1[2] - dp0[2] * dp1[1];
                 norm[1] = dp0[2] * dp1[0] - dp0[0] * dp1[2];
                 norm[2] = dp0[0] * dp1[1] - dp0[1] * dp1[0];
                 tmp = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
                 tmp = 1.0/sqrt(tmp);
                 for (j = 0; j < dim; j++) {
                     norm[j] *= tmp;
                 }
                 memcpy(normal_of_face + (k * dim), norm, (size_t)szdim);
 
                 offset += nn;
             }
             nnode_grp_oso[ifile][g] = nnode;
             nedge_grp_oso[ifile][g] = nedge;
             nface_grp_oso[ifile][g] = nface;
             memcpy(ctr_grp_oso[ifile] + 3 * g, ctr, (size_t)szdim);
 
             nodelist_for_edge_oso[ifile][g] = nodelist_for_edge;
             edgelist_for_face_oso[ifile][g] = edgelist_for_face;
             nodelist_for_face_oso[ifile][g] = nodelist_for_face;
             nnode_for_face_oso[ifile][g]    = nnode_for_face;
             coords_grp_oso[ifile][g]        = coords;
             normal_grp[g]                   = normal_of_face;
 
             free(xyz1);
             free(xyz2);
 
 
         }   // end of grp
 
     }    // end of ngrp > 1
 
     return;
 }
void check_plane(int dim, double *xyz1, double *xyz2, int nn0, int *ifcyl,
                 int *nnode, int *nedge, int *nface,
                 int *nodelist_for_edge, int *edgelist_for_face,
                 int *nedge_for_face, int *nodelist_for_face,
                 double *coords, double *ctr, double *normal_of_face)
{
//   ifcyl:   input and output
//            = 0 in input: check whether the ending faces are planes or not.
//           != 0 in input: don't check whether they are planes or not.
//            = 0 in output: results are not planes.
//            > 0 in output, planes and cylinder along ifcyl direction (1-based)
//            < 0 in output, planes, but conic along -ifcyl direction (1-based).
//
     int    i, i1, i2, j, k, k1, plane, isplane, same, cyl, sz, szdim;
     int    nn, n, n1, n2, offset;
     int    *edgelist, *nodelist;
     int    *facelist_for_zone;
     double norm[3], dp0[3], dp1[3], *ptrs[2];
     double tmp, dot, vol, dxmax, dx, small;
     double *p0, *p1, *p2;
 
//   check xyz1 and xyz2 are two planes
 
     szdim = dim * sizeof(double);
 
     ptrs[0] = xyz1;
     ptrs[1] = xyz2;
 
     if (*ifcyl == 0) {
 
         for (plane = 0; plane < 2; plane++) {
             p0 = ptrs[plane];
             p1 = p0 + dim;
             p2 = p1 + dim;
             for (i = 0; i < dim; i++) {
                 dp0[i] = p1[i] - p0[i];
                 dp1[i] = p2[i] - p1[i];
             }
             norm[0] = dp0[1] * dp1[2] - dp0[2] * dp1[1];
             norm[1] = dp0[2] * dp1[0] - dp0[0] * dp1[2];
             norm[2] = dp0[0] * dp1[1] - dp0[1] * dp1[0];
             tmp = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
             tmp = 1.0/sqrt(tmp);
             for (i = 0; i < dim; i++) {
                 norm[i] *= tmp;
             }
             dxmax = 0.0;
             for (i = 0; i < nn0; i++) {
                 i1 = (i + 1) % dim;
                 p0 = ptrs[plane] + (i  * dim);
                 p1 = ptrs[plane] + (i1 * dim);
                 for (k = 0; k < dim; k++) {
                     dx = fabs(p1[k] - p0[k]);
                     if (dx > dxmax) dxmax = dx;
                 }
             }
             small = 1.0e-12 * dxmax;
 
             isplane = 1;
             for (i = 0; i < nn0; i++) {
                 i1 = (i + 1) % dim;
                 p0 = ptrs[plane] + (i  * dim);
                 p1 = ptrs[plane] + (i1 * dim);
                 dot = 0.0;
                 for (k = 0; k < dim; k++) {
                      dot += ((p1[k] - p0[k]) * norm[k]); ;
                 }
                 if (fabs(dot) > small) {
                     isplane = 0;
                     break;
                 }
             }
             assert(isplane);
         }
         for (k = 0; k < dim; k++) {
             same = 1;
             for (i = 1; i < nn0; i++) {
                 p1 = xyz1 + (i * dim);
                 p2 = xyz2 + (i * dim);
                 if ((p1[k] != xyz1[k]) || (p2[k] != xyz2[k])) {
                     same = 0;
                     break;
                 }
             }
             if (same) {
                 cyl = 1;
                 for (i = 0; i < nn0; i++) {
                     p1 = xyz1 + (i * dim);
                     p2 = xyz2 + (i * dim);
                     for (k1 = 0; k1 < dim; k1++) {
                         if (k1 == k) continue;
                         if (p1[k1] != p2[k1]) {
                             cyl = 0;
                             break;
                         }
                     }
                 }
                 if (cyl) {
                     *ifcyl = k + 1;
                 }
                 else {
                     *ifcyl = -(k + 1);
                 }
                 break;
             }
             if (*ifcyl != 0) break;
         }
     }
     *nnode = nn0 + nn0;
     *nedge = 3 * nn0;
     *nface = nn0 + 2;
 
//   bottom edges
     nodelist = nodelist_for_edge;
     for (i = 0; i < nn0; i++) {
         i1 = (i + 1) % nn0;
         nodelist[0] = i;
         nodelist[1] = i1;
         nodelist += 2;
     }
//   top edge
     for (i = 0; i < nn0; i++) {
         i1 = (i + 1) % nn0;
         nodelist[0] = nn0 + i;
         nodelist[1] = nn0 + i1;
         nodelist += 2;
     }
//   edge at side
     for (i = 0; i < nn0; i++) {
         nodelist[0] = i;
         nodelist[1] = i + nn0;
         nodelist += 2;
     }
//   faces at side
 
     edgelist = edgelist_for_face;
     nodelist = nodelist_for_face;
     for (i = 0; i < nn0; i++) {
         i1 = (i + 1) % nn0;
         edgelist[0] = i;                   // bottom edge
         edgelist[1] = nn0 + nn0 + i1;  // side edge
         edgelist[2] = nn0 + i;           // top edge
         edgelist[3] = nn0 + nn0 + i;   // side edge
         nedge_for_face[i] = 4;
         edgelist += 4;
 
         nodelist[0] = i;
         nodelist[1] = i1;
         nodelist[2] = nn0 + i1;
         nodelist[3] = nn0 + i;
         nodelist += 4;
     }
//   bottom face
     for (i = 0; i < nn0; i++) {
         i1 = nn0 + i;
         edgelist[i]  = i;
         edgelist[i1] = i1;
 
         nodelist[i]  = i;
         nodelist[i1] = i1;
     }
     nedge_for_face[nn0]   = nn0;
     nedge_for_face[nn0+1] = nn0;
 
     sz = nn0 * dim * sizeof(double);
     memcpy(coords, xyz1, (size_t)sz);
     memcpy(coords + (dim * nn0), xyz2, (size_t)sz);
 
     facelist_for_zone = (int *) malloc((*nface) * sizeof(int));
     for (i = 0; i < *nface; i++) {
         facelist_for_zone[i] = i;
     }
     cal_ctr0(*nface, *nnode, facelist_for_zone, nedge_for_face,
              nodelist_for_face, coords, ctr, &vol);
     free(facelist_for_zone);
 
//   calculate normals of each face
 
     offset = 0;
     for (k = 0; k < *nface; k++) {
         nn = nedge_for_face[k];
         nodelist = nodelist_for_face + offset;
         i  = 0;
         i1 = (i + 1) % nn;
         i2 = (i1+ 1) % nn;
         n = nodelist[i];
         n1 = nodelist[i1];
         n2 = nodelist[i2];
         p0 = coords + (n  * dim);
         p1 = coords + (n1 * dim);
         p2 = coords + (n2 * dim);
         for (i = 0; i < dim; i++) {
             dp0[i] = p1[i] - p0[i];
             dp1[i] = p2[i] - p1[i];
         }
         norm[0] = dp0[1] * dp1[2] - dp0[2] * dp1[1];
         norm[1] = dp0[2] * dp1[0] - dp0[0] * dp1[2];
         norm[2] = dp0[0] * dp1[1] - dp0[1] * dp1[0];
         tmp = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
         tmp = 1.0/sqrt(tmp);
         for (j = 0; j < dim; j++) {
             norm[j] *= tmp;
         }
         memcpy(normal_of_face + (k * dim), norm, (size_t)szdim);
 
         offset += nn;
     }
 
     return;
 }
 
void gconic_rec(int ifinquiry, int *mixed, int ifcyl,
             int dim, double r0, double r1, double *ctr0, double *ctr1,
             double *xl, double *dx, double *vol)
{
     int myifcyl, i, i1, szdim;
     double tmp, myctr0[3], myctr1[3], myxl[3], mydx[3];
 
     if (dim == 3) {
         if ((ifcyl == 2) || (ifcyl == -2)) {
             for (i = 0; i < dim; i++) {
                 i1 = (i + 1) % dim;
                 myctr0[i] = ctr0[i1];
                 myctr1[i] = ctr1[i1];
                 myxl[i] = xl[i1];
                 mydx[i] = dx[i1];
             }
         }
         else if ((ifcyl == 3) || (ifcyl == -3)) {
             for (i = 0; i < dim; i++) {
                 i1 = (i + 2) % dim;
                 myctr0[i] = ctr0[i1];
                 myctr1[i] = ctr1[i1];
                 myxl[i] = xl[i1];
                 mydx[i] = dx[i1];
             }
         }
         else {
             szdim = dim * sizeof(double);
             memcpy(myctr0, ctr0, (size_t)szdim);
             memcpy(myctr1, ctr1, (size_t)szdim);
             memcpy(myxl,   xl,   (size_t)szdim);
             memcpy(mydx,   dx,   (size_t)szdim);
         }
     }
     else {
         szdim = dim * sizeof(double);
         memcpy(myctr0, ctr0, (size_t)szdim);
         memcpy(myctr1, ctr1, (size_t)szdim);
         memcpy(myxl,   xl,   (size_t)szdim);
         memcpy(mydx,   dx,   (size_t)szdim);
     }
     if ((dim == 3) && (ifcyl == 0)) {
         gconic_rec_rotate(ifinquiry, mixed, dim, r0, r1, myctr0, myctr1, myxl, mydx, vol);
     }
     else {
         if (ifcyl < 0) {
             myifcyl = -1;
         }
         else {
             myifcyl = 1;
             if (myctr0[0] > myctr1[0]) {
                 tmp = myctr1[0];
                 myctr1[0] = myctr0[0];
                 myctr0[0] = tmp;
             }
         }
         gconic_rec_cyl(ifinquiry, mixed, myifcyl, dim, r0, r1, myctr0, myctr1, myxl, mydx, vol);
     }
     return;
 }
 
//   gellipse_rec_new
void gellipse_rec(int ifinquiry, int *mixed, int geop, double vcell, int dim,
                  double *ctr, double *r, double *xl, double *dx, double *vol)
{
//   Without 3D shrinking.
//   2D shrink
//   dr is the direction of the ellipse.
 
     int    i0,  i, k, n0, szdim, algned;
     double s2, r2, rxy, rxyz, r0, factor, tmp, dxcell;
     double theta, phi, sint, sinp, cost, cosp;
     double xr[3], myxl[3], myxr[3], ctrp[3], myctr[3], mydx[3];
     double c4[4][2], cp4[4][2], c8[8][3], cp8[8][3];
     double ctroid[2];
     double *cp, *c;
 
     int    nface, nedge, nnode;
     int    nedge_for_face[6], edgelist_for_face[24], nodelist_for_edge[48];
 
     szdim  = dim * sizeof(double);
     dxcell = dx[0];
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     if (!aligned_sav) {
         sinp = sinphi_sav;
         cosp = cosphi_sav;
         sint = sintheta_sav;
         cost = costheta_sav;
 
         for (i = 0; i < dim; i++) {
             myxl[i] = xl[i] - ctr[i];
             myxr[i] = xr[i] - ctr[i];
         }
         if (dim == 2) {
 
             c4[0][0] = myxl[0];
             c4[0][1] = myxl[1];
 
             c4[1][0] = myxr[0];
             c4[1][1] = myxl[1];
 
             c4[2][0] = myxr[0];
             c4[2][1] = myxr[1];
 
             c4[3][0] = myxl[0];
             c4[3][1] = myxr[1];
 
             for (i = 0; i < 4; i++) {
                 c  = c4[i];
                 cp = cp4[i];
                 cp[0] =  cosp * c[0] + sinp * c[1];
                 cp[1] = -sinp * c[0] + cosp * c[1];
             }
//             ctrp[0] =  cosp * ctr[0] + sinp * ctr[1];
//             ctrp[1] = -sinp * ctr[0] + cosp * ctr[1];
 
//           scaling
 
             r0 = r[0];
             for (i = 1; i < dim; i++) {
                 if (r[i] < r0) {
                     r0 = r[i];
                 }
             }
             factor = 1.0;
             for (i = 0; i < dim; i++) {
                 tmp     = r[i]/r0;
                 factor *= tmp;
                 tmp     = 1.0/tmp;
                 for (k = 0; k < 4; k++) {
                     cp4[k][i] *= tmp;
                 }
//                 ctrp[i] *= tmp;
             }
//           rotate back
             sinp = - sinp;
 
             for (i = 0; i < 4; i++) {
                 c  = c4[i];
                 cp = cp4[i];
                 c[0] =  cosp * cp[0] + sinp * cp[1];
                 c[1] = -sinp * cp[0] + cosp * cp[1];
             }
//             myctr[0] =  cosp * ctrp[0] + sinp * ctrp[1];
//             myctr[1] = -sinp * ctrp[0] + cosp * ctrp[1];
 
             myctr[0] = 0.0;
             myctr[1] = 0.0;
 
             sph_poly2d(ifinquiry, geop, dim, myctr, r0, c4[0], 4, mixed, vol, ctroid);
 
             if (!ifinquiry) {
                 if (*mixed) {
                     *vol *= factor;
                 }
                 else {
                     *vol = vcell;
                 }
                 assert(*vol >= 0.0);
                 assert(*vol <= vcell + small);
             }
             else if (*vol != 0.0) {
                 *vol = vcell;
             }
         }
         else if (dim == 3) {
 
             c8[0][0] = myxl[0];
             c8[0][1] = myxl[1];
             c8[0][2] = myxl[2];
 
             c8[1][0] = myxr[0];
             c8[1][1] = myxl[1];
             c8[1][2] = myxl[2];
 
             c8[2][0] = myxr[0];
             c8[2][1] = myxr[1];
             c8[2][2] = myxl[2];
 
             c8[3][0] = myxl[0];
             c8[3][1] = myxr[1];
             c8[3][2] = myxl[2];
 
             c8[4][0] = myxl[0];
             c8[4][1] = myxl[1];
             c8[4][2] = myxr[2];
 
             c8[5][0] = myxr[0];
             c8[5][1] = myxl[1];
             c8[5][2] = myxr[2];
 
             c8[6][0] = myxr[0];
             c8[6][1] = myxr[1];
             c8[6][2] = myxr[2];
 
             c8[7][0] = myxl[0];
             c8[7][1] = myxr[1];
             c8[7][2] = myxr[2];
 
             for (i = 0; i < 8; i++) {
                 c  = c8[i];
                 cp = cp8[i];
                 cp[0] =  cost * cosp * c[0] + cost * sinp * c[1] - sint * c[2];
                 cp[1] = -       sinp * c[0] +        cosp * c[1];
                 cp[2] =  sint * cosp * c[0] + sint * sinp * c[1] + cost * c[2];
             }
//             ctrp[0] =  cost * cosp * ctr[0] + cost * sinp * ctr[1] - sint * ctr[2];
//             ctrp[1] = -       sinp * ctr[0] +        cosp * ctr[1];
//             ctrp[2] =  sint * cosp * ctr[0] + sint * sinp * ctr[1] + cost * ctr[2];
 
//             tmp = (fabs(ctrp[0]) + fabs(ctrp[1]) + fabs(ctrp[2]))/dx[0];
//             assert(tmp < small);
//
             ctrp[0] = 0.0;
             ctrp[1] = 0.0;
             ctrp[2] = 0.0;
 
             nface = 6;
             nedge = 12;
             nnode = 8;
             for (i = 0; i < nface; i++) {
                 nedge_for_face[i] = 4;
             }
             // x-edge
             nodelist_for_edge[0] = 0;
             nodelist_for_edge[1] = 1;
 
             nodelist_for_edge[2] = 2;
             nodelist_for_edge[3] = 3;
 
             nodelist_for_edge[4] = 6;
             nodelist_for_edge[5] = 7;
 
             nodelist_for_edge[6] = 4;
             nodelist_for_edge[7] = 5;
             // y-edges
             nodelist_for_edge[8]  = 0;
             nodelist_for_edge[9]  = 3;
 
             nodelist_for_edge[10] = 1;
             nodelist_for_edge[11] = 2;
 
             nodelist_for_edge[12] = 5;
             nodelist_for_edge[13] = 6;
 
             nodelist_for_edge[14] = 4;
             nodelist_for_edge[15] = 7;
             // z-edges
             nodelist_for_edge[16] = 0;
             nodelist_for_edge[17] = 4;
 
             nodelist_for_edge[18] = 1;
             nodelist_for_edge[19] = 5;
 
             nodelist_for_edge[20] = 2;
             nodelist_for_edge[21] = 6;
 
             nodelist_for_edge[22] = 3;
             nodelist_for_edge[23] = 7;
 
             // x-faces
             edgelist_for_face[0] = 4;
             edgelist_for_face[1] = 11;
             edgelist_for_face[2] = 7;
             edgelist_for_face[3] = 8;
 
             edgelist_for_face[4] = 5;
             edgelist_for_face[5] = 9;
             edgelist_for_face[6] = 6;
             edgelist_for_face[7] = 10;
 
             // y-faces
             edgelist_for_face[8]  = 0;
             edgelist_for_face[9]  = 8;
             edgelist_for_face[10] = 3;
             edgelist_for_face[11] = 9;
 
             edgelist_for_face[12] = 1;
             edgelist_for_face[13] = 10;
             edgelist_for_face[14] = 2;
             edgelist_for_face[15] = 11;
 
             // z-faces
             edgelist_for_face[16] = 0;
             edgelist_for_face[17] = 5;
             edgelist_for_face[18] = 1;
             edgelist_for_face[19] = 4;
 
             edgelist_for_face[20] = 2;
             edgelist_for_face[21] = 6;
             edgelist_for_face[22] = 3;
             edgelist_for_face[23] = 7;
 
             ellipse_poly3d(ifinquiry, dim, ctrp, r, cp8[0], nnode, nedge, nface,
                            nedge_for_face, nodelist_for_edge, edgelist_for_face,
                            vcell, mixed, vol);
 
             assert(*vol >= 0.0);
             assert(*vol <= vcell + small);
         }
     }
     else {
         memcpy(ctrp, ctr, (size_t)szdim);
 
         r0 = r[0];
         for (i = 1; i < dim; i++) {
             if (r[i] < r0) {
                 r0 = r[i];
             }
         }
         factor = 1.0;
         for (i = 0; i < dim; i++) {
             tmp      = r[i]/r0;
             factor  *= tmp;
             tmp      = 1.0/tmp;
             myxl[i]  = xl[i]  * tmp;
             mydx[i]  = dx[i]  * tmp;
             myctr[i] = ctr[i] * tmp;
         }
         gsph_rec(ifinquiry, mixed, geop, dim, myctr, r0, myxl, mydx, vol);
 
         if (!ifinquiry) {
             *vol *= factor;
 
             assert(*vol >= 0.0);
             assert(*vol <= vcell + small);
         }
         else if (*vol != 0.0) {
             *vol = vcell;
         }
     }
 
     return;
 }
 
 
 
 
// gellipse_rec_old
void gellipse_rec_old(int ifinquiry, int *mixed, int geop, double vcell, int dim,
                  double *ctr, double *r, double *xl, double *dx, double *vol)
{
//   With 3D shrinking
//   3D shrink
//   dr is the direction of the ellipse.
 
     int    i0,  i, k, n0, szdim, algned;
     double s2, r2, rxy, rxyz, r0, factor, tmp;
     double theta, phi, sint, sinp, cost, cosp;
     double xr[3], ctrp[3], myctr[3], myxl[3], mydx[3];
     double c4[4][2], cp4[4][2], c8[8][3], cp8[8][3];
     double ctroid[2];
     double *cp, *c;
 
     int    nface, nedge, nnode;
     int    nedge_for_face[6], edgelist_for_face[24], nodelist_for_edge[48];
 
     szdim = dim * sizeof(double);
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     if (!aligned_sav) {
         sinp = sinphi_sav;
         cosp = cosphi_sav;
         sint = sintheta_sav;
         cost = costheta_sav;
 
         if (dim == 2) {
             c4[0][0] = xl[0];
             c4[0][1] = xl[1];
 
             c4[1][0] = xr[0];
             c4[1][1] = xl[1];
 
             c4[2][0] = xr[0];
             c4[2][1] = xr[1];
 
             c4[3][0] = xl[0];
             c4[3][1] = xr[1];
 
             for (i = 0; i < 4; i++) {
                 c  = c4[i];
                 cp = cp4[i];
                 cp[0] =  cosp * c[0] + sinp * c[1];
                 cp[1] = -sinp * c[0] + cosp * c[1];
             }
             ctrp[0] =  cosp * ctr[0] + sinp * ctr[1];
             ctrp[1] = -sinp * ctr[0] + cosp * ctr[1];
 
//           scaling
 
             r0 = r[0];
             for (i = 1; i < dim; i++) {
                 if (r[i] < r0) {
                     r0 = r[i];
                 }
             }
             factor = 1.0;
             for (i = 0; i < dim; i++) {
                 tmp     = r[i]/r0;
                 factor *= tmp;
                 tmp     = 1.0/tmp;
                 for (k = 0; k < 4; k++) {
                     cp4[k][i] *= tmp;
                 }
                 ctrp[i] *= tmp;
             }
//           rotate back
             sinp = - sinp;
             for (i = 0; i < 4; i++) {
                 c  = c4[i];
                 cp = cp4[i];
                 c[0] =  cosp * cp[0] + sinp * cp[1];
                 c[1] = -sinp * cp[0] + cosp * cp[1];
             }
             myctr[0] =  cosp * ctrp[0] + sinp * ctrp[1];
             myctr[1] = -sinp * ctrp[0] + cosp * ctrp[1];
 
             sph_poly2d(ifinquiry, geop, dim, myctr, r0, c4[0], 4, mixed, vol, ctroid);
 
             if (!ifinquiry) {
                 if (*mixed) {
                     *vol *= factor;
                 }
                 else {
                     *vol = vcell;
                 }
                 assert(*vol >= 0.0);
                 assert(*vol <= vcell + small);
             }
             else if (*vol != 0.0) {
                 *vol = vcell;
             }
         }
         else if (dim == 3) {
 
             c8[0][0] = xl[0];
             c8[0][1] = xl[1];
             c8[0][2] = xl[2];
 
             c8[1][0] = xr[0];
             c8[1][1] = xl[1];
             c8[1][2] = xl[2];
 
             c8[2][0] = xr[0];
             c8[2][1] = xr[1];
             c8[2][2] = xl[2];
 
             c8[3][0] = xl[0];
             c8[3][1] = xr[1];
             c8[3][2] = xl[2];
 
             c8[4][0] = xl[0];
             c8[4][1] = xl[1];
             c8[4][2] = xr[2];
 
             c8[5][0] = xr[0];
             c8[5][1] = xl[1];
             c8[5][2] = xr[2];
 
             c8[6][0] = xr[0];
             c8[6][1] = xr[1];
             c8[6][2] = xr[2];
 
             c8[7][0] = xl[0];
             c8[7][1] = xr[1];
             c8[7][2] = xr[2];
 
             for (i = 0; i < 8; i++) {
                 c  = c8[i];
                 cp = cp8[i];
                 cp[0] =  cost * cosp * c[0] + cost * sinp * c[1] - sint * c[2];
                 cp[1] = -       sinp * c[0] +        cosp * c[1];
                 cp[2] =  sint * cosp * c[0] + sint * sinp * c[1] + cost * c[2];
             }
             ctrp[0] =  cost * cosp * ctr[0] + cost * sinp * ctr[1] - sint * ctr[2];
             ctrp[1] = -       sinp * ctr[0] +        cosp * ctr[1];
             ctrp[2] =  sint * cosp * ctr[0] + sint * sinp * ctr[1] + cost * ctr[2];
             tmp     = (fabs(ctrp[1]) + fabs(ctrp[2]))/dx[0];
             assert(tmp < small);
             ctrp[1] = 0.0;
             ctrp[2] = 0.0;
 
             r0 = r[0];
             for (i = 1; i < dim; i++) {
                 if (r[i] < r0) {
                     r0 = r[i];
                 }
             }
//           scaling
 
             factor = 1.0;
             for (i = 0; i < dim; i++) {
                 tmp     = r[i]/r0;
                 factor *= tmp;
                 tmp     = 1.0/tmp;
                 for (k = 0; k < 8; k++) {
                     cp8[k][i] *= tmp;
                 }
                 ctrp[i] *= tmp;
             }
 
/****
//           rotate back
//           it should be a different angle
 
             sint = -sint;
             sinp = -sinp;
             for (i = 0; i < 8; i++) {
                 c  =  c8[i];
                 cp = cp8[i];
                 c[0] =  cost * cosp * cp[0] + cost * sinp * cp[1] - sint * cp[2];
                 c[1] = -       sinp * cp[0] +        cosp * cp[1];
                 c[2] =  sint * cosp * cp[0] + sint * sinp * cp[1] + cost * cp[2];
             }
             myctr[0] =  cost * cosp * ctrp[0] + cost * sinp * ctrp[1] - sint * ctrp[2];
             myctr[1] = -       sinp * ctrp[0] +        cosp * ctrp[1];
             myctr[2] =  sint * cosp * ctrp[0] + sint * sinp * ctrp[1] + cost * ctrp[2];
***/
 
             nface = 6;
             nedge = 12;
             nnode = 8;
             for (i = 0; i < nface; i++) {
                 nedge_for_face[i] = 4;
             }
             // x-edge
             nodelist_for_edge[0] = 0;
             nodelist_for_edge[1] = 1;
 
             nodelist_for_edge[2] = 2;
             nodelist_for_edge[3] = 3;
 
             nodelist_for_edge[4] = 6;
             nodelist_for_edge[5] = 7;
 
             nodelist_for_edge[6] = 4;
             nodelist_for_edge[7] = 5;
             // y-edges
             nodelist_for_edge[8]  = 0;
             nodelist_for_edge[9]  = 3;
 
             nodelist_for_edge[10] = 1;
             nodelist_for_edge[11] = 2;
 
             nodelist_for_edge[12] = 5;
             nodelist_for_edge[13] = 6;
 
             nodelist_for_edge[14] = 4;
             nodelist_for_edge[15] = 7;
             // z-edges
             nodelist_for_edge[16] = 0;
             nodelist_for_edge[17] = 4;
 
             nodelist_for_edge[18] = 1;
             nodelist_for_edge[19] = 5;
 
             nodelist_for_edge[20] = 2;
             nodelist_for_edge[21] = 6;
 
             nodelist_for_edge[22] = 3;
             nodelist_for_edge[23] = 7;
 
             // x-faces
             edgelist_for_face[0] = 4;
             edgelist_for_face[1] = 11;
             edgelist_for_face[2] = 7;
             edgelist_for_face[3] = 8;
 
             edgelist_for_face[4] = 5;
             edgelist_for_face[5] = 9;
             edgelist_for_face[6] = 6;
             edgelist_for_face[7] = 10;
 
             // y-faces
             edgelist_for_face[8]  = 0;
             edgelist_for_face[9]  = 8;
             edgelist_for_face[10] = 3;
             edgelist_for_face[11] = 9;
 
             edgelist_for_face[12] = 1;
             edgelist_for_face[13] = 10;
             edgelist_for_face[14] = 2;
             edgelist_for_face[15] = 11;
 
             // z-faces
             edgelist_for_face[16] = 0;
             edgelist_for_face[17] = 5;
             edgelist_for_face[18] = 1;
             edgelist_for_face[19] = 4;
 
             edgelist_for_face[20] = 2;
             edgelist_for_face[21] = 6;
             edgelist_for_face[22] = 3;
             edgelist_for_face[23] = 7;
 
             sph_poly3d(ifinquiry, dim, ctrp, r0, cp8[0], nnode, nedge, nface,
                        nedge_for_face, nodelist_for_edge, edgelist_for_face,
                        vcell/factor, mixed, vol);
 
             *vol *= factor;
 
             assert(*vol >= 0.0);
             assert(*vol <= vcell + small);
         }
     }
     else {
         memcpy(ctrp, ctr, (size_t)szdim);
 
         r0 = r[0];
         for (i = 1; i < dim; i++) {
             if (r[i] < r0) {
                 r0 = r[i];
             }
         }
         factor = 1.0;
         for (i = 0; i < dim; i++) {
             tmp      = r[i]/r0;
             factor  *= tmp;
             tmp      = 1.0/tmp;
             myxl[i]  = xl[i]  * tmp;
             mydx[i]  = dx[i]  * tmp;
             myctr[i] = ctr[i] * tmp;
         }
         gsph_rec(ifinquiry, mixed, geop, dim, myctr, r0, myxl, mydx, vol);
 
         if (!ifinquiry) {
             if (*mixed) {
                 *vol *= factor;
             }
             else {
                 *vol = vcell;
             }
             assert(*vol >= 0.0);
             assert(*vol <= vcell + small);
         }
         else if (*vol != 0.0) {
             *vol = vcell;
         }
     }
 
     return;
 }
 
void gconic_rec_rotate(int ifinquiry, int *mixed,
                int dim, double r0, double r1, double *ctr0, double *ctr1,
                double *xl, double *dx, double *vol)
{
     int    i, i1, k, nnode, nedge, nface;
     int    intersected, anyinside, alloutside, allinside, ifinside[8];
     int    nedge_for_face[6], nodelist_for_edge[24], edgelist_for_face[24];
     double vcell, sinp, cosp, sint, cost, drdx, tmp;
     double x0, x1, r, r2, ds2, ddx;
     double xpl[3], xpr[3], c8[8][3], cp8[8][3], ctrp[3], ctr1p[3];
     double *c, *cp;
 
     vcell = 1.0;
     for (i = 0; i < dim; i++) {
         xpl[i]  = xl[i] - ctr0[i];
         xpr[i]  = xpl[i] + dx[i];
         ctrp[i] = ctr1[i] - ctr0[i];
         vcell  *= dx[i];
     }
     sinp = sinphi_conic;
     cosp = cosphi_conic;
     sint = sintheta_conic;
     cost = costheta_conic;
 
     c8[0][0] = xpl[0];
     c8[0][1] = xpl[1];
     c8[0][2] = xpl[2];
 
     c8[1][0] = xpr[0];
     c8[1][1] = xpl[1];
     c8[1][2] = xpl[2];
 
     c8[2][0] = xpr[0];
     c8[2][1] = xpr[1];
     c8[2][2] = xpl[2];
 
     c8[3][0] = xpl[0];
     c8[3][1] = xpr[1];
     c8[3][2] = xpl[2];
 
     c8[4][0] = xpl[0];
     c8[4][1] = xpl[1];
     c8[4][2] = xpr[2];
 
     c8[5][0] = xpr[0];
     c8[5][1] = xpl[1];
     c8[5][2] = xpr[2];
 
     c8[6][0] = xpr[0];
     c8[6][1] = xpr[1];
     c8[6][2] = xpr[2];
 
     c8[7][0] = xpl[0];
     c8[7][1] = xpr[1];
     c8[7][2] = xpr[2];
 
     nnode = 8;
     for (i = 0; i < nnode; i++) {
         c  = c8[i];
         cp = cp8[i];
         cp[0] =  cost * cosp * c[0] + cost * sinp * c[1] - sint * c[2];
         cp[1] = -       sinp * c[0] +        cosp * c[1];
         cp[2] =  sint * cosp * c[0] + sint * sinp * c[1] + cost * c[2];
     }
     ctr1p[0] =  cost * cosp * ctrp[0] + cost * sinp * ctrp[1] - sint * ctrp[2];
     ctr1p[1] = -       sinp * ctrp[0] +        cosp * ctrp[1];
     ctr1p[2] =  sint * cosp * ctrp[0] + sint * sinp * ctrp[1] + cost * ctrp[2];
 
     tmp = (fabs(ctr1p[1]) + fabs(ctr1p[2]))/dx[0];
     assert(tmp < small);
 
     ctr1p[1] = 0.0;
     ctr1p[2] = 0.0;
 
     x0 = 0.0;         //  ctr0p[0];
     x1 = ctr1p[0];
 
     drdx = (r1 - r0)/(x1 - x0);
 
//   check whether all vertices are inside the sphere.
 
     for (i = 0; i < nnode; i++) {
         cp  = cp8[i];
         ddx = cp[0] - x0;
         r   = r0 + drdx * ddx;
         r2  = r * r;
         ds2 = 0.0;
         for (k = 1; k < dim; k++) {
             ds2 += (cp[k] * cp[k]);
         }
         if (ds2 < r2) {
             if ((x0 < cp[0]) && (cp[0] < x1)) {
                 ifinside[i] = 1;  // inside
             }
             else if ((cp[0] < x0) || (cp[0] > x1)) {
                 ifinside[i] = 0;  // outside
             }
             else {
                 ifinside[i] = -1; // touch
             }
         }
         else if (ds2 > r2) {
             ifinside[i] = 0;  // outside
         }
         else { // (ds2 == r2)
             if ((x0 <= cp[0]) && (cp[0] <= x1)) {
                 ifinside[i] = -1; // touch
             }
             else {
                 ifinside[i] = 0;  // outside
             }
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < nnode; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < nnode; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < nnode; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *vol   = vcell;
         *mixed = 0;
         return;
     }
     else if (alloutside) {
 
//       conservetivatively guarantee outside
 
/***
         intersected = 0;
         for (i = 0; i < nnode; i++) {
             cp  = cp8[i];
             ddx = cp[0] - x0;
             r   = r0 + drdx * ddx;
 
             for (k = 1; k < dim; k++) {
                 if ((-r < cp[k]) && (cp[k] < r)) {
                     intersected = 1;       // possiblly
                     break;
                 }
             }
             if (intersected) {
                 break;
             }
         }
*****/
         *vol   = 0.0;
         *mixed = 0;
         return;
     }
     else if (anyinside) {
         if (ifinquiry) {
            *mixed = 1;
            *vol = small;
            return;
         }
     }
     else if (ifinquiry) {
         *mixed = 0;
         *vol = 0.0;
         return;
     }
     nface = 6;
     nedge = 12;
     nnode = 8;
     for (i = 0; i < nface; i++) {
         nedge_for_face[i] = 4;
     }
     // x-edge
     nodelist_for_edge[0] = 0;
     nodelist_for_edge[1] = 1;
 
     nodelist_for_edge[2] = 2;
     nodelist_for_edge[3] = 3;
 
     nodelist_for_edge[4] = 6;
     nodelist_for_edge[5] = 7;
 
     nodelist_for_edge[6] = 4;
     nodelist_for_edge[7] = 5;
     // y-edges
     nodelist_for_edge[8]  = 0;
     nodelist_for_edge[9]  = 3;
 
     nodelist_for_edge[10] = 1;
     nodelist_for_edge[11] = 2;
 
     nodelist_for_edge[12] = 5;
     nodelist_for_edge[13] = 6;
 
     nodelist_for_edge[14] = 4;
     nodelist_for_edge[15] = 7;
     // z-edges
     nodelist_for_edge[16] = 0;
     nodelist_for_edge[17] = 4;
 
     nodelist_for_edge[18] = 1;
     nodelist_for_edge[19] = 5;
 
     nodelist_for_edge[20] = 2;
     nodelist_for_edge[21] = 6;
 
     nodelist_for_edge[22] = 3;
     nodelist_for_edge[23] = 7;
 
     // x-faces
     edgelist_for_face[0] = 4;
     edgelist_for_face[1] = 11;
     edgelist_for_face[2] = 7;
     edgelist_for_face[3] = 8;
 
     edgelist_for_face[4] = 5;
     edgelist_for_face[5] = 9;
     edgelist_for_face[6] = 6;
     edgelist_for_face[7] = 10;
 
     // y-faces
     edgelist_for_face[8]  = 0;
     edgelist_for_face[9]  = 8;
     edgelist_for_face[10] = 3;
     edgelist_for_face[11] = 9;
 
     edgelist_for_face[12] = 1;
     edgelist_for_face[13] = 10;
     edgelist_for_face[14] = 2;
     edgelist_for_face[15] = 11;
 
     // z-faces
     edgelist_for_face[16] = 0;
     edgelist_for_face[17] = 5;
     edgelist_for_face[18] = 1;
     edgelist_for_face[19] = 4;
 
     edgelist_for_face[20] = 2;
     edgelist_for_face[21] = 6;
     edgelist_for_face[22] = 3;
     edgelist_for_face[23] = 7;
 
     conic_poly3d(ifinquiry, dim, x0, r0, x1, r1, cp8[0], nnode, nedge, nface,
                  nedge_for_face, nodelist_for_edge, edgelist_for_face,
                  vcell, mixed, vol);
     *mixed = 1;
     assert(*vol >= 0.0);
     assert(*vol <= vcell + small);
 
     return;
 }
 
void gconic_rec_cyl(int ifinquiry, int *mixed, int ifcyl,
                int dim, double r0, double r1, double *cl, double *cr,
                double *xl, double *dx, double *vol)
{
//   The top and bottom are assumed to be on the y-z plane.
//   If ifinquiry = 1, only determine the cell is mixed or npt through output mixed.
//   If ifinquiry = 0, will calculate vol.
//   If mixed = 0, and vol != 0.0, the cell is fully within the structure.
//   If mixed = 0, and vol == 0.0, the cell is outside of the structure.
 
     int    ifinside[8], intersected, allinside, anyinside, alloutside;
     int    i, k;
 
     double vcell, dxinv, drdx, dcydx, dczdx, ddx, r, r2, ds2;
     double xr[3], ctr[2], dc[2], *p, c8[8][3];
 
     vcell  = 1.0;
     for (i = 0; i < dim; i++) {
         xr[i]  = xl[i] + dx[i];
         vcell *= dx[i];
     }
     c8[0][0] = xl[0];
     c8[0][1] = xl[1];
     c8[0][2] = xl[2];
 
     c8[1][0] = xr[0];
     c8[1][1] = xl[1];
     c8[1][2] = xl[2];
 
     c8[2][0] = xr[0];
     c8[2][1] = xr[1];
     c8[2][2] = xl[2];
 
     c8[3][0] = xl[0];
     c8[3][1] = xr[1];
     c8[3][2] = xl[2];
 
     c8[4][0] = xl[0];
     c8[4][1] = xl[1];
     c8[4][2] = xr[2];
 
     c8[5][0] = xr[0];
     c8[5][1] = xl[1];
     c8[5][2] = xr[2];
 
     c8[6][0] = xr[0];
     c8[6][1] = xr[1];
     c8[6][2] = xr[2];
 
     c8[7][0] = xl[0];
     c8[7][1] = xr[1];
     c8[7][2] = xr[2];
 
     dxinv = 1.0/(cr[0] - cl[0]);
     drdx  = (r1 - r0) * dxinv;
     dcydx = (cr[1] - cl[1]) * dxinv;
     dczdx = (cr[2] - cl[2]) * dxinv;
 
     for (i = 0; i < 8; i++) {
         p   = c8[i];
         ddx = p[0] - cl[0];
         r   = r0 + drdx * ddx;
         r2  = r * r;
         ctr[0] = cl[1] + dcydx * ddx;
         ctr[1] = cl[2] + dczdx * ddx;
         dc[0]  = p[1] - ctr[0];
         dc[1]  = p[2] - ctr[1];
         ds2    = 0.0;
         for (k = 0; k < 2; k++) {
             ds2 += (dc[k] * dc[k]);
         }
         if (ds2 < r2) {
             if ((cl[0] < p[0]) && (p[0] < cr[0])) {
                 ifinside[i] = 1;  // inside
             }
             else if ((p[0] < cl[0]) || (p[0] > cr[0])) {
                 ifinside[i] = 0;  // outside
             }
             else {
                 ifinside[i] = -1; // touch
             }
         }
         else if (ds2 > r2) {
             ifinside[i] = 0;  // outside
         }
         else { // (ds2 == r2)
             if ((cl[0] <= p[0]) && (p[0] <= cr[0])) {
                 ifinside[i] = -1; // touch
             }
             else {
                 ifinside[i] = 0;  // outside
             }
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < 8; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < 8; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < 8; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *vol   = vcell;
         *mixed = 0;
         return;
     }
     else if (alloutside) {
 
/******
//       conservetivatively guarantee outside
 
         intersected = 0;
         for (i = 0; i < 8; i++) {
             p   = c8[i];
             ddx = p[0] - cl[0];
             r   = r0 + drdx * ddx;
             r2  = r * r;
             ctr[0] = cl[1] + dcydx * ddx;
             ctr[1] = cl[2] + dczdx * ddx;
             dc[0]  = p[1] - ctr[0];
             dc[1]  = p[2] - ctr[1];
             ds2    = 0.0;
 
             for (k = 0; k < 2; k++) {
                 if ( (-r < dc[k]) && (dc[k] < r)) {
                     intersected = 1;
                     break;
                 }
             }
             if (intersected) {
                 break;
             }
         }
***********/
         *vol   = 0.0;
         *mixed = 0;
         return;
     }
     else if (anyinside) {
         if (ifinquiry) {
            *mixed = 1;
            *vol = small;
            return;
         }
     }
     if (ifcyl > 0) {
         conic_rec0(dim, r0, cl, cr, xl, dx, vol);
     }
     else if (ifcyl < 0) {
         conic_rec(dim, r0, r1, cl, cr, xl, dx, vol);
     }
     else {
          printf("ERROR: gconic_rec_cyl\n");
     }
     *mixed = 1;
 
     return;
 }
 
void conic_rec(int dim, double r0, double r1, double *cl, double *cr,
             double *xl, double *dx, double *vol)
{
//   The top and bottom are assumed to be on the y-z plane.
 
     int    ifinquiry, mixed;
 
     int    i, i2, npart, npart2, nit;
     double sixth, vcell, xlow, xhgh, dxinv, drdx, dcydx, dczdx;
     double r, ddx, ddx0, dx6, x, area, dvol, vol_old, err;
     double *aend, *am, *tmp;
     double xr[3], c0[2], dc[2], ctr[2];
 
     sixth    = 0.1666666666666666666667;
     ifinquiry = 0;
 
     vcell = 1.0;
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
         vcell *= dx[i];
     }
     xlow  = MAX(cl[0], xl[0]);
     xhgh  = MIN(cr[0], xr[0]);
 
     if (xlow >= xhgh) {
         *vol = 0.0;
         return;
     }
     dxinv = 1.0/(cr[0] - cl[0]);
     drdx  = (r1 - r0) * dxinv;
     dcydx = (cr[1] - cl[1]) * dxinv;
     dczdx = (cr[2] - cl[2]) * dxinv;
 
     npart = 4;
     aend  = (double *) malloc((npart + npart + 1) * sizeof(double));
     am    = aend + (npart + 1);
     c0[0] = xl[1];
     c0[1] = xl[2];
     dc[0] = dx[1];
     dc[1] = dx[2];
 
     ddx0   = (xhgh - xlow)/(double)npart;
     x  = xlow - ddx0;
     for (i = 0; i <= npart; i++) {
         x += ddx0;
         ddx = x - cl[0];
         r      = r0    + drdx  * ddx;
         ctr[0] = cl[1] + dcydx * ddx;
         ctr[1] = cl[2] + dczdx * ddx;
 
         gsph_rec(0, &mixed, 1, 2, ctr, r, c0, dc, &area);
         aend[i] = area;
     }
     vol_old = 1.0e+30;
     err = 1.0;
     nit = 0;
     while (err > tol) {
         ddx0 = (xhgh - xlow)/(double)npart;
         dx6  = sixth * ddx0;
         *vol = 0.0;
         x    = xlow - 0.5 * ddx0;
         for (i = 0; i < npart; i++) {
             x += ddx0;
             ddx    = x     - cl[0];
             r      = r0    + drdx  * ddx;
             ctr[0] = cl[1] + dcydx * ddx;
             ctr[1] = cl[2] + dczdx * ddx;
             gsph_rec(0, &mixed, 1, 2, ctr, r, c0, dc, &area);
             am[i] = area;
             dvol  = dx6 *(aend[i] + 4.0 * am[i] + aend[i+1]);
             *vol += dvol;
         }
         err = fabs(*vol - vol_old)/vcell;
 
         if (err > tol) {
             npart2 = npart + npart;
             tmp = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
             assert(tmp);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 tmp[i2]   = aend[i];
                 tmp[i2+1] = am[i];
             }
             tmp[npart2] = aend[npart];
             free(aend);
             aend = tmp;
             am   = aend + (npart2 + 1);
 
             npart = npart2;
             nit++;
         }
         vol_old = *vol;
     }
     free(aend);
     assert(*vol >= 0.0);
     assert(*vol <= vcell + small);
 
     return;
 }
 
 
void conic_rec0(int dim, double r, double *cl, double *cr,
             double *xl, double *dx, double *vol)
{
//   The top and bottom of the cylinder are assumed to be on the y-z plane.
 
     int    mixed;
     int    i;
     double xlow, xhgh, area;
 
     assert((cl[1] == cr[1]) && (cl[2] == cr[2]));
 
     xlow  = MAX(cl[0], xl[0]);
     xhgh  = MIN(cr[0], xl[0] + dx[0]);
 
     if (xlow >= xhgh) {
         *vol = 0.0;
         return;
     }
     gsph_rec(0, &mixed, 1, 2, cl + 1, r, xl + 1, dx + 1, &area);
     *vol = area * (xhgh - xlow);
 
     return;
 }
 
void rec_rec(int ifinquiry, int *ifmixed, int geop,
             double vcell, int dim, double *xxl, double *xxr,
             double *xl, double *dx, double *vol)
{
     int i, k, nc, outside, inside;
     double z0, z1, r0, r1, h, dr, dvol, small;
     double xr[3];
     double myxl[2][2], myxr[2][2], myxxl[2][2], myxxr[2][2];
     double *cl, *cr, *ccl, *ccr;
 
     small = 1.0e-10 * dx[0];
 
     for (k = 0; k < dim; k++) {
         xr[k] = xl[k] + dx[k];
     }
     outside = 0;
     for (k = 0; k < dim; k++) {
         if ((xxl[k] + small > xr[k]) || (xxr[k] < xl[k] + small)) {
             outside = 1;
             break;
         }
     }
     if (outside) {
         *vol     = 0.0;
         *ifmixed = 0;
         return;
     }
     inside = 1;
     for (k = 0; k < dim; k++) {
         if ( !((xl[k] + small > xxl[k]) && (xr[k] < xxr[k] + small)) ) {
            inside = 0;
            break;
         }
     }
     if (inside) {
         *vol = vcell;
         *ifmixed = 0;
         return;
     }
     if (ifinquiry == 1) {
         *vol = small;
         *ifmixed = 1;
         return;
     }
     *ifmixed = 1;
     if ((dim == 1) || (geop == 1)) {
         *vol = 1.0;
         for (k = 0; k < dim; k++) {
             z0 = MAX(xl[k], xxl[k]);
             z1 = MIN(xr[k], xxr[k]);
             *vol *= (z1 - z0);
         }
     }
     else if ((dim == 2) && (geop == 2)) { 
         z0  = MAX(xl[1], xxl[1]);
         z1  = MIN(xr[1], xxr[1]); 
         h   = MAX(0.0, z1 - z0);
         r0  = MAX(xl[0], xxl[0]); 
         r1  = MIN(xr[0], xxr[0]);
         dr  = MAX(0.0, r1 - r0);
         *vol = h * (r1 + r0) * dr;  
     } 
     else if ((dim == 1) && (geop == 2)) {  // ????  
         *vol = 1.0;
         for (k = 0; k < dim; k++) {
             z0 = MAX(xl[k], xxl[k]);
             z1 = MIN(xr[k], xxr[k]);
             *vol *= (z1*z1 - z0*z0);
         }
     }
     else if ((dim == 1) || (geop == 3)) {
         *vol = (4.0/3.0);
         for (k = 0; k < dim; k++) {
             z0 = MAX(xl[k], xxl[k]);
             z1 = MIN(xr[k], xxr[k]);
             *vol *= (z1*z1*z1 - z0*z0*z0);
         }
     }
     else if (dim == 2) {
         if ((xl[0] < 0.0) && (xr[0] > 0.0)) {
             nc = 2;
             myxl[0][0] = 0.0;
             myxr[0][0] = -xl[0];
 
             myxl[1][0] = 0.0;
             myxr[1][0] =  xr[0];
 
             myxxl[0][0] = 0.0;
             myxxr[0][0] = MIN(0.0, -xxl[0]);
             myxxl[1][0] = MAX(0.0, xxl[0]);
             myxxr[1][0] = xxr[0];
         }
         else if (xl[0] < 0.0) {
             nc = 1;
             myxl[0][0]  = -xr[0];
             myxr[0][0]  = -xl[0];
             myxxl[0][0] = -xxr[0];
             myxxr[0][0] = -xxl[0];
         }
         else {
             nc = 1;
             myxl[0][0]  = xl[0];
             myxr[0][0]  = xr[0];
             myxxl[0][0] = xxl[0];
             myxxr[0][0] = xxr[0];
         }
         for (k = 0; k < nc; k++) {
             myxl[k][1]  = xl[1];
             myxr[k][1]  = xr[1];
             myxxl[k][1] = xxl[1];
             myxxr[k][1] = xxr[1];
         }
         *vol = 0.0;
         for (k = 0; k < nc; k++) {
             cl = myxl[k];
             cr = myxr[k];
             ccl = myxxl[k];
             ccr = myxxr[k];
             r0  = MAX(cl[0], ccl[0]);
             r1  = MIN(cr[0], ccr[0]);
             z0  = MAX(cl[1], ccl[1]);
             z1  = MIN(cr[1], ccr[1]);
             dvol = (r1 + r0)*(r1 - r0)*(z1 - z0);  // * PI
             *vol += dvol;
         }
     }
     return;
 }
 
void gpoly3d_rec_oso(int nreg, int ifinquiry, int *ifmixed, double vcell,
                         int ifcyl, double *xl, double *dx, double *vol)
{
//   ifcyl  : input
//            1 : along x
//            2 : along y
//            3 : along z
 
     int  inquiry, mixed;
     int  dim, ifile, i, k, g, nn, option;
     int  nnode, nedge, nface, ngrp;
     int  *nnode_grp, *nedge_grp, *nface_grp;
     int  *nodelist_for_edge, *edgelist_for_face, *nodelist_for_face;
     int  *nnode_for_face;
     int  **nodelist_for_edge_grp, **edgelist_for_face_grp;
     int  **nodelist_for_face_grp, **nnode_for_face_grp;
     double dvol;
 
     double *xyz, *factor_grp, *normal_of_face, *ctr, *ctr_grp;
     double **normal_of_face_grp, **coords_grp;
 
//   Since I do not have the code working for the intersection between
//   a concave 3D polyhedron and a plane yet, I have to rely on my old
//   approach until I add the capability.
 
     dim = 3;
 
     ifile  = -1;
     option = 0;
     for (i = 0; i < nfile_oso; i++) {
         if (reg_oso[i] == nreg) {
             ifile = i;
             break;
         }
     }
     assert(ifile >= 0);
     ngrp       = npart_oso[ifile];
     factor_grp = factor_part_oso[ifile];
 
     nnode_grp = nnode_grp_oso[ifile];
     nedge_grp = nedge_grp_oso[ifile];
     nface_grp = nface_grp_oso[ifile];
 
     nodelist_for_edge_grp = nodelist_for_edge_oso[ifile];
     edgelist_for_face_grp = edgelist_for_face_oso[ifile];
     nodelist_for_face_grp = nodelist_for_face_oso[ifile];
     nnode_for_face_grp    = nnode_for_face_oso[ifile];
     normal_of_face_grp    = normal_of_face_oso[ifile];
     coords_grp            = coords_grp_oso[ifile];
     ctr_grp               = ctr_grp_oso[ifile];
 
     if (ngrp == 1) {
             nnode = nnode_grp[0];
             nedge = nedge_grp[0];
             nface = nface_grp[0];
             nodelist_for_edge = nodelist_for_edge_grp[0];
             edgelist_for_face = edgelist_for_face_grp[0];
             nodelist_for_face = nodelist_for_face_grp[0];
             nnode_for_face    = nnode_for_face_grp[0];
             normal_of_face    = normal_of_face_grp[0];
             xyz               = coords_grp[0];
             ctr               = ctr_grp;
 
             gpoly3d_rec(nreg, ifinquiry, ifmixed, ifcyl, dim, xyz, nnode,
                         nedge, nface, nodelist_for_edge, edgelist_for_face,
                         nnode_for_face, nodelist_for_face, nnode_for_face,
                         ctr, normal_of_face,
                         xl, dx, vol);
     }
     else if (ngrp > 1)  {
         inquiry = ifinquiry;;
//       inquiry = 0;
 
         *ifmixed = 0;
         *vol = 0.0;
         for (g = 0; g < ngrp; g++) {
 
             nnode = nnode_grp[g];
             nedge = nedge_grp[g];
             nface = nface_grp[g];
             nodelist_for_edge = nodelist_for_edge_grp[g];
             edgelist_for_face = edgelist_for_face_grp[g];
             nodelist_for_face = nodelist_for_face_grp[g];
             nnode_for_face    = nnode_for_face_grp[g];
             normal_of_face    = normal_of_face_grp[g];
             xyz               = coords_grp[g];
             ctr               = ctr_grp + (dim * g);
 
             gpoly3d_rec(nreg, inquiry, &mixed, ifcyl, dim, xyz, nnode,
                         nedge, nface, nodelist_for_edge, edgelist_for_face,
                         nnode_for_face, nodelist_for_face, nnode_for_face,
                         ctr, normal_of_face,
                         xl, dx, &dvol);
 
             *vol += (factor_grp[g] * dvol);
             *ifmixed += mixed;
         }
     }
 
     return;
 }
 
void gpoly_rec(int mype,
               int nreg, int ifinquiry, int *ifmixed, int geop, double vcell,
               int ifcyl, int dim, double *coords, int nnode,
               double *x2min, double *x2max,
               int nedge, int nface, int *nodelist_for_edge,
               int *edgelist_for_face, int *nedge_ea_face,
               int *nodelist_for_face, int *num_nodes_for_face,
               double *ctr, double *norm_of_face,
               double *xl, double *dx, double *vol)
{
     int  inquiry, mixed;
     int  ifile, i, k, g, ngrp, nn;
     double dvol, x2mn[2], x2mx[2];
 
     int    *nng;
     double *factor_grp, *p;
     double *xys0, **xysg, *xys;
 
     ifile  = -1;
     ngrp   = -1;
     for (i = 0; i < nfile_oso; i++) {
         if (reg_oso[i] == nreg) {
             ifile = i;
             ngrp  = npart_oso[ifile];
             break;
         }
     }
     if ((ifile < 0) || (ngrp <= 0) || (if_decompose_concave == 0)) {
         if (dim == 2) {
             poly2d_rec(ifinquiry, ifmixed, geop, vcell, dim, nnode, coords,
                        x2min, x2max, xl, dx, vol);
         }
         else if (dim == 3) {
             gpoly3d_rec(nreg, ifinquiry, ifmixed, ifcyl, dim, coords, nnode,
                         nedge, nface, nodelist_for_edge, edgelist_for_face,
                         nedge_ea_face, nodelist_for_face, num_nodes_for_face,
                         ctr, norm_of_face,
                         xl, dx, vol);
         }
     }
     else {
         *ifmixed = 0;
 
         ngrp       = npart_oso[ifile];
         xysg       = coords_part_oso[ifile];
         nng        = nnode_part_oso[ifile];
         factor_grp = factor_part_oso[ifile];
 
         if (ngrp > 1) {
             inquiry = 0;
         }
         else {
             inquiry = ifinquiry;
         }
         *vol = 0.0;
         for (g = 0; g < ngrp; g++) {
             xys = xysg[g];
             nn  = nng[g];
 
             if (dim == 2) {
                 x2mn[0] = xys[0];
                 x2mn[1] = xys[1];
                 memcpy(x2mx, x2mn, (size_t)(2 * sizeof(double)));
                 for (i = 1; i < nn; i++) {
                     p = xys + (i + i);
                     for (k = 0; k < 2; k++) {
                         if (p[k] < x2mn[k]) {
                             x2mn[k] = p[k];
                         }
                         else if (p[k] > x2mx[k]) {
                             x2mx[k] = p[k];
                         }
                     }
                 }
                 poly2d_rec(inquiry, &mixed, geop, vcell, dim, nn, xys,
                            x2mn, x2mx, xl, dx, &dvol);
 
                 *vol += (factor_grp[g] * dvol);
             }
             else if (dim == 3) {
                 gpoly3d_rec(nreg, inquiry, &mixed, ifcyl, dim, xys, nn,
                             nedge, nface, nodelist_for_edge, edgelist_for_face,
                             nedge_ea_face, nodelist_for_face, num_nodes_for_face,
                             ctr, norm_of_face,
                             xl, dx, &dvol);
             }
             *vol += (factor_grp[g] * dvol);
             *ifmixed += mixed;
         }
     }
     return;
 }
 
void oso_rec(int mype,
             int nreg, int ifinquiry, int *ifmixed, int geop, double vcell,
             int dim, double *xl, double *dx, double *vol)
{
//   option: 1 or 2
 
     char name[128];
     int  inquiry, mixed, found, szdim, fileid;
     int  ifile, nnode, i, k, g, ngrp;
     double factor, dvol;
     double x2min[2], x2max[2];
     double *p, *coords, *factor_grp;
 
     szdim  = 2 * sizeof(double);
 
     ifile  = -1;
     for (i = 0; i < nfile_oso; i++) {
         if (reg_oso[i] == nreg) {
             ifile = i;
             break;
         }
     }
     assert(ifile >= 0);
 
     ngrp = npart_oso[ifile];
 
     if (ngrp == 1) {
         nnode = nnode_part_oso[ifile][0];
         coords = coords_part_oso[ifile][0];
         x2min[0] = coords[0];
         x2min[1] = coords[1];
         memcpy(x2max, x2min, (size_t)szdim);
 
         for (i = 1; i < nnode; i++) {
             p = coords + (i + i);
             for (k = 0; k < 2; k++) {
                 if (p[k] < x2min[k]) {
                     x2min[k] = p[k];
                 }
                 else if (p[k] > x2max[k]) {
                     x2max[k] = p[k];
                 }
             }
         }
         poly2d_rec(ifinquiry, ifmixed, geop, vcell, dim, nnode, coords,
                    x2min, x2max, xl, dx, vol);
     }
     else  {
         inquiry = 0;
         *vol = 0.0;
 
         for (g = 0; g < ngrp; g++) {
             nnode = nnode_part_oso[ifile][g];
             coords = coords_part_oso[ifile][g];
             factor = factor_part_oso[ifile][g];
 
             x2min[0] = coords[0];
             x2min[1] = coords[1];
             memcpy(x2max, x2min, (size_t)szdim);
             for (i = 1; i < nnode; i++) {
                 p = coords + (i + i);
                 for (k = 0; k < 2; k++) {
                     if (p[k] < x2min[k]) {
                         x2min[k] = p[k];
                     }
                     else if (p[k] > x2max[k]) {
                         x2max[k] = p[k];
                     }
                 }
             }
             poly2d_rec(inquiry, &mixed, geop, vcell, dim, nnode, coords,
                        x2min, x2max, xl, dx, &dvol);
 
             *vol     += (factor * dvol);
             *ifmixed += mixed;
         }
     }
     if ((*vol < vcell) && (*vol > 0.0)) {
         *ifmixed = 1;
     }
     else {
         *ifmixed = 0;
     }
     return;
 }
 
 
 
 
 
 
 
 
 
 
void poly2d_rec(int ifinquiry, int *mixed, int geop, double vcell, int dim,
                int nnode, double *xys, double *x2min, double *x2max,
                double *xl, double *dx, double *vol)
{
     point_position  pos, pos_old;
     const  int mxnn = 128;
     int    err, allocated, anyinside, alloutside, allinside, intersected;
     int    i, i1, k, nxing, nint, szdim, touch_only;
     int    myifinside[mxnn], *ifinside;
//     double x2min[2], x2max[2];
     double xr[2], rint[4];
     double myxyint[mxnn + mxnn], c5[5][2];
     double *xyint;
     double *p0, *p1;
 
     szdim  = dim * sizeof(double);
 
     if (nnode < mxnn) {
         ifinside = myifinside;
//         xyint  = myxyint;
           allocated = 0;
     }
     else {
         ifinside = (int *) malloc((nnode + 1) * sizeof(int));
//         xyint  = (double *) malloc((nnode + 1) * szdim);
           allocated = 1;
     }
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     c5[0][0] = xl[0];
     c5[0][1] = xl[1];
 
     c5[1][0] = xr[0];
     c5[1][1] = xl[1];
 
     c5[2][0] = xr[0];
     c5[2][1] = xr[1];
 
     c5[3][0] = xl[0];
     c5[3][1] = xr[1];
     memcpy(c5[4], c5[0], (size_t)szdim);
 
/****
     for (k = 0; k < dim; k++) {
         x2min[k] = xys[k];
         x2max[k] = x2min[k];
     }
     for (i = 1; i < nnode; i++) {
         p0 = xys + (i + i);
         for (k = 0; k < dim; k++) {
             if (x2min[k] > p0[k]) {
                 x2min[k] = p0[k];
             }
             else if (x2max[k] < p0[k]) {
                 x2max[k] = p0[k];
             }
         }
     }
****/
     for (k = 0; k < dim; k++) {
         if ((x2min[k] >= xr[k]) || (x2max[k] < xl[k])) {
             *vol = 0.0;
             *mixed = 0;
             if (allocated) {
                 free(ifinside);
//               free(xyint);
             }
             return;
         }
     }
//   check each vertex of polygon is inside the polygon or not
 
     for (i = 0; i < 4; i++) {
         p0 = c5[i];
         pos = point_in_polygon(xys, nnode, x2min, x2max, p0, &touch_only);
//       pos_old = point_in_polygon_old(xys, nnode, p0, &touch_only);
//       assert(pos == pos_old);
 
         if (pos == inside) {
             ifinside[i] = 1;
         }
         else if (pos == outside) {
             ifinside[i] = 0;
         }
         else if (pos == touching) {
             ifinside[i] =  -1;
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < 4; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < 4; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < 4; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *vol = vcell;
         *mixed = 0;
         if (allocated) {
             free(ifinside);
//             free(xyint);
         }
         return;
     }
     else if (anyinside) {
         *mixed = 1;
         *vol   = small;
         if (ifinquiry) {
             if (allocated) {
                 free(ifinside);
//                 free(xyint);
             }
             return;
         }
     }
     else if (alloutside) {
 
//       check whether there is any intersection of edges of two polygons.
 
         intersected = 0;
         for (i = 0; i < nnode; i++) {
             p0 = xys + (i + i);
             for (k = 0; k < 4; k++) {
                 p1 = c5[k];
                 nxing = line_intersect2(p0, p1, rint, &touch_only);
                 if (nxing) {
                     intersected = 1;
                     break;
                 }
             }
             if (intersected) {
                 break;
             }
         }
         if (intersected) {
             *mixed = 1;
             if (ifinquiry) {
                 *vol = small;
                 if (allocated) {
                     free(ifinside);
//                     free(xyint);
                 }
                 return;
             }
         }
         else {
             *mixed = 0;
             *vol = 0.0;
             if (allocated) {
                 free(ifinside);
//                 free(xyint);
             }
             return;
         }
     }
     for (i = 0; i < nnode; i++) {
         p0 = xys + (i + i);
         if ( (p0[0] > xl[0]) && (p0[0] < xr[0]) &&
              (p0[1] > xl[1]) && (p0[1] < xr[1])) {
              ifinside[i] = 1;   // inside
         }
         else if ((p0[0] < xl[0]) || (p0[0] > xr[0]) ||
                  (p0[1] < xl[1]) || (p0[1] > xr[1])) {
              ifinside[i] = 0; // outside
         }
         else if ((p0[0] == xl[0]) || (p0[0] == xr[0]) ||
                  (p0[1] == xl[1]) || (p0[1] == xr[1]) ) {
 
              ifinside[i] = -1; // touch
         }
         else {
              printf("ERROR: ifinside in poly2d_rec\n");
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < nnode; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < nnode; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < nnode; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         if (geop == 1) {
             for (i = 0; i < nnode; i++) {
                 i1 = (i + 1) % nnode;
                 p0 = xys + (i + i);
                 p1 = xys + (i1 + i1);
                 *vol += (p0[0] * p1[1] - p0[1] * p1[0]);
             }
             *vol *= 0.5;
             if (*vol < 0.0) *vol = - (*vol);
         }
         else if (geop == 2) {
             rz_area(nnode, xys, vol);
         }
         if (*vol >= vcell) {
             *mixed = 0;
         }
         else {
             *mixed = 1;
         }
     }
     else if (anyinside) {
         *mixed = 1;
         *vol = small;
         if (ifinquiry) {
             if (allocated) {
                 free(ifinside);
//                 free(xyint);
             }
             return;
         }
     }
     xyint = NULL;
     err = hull_test(xys, nnode, c5[0], 4, &xyint, &nint);
 
//     err = myihull(xys, nnode, x2min, x2max, c5[0], 4, dx[0], xyint, &nint);
 
//     if (err == -1) {
//         ihull2(xys, nnode, c5[0], 4, dx[0], xyint, &nint, NULL);
////         write_hull("test1", 2, NULL, 0, 1, xys, nnode, c5[0], 4, xyint, nint);
//     }
     *vol = 0.0;
     if (nint >= 3) {
         if (geop == 1) {
             for (i = 0; i < nint; i++) {
                 i1 = (i + 1) % nint;
                 p0 = xyint + (i + i);
                 p1 = xyint + (i1 + i1);
                 *vol += (p0[0] * p1[1] - p0[1] * p1[0]);
             }
             *vol *= 0.5;
             if (*vol < 0.0) *vol = - (*vol);
         }
         else if (geop == 2) {
             rz_area(nint, xyint, vol);
         }
     }
     if (xyint) free(xyint);
 
     assert(*vol >= 0.0);
 
/***
     if ((*vol > vcell + small) && !err) {
         write_hull("test1", 2, NULL, 0, 1, xys, nnode, c5[0], 4, xyint, nint);
         ihull2(xys, nnode, c5[0], 4, dx[0], xyint, &nint, NULL);
         write_hull("test2", 2, NULL, 0, 1, xys, nnode, c5[0], 4, xyint, nint);
         *vol = 0.0;
         if (nint >= 3) {
             if (geop == 1) {
                 for (i = 0; i < nint; i++) {
                     i1 = (i + 1) % nint;
                     p0 = xyint + (i + i);
                     p1 = xyint + (i1 + i1);
                     *vol += (p0[0] * p1[1] - p0[1] * p1[0]);
                 }
                 *vol *= 0.5;
                 if (*vol < 0.0) *vol = - (*vol);
             }
             else if (geop == 2) {
                 rz_area(nint, xyint, vol);
             }
         }
     }
******/
     assert(*vol <= vcell + small);
 
     *vol = MIN(*vol, vcell);
     if ((*vol == 0.0) || (*vol == vcell)) {
         *mixed = 0;
     }
     else {
         *mixed = 1;
     }
 
     if (allocated) {
         free(ifinside);
//         free(xyint);
     }
     return;
 }
 
void gpoly3d_rec_old(int nreg, int ifinquiry, int *mixed, int ifcyl, int dim,
                 double *coords, int nnode,
                 int nedge, int nface, int *nodelist_for_edge,
                 int *edgelist_for_face, int *nedge_ea_face,
                 int *nodelist_for_face, int *num_nodes_for_face,
                 double *ctr, double *norm_of_face,
                 double *xl, double *dx, double *vol)
{
     int    nn, i, i2, i3, k;
     double z0, z1;
     double myxl[3], mydx[3];
     double xy1[66], xy2[66];
     double *xyz1, *xyz2, *pt1, *pt2, *ps1, *ps2;
 
     assert(nnode < 32);
 
     nn = nnode / 2;
     assert(nn + nn == nnode);
     xyz1 = coords;
     xyz2 = coords + (nn * dim);
 
 
/*********** For DEBUG
double t = 0.25 * PI;
double p = 0.25 * PI;
double sinp, cosp, sint, cost;
double pxyz1[100], pxyz2[100], *c, *cp;
 
sint = sin(t);
cost = cos(t);
sinp = sin(p);
cosp = cos(p);
 
for (i = 0; i < nn; i++) {
    c  = xyz1 + dim * i;
    cp = pxyz1 + dim * i;
    cp[0] =  cost * cosp * c[0] + cost * sinp * c[1] - sint * c[2];
    cp[1] = -       sinp * c[0] +        cosp * c[1];
    cp[2] =  sint * cosp * c[0] + sint * sinp * c[1] + cost * c[2];
}
for (i = 0; i < nn; i++) {
    c  = xyz2 + dim * i;
    cp = pxyz2 + dim * i;
    cp[0] =  cost * cosp * c[0] + cost * sinp * c[1] - sint * c[2];
    cp[1] = -       sinp * c[0] +        cosp * c[1];
    cp[2] =  sint * cosp * c[0] + sint * sinp * c[1] + cost * c[2];
}
**********/
 
     if ((ifcyl == 3) || (ifcyl == -3)) {
         z0 = xyz1[2];
         z1 = xyz2[2];
         for (i = 0; i < nn; i++) {
             i2 = i + i;
             i3 = i2 + i;
             pt1 = xy1 + i2;
             pt2 = xy2 + i2;
             ps1 = xyz1 + i3;
             ps2 = xyz2 + i3;
             for (k = 0; k < 2; k++) {
                 pt1[k] = ps1[k];
                 pt2[k] = ps2[k];
             }
          }
          for (k = 0; k < dim; k++) {
              myxl[k] = xl[k];
              mydx[k] = dx[k];
          }
      }
      else if ((ifcyl == 2) || (ifcyl == -2)) {
          z0 = xyz1[1];
          z1 = xyz2[1];
          for (i = 0; i < nn; i++) {
             i2 = i + i;
             i3 = i2 + i;
             pt1 = xy1 + i2;
             pt2 = xy2 + i2;
             ps1 = xyz1 + i3;
             ps2 = xyz2 + i3;
 
             pt1[0] = ps1[2];
             pt1[1] = ps1[0];
             pt2[0] = ps2[2];
             pt2[1] = ps2[0];
          }
          myxl[0] = xl[2];
          myxl[1] = xl[0];
          myxl[2] = xl[1];
          mydx[0] = dx[2];
          mydx[1] = dx[0];
          mydx[2] = dx[1];
      }
      else if ((ifcyl == 1) || (ifcyl == -1)) {
          z0 = xyz1[0];
          z1 = xyz2[0];
          for (i = 0; i < nn; i++) {
             i2 = i + i;
             i3 = i2 + i;
             pt1 = xy1 + i2;
             pt2 = xy2 + i2;
             ps1 = xyz1 + i3;
             ps2 = xyz2 + i3;
 
             pt1[0] = ps1[1];
             pt1[1] = ps1[2];
             pt2[0] = ps2[1];
             pt2[1] = ps2[2];
          }
          myxl[0] = xl[1];
          myxl[1] = xl[2];
          myxl[2] = xl[0];
          mydx[0] = dx[1];
          mydx[1] = dx[2];
          mydx[2] = dx[0];
      }
      if (ifcyl == 0) {
          xy1[nn + nn]     = xy1[0];
          xy1[nn + nn + 1] = xy1[1];
          xy2[nn + nn]     = xy2[0];
          xy2[nn + nn + 1] = xy2[1];
 
          poly3dg_rec(ifinquiry, mixed, ifcyl, dim, coords, nnode,
                      nedge, nface,
                      nodelist_for_edge, edgelist_for_face, nedge_ea_face,
                      nodelist_for_face, num_nodes_for_face,
                      ctr, norm_of_face, xl, dx, vol);
      }
      else {
          memcpy(xy1 + (nn + nn), xy1, (size_t)(2 * sizeof(double)));
          poly3d_rec(ifinquiry, mixed, ifcyl, dim,
                     nn, xy1, xy2, z0, z1, myxl, mydx, vol);
      }
      return;
}
 
void poly3dg_rec(int ifinquiry, int *mixed, int ifcyl, int dim,
                 double *coords, int nnode,
                 int nedge, int nface, int *nodelist_for_edge,
                 int *edgelist_for_face, int *nedge_ea_face,
                 int *nodelist_for_face, int *num_nodes_for_face,
                 double *ctr, double *norm_of_face,
                 double *xl, double *dx, double *vol)
{
//   This is for an arbitrary polyhedrons to intersect a cube.
 
     int    szdim, szdim2;
     int    i, i0, i1, i2, k, k1, f, n0, ns, np, s, offset;
     int    anyinside, alloutside, allinside, found, done;
     int    nit, npart, npart2, nint, nxyint;
     int    *edgelist, *nodelist;
     int    ifinside[8];
     double sixth, dxinv, vcell, x0, x1, xlow, xhgh, tmp;
     double dctr, dcube, dot_ctr, dot_cube;
     double area, areamx, dxmin, x, ddx, dx6, dvol, vol_old, err;
     double *norm, *c, *p0, *p1, *p2, *fend, *fm, *array;
     double xr[3], x3min[3], x3max[3], c8[8][3], c5[5][2];
     double xs[32], vols[32], tmp1d[64];
     double pts_int[64], *xyint; // xyint[64];
     int    ngrp, g;
     double dxgrp, xlow_grp, xhgh_grp, vol_grp;
 
//////// For debug
 
     char name_cube[] = "cube";
     char name_poly[] = "poly";
     char name_ctr[]  = "ctr";
     int fid;
     int mynodelist[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
     double myx[32], myy[32], myz[32];
 
//     mio_Unstructured_Mesh mesh;
//     mio_Coord *coord;
//////////////////////
 
 
     assert(nnode < 32);
     assert(nface < 32);
 
     sixth = 0.1666666666666667;
     vcell = 1.0;
     for (i = 0; i < dim; i++) {
         vcell *= dx[i];
         xr[i] = xl[i] + dx[i];
     }
//   calculate xmin, xmax, etc for the polyhedron.
 
     szdim2 = sizeof(double);
     szdim  = dim * szdim2;
     szdim2 += szdim2;
 
     memcpy(x3min, coords, (size_t)szdim);
     memcpy(x3max, coords, (size_t)szdim);
     for (i = 1; i < nnode; i++) {
         p0 = coords + (i * dim);
         for (k = 0; k < dim; k++) {
             if (p0[k] < x3min[k]) {
                 x3min[k] = p0[k];
             }
             else if (p0[k] > x3max[k]) {
                 x3max[k] = p0[k];
             }
         }
     }
     alloutside = 0;
     for (k = 0; k < dim; k++) {
         if ((xl[k] >= x3max[k]) || (xr[k] <= x3min[k])) {
             alloutside = 1;
             break;
         }
     }
     if (alloutside) {
         *mixed = 0;
         *vol   = 0.0;
         return;
     }
     c8[0][0] = xl[0];
     c8[0][1] = xl[1];
     c8[0][2] = xl[2];
 
     c8[1][0] = xr[0];
     c8[1][1] = xl[1];
     c8[1][2] = xl[2];
 
     c8[2][0] = xr[0];
     c8[2][1] = xr[1];
     c8[2][2] = xl[2];
 
     c8[3][0] = xl[0];
     c8[3][1] = xr[1];
     c8[3][2] = xl[2];
 
     c8[4][0] = xl[0];
     c8[4][1] = xl[1];
     c8[4][2] = xr[2];
 
     c8[5][0] = xr[0];
     c8[5][1] = xl[1];
     c8[5][2] = xr[2];
 
     c8[6][0] = xr[0];
     c8[6][1] = xr[1];
     c8[6][2] = xr[2];
 
     c8[7][0] = xl[0];
     c8[7][1] = xr[1];
     c8[7][2] = xr[2];
 
//   check whether all vertices are inside the sphere.
 
     for (i = 0; i < 8; i++) {
         c  = c8[i];
         ifinside[i] = 1;   // inside
 
         offset = 0;
         for (f = 0; f < nface; f++) {
             norm = norm_of_face + (f * dim);
             nodelist = nodelist_for_face + offset;
             n0 = nodelist[0];
             p0 = coords + (n0 * dim);
             dot_ctr  = 0.0;
             dot_cube = 0.0;
             for (k = 0; k < dim; k++) {
                 dctr  = ctr[k] - p0[k];
                 dcube = c[k]   - p0[k];
                 dot_ctr  += (norm[k] * dctr);
                 dot_cube += (norm[k] * dcube);
             }
             if (dot_ctr * dot_cube < 0.0) {
                 ifinside[i] = 0;  // outside
                 break;
             }
             else if (dot_ctr * dot_cube == 0.0) {
                 ifinside[i] = -1; // touch
             }
             offset += num_nodes_for_face[f];
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < 8; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < 8; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < 8; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *vol   = vcell;
         *mixed = 0;
         return;
     }
     else if (alloutside) {
         *vol   = 0.0;
         *mixed = 0;
         return;
     }
     else if (anyinside) {
         if (ifinquiry) {
            *mixed = 1;
            *vol = small;
            return;
         }
     }
     else if (ifinquiry) {
         *mixed = 0;
         *vol = 0.0;
         return;
     }
 
 
 
//   FOR DEBUG
 
/***
 
     mio_create_file("test", mio_file_create, mio_independent, &fid);
 
//      mio_init(mio_umesh, -1, &mesh);
//      coord = &(mesh.coord);
//      coord->datatype = mio_double;
//      coord->coord[0] = (void *)ctr;
//      coord->coord[1] = (void *)(ctr + 1);
//      coord->coord[2] = (void *)(ctr + 2);
 
//      mesh.name     = name_ctr;
//      mesh.dims     = dim;
//      mesh.sizes[3] = 1;
//      mesh.type     = mio_point;
 
//      mio_write(mio_umesh, fid, &mesh);
 
      mio_init(mio_umesh, -1, &mesh);
      coord = &(mesh.coord);
      coord->datatype = mio_double;
 
      mesh.name     = name_ctr;
      mesh.idmin    = 0;
      mesh.dims     = dim;
      mesh.datatype = mio_int;
      mesh.sizes[3] = 8;
      mesh.sizes[0] = 1;
      mesh.type     = mio_hex;
      mesh.nodelist_for_zone = mynodelist;
      for (i = 0; i < 8; i++) {
          myx[i] = ctr[0] + 0.1 * c8[i][0];
          myy[i] = ctr[1] + 0.1 * c8[i][1];
          myz[i] = ctr[2] + 0.1 * c8[i][2];
      }
      coord->coord[0] = (void *)myx;
      coord->coord[1] = (void *)myy;
      coord->coord[2] = (void *)myz;
      mio_write(mio_umesh, fid, &mesh);
 
      mio_init(mio_umesh, -1, &mesh);
      coord = &(mesh.coord);
      coord->datatype = mio_double;
 
      mesh.name = name_cube;
      mesh.idmin    = 0;
      mesh.dims     = dim;
      mesh.datatype = mio_int;
      mesh.sizes[3] = 8;
      mesh.sizes[0] = 1;
      mesh.type = mio_hex;
      mesh.nodelist_for_zone = mynodelist;
      for (i = 0; i < 8; i++) {
          myx[i] = c8[i][0];
          myy[i] = c8[i][1];
          myz[i] = c8[i][2];
      }
 
      coord->coord[0] = (void *)myx;
      coord->coord[1] = (void *)myy;
      coord->coord[2] = (void *)myz;
      mio_write(mio_umesh, fid, &mesh);
 
      mio_init(mio_umesh, -1, &mesh);
      coord = &(mesh.coord);
      coord->datatype = mio_double;
 
      mesh.name     = name_poly;
      mesh.idmin    = 0;
      mesh.dims     = dim;
      mesh.datatype = mio_int;
      mesh.sizes[3] = nnode;
      mesh.sizes[1] = nface;
      mesh.sizes[0] = 1;
      mesh.type = mio_general_mesh;
      mesh.nodelist_for_face = nodelist_for_face;
      mesh.facelist_for_zone = mynodelist;
      mesh.num_nodes_for_face = num_nodes_for_face;
      mesh.num_faces_for_zone = &nface;
 
      for (i = 0; i < nnode; i++) {
          p0   = coords + (i * dim);
          myx[i] = p0[0];
          myy[i] = p0[1];
          myz[i] = p0[2];
      }
      coord->coord[0] = (void *)myx;
      coord->coord[1] = (void *)myy;
      coord->coord[2] = (void *)myz;
      mio_write(mio_umesh, fid, &mesh);
 
      mio_close_file(fid);
 
///////////////////////
****/
 
     *mixed = 1;
     c5[0][0] = c8[0][1];
     c5[0][1] = c8[0][2];
 
     c5[1][0] = c8[3][1];
     c5[1][1] = c8[3][2];
 
     c5[2][0] = c8[7][1];
     c5[2][1] = c8[7][2];
 
     c5[3][0] = c8[4][1];
     c5[3][1] = c8[4][2];
 
     memcpy(c5[4], c5[0], (size_t)szdim2);
 
//   determine the segments
 
     xlow_grp = xl[0];
     xhgh_grp = xr[0];
 
     dxgrp = (xhgh_grp - xlow_grp)/(double)ngrp;
     xlow  = xlow_grp - dxgrp;
 
     for (g = 0; g < ngrp; g++) {
         xlow += dxgrp;
         xhgh  = xlow + dxgrp;
 
         vol_grp = 0.0;
         npart = 8;
 
    //   Simpson integral
 
         fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
         assert(fend);
 
         fm    = fend + (npart + 1);
         ddx   = (xhgh - xlow)/(double)npart;
         x     = xlow - ddx;
         for (i = 0; i <= npart; i++) {
             x += ddx;
 
    //       find the polygon resulted from the intersection between
    //       the plane x = x and the polyhedron.
 
             plane_poly3d(dim, x, coords, nnode, nedge, nface, nedge_ea_face,
                          nodelist_for_edge, edgelist_for_face,
                          &nint, pts_int);
             area = 0.0;
             if (nint >= 3) {
                 memcpy(pts_int + (nint + nint), pts_int, (size_t)szdim2);
 
                 xyint = NULL;
                 err = hull_test(pts_int, nint, c5[0], 4, &xyint, &nxyint);
 
//                 err = myihull(pts_int, nint, NULL, NULL, c5[0], 4, dx[0], xyint, &nxyint);
//                 if (err == -1) {
//                     ihull2(pts_int, nint, c5[0], 4, dx[0], xyint, &nxyint, tmp1d);
//                 }
 
                 if (nxyint >= 3) {
                     for (k = 0; k < nxyint; k++) {
                         k1 = (k + 1) % nxyint;
                         p0 = xyint + (k  + k);
                         p1 = xyint + (k1 + k1);
                         area += (p0[0] * p1[1] - p0[1] * p1[0]);
                     }
                     area *= 0.5;
                     if (area < 0.0) area = -area;
                 }
                 if (xyint) free(xyint);
             }
             fend[i] = area;
             areamx = MAX(areamx, area);
         }
         vol_old = 1.0e+30;
         done = 0;
         nit = 0;
         while (!done) {
              ddx = (xhgh - xlow)/(double)npart;
              dx6 = sixth * ddx;
              vol_grp = 0.0;
              x = xlow - 0.5 * ddx;
              for (i = 0; i < npart; i++) {
                  x += ddx;
 
                  plane_poly3d(dim, x, coords, nnode, nedge, nface, nedge_ea_face,
                               nodelist_for_edge, edgelist_for_face,
                               &nint, pts_int);
                  area = 0.0;
                  if (nint >= 3) {
                      memcpy(pts_int + (nint + nint), pts_int, (size_t)szdim2);
 
                      xyint = NULL;
                      err = hull_test(pts_int, nint, c5[0], 4, &xyint, &nxyint);
 
//                      err = myihull(pts_int, nint, NULL, NULL, c5[0], 4, dx[0], xyint, &nxyint);
//                      if (err == -1) {
//                          ihull2(pts_int, nint, c5[0], 4, dx[0], xyint, &nxyint, tmp1d);
//                      }
 
                      if (nxyint >= 3) {
                          for (k = 0; k < nxyint; k++) {
                              k1 = (k + 1) % nxyint;
                              p0 = xyint + (k  + k);
                              p1 = xyint + (k1 + k1);
                              area += (p0[0] * p1[1] - p0[1] * p1[0]);
                          }
                          area *= 0.5;
                          if (area < 0.0) area = -area;
                      }
                      if (xyint) free(xyint);
                  }
                  fm[i] = area;
                  dvol  = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
                  vol_grp += dvol;
                  areamx = MAX(areamx, area);
              }
              err   = fabs(vol_grp - vol_old)/vcell;
              dxmin = tol * (vcell/areamx);
              dxmin = 0.1 * sqrt(dxmin);
 
    //        if ((err > tol) && (ddx > dxmin)) {
              if (err > tol) {
                  npart2 = npart + npart;
                  array  = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
                  assert(array);
 
                  for (i = 0; i < npart; i++) {
                      i2 = i + i;
                      array[i2]   = fend[i];
                      array[i2+1] = fm[i];
                  }
                  array[npart2] = fend[npart];
                  free(fend);
                  fend = array;
                  fm   = fend + (npart2 + 1);
 
                  npart = npart2;
                  nit++;
              }
              else {
                  done = 1;
              }
              vol_old = vol_grp;
         }
         free(fend);
 
         *vol += vol_grp;
 
     }
 
     return;
 }
 
void poly3d_rec(int ifinquiry, int *mixed, int ifcyl, int dim,
                 int nnode, double *xys1,
                 double *xys2,
                 double z0, double z1, double *xl, double *dx, double *vol)
{
//   The top and bottom are assumed to be on the y-z plane.
//   If ifinquiry = 1, only determine the cell is mixed or npt through output mixed.
//   If ifinquiry = 0, will calculate vol.
//   If mixed = 0, and vol != 0.0, the cell is fully within the structure.
//   If mixed = 0, and vol == 0.0, the cell is outside of the structure.
 
     point_position pos, pos_old;
 
     int    ifinside[100], allinside, anyinside, alloutside;
     int    i, j, k, nn, offset, k_touch;
     double ptr[ 2 * 100], dpsdz[100][2];
     double vcell, dzinv, drdx, dcydx, dczdx, ddz, zlow, zhgh;
     double *p1, *p2, *p, *dp, *ps;
     double xr[3], ctr[2], dc[2], c8[8][3];
 
     assert(nnode < 100);
 
     vcell  = 1.0;
     for (i = 0; i < dim; i++) {
         xr[i]  = xl[i] + dx[i];
         vcell *= dx[i];
     }
     if ((xl[2] >= z1) || (xr[2] <= z0)) {
         *mixed = 0;
         *vol   = 0.0;
         return;
     }
     c8[0][0] = xl[0];
     c8[0][1] = xl[1];
     c8[0][2] = xl[2];
 
     c8[1][0] = xr[0];
     c8[1][1] = xl[1];
     c8[1][2] = xl[2];
 
     c8[2][0] = xr[0];
     c8[2][1] = xr[1];
     c8[2][2] = xl[2];
 
     c8[3][0] = xl[0];
     c8[3][1] = xr[1];
     c8[3][2] = xl[2];
 
     if (ifcyl > 0) {
 
         nn = 4;
         for (k = 0; k < nn; k++) {
             pos = point_in_polygon(xys1, nnode, NULL, NULL, c8[k], &k_touch);
//           pos_old = point_in_polygon_old(xys1, nnode, c8[k], &k_touch);
//           assert(pos_old == pos);
 
             if (pos == inside) {
                 ifinside[k] = 1;  // inside
             }
             else if (pos == outside) {
                 ifinside[k] = 0;   // outside
             }
             else {
                 ifinside[k] = -1;  // touch
             }
         }
     }
     else {
         nn = 8;
 
         c8[4][0] = xl[0];
         c8[4][1] = xl[1];
         c8[4][2] = xr[2];
 
         c8[5][0] = xr[0];
         c8[5][1] = xl[1];
         c8[5][2] = xr[2];
 
         c8[6][0] = xr[0];
         c8[6][1] = xr[1];
         c8[6][2] = xr[2];
 
         c8[7][0] = xl[0];
         c8[7][1] = xr[1];
         c8[7][2] = xr[2];
 
         dzinv = 1.0/(z1 - z0);
 
         offset = 0;
         for (i = 0; i < nnode; i++) {
             p1 = xys1 + offset;
             p2 = xys2 + offset;
             dp = dpsdz[i];
             for (k = 0; k < 2; k++) {
                 dp[k] = (p2[k] - p1[k]) * dzinv;
             }
             offset += 2;
         }
 
         for (k = 0; k < nn; k++) {
             p   = c8[k];
             ddz = p[2] - z0;
 
             offset = 0;
             for (j = 0; j < nnode; j++) {
                 p1 = xys1 + offset;
                 ps = ptr  + offset;
                 dp = dpsdz[j];
                 for (i = 0; i < 2; i++) {
                     ps[i] = p1[i] + dp[i] * ddz;
                 }
                 offset += 2;
             }
             pos = point_in_polygon(ptr, nnode, NULL, NULL, p, &k_touch);
//           pos_old = point_in_polygon_old(ptr, nnode, p, &k_touch);
//           assert(pos == pos_old);
 
             if (pos == inside) {
                 ifinside[k] = 1;  // inside
             }
             else if (pos == outside) {
                 ifinside[k] = 0;   // outside
             }
             else {
                 ifinside[k] = -1;  // touch
             }
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < nn; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < nn; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < nn; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         zlow = MAX(z0, xl[2]);
         zhgh = MIN(z1, xr[2]);
         *vol = (zhgh - zlow) * dx[0] * dx[1];
         *mixed = 1;
         return;
     }
     else if (alloutside) {
         *vol   = 0.0;
         *mixed = 0;
         return;
     }
     else if (anyinside) {
         if (ifinquiry) {
            *mixed = 1;
            *vol   = small;
            return;
         }
         else {
            if (ifcyl > 0) {
                poly_rec3d0(dim, nnode, xys1, NULL, NULL, z0, z1, xl, dx, vol);
            }
            else {
                poly_rec3d(dim, nnode, xys1, xys2, z0, z1, xl, dx, vol);
            }
            *mixed = 1;
         }
     }
     else {
         *vol   = 0.0;
         *mixed = 0;
     }
     assert(*vol >= 0.0);
     assert(*vol <= vcell + small);
 
     return;
 }
 
void poly_rec3d(int dim, int nnode, double *xys1,
                double *xys2, double z0, double z1,
                double *xl, double *dx, double *vol)
{
 
//   This is very slow, need to be investigated.
 
     int    szdim2, done, ifinquiry, mixed;
     int    i, i1, i2, k, offset, nint, npart, npart2, nit;
     double sixth, dzmin, vcell, dzinv, zlow, zhgh, areamx;
     double ddz, ddz0, dz6, z, dstot, area, dvol, vol_old, err;
     double *aend, *am, *tmp1d, *p0, *p1, *p2, *dp;
     double xr[3], c0[2], dc[2];
 
     double rp[10];
     double dpdz[200];
//     double xyint[200], xyin[100];
     double *xyint;
     double xym[200], xytmp[200];
     double *poly;
     double pr[10];
 
     szdim2 = 2 * sizeof(double);
     sixth     = 0.1666666666666666666667;
     ifinquiry = 0;
 
     vcell = 1.0;
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
         vcell *= dx[i];
     }
     pr[0] = xl[0];
     pr[1] = xl[1];
 
     pr[2] = xr[0];
     pr[3] = xl[1];
 
     pr[4] = xr[0];
     pr[5] = xr[1];
 
     pr[6] = xl[0];
     pr[7] = xr[1];
     memcpy(pr + 8, pr, (size_t)szdim2);
 
     dzinv = 1.0/(z1 - z0);
     offset = 0;
     for (i = 0; i < nnode; i++) {
         p1 = xys1 + offset;
         p2 = xys2 + offset;
         dp = dpdz + offset;
         for (k = 0; k < 2; k++) {
             dp[k] = (p2[k] - p1[k]) * dzinv;
         }
         offset += 2;
     }
     npart = 4;
     aend  = (double *) malloc((npart + npart + 1) * sizeof(double));
     am    = aend + (npart + 1);
     c0[0] = xl[0];
     c0[1] = xl[1];
     dc[0] = dx[0];
     dc[1] = dx[1];
 
     zlow = MAX(z0, xl[2]);
     zhgh = MIN(z1, xr[2]);
     ddz0  = (zhgh - zlow)/(double)npart;
     z  = zlow - ddz0;
     for (i = 0; i <= npart; i++) {
         z += ddz0;
         ddz    = z - z0;
         offset = 0;
         for (i1 = 0; i1 < nnode; i1++) {
             p1 = xys1 + offset;
             p2 = xys2 + offset;
             p0 = xym  + offset;
             dp = dpdz + offset;
 
             for (k = 0; k < 2; k++) {
                 p0[k] = p1[k] + dp[k] * ddz;
             }
             offset += 2;
         }
/***
         dstot = 0.0;
         for (i1 = 0; i1 < nnode; i1++) {
             i2 = (i1 + 1) % nnode;
             p1 = xym + (i1 + i1);
             p2 = xym + (i2 + i2);
             for (k = 0; k < 2; k++) {
                 dstot += fabs(p2[k] - p1[k]);
             }
         }
***/
         memcpy(xym + (nnode + nnode), xym, (size_t)szdim2);
         area= 0.0;
 
         xyint = NULL;
         err = hull_test(xym, nnode, pr, 4, &xyint, &nint);
 
//         err = myihull(xym, nnode, NULL, NULL, pr, 4, dx[0], xyint, &nint);
//         if (err == -1) {
//             ihull2(xym, nnode, pr, 4, dx[0], xyint, &nint, xytmp);
//         }
 
         if (nint >= 3) {
             for (i1 = 0; i1 < nint; i1++) {
                 i2 = (i1 + 1) % nint;
                 p1 = xyint + (i1 + i1);
                 p2 = xyint + (i2 + i2);
                 area += (p1[0] * p2[1] - p1[1] * p2[0]);
             }
             area *= 0.5;
             if (area < 0.0) area = - area;
         }
         aend[i] = area;
 
         if (xyint) free(xyint);
     }
     done = 0;
     dzmin = 0.1 * tol * dx[2];
     ddz0 = dx[2];
     vol_old = 1.0e+30;
     err = 1.0;
     nit = 0;
     while (!done)  {
         ddz0 = (zhgh - zlow)/(double)npart;
         dz6  = sixth * ddz0;
         *vol = 0.0;
         z    = zlow - 0.5 * ddz0;
         areamx = small;
         for (i = 0; i < npart; i++) {
             z += ddz0;
             ddz    = z - z0;
             offset = 0;
             for (i1 = 0; i1 < nnode; i1++) {
                 p1 = xys1 + offset;
                 p2 = xys2 + offset;
                 p0 = xym  + offset;
                 dp = dpdz + offset;
                 for (k = 0; k < 2; k++) {
                     p0[k] = p1[k] + dp[k] * ddz;
                 }
                 offset += 2;
             }
             memcpy(xym + (nnode + nnode), xym, (size_t)szdim2);
 
             xyint = NULL;
             err = hull_test(xym, nnode, pr, 4, &xyint, &nint);
 
//             err = myihull(xym, nnode, NULL, NULL, pr, 4, dx[0], xyint, &nint);
//             if (err == -1) {
//                 ihull2(xym, nnode, pr, 4, dx[0], xyint, &nint, xytmp);
//             }
             area = 0.0;
             if (nint >= 3) {
                 for (i1 = 0; i1 < nint; i1++) {
                     i2 = (i1 + 1) % nint;
                     p1 = xyint + (i1 + i1);
                     p2 = xyint + (i2 + i2);
                     area += (p1[0] * p2[1] - p1[1] * p2[0]);
                 }
                 area *= 0.5;
                 if (area < 0.0) area = - area;
                 if (area > areamx) areamx = area;
             }
             am[i] = area;
             dvol  = dz6 *(aend[i] + 4.0 * am[i] + aend[i+1]);
             *vol += dvol;
 
             if (xyint) free(xyint);
         }
         err = fabs(*vol - vol_old)/vcell;
         dzmin = 0.1 * sqrt(tol * vcell/areamx);
 
//       if ((err > tol) && (ddz0 > dzmin)) {
         if (err > tol) {
             npart2 = npart + npart;
             tmp1d = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
             assert(tmp1d);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 tmp1d[i2]   = aend[i];
                 tmp1d[i2+1] = am[i];
             }
             tmp1d[npart2] = aend[npart];
             free(aend);
             aend = tmp1d;
             am   = aend + (npart2 + 1);
 
             npart = npart2;
             nit++;
         }
         else {
             done = 1;
         }
         vol_old = *vol;
     }
     free(aend);
     assert(*vol >= 0.0);
     assert(*vol <= vcell + small);
 
     return;
 }
 
 
void poly_rec3d0(int dim, int nnode, double *xys, double *x2min, double *x2max,
                 double z0, double z1,
                 double *xl, double *dx, double *vol)
{
//   The top and bottom of the cylinder are assumed to be on the y-z plane.
 
     int    err, i, i1, nint;
     double area, zlow, zhgh, dz;
     double xr[3], c5[5][2], tmp1d[64], *xyint; // xyint[64];
     double *p0, *p1;
 
     assert(nnode < 32);
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     c5[0][0] = xl[0];
     c5[0][1] = xl[1];
 
     c5[1][0] = xr[0];
     c5[1][1] = xl[1];
 
     c5[2][0] = xr[0];
     c5[2][1] = xr[1];
 
     c5[3][0] = xl[0];
     c5[3][1] = xr[1];
     c5[4][0] = c5[0][0];
     c5[4][1] = c5[0][1];
 
     xyint = NULL;
     err = hull_test(xys, nnode, c5[0], 4, &xyint, &nint);
 
//     err = myihull(xys, nnode, x2min, x2max, c5[0], 4, dx[0], xyint, &nint);
//     if (err == -1) {
//         ihull2(xys, nnode, c5[0], 4, dx[0], xyint, &nint, tmp1d);
//     }
     *vol = 0.0;
     area = 0.0;
     if (nint >= 3) {
         for (i = 0; i < nint; i++) {
             i1 = (i + 1) % nint;
             p0 = xyint + (i + i);
             p1 = xyint+ (i1 + i1);
             area += (p0[0] * p1[1] - p0[1] * p1[0]);
         }
         area *= 0.5;
         if (area < 0.0) area = - area;
 
         zlow = MAX(z0, xl[2]);
         zhgh = MIN(z1, xr[2]);
         dz   = MAX(0.0, zhgh - zlow);
         *vol = dz * area;
     }
     if (xyint) free(xyint);
 
     return;
 }
 
void psph_rec(int geop, int dim, double *ctr, double radius,
              double ang_start, double ang_end,
              double *xl, double *dx, double *vol)
{
     int szdim, i, i1, k, n, n1, nn, skip, mixed, intsected[2];
     double small, ang, tmp, t, t2[2];
     double xr[2], c4[4][2], r2[4], *p0, *p1;
     double ang2[2], tana[2], c[10][2];
     double ctr0[3], ctroid[3];
 
     small = 1.0e-10 * dx[0];
 
     if (dim != 2) {
         printf("ERROR: dim !=2 in psph_rec\n");
         return;
     }
     szdim = dim * sizeof(double);
 
     ang2[0] = ang_start;
     ang2[1] = ang_end;
     t = PI / 180.0;
     for (k = 0; k < 2; k++) {
         if (ang2[k] == 0.0) {
             tana[k] = 0.0;
         }
         else if (ang2[k] == 90.0) {
             tana[k] = 1.0e+10;
         }
         else {
             tana[k] = tan(t * ang2[k]);
         }
     }
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     c4[0][0] = xl[0] - ctr[0];
     c4[0][1] = xl[1] - ctr[1];
     c4[1][0] = xr[0] - ctr[0];
     c4[1][1] = xl[1] - ctr[1];
     c4[2][0] = xr[0] - ctr[0];
     c4[2][1] = xr[1] - ctr[1];
     c4[3][0] = xl[0] - ctr[0];
     c4[3][1] = xr[1] - ctr[1];
 
     n = 0;
     for (i = 0; i < 4; i++) {
         i1 = (i + 1) % 4;
         for (k = 0; k < 2; k++) {
             intsected[k] = 0;
         }
         skip = 0;
         for (k = 0; k < 2; k++) {
             if (skip) continue;
 
             ang = ang2[k];
             if (ang == 0.0) {
                 if ((i == 0) || (i == 2)) {
                     tmp = xl[1] - ctr[1];
                     if (fabs(tmp) < tiny) {
                         intsected[0] = 1;
                         intsected[1] = 1;
                         memcpy(r2,     c4[i],  (size_t)szdim);
                         memcpy(r2 + 2, c4[i1], (size_t)szdim);
                         skip = 1;
                         continue;
                     }
                 }
             }
             else if (ang == 90.0) {
                 if ((i == 1) || (i == 3)) {
                     tmp = xl[0] - ctr[0];
                     if (fabs(tmp) < tiny) {
                         intsected[0] = 1;
                         intsected[1] = 1;
                         memcpy(r2,     c4[i],  (size_t)szdim);
                         memcpy(r2 + 2, c4[i1], (size_t)szdim);
                         skip = 1;
                         continue;
                     }
                 }
             }
//           find the intersection between the ray and segment
 
             t = tana[k] *(c4[i1][0] - c4[i][0]) - (c4[i1][1] - c4[i][1]);
             t = (c4[i][1] - tana[k] * c4[i][0]) / t;
             t2[k] = t;
             if ((t >= 0.0) && (t <= 1.0)) {
                 intsected[k]  = 1;
                 r2[k + k]     = c4[i][0] + t * (c4[i1][0] - c4[i][0]);
                 r2[k + k + 1] = c4[i][1] + t * (c4[i1][1] - c4[i][1]);
             }
        }
        if (intsected[0] && intsected[1]) {
            if ((i == 0) || (i == 3)) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2 + 2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2 + 2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2 + 2, (size_t)szdim);
                    n++;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2, (size_t)szdim);
                    n++;
                }
            }
            else if ((i == 1) || (i == 2)) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2, (size_t)(szdim + szdim));
                        n += 2;
                    }
                }
                else {
                    memcpy(c[n], r2, (size_t)(szdim + szdim));
                    n += 2;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2 + 2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2 + 2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2 + 2, (size_t)szdim);
                    n++;
                }
            }
        }
        else if (intsected[0]) {
            if (i == 0) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = c4[i];
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], c4[i], (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], c4[i], (size_t)szdim);
                    n++;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2, (size_t)szdim);
                    n++;
                }
            }
            else if (i == 1) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2, (size_t)szdim);
                    n++;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = c4[i1];
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], c4[i1], (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], c4[i1], (size_t)szdim);
                    n++;
                }
            }
            else if (i == 2) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2, (size_t)szdim);
                    n++;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = c4[i1];
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], c4[i1], (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], c4[i1], (size_t)szdim);
                    n++;
                }
            }
            else {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2, (size_t)szdim);
                    n++;
                }
            }
        }
        else if (intsected[1]) {
            if (i == 0) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2 + 2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2+2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2+2, (size_t)szdim);
                    n++;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = c4[i1];
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], c4[i1], (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], c4[i1], (size_t)szdim);
                    n++;
                }
            }
            else if (i == 1) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = c4[i];
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], c4[i], (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], c4[i], (size_t)szdim);
                    n++;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2 + 2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2+2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2+2, (size_t)szdim);
                    n++;
                }
            }
            else if (i == 2) {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = c4[i];
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], c4[i], (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], c4[i], (size_t)szdim);
                    n++;
                }
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2 + 2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2+2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2+2, (size_t)szdim);
                    n++;
                }
            }
            else if (i == 3)  {
                if (n > 0) {
                    n1 = n - 1;
                    p0 = c[n1];
                    p1 = r2 + 2;
                    tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                    if (tmp > small) {
                        memcpy(c[n], r2+2, (size_t)szdim);
                        n++;
                    }
                }
                else {
                    memcpy(c[n], r2+2, (size_t)szdim);
                    n++;
                }
            }
        }
        else if (((t2[0] < 0.0)&&(t2[1] > 1.0))||((t2[0] > 1.0)&&(t2[1] < 0.0))) {
            if (n > 0) {
                n1 = n - 1;
                p0 = c[n1];
                p1 = c4[i];
                tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                if (tmp > small) {
                    memcpy(c[n], c4[i], (size_t)szdim);
                    n++;
                }
            }
            else {
                memcpy(c[n], c4[i], (size_t)szdim);
                n++;
            }
            if (n > 0) {
                n1 = n - 1;
                p0 = c[n1];
                p1 = c4[i1];
                tmp = fabs(p0[0] - p1[0]) + fabs(p0[1] - p1[1]);
                if (tmp > small) {
                    memcpy(c[n], c4[i1], (size_t)szdim);
                    n++;
                }
            }
            else {
                memcpy(c[n], c4[i1], (size_t)szdim);
                n++;
            }
        }
    }
    memcpy(c[n], c[0], (size_t)szdim);
    for (k = 0; k < dim; k++) {
        ctr0[k] = 0.0;
    }
    sph_poly2d(0, geop, dim, ctr0, radius, c[0], n, &mixed, vol, ctroid);
 
    return;
 }
 
 
void gsph_rec(int ifinquiry, int *mixed,
              int geop, int dim, double *ctr, double r,
              double *xl, double *dx, double *vol)
{
/**
     allow x-plane x = ctr[0], or y-plane, y = cr0[1], or z plane z = ctr[2],
     cuts the rectangular.
**/
/**
 *           7-----------6
 *          /.         / |
 *         / .        /  |         z
 *        4----------5   |         ^
 *        |  3.......|...2         |     . y
 *        | .        |  /          |   .
 *        |.         | /           | .
 *        0----------1             --------->  x
 *
 ***/
     int    i, i1, k, nr, szdim, skip;
     int    ifinside[8], cut[3];
     int    allinside, anyinside, alloutside, intersected;
 
     double xr[3];
     double tmp, r2, ds2, vcell, dvol;
     double *c0, *dc0, cl[8][3], dc[8][3], ddx[3];
     double *p0, *p1, delx[2], a, b, c, delt2, rt, t0, t1;
 
     szdim = dim * sizeof(double);
 
//   determine whether the x- or y- or z-plane cuts the rectangular
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     if (dim == 2) {
         nr = 4;
 
         cl[0][0] = xl[0];
         cl[0][1] = xl[1];
 
         cl[1][0] = xr[0];
         cl[1][1] = xl[1];
 
         cl[2][0] = xr[0];
         cl[2][1] = xr[1];
 
         cl[3][0] = xl[0];
         cl[3][1] = xr[1];
     }
     else if (dim == 3) {
         nr = 8;
 
         cl[0][0] = xl[0];
         cl[0][1] = xl[1];
         cl[0][2] = xl[2];
 
         cl[1][0] = xr[0];
         cl[1][1] = xl[1];
         cl[1][2] = xl[2];
 
         cl[2][0] = xr[0];
         cl[2][1] = xr[1];
         cl[2][2] = xl[2];
 
         cl[3][0] = xl[0];
         cl[3][1] = xr[1];
         cl[3][2] = xl[2];
 
         cl[4][0] = xl[0];
         cl[4][1] = xl[1];
         cl[4][2] = xr[2];
 
         cl[5][0] = xr[0];
         cl[5][1] = xl[1];
         cl[5][2] = xr[2];
 
         cl[6][0] = xr[0];
         cl[6][1] = xr[1];
         cl[6][2] = xr[2];
 
         cl[7][0] = xl[0];
         cl[7][1] = xr[1];
         cl[7][2] = xr[2];
     }
     r2 = r * r;
     *mixed = 0;
     for (k = 0; k < nr; k++) {
         c0 = cl[k];
         ds2 = 0.0;
         for (i = 0; i < dim; i++) {
             tmp = c0[i] - ctr[i];
             ds2 += (tmp * tmp);
         }
         if (ds2 < r2) {
             ifinside[k] = 1;  // inside
         }
         else if (ds2 > r2) {
             ifinside[k] = 0; // outside
         }
         else {
             ifinside[k] = -1;
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < nr; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < nr; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < nr; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
 
         vcell  = 1.0;
         for (i = 0; i < dim; i++) {
             vcell *= dx[i];
         }
         if ((dim == 2) && (geop == 2)) {
             vcell *= fabs(xl[0] + xr[0]);
         }
         *vol   = vcell;
         *mixed = 0;
         return;
     }
     else if (alloutside) {
 
         intersected = 0;
         if (dim == 2) {
             for (i = 0; i < nr; i++) {
                 i1 = (i + 1) % nr;
                 p0 = cl[i];
                 p1 = cl[i1];
                 for (k = 0; k < dim; k++) {
                     ddx[k]   = p1[k] - p0[k];
                     delx[k] = p0[k] - ctr[k];
                 }
                 a = 0.0;
                 b = 0.0;
                 c = 0.0;
                 for (k = 0; k < dim; k++) {
                     a += (ddx[k] * ddx[k]);
                     b += (ddx[k] * delx[k]);
                     c += (delx[k] * delx[k]);
                 }
                 b  = b + b;
                 c -= r2;
                 delt2 = b * b - 4.0 * a * c;
                 if ((delt2 < 0.0) || (a == 0.0)) continue;
 
                 rt = sqrt(delt2);
                 tmp = 0.5 / a;
                 t0  = (-b - rt) * tmp;
                 t1  = (-b + rt) * tmp;
                 if (((t0 >= 0.0) && (t0 <= 1.0)) || ((t1 >= 0.0) && (t1 <= 1.0))) {
                     intersected = 1;
                     break;
                 }
             }
             if (!intersected) {
                 *vol   = 0.0;
                 *mixed = 0;
                 return;
             }
         }
         else if (dim == 3) {
             check_intersected_sphere(dim, xl, xr, ctr, r, &intersected);
         }
         if (!intersected) {
             *vol   = 0.0;
             *mixed = 0;
             return;
         }
     }
     else if (anyinside) {
         if (ifinquiry) {
            *mixed = 1;
            *vol = small;
            return;
         }
         else {
            *mixed = 1;
         }
     }
     else {
         *mixed = 0;
         *vol = 0.0;
     }
 
     for (i = 0; i < dim; i++) {
         if ((xl[i] < ctr[i]) && (xr[i] > ctr[i])) {
             cut[i] = 1;
         }
         else {
             cut[i] = 0;
         }
     }
     if (dim == 1) {
 
     }
     else if (dim == 2) {
         if (!cut[0] && !cut[1]) {
             memcpy(cl[0], xl, (size_t)szdim);
             memcpy(dc[0], dx, (size_t)szdim);
             nr = 1;
         }
         else if (cut[0] && !cut[1]) {
             ddx[0] = ctr[0] - xl[0];
             memcpy(cl[0], xl, (size_t)szdim);
             dc[0][0] = ddx[0];
             dc[0][1] = dx[1];
             cl[1][0] = ctr[0];
             cl[1][1] = xl[1];
             dc[1][0] = dx[0] - ddx[0];
             dc[1][1] = dx[1];
             nr = 2;
         }
         else if (!cut[0] && cut[1]) {
             ddx[1] = ctr[1] - xl[1];
             memcpy(cl[0], xl, (size_t)szdim);
             dc[0][0] = dx[0];
             dc[0][1] = ddx[1];
             cl[1][0] = xl[0];
             cl[1][1] = ctr[1];
             dc[1][0] = dx[0];
             dc[1][1] = dx[1] - ddx[1];
             nr = 2;
         }
         else {  // (cut[0] && cut[1])
             for (i = 0; i < dim; i++) {
                 ddx[i] = ctr[i] - xl[i];
             }
             cl[0][0] = xl[0];
             cl[0][1] = xl[1];
             dc[0][0] = ddx[0];
             dc[0][1] = ddx[1];
             cl[1][0] = ctr[0];
             cl[1][1] = xl[1];
             dc[1][0] = dx[0] - ddx[0];
             dc[1][1] = ddx[1];
             cl[2][0] = xl[0];
             cl[2][1] = ctr[1];
             dc[2][0] = ddx[0];
             dc[2][1] = dx[1] - ddx[1];
             cl[3][0] = ctr[0];
             cl[3][1] = ctr[1];
             dc[3][0] = dx[0] - ddx[0];
             dc[3][1] = dx[1] - ddx[1];
             nr = 4;
         }
     }
     else if (dim == 3) {
         if (!cut[0] && !cut[1] && !cut[2]) {
             memcpy(cl[0], xl, (size_t)szdim);
             memcpy(dc[0], dx, (size_t)szdim);
             nr = 1;
         }
         else if (cut[0] && !cut[1] && !cut[2]) {
             ddx[0] = ctr[0] - xl[0];
 
             for (i = 0; i < 2; i++) {
                 memcpy(cl[i], xl, (size_t)szdim);
                 memcpy(dc[i], dx, (size_t)szdim);
             }
             dc[0][0] = ddx[0];
             cl[0][0] = ctr[0];
             dc[1][0] = dx[0] - ddx[0];
 
             nr = 2;
         }
         else if (!cut[0] && cut[1] && !cut[2]) {
             ddx[1] = ctr[1] - xl[1];
 
             for (i = 0; i < 2; i++) {
                 memcpy(cl[i], xl, (size_t)szdim);
                 memcpy(dc[i], dx, (size_t)szdim);
             }
             dc[0][1] = ddx[1];
             cl[1][1] = ctr[1];
             dc[1][1] = dx[1] - ddx[1];
 
             nr = 2;
         }
         else if (!cut[0] && !cut[1] && cut[2]) {
             ddx[2] = ctr[2] - xl[2];
 
             for (i = 0; i < 2; i++) {
                 memcpy(cl[i], xl, (size_t)szdim);
                 memcpy(dc[i], dx, (size_t)szdim);
             }
             dc[0][2] = ddx[2];
             cl[1][2] = ctr[2];
             dc[1][2] = dx[2] - ddx[2];
 
             nr = 2;
         }
         else if (cut[0] && cut[1] && !cut[2]) {
             ddx[0] = ctr[0] - xl[0];
             ddx[1] = ctr[1] - xl[1];
 
             memcpy(cl[0], xl, (size_t)szdim);
             dc[0][0] = ddx[0];
             dc[0][1] = ddx[0];
             dc[0][2] = dx[2];
 
             cl[1][0] = ctr[0];
             cl[1][1] = xl[1];
             cl[1][2] = xl[2];
             dc[1][0] = dx[0] - ddx[0];
             dc[1][1] = ddx[1];
             dc[1][2] = dx[2];
 
             cl[2][0] = xl[0];
             cl[2][1] = ctr[1];
             cl[2][2] = xl[2];
             dc[2][0] = ddx[0];
             dc[2][1] = dx[1] - ddx[1];
             dc[2][2] = dx[2];
 
             cl[3][0] = ctr[0];
             cl[3][1] = ctr[1];
             cl[3][2] = xl[2];
             dc[3][0] = dx[0] - ddx[0];
             dc[3][1] = dx[1] - ddx[1];
             dc[3][2] = dx[2];
 
             nr = 4;
         }
         else if (cut[0] && !cut[1] && cut[2]) {
             ddx[0] = ctr[0] - xl[0];
             ddx[2] = ctr[2] - xl[2];
 
             memcpy(cl[0], xl, (size_t)szdim);
             dc[0][0] = ddx[0];
             dc[0][1] = dx[1];
             dc[0][2] = ddx[2];
 
             cl[1][0] = ctr[0];
             cl[1][1] = xl[1];
             cl[1][2] = xl[2];
             dc[1][0] = dx[0] - ddx[0];
             dc[1][1] = dx[1];
             dc[1][2] = ddx[2];
 
             cl[2][0] = xl[0];
             cl[2][1] = xl[1];
             cl[2][2] = ctr[2];
             dc[2][0] = ddx[0];
             dc[2][1] = dx[1];
             dc[2][2] = dx[2] - ddx[2];
 
             cl[3][0] = ctr[0];
             cl[3][1] = xl[1];
             cl[3][2] = ctr[2];
             dc[3][0] = dx[0] - ddx[0];
             dc[3][1] = dx[1];
             dc[3][2] = dx[2] - ddx[2];
 
             nr = 4;
         }
         else if (!cut[0] && cut[1] && cut[2]) {
 
             ddx[1] = ctr[1] - xl[1];
             ddx[2] = ctr[2] - xl[2];
 
             memcpy(cl[0], xl, (size_t)szdim);
             dc[0][0] = dx[0];
             dc[0][1] = ddx[1];
             dc[0][2] = ddx[2];
 
             cl[1][0] = xl[0];
             cl[1][1] = ctr[1];
             cl[1][2] = xl[2];
             dc[1][0] = dx[0];
             dc[1][1] = dx[1] - ddx[1];
             dc[1][2] = ddx[2];
 
             cl[2][0] = xl[0];
             cl[2][1] = xl[1];
             cl[2][2] = ctr[2];
             dc[2][0] = dx[0];
             dc[2][1] = ddx[1];
             dc[2][2] = dx[2] - ddx[2];
 
             cl[3][0] = xl[0];
             cl[3][1] = ctr[1];
             cl[3][2] = ctr[2];
             dc[3][0] = dx[0];
             dc[3][1] = dx[1] - ddx[1];
             dc[3][2] = dx[2] - ddx[2];
 
             nr = 4;
         }
         else {   // all cut
             for (i = 0; i < dim; i++) {
                 ddx[i] = ctr[i] - xl[i];
             }
             memcpy(cl[0], xl, (size_t)szdim);
             for (i = 0; i < dim; i++) {
                 dc[0][i] = ddx[i];
             }
             cl[1][0] = ctr[0];
             cl[1][1] = xl[1];
             cl[1][2] = xl[2];
             dc[1][0] = dx[0] - ddx[0];
             dc[1][1] = ddx[1];
             dc[1][2] = ddx[2];
 
             cl[2][0] = xl[0];
             cl[2][1] = ctr[1];
             cl[2][2] = xl[2];
             dc[2][0] = ddx[0];
             dc[2][1] = dx[1] - ddx[1];
             dc[2][2] = ddx[2];
 
             cl[3][0] = ctr[0];
             cl[3][1] = ctr[1];
             cl[3][2] = xl[2];
             dc[3][0] = dx[0] - ddx[0];
             dc[3][1] = dx[1] - ddx[1];
             dc[3][2] = ddx[2];
 
             cl[4][0] = xl[0];
             cl[4][1] = xl[1];
             cl[4][2] = ctr[2];
             dc[4][0] = ddx[0];
             dc[4][1] = ddx[1];
             dc[4][2] = dx[2] - ddx[2];
 
             cl[5][0] = ctr[0];
             cl[5][1] = xl[1];
             cl[5][2] = ctr[2];
             dc[5][0] = dx[0] - ddx[0];
             dc[5][1] = ddx[1];
             dc[5][2] = dx[2] - ddx[2];
 
             cl[6][0] = xl[0];
             cl[6][1] = ctr[1];
             cl[6][2] = ctr[2];
             dc[6][0] = ddx[0];
             dc[6][1] = dx[1] - ddx[1];
             dc[6][2] = dx[2] - ddx[2];
 
             cl[7][0] = ctr[0];
             cl[7][1] = ctr[1];
             cl[7][2] = ctr[2];
             dc[7][0] = dx[0] - ddx[0];
             dc[7][1] = dx[1] - ddx[1];
             dc[7][2] = dx[2] - ddx[2];
 
             nr = 8;
         }
     }
     vcell = 0.0;
     for (i = 0; i < nr; i++) {
         c0  = cl[i];
         dc0 = dc[i];
         dvol = 1.0;
         for (k = 0; k < dim; k++) {
             dvol *= dc0[k];
         }
         if ((dim == 2) && (geop == 2)) { 
             dvol *= (c0[0] + c0[0] + dc0[0]);
         } 
         vcell += fabs(dvol);
     }
     *vol = 0.0;
     for (i = 0; i < nr; i++) {
         c0  = cl[i];
         dc0 = dc[i];
         dvol = 0.0;
         skip = 0;
         for (k = 0; k < dim; k++) {
             tmp = dc0[k]/dx[k];
             if (fabs(tmp) < small) {
                 skip = 1;
                 break;
             }
         }
         if (!skip) {
             sph_rec(geop, dim, ctr, r, c0, dc0, vcell, &dvol);
         }
         *vol += dvol;
     }
     assert((*vol >= 0.0) && (*vol <= vcell + small));
 
     return;
 }
 
void check_intersected_sphere(int dim, double *xl, double *xr,
                              double *ctr, double r, int *intersected)
{
     int szdim, valid;
     int lh, i, i1, i2;
     double xlr[2][3];
     double dxmin2[3], r2, tmp, *zl, z, dz, dz2, xy2, dx2i1, dxi1, xi1;
 
     szdim = dim * sizeof(double);
     memcpy(xlr[0], xl, (size_t)szdim);
     memcpy(xlr[1], xr, (size_t)szdim);
 
     r2 = r * r;
     for (i = 0; i < dim; i++) {
         dxmin2[i] = xlr[0][i] - ctr[i];
         dxmin2[i] *= dxmin2[i];
         tmp = xlr[1][i] - ctr[i];
         tmp *= tmp;
         if (tmp < dxmin2[i]) dxmin2[i] = tmp;
     }
     *intersected = 0;
 
     for (lh = 0; lh < 2; lh++) {
         zl = xlr[lh];
         for (i = 0; i < dim; i++) {
             dz = zl[i] - ctr[i];
             dz2 = dz * dz;
             xy2 = r2 - dz2;
             if (xy2 < 0.0) continue;
 
             valid = 1;
             for (i1 = 0; i1 < dim; i1++) {
                 if (i1 == i) continue;
                 for (i2 = 0; i2 < dim; i2++) {
                     if ((i2 == i) || (i2 == i1)) continue;
 
                     dx2i1 = xy2 - dxmin2[i2];
                     if (dx2i1 < 0.0) {
                         valid = 0;
                         break;
                     }
                     dxi1 = sqrt(dx2i1);
                     xi1  = ctr[i1] + dxi1;
                     if ((xi1 < xl[i1]) || (xi1 > xr[i1])) {
                         valid = 0;
                         break;
                     }
                 }
             }
             if (valid) {
                 *intersected = 1;   // possiblly intersects
                 break;
             }
         }
         if (*intersected) {
             break;
         }
    }
    return;
 }
 
void poly_rec2d_overlap(int dim, int nnode, double *pts, double *xl, double *dx,
                        int *overlapped)
{
     int    i, k, in, offset, k_touch;
     double *p;
     double xr[2], c4[4][2];
     point_position pos, pos_old;
 
     for (k = 0; k < 2; k++) {
         xr[k] = xl[k] + dx[k];
     }
//   check whether any vertex is in the rectangular
 
     in = 0;
     for (i = 0; i < nnode; i++) {
         p = pts + offset;
         in = 1;
         for (k = 0; k < dim; k++) {
             if ((p[k] <= xl[k]) || (p[k] >= xr[k])) {
                 in = 0;
                 break;
             }
         }
         if (in) break;
     }
     if (in) {
         *overlapped = 1;
         return;
     }
//   check any corner of the rectangular is in the polygon.
 
     for (i = 0; i < 4; i++) {
         p = c4[i];
         in = 0;
         pos = point_in_polygon(pts, nnode, NULL, NULL, p, &k_touch);
//       pos_old = point_in_polygon_old(pts, nnode, p, &k_touch);
//       assert(pos == pos_old);
 
         if (pos == inside) {
             in = 1;;
             break;
         }
     }
     if (in) {
         *overlapped = 1;
         return;
     }
     else {
         *overlapped = 0;
     }
     return;
}
 
 
void poly_rec2d(int geop,
                int dim, int nnode, double *pts, double *xl, double *dx,
                double *vol)
{
//   (pts[0], pts[1]) is the 1st vertex, and so on.
//   pts is assumed is (nnode + 1) * 2 long in memory.
 
     int    err,  i, i1, offset, n, nint;
     double *p0, *p1;
     double xr[3], c5[5][2];
 
     double *xyint;
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
     c5[0][0] = xl[0];
     c5[0][1] = xl[1];
 
     c5[1][0] = xr[0];
     c5[1][1] = xl[1];
 
     c5[2][0] = xr[0];
     c5[2][1] = xr[1];
 
     c5[3][0] = xl[0];
     c5[3][1] = xr[1];
     c5[4][0] = c5[0][0];
     c5[4][1] = c5[0][1];
 
     xyint = NULL;
     err = hull_test(pts, nnode, &(c5[0][0]), 4, &xyint, &nint);
 
     *vol = 0.0;
     if (nint >= 3) {
         if (geop == 1) {
             for (i = 0; i < nint; i++) {
                 i1 = (i + 1) % nint;
                 p0 = xyint + (i + i);
                 p1 = xyint + (i1 + i1);
                 *vol += (p0[0] * p1[1] - p0[1] * p1[0]);
             }
             *vol *= 0.5;
             if (*vol < 0.0) *vol = - (*vol);
         }
         else if (geop == 2) {
             rz_area(nint, xyint, vol);
         }
     }
     if (xyint) free(xyint);
 
     return;
 }
 
void sph_poly3d(int ifinquiry, int dim, double *ctr, double r,
           double *ptrs, int nnode, int nedge, int nface,
           int *nedge_for_face, int *nodelist_for_edge, int *edgelist_for_face,
           double vcell, int *mixed, double *vol)
{
     int    szdim, allinside, alloutside, intersected, npart, npart2;
     int    i, i2, k, offset, nint, nit;
     int    isinside[100], mxed;
     double sixth, r2, tmp, ds2, xlow, xhgh, scale, ddx, x, rx, area;
     double yint, zint, err, dx6, dvol, vol_old;
     double *p0, *p1, *fend, *fm, *array;
     double ctroid[2];
     double ptrs_int[200];
 
     assert(nnode < 100);
     szdim = dim * sizeof(double);
 
     sixth = 0.16666666666666666667;
     r2 = r * r;
 
//   check whether all vertices are in inside the cricle.
 
     offset = 0;
     for (i = 0; i < nnode; i++) {
         p0 = ptrs + offset;
         ds2 = 0.0;
         for (k = 0; k < dim; k++) {
             tmp = p0[k] - ctr[k];
             ds2 += (tmp * tmp);
         }
         if (ds2 < r2) {
             isinside[i] = 1;
         }
         else if (ds2 > r2) {
             isinside[i] = 0;
         }
         else {  // touch
             isinside[i] = -1;
         }
         offset += dim;
     }
     allinside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *mixed = 0;
         *vol   = vcell;
         return;
     }
     alloutside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i]) {
             alloutside = 0;
             break;
         }
     }
     if (alloutside) {
 
/**********
//       conservatively quarantee
 
         intersected = 0;
         offset = 0;
         for (i = 0; i < nnode; i++) {
             p0 = ptrs + offset;
             ds2 = 0.0;
             for (k = 0; k < dim; k++) {
                 if ((-ctr[k] < p0[k]) && (p0[k] < ctr[k])) {
                     intersected = 1;
                     break;
                 }
             }
             if (intersected) {
                 break;
             }
             offset += dim;
         }
*********/
         *mixed = 0;
         *vol   = 0.0;
         return;
     }
     *mixed = 1;
     if (ifinquiry) {
         *vol = small;
         return;
     }
     xlow = ptrs[0];
     xhgh = xlow;
     offset = 0;
     for (i = 0; i < nnode; i++) {
         p0 = ptrs + offset;
         if (p0[0] < xlow) {
             xlow = p0[0];
         }
         else if (p0[0] > xhgh) {
             xhgh = p0[0];
         }
         offset += dim;
     }
     scale = xhgh - xlow;
 
     ddx = scale;
     for (i = 0; i < nnode; i++) {
         p0 = ptrs + (i * dim);
         for (k = i + 1; k < nnode; k++) {
             p1 = ptrs + (k * dim);
             tmp = fabs(p1[0] - p0[0]);
             ddx = MIN(ddx, tmp);
         }
     }
     ddx = MAX(ddx, 0.001 * scale);
     npart = scale/ddx;
 
//   Simpson integral
 
    fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
    assert(fend);
    fm    = fend + (npart + 1);
    ddx   = (xhgh - xlow)/(double)npart;
    x     = xlow - ddx;
    for (i = 0; i <= npart; i++) {
        x += ddx;
 
//      find the polygon resulted from the intersection between
//      the plane x = x and the polyhedron.
 
        plane_poly3d(dim, x, ptrs, nnode, nedge, nface, nedge_for_face,
                     nodelist_for_edge, edgelist_for_face,
                     &nint, ptrs_int);
 
        if (nint < 3) {
            area = 0.0;
        }
        else {
            tmp = x - ctr[0];
            rx  = sqrt(MAX(r2 - tmp * tmp, 0.0));
            sph_poly2d(0, 1, 2, ctr + 1, rx, ptrs_int, nint, &mxed, &area, ctroid);
        }
        fend[i] = area;
    }
    err = 1.0;
    nit = 0;
    while (err > tol) {
         ddx = (xhgh - xlow)/(double)npart;
         dx6 = sixth * ddx;
         *vol = 0.0;
         x = xlow - 0.5 * ddx;
         for (i = 0; i < npart; i++) {
             x += ddx;
 
             plane_poly3d(dim, x, ptrs, nnode, nedge, nface, nedge_for_face,
                          nodelist_for_edge, edgelist_for_face,
                          &nint, ptrs_int);
             if (nint < 3) {
                 area = 0.0;
             }
             else {
                 tmp = x - ctr[0];
                 rx  = sqrt(MAX(r2 - tmp * tmp, 0.0));
                 sph_poly2d(0, 1, 2, ctr + 1, rx, ptrs_int, nint, &mxed, &area, ctroid);
             }
             fm[i] = area;
             dvol  = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
             *vol += dvol;
         }
         if (nit) {
             err = fabs(*vol - vol_old)/vcell;
         }
         if (err > tol) {
             npart2 = npart + npart;
             array = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
             assert(array);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 array[i2]   = fend[i];
                 array[i2+1] = fm[i];
             }
             array[npart2] = fend[npart];
             free(fend);
             fend = array;
             fm   = fend + (npart2 + 1);
 
             npart = npart2;
             nit++;
         }
         vol_old = *vol;
     }
     free(fend);
 
    return;
 }
 
void ellipse_poly3d(int ifinquiry, int dim, double *ctr, double *r,
           double *ptrs, int nnode, int nedge, int nface,
           int *nedge_for_face, int *nodelist_for_edge, int *edgelist_for_face,
           double vcell, int *mixed, double *vol)
{
//   Simpson formula
 
     int    isinside[100], intersected, mxed, found, done;
     int    szdim, allinside, alloutside, npart, npart2;
     int    dim0, dim1, i0, i, i2, k, s, offset, nint, nit;
     int    np, ns;
     double dxmin, sixth, ds2, xs[50];
     double r2, r0, tmp, xlow, xhgh, scale, ddx, x, rx, area, areamx;
     double factors[2], factor, yint, zint, err, dx6, dvol, vol_old;
     double *p0, *p1, *fend, *fm, *array;
     double ctroid[2], rp[2], ctrp[3], ptrs_int[300];
     int    ngrp, g;
     double dxgrp, xlow_grp, xhgh_grp, vol_grp;
 
     assert(nnode < 100);
     szdim = dim * sizeof(double);
     sixth = 0.16666666666666666667;
 
     if (major_axis_sav == 0) {
         dim0 = 1;
         dim1 = 2;
     }
     else if (major_axis_sav == 2) {
         dim0 = 0;
         dim1 = 1;
     }
//   shrink only to check whether inside or outside
 
     r0 = r[0];
     for (i = 1; i < dim; i++) {
         if (r[i] < r0) {
             r0 = r[i];
         }
     }
//   scale
 
     r2     = r0 * r0;
     factor = 1.0;
     for (i = 0; i < dim; i++) {
         tmp     = r[i]/r0;
         factor *= tmp;
         tmp     = 1.0/tmp;
         for (k = 0; k < nnode; k++) {
             p0 = ptrs     + (k * dim);
             p1 = ptrs_int + (k * dim);
             p1[i] = tmp * p0[i];
         }
         ctrp[i] = tmp * ctr[i];
     }
//   check whether all vertices are inside the sphere.
 
     offset = 0;
     for (i = 0; i < nnode; i++) {
         p0 = ptrs_int + offset;
         ds2 = 0.0;
         for (k = 0; k < dim; k++) {
             tmp = p0[k] - ctrp[k];
             ds2 += (tmp * tmp);
         }
         if (ds2 < r2) {
             isinside[i] = 1;
         }
         else if (ds2 > r2) {
             isinside[i] = 0;
         }
         else {  // touch
             isinside[i] = -1;
         }
         offset += dim;
     }
     allinside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *mixed = 0;
         *vol   = vcell;
         return;
     }
     alloutside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i]) {
             alloutside = 0;
             break;
         }
     }
     if (alloutside) {
 
/*******
//       conservatively guarantee
 
         intersected = 0;
         offset = 0;
         for (i = 0; i < nnode; i++) {
             p0 = ptrs_int + offset;
             ds2 = 0.0;
             for (k = 0; k < dim; k++) {
                 if ((-ctrp[k] < p0[k]) && (p0[k] < ctrp[k])) {
                     intersected = 1;
                     break;
                 }
             }
             if (intersected) {
                 break;
             }
             offset += dim;
         }
**********/
         *mixed = 0;
         *vol   = 0.0;
         return;
     }
     *mixed = 1;
     if (ifinquiry) {
         *vol = small;
         return;
     }
     r0         = r[dim0];
     factors[0] = r0/r[dim0];
     factors[1] = r0/r[dim1];
     factor     = 1.0/(factors[0] * factors[1]);
 
     np = 0;
     for (i = 0; i < nnode; i++) {
         p0 = ptrs + (i * dim);
         found = 0;
         for (k = 0; k < np; k++) {
              if (xs[k] == p0[major_axis_sav]) {
                 found = 1;
                 break;
             }
         }
         if (!found) {
             xs[np] = p0[major_axis_sav];
             np++;
         }
     }
     for (i = 0; i < np; i++) {
         xlow = xs[i];
         i0 = i;
         for (k = i + 1; k < np; k++) {
             if (xlow > xs[k]) {
                 xlow = xs[k];
                 i0 = k;
             }
         }
         if (i0 != i) {
             tmp = xs[i];
             xs[i] = xs[i0];
             xs[i0] = tmp;
         }
     }
     ns = np - 1;
     areamx = 0.0;
     dxmin = small;
 
     *vol = 0.0;
     ngrp = 10;
     for (s = 0; s < ns; s++) {
         xlow_grp = xs[s];
         xhgh_grp = xs[s+1];
         if (xhgh_grp <= xlow_grp) continue;
 
         dxgrp = (xhgh_grp - xlow_grp)/(double)ngrp;
         xlow  = xlow_grp - dxgrp;
 
         for (g = 0; g < ngrp; g++) {
             xlow += dxgrp;
             xhgh  = xlow + dxgrp;
 
             npart = 8;
 
         //  Simpson integral
 
             fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
             assert(fend);
 
             fm    = fend + (npart + 1);
             ddx   = (xhgh - xlow)/(double)npart;
             x     = xlow - ddx;
             for (i = 0; i <= npart; i++) {
                 x += ddx;
 
                 tmp = (x - ctr[major_axis_sav])/r[major_axis_sav];
                 tmp = 1.0 - tmp * tmp;
                 if (tmp <= 0.0) {
                     area = 0.0;
                 }
                 else {
                     tmp   = sqrt(tmp);
                     rp[0] = tmp * r[dim0];
                     rp[1] = tmp * r[dim1];
 
             //      find the polygon resulted from the intersection between
             //      the plane x = x and the polyhedron.
 
                     plane_poly3d(dim, x, ptrs, nnode, nedge, nface, nedge_for_face,
                                  nodelist_for_edge, edgelist_for_face,
                                  &nint, ptrs_int);
 
                     if (nint < 3) {
                         area = 0.0;
                     }
                     else {
             //          scaling
             //
                         r0 = rp[0];
                         for (k = 0; k < nint; k++) {
                             p0 = ptrs_int + (k + k);
                             for (i0 = 0; i0 < 2; i0++) {
                                 p0[i0] *= factors[i0];
                             }
                         }
                         ctrp[0] = factors[0] * ctr[dim0];
                         ctrp[1] = factors[1] * ctr[dim1];
                         sph_poly2d(0, 1, 2, ctrp, r0, ptrs_int, nint, &mxed, &area, ctroid);
                         area *= factor;
                     }
                 }
                 fend[i] = area;
                 areamx = MAX(areamx, area);
             }
             vol_old = 1.0e+30;
             done = 0;
             nit = 0;
             while (!done) {
                  ddx = (xhgh - xlow)/(double)npart;
                  dx6 = sixth * ddx;
                  vol_grp = 0.0;
                  x = xlow - 0.5 * ddx;
                  for (i = 0; i < npart; i++) {
                      x += ddx;
 
                      tmp = (x - ctr[major_axis_sav])/r[major_axis_sav];
                      tmp = 1.0 - tmp * tmp;
                      if (tmp <= 0.0) {
                          area = 0.0;
                      }
                      else {
                          tmp   = sqrt(tmp);
                          rp[0] = tmp * r[dim0];
                          rp[1] = tmp * r[dim1];
 
                          plane_poly3d(dim, x, ptrs, nnode, nedge, nface, nedge_for_face,
                                       nodelist_for_edge, edgelist_for_face,
                                       &nint, ptrs_int);
 
                          if (nint < 3) {
                              area = 0.0;
                          }
                          else {
             //               scaling
             //
                              r0 = rp[0];
                              for (k = 0; k < nint; k++) {
                                  p0 = ptrs_int + (k + k);
                                  for (i0 = 0; i0 < 2; i0++) {
                                      p0[i0] *= factors[i0];
                                  }
                              }
                              ctrp[0] = factors[0] * ctr[dim0];
                              ctrp[1] = factors[1] * ctr[dim1];
                              sph_poly2d(0, 1, 2, ctrp, r0, ptrs_int, nint, &mxed, &area, ctroid);
                              area *= factor;
                          }
                      }
                      fm[i] = area;
                      dvol  = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
                      vol_grp += dvol;
                      areamx = MAX(areamx, area);
                  }
                  err   = fabs(vol_grp - vol_old)/vcell;
                  dxmin = tol * (vcell/areamx);
                  dxmin = 0.1 * sqrt(dxmin);
 
                  if ((err > tol) && (ddx > dxmin)) {
                      npart2 = npart + npart;
                      array  = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
                      assert(array);
                      for (i = 0; i < npart; i++) {
                          i2 = i + i;
                          array[i2]   = fend[i];
                          array[i2+1] = fm[i];
                      }
                      array[npart2] = fend[npart];
                      free(fend);
                      fend = array;
                      fm   = fend + (npart2 + 1);
 
                      npart = npart2;
                      nit++;
                  }
                  else {
                      done = 1;
                  }
                  vol_old = vol_grp;
             }
             free(fend);
 
             *vol += vol_grp;
         }
     }
     return;
 }
 
 
void conic_poly3d(int ifinquiry, int dim,
           double x0, double r0, double x1, double r1,
           double *pts, int nnode, int nedge, int nface,
           int *nedge_for_face, int *nodelist_for_edge, int *edgelist_for_face,
           double vcell, int *mixed, double *vol)
{
//   Simpson formula
 
     int    mxed, found, done;
     int    szdim, npart, npart2;
     int    i0, i, i2, k, s, offset, nint, nit;
     int    np, ns;
     double mytol, dxinv, dxmin, sixth, ds2, xs[50], vols[50];
     double tmp, drdx, r, xlow, xhgh, scale, ddx, x, area, areamx;
     double yint, zint, err, dx6, dvol, vol_old;
     double *p0, *p1, *fend, *fm, *array;
     double ctr[2], ctroid[2], pts_int[300];
     int    ngrp, g;
     double dxgrp, xlow_grp, xhgh_grp, vol_grp;
 
     assert(nnode < 100);
     szdim = dim * sizeof(double);
     sixth = 0.16666666666666666667;
 
//   for some reason, Simpson integral convers too slowly. I temporarily
//   set a larger tolerance.
//
     mytol = 10.0 * tol;
 
 
     ctr[0] = 0.0;
     ctr[1] = 0.0;
     dxinv  = 1.0/(x1 - x0);
 
     np = 0;
     for (i = 0; i < nnode; i++) {
         p0 = pts + (i * dim);
         x  = p0[0];
         if (x < x0) {
             x = x0;
         }
         else if (x > x1) {
             x = x1;
         }
         found = 0;
         for (k = 0; k < np; k++) {
              if (fabs(xs[k] - x)*dxinv < small) {
                 found = 1;
                 break;
             }
         }
         if (!found) {
             xs[np] = x;
             np++;
         }
     }
     for (i = 0; i < np; i++) {
         xlow = xs[i];
         i0 = i;
         for (k = i + 1; k < np; k++) {
             if (xlow > xs[k]) {
                 xlow = xs[k];
                 i0 = k;
             }
         }
         if (i0 != i) {
             tmp = xs[i];
             xs[i] = xs[i0];
             xs[i0] = tmp;
         }
     }
     ns    = np - 1;
     dxmin = small;
     drdx  = (r1 - r0)/(x1 - x0);
 
     areamx = 0.0;
     *vol = 0.0;
     for (s = 0; s < ns; s++) {
         xlow_grp  = xs[s];
         xhgh_grp  = xs[s+1];
         vols[s]   = 0.0;
         if (xhgh_grp <= xlow_grp) continue;
 
         ngrp = 10;
         dxgrp = (xhgh_grp - xlow_grp)/(double)ngrp;
 
         xlow = xlow_grp - dxgrp;
         for (g = 0; g < ngrp; g++) {
             xlow += dxgrp;
             xhgh  = xlow + dxgrp;
 
             vol_grp = 0.0;
             npart = 8;
 
    //       Simpson integral
 
             fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
             assert(fend);
 
             fm  = fend  + (npart + 1);
             ddx = (xhgh - xlow)/(double)npart;
             x   = xlow  - ddx;
             for (i = 0; i <= npart; i++) {
                 x += ddx;
                 r  = r0 + drdx *(x - x0);
 
    //           find the polygon resulted from the intersection between
    //           the plane x = x and the polyhedron.
 
                 plane_poly3d(dim, x, pts, nnode, nedge, nface, nedge_for_face,
                              nodelist_for_edge, edgelist_for_face,
                              &nint, pts_int);
 
                 if (nint < 3) {
                     area = 0.0;
                 }
                 else {
                     sph_poly2d(0, 1, 2, ctr, r, pts_int, nint, &mxed, &area, ctroid);
                 }
                 fend[i] = area;
                 areamx  = MAX(areamx, area);
             }
             vol_old = 1.0e+30;
             done = 0;
             nit = 0;
             while (!done) {
                  ddx = (xhgh - xlow)/(double)npart;
                  dx6 = sixth * ddx;
                  vol_grp = 0.0;
                  x = xlow - 0.5 * ddx;
                  for (i = 0; i < npart; i++) {
                      x += ddx;
                      r  = r0 + drdx *(x - x0);
 
                      plane_poly3d(dim, x, pts, nnode, nedge, nface, nedge_for_face,
                                   nodelist_for_edge, edgelist_for_face,
                                   &nint, pts_int);
 
                      if (nint < 3) {
                          area = 0.0;
                      }
                      else {
                          sph_poly2d(0, 1, 2, ctr, r, pts_int, nint, &mxed, &area, ctroid);
                      }
                      fm[i] = area;
                      dvol  = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
                      vol_grp += dvol;
                      areamx  = MAX(areamx, area);
                  }
                  err   = fabs(vol_grp - vol_old)/vcell;
                  dxmin = tol * (vcell/areamx);
                  dxmin = sqrt(dxmin);
 
    //            if ((err > tol) && (ddx > dxmin)) {
 
                  if (err > mytol) {
                      npart2 = npart + npart;
                      array  = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
                      assert(array);
                      for (i = 0; i < npart; i++) {
                          i2 = i + i;
                          array[i2]   = fend[i];
                          array[i2+1] = fm[i];
                      }
                      array[npart2] = fend[npart];
                      free(fend);
                      fend = array;
                      fm   = fend + (npart2 + 1);
 
                      npart = npart2;
                      nit++;
                  }
                  else {
                      done = 1;
                  }
                  vol_old = vol_grp;
              }
              free(fend);
 
              *vol += vol_grp;
          }
     }
     return;
 }
 

void plane_poly3d(int dim, double xint, double *ptrs,
                  int nnode, int nedge, int nface,
                  int *nedge_for_face, int *nodelist_for_edge,
                  int *edgelist_for_face,
                  int *nint, double *ptrs_int)
{
     int    szdim2, skip, failed;
     int    f, e, n0, n1, ne, i0, i, i1, i2, k, offset;
     int    included[20], edge_intersected[20];
     int    *edgelist, *nodelist;
 
     double t, yint, zint, normz, mynorm;
     double dp3[3], dr2[20][2];
     double *dp, *p0, *p1, *p2, *pt, *pi;
     double ptrs_edge_int[20][2];
     double dxy0[2], dxy1[2], tmpxy[2];
     double dxmxinv, a, sina, cosa;
 
     szdim2 = 2 * sizeof(double);
 
 
//   intersection between each edge and the plane x = xint.
 
     for (i = 0; i < nnode; i++) {
         included[i] = 0;
     }
     offset = 0;
     dxmxinv = 0.0;
     for (e = 0; e < nedge; e++) {
         edge_intersected[e] = 0;
         pt = ptrs_edge_int[e];
 
         nodelist = nodelist_for_edge + offset;
         n0 = nodelist[0];
         n1 = nodelist[1];
         p0 = ptrs + (dim * n0);
         p1 = ptrs + (dim * n1);
         for (k = 0; k < dim; k++) {
             dp3[k] = p1[k] - p0[k];
             if (dxmxinv < fabs(dp3[k])) {
                 dxmxinv = fabs(dp3[k]);
             }
         }
         if (dp3[0] == 0.0) {
             skip = 1;
         }
         else {
             skip = 0;
         }
         if (!skip) {
             t = (xint - p0[0])/dp3[0];
             if ((t < 0.0) || (t > 1.0)) {
                 skip = 1;
             }
         }
         if (!skip) {
             if (t == 0.0) {
                 if (!included[n0]) {
                     memcpy(pt, p0 + 1, (size_t)szdim2);
                     included[n0] = 1;
                     edge_intersected[e] = 1;
                 }
             }
             else if (t == 1.0) {
                 if (!included[n1]) {
                     memcpy(pt, p1 + 1, (size_t)szdim2);
                     included[n1] = 1;
                     edge_intersected[e] = 1;
                 }
             }
             else {
                 yint  = p0[1] + dp3[1] * t;
                 zint  = p0[2] + dp3[2] * t;
                 pt[0] = yint;
                 pt[1] = zint;
                 edge_intersected[e] = 1;
             }
         }
         offset += 2;
     }
     dxmxinv = 1.0/dxmxinv;
     for (i = 0; i < nedge; i++) {
         included[i] = 0;
     }
     *nint = 0;
     offset = 0;
     for (f = 0; f < nface; f++) {
 
         ne       = nedge_for_face[f];
         edgelist = edgelist_for_face + offset;
 
         for (i = 0; i < ne; i++) {
             e  = edgelist[i];
             if (!edge_intersected[e] || included[e]) continue;
 
             pi  = ptrs_edge_int[e];
 
//           check whether pi is cery close to the one in the list.
//
             skip = 0;
             for (i1 = 0; i1 < *nint; i1++) {
                 pt  = ptrs_int + (i1 + i1);
                 a  = (fabs(pt[0]-pi[0]) + fabs(pt[1]-pi[1])) * dxmxinv;
                 if (a < small) {
                     skip = 1;
                     break;
                 }
             }
             if (!skip) {
                 pt  = ptrs_int + (*nint + *nint);
                 memcpy(pt, pi, (size_t)szdim2);
                 (*nint)++;
             }
             included[e] = 1;
         }
         offset += ne;
     }
     if (*nint < 4) {
         return;
     }
//   put nodes into the correct order
 
     nodelist = (int *) malloc((*nint) * sizeof(int));
     for (i = 0; i < *nint; i++) { 
         nodelist[i] = i;
     } 
     mychull_sort(ptrs_int, nodelist, *nint, &failed);
 
     free(nodelist);

     return;
 }
 

void sph_poly2d(int ifinquiry, int geop, int dim, double *ctr, double r,
                double *ptrs, int nnode,
                int *mixed, double *vol, double *ctroid)
{
//   ctroid[2] is the centroid of the overlap if there is an overlap.
 
     int   intersected;
     int   allinside, alloutside, ifescape, already, included[100];
     int   n0, i, i0, i1, k, k1, offset, ncr_tot, ncr, i_inside, nint, szdim;
     int   nn_poly, nn_cir;
     int   isinside[100], ncross[100], i0_cross[100], i1_cross[100];
     int   first_cross[100];
     double third, two3rd;
     double r2, r3, ds2, tmp, a, b, c, delt2, rt, t0, t1;
     double sina, vol_poly, vol_cord;
     double myctr[2], dx[2], delx[2], dp[2];
     double cross[100], ptrs_res[200], ptrs_cir[200];
     double *p0, *p1, *pt, *ptc;
     double dv, r0, z0, r1, z1, dz, rc, zc, r01, r02, r12, r0r1;
     double area_pie, area_tri, ha, sinha, cosha, xcarea, dxc_pie;
     double ds, factor, pm, rc_pie, rc_tri;
 
     double coords_tri[6], *pt0, *pt1, myarea;
 
     assert(nnode < 100);
     szdim = dim * sizeof(double);
 
     third  = 0.33333333333333333333333;
     two3rd = third + third;
 
     r2 = r * r;
     r3 = r * r2;
 
     *vol = 0.0;
 
//   ctroid will be the centroid of the overlap if there is a overlap.
     memcpy(ctroid, ctr, (size_t)szdim);
 
//   check whether all vertices are in inside the cricle.
 
     offset = 0;
     for (i = 0; i < nnode; i++) {
         p0 = ptrs + offset;
         ds2 = 0.0;
         for (k = 0; k < dim; k++) {
             tmp = p0[k] - ctr[k];
             ds2 += (tmp * tmp);
         }
         if (ds2 < r2) {
             isinside[i] = 1;
         }
         else if (ds2 > r2) {
             isinside[i] = 0;
         }
         else {  // touch
             isinside[i] = -1;
         }
         offset += dim;
     }
     allinside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *mixed = 0;
         *vol   = 0.0;
         ctroid[0] = 0.0;
         ctroid[1] = 0.0;
 
         if (geop == 1) {
             for (i = 0; i < nnode; i++) {
                 i1 = (i + 1) % nnode;
                 p0 = ptrs + (dim * i);
                 p1 = ptrs + (dim * i1);
                 dv = p0[0] * p1[1] - p1[0] * p0[1];
                 ctroid[0] += ((p0[0] + p1[0]) * dv);
                 ctroid[1] += ((p0[1] + p1[1]) * dv);
                 *vol += dv;
             }
             *vol *= 0.5;
             if (*vol != 0.0) {
                 factor = 1.0/(6.0 * (*vol));
                 ctroid[0] *= factor;
                 ctroid[1] *= factor;
             }
             if (*vol < 0.0) *vol = - (*vol);
         }
         else if (geop == 2) {
             for (i = 0; i < nnode; i++) {
                 i1 = (i + 1) % nnode;
                 p0 = ptrs + (dim * i);
                 p1 = ptrs + (dim * i1);
 
                 r0 = p0[0];
                 z0 = p0[1];
                 r1 = p1[0];
                 z1 = p1[1];
 
                 dz   = z1 - z0;
                 r02  = r0 * r0;
                 r12  = r1 * r1;
                 r0r1 = r0 * r1;
 
//               dv = 2 * integral (r dr dz)
 
                 dv = dz *(r12 + r0r1 + r02);
                 *vol += dv;
 
                 rc = dz *(r02 * r0 + r02 * r1 + r0 * r12 + r12 * r1);
                 ctroid[0] += rc;
 
                 r01 = r0 + r1;
                 zc = dz *(r01 * r01 *(z0 + z1) + 2.0 *(r02 * z0 + r12 * z1));
                 ctroid[1] += zc;
             }
             *vol /= 3.0;
 
             if (*vol < 0.0) {
                 *vol = - (*vol);
                 ctroid[0] = -ctroid[0];
                 ctroid[1] = -ctroid[1];
             }
         }
         else {
             printf("ERROR: geop != 1, or 2, not yet.\n");
         }
         return;
     }
     alloutside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i]) {
             alloutside = 0;
             break;
         }
     }
     if (alloutside) {
 
//       check whether the circle intersects with any edge
 
         intersected = 0;
         for (i = 0; i < nnode; i++) {
             i1 = (i + 1) % nnode;
             p0 = ptrs + (dim * i);
             p1 = ptrs + (dim * i1);
 
             for (k = 0; k < dim; k++) {
                 dx[k]   = p1[k] - p0[k];
                 delx[k] = p0[k] - ctr[k];
             }
             a = 0.0;
             b = 0.0;
             c = 0.0;
             for (k = 0; k < dim; k++) {
                 a += (dx[k]   * dx[k]);
                 b += (dx[k]   * delx[k]);
                 c += (delx[k] * delx[k]);
             }
             b  = b + b;
             c -= r2;
             delt2 = b * b - 4.0 * a * c;
             if ((delt2 < 0.0) || (a == 0.0)) continue;
 
             rt = sqrt(delt2);
             tmp = 0.5 /a;
             t0  = (-b - rt) * tmp;
             t1  = (-b + rt) * tmp;
 
             if (((t0 >= 0.0) && (t0 <= 1.0)) || ((t1 >= 0.0) && (t1 <= 1.0))) {
                 intersected = 1;
                 break;
             }
         }
         if (!intersected) {
             *mixed = 0;
             *vol = 0.0;
             return;
         }
     }
//   find intersection between circle and line segement
 
     ncr_tot = 0;
     for (i = 0; i < nnode; i++) {
         i1 = (i + 1) % nnode;
 
//          if ((isinside[i] == 1) && (isinside[i1] == 1)) continue;
 
         p0 = ptrs + (dim * i);
         p1 = ptrs + (dim * i1);
 
         for (k = 0; k < dim; k++) {
             dx[k]   = p1[k] - p0[k];
             delx[k] = p0[k] - ctr[k];
         }
         a = 0.0;
         b = 0.0;
         c = 0.0;
         for (k = 0; k < dim; k++) {
             a += (dx[k] * dx[k]);
             b += (dx[k] * delx[k]);
             c += (delx[k] * delx[k]);
         }
         b  = b + b;
         c -= r2;
         delt2 = b * b - 4.0 * a * c;
 
         first_cross[i] = -1;
         if ((delt2 < 0.0) || (a == 0.0)) {
             ncross[i] = 0;
             continue;
         }
         rt = sqrt(delt2);
         tmp = 0.5 / a;
         t0  = (-b - rt) * tmp;
         t1  = (-b + rt) * tmp;
 
         if (fabs(t0 - t1) < small) {
             ncross[i] = 0;
             continue;
         }
         ncross[i] = 0;
         if ((t0 >= 0.0) && (t0 <= 1.0)) {
             ncross[i]++;
             cross[ncr_tot]    = t0;
             i0_cross[ncr_tot] = i;
             i1_cross[ncr_tot] = i1;
             first_cross[i]    = ncr_tot;
             ncr_tot++;
         }
         if ((t1 >= 0.0) && (t1 <= 1.0)) {
             if (ncross[i] == 0) {
                 first_cross[i] = ncr_tot;
             }
             ncross[i]++;
             cross[ncr_tot]    = t1;
             i0_cross[ncr_tot] = i;
             i1_cross[ncr_tot] = i1;
             ncr_tot++;
         }
     }
     if (ncr_tot < 2) {
         *mixed = 0;
         return;
     }
     if (ifinquiry) {
         *mixed = 1;
         return;
     }
     i_inside = -1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i] == 1) {
             i_inside = i;
             break;
         }
     }
     for (i = 0; i < nnode; i++) {
         included[i] = 0;
     }
     ctroid[0] = 0.0;
     ctroid[1] = 0.0;
     nn_cir    = 0.0;
 
     *mixed = 1;
     if (i_inside >= 0) {
         p0 = ptrs + (i_inside * dim);
         memcpy(ptrs_res, p0, (size_t)szdim);
         included[i_inside] = 1;
         nn_poly = 1;
         pt     = ptrs_res + dim;
         ptc    = ptrs_cir;
 
         for (i = 0; i < nnode; i++) {
             i0 = (i_inside + i) % nnode;
             i1 = (i0 + 1) % nnode;
////             if ((ncross[i0] == 0) && (isinside[i1] == 1)) {
             if ((ncross[i0] == 0) && isinside[i1]) {
                 if (!included[i1]) {
                     p1 = ptrs + (i1 * dim);
                     memcpy(pt, p1, (size_t)szdim);
                     nn_poly++;
                     pt += dim;
                 }
                 included[i1] = 1;
             }
             else if (ncross[i0] == 1) {
                 ncr = first_cross[i0];
                 t0  = cross[ncr];
                 p0  = ptrs + (i0 * dim);
                 p1  = ptrs + (i1 * dim);
                     pt[0] = p0[0] + (p1[0] - p0[0]) * t0;
                     pt[1] = p0[1] + (p1[1] - p0[1]) * t0;
                     memcpy(ptc, pt, (size_t)szdim);
                     nn_poly++;
                     nn_cir++;
                     ptc += dim;
                     pt  += dim;
 
/////              if ((isinside[i1] == 1) && !included[i1]) {
                   if (isinside[i1] && !included[i1]) {
                     memcpy(pt, p1, (size_t)szdim);
                     included[i1] = 1;
                     nn_poly++;
                     pt += dim;
                 }
             }
             else if (ncross[i0] == 2) {
                 ncr = first_cross[i0];
                 t0  = cross[ncr];
                 p0  = ptrs + (i0 * dim);
                 p1  = ptrs + (i1 * dim);
 
                     pt[0] = p0[0] + (p1[0] - p0[0]) * t0;
                     pt[1] = p0[1] + (p1[1] - p0[1]) * t0;
                     memcpy(ptc, pt, (size_t)szdim);
                     nn_poly++;
                     nn_cir++;
                     ptc += dim;
                     pt  += dim;
 
                 t0 = cross[ncr+1];
                     pt[0] = p0[0] + (p1[0] - p0[0]) * t0;
                     pt[1] = p0[1] + (p1[1] - p0[1]) * t0;
                     memcpy(ptc, pt, (size_t)szdim);
                     nn_poly++;
                     nn_cir++;
                     ptc += dim;
                     pt  += dim;
             }
         }
         nint = (pt - ptrs_res)/dim;
         vol_poly = 0.0;
 
         if (geop == 1) {
             myctr[0] = 0.0;
             myctr[1] = 0.0;
             for (i = 0; i < nn_poly; i++) {
                 i1 = (i + 1) % nn_poly;
                 p0 = ptrs_res + (dim * i);
                 p1 = ptrs_res + (dim * i1);
                 dv = p0[0] * p1[1] - p1[0] * p0[1];
                 myctr[0] += ((p0[0] + p1[0]) * dv);
                 myctr[1] += ((p0[1] + p1[1]) * dv);
                 vol_poly += dv;
             }
             vol_poly *= 0.5;
             if (vol_poly != 0.0) {
                 factor = 1.0/(6.0 * vol_poly);
                 myctr[0] *= factor;
                 myctr[1] *= factor;
             }
             if (vol_poly < 0.0) vol_poly = - vol_poly;
             for (i = 0; i < dim; i++) {
                 ctroid[i] += (vol_poly * myctr[i]);
             }
         }
         else if (geop == 2) {
             for (i = 0; i < nn_poly; i++) {
                 i1 = (i + 1) % nn_poly;
                 p0 = ptrs_res + (dim * i);
                 p1 = ptrs_res + (dim * i1);
                 r0 = p0[0];
                 z0 = p0[1];
                 r1 = p1[0];
                 z1 = p1[1];
 
                 dz   = z1 - z0;
                 r02  = r0 * r0;
                 r12  = r1 * r1;
                 r0r1 = r0 * r1;
 
//               dv = 2 * integral (r dr dz)
 
                 dv = dz *(r12 + r0r1 + r02);
                 vol_poly += dv;
 
                 rc = dz *(r02 * r0 + r02 * r1 + r0 * r12 + r12 * r1);
                 myctr[0] += rc;
 
                 r01 = r0 + r1;
                 zc = dz *(r01 * r01 *(z0 + z1) + 2.0 *(r02 * z0 + r12 * z1));
                 myctr[1] += zc;
             }
             vol_poly /= 3.0;
 
             for (i = 0; i < dim; i++) {
                 ctroid[i] += (vol_poly * myctr[i]);
             }
             if (vol_poly < 0.0) {
                 vol_poly = - vol_poly;
             }
         }
     }
     else {
         vol_poly = 0.0;
     }
//   calculate the volume of circular sections
 
     nint = (ptc - ptrs_cir)/dim;
 
     ifescape = 0;
     *vol = 0.0;
 
     for (i = 0; i < nn_cir-1; i++) {
 
         if (!ifescape) {
             i1 = i + 1;
             p0 = ptrs_cir + (dim * i);
             p1 = ptrs_cir + (dim * i1);
             ds = fabs(p1[0] - p0[0]) + fabs(p1[1] - p0[1]);
             if (ds < small * r) {
//           if ((p0[0] == p1[0]) && (p0[1] == p1[1])) {
                 ifescape = 1 - ifescape;
                 continue;
             }
             if (geop == 1) {
                 ds2 = 0.0;
                 for (k = 0; k < dim; k++) {
                     tmp  = p1[k] - p0[k];
                     ds2 += (tmp * tmp);
                 }
                 c = sqrt(ds2);
                 a = 2.0 * asin(0.5 * c / r);
                 sina = sin(a);
                 assert(sina >= 0.0);
                 area_pie = 0.5 * r2 * a;
                 area_tri = 0.5 * r2 * sina;
                 vol_cord = area_pie - area_tri;
 
                 *vol += vol_cord;
 
//               center of the "pie"
 
                 sinha   = 0.5 * c / r;
                 cosha   = sqrt(MAX(1.0 - sinha * sinha, 0.0));
                 xcarea  = two3rd * r3 * sinha;
                 dxc_pie = xcarea  / area_pie;
                 factor  = dxc_pie / (r * cosha);
 
/*********
//               centroid of the triangle (ctr[k], p0, p1);
                 memcpy(coords_tri,     ctr, (size_t)szdim);
                 memcpy(coords_tri + 2, p0,  (size_t)szdim);
                 memcpy(coords_tri + 4, p1,  (size_t)szdim);
                 myarea   = 0.0;
                 myctr[0] = 0.0;
                 myctr[1] = 0.0;
                 for (k = 0; k < 3; k++) {
                     k1 = (k + 1) % 3;
                     pt0 = coords_tri + (k  + k);
                     pt1 = coords_tri + (k1 + k1);
                     dv = pt0[0] * pt1[1] - pt1[0] * pt0[1];
                     myctr[0] += ((pt0[0] + pt1[0]) * dv);
                     myctr[1] += ((pt0[1] + pt1[1]) * dv);
                     myarea += dv;
                 }
                 myarea *= 0.5;
                 if (myarea != 0.0) {
                     tmp = 1.0/(6.0 * myarea);
                     myctr[0] *= tmp;
                     myctr[1] *= tmp;
                 }
                 if (myarea < 0.0) myarea = - myarea;
 
//               myarea should be the same as area_tri
//               myctr should be the same as rc_tri later.
*************/
 
 
                 for (k = 0; k < dim; k++) {
                     pm = 0.5 *(p0[k] + p1[k]);
                     rc_pie = ctr[k] + factor * (pm - ctr[k]);
                     rc_pie *= area_pie;
 
//                   center of the "triangle", (ctr, p0, p1)
 
                     rc_tri = third *(ctr[k] + p0[k] + p1[k]);
 
                     rc_tri *= area_tri;
 
//                   add the product of centroid and area for the "cord"
 
                     ctroid[k] += (rc_pie - rc_tri);
                 }
             }
             else if (geop == 2) {
                 cal_vf_cyl_cord(r, ctr, p0, p1, &vol_cord);
                 *vol += vol_cord;
             }
         }
         ifescape = 1 - ifescape;
     }
     *vol += vol_poly;
 
     if (*vol  != 0.0) {
         factor = 1.0/(*vol);
         for (k = 0; k < dim; k++) {
             ctroid[k] *= factor;
         }
     }
     return;
 }
 
void cal_vf_cyl_cord(double radius, double *ctr, double *p0, double *p1,
                     double *vol)
{
//  This is for geop = 2, cylindral coordinate
//
//   vol = int_y0^y1 [  (R2 - y * y) - (x0 + slope (y - y0) )^2 ]  dy
 
    int i;
    double mypi, third;
    double r0[2], r1[2], slope, s2, x2, y2, xy, r2, c0, c1, c2, tmp;
 
    mypi = 1.0;
    third = 0.3333333333333333333333;
 
    if (p0[1] == p1[1]) {
        printf("ERROR: y0 == y1 in cal_vf_cyl_cord\n");
        return;
    }
    for (i = 0; i < 2; i++) {
        r0[i] = p0[i] - ctr[i];
        r1[i] = p1[i] - ctr[i];
    }
    slope = (r1[0] - r0[0])/(r1[1] - r0[1]);
 
    s2 = slope * slope;
    x2 = r0[0] * r0[0];
    y2 = r0[1] * r0[1];
    xy = r0[0] * r0[1];
    r2 = radius * radius;
    c0 = r2 - x2 + 2.0 * slope * xy - s2 * y2;
    c1 = - 2.0 * r0[0] * slope + s2 * 2.0 * r0[1];
    c2 = - (1.0 + s2);
 
    tmp  = r1[1] * r1[1];
    *vol = third * c2 *(tmp * r1[1] - y2 * r0[1])
         + 0.5 * c1 *(tmp - y2) + c0 *(r1[1] - r0[1]);
 
    if (r1[1] < r0[1]) {
        *vol = -(*vol);
    }
    assert(*vol >= 0.0);
 
    *vol *= mypi;
 
    return;
 }
 
 
void sph_poly2d_new(int ifinquiry, int geop, int dim, double *ctr, double r,
                double *ptrs, int nnode,
                int *mixed, double *vol, double *ctroid)
{
//   ctroid[2] is the centroid of the overlap if there is an overlap.
 
     int   intersected;
     int   allinside, alloutside, ifescape, already, included[100];
     int   n0, i0, i1, i, k, offset, ncr_tot, ncr, i_inside, nint, szdim;
     int   nn_poly, nn_cir;
     int   isinside[100], ncross[100], i0_cross[100], i1_cross[100];
     int   first_cross[100], node_exact[100];
     double third, two3rd;
     double r2, r3, ds2, tmp, a, b, c, delt2, rt, t0, t1;
     double sina, vol_poly, vol_cord;
     double myctr[2], dx[2], delx[2], dp[2];
     double cross[100], ptrs_res[200], ptrs_cir[200];
     double *p0, *p1, *pt, *ptc;
     double dv, r0, z0, r1, z1, dz, rc, zc, r01, r02, r12, r0r1;
     double area_pie, area_tri, ha, sinha, cosha, xcarea, dxc_pie;
     double factor, pm, rc_pie, rc_tri;
 
     assert(nnode < 100);
     szdim = dim * sizeof(double);
 
     third  = 0.33333333333333333333333;
     two3rd = third + third;
 
     r2 = r * r;
     r3 = r * r2;
 
     *vol = 0.0;
 
//   ctroid will be the centroid of the overlap if there is a overlap.
     memcpy(ctroid, ctr, (size_t)szdim);
 
//   check whether all vertices are in inside the cricle.
 
     offset = 0;
     for (i = 0; i < nnode; i++) {
         p0 = ptrs + offset;
         ds2 = 0.0;
         for (k = 0; k < dim; k++) {
             tmp = p0[k] - ctr[k];
             ds2 += (tmp * tmp);
         }
         if (ds2 < r2) {
             isinside[i] = 1;
         }
         else if (ds2 > r2) {
             isinside[i] = 0;
         }
         else {  // touch
             isinside[i] = -1;
         }
         offset += dim;
     }
     allinside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *mixed = 0;
         *vol   = 0.0;
         if (geop != 1) {
             printf("ERROR: geop in sph_poly2d");
             exit(1);
         }
         for (i = 0; i < nnode; i++) {
             i1 = (i + 1) % nnode;
             p0 = ptrs + (dim * i);
             p1 = ptrs + (dim * i1);
             dv = p0[0] * p1[1] - p1[0] * p0[1];
             *vol += dv;
         }
         *vol *= 0.5;
         if (*vol < 0.0) *vol = - (*vol);
 
         return;
     }
     alloutside = 1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i]) {
             alloutside = 0;
             break;
         }
     }
     if (alloutside) {
//       check whether the circle intersects with any edge
         intersected = 0;
         for (i = 0; i < nnode; i++) {
             i1 = (i + 1) % nnode;
             p0 = ptrs + (dim * i);
             p1 = ptrs + (dim * i1);
 
             for (k = 0; k < dim; k++) {
                 dx[k]   = p1[k] - p0[k];
                 delx[k] = p0[k] - ctr[k];
             }
             a = 0.0;
             b = 0.0;
             c = 0.0;
             for (k = 0; k < dim; k++) {
                 a += (dx[k]   * dx[k]);
                 b += (dx[k]   * delx[k]);
                 c += (delx[k] * delx[k]);
             }
             b  = b + b;
             c -= r2;
             delt2 = b * b - 4.0 * a * c;
             if ((delt2 < 0.0) || (a == 0.0)) continue;
 
             rt = sqrt(delt2);
             tmp = 0.5 / a;
             t0  = (-b - rt) * tmp;
             t1  = (-b + rt) * tmp;
 
             if (((t0 >= 0.0) && (t0 <= 1.0)) || ((t1 >= 0.0) && (t1 <= 1.0))) {
                 intersected = 1;
                 break;
             }
         }
         if (!intersected) {
             *mixed = 0;
             *vol = 0.0;
             return;
         }
     }
//   find intersection between circle and line segement
 
     ncr_tot = 0;
     for (i = 0; i < nnode; i++) {
         i1 = (i + 1) % nnode;
 
//          if ((isinside[i] == 1) && (isinside[i1] == 1)) continue;
 
         p0 = ptrs + (dim * i);
         p1 = ptrs + (dim * i1);
 
         for (k = 0; k < dim; k++) {
             dx[k]   = p1[k] - p0[k];
             delx[k] = p0[k] - ctr[k];
         }
         a = 0.0;
         b = 0.0;
         c = 0.0;
         for (k = 0; k < dim; k++) {
             a += (dx[k] * dx[k]);
             b += (dx[k] * delx[k]);
             c += (delx[k] * delx[k]);
         }
         b  = b + b;
         c -= r2;
         delt2 = b * b - 4.0 * a * c;
 
         first_cross[i] = -1;
         if ((delt2 < 0.0) || (a == 0.0)) {
             ncross[i] = 0;
             continue;
         }
         rt = sqrt(delt2);
         tmp = 0.5 / a;
         t0  = (-b - rt) * tmp;
         t1  = (-b + rt) * tmp;
 
         if (fabs(t0 - t1) < small) {
             ncross[i] = 0;
             continue;
         }
         ncross[i] = 0;
         if ((t0 >= 0.0) && (t0 <= 1.0)) {
             ncross[i]++;
             cross[ncr_tot]    = t0;
             i0_cross[ncr_tot] = i;
             i1_cross[ncr_tot] = i1;
             first_cross[i]    = ncr_tot;
             if (t0 == 0.0) {
                 node_exact[ncr_tot] = i;
             }
             else if (t0 == 1.0) {
                 node_exact[ncr_tot] = i1;
             }
             else {
                 node_exact[ncr_tot] = -1;
             }
             ncr_tot++;
         }
         if ((t1 >= 0.0) && (t1 <= 1.0)) {
             if (ncross[i] == 0) {
                 first_cross[i] = ncr_tot;
             }
             ncross[i]++;
             cross[ncr_tot]    = t1;
             i0_cross[ncr_tot] = i;
             i1_cross[ncr_tot] = i1;
             if (t1 == 0.0) {
                 node_exact[ncr_tot] = i;
             }
             else if (t1 == 1.0) {
                 node_exact[ncr_tot] = i1;
             }
             else {
                 node_exact[ncr_tot] = -1;
             }
             ncr_tot++;
         }
     }
     if (ncr_tot < 2) {
         *mixed = 0;
         return;
     }
     if (ifinquiry) {
         *mixed = 1;
         return;
     }
     i_inside = -1;
     for (i = 0; i < nnode; i++) {
         if (isinside[i] == 1) {
             i_inside = i;
             break;
         }
     }
     for (i = 0; i < nnode; i++) {
         included[i] = 0;
     }
     ctroid[0] = 0.0;
     ctroid[1] = 0.0;
 
     *mixed = 1;
     if (i_inside >= 0) {
         p0 = ptrs + (i_inside * dim);
         memcpy(ptrs_res, p0, (size_t)szdim);
         included[i_inside] = 1;
         nn_poly = 1;
         pt     = ptrs_res + dim;
         ptc    = ptrs_cir;
         nn_cir = 0.0;
 
         for (i = 0; i < nnode; i++) {
             i0 = (i_inside + i) % nnode;
             i1 = (i0 + 1) % nnode;
             if ((ncross[i0] == 0) && (isinside[i1] == 1)) {
                 if (!included[i1]) {
                     p1 = ptrs + (i1 * dim);
                     memcpy(pt, p1, (size_t)szdim);
                     nn_poly++;
                     pt += dim;
                 }
                 included[i1] = 1;
             }
             else if (ncross[i0] == 1) {
                 ncr = first_cross[i0];
                 t0  = cross[ncr];
                 p0  = ptrs + (i0 * dim);
                 p1  = ptrs + (i1 * dim);
                 n0  = node_exact[ncr];
                 already = 0;
                 if (n0 >= 0) {
                     already = included[n0];
                 }
                 if (!already)  {
                     pt[0] = p0[0] + (p1[0] - p0[0]) * t0;
                     pt[1] = p0[1] + (p1[1] - p0[1]) * t0;
                     memcpy(ptc, pt, (size_t)szdim);
                     nn_poly++;
                     nn_cir++;
                     ptc += dim;
                     pt  += dim;
                     if (n0 >= 0) {
                         included[n0] = 1;
                     }
                 }
                 if ((isinside[i1] == 1) && !included[i1]) {
                     memcpy(pt, p1, (size_t)szdim);
                     included[i1] = 1;
                     nn_poly++;
                     pt += dim;
                 }
             }
             else if (ncross[i0] == 2) {
                 ncr = first_cross[i0];
                 t0  = cross[ncr];
                 p0  = ptrs + (i0 * dim);
                 p1  = ptrs + (i1 * dim);
 
                 n0  = node_exact[ncr];
                 already = 0;
                 if (n0 >= 0) {
                     already = included[n0];
                 }
                 if (!already)  {
                     pt[0] = p0[0] + (p1[0] - p0[0]) * t0;
                     pt[1] = p0[1] + (p1[1] - p0[1]) * t0;
                     memcpy(ptc, pt, (size_t)szdim);
                     nn_poly++;
                     nn_cir++;
                     ptc += dim;
                     pt  += dim;
                     if (n0 >= 0) {
                         included[n0] = 1;
                     }
                 }
                 t0 = cross[ncr+1];
                 n0 = node_exact[ncr+1];
                 already = 0;
                 if (n0 >= 0) {
                     already = included[n0];
                 }
                 if (!already)  {
                     pt[0] = p0[0] + (p1[0] - p0[0]) * t0;
                     pt[1] = p0[1] + (p1[1] - p0[1]) * t0;
                     memcpy(ptc, pt, (size_t)szdim);
                     nn_poly++;
                     nn_cir++;
                     ptc += dim;
                     pt  += dim;
                     if (n0 >= 0) {
                         included[n0] = 1;
                     }
                 }
             }
         }
         nint = (pt - ptrs_res)/dim;
         vol_poly = 0.0;
 
         if (geop == 1) {
             myctr[0] = 0.0;
             myctr[1] = 0.0;
             for (i = 0; i < nn_poly; i++) {
                 i1 = (i + 1) % nn_poly;
                 p0 = ptrs_res + (dim * i);
                 p1 = ptrs_res + (dim * i1);
                 dv = p0[0] * p1[1] - p1[0] * p0[1];
                 myctr[0] += ((p0[0] + p1[0]) * dv);
                 myctr[1] += ((p0[1] + p1[1]) * dv);
                 vol_poly += dv;
             }
             vol_poly *= 0.5;
             if (vol_poly != 0.0) {
                 factor = 1.0/(6.0 * vol_poly);
                 myctr[0] *= factor;
                 myctr[1] *= factor;
             }
             if (vol_poly < 0.0) vol_poly = - vol_poly;
             for (i = 0; i < dim; i++) {
                 ctroid[i] += (vol_poly * myctr[i]);
             }
         }
         else if (geop == 2) {
             for (i = 0; i < nn_poly; i++) {
                 i1 = (i + 1) % nn_poly;
                 p0 = ptrs_res + (dim * i);
                 p1 = ptrs_res + (dim * i1);
                 r0 = p0[0];
                 z0 = p0[1];
                 r1 = p1[0];
                 z1 = p1[1];
 
                 dz   = z1 - z0;
                 r02  = r0 * r0;
                 r12  = r1 * r1;
                 r0r1 = r0 * r1;
 
//               dv = 2 * integral (r dr dz)
 
                 dv = dz *(r12 + r0r1 + r02);
                 vol_poly += dv;
 
                 rc = dz *(r02 * r0 + r02 * r1 + r0 * r12 + r12 * r1);
                 myctr[0] += rc;
 
                 r01 = r0 + r1;
                 zc = dz *(r01 * r01 *(z0 + z1) + 2.0 *(r02 * z0 + r12 * z1));
                 myctr[1] += zc;
             }
             for (i = 0; i < dim; i++) {
                 ctroid[i] += (vol_poly * myctr[i]);
             }
         }
     }
     else {
         vol_poly = 0.0;
     }
//   calculate the volume of circular sections
 
     nint = (ptc - ptrs_cir)/dim;
 
     ifescape = 0;
     *vol = 0.0;
     for (i = 0; i < nn_cir-1; i++) {
 
         if (!ifescape) {
             i1 = i + 1;
             p0 = ptrs_cir + (dim * i);
             p1 = ptrs_cir + (dim * i1);
             if ((p0[0] == p1[0]) && (p0[1] == p1[1])) continue;
 
             if (geop == 1) {
                 ds2 = 0.0;
                 for (k = 0; k < dim; k++) {
                     tmp  = p1[k] - p0[k];
                     ds2 += (tmp * tmp);
                 }
                 c = sqrt(ds2);
                 a = 2.0 * asin(0.5 * c / r);
                 sina = sin(a);
                 area_pie = 0.5 * r2 * a;
                 area_tri = 0.5 * r2 * sina;
                 vol_cord = area_pie - area_tri;
 
                 *vol += vol_cord;
 
//               center of the "pie"
 
                 ha      = 0.5 * a;
                 sinha   = sin(ha);
                 cosha   = sqrt(MAX(1.0 - sinha * sinha, 0.0));
                 xcarea  = two3rd * r3 * sinha;
                 dxc_pie = xcarea / area_pie;
                 factor  = dxc_pie / (r * cosha);
 
                 for (k = 0; k < dim; k++) {
                     pm = 0.5 *(p0[k] + p1[k]);
                     rc_pie = ctr[k] + factor * (pm - ctr[k]);
                     rc_pie *= area_pie;
 
//                   center of the "triangle", (ctr, p0, p1)
 
                     rc_tri = third *(ctr[k] + p0[k] + p1[k]);
                     rc_tri *= area_tri;
 
//                   add the product of centroid and area for the "cord"
 
                     ctroid[k] += (rc_pie - rc_tri);
                 }
             }
             else if (geop == 2) {
                 printf("ERROR: geop = 2 not done yet for sph_poly2d.\n");
             }
         }
         ifescape = 1 - ifescape;
     }
     *vol += vol_poly;
 
     if (*vol  != 0.0) {
         factor = 1.0/(*vol);
         for (k = 0; k < dim; k++) {
             ctroid[k] *= factor;
         }
     }
     return;
 }
 
 
void sph_rec(int geop, int dim, double *ctr, double r, double *xl, double *dx, double vcell, double *vol)
{
/**
     Assume that x-plane x = ctr[0], or y-plane, y = cr0[1], or z plane z = ctr[2],
     does't cut the rectangular.
 **/
     int    i;
     double xc[3], cl[3], cr[3], pcl[3], pcr[3];
 
     if (dim == 2) { 
         sph_rec2d_ctr(geop, dim, r, ctr, xl, dx, vol); 
     } 
     else if (dim == 3) {  
         for (i = 0; i < dim; i++) {
             cl[i] = xl[i] - ctr[i];
             cr[i] = cl[i] + dx[i];
             xc[i] = 0.5 *(cl[i] + cr[i]);
         }
         for (i = 0; i < dim; i++) {
             if (xc[i] < 0.0) {
                 pcl[i] = - cr[i];
             }
             else {
                 pcl[i] =   cl[i];
             }
         }
         sph_rec3d(dim, r, pcl, dx, vcell, vol);
     }
     return;
 }
 
void sph_rec2d_ctr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/***
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0).
 
     geop:  = 1   for Cartesian
            = 2   for Cylindrical
            = 3   for Spherical, not yet
 
     dim        : input
     r          : input, radius of the sphere
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assume all vertices of the rectangular,
     x0[0] >= 0. x0[1] >= 0
****/
     int i; 
     double small, r2, x1[2], cl[2], cr[2], r2mn, r2mx, r2x, r2y;
 
     small  = 1.0e-08 * dx[0];

     r2 = r * r;

     for (i = 0; i < dim; i++) { 
         x1[i] = x0[i] + dx[i];
         cl[i] = x0[i] - ctr[i];
         cr[i] = cl[i] + dx[i];
     } 
     if ((ctr[0] <= x0[0] + small) && (ctr[1] <= x0[1] + small)) { 
//       ctr is at the lower left of the rectangular, ll  
         r2mn = cl[0] * cl[0] + cl[1] * cl[1];
         r2mx = cr[0] * cr[0] + cr[1] * cr[1];
         r2x  = cr[0] * cr[0] + cl[1] * cl[1]; 
         r2y  = cr[1] * cr[1] + cl[0] * cl[0]; 
     } 
     else if ((ctr[0] + small >= x1[0]) && (ctr[1] + small >= x1[1])) { 
//       ctr is at the upper right of the rectangular, ur  
         r2mx = cl[0] * cl[0] + cl[1] * cl[1]; 
         r2mn = cr[0] * cr[0] + cr[1] * cr[1]; 
         r2x  = cl[0] * cl[0] + cr[1] * cr[1];
         r2y  = cr[0] * cr[0] + cl[1] * cl[1];  
     }
     else if ((ctr[0] <= x0[0] + small) && (ctr[1] + small >= x1[1])) {
//       ctr is at the upper left of the reactangular, ul   
         r2mn = cl[0] * cl[0] + cr[1] * cr[1];
         r2mx = cr[0] * cr[0] + cl[1] * cl[1]; 
         r2x  = cr[0] * cr[0] + cr[1] * cr[1];
         r2y  = cl[0] * cl[0] + cl[1] * cl[1];
     }
     else if ((ctr[0] + small >= x1[0]) && (ctr[1] <= x0[1] + small)) {  
//       ctr is at the lower right of the rectangular, lr  
         r2mn = cr[0] * cr[0] + cl[1] * cl[1]; 
         r2mx = cl[0] * cl[0] + cr[1] * cr[1]; 
         r2x  = cl[0] * cl[0] + cl[1] * cl[1];
         r2y  = cr[0] * cr[0] + cr[1] * cr[1];  
     }  
     else { 
         printf("ERROR: a rectangle is cut by a plane through the ctr of sphere.\n");
         assert(0); 
     } 
     if (r2mn >= r2) {
         *vol = 0.0;
     }
     else if (r2mx <= r2) { 
         *vol = dx[0] * dx[1]; 
         if (geop == 2) { 
             *vol *= (x0[0] + x0[0] + dx[0]); 
         }
     } 
     else if ((ctr[0] <= x0[0] + small ) && (ctr[1] <= x0[1] + small)) { 
//       The center of sphere is at the lower left of the rectangular, ll   

         if ((r2 <= r2x) && (r2 <= r2y)) { 
             sph_rec2d0_ll(geop, dim, r, ctr, x0, dx, vol); 
         } 
         else if ((r2 >= r2x) && (r2 >= r2y)) {
             sph_rec2dxy_ll(geop, dim, r, ctr, x0, dx, vol);
         } 
         else if ((r2 > r2x) && (r2 <= r2y)) {
             sph_rec2dx_ll(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 <= r2x) && (r2 > r2y)) {
             sph_rec2dy_ll(geop, dim, r, ctr, x0, dx, vol);
         }
     }  
     else if ((ctr[0] + small >= x1[0]) && (ctr[1] + small >= x1[1])) { 
//       The center of sphere is at the upper right of the rectangular, ur  
         if ((r2 <= r2x) && (r2 <= r2y)) {
             sph_rec2d0_ur(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 >= r2x) && (r2 >= r2y)) {
             sph_rec2dxy_ur(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 > r2x) && (r2 <= r2y)) {
             sph_rec2dx_ur(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 <= r2x) && (r2 > r2y)) {
             sph_rec2dy_ur(geop, dim, r, ctr, x0, dx, vol);
         }
     }
     else if ((ctr[0] <= x0[0] + small) && (ctr[1] + small >= x1[1])) { 

//       The center of sphere is at the upper left of the rectangular, ul  
         if ((r2 <= r2x) && (r2 <= r2y)) {
             sph_rec2d0_ul(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 >= r2x) && (r2 >= r2y)) {
             sph_rec2dxy_ul(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 > r2x) && (r2 <= r2y)) {
             sph_rec2dx_ul(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 <= r2x) && (r2 > r2y)) {
             sph_rec2dy_ul(geop, dim, r, ctr, x0, dx, vol);
         }
     }
     else if ((ctr[0] + small >= x1[0]) && (ctr[1] <= x0[1] + small)) {
          
//       The center of sphere is at the lower right of the rectangular, lr 
         if ((r2 <= r2x) && (r2 <= r2y)) {
             sph_rec2d0_lr(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 >= r2x) && (r2 >= r2y)) {
             sph_rec2dxy_lr(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 > r2x) && (r2 <= r2y)) {
             sph_rec2dx_lr(geop, dim, r, ctr, x0, dx, vol);
         }
         else if ((r2 <= r2x) && (r2 > r2y)) {
             sph_rec2dy_lr(geop, dim, r, ctr, x0, dx, vol);
         }
     } 
     return;
 }

double integral_l(double r, double *ctr, double *x0, double *dx, double y0, double y1)
{   
/***
                -------------
                |           |
      .         |           |
           .    |           |
            y1 .| 1         |
                | .         |
            y0  ----.--------
                   2  .
or

                        o 
            y1  -------------
                |     o     |
                |    o      |
                |  o        |
            y0  o           |
                |           |
           o    ----.--------

*****/
//   vol = integral from y0 to y1  pi * (x^2 - x0[0]^2) dy 

     int niter_mx, npart_mx, niter, npart, hpart, i;
     double sixth;
     double r2, x02, dy, vol_inv, res_previous, res, ddy, y, yyc, tmp, x, err;
     double fp[2049], fm[2048]; 

//     double third;
//     double xc2, y0yc, y1yc, a, tmp1, tmp2, tmp3, rt, v1, v0;

     sixth = 0.1666666666666666666666667;
     dy = y1 - y0;
     r2  = r * r;
     x02 = x0[0] * x0[0];

//     third = 0.33333333333333333333333;
//     xc2 = ctr[0] * ctr[0];
//     tmp1 = (xc2 - x02 + r2) * dy;
//     y0yc = y0 - ctr[1];
//     y1yc = y1 - ctr[1];
//     tmp2 = third *(y1yc * y1yc + y0yc * y0yc + y0yc * y1yc) * dy;
// 
//     if (ctr[0] == 0.0) { 
//         tmp3 = 0.0;
//     }
//     else { 
//         tmp = r2 - y1yc * y1yc;
//         tmp = MAX(tmp, 0.0);
//         rt  = sqrt(tmp);
//
//         tmp = r /rt;
//         a   = atan(tmp);
//         tmp = r2 * a;
//         rt *= y1yc;
//         v1  = rt + tmp; 
//
//         tmp = r2 - y0yc * y0yc;
//         tmp = MAX(tmp, 0.0);
//         rt  = sqrt(tmp);
//
//         tmp = r /rt; 
//         a   = atan(tmp);
//         tmp = r2 * a; 
//         rt *= y0yc;  
//         v0  = rt + tmp;
// 
//         tmp3 = ctr[0] *(v1 - v0);  
//     } 
//     res = tmp1 - tmp2 + tmp3;

//   Simpson Integral   intgeral_a^b f(x) dx = (h/3) [f(a) + 4 f((a+b)/2) + f(b)]; h = (b-a)/2.

     niter_mx = 10;
     npart_mx = 2048; 
      
     niter = 0;
     npart = 1;
     ddy   = dy;

     y     = y0;  
     yyc   = y - ctr[1]; 
     tmp   = MAX(r2 - yyc * yyc, 0.0);
     x     = ctr[0] + sqrt(tmp);   
     fp[0] = x * x - x02; 

     y     = y1;
     yyc   = y - ctr[1];
     tmp   = MAX(r2 - yyc * yyc, 0.0);
     x     = ctr[0] + sqrt(tmp);
     fp[1] = x * x - x02; 

     err = 1.0;
     vol_inv = dx[0] * dx[1] * (x0[0] + x0[0] + dx[0]);
     res_previous = vol_inv;
     vol_inv = 1.0/vol_inv;

     while ((err > tol) && (niter < niter_mx)) { 

           hpart = npart/2;
           for (i = hpart; i > 0; i--) { 
               fp[i+i] = fp[i];     
           } 
           for (i = 0; i < hpart; i++) {
               fp[i+i+1] = fm[i];
           }
           res = 0.0;
           y = y0 + 0.5 * ddy;
           for (i = 0; i < npart; i++) {  
               yyc   = y - ctr[1];
               tmp   = MAX(r2 - yyc * yyc, 0.0); 
               x     = ctr[0] + sqrt(tmp); 
               fm[i] = x * x - x02;
               res += (4.0 * fm[i]);

               res += (fp[i] + fp[i+1]); 

               y += ddy;
           }  
           res *= (sixth * ddy);

           err = fabs(res - res_previous) * vol_inv; 
           res_previous   = res;

           niter++;
           ddy *= 0.5;
           npart = npart + npart;
     } 

     return res; 
   } 


double integral_r(double r, double *ctr, double *x0, double *dx, double y1, double y2)
{   
/***
                                        
                    o    
                      o  
             y2 -------2-----
                |        o  |
                |          o| 1  y1  
                |           | o  
                |           |    o 
                |           |
                -------------
                           x1  
or 

 
                -------------        o  
                |           |   o  
                |           o  y2 
                |       o   |  
                |           |   
                |    o      |    
                ---o---------  y1 
                   2  
                 o
**/
//     vol = integral from y1 to y2  pi * (x1^2 - x^2) dy 

       int niter_mx, npart_mx, niter, npart, hpart, i;
       double sixth;
       double r2, x1, x12, dy, vol_inv, res_previous, res, ddy, y, yyc, tmp, x, err;
       double fp[2049], fm[2048];

       sixth = 0.1666666666666666666666667;
       dy = y2 - y1;
       r2  = r * r;
       x1  = x0[0] + dx[0];
       x12 = x1 * x1;

     niter_mx = 10;
     npart_mx = 2048;

     niter = 0;
     npart = 1;
     ddy   = dy;

     y     = y1;
     yyc   = y - ctr[1];
     tmp   = MAX(r2 - yyc * yyc, 0.0);
     x     = ctr[0] - sqrt(tmp);
     fp[0] = x12 - x * x;

     y     = y2;
     yyc   = y - ctr[1];
     tmp   = MAX(r2 - yyc * yyc, 0.0);
     x     = ctr[0] - sqrt(tmp);
     fp[1] = x12 - x * x;

     err = 1.0;
     vol_inv = dx[0] * dx[1] * (x0[0] + x0[0] + dx[0]);
     res_previous = vol_inv;
     vol_inv = 1.0/vol_inv;

     while ((err > tol) && (niter < niter_mx)) {

           hpart = npart/2;
           for (i = hpart; i > 0; i--) {
               fp[i+i] = fp[i];
           }
           for (i = 0; i < hpart; i++) {
               fp[i+i+1] = fm[i];
           }
           res = 0.0;
           y = y1 + 0.5 * ddy;
           for (i = 0; i < npart; i++) {
               yyc   = y - ctr[1];
               tmp   = MAX(r2 - yyc * yyc, 0.0);
               x     = ctr[0] - sqrt(tmp);
               fm[i] = x12 - x * x;
               res += (4.0 * fm[i]);

               res += (fp[i] + fp[i+1]);

               y += ddy;
           }
           res *= (sixth * ddy);

           err = fabs(res - res_previous) * vol_inv;
           res_previous   = res;

           niter++;
           ddy *= 0.5;
           npart = npart + npart;
     }
     return res; 
 } 


void sph_rec2d0_ll(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only one vertex, x0, of the the rectangular.
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : ctr of sphere  
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     The center, ctr, is assumed at the lower left of the rectangular. 

                -------------
                |           |
      .         |           |
           .    |           |
               .| 1         |
                | .         |
                ----.--------
                   2  .
                        .
                         .
**/
     double third, r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
 
     third = 0.33333333333333333333333;
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] and the circle
 
     xi1 = x0[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] + sqrt(y2);
 
     if (geop == 1) {
 
//       intersection of the line y = x0[1] and the circle
 
         yi2 = x0[1];
         y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
         x2  = MAX(r2 - y2, 0.0);
         xi2 = ctr[0] + sqrt(x2);
 
         x2  = xi2 - xi1;
         y2  = yi1 - yi2;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   vol of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi2 - x0[0])*(yi1 - x0[1]);
         *vol = volt + vols;
     }
     else if (geop == 2) {

//       vol = integral from y0 to yi1  pi * (x^2 - x0[0]^2) dy 

         *vol = integral_l(r, ctr, x0, dx, x0[1], yi1);
     }
     return;
 }
 
void sph_rec2dx_ll(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim} : input, the center of the sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assumed that the center is at the lower left of the rectangular,
 
      o         -------------
           o    |           |
                o           |
                |    o      |
                |        o  |
                |           i 
                |           |  o
                -------------     o
 
**/
     double r2, xi, x2, y2, yi, volr, x00[2], dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] + dx[0] and the circle
 
     xi = x0[0] + dx[0];
     x2 = (xi - ctr[0]) * (xi - ctr[0]);
     y2 = MAX(r2 - x2, 0.0);
     yi = ctr[1] + sqrt(y2);
     volr = (yi - x0[1]) * dx[0];
 
     if (geop == 2) {
         volr *= (xi + x0[0]);
     }
     x00[0]  = x0[0];
     x00[1]  = yi;
     dx00[0] = dx[0];
     dx00[1] = MAX(0.0, x0[1] + dx[1] - yi);
 
     sph_rec2d0_ll(geop, dim, r, ctr, x00, dx00,  vol);
 
     *vol += volr;
 
     return;
 }
 
void sph_rec2dy_ll(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0], x0[1] + dx[1])
 
     dim        : input
     r          : input, radius of the sphere
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assumed that the center ctr is at the lower left of the rectangular,
     x0[0] >= 0. x0[1] >= 0, x0[2] > 0
 
          o
                o  i
                ---o-------------
                |     o         |
                |        o      |
                |          o    |
                |               |
                |            o  |
                --------------o--
                               o
**/
     double r2, xi, x2, y2, yi, volr, x00[2], dx00[2];
 
     r2 = r * r;
 
//   intersection of the line y = x0[1] + dx[1] and the circle
 
     yi = x0[1] + dx[1];
     y2 = (yi - ctr[1]) * (yi - ctr[1]);
     x2 = MAX(r2 - y2, 0.0);
     xi = ctr[0] + sqrt(x2);
     volr = (xi - x0[0]) * dx[1];
 
     if (geop == 2) {
         volr *= (xi + x0[0]);
     }
     x00[0]  = xi;
     x00[1]  = x0[1];
     dx00[0] = MAX(0.0, x0[0] + dx[0] - xi);
     dx00[1] = dx[1];
 
     sph_rec2d0_ll(geop, dim, r, ctr, x00, dx00,  vol);
 
     *vol += volr;
 
     return;
 }
 
void sph_rec2dxy_ll(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : input, the center of sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assumed that the center is at the lower left of the rectangular,
 
                o
                    o 2
                -------o-----
                |         o |
                |           o 1
                |           | o
                |           |  o
                |           |   o
                -------------
 
**/
     double r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
     double x00[2], dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] + dx[0] and the circle
 
     xi1 = x0[0] + dx[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] + sqrt(y2);
 
//   intersection of the line y = x0[1] + dx[1] and the circle
 
     yi2 = x0[1] + dx[1];
     y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
     x2  = MAX(r2 - y2, 0.0);
     xi2 = ctr[0] + sqrt(x2);
 
     if (geop == 1) {
         x2  = xi1 - xi2;
         y2  = yi2 - yi1;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   area of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi1 - xi2) * (yi2 - yi1);
 
         *vol = dx[0] * dx[1] - volt + vols;
     }
     else if (geop == 2) {

         vols = dx[0] *(yi1 - x0[1]) *(xi1 + x0[0]);
         volt = (yi2 - yi1) * (xi2 - x0[0]) * (xi2 + x0[0]);
 
         x00[0]  = xi2;
         x00[1]  = yi1;
         dx00[0] = MAX(0.0, xi1 - xi2);
         dx00[1] = MAX(0.0, yi2 - yi1);
         sph_rec2d0_ll(geop, dim, r, ctr, x00, dx00, vol);
 
         *vol += (vols + volt);
     }
     return;
 }

void sph_rec2d0_ul(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only one vertex, x0, of the the rectangular.
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : ctr of sphere  
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assumed that the center is the upper left of the rectangular,
 
                         o
                        o 
                -------2-----
                |     o     |
                |   o       |
                | o         |
              1 o           |
              o |           |
           o    ----.--------
                          
**/
     double third, r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
 
     third = 0.33333333333333333333333;
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] and the circle
 
     xi1 = x0[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] - sqrt(y2);
 
     if (geop == 1) {
 
//       intersection of the line y = x0[1] + dx[1] and the circle
 
         yi2 = x0[1] + dx[1];
         y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
         x2  = MAX(r2 - y2, 0.0);
         xi2 = ctr[0] + sqrt(x2);
 
         x2  = xi2 - xi1;
         y2  = yi2 - yi1;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   vol of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi2 - xi1)*(yi2 - yi1);
         *vol = volt + vols;
     }
     else if (geop == 2) {

//       vol = integral from yi1 to yi2  pi * (x^2 - x0[0]^2) dy 

         *vol = integral_l(r, ctr, x0, dx, yi1, x0[1] + dx[1]); 
     }
     return;
 }
 
void sph_rec2dx_ul(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim} : input, the center of the sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 


                                  o 
                -------------   o 
                |           | o   
                |           o i  
                |        o  |
                |     o     |
                |  o        | 
              o |           |  
          o     -------------     

     It is assumed that the center is the upper left of the rectangular,
**/
     double r2, xi, x2, y2, yi, volr, dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] + dx[0] and the circle
 
     xi = x0[0] + dx[0];
     x2 = (xi - ctr[0]) * (xi - ctr[0]);
     y2 = MAX(r2 - x2, 0.0);
     yi = ctr[1] - sqrt(y2);
     volr = (x0[1] + dx[1] - yi) * dx[0];
 
     if (geop == 2) {
         volr *= (xi + x0[0]);
     }
     dx00[0] = dx[0];
     dx00[1] = MAX(0.0, yi - x0[1]);
 
     sph_rec2d0_ul(geop, dim, r, ctr, x0, dx00,  vol);
 
     *vol += volr;
 
     return;
 }
 

void sph_rec2dy_ul(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0], x0[1] + dx[1])
 
     dim        : input
     r          : input, radius of the sphere
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assumed that the ctr is at the upper left of the rectangle
 
                               o .
                .             o  
                -------------o---
                |           o   |
                |          o    |  
                |         o     |
                |       o       |
                |    o          |
                ---o-------------
                 o  i 
              o
**/
     double r2, xi, x2, y2, yi, volr, x00[2], dx00[2];
 
     r2 = r * r;
 
//   intersection of the line y = x0[1] and the circle
 
     yi = x0[1];
     y2 = (yi - ctr[1]) * (yi - ctr[1]);
     x2 = MAX(r2 - y2, 0.0);
     xi = ctr[0] + sqrt(x2);

     volr = (xi - x0[0]) * dx[1];
 
     if (geop == 2) {
         volr *= (xi + x0[0]);
     }
     x00[0]  = xi;
     x00[1]  = x0[1];
     dx00[0] = MAX(0.0, x0[0] + dx[0] - xi);
     dx00[1] = dx[1];
 
     sph_rec2d0_ul(geop, dim, r, ctr, x00, dx00,  vol);
 
     *vol += volr;
 
     return;
 }
 
void sph_rec2dxy_ul(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : input, the center of sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assumed that ctr is at the upper left of the rectangular
 
                -------------  o       
                |           |      
                |           o  1
                |         o |  
                |           |   
                |      o    |    
                ---o---------
                o  2  
            o
**/
     double r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
     double x00[2], dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] + dx[0] and the circle
 
     xi1 = x0[0] + dx[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] - sqrt(y2);
 
//   intersection of the line y = x0[1] + dx[1] and the circle
 
     yi2 = x0[1];
     y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
     x2  = MAX(r2 - y2, 0.0);
     xi2 = ctr[0] + sqrt(x2);
 
     if (geop == 1) {
         x2  = xi1 - xi2;
         y2  = yi1 - yi2;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   area of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi1 - xi2) * (yi1 - yi2);
 
         *vol = dx[0] * dx[1] - volt + vols;
     }
     else if (geop == 2) {
         vols = dx[0] * (x0[1] + dx[1] - yi1) * (x0[0] + xi1); 
         volt = (yi1 - x0[1]) * (xi2 - x0[0]) * (xi2 + x0[0]);
 
         x00[0]  = xi2;
         x00[1]  = x0[1];
         dx00[0] = MAX(0.0, xi1 - xi2);
         dx00[1] = MAX(0.0, yi1 - yi2);

         sph_rec2d0_ul(geop, dim, r, ctr, x00, dx00, vol);
 
         *vol += (vols + volt);
     }
     return;
 }

void sph_rec2d0_ur(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only one vertex, x0, of the the rectangular.
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : ctr of sphere  
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     The ctr is assumed to be at the upper right of the rectangular 
 
                    o    
                      o  
                -------2-----
                |        o  |
                |          o| 1 
                |           | o  
                |           |    o 
                |           |
                ----.--------
**/
     double third, r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
 
     third = 0.33333333333333333333333;
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] and the circle
 
     xi1 = x0[0] + dx[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] - sqrt(y2);
 
     if (geop == 1) {
 
//       intersection of the line y = x0[1] + dx[1] and the circle
 
         yi2 = x0[1] + dx[1];
         y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
         x2  = MAX(r2 - y2, 0.0);
         xi2 = ctr[0] - sqrt(x2);
 
         x2  = xi1 - xi2;
         y2  = yi2 - yi1;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   vol of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi1 - xi2) * (yi2 - yi1);
         *vol = volt + vols;
     }
     else if (geop == 2) {

//       vol = integral from yi1 to yi2  pi * (x^2 - x0[0]^2) dy 

         *vol = integral_r(r, ctr, x0, dx, yi1, x0[1] + dx[1]); 
     }
     return;
 }
 
void sph_rec2dx_ur(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim} : input, the center of the sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.

     The ctr is assumed to be at the upper right of the rectangular 

             o                        
                -------------      
                |           |     
                o i         |    
                |           |
                |     o     |  
                |           o      
                |           |      o 
                -------------     
              x0
 
**/
     double r2, xi, x2, y2, yi, volr, dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] and the circle
 
     xi = x0[0];
     x2 = (xi - ctr[0]) * (xi - ctr[0]);
     y2 = MAX(r2 - x2, 0.0);
     yi = ctr[1] - sqrt(y2);
     volr = (x0[1] + dx[1] - yi) * dx[0];
 
     if (geop == 2) {
         volr *= (x0[0] + x0[0] + dx[0]);
     }
     dx00[0] = dx[0];
     dx00[1] = MAX(0.0, yi - x0[1]);
 
     sph_rec2d0_ur(geop, dim, r, ctr, x0, dx00, vol);
 
     *vol += volr;
 
     return;
 }
 
void sph_rec2dy_ur(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0], x0[1] + dx[1])
 
     dim        : input
     r          : input, radius of the sphere
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     The ctr is assumed to be at the upper right of the rectangular 
 
                  o         
                   o            
                ----o------------
                |    o          |
                |     o         |  
                |       o       |
                |        o      |
                |          o    |
                -------------o---
              x0            i  o  
                                 o
**/
     double r2, xi, x2, y2, yi, volr, dx00[2];
 
     r2 = r * r;
 
//   intersection of the line y = x0[1] and the circle
 
     yi = x0[1];
     y2 = (yi - ctr[1]) * (yi - ctr[1]);
     x2 = MAX(r2 - y2, 0.0);
     xi = ctr[0] - sqrt(x2);
     volr = (x0[0] + dx[0] - xi) * dx[1];
 
     if (geop == 2) {
         volr *= (x0[0] + dx[0] + xi);
     }
     dx00[0] = MAX(0.0, xi - x0[0]);
     dx00[1] = dx[1];
 
     sph_rec2d0_ur(geop, dim, r, ctr, x0, dx00,  vol);
 
     *vol += volr;
 
     return;
 }
 
void sph_rec2dxy_ur(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : input, the center of sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     The ctr is assumed to be at the upper right of the rectangular 
 
          o     -------------      
           o    |           |    
                |           |   
              o |           |  
                o 1         |   
                | o         |    
                ------o------
                      2    o    
                                  o
**/
     double r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
     double dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] and the circle
 
     xi1 = x0[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] - sqrt(y2);
 
//   intersection of the line y = x0[1] and the circle
 
     yi2 = x0[1];
     y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
     x2  = MAX(r2 - y2, 0.0);
     xi2 = ctr[0] - sqrt(x2);
 
     if (geop == 1) {
         x2  = xi2 - xi1;
         y2  = yi1 - yi2;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   area of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi2 - xi1) * (yi1 - yi2);
 
         *vol = dx[0] * dx[1] - volt + vols;
     }
     else if (geop == 2) {

         vols = dx[0] * (x0[1] + dx[1] - yi1) * (x0[0] + x0[0] + dx[0]); 
         volt = (yi1 - x0[1]) * (x0[0] + dx[0] - xi2) * (x0[0] + dx[0] + xi2);
 
         dx00[0] = MAX(0.0, xi2 - x0[0]);
         dx00[1] = MAX(0.0, yi1 - x0[1]);
         sph_rec2d0_ur(geop, dim, r, ctr, x0, dx00, vol);
 
         *vol += (vols + volt);
     }
     return;
 }

void sph_rec2d0_lr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only one vertex, x0, of the the rectangular.
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : ctr of sphere  
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     The ctr is assumed to be at the lower right of the rectangular 
 
                -------------
                |           |   
                |           |   o  
                |           1    
                |        o  |
                ----.--2-----
                        
                     o
**/
     double third, r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
 
     third = 0.33333333333333333333333;
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] + dx[0] and the circle
 
     xi1 = x0[0] + dx[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] + sqrt(y2);
 
     if (geop == 1) {
 
//       intersection of the line y = x0[1] and the circle
 
         yi2 = x0[1];
         y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
         x2  = MAX(r2 - y2, 0.0);
         xi2 = ctr[0] - sqrt(x2);
 
         x2  = xi1 - xi2;
         y2  = yi1 - yi2;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   vol of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi1 - xi2) * (yi1 - yi2);
         *vol = volt + vols;
     }
     else if (geop == 2) {

//       vol = integral from yi2 to yi1  pi * (xi1^2 - x^2) dy 

         *vol = integral_r(r, ctr, x0, dx, x0[1], yi1); 
     }
     return;
 }
 
void sph_rec2dx_lr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim} : input, the center of the sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.

     The ctr is assumed to be at the lower right of the rectangular 

                -------------      
                |           |     
                |           |      o  
                |           o  
                |    o      |  
                o i         |      
                |           |        
            o   -------------     
               x0
**/
     double r2, xi, x2, y2, yi, volr, x00[2], dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] and the circle
 
     xi = x0[0];
     x2 = (xi - ctr[0]) * (xi - ctr[0]);
     y2 = MAX(r2 - x2, 0.0);
     yi = ctr[1] + sqrt(y2);
     volr = (yi - x0[1]) * dx[0];
 
     if (geop == 2) {
         volr *= (x0[0] + x0[0] + dx[0]);
     }
     x00[0] = xi;
     x00[1] = yi;  
     dx00[0] = dx[0];
     dx00[1] = MAX(0.0, x0[1] + dx[1] - yi);
 
     sph_rec2d0_lr(geop, dim, r, ctr, x00, dx00, vol);
 
     *vol += volr;
 
     return;
 }
 
void sph_rec2dy_lr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0], x0[1] + dx[1])
 
     dim        : input
     r          : input, radius of the sphere
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     The ctr is assumed to be at the lower right of the rectangular 
 
                            
                              i     o 
                --------------o--
                |               |
                |        o      |  
                |               |
                |     o         |
                |               |
                ----o------------
              x0                  
                   o              
**/
     double r2, xi, x2, y2, yi, volr, dx00[2];
 
     r2 = r * r;
 
//   intersection of the line y = x0[1] + dx[1] and the circle
 
     yi = x0[1] + dx[1];
     y2 = (yi - ctr[1]) * (yi - ctr[1]);
     x2 = MAX(r2 - y2, 0.0);
     xi = ctr[0] - sqrt(x2);
     volr = (x0[0] + dx[0] - xi) * dx[1];
 
     if (geop == 2) {
         volr *= (x0[0] + dx[0] + xi);
     }
     dx00[0] = MAX(0.0, xi - x0[0]);
     dx00[1] = dx[1];
 
     sph_rec2d0_lr(geop, dim, r, ctr, x0, dx00, vol);
 
     *vol += volr;
 
     return;
 }
 
void sph_rec2dxy_lr(int geop, int dim, double r, double *ctr, double *x0, double *dx, double *vol)
{
/**
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only two vertices,
       (x0[0], x0[1])
       (x[0] + dx[0], x0[1])
 
     dim        : input
     r          : input, radius of the sphere
     ctr[0:dim) : input, the center of sphere 
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     The ctr is assumed to be at the lower right of the rectangular 
                              
                             o 
                    2   o 
                -----o-------      
                | o         |    
              1 o           |   
              o |           |  
                |           |   
            o   |           |    
                -------------
**/

     double r2, xi1, yi1, xi2, yi2, x2, y2, c, a, sina, vols, volt;
     double x00[2], dx00[2];
 
     r2 = r * r;
 
//   intersection of the line x = x0[0] and the circle
 
     xi1 = x0[0];
     x2  = (xi1 - ctr[0]) * (xi1 - ctr[0]);
     y2  = MAX(r2 - x2, 0.0);
     yi1 = ctr[1] + sqrt(y2);
 
//   intersection of the line y = x0[1] + dx[1] and the circle
 
     yi2 = x0[1] + dx[1];
     y2  = (yi2 - ctr[1]) * (yi2 - ctr[1]);
     x2  = MAX(r2 - y2, 0.0);
     xi2 = ctr[0] - sqrt(x2);
 
     if (geop == 1) {
         x2  = xi2 - xi1;
         y2  = yi2 - yi1;
         c   = sqrt(y2 * y2 + x2 * x2);
 
    //   area of the segment
         a   = 2.0 * asin(0.5 * c/r);
         sina  = sin(a);
         vols = 0.5 * r2 *(a - sina);
 
    //   area of the right triamgle
 
         volt = 0.5 *(xi2 - xi1) * (yi2 - yi1);
 
         *vol = dx[0] * dx[1] - volt + vols;
     }
     else if (geop == 2) {

         vols = dx[0] * (yi1 - x0[1]) * (x0[0] + x0[0] + dx[0]); 
         volt = (yi2 - yi1) * (x0[0] + dx[0] - xi2) * (x0[0] + dx[0] + xi2);
 
         x00[0]  = x0[0];
         x00[1]  = yi1;
         dx00[0] = MAX(0.0, xi2 - xi1);
         dx00[1] = MAX(0.0, yi2 - yi1);
         sph_rec2d0_lr(geop, dim, r, ctr, x00, dx00, vol);
 
         *vol += (vols + volt);
     }
     return;
 }

 
void sph_rec3d(int dim, double r, double *x0, double *dx, double vcell, double *vol)
{
/**
 *   It is assume all vertices of the rectangular,
 *   x0[0] >= 0. x0[1] >= 0, x0[2] > 0
 *
 *           7-----------6
 *          /.         / |
 *         / .        /  |         z
 *        4----------5   |         ^
 *        |  3.......|...2         |     . y
 *        | .        |  /          |   .
 *        |.         | /           | .
 *        0----------1             --------->  x
 */
     int i, d0, d1, d2;
     double xl[3], dxe[3], volc;
 
     d0 = 0;
     for (i = 1; i < dim; i++) {
         if (dx[i] < dx[d0]) {
             d0 = i;
         }
     }
     d1 = (d0 + 1) % dim;
     d2 = (d1 + 1) % dim;
     xl[0] = x0[d0];
     xl[1] = x0[d1];
     xl[2] = x0[d2];
     dxe[0] = dx[d0];
     dxe[1] = dx[d1];
     dxe[2] = dx[d2];
 
     sph_rec3dx(dim, r, xl, dxe, vcell, vol);
 
     return;
}
 
void sph_rec3dx(int dim, double r, double *x0, double *dx, double vcell, double *vol)
{
/**
 *   It is assume all vertices of the rectangular,
 *   x0[0] >= 0. x0[1] >= 0, x0[2] > 0
 *
 *           7-----------6
 *          /.         / |
 *         / .        /  |         z
 *        4----------5   |         ^
 *        |  3.......|...2         |     . y
 *        | .        |  /          |   .
 *        |.         | /           | .
 *        0----------1             --------->  x
 */
 
     int i, idir, szdim;
     double r2, vol0, ds2, xr[3], xl2[3], xr2[3];
     double dsmn, dsmx, x0tmp, x1tmp;
     double c8[8][3], ds8[8];
     double xl[3], dxe[3];
     double borc, x0s2, xs2, xrr2;
     double x0s, xs, xrr, v1, v2, v3;
 
     szdim = dim * sizeof(double);
 
     r2 = r * r;
     vol0 = 1.0;
     for (i = 0; i < dim; i++) {
         xr[i] = x0[i] + dx[i];
         xr2[i] = xr[i] * xr[i];
         xl2[i] = x0[i] * x0[i];
 
         vol0 *= dx[i];
     }
     c8[0][0] = x0[0];
     c8[0][1] = x0[1];
     c8[0][2] = x0[2];
 
     c8[1][0] = xr[0];
     c8[1][1] = x0[1];
     c8[1][2] = x0[2];
 
     c8[2][0] = xr[0];
     c8[2][1] = xr[1];
     c8[2][2] = x0[2];
 
     c8[3][0] = x0[0];
     c8[3][1] = xr[1];
     c8[3][2] = x0[2];
 
     c8[4][0] = x0[0];
     c8[4][1] = x0[1];
     c8[4][2] = xr[2];
 
     c8[5][0] = xr[0];
     c8[5][1] = x0[1];
     c8[5][2] = xr[2];
 
     c8[6][0] = xr[0];
     c8[6][1] = xr[1];
     c8[6][2] = xr[2];
 
     c8[7][0] = x0[0];
     c8[7][1] = xr[1];
     c8[7][2] = xr[2];
 
     for (i = 0; i < 8; i++) {
         ds2 = 0.0;
         for (idir = 0; idir < 3; idir++) {
             ds2 += (c8[i][idir] * c8[i][idir]);
         }
         ds8[i] = ds2;
     }
     dsmn = ds8[0];
     dsmx = ds8[0];
     for (i = 0; i < 8; i++) {
         if (ds8[i] < dsmn) {
             dsmn = ds8[i];
         }
         else if (dsmx < ds8[i]) {
             dsmx = ds8[i];
         }
     }
     if (dsmx <= r2) {
         *vol = vol0;
         return;
     }
     else if (dsmn >= r2) {
         *vol = 0.0;
         return;
     }
//   ssph_rec3d_general(dim, r, x0, dx, vcell, vol);
     ssph_rec3d_s(dim, r, x0, dx, vcell, vol);
     return;
 }
 
void sph_rec3d_b(int dim, double r, double *xl, double *dx, double vcell, double *vol)
{
/**
 *   It is assume all vertices of the rectangular,
 *   xl[0] >= 0. x0[1] >= 0, xl[2] > 0
 *
 *           7-----------6
 *          /.         / |
 *         / .        /  |         z
 *        4----------5   |         ^
 *        |  3.......|...2         |     . y
 *        | .        |  /          |   .
 *        |.         | /           | .
 *        0----------1             --------->  x
 */
 
     int    i;
     double r2, xr[3], xl2[3], xr2[3];
     double xlower, xupper;
     double borc, x1ss2, x1s2, xrr2;
     double x1ss, x1s, xrr, v1, v2;
 
     r2 = r * r;
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
         xr2[i] = xr[i] * xr[i];
         xl2[i] = xl[i] * xl[i];
     }
     x1s2 = r2 - xr2[1] - xl2[2];
     xrr2 = r2 - xl2[1] - xl2[2];
 
     x1ss2 = r2 - xl2[1] - xl2[2];
     x1ss  = sqrt(MAX(x1ss2, 0.0));
     borc  = r2 - xl2[0] - xr2[2] - xl2[1];
     assert(borc <= small);
 
     if (xr2[0] <= xrr2) {
         if (x1s2 < xl2[0]) {
             xupper = MIN(x1ss, xr[0]);
             ssph_rec3d_a2(dim, r, xl, dx, xl[0], xupper, vcell, vol);
         }
         else {
             x1s = sqrt(MAX(x1s2, 0.0));
 
             if (xr[0] <= x1s) {
                 xupper = MIN(x1ss, xr[0]);
                 ssph_rec3d_a1(dim, r, xl, dx, xl[0], xupper, vcell, vol);
             }
             else {
                 ssph_rec3d_a1(dim, r, xl, dx, xl[0], x1s,  vcell, &v1);
                 xupper = MIN(x1ss, xr[0]);
                 ssph_rec3d_a2(dim, r, xl, dx, x1s, xupper, vcell, &v2);
                 *vol = v1 + v2;
             }
         }
     }
     else if (xr2[0] >= xrr2) {
         xrr = sqrt(MAX(xrr2, 0.0));
         if (x1s2 < xl2[0]) {
             xupper = MIN(x1ss, xrr);
             ssph_rec3d_a2(dim, r, xl, dx, xl[0], xupper, vcell, vol);
         }
         else {
             x1s = sqrt(MAX(x1s2, 0.0));
             if (xrr <= x1s) {
                 xupper = MIN(x1ss, xrr);
                 ssph_rec3d_a1(dim, r, xl, dx, xl[0], xrr, vcell, vol);
             }
             else {
                 ssph_rec3d_a1(dim, r, xl, dx, xl[0], x1s, vcell, &v1);
                 xupper = MIN(x1ss, xrr);
                 ssph_rec3d_a2(dim, r, xl, dx, x1s, xupper, vcell, &v2);
                 *vol = v1 + v2;
             }
         }
     }
     return;
 }
 
void bintegral(double rsphere, double *d3, double *val, int *np)
{
//   This Boole's formula converges much more slowly than Simpson's.
 
     int   i, i2, npart, npart2, nit;
     double c7, c32, c12;
     double rsphere2, a2, b2, x0, x1, dx, hdx, qdx, x, x2, r2, u, v;
     double ratio, arctan, dval, val_old, err;
     double *fend, *fm, *fl, *fr, *tmp;
 
 
     x   = 1.0/90.0;
     c7  = x * 7.0;
     c12 = x * 12.0;
     c32 = x * 32.0;
 
     rsphere2 = rsphere * rsphere;
     b2       = rsphere2 - d3[1] * d3[1];
     a2       = rsphere2 - d3[2] * d3[2];
     x0       = d3[0];
     x1       = sqrt(MAX(rsphere2 - d3[1] * d3[1] - d3[2] * d3[2], 0.0));
 
     npart = 4;
     fend  = (double *) malloc((4*npart + 1) * sizeof(double));
     fm    = fend + (npart + 1);
     fl    = fm + npart;
     fr    = fl + npart;
 
     dx = (x1 - x0)/(double)npart;
     x  = x0;
     for (i = 0; i <= npart; i++) {
         x2 = x * x;
         r2 = rsphere2 - x2;
         u  = sqrt(MAX(a2 - x2, 0.0)) /(d3[2] + tiny);
         v  = sqrt(MAX(b2 - x2, 0.0)) /(d3[1] + tiny);
/***
         ratio = 1.0 - u * v;
         ratio = (u + v) / ratio;
         arctan = atan(ratio);
***/
         u = atan(u);
         v = atan(v);
         arctan = u + v;
 
         fend[i] = r2 * arctan;
         x += dx;
     }
     x = x0 + hdx;
     for (i = 0; i < npart; i++) {
         x2 = x * x;
         r2 = rsphere2 - x2;
         u  = sqrt(MAX(a2 - x2, 0.0)) /(d3[2] + tiny);
         v  = sqrt(MAX(b2 - x2, 0.0)) /(d3[1] + tiny);
/***
         ratio = 1.0 - u * v;
         ratio = (u + v) / ratio;
         arctan = atan(ratio);
***/
         u = atan(u);
         v = atan(v);
         arctan = u + v;
         fm[i] = r2 * arctan;
         x += dx;
     }
     err = 1.0;
     nit = 0;
     while (err > tol) {
         dx = (x1 - x0)/(double)npart;
         hdx = 0.5  * dx;
         qdx = 0.25 * dx;
         *val = 0.0;
         x = x0 + qdx;
         for (i = 0; i < npart; i++) {
             x2 = x * x;
             r2 = rsphere2 - x2;
             u  = sqrt(MAX(a2 - x2, 0.0)) /(d3[2] + tiny);
             v  = sqrt(MAX(b2 - x2, 0.0)) /(d3[1] + tiny);
/****
             ratio = 1.0 - u * v;
             ratio = (u + v) / ratio;
             arctan = atan(ratio);
***/
             u = atan(u);
             v = atan(v);
             arctan = u + v;
 
             fl[i]  = r2 * arctan;
 
             x += hdx;
             x2 = x * x;
             r2 = rsphere2 - x2;
             u  = sqrt(MAX(a2 - x2, 0.0)) /(d3[2] + tiny);
             v  = sqrt(MAX(b2 - x2, 0.0)) /(d3[1] + tiny);
/****
             ratio = 1.0 - u * v;
             ratio = (u + v) / ratio;
             arctan = atan(ratio);
***/
             u = atan(u);
             v = atan(v);
             arctan = u + v;
 
             fr[i]  = r2 * arctan;
             dval = dx*(c7 *fend[i] + c32*fl[i] + c12*fm[i] + c32*fr[i] + c7*fend[i+1]);
             *val  += dval;
 
             x += hdx;
         }
         if (nit) {
             err = fabs(*val - val_old)/val_old;
         }
         if (err > tol) {
             npart2 = npart + npart;
             tmp = (double *) malloc((4*npart2 + 1) * sizeof(double));
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 tmp[i2]   = fend[i];
                 tmp[i2+1] = fm[i];
             }
             tmp[npart2] = fend[npart];
             fm = tmp + (npart2 + 1);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 fm[i2]   = fl[i];
                 fm[i2+1] = fr[i];
             }
             fl = fm + npart2;
             fr = fl + npart2;
 
             free(fend);
             fend = tmp;
 
             npart = npart2;
             nit++;
         }
         val_old = *val;
     }
     free(fend);
 
     *np = npart;
 
     return;
 }
 
void sintegral(double rsphere, double *d3, double x1, double *val, int *np)
{
//   Simpson rule
 
     int   i, i2, npart, npart2, nit;
     double sixth, rsphere2, a2, b2, x0,  dx, dx6, x, x2, r2, u, v;
     double pi2, ratio, arctan, dval, val_old, err;
     double *fend, *fm, *tmp;
 
     sixth    = 0.1666666666666666666667;
     pi2      = 0.5 * PI;
 
     rsphere2 = rsphere * rsphere;
     b2       = rsphere2 - d3[1] * d3[1];
     a2       = rsphere2 - d3[2] * d3[2];
     x0       = d3[0];
 
     npart = 4;
     fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
     fm    = fend + (npart + 1);
 
     assert(x1 >= x0);
 
     dx = (x1 - x0)/(double)npart;
     x  = x0 - dx;
     for (i = 0; i <= npart; i++) {
         x += dx;
         x2 = x * x;
         r2 = rsphere2 - x2;
         u  = sqrt(MAX(a2 - x2, 0.0)) /(d3[2] + tiny);
         v  = sqrt(MAX(b2 - x2, 0.0)) /(d3[1] + tiny);
 
/***
         ratio = 1.0 - u * v;
         ratio = (u + v) / ratio;
         arctan = atan(ratio);
***/
 
         u = atan(u);
         v = atan(v);
         arctan = u + v;
 
         fend[i] = r2 *(arctan - pi2);
     }
     err = 1.0;
     nit = 0;
     while (err > tol) {
         dx = (x1 - x0)/(double)npart;
         dx6 = sixth * dx;
         *val = 0.0;
         x = x0 - 0.5 * dx;
         for (i = 0; i < npart; i++) {
             x += dx;
             x2 = x * x;
             r2 = rsphere2 - x2;
             u  = sqrt(MAX(a2 - x2, 0.0)) /(d3[2] + tiny);
             v  = sqrt(MAX(b2 - x2, 0.0)) /(d3[1] + tiny);
 
/***
             ratio = 1.0 - u * v;
             ratio = (u + v) / ratio;
             arctan = atan(ratio);
***/
 
             u = atan(u);
             v = atan(v);
             arctan = u + v;
 
 
             fm[i]  = r2 *(arctan - pi2);
             dval   = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
             *val  += dval;
         }
         if (nit) {
             err = fabs(*val - val_old)/val_old;
         }
         if (err > tol) {
             npart2 = npart + npart;
             tmp = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
             assert(tmp);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 tmp[i2]   = fend[i];
                 tmp[i2+1] = fm[i];
             }
             tmp[npart2] = fend[npart];
             free(fend);
             fend = tmp;
             fm   = fend + (npart2 + 1);
 
             npart = npart2;
             nit++;
         }
         val_old = *val;
     }
     free(fend);
 
     *np = npart;
     return;
 }
 
void ssphere_int(double x0, double x1, double rsphere, double d, double *vol)
{
     int npart, npart2, nit, i, i2;
     double rsphere2, err, dx, dx6, x, r2, r, dbyr, acosa, ang, sina, area;
     double sixth, dvol, vol_old;
     double *fend, *fm, *tmp;
 
     sixth    = 0.1666666666666666666667;
 
     rsphere2 = rsphere * rsphere;
 
     npart = 4;
     fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
     fm    = fend + (npart + 1);
     dx = (x1 - x0)/(double)npart;
     x = x0 - dx;
     for (i = 0; i <= npart; i++) {
         x += dx;
         r2 = MAX(rsphere2 - x * x, 0.0);
         r  = sqrt(r2);
         dbyr = d / r;
         acosa = acos(dbyr);
         ang   = acosa + acosa;
         sina  = sin(ang);
         area  = 0.25 * r2 *(ang - sina);
         fend[i] = area;
     }
     err = 1.0;
     nit = 0;
     while (err > tol) {
         dx = (x1 - x0)/(double)npart;
         dx6 = sixth * dx;
         x   = x0 - 0.5 * dx;
         *vol = 0.0;
         for (i = 0; i < npart; i++) {
             x += dx;
             r2 = MAX(rsphere2 - x * x, 0.0);
             r  = sqrt(r2);
             dbyr = d / r;
             acosa = acos(dbyr);
             ang   = acosa + acosa;
             sina  = sin(ang);
             area  = 0.25 * r2 *(ang - sina);
             area = MAX(area, 0.0);
             fm[i] = area;
             dvol  = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
             *vol  += dvol;
         }
         if (nit) {
             err = fabs(*vol - vol_old)/vol_old;
         }
         if (err > tol) {
             npart2 = npart + npart;
             tmp = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
             assert(tmp);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 tmp[i2]   = fend[i];
                 tmp[i2+1] = fm[i];
             }
             tmp[npart2] = fend[npart];
             free(fend);
             fend = tmp;
             fm   = fend + (npart2 + 1);
 
             npart = npart2;
             nit++;
         }
         vol_old = *vol;
     }
     free(fend);
 
     return;
 }
 
 
void sphere_int(double x0, double x1, double rsphere, double d, double *vol)
{
     int npart, nit, i;
     double rsphere2, err, dx, x, r2, r, dbyr, acosa, ang, sina, area, dvol;
     double vol_old;
 
     rsphere2 = rsphere * rsphere;
 
     npart = 4;
     err = 1.0;
 
     nit = 0;
     while (err > tol) {
         dx = (x1 - x0)/(double)npart;
         *vol = 0.0;
         x = x0 - 0.5 * dx;
         for (i = 0; i < npart; i++) {
             x += dx;
             r2 = MAX(rsphere2 - x * x, 0.0);
             r  = sqrt(r2);
             dbyr = d / r;
             acosa = acos(dbyr);
             ang   = acosa + acosa;
             sina  = sin(ang);
             area  = 0.25 * r2 *(ang - sina);
             dvol  = area * dx;
             *vol  += dvol;
         }
         if (nit) {
             err = fabs(*vol - vol_old)/vol_old;
         }
         npart += npart;
         nit++;
         vol_old = *vol;
     }
     return;
 }
 
void dereax(double rsphere, double d2, double fixedx, double y, double *f)
{
     double rsphere2, x2, y2, r2, ry2, rt, tmp, atantmp;
 
     rsphere2 = rsphere * rsphere;
     x2       = fixedx  * fixedx;
     r2       = rsphere2 - x2;
     y2       = y * y;
 
     if (y2 != r2) {
         ry2      = MAX(r2 - y * y, 0.0);
         rt       = sqrt(ry2);
         tmp      = y/(rt + tiny);
         atantmp  = atan(tmp);
         *f       = 0.5 *(y * rt + r2 * atantmp) - d2 * y;
     }
     else {
         *f       = 0.25 * PI * r2;
     }
     return;
  }
 
 
 
void sph_rec3d0_old(int dim, double r, double *x0, double *dx, double *vol)
{
/**
     This has a big problem:
 
     solution (a small number) = big_number_1 - big_number_2.
 
 
     This is to calculate the volume of the inner part covered by the sphere.
     The sphere is center at (0,0). This is for a special case in which
     the sphere covers only one vertex, x0, of the the rectangular.
 
     dim        : input
     r          : input, radius of the sphere
     x0[0:dim)  : input, the coordinate of the rectangular of the lower coarner.
     dx[0:dim)  : input, the width of the rectangular.
     vol        : output, the volume of the part of the rectangular covered
                  by the sphere.
 
     It is assume all vertices of the rectangular,
     x0[0] >= 0. x0[1] >= 0, x0[2] > 0
 
                -------------
                |           |
      .         |           |
           .    |           |
               .|           |
                | .         |
                ----.--------
                      .
                         .
                          .
**/
     int    i;
     double pi6, pi24;
     double r2, r0, h, a2, a, d, vxcap, vycap, vzcap, vxy, vxz, vyz, vsphere, vr;
 
     pi6  = 0.16666666666666666667 * PI;
     pi24 = 0.25 * pi6;
 
     r2 = r * r;
 
     r0 = 0.0;
     for (i = 0; i < dim; i++) {
         r0 += (x0[i] * x0[i]);
     }
     r0 = sqrt(r0);
     if (r0 >= r) {
         *vol = 0.0;
         return;
     }
//   calculate the spherical cap cut by the plane x = x0[0].
 
     d  = x0[0];
     h  = r - d;
     a2 = r2 - d * d;
     vxcap = pi24 * h *(3.0 * a2 + h * h);
 
//   calculate the spherical cap cut by the plane y = x0[1].
 
     d  = x0[1];
     h  = r - d;
     a2 = r2 - d * d;
     vycap = pi24 * h *(3.0 * a2 + h * h);
 
//   calculate the spherical cap cut by the plane z = x0[2];
 
     d  = x0[2];
     h  = r - d;
     a2 = r2 - d * d;
     vzcap = pi24 * h *(3.0 * a2 + h * h);
 
//   the volume of the rectangular
 
     vr = 1.0;
     for (i = 0; i < dim; i++) {
         vr *= x0[i];
     }
//   the volume of the 8th sphere
     vsphere = pi6 * r2 * r;
 
//   calculate vxy
 
     d  = x0[0];
     a2 = MAX(r2 - d * d, 0.0);
     a  = sqrt(a2);
//     sphere_int(x0[1], a, r, d, &vxy);
     ssphere_int(x0[1], a, r, d, &vxy);
 
//   calculate vxz
 
     d  = x0[0];
     a2 = MAX(r2 - d * d, 0.0);
     a  = sqrt(a2);
//     sphere_int(x0[2], a, r, d, &vxz);
     ssphere_int(x0[2], a, r, d, &vxz);
 
     //   calculate vyz
 
     d  = x0[1];
     a2 = MAX(r2 - d * d, 0.0);
     a  = sqrt(a2);
//     sphere_int(x0[2], a, r, d, &vyz);
     ssphere_int(x0[2], a, r, d, &vyz);
 
     *vol = (vsphere - vr) - (vxcap + vycap + vzcap - (vxy + vxz + vyz));
 
     return;
 }
 
 
/**********
void plane_rec3d(int dim, double *p, double *norm, double *xl, double *dx,
                 double *vol)
{
//   The plane (p0, norm) cut the rectangular into two parts
//   vol[0] refers to the lower side of the norm direction.
//   vol[1] refers to the higher side of the norm difrection.
 
//   vol[2]: output
 
     int  i, nnorm;
     double ptr[3], normal[3], v, tmp;
 
     v = 1.0;
     for (i = 0; i < dim; i++) {
         ptr[i]    = (p[i] - xl[i])/dx[i];
         normal[i] = norm[i]/dx[i];
         v *= dx[i];
     }
     tmp = 0.0;
     for (i = 0; i < dim; i++) {
         tmp += (normal[i] * normal[i]);
     }
     tmp = 1.0/tmp;
     for (i = 0; i < dim; i++) {
         normal[i] *= tmp;
     }
     nnorm = 0;
     for (i = 0; i < dim; i++) {
        if (normal[i] != 0.0) {
            nnorm++;
        }
     }
     if (nnorm == 3) {
         normal[2] = -normal[2];
         cube_plane_3(p, normal, vol);
     }
     else if (nnorm == 2) {
         cube_plane_2(p, norm, vol);
     }
     else if (nnorm == 1) {
          *vol = p[0] * norm[0] + p[1] * norm[1] + p[2] * norm[2];
     }
     *vol *= v;
 
     return;
 }
 
void cube_plane_3(double *p, double *norm, double *vol)
{
     double v, normal[3];
 
     if ((norm[0] > 0.0) && (norm[1] > 0.0) && (norm[2] > 0.0)) {
         cube_plane_ppp(p, norm, vol);
     }
     else if ((norm[0] < 0.0) && (norm[1] > 0.0) && (norm[2] > 0.0)) {
         normal[0] = -norm[0];
         normal[1] =  norm[2];
         normal[2] =  norm[1];
 
         cube_plane_ppp(p, normal, &v);
     }
     else if ((norm[0] > 0.0) && (norm[1] < 0.0) && (norm[2] > 0.0)) {
          normal[0] =  norm[2];
          normal[1] = -norm[1];
          normal[2] =  norm[0];
 
          cube_plane_ppp(p, normal, &v);
     }
     else if ((norm[0] > 0.0) && (norm[1] > 0.0) && (norm[2] < 0.0)) {
          normal[0] =  norm[1];
          normal[1] =  norm[0];
          normal[2] = -norm[2];
 
          cube_plane_ppp(p, normal, &v);
     }
     else if ((norm[0] < 0.0) && (norm[1] < 0.0) && (norm[2] > 0.0)) {
          normal[0] = - norm[0];
          normal[1] = - norm[1];
          normal[2] =   norm[2];
 
          cube_plane_ppp(p, normal, &v);
     }
     else if ((norm[0] < 0.0) && (norm[1] > 0.0) && (norm[2] < 0.0)) {
          normal[0] = - norm[0];
          normal[1] =   norm[1];
          normal[2] = - norm[2];
 
          cube_plane_ppp(p, normal, &v);
     }
     else if ((norm[0] > 0.0) && (norm[1] < 0.0) && (norm[2] < 0.0)) {
          normal[0] =   norm[0];
          normal[1] = - norm[1];
          normal[2] = - norm[2];
 
          cube_plane_ppp(p, normal, &v);
     }
     else if ((norm[0] < 0.0) && (norm[1] < 0.0) && (norm[2] < 0.0)) {
 
          normal[0] = - norm[1];
          normal[1] = - norm[0];
          normal[2] = - norm[2];
 
          cube_plane_ppp(p, normal, &v);
     }
     return;
 }
 
void cube_plane_ppp(double *p, double *norm, double *vol)
{
     double v, normal[3];
 
     if ((norm[2] >= norm[1]) && (norm[1] >= norm[0])) {
         cube_plane_n2n1n0(p, norm, vol);
     }
     else if ((norm[2] >= norm[0]) && (norm[0] >= norm[1])) {
         cube_plane_n2n0n1(p, norm, vol);
     }
     else if ((norm[0] >= norm[2])&&(norm[2] >= norm[1])) {
         normal[0] = norm[1];
         normal[1] = norm[2];
         normal[2] = norm[0];
         cube_plane_n2n1n0(p, normal, &v);
     }
     else if ((norm[0] >= norm[1])&&(norm[1] >= norm[2])) {
         normal[0] = norm[1];
         normal[1] = norm[2];
         normal[2] = norm[0];
         get_ds_n2n0n1(p, normal, &v);
     }
     else if ((norm[1] >= norm[0])&&(norm[0] >= norm[2])) {
         normal[0] = norm[2];
         normal[1] = norm[0];
         normal[2] = norm[1];
         get_ds_n2n1n0(p, normal, &v);
     }
     else if ((norm[1] >= norm[2])&&(norm[2] >= norm[0])) {
         normal[0] = norm[2];
         normal[1] = norm[0];
         normal[2] = norm[1];
         get_ds_n2n0n1(p, normal, &v);
     }
     return;
 }
 
void cube_plane_n2n1n0(double *p, double *norm, double *vol)
{
     int dim;
     double sixth, ds, la[3];
 
     sixth = 0.6666666666666666667;
     dim = 3;
     ds = norm[0] * p[0] + norm[1] * p[1] + norm[2] * p[2];
     if (ds <= norm[0]) {
//       Fig.1
//       cube_1_inv
 
         la[0] = ds / norm[0];   // on e1
         la[1] = ds / norm[1];   // on e3
         la[2] = ds / norm[2];   // on e8
         *vol = sixth * la[0] * la[1] * la[2];
 
//         vol = cube_1(1, dim, norm,  ds, NULL, NULL, NULL,
//                      NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
//                      NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
     }
     if (ds <= norm[1]) {
//       Fig.2
//       cube_2_inv
 
         *vol = cube_2(1, dim, norm, ds, NULL, NULL, NULL,
                      NULL, NULL, NULL, NULL,  NULL, NULL, 0, NULL,
                      0.0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
     }
     if (norm[2] <= norm[0] + norm[1]) {
         if (ds <= norm[2]) {
//           Fig.3a
//          cube_3a_inv
 
             *vol = cube_3a(1, dim, norm, ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[1]) {
//           Fig.4a
//           cube_4a_inv
 
             *vol = cube_4a(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[2]) {
//           Fig.5a
//           cube_5a_inv
 
             *vol = cube_5a(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[1] + norm[2]) {
//           Fig.6a
//           cube_6a_inv
 
             *vol = cube_6a(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[1] + norm[2]) {
//           Fig.7a
//           cube_7a or cube_7a0
 
             *vol = 1.0 - sixth * ds * ds * ds /(norm[0] * norm[1] * norm[2]);
 
//             vol = cube_7a(dim, norm, volume, ds, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
 
         }
         else {
             printf("ERROR: max ds in cube_int_n2n1n0 at location 1\n");
         }
     }
     else {
         if (ds <= norm[0] + norm[1]) {
//           Fig.3b
//           cube_3b_inv
 
             *vol = cube_3b(1, dim, norm, ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
 
         }
         else if (ds <= norm[2]) {
//           Fig.4b
//           cube_4b_inv
 
             *vol = cube_4b(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0, NULL, NULL, NULL, NULL, NULL);
 
         }
         else if (ds <= norm[0] + norm[2]) {
//           Fig.5b
//           cube_5b_inv
 
             *vol = cube_5b(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[1] + norm[2]) {
//           Fig.6b
//           cube_6b_inv
 
             *vol = cube_6b(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[1] + norm[2]) {
//           Fig.7b
//           cube_7b or cube_7b0
 
 
             *vol = 1.0 - sixth * ds * ds * ds /(norm[0] * norm[1] * norm[2]);
 
//             vol = cube_7b(dim, norm, volume, ds, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
         }
         else {
             printf("ERROR: max ds in cube_int_n2n1n0 at location 2\n");
 
         }
     }
     return;
 }
 
void cube_plane_n2n0n1(double *p, double *norm, double *vol)
{
     int dim;
     double sixth, ds, la[3];
 
     sixth = 0.66666666666667;
     dim = 3;
     ds = norm[0] * p[0] + norm[1] * p[1] + norm[2] * p[2];
     if (ds <= norm[1]) {
//       Fig.1
//       cube_11
 
         la[0] = ds / norm[0];   // on e1
         la[1] = ds / norm[1];   // on e3
         la[2] = ds / norm[2];   // on e8
         *vol = sixth * la[0] * la[1] * la[2];
 
//         vol = cube_1(1, dim, norm,  ds, NULL, NULL, NULL,
//                      NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
//                      NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
     }
     if (ds <= norm[0]) {
//       Fig.12
//       cube_12
 
         *vol = cube_12(1, dim, norm, ds, NULL, NULL, NULL,
                      NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                      0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
     }
     if (norm[2] <= norm[0] + norm[1]) {
         if (ds <= norm[2]) {
//           Fig.13a
//          cube_13a
 
             *vol = cube_13a(1, dim, norm, ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[1]) {
//           Fig.14a
//           cube_14a
 
             *vol = cube_14a(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[1] + norm[2]) {
//           Fig.15a
//           cube_15a
 
             *vol = cube_15a(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[2]) {
//           Fig.16a
//           cube_16a
 
             *vol = cube_16a(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[1] + norm[2]) {
//           Fig.17a
//           cube_17a
 
             *vol = 1.0 - sixth * ds * ds * ds /(norm[0] * norm[1] * norm[2]);
 
//             vol = cube_7a(dim, norm, volume, ds, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
 
         }
         else {
             printf("ERROR: max ds in cube_int_n2n1n0 at location 1\n");
         }
     }
     else {
         if (ds <= norm[0] + norm[1]) {
//           Fig.13b
//           cube_13b
 
             *vol = cube_13b(1, dim, norm, ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
 
         }
         else if (ds <= norm[2]) {
//           Fig.14b
//           cube_14b
 
             *vol = cube_14b(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
 
         }
         else if (ds <= norm[1] + norm[2]) {
//           Fig.15b
//           cube_5b
 
             *vol = cube_15b(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[2]) {
//           Fig.16b
//           cube_16b
 
             *vol = cube_16b(1, dim, norm,  ds, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
                           0, NULL, (double)0.0, NULL, NULL, NULL, NULL, NULL);
         }
         else if (ds <= norm[0] + norm[1] + norm[2]) {
//           Fig.17b
//           cube_17b
 
 
             *vol = 1.0 - sixth * ds * ds * ds /(norm[0] * norm[1] * norm[2]);
 
//             vol = cube_17b(dim, norm, volume, ds, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
//                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
         }
         else {
             printf("ERROR: max ds in cube_int_n2n0n1 at location 2\n");
 
         }
     }
     return;
 }
 
void cube_plane_2(double * p, double *norm, double *vol)
{
     double ptr[2], normal[2], mynormal[2];
 
     if (norm[0] == 0.0) {
         normal[0] = norm[1];
         normal[1] = norm[2];
         ptr[0] = p[1];
         ptr[1] = p[2];
     }
     else if (norm[1] == 0.0) {
         normal[0] = norm[2];
         normal[1] = norm[0];
         ptr[0] = p[2];
         ptr[1] = p[0];
     }
     else if (norm[2] == 0.0) {
         normal[0] = norm[0];
         normal[1] = norm[1];
         ptr[0] = p[0];
         ptr[1] = p[1];
     }
     if ((normal[0] > 0.0) && (normal[1] > 0.0)) {
         square_line_pp(ptr, normal, vol);
     }
     else if ((normal[0] < 0.0) && (normal[1] > 0.0)) {
         mynormal[0] =  normal[1];
         mynormal[1] = -normal[0];
         square_line_pp(ptr, mynormal, vol);
     }
     else if ((normal[0] > 0.0) && (normal[1] < 0.0)) {
         mynormal[0] = -normal[1];
         mynormal[1] =  normal[0];
         square_line_pp(ptr, mynormal, vol);
     }
     else if ((normal[0] < 0.0) && (normal[1] < 0.0)) {
         mynormal[0] = -normal[0];
         mynormal[1] = -normal[1];
         square_line_pp(ptr, mynormal, vol);
     }
     return;
 }
************/
 
void square_line_pp(double *p, double *norm, double *area)
{
     if (norm[1] >= norm[0]) {
         square_line_n1n0(p, norm, area);
     }
     else {
         square_line_n0n1(p, norm, area);
     }
     return;
 }
 
void square_line_n1n0(double *p, double *norm, double *area)
{
//   area is the lower part along the norm direction.
 
     double s, ds;
 
     s = p[0] * norm[0] + p[1] * norm[1];
     if (s <= norm[0]) {
         *area = 0.5 * s * s /(norm[0] * norm[1]);
     }
     else if (s <= norm[1]) {
         *area = 0.5 * norm[0] / norm[1] + (s - norm[0])/norm[1];
     }
     else if (s <= norm[0] + norm[1]) {
         ds = norm[0] + norm[1] - s;
         *area = 1.0 - 0.5 * ds * ds / (norm[0] * norm[1]);
     }
     return;
  }
 
void square_line_n0n1(double *p, double *norm, double *area)
{
//   area is the lower part along the norm direction.
 
     double s, ds;
 
     s = p[0] * norm[0] + p[1] * norm[1];
     if (s <= norm[1]) {
         *area = 0.5 * s * s /(norm[0] * norm[1]);
     }
     else if (s <= norm[0]) {
         *area = 0.5 * norm[1] / norm[0] + (s - norm[1])/norm[0];
     }
     else if (s <= norm[0] + norm[1]) {
         ds = norm[0] + norm[1] - s;
         *area = 1.0 - 0.5 * ds * ds / (norm[0] * norm[1]);
     }
     return;
  }
 
 
 
 
 
 
 
 
////////////////////////////////////////////////////////////////////////////
 
void aintegral(double rsphere, double x0, double x1, double p0,
               double *val)
{
//   integral_from_x0_to_x1 sqrt(R^2 - p0^2 - x^2)] dx
 
     double rsphere2, p02, r2, x02, x12, rt1, rt0, u, v, arctan;
     double tmp;
 
     rsphere2 = rsphere * rsphere;
     p02      = p0 * p0;
     r2       = rsphere2 - p02;
     x02      = x0 * x0;
     x12      = x1 * x1;
 
     tmp      = r2 - x12;
     assert(tmp >= -small);
     tmp      = MAX(0.0, tmp);
     rt1      = sqrt(tmp);
     tmp      = r2 - x02 + tiny;
     assert(tmp >= -small);
     tmp      = MAX(0.0, tmp);
     rt0      = sqrt(tmp);
 
     *val  = x1 * rt1 - x0 * rt0;
     u     = x1 / rt1;
     v     = x0 / rt0;
     u     = atan(u);
     v     = atan(v);
     arctan = u - v;
 
     *val  += r2 * arctan;
     *val  *= 0.5;
 
     return;
 }
 
void sintegral_minus(double rsphere, double x0, double x1,
                    double p0, double m0,
                    double vcell, double *val, int *np)
{
/***
   Simpson rule for
 
   integral_from_x0_to_x1
   (R^2 - x^2){ atan [sqrt(R^2 - p0^2 - x^2)]/p0 -
                atan [sqrt(R^2 - m0^2 - x^2)]/m0
              } dx
***/
 
     int    i, i2, npart, npart2, nit;
     double sixth, rsphere2, rp2, rm2;
     double dx, dx6, x, x2, r2, p02, m02, tmpu, tmpv, u, v;
     double pi2, ratio, arctan, dval, val_old, err;
     double *fend, *fm, *tmp;
 
     sixth    = 0.1666666666666666666667;
 
     rsphere2 = rsphere * rsphere;
     p02      = p0 * p0;
     m02      = m0 * m0;
     rp2      = rsphere2 - p02;
     rm2      = rsphere2 - m02;
 
     npart = 4;
     fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
     fm    = fend + (npart + 1);
 
     assert(x1 >= x0);
 
     dx = (x1 - x0)/(double)npart;
     x  = x0 - dx;
     for (i = 0; i <= npart; i++) {
         x += dx;
         x2 = x * x;
         r2 = rsphere2  - x2;
         tmpu = rp2 - x2;
         tmpu = MAX(tmpu, tiny);
//         assert(tmpu >= 0.0);
         u  = (sqrt(tmpu)) /(p0 + tiny);
         tmpv= rm2 - x2;
         tmpv = MAX(tmpv, tiny);
//         assert(tmpv >= 0.0);
         v  = (sqrt(tmpv)) /(m0 + tiny);
 
         u = atan(u);
         v = atan(v);
         arctan = u - v;
 
         fend[i] = r2 * arctan;
     }
     err = 1.0;
     nit = 0;
     while (err > tol) {
         dx = (x1 - x0)/(double)npart;
         dx6 = sixth * dx;
         *val = 0.0;
         x = x0 - 0.5 * dx;
         for (i = 0; i < npart; i++) {
             x += dx;
             x2 = x * x;
             r2 = rsphere2 - x2;
 
             tmpu = rp2 - x2;
             tmpu = MAX(tmpu, tiny);
//             assert(tmpu >= 0.0);
             u  = (sqrt(tmpu)) /(p0 + tiny);
             tmpv = rm2 - x2;
             tmpv = MAX(tmpv, tiny);
//             assert(tmpv >= 0.0);
             v  = (sqrt(tmpv)) /(m0 + tiny);
 
             u = atan(u);
             v = atan(v);
             arctan = u - v;
 
             fm[i]  = r2 * arctan;
             dval   = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
             *val  += dval;
         }
         if (nit) {
             err = fabs(*val - val_old)/vcell;
         }
         if (err > tol) {
             npart2 = npart + npart;
             tmp = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
             assert(tmp);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 tmp[i2]   = fend[i];
                 tmp[i2+1] = fm[i];
             }
             tmp[npart2] = fend[npart];
             free(fend);
             fend = tmp;
             fm   = fend + (npart2 + 1);
 
             npart = npart2;
             nit++;
         }
         val_old = *val;
     }
     free(fend);
 
     *np = npart;
     return;
 }
 
void sintegral_plus(double rsphere, double x0, double x1,
                    double z0, double y0, double f0,
                    double vcell, double *val, int *np)
{
/***
   Simpson rule for
 
   integral_from_x0_to_x1
   (R^2 - x^2){ atan [sqrt(R^2 - z0^2 - x^2)]/z0 +
                atan [sqrt(R^2 - y0^2 - x^2)]/y0 -
                f0
              } dx
***/
 
     int    i, i2, npart, npart2, nit;
     double sixth, rsphere2, rz2, ry2, dx, dx6;
     double x, x2, r2, z02, y02, tmpu, tmpv, u, v;
     double pi2, ratio, arctan, dval, val_old, err;
     double *fend, *fm, *tmp;
 
     sixth    = 0.1666666666666666666667;
 
     rsphere2 = rsphere * rsphere;
     z02      = z0 * z0;
     y02      = y0 * y0;
     rz2      = rsphere2 - z02;
     ry2      = rsphere2 - y02;
 
     npart = 4;
     fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
     fm    = fend + (npart + 1);
 
     assert(x1 >= x0);
 
     dx = (x1 - x0)/(double)npart;
     x  = x0 - dx;
     for (i = 0; i <= npart; i++) {
         x += dx;
         x2 = x * x;
         r2 = rsphere2  - x2;
         tmpu = rz2 - x2;
         tmpu = MAX(tmpu, tiny);
//         assert(tmpu >= 0);
         u  = (sqrt(tmpu)) /(z0 + tiny);
         tmpv = ry2 - x2;
         tmpv = MAX(tmpv, tiny);
//         assert(tmpv >= 0.0);
         v  = (sqrt(tmpv)) /(y0 + tiny);
 
/***
         ratio = 1.0 - u * v;
         ratio = (u + v) / ratio;
         arctan = atan(ratio);
***/
         u = atan(u);
         v = atan(v);
         arctan = u + v;
 
         fend[i] = r2 *(arctan - f0);
     }
     err = 1.0;
     nit = 0;
     while (err > tol) {
         dx = (x1 - x0)/(double)npart;
         dx6 = sixth * dx;
         *val = 0.0;
         x = x0 - 0.5 * dx;
         for (i = 0; i < npart; i++) {
             x += dx;
             x2 = x * x;
             r2 = rsphere2 - x2;
             tmpu = rz2 - x2;
             tmpu = MAX(tmpu, tiny);
//             assert(tmpu >= 0.0);
             u  = (sqrt(tmpu)) /(z0 + tiny);
             tmpv = ry2 - x2;
             tmpv = MAX(tmpv, tiny);
//             assert(tmpv >= 0.0);
             v  = (sqrt(tmpv)) /(y0 + tiny);
 
/***
             ratio = 1.0 - u * v;
             ratio = (u + v) / ratio;
             arctan = atan(ratio);
***/
 
             u = atan(u);
             v = atan(v);
             arctan = u + v;
 
 
             fm[i]  = r2 *(arctan - f0);
             dval   = dx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
             *val  += dval;
         }
         if (nit) {
             err = fabs(*val - val_old)/vcell;
         }
         if (err > tol) {
             npart2 = npart + npart;
             tmp = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
             assert(tmp);
             for (i = 0; i < npart; i++) {
                 i2 = i + i;
                 tmp[i2]   = fend[i];
                 tmp[i2+1] = fm[i];
             }
             tmp[npart2] = fend[npart];
             free(fend);
             fend = tmp;
             fm   = fend + (npart2 + 1);
 
             npart = npart2;
             nit++;
         }
         val_old = *val;
     }
     free(fend);
 
     *np = npart;
     return;
 }
 
 
void ssph_rec3d_a1(int dim, double rsphere, double *xl, double *dx,
                   double xlower, double xupper, double vcell, double *vol)
{
     int    i, npart;
     double pi2, xr[3], xr2[3], xl2[3];
     double rsphere2, xrr, volc;
     double aint, bint, cint;
 
     if (xupper <= xlower) {
         *vol = 0.0;
         return;
     }
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
         xr2[i] = xr[i] * xr[i];
         xl2[i] = xl[i] * xl[i];
     }
     rsphere2 = rsphere * rsphere;
 
     volc = - xl[2] * dx[1] * (xupper- xlower);
     aintegral(rsphere, xlower, xupper, xr[1], &aint);
     aintegral(rsphere, xlower, xupper, xl[1], &bint);
     sintegral_minus(rsphere, xlower, xupper, xl[1], xr[1], vcell, &cint, &npart);
 
     *vol = volc + 0.5 *(xr[1] * aint - xl[1] * bint + cint);
 
     return;
 }
 
void ssph_rec3d_a2(int dim, double rsphere, double *xl, double *dx,
                   double xlower, double xupper, double vcell, double *vol)
{
     int    i, npart;
     double pi2, xr[3], xr2[3], xl2[3];
     double rsphere2, xrr, volc;
     double aint, bint, cint;
 
     if (xupper <= xlower) {
         *vol = 0.0;
         return;
     }
     pi2 = 0.5 * PI;
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
         xr2[i] = xr[i] * xr[i];
         xl2[i] = xl[i] * xl[i];
     }
     rsphere2 = rsphere * rsphere;
 
     volc = xl[1] * xl[2] *(xupper- xlower);
     aintegral(rsphere, xlower, xupper, xl[2], &aint);
     aintegral(rsphere, xlower, xupper, xl[1], &bint);
     sintegral_plus(rsphere, xlower, xupper, xl[2], xl[1], pi2, vcell, &cint, &npart);
 
     *vol = volc + 0.5 *(-xl[2] * aint - xl[1] * bint + cint);
 
     return;
 }
 
void ssph_rec3d_a3(int dim, double rsphere, double *xl, double *dx, double *vol)
{
     int    i;
     double xr[3], xr2[3], xl2[3];
     double rsphere2, r2, x0s2, x0s, x1s2;
     double tmp, rt0, rts, u, v;
 
     if (dx[0] <= 0.0) {
         *vol = 0.0;
         return;
     }
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
         xr2[i] = xr[i] * xr[i];
         xl2[i] = xl[i] * xl[i];
     }
     rsphere2 = rsphere * rsphere;
     r2       = rsphere2 - xl2[2];
     x0s2 = MAX(rsphere2 - xr2[1] - xr2[2], 0.0);
     x0s  = sqrt(x0s2);
     x1s2 = r2 - xr2[1];
     if (x1s2 <= xl2[0]) {
         tmp = MAX(r2 - xl2[0], 0.0);
         rt0 = sqrt(tmp);
         tmp = MAX(r2 - x0s2, 0.0);
         rts = sqrt(tmp);
         u   = rt0/xl[0];
         v   = rts/x0s;
         u   = atan(u);
         v   = atan(v);
         *vol = 0.5 * dx[2] *(x0s * rts - xl[0] * rt0  +
                              r2 *(u - v) - 2.0 * xl[1] * (x0s - xl[0]));
     }
     else {
         *vol = dx[2] * dx[1] *(x0s - xl[0]);
     }
     return;
 }
 
 
void ssph_rec3d_general(int dim, double rsphere, double *xl, double *dx,
                        double volc, double *vol)
 
{
     int i, npart, it, itmx, climited;;
     double xr[3], xr2[3], rsphere2, r2, sol, err;
     double x2int, xint, ddx, xm, xm2, rm2, cm2, cm, y1, ys2, ys;
     double tmp, rt0, rt1, u, v, a, dy, dvol, ylow, yhgh;
 
     *vol = 0.0;
 
     for (i = 0; i < dim; i++) {
         xr[i]  = xl[i] + dx[i];
         xr2[i] = xr[i] * xr[i];
     }
     rsphere2 = rsphere * rsphere;
     r2       = rsphere2 - xl[2] * xl[2];
 
//   choose the initial npart
 
//   intersect between sphere r = R and line (y = y0, z = z0)
     x2int = r2 - xl[1] * xl[1];
     if (x2int <= 0) {
         return;
     }
     xint  = sqrt(x2int);
     if ((xint <= xl[0]) || (xint >= xr[0])) {
         ddx = dx[0];
     }
     else {
         ddx = MIN(xint - xl[0], xr[0] - xint);
     }
//   intersect between sphere r = R and line (y = y1, z = z1)
     x2int = rsphere2 - xr2[1] - xr2[2];
     if (x2int > 0.0) {
         xint  = sqrt(x2int);
         if ((xl[0] < xint) && (xint < xr[0])) {
             ddx = MIN(MIN(ddx, xr[0] - xint), xint - xl[0]);
         }
     }
     npart = (int)(dx[0]/ddx) + 32;
 
     itmx = 20;
     sol  = 1000.0 * volc;
     err  = 1.0;
     it   = 0;
 
     while ((err > tol) && (it < itmx)) {
 
         ddx = dx[0] / (double) npart;
         xm  = xl[0] - 0.5 * ddx;
         *vol = 0.0;
         for (i = 0; i < npart; i++) {
             xm += ddx;
             xm2 = xm * xm;
             rm2 = rsphere2 - xm2;
             cm2 = r2 - xm2;
             if (cm2 < 0.0) continue;
             cm  = sqrt(cm2);
             if (cm <= xl[1]) continue;
 
             y1  = MIN(cm, xr[1]);
             ys2 = rm2 - xr2[2];
             if (ys2 < 0.0) {
                 climited = 1;
             }
             else {
                ys = sqrt(ys2);
                if (ys <= xl[1]) {
                    climited = 1;
                }
                else {
                    climited = 0;
                }
             }
             if (climited) {
 
                 tmp = MAX(rm2 - xl[1] * xl[1], 0.0);
                 rt0 = sqrt(tmp);
                 tmp = rm2 - y1 * y1;
                 assert(tmp > -small);
                 tmp = MAX(0.0, tmp);
                 rt1 = sqrt(tmp);
                 u   = y1/rt1;
                 v   = xl[1]/rt0;
                 u   = atan(u);
                 v   = atan(v);
                 a   = 0.5 *(y1 * rt1 - xl[1] * rt0 + rm2 *(u - v))
                     - xl[2] *(y1 - xl[1]);
                 dvol = (ddx * a);
             }
             else {
                 yhgh = MIN(y1, ys);
                 dy = MAX(0.0, yhgh - xl[1]);
                 dvol = ddx * dx[2] * dy;
 
                 ylow = MAX(ys, xl[1]);
                 yhgh = y1;
                 if (yhgh > ylow) {
                     tmp = rm2 - ylow * ylow;
                     tmp = MAX(tmp, tiny);
                     rt0 = sqrt(tmp);
                     tmp = rm2 - yhgh * yhgh;
                     tmp = MAX(tmp, tiny);
                     rt1 = sqrt(tmp);
                     u   = yhgh/(rt1 + tiny);
                     v   = ylow/(rt0 + tiny);
                     u   = atan(u);
                     v   = atan(v);
                     a   = 0.5 *(yhgh * rt1 - ylow * rt0 + rm2 *(u - v))
                         - xl[2] *(yhgh - ylow);
                     dvol += (ddx * a);
                 }
             }
             dvol = MAX(MIN(dvol, volc), 0.0);
 
//             assert(dvol > 0.0);
             *vol += dvol;
         }
         err = fabs(*vol - sol)/volc;
         sol = *vol;
         npart += npart;
         it++;
     }
     return;
 }
 
 
void ssph_rec3d_s(int dim, double rsphere, double *xl, double *dx,
                        double volc, double *vol)
 
{
     int    ifsimpson;
     int    i, i2, npart, npart2, it, itmx, climited;;
     double xr[3], xr2[3], rsphere2, r2, sol, err;
     double x2int, xint, ddx, xm, xm2, rm2, cm2, cm, y1, ys2, ys;
     double fmi, tmp, rt0, rt1, u, v, a, dy, dvol, ylow, yhgh;
     double *fend, *fm, *array, sixth, ddx6;
 
     *vol = 0.0;
 
     sixth = 0.1666666666666666666667;
 
     for (i = 0; i < dim; i++) {
         xr[i]  = xl[i] + dx[i];
         xr2[i] = xr[i] * xr[i];
     }
     rsphere2 = rsphere * rsphere;
     r2       = rsphere2 - xl[2] * xl[2];
 
//   choose the initial npart
 
//   intersect between sphere r = R and line (y = y0, z = z0)
     x2int = r2 - xl[1] * xl[1];
     if (x2int <= 0) {
         return;
     }
     xint  = sqrt(x2int);
     if ((xint <= xl[0]) || (xint >= xr[0])) {
         ddx = dx[0];
     }
     else {
         ddx = MIN(xint - xl[0], xr[0] - xint);
         ddx = MAX(ddx, 0.001 * dx[0]);
     }
//   intersect between sphere r = R and line (y = y1, z = z1)
     x2int = rsphere2 - xr2[1] - xr2[2];
     if (x2int > 0.0) {
         xint  = sqrt(x2int);
         if ((xl[0] < xint) && (xint < xr[0])) {
             ddx = MIN(MIN(ddx, xr[0] - xint), xint - xl[0]);
             ddx = MAX(ddx, 0.001 * dx[0]);
         }
     }
     ifsimpson = 1;
 
     npart = (int)(dx[0]/ddx) + 32;
     fend  = NULL;
 
     if (ifsimpson) {
         fend  = (double *) malloc((npart + npart + 1) * sizeof(double));
         fm    = fend + (npart + 1);
 
         ddx = dx[0] / (double) npart;
         xm  = xl[0] - ddx;
         for (i = 0; i <= npart; i++) {
             xm += ddx;
 
             xm2 = xm * xm;
             rm2 = rsphere2 - xm2;
             cm2 = r2 - xm2;
             if (cm2 < 0.0) {
                 fend[i] = 0.0;
                 continue;
             }
             cm  = sqrt(cm2);
             if (cm <= xl[1]) {
                 fend[i] = 0.0;
                 continue;
             }
             y1  = MIN(cm, xr[1]);
             ys2 = rm2 - xr2[2];
             if (ys2 < 0.0) {
                 climited = 1;
             }
             else {
                ys = sqrt(ys2);
                if (ys <= xl[1]) {
                    climited = 1;
                }
                else {
                    climited = 0;
                }
             }
             if (climited) {
 
                 tmp = MAX(rm2 - xl[1] * xl[1], 0.0);
                 rt0 = sqrt(tmp);
                 tmp = rm2 - y1 * y1;
                 assert(tmp > -small);
                 tmp = MAX(0.0, tmp);
                 rt1 = sqrt(tmp);
                 u   = y1/rt1;
                 v   = xl[1]/rt0;
                 u   = atan(u);
                 v   = atan(v);
                 a   = 0.5 *(y1 * rt1 - xl[1] * rt0 + rm2 *(u - v))
                     - xl[2] *(y1 - xl[1]);
                 fend[i] = a;
             }
             else {
                 yhgh = MIN(y1, ys);
                 dy = MAX(0.0, yhgh - xl[1]);
                 fend[i] = dx[2] * dy;
 
                 ylow = MAX(ys, xl[1]);
                 yhgh = y1;
                 if (yhgh > ylow) {
                     tmp = rm2 - ylow * ylow;
                     tmp = MAX(tmp, tiny);
                     rt0 = sqrt(tmp);
                     tmp = rm2 - yhgh * yhgh;
                     tmp = MAX(tmp, tiny);
                     rt1 = sqrt(tmp);
                     u   = yhgh/(rt1 + tiny);
                     v   = ylow/(rt0 + tiny);
                     u   = atan(u);
                     v   = atan(v);
                     a   = 0.5 *(yhgh * rt1 - ylow * rt0 + rm2 *(u - v))
                         - xl[2] *(yhgh - ylow);
                     fend[i] += a;
                 }
             }
         }
     }
     itmx = 20;
     sol  = 1000.0 * volc;
     err  = 1.0;
     it   = 0;
     while ((err > tol) && (it < itmx)) {
 
         ddx  = dx[0] / (double) npart;
         ddx6 = sixth * ddx;
 
         xm   = xl[0] - 0.5 * ddx;
         *vol = 0.0;
 
         if (ifsimpson) {
 
             for (i = 0; i < npart; i++) {
                 xm += ddx;
                 xm2 = xm * xm;
                 rm2 = rsphere2 - xm2;
                 cm2 = r2 - xm2;
                 if (cm2 < 0.0) {
                     fm[i] = 0.0;
                     dvol = ddx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
                     dvol = MAX(MIN(dvol, volc), 0.0);
                     *vol += dvol;
                     continue;
                 }
                 cm  = sqrt(cm2);
                 if (cm <= xl[1]) {
                     fm[i] = 0.0;
                     dvol = ddx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
                     dvol = MAX(MIN(dvol, volc), 0.0);
                     *vol += dvol;
                     continue;
                 }
                 y1  = MIN(cm, xr[1]);
                 ys2 = rm2 - xr2[2];
                 if (ys2 < 0.0) {
                     climited = 1;
                 }
                 else {
                    ys = sqrt(ys2);
                    if (ys <= xl[1]) {
                        climited = 1;
                    }
                    else {
                        climited = 0;
                    }
                 }
                 if (climited) {
 
                     tmp = MAX(rm2 - xl[1] * xl[1], 0.0);
                     rt0 = sqrt(tmp);
                     tmp = rm2 - y1 * y1;
                     assert(tmp > -small);
                     tmp = MAX(0.0, tmp);
                     rt1 = sqrt(tmp);
                     u   = y1/rt1;
                     v   = xl[1]/rt0;
                     u   = atan(u);
                     v   = atan(v);
                     a   = 0.5 *(y1 * rt1 - xl[1] * rt0 + rm2 *(u - v))
                         - xl[2] *(y1 - xl[1]);
                     fm[i] = a;
                 }
                 else {
                     yhgh = MIN(y1, ys);
                     dy = MAX(0.0, yhgh - xl[1]);
                     fm[i]= dx[2] * dy;
 
                     ylow = MAX(ys, xl[1]);
                     yhgh = y1;
                     if (yhgh > ylow) {
                         tmp = rm2 - ylow * ylow;
                         tmp = MAX(tmp, tiny);
                         rt0 = sqrt(tmp);
                         tmp = rm2 - yhgh * yhgh;
                         tmp = MAX(tmp, tiny);
                         rt1 = sqrt(tmp);
                         u   = yhgh/(rt1 + tiny);
                         v   = ylow/(rt0 + tiny);
                         u   = atan(u);
                         v   = atan(v);
                         a   = 0.5 *(yhgh * rt1 - ylow * rt0 + rm2 *(u - v))
                             - xl[2] *(yhgh - ylow);
                         fm[i] += a;
                     }
                 }
                 dvol = ddx6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
                 dvol = MAX(MIN(dvol, volc), 0.0);
 
    //            assert(dvol > 0.0);
                 *vol += dvol;
             }
         }
         else {
             for (i = 0; i < npart; i++) {
                 xm += ddx;
                 xm2 = xm * xm;
                 rm2 = rsphere2 - xm2;
                 cm2 = r2 - xm2;
                 if (cm2 < 0.0) {
                     fmi = 0.0;
                     continue;
                 }
                 cm  = sqrt(cm2);
                 if (cm <= xl[1]) {
                     fmi = 0.0;
                     continue;
                 }
                 y1  = MIN(cm, xr[1]);
                 ys2 = rm2 - xr2[2];
                 if (ys2 < 0.0) {
                     climited = 1;
                 }
                 else {
                    ys = sqrt(ys2);
                    if (ys <= xl[1]) {
                        climited = 1;
                    }
                    else {
                        climited = 0;
                    }
                 }
                 if (climited) {
 
                     tmp = MAX(rm2 - xl[1] * xl[1], 0.0);
                     rt0 = sqrt(tmp);
                     tmp = rm2 - y1 * y1;
                     assert(tmp > -small);
                     tmp = MAX(0.0, tmp);
                     rt1 = sqrt(tmp);
                     u   = y1/rt1;
                     v   = xl[1]/rt0;
                     u   = atan(u);
                     v   = atan(v);
                     a   = 0.5 *(y1 * rt1 - xl[1] * rt0 + rm2 *(u - v))
                         - xl[2] *(y1 - xl[1]);
                     fmi = a;
                 }
                 else {
                     yhgh = MIN(y1, ys);
                     dy   = MAX(0.0, yhgh - xl[1]);
                     fmi  = dx[2] * dy;
 
                     ylow = MAX(ys, xl[1]);
                     yhgh = y1;
                     if (yhgh > ylow) {
                         tmp = rm2 - ylow * ylow;
                         tmp = MAX(tmp, tiny);
                         rt0 = sqrt(tmp);
                         tmp = rm2 - yhgh * yhgh;
                         tmp = MAX(tmp, tiny);
                         rt1 = sqrt(tmp);
                         u   = yhgh/(rt1 + tiny);
                         v   = ylow/(rt0 + tiny);
                         u   = atan(u);
                         v   = atan(v);
                         a   = 0.5 *(yhgh * rt1 - ylow * rt0 + rm2 *(u - v))
                             - xl[2] *(yhgh - ylow);
                         fmi += a;
                     }
                 }
                 dvol = ddx * fmi;
                 dvol = MAX(MIN(dvol, volc), 0.0);
 
    //            assert(dvol > 0.0);
                 *vol += dvol;
             }
         }
         if (it) {
             err = fabs(*vol - sol)/volc;
         }
         if (err > tol) {
             npart2 = npart + npart;
             if (ifsimpson) {
                 array = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
                 if (array) {   // successfully allocated
                     for (i = 0; i < npart; i++) {
                         i2 = i + i;
                         array[i2]   = fend[i];
                         array[i2+1] = fm[i];
                     }
                     array[npart2] = fend[npart];
                     free(fend);
                     fend = array;
                     fm   = fend + (npart2 + 1);
                 }
                 else {
                     if (fend) free(fend);
                     fend = NULL;
                     ifsimpson = 0;
                 }
             }
             npart = npart2;
             it++;
         }
         sol = *vol;
     }
     if (fend) free(fend);
 
     return;
 }
 
void write_hull(char *fname, double *r1, int n1, double *r2, int n2, double *ri, int ni)
{
/****
         char   name[16];
         int    i, n, fid; 
         int    *list_default;
         double *x, *y, *p1;
         mio_Coord *coord;
         mio_Unstructured_Mesh mesh;
         mio_Group grp;

         n = ni + n1 + n2; 
         list_default = (int *) malloc((n + n) * sizeof(int));
         for (i = 0; i < n+n; i++) { 
             list_default[i] = i;
         } 
         x = (double *) malloc((n + n) * sizeof(double));
         y = x + n;
      
//         mio_open_file(fname, mio_file_create, &fid);

         for (i = 0; i < n1; i++) {   
             p1   = r1 + (i + i);   
             x[i] = p1[0];
             y[i] = p1[1];
         }  
         strcpy(name, "poly1");
         mio_init(fid, mio_umesh, -1, &mesh);
         coord = &(mesh.coord);
         mesh.name       = name;
         mesh.dims       = 2;
         mesh.idmin      = 0;
         mesh.datatype   = mio_int;
         coord->datatype = mio_double;
         coord->coord[0] = x;
         coord->coord[1] = y;

//       polygon with 1441 or more nodes seems not working in Ensight.

         mesh.type     = mio_general_mesh;
         mesh.sizes[1] = 1;
         mesh.sizes[3] = n1;
         mesh.num_nodes_for_face = &n1;
         mesh.nodelist_for_face  = list_default;

         mio_write(fid, mio_umesh, fid, &mesh);

         for (i = 0; i < n2; i++) {
             list_default[i] = i;
         }
         for (i = 0; i < n2; i++) {
             p1   = r2 + (i + i);
             x[i] = p1[0];
             y[i] = p1[1];
         }
         strcpy(name, "poly2");
         mio_init(fid, mio_umesh, -1, &mesh);
         coord = &(mesh.coord);
         mesh.name       = name;
         mesh.dims       = 2;
         mesh.idmin      = 0;

         mesh.type       = mio_general_mesh;
         mesh.datatype   = mio_int;
         coord->datatype = mio_double;
         mesh.sizes[1] = 1;
         mesh.sizes[3] = n2;

         coord->coord[0] = x; 
         coord->coord[1] = y;
         mesh.num_nodes_for_face = &n2;
         mesh.nodelist_for_face  = list_default;
         mio_write(fid, mio_umesh, fid, &mesh);
          

         for (i = 0; i < ni; i++) {
             list_default[i] = i;
         }
         for (i = 0; i < ni; i++) {
             p1   = ri + (i + i);
             x[i] = p1[0];
             y[i] = p1[1];
         }
         strcpy(name, "poly_int");
         mio_init(fid, mio_umesh, -1, &mesh);
         coord = &(mesh.coord);
         mesh.name       = name;
         mesh.dims       = 2;
         mesh.idmin      = 0;
         mesh.type       = mio_general_mesh;
         mesh.datatype   = mio_int;
         coord->datatype = mio_double;

         mesh.sizes[1] = 1;
         mesh.sizes[3] = ni;

         coord->coord[0] = x;
         coord->coord[1] = y;
         mesh.num_nodes_for_face = &ni;
         mesh.nodelist_for_face  = list_default;
         mio_write(fid, mio_umesh, fid, &mesh);

         mio_close_file(fid);

         free(x);
         free(list_default);
****/

     return;
 }
 
void gpoly3d_rec(int nreg, int ifinquiry, int *mixed, int ifcyl, int dim,
                 double *coords, int nnode,
                 int nedge, int nface, int *nodelist_for_edge,
                 int *edgelist_for_face, int *nedge_ea_face,
                 int *nodelist_for_face, int *num_nodes_for_face,
                 double *ctr, double *norm_of_face,
                 double *xl, double *dx, double *vol)
{
//   ctr and ctr norm_of_face are not used.
 
     int sz_nlist, sz_nlist_mx, nnode_out_mx, nface_out_mx;
     int sz, i, nnode_out, nface_out;
     int *list_default;
     double *coords1d_out;
     int    *nnode_for_face_out;
     int    *nodelist_for_face_out;
     double xr[3], ctrp[3], *c;
 
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
     }
 
/*******
 
     int i, p, npoly, nf, nn;
     int  *nnode_ea_face, *nodelist_ea_face;
 
     amr_Loc poly_loc[27], loc;
 
     int list_default[16], nnode_out[27];
     double my_coords1d_out[27][128], *coords1d_out[27];
     int    nface_for_zone_out[27];
     int    my_nnode_for_face_out[27][16], *nnode_for_face_out[27];
     int    my_nodelist_for_face_out[27][128], *nodelist_for_face_out[27];
     double xr[3], ctrp[3], *c;
     double *coords_p;
 
     for (i = 0; i < 27; i++) {
         coords1d_out[i] = my_coords1d_out[i];
         nnode_for_face_out[i] = my_nnode_for_face_out[i];
         nodelist_for_face_out[i] = my_nodelist_for_face_out[i];
     }
     poly_rec3d_unsplit(nz, ifinquiry, 3, xl, xr, nface, nedge, nnode, coords,
                    num_nodes_for_face, nodelist_for_face,
                    edgelist_for_face, nodelist_for_edge,
                    &npoly, poly_loc,
                    nnode_out, coords1d_out,
                    nface_for_zone_out, nnode_for_face_out,
                    nodelist_for_face_out);
     if (ifinquiry) {
         if (npoly == 0) {
             *vol   = 0.0;
             *mixed = 0;
         }
         else {
             *vol   = small;
             *mixed = 1;
         }
     }
     else if (npoly == 0) {
         *vol   = 0.0;
         *mixed = 0;
     }
     else {
         *vol = 0.0;
         for (p = 0; p < npoly; p++) {
             loc = poly_loc[p];
             if (loc != inner) continue;
 
             nn = nnode_out[p];
             coords_p = coords1d_out[p];
             nf = nface_for_zone_out[p];
             nnode_ea_face = nnode_for_face_out[p];
             nodelist_ea_face = nodelist_for_face_out[p];
 
             for (i = 0; i < nf; i++) {
                 list_default[i] = i;
             }
             cal_ctr0(nf, nn, list_default, nnode_ea_face,
                      nodelist_ea_face, coords_p, ctrp, vol);
             break;
         }
     }
************/
 
     sz_nlist = 0;
     for (i = 0; i < nface; i++) {
         sz_nlist += num_nodes_for_face[i];
     }
     sz_nlist_mx  = sz_nlist + sz_nlist + 256;
     nnode_out_mx = nnode + nnode + 32;
     nface_out_mx = nface + nface + 32;
 
     sz = nface_out_mx + nface_out_mx + sz_nlist_mx;
     list_default          = (int    *) malloc(sz * sizeof(int));
     nnode_for_face_out    = list_default + nface_out_mx;
     nodelist_for_face_out = nnode_for_face_out + nface_out_mx;
 
     coords1d_out          = (double *) malloc(nnode_out_mx * 3 * sizeof(double));
 
     poly_rec3d_split(nreg, ifinquiry, 3, xl, xr, nface, nedge, nnode, coords,
                    num_nodes_for_face, nodelist_for_face,
                    edgelist_for_face, nodelist_for_edge,
                    &nnode_out, coords1d_out,
                    &nface_out, nnode_for_face_out,
                    nodelist_for_face_out);
 
     sz_nlist = 0;
     for (i = 0; i < nface_out; i++) {
         sz_nlist += nnode_for_face_out[i];
     }
     if (!ifinquiry) {
         assert(sz_nlist  <= sz_nlist_mx);
         assert(nnode_out <= nnode_out_mx);
         assert(nface_out <= nface_out_mx);
     }
     if (ifinquiry) {
         if (nnode_out == 0) {
             *vol   = 0.0;
             *mixed = 0;
         }
         else {
             *vol   = small;
             *mixed = 1;
         }
     }
     else if (nnode_out == 0) {
         *vol   = 0.0;
         *mixed = 0;
     }
     else {
         *vol = 0.0;
         for (i = 0; i < nface_out; i++) {
             list_default[i] = i;
         }
         cal_ctr0(nface_out, nnode_out, list_default, nnode_for_face_out,
                  nodelist_for_face_out, coords1d_out, ctrp, vol);
     }
     free(list_default);
     free(coords1d_out);
 
     return;
 }
 
//////////////////////////////////////////////////////////////////////////////////////
 
void poly_rec3d_split(int nreg, int ifinquiry, int dim,
                double *xl, double *xr,
                int nface, int nedge,  int nnode,  double *coords1d,
                int *num_node_for_face, int *nodelist_for_face,
                int *edgelist_for_face, int *nodelist_for_edge,
                int *nnode_out, double *coords1d_out,
                int *nface_out, int *nnode_for_face_out,
                int *nodelist_for_face_out)
{
     int  nnode_mx, nedge_mx, nface_mx, sz_nlist_mx;
 
     int  szint, szdim, inside, *inside_ea_node;
     int  anyinside, alloutside, allinside;
     int  i, i1, i2, j, j1, k, idim, lr, lh, copy_read, copy_write;
     int  n, n0, n1, n2, nn, m0, m1, nf, f, f1, e, offset;
     int  sz, sz_nlist;
     int  my_nnode, my_nedge, my_nface;
     int  nnode_this_p, nedge_this_p;
     int  *old2new_nodeid, ifinside[8];
     int  nnode_p[2], nedge_p[2], nface_p[2];
     int  *nnode_for_face_p[2];
     int  *nodelist_for_face_p[2];
     int  *edgelist_for_face_p[2];
     int  *nodelist_for_edge_p[2];
 
     int  nnode_out2, nface_for_zone_out2[2];
     int  *nodelist_for_face_out1d, *nodelist_for_face_out2[2];
     int  *nnode_for_face_out1d,      *nnode_for_face_out2[2];
     int  nnode_intface, *nodelist_intface;
     int  *nodelist_for_edge_work;
     int  *edgelist_for_face_work;
     int  *list_default;
 
     int  *nodelist_s, *edgelist_s, *nnode_for_face_s;
     int  *nodelist_t, *edgelist_t;
     int  *nodelist_exist;
 
     int  *my_num_node_for_face, *my_nodelist_for_face;
     int  *my_edgelist_for_face, *my_nodelist_for_edge;
 
     double *planes, x, *my_coords1d;
     double *coords1d_p[2];
     double *coords1d_out2;
     double xmin[3], xmax[3], dx[3], ctr[3], vol, c3[3];
     double *norm_of_face;
     double *norm, *p0, *p1, *p2, dp0[3], dp1[3], tmp;
     double dot_ctr, dot_cube, dctr, dcube;
     double c8[8][3], a[3][3], factor, cost, sint, cosp, sinp;
 
     double *xs, *ys, *zs;
     double *coord1ds_t;
     double *cs, *ct, *my_coord1d, xplane, dl;
 
     szint  = sizeof(int);
     szdim  = sizeof(double) * dim;
 
//   check whether all nodes of poly are obviously outside the rectangular.
 
     cs = coords1d;
     memcpy(xmin, cs, (size_t)szdim);
     memcpy(xmax, cs, (size_t)szdim);
     cs += dim;
     for (i = 1; i < nnode; i++) {
         for (k = 0; k < dim; k++) {
             if (cs[k] < xmin[k]) {
                 xmin[k] = cs[k];
             }
             else if (cs[k] > xmax[k]) {
                 xmax[k] = cs[k];
             }
         }
         cs += dim;
     }
     alloutside = 0;
     for (k = 0; k < dim; k++) {
         if ((xmax[k] <= xl[k]) || (xmin[k] >= xr[k]))  {
             alloutside = 1;
             break;
         }
     }
     if (alloutside) {
         *nnode_out = 0;
         *nface_out = 0;
         return;
     }
     norm_of_face   = (double *) malloc(nface * 3 * sizeof(double));
 
     nnode_mx    = nnode + nnode + 32;
     nedge_mx    = nedge + nedge + 32;
     nface_mx    = nface + nface + 32;
 
     sz_nlist = 0;
     for (i = 0; i < nface; i++) {
         sz_nlist += num_node_for_face[i];
     }
     sz_nlist_mx = sz_nlist + sz_nlist + 256;
 
     sz = nnode_mx + nface + nnode + 8;
     old2new_nodeid = (int *) malloc(sz * sizeof(int));
     list_default   = old2new_nodeid + nnode_mx;
     inside_ea_node = list_default + nface;
 
     for (i = 0; i < 2; i++) {
         sz = nface_mx + sz_nlist_mx + sz_nlist_mx + nedge_mx + nedge_mx;
         nnode_for_face_p[i]    = (int *) malloc(sz * sizeof(int));
         nodelist_for_face_p[i] = nnode_for_face_p[i]    + nface_mx;
         edgelist_for_face_p[i] = nodelist_for_face_p[i] + sz_nlist_mx;
         nodelist_for_edge_p[i] = edgelist_for_face_p[i] + sz_nlist_mx;
 
         coords1d_p[i] = (double *) malloc(3 * nnode_mx * sizeof(double));
     }
     sz = sz_nlist_mx + sz_nlist_mx + nface_mx + nface_mx + nnode_mx + nedge_mx + nedge_mx + sz_nlist_mx;
 
     nodelist_for_face_out1d = (int *) malloc(sz * sizeof(int));
     nnode_for_face_out1d    = nodelist_for_face_out1d + (sz_nlist_mx + sz_nlist_mx);
     nodelist_intface        = nnode_for_face_out1d + (nface_mx + nface_mx);
     nodelist_for_edge_work  = nodelist_intface     + nnode_mx;
     edgelist_for_face_work  = nodelist_for_edge_work  + (nedge_mx + nedge_mx);
 
 
     coords1d_out2 = (double *) malloc(6 * nnode_mx * sizeof(double));
     xs            = coords1d_out2 +(nnode_mx * 3);
     ys            = xs            + nnode_mx;
     zs            = ys            + nnode_mx;
 
     for (i = 0; i < 2; i++) {
         nodelist_for_face_out2[i] = nodelist_for_face_out1d + (i * sz_nlist_mx);
         nnode_for_face_out2[i]    = nnode_for_face_out1d    + (i * nface_mx);
     }
     allinside = 1;
     for (n = 0; n < nnode; n++) {
         inside_ea_node[n] = 1;
         cs = coords1d + (n * dim);
         for (i = 0; i < dim; i++) {
             if ((cs[i] < xl[i]) || (cs[i] > xr[i])) {
                 allinside = 0;
                 inside_ea_node[n] = 0;
             }
         }
     }
//   check whether the polyhedron is inside the rectangular.
 
     if (allinside) {
         *nnode_out = nnode;
         *nface_out = nface;
 
         memcpy(coords1d_out,          coords1d,          (size_t)(nnode    * szdim));
         memcpy(nnode_for_face_out,    num_node_for_face, (size_t)(nface    * szint));
         memcpy(nodelist_for_face_out, nodelist_for_face, (size_t)(sz_nlist * szint));
 
          free(old2new_nodeid);
 
          for (i = 0; i < 2; i++) {
              free(nnode_for_face_p[i]);
 
              free(coords1d_p[i]);
          }
          free(nodelist_for_face_out1d);
 
          free(coords1d_out2);
 
          free(norm_of_face);
 
         return;
     }
//   check whether all nodes of poly are obviously outside the rectangular.
 
     offset = 0;
     for (k = 0; k < nface; k++) {
         nn = num_node_for_face[k];
         nodelist_s = nodelist_for_face + offset;
         i = 0;
         i1 = (i + 1) % nn;
         i2 = (i1+ 1) % nn;
         n  = nodelist_s[i];
         n1 = nodelist_s[i1];
         n2 = nodelist_s[i2];
         p0 = coords1d + (n  * dim);
         p1 = coords1d + (n1 * dim);
         p2 = coords1d + (n2 * dim);
         for (i = 0; i < dim; i++) {
             dp0[i] = p1[i] - p0[i];
             dp1[i] = p2[i] - p1[i];
         }
         norm = norm_of_face + (3 * k);
         norm[0] = dp0[1] * dp1[2] - dp0[2] * dp1[1];
         norm[1] = dp0[2] * dp1[0] - dp0[0] * dp1[2];
         norm[2] = dp0[0] * dp1[1] - dp0[1] * dp1[0];
         tmp = norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2];
         tmp = 1.0/sqrt(tmp);
         for (i1 = 0; i1 < dim; i1++) {
             norm[i1] *= tmp;
         }
         offset += nn;
     }
     for (i = 0; i < nface; i++) {
         list_default[i] = i;
     }
     cal_ctr0(nface, nnode, list_default, num_node_for_face,
              nodelist_for_face, coords1d, ctr, &vol);
 
     c8[0][0] = xl[0];
     c8[0][1] = xl[1];
     c8[0][2] = xl[2];
 
     c8[1][0] = xr[0];
     c8[1][1] = xl[1];
     c8[1][2] = xl[2];
 
     c8[2][0] = xr[0];
     c8[2][1] = xr[1];
     c8[2][2] = xl[2];
 
     c8[3][0] = xl[0];
     c8[3][1] = xr[1];
     c8[3][2] = xl[2];
 
     c8[4][0] = xl[0];
     c8[4][1] = xl[1];
     c8[4][2] = xr[2];
 
     c8[5][0] = xr[0];
     c8[5][1] = xl[1];
     c8[5][2] = xr[2];
 
     c8[6][0] = xr[0];
     c8[6][1] = xr[1];
     c8[6][2] = xr[2];
 
     c8[7][0] = xl[0];
     c8[7][1] = xr[1];
     c8[7][2] = xr[2];
 
//   check whether all vertices are inside the sphere.
 
     for (i = 0; i < 8; i++) {
         cs  = c8[i];
         ifinside[i] = 1;   // inside
 
         offset = 0;
         for (f = 0; f < nface; f++) {
             norm = norm_of_face + (3 * f);
             nodelist_s = nodelist_for_face + offset;
             n0 = nodelist_s[0];
             p0 = coords1d + (n0 * dim);
             dot_ctr  = 0.0;
             dot_cube = 0.0;
             for (k = 0; k < dim; k++) {
                 dctr  = ctr[k] - p0[k];
                 dcube = cs[k]  - p0[k];
                 dot_ctr  += (norm[k] * dctr);
                 dot_cube += (norm[k] * dcube);
             }
             if (dot_ctr * dot_cube < 0.0) {
                 ifinside[i] = 0;  // outside
                 break;
             }
             else if (dot_ctr * dot_cube == 0.0) {
                 ifinside[i] = -1; // touch
             }
             offset += num_node_for_face[f];
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < 8; i++) {
         if (ifinside[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < 8; i++) {
         if (ifinside[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < 8; i++) {
         if (ifinside[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
         *nnode_out = 8;
         *nface_out = 6;
 
         memcpy(coords1d_out, c8[0], (size_t)(8 * szdim));
         for (i = 0; i < 6; i++) {
             nnode_for_face_out[i] = 4;
         }
         nodelist_t = nodelist_for_face_out;
         nodelist_t[0] = 0;
         nodelist_t[1] = 3;
         nodelist_t[2] = 7;
         nodelist_t[3] = 4;
         nodelist_t   += 4;
 
         nodelist_t[0] = 1;
         nodelist_t[1] = 5;
         nodelist_t[2] = 6;
         nodelist_t[3] = 2;
         nodelist_t   += 4;
 
         nodelist_t[0] = 0;
         nodelist_t[1] = 4;
         nodelist_t[2] = 5;
         nodelist_t[3] = 1;
         nodelist_t   += 4;
 
         nodelist_t[0] = 2;
         nodelist_t[1] = 6;
         nodelist_t[2] = 7;
         nodelist_t[3] = 3;
         nodelist_t   += 4;
 
         nodelist_t[0] = 0;
         nodelist_t[1] = 1;
         nodelist_t[2] = 2;
         nodelist_t[3] = 3;
         nodelist_t   += 4;
 
         nodelist_t[0] = 4;
         nodelist_t[1] = 7;
         nodelist_t[2] = 6;
         nodelist_t[3] = 5;
 
          free(old2new_nodeid);
 
          for (i = 0; i < 2; i++) {
              free(nnode_for_face_p[i]);
 
              free(coords1d_p[i]);
          }
          free(nodelist_for_face_out1d);
 
          free(coords1d_out2);
 
          free(norm_of_face);
 
         return;
     }
     for (i = 0; i < nnode; i++) {
         anyinside += inside_ea_node[i];
     }
     if (ifinquiry && anyinside) {
         *nnode_out = 1;
         *nface_out = 1;
 
          free(old2new_nodeid);
 
          for (i = 0; i < 2; i++) {
              free(nnode_for_face_p[i]);
 
              free(coords1d_p[i]);
          }
          free(nodelist_for_face_out1d);
 
          free(coords1d_out2);
 
          free(norm_of_face);
 
         return;
     }
     if (ifinquiry) {
         *nnode_out = 1;
         *nface_out = 1;
 
         free(old2new_nodeid);
 
          for (i = 0; i < 2; i++) {
              free(nnode_for_face_p[i]);
 
              free(coords1d_p[i]);
          }
          free(nodelist_for_face_out1d);
 
          free(coords1d_out2);
 
          free(norm_of_face);
 
         return;
     }
//   find intersecting
 
     for (i = 0; i < dim; i++) {
         dx[i] = xr[i] - xl[i];
     }
//   find intersecting
 
     copy_read  = 0;
     memcpy(nnode_for_face_p[copy_read],    num_node_for_face, (size_t)(nface     *szint));
     memcpy(nodelist_for_face_p[copy_read], nodelist_for_face, (size_t)(sz_nlist  *szint));
     memcpy(edgelist_for_face_p[copy_read], edgelist_for_face, (size_t)(sz_nlist  *szint));
     memcpy(nodelist_for_edge_p[copy_read], nodelist_for_edge, (size_t)((2*nedge) *szint));
     memcpy(coords1d_p[copy_read],          coords1d,          (size_t)(nnode     *szdim));
 
     nnode_p[copy_read] = nnode;
     nedge_p[copy_read] = nedge;
     nface_p[copy_read] = nface;
 
     for (lr = 0; lr < 2; lr++) {
         if (lr == 0) {
             planes = xl;
             lh = 1;
         }
         else if (lr == 1) {
             planes = xr;
             lh = 0;
         }
         for (idim = 0; idim < dim; idim++) {
             x  = planes[idim];
             copy_write = 1 - copy_read;
 
             my_nnode             = nnode_p[copy_read];
             my_nedge             = nedge_p[copy_read];
             my_nface             = nface_p[copy_read];
             my_coords1d          = coords1d_p[copy_read];
             my_num_node_for_face = nnode_for_face_p[copy_read];
             my_nodelist_for_face = nodelist_for_face_p[copy_read];
             my_edgelist_for_face = edgelist_for_face_p[copy_read];
             my_nodelist_for_edge = nodelist_for_edge_p[copy_read];
 
             poly3d_cut_by_plane(dim, idim, x, xl, dx[0],
                            my_nface, my_nedge, my_nnode, my_coords1d,
                            my_num_node_for_face, my_nodelist_for_face,
                            my_edgelist_for_face, my_nodelist_for_edge,
                            &nnode_out2, coords1d_out2,
                            nface_for_zone_out2, nnode_for_face_out2,
                            nodelist_for_face_out2,
                            &nnode_intface, nodelist_intface);
 
 
 
             nf = nface_for_zone_out2[lh];
             if (nf == 0) {
                 *nnode_out = 0;
                 *nface_out = 0;
 
                 free(old2new_nodeid);
 
                 for (i = 0; i < 2; i++) {
                     free(nnode_for_face_p[i]);
 
                     free(coords1d_p[i]);
                 }
                 free(nodelist_for_face_out1d);
 
                 free(coords1d_out2);
 
                 free(norm_of_face);
 
                 return;
             }
             for (n = 0; n < nnode_out2; n++) {
                 old2new_nodeid[n] = -1;
             }
    //       fine the net number of nodes
 
             nodelist_s       = nodelist_for_face_out2[lh];
             nnode_for_face_s = nnode_for_face_out2[lh];
 
             nnode_this_p = 0;
             for (f = 0; f < nf; f++) {
                 nn = nnode_for_face_s[f];
                 for (j = 0; j < nn; j++) {
                     n = nodelist_s[j];
                     if (old2new_nodeid[n] < 0) {
                         old2new_nodeid[n] = nnode_this_p;
                         nnode_this_p++;
                     }
                 }
                 nodelist_s += nn;
             }
    //       find nedge, num_edge_for_face, edgelist_for_face, and nodelist_for_edge
 
             nodelist_t = nodelist_for_edge_work;
             edgelist_t = edgelist_for_face_work;
             nodelist_s = nodelist_for_face_out2[lh];
 
             nedge_this_p = 0;
 
             for (f = 0; f < nf; f++) {
                 nn = nnode_for_face_s[f];
                 for (j = 0; j < nn; j++) {
                     j1 = (j + 1) % nn;
                     n0 = nodelist_s[j];
                     n1 = nodelist_s[j1];
 
    //               check whether this edge is already exist
 
                     e = -1;
                     for (k = 0; k < nedge_this_p; k++) {
                         nodelist_exist = nodelist_for_edge_work + (k + k);
                         m0 = nodelist_exist[0];
                         m1 = nodelist_exist[1];
                         if (((m0 == n0)&&(m1 == n1)) || ((m0 == n1)&&(m1 == n0))) {
                             e = k;
                             break;
                         }
                     }
                     if (e >= 0) {
                         edgelist_t[j] = e;
                     }
                     else {
                         edgelist_t[j] = nedge_this_p;
                         nodelist_t[0] = n0;
                         nodelist_t[1] = n1;
                         nodelist_t   += 2;
                         nedge_this_p++;
                     }
                 }
                 nodelist_s += nn;
                 edgelist_t += nn;
             }
             nnode_p[copy_write] = nnode_this_p;
             nedge_p[copy_write] = nedge_this_p;
             nface_p[copy_write] = nf;
 
             memcpy(nnode_for_face_p[copy_write], nnode_for_face_out2[lh], (size_t)(nf * szint));
 
             nnode_for_face_s = nnode_for_face_out2[lh];
             sz_nlist = 0;
             for (f = 0; f < nf; f++) {
                 sz_nlist += nnode_for_face_s[f];
             }
             edgelist_s = edgelist_for_face_work;
             edgelist_t = edgelist_for_face_p[copy_write];
             memcpy(edgelist_t, edgelist_s, (size_t)(sz_nlist * szint));
 
             nodelist_t = nodelist_for_face_p[copy_write];
             nodelist_s = nodelist_for_face_out2[lh];
 
             for (j = 0; j < sz_nlist; j++) {
                 n0 = nodelist_s[j];
                 n1 = old2new_nodeid[n0];
                 assert(n1 >= 0);
                 nodelist_t[j] = n1;
             }
             nodelist_s = nodelist_for_edge_work;
             nodelist_t = nodelist_for_edge_p[copy_write];
 
             for (j = 0; j < nedge_this_p + nedge_this_p; j++) {
                 n0 = nodelist_s[j];
                 n1 = old2new_nodeid[n0];
                 assert(n1 >= 0);
                 nodelist_t[j] = n1;
             }
             coord1ds_t = coords1d_p[copy_write];
             for (j = 0; j < nnode_out2; j++) {
                 n  = old2new_nodeid[j];
                 if (n >= 0) {
                     cs = coords1d_out2 + (dim * j);
                     ct = coord1ds_t    + (dim * n);
                     memcpy(ct, cs, (size_t)szdim);
                 }
             }
/******
//////////// FOR DEBUG   ///////////////////////////////////////////////////////////////
             {
                char name[64];
                char name1[] = "poly";
                char name2[] = "rect";
                char nameo[] = "poly_int";
                int  fileid;
                int  default_list[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
                double xarray[128], yarray[128], zarray[128];
                mio_Coord *coord;
                mio_Unstructured_Mesh mesh;
 
                sprintf(name, "poly_cut");
 
                mio_create_file(name, mio_file_create, mio_independent, &fileid);
                mio_init(fileid, mio_umesh, -1, &mesh);
                coord = &(mesh.coord);
                coord->datatype = mio_double;
                mesh.name       = name1;
                mesh.dims       = 3;
                mesh.idmin      = 0;
                mesh.type       = mio_general_mesh;
                mesh.datatype   = mio_int;
                coord->datatype = mio_double;
 
                mesh.sizes[0]           = 1;
                mesh.sizes[1]           = nface;
                mesh.sizes[3]           = nnode;
                mesh.num_faces_for_zone = &nface;
                mesh.facelist_for_zone  = default_list;
                mesh.num_nodes_for_face = num_node_for_face;
                mesh.nodelist_for_face  = nodelist_for_face;
 
                for (j = 0; j < nnode; j++) {
                    cs = coords1d + (dim * j);
                    xarray[j] = cs[0];
                    yarray[j] = cs[1];
                    zarray[j] = cs[2];
                }
                coord->coord[0] = (void *)xarray;
                coord->coord[1] = (void *)yarray;
                coord->coord[2] = (void *)zarray;
 
                mio_write(fileid, mio_umesh, fileid, &mesh);
 
                mio_init(fileid, mio_umesh, -1, &mesh);
                coord = &(mesh.coord);
                coord->datatype = mio_double;
                mesh.name       = name2;
                mesh.dims       = 3;
                mesh.idmin      = 0;
                mesh.type       = mio_hex;
                mesh.datatype   = mio_int;
                coord->datatype = mio_double;
 
                mesh.sizes[0]           = 1;
                mesh.sizes[3]           = 8;
                mesh.nodelist_for_zone  = default_list;
 
                for (j = 0; j < 8; j++) {
                    xarray[j] = c8[j][0];
                    yarray[j] = c8[j][1];
                    zarray[j] = c8[j][2];
                }
                coord->coord[0] = (void *)xarray;
                coord->coord[1] = (void *)yarray;
                coord->coord[2] = (void *)zarray;
 
                mio_write(fileid, mio_umesh, fileid, &mesh);
 
                mio_init(fileid, mio_umesh, -1, &mesh);
                coord = &(mesh.coord);
                coord->datatype = mio_double;
                mesh.name       = nameo;
                mesh.dims       = 3;
                mesh.idmin      = 0;
                mesh.type       = mio_general_mesh;
                mesh.datatype   = mio_int;
                coord->datatype = mio_double;
 
                mesh.sizes[0]           = 1;
                mesh.sizes[1]           = nface_p[copy_write];
                mesh.sizes[3]           = nnode_p[copy_write];
                mesh.num_faces_for_zone = &(nface_p[copy_write]);
                mesh.facelist_for_zone  = default_list;
                mesh.num_nodes_for_face = nnode_for_face_p[copy_write];
                mesh.nodelist_for_face  = nodelist_for_face_p[copy_write];
 
                for (j = 0; j < nnode_p[copy_write]; j++) {
                    cs = coords1d_p[copy_write] + (dim * j);
                    xarray[j] = cs[0];
                    yarray[j] = cs[1];
                    zarray[j] = cs[2];
                }
                coord->coord[0] = (void *)xarray;
                coord->coord[1] = (void *)yarray;
                coord->coord[2] = (void *)zarray;
 
                mio_write(fileid, mio_umesh, fileid, &mesh);
 
                mio_close_file(fileid);
             }
******/
///////////////////////////////////////////////////////////////////////////////////////////`
 
             copy_read = 1 - copy_read;
         }   // idim
     }       // lr
     *nnode_out = nnode_p[copy_write];
     memcpy(coords1d_out, coords1d_p[copy_write], (size_t)((*nnode_out) * szdim));
     *nface_out = nface_p[copy_write];
     sz_nlist = 0;
     for (f = 0; f < *nface_out; f++) {
         sz_nlist += nnode_for_face_p[copy_write][f];
     }
     memcpy(nnode_for_face_out, nnode_for_face_p[copy_write], (size_t)((*nface_out) * szint));
     memcpy(nodelist_for_face_out, nodelist_for_face_p[copy_write], (size_t)(sz_nlist * szint));
 
/***********
//////  FOR  DEBUG    //////////////////////////////////////////////////////////////////////
     {
        char name[64];
        char name1[] = "poly";
        char name2[] = "rect";
        char nameo[] = "poly_int";
        int  fileid;
        int  default_list[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
        double xarray[128], yarray[128], zarray[128];
        mio_Coord *coord;
        mio_Unstructured_Mesh mesh;
 
        sprintf(name, "poly_cut");
 
        mio_create_file(name, mio_file_create, mio_independent, &fileid);
        mio_init(fileid, mio_umesh, -1, &mesh);
        coord = &(mesh.coord);
        coord->datatype = mio_double;
        mesh.name       = name1;
        mesh.dims       = 3;
        mesh.idmin      = 0;
        mesh.type       = mio_general_mesh;
        mesh.datatype   = mio_int;
        coord->datatype = mio_double;
 
        mesh.sizes[0]           = 1;
        mesh.sizes[1]           = nface;
        mesh.sizes[3]           = nnode;
        mesh.num_faces_for_zone = &nface;
        mesh.facelist_for_zone  = default_list;
        mesh.num_nodes_for_face = num_node_for_face;
        mesh.nodelist_for_face  = nodelist_for_face;
 
        for (j = 0; j < nnode; j++) {
            cs = coords1d + (dim * j);
            xarray[j] = cs[0];
            yarray[j] = cs[1];
            zarray[j] = cs[2];
        }
        coord->coord[0] = (void *)xarray;
        coord->coord[1] = (void *)yarray;
        coord->coord[2] = (void *)zarray;
 
        mio_write(fileid, mio_umesh, fileid, &mesh);
 
        mio_init(fileid, mio_umesh, -1, &mesh);
        coord = &(mesh.coord);
        coord->datatype = mio_double;
        mesh.name       = name2;
        mesh.dims       = 3;
        mesh.idmin      = 0;
        mesh.type       = mio_hex;
        mesh.datatype   = mio_int;
        coord->datatype = mio_double;
 
        mesh.sizes[0]           = 1;
        mesh.sizes[3]           = 8;
        mesh.nodelist_for_zone  = default_list;
 
        for (j = 0; j < 8; j++) {
            xarray[j] = c8[j][0];
            yarray[j] = c8[j][1];
            zarray[j] = c8[j][2];
        }
        coord->coord[0] = (void *)xarray;
        coord->coord[1] = (void *)yarray;
        coord->coord[2] = (void *)zarray;
 
        mio_write(fileid, mio_umesh, fileid, &mesh);
 
        mio_init(fileid, mio_umesh, -1, &mesh);
        coord = &(mesh.coord);
        coord->datatype = mio_double;
        mesh.name       = nameo;
        mesh.dims       = 3;
        mesh.idmin      = 0;
        mesh.type       = mio_general_mesh;
        mesh.datatype   = mio_int;
        coord->datatype = mio_double;
 
        mesh.sizes[0]           = 1;
        mesh.sizes[1]           = *nface_out;
        mesh.sizes[3]           = *nnode_out;
        mesh.num_faces_for_zone = nface_out;
        mesh.facelist_for_zone  = default_list;
        mesh.num_nodes_for_face = nnode_for_face_out;
        mesh.nodelist_for_face  = nodelist_for_face_out;
 
        for (j = 0; j < *nnode_out; j++) {
            cs = coords1d_out + (dim * j);
            xarray[j] = cs[0];
            yarray[j] = cs[1];
            zarray[j] = cs[2];
        }
        coord->coord[0] = (void *)xarray;
        coord->coord[1] = (void *)yarray;
        coord->coord[2] = (void *)zarray;
 
        mio_write(fileid, mio_umesh, fileid, &mesh);
 
        mio_close_file(fileid);
     }
*********/
///////////////////////////////////////////////////////////////////////////////////////////
 
     free(old2new_nodeid);
 
     for (i = 0; i < 2; i++) {
         free(nnode_for_face_p[i]);
 
         free(coords1d_p[i]);
     }
     free(nodelist_for_face_out1d);
 
     free(coords1d_out2);
 
     free(norm_of_face);
 
     return;
 }
 
 
 
void poly3d_cut_by_plane(int dim, int norm_direction, double plane_location,
                 double *x0_scale, double dx_scale,
                 int nface, int nedge0, int nnode,  double *coords1d,
                 int *num_node_for_face, int *nodelist_for_face,
                 int *edgelist_for_face0, int *nodelist_for_edge0,
                 int *nnode_out, double *coords1d_out,
                 int *nface_for_zone_out, int **nnode_for_face_out,
                 int **nodelist_for_face_out,
                 int *nnode_intface, int *nodelist_intface)
{
//   If nedge = 0 in input, there are no edges involed and therefore
//   edgelist_for_face and nodelist_for_edge are not vlaid.
 
     int  szdim;
     int  i, k, i1, f, n0, n1, nn, nf, e, m0, m1, offset;
     int  sz_nlist, sz_nlist_mx, nedge_mx, nnode_mx;
 
     int  *nodelist_s, *nodelist_t, *edgelist_t;
 
     double dxinv, xplane, xmax, xmin;
     double c3[3], *ct, *cs;
     double *my_coords1d;
 
     int mynnode[64], mynodes[64], include_face[64];
     int new2old_face[64], old2new_face[64];
     int offset_new, nnode_new, nf_new, nn_new, n, idim, part;
     int *nodes, *nodes_new, *nodes_old;
     int *nnode_ea_face, *nodelist_ea_face;
     int *include_node, *old2new_node, *new2old_node;
     double small, small2, dx, ds2;
 
     int nedge, *edgelist_for_face, *nodelist_for_edge;
 
     small = 1.0e-10;
 
     szdim = dim * sizeof(double);
     dxinv = 1.0/dx_scale;
 
     xmax = dxinv * coords1d[norm_direction];
     xmin = xmax;
 
     for (i = 1; i < nnode; i++) {
         ct = coords1d + (dim * i + norm_direction);
         dx = (*ct) * dxinv;
         if (xmax < dx) {
             xmax = dx;
         }
         else if (xmin > dx) {
             xmin = dx;
         }
     }
     xplane = plane_location * dxinv;
     if (xmin + small >= xplane) {
         n0 = 1;
     }
     else if (xmax - small <= xplane) {
         n0 = 0;
     }
     else {
         n0 = -1;
     }
     sz_nlist = 0;
     for (i = 0; i < nface; i++) {
         sz_nlist += num_node_for_face[i];
     }
     if (n0 >= 0) {
         *nnode_out = nnode;
         memcpy(coords1d_out, coords1d, (size_t)(nnode * szdim));
         nface_for_zone_out[n0] = nface;
         nface_for_zone_out[1-n0] = 0;
         memcpy(nnode_for_face_out[n0], num_node_for_face,
               (size_t)(nface * sizeof(int)));
         memcpy(nodelist_for_face_out[n0], nodelist_for_face,
               (size_t)(sz_nlist * sizeof(int)));
         *nnode_intface = 0;
 
         return;
     }
     sz_nlist_mx = sz_nlist + 256;
     nedge_mx    = nnode + nnode + 64;
     nnode_mx    = nnode + nnode + 64;
 
     my_coords1d = (double *) malloc(dim * nnode_mx * sizeof(double));
 
     if ((nedge == 0) && (nface != 0)) {
 
//       find nedge, num_edge_for_face, edgelist_for_face, and nodelist_for_edge,
 
         get_edges(nface, nnode, num_node_for_face,
                   nodelist_for_face, &nedge,
                   &edgelist_for_face, &nodelist_for_edge);
     }
     else {
         nedge = nedge0;
         edgelist_for_face = edgelist_for_face0;
         nodelist_for_edge = nodelist_for_edge0;
     }
     if (norm_direction == 0) {
         for (i = 0; i < nnode; i++) {
             ct = my_coords1d + (dim * i);
             cs = coords1d    + (dim * i);
             ct[0] = dxinv *(cs[0] - plane_location);
             ct[1] = dxinv *(cs[1] - x0_scale[1]);
             ct[2] = dxinv *(cs[2] - x0_scale[2]);
         }
     }
     else if (norm_direction == 1) {
         for (i = 0; i < nnode; i++) {
             ct = my_coords1d + (dim * i);
             cs = coords1d    + (dim * i);
             ct[0] = dxinv *(cs[1] - plane_location);
             ct[1] = dxinv *(cs[2] - x0_scale[2]);
             ct[2] = dxinv *(cs[0] - x0_scale[0]);
         }
     }
     else if (norm_direction == 2) {
         for (i = 0; i < nnode; i++) {
             ct = my_coords1d + (dim * i);
             cs = coords1d    + (dim * i);
             ct[0] = dxinv *(cs[2] - plane_location);
             ct[1] = dxinv *(cs[0] - x0_scale[0]);
             ct[2] = dxinv *(cs[1] - x0_scale[1]);
         }
     }
     xplane = 0.0;
 
     poly3d_cut_by_xplane(dim, xplane, nface, nedge, nnode, my_coords1d,
                          num_node_for_face, edgelist_for_face,
                          nodelist_for_edge, nodelist_for_face,
                          nnode_out, coords1d_out,
                          nface_for_zone_out, nnode_for_face_out,
                          nodelist_for_face_out,
                          nnode_intface, nodelist_intface);
 
//////////////////////////////////////////////////////////////////////
//   remove possible redundant nodes
 
     small2 = small * small;
     old2new_node = (int *) malloc(3 * (*nnode_out) * sizeof(int));
     new2old_node = old2new_node + (*nnode_out);
     include_node = new2old_node + (*nnode_out);
 
     for (i = 0; i < *nnode_out; i++) {
         old2new_node[i] = -1;
     }
     nnode_new       = 1;
     old2new_node[0] = 0;
     new2old_node[0] = 0;
     include_node[0] = 1;
 
     for (i = 1; i < *nnode_out; i++) {
         ct = coords1d_out + (i * dim);
         for (k = 0; k < nnode_new; k++) {
             n  = new2old_node[k];
             cs = coords1d_out + (n * dim);
             ds2 = 0.0;
             for (idim = 0; idim < dim; idim++) {
                 dx = ct[idim] - cs[idim];
                 ds2 += (dx * dx);
             }
             if (ds2 < small2) {
                 old2new_node[i] = k;
                 include_node[i] = 0;
                 break;
             }
         }
         if (old2new_node[i] < 0) {
             include_node[i] = 1;
             old2new_node[i] = nnode_new;
             new2old_node[nnode_new] = i;
 
             nnode_new++;
         }
     }
     if (nnode_new < nnode) {
 
         for (n = 0; n < nnode_new; n++) {
             i  = new2old_node[n];
             ct = coords1d_out + (n * dim);
             cs = coords1d_out + (i * dim);
             memcpy(ct, cs, (size_t)(dim * sizeof(double)));
         }
         for (part = 0; part < 2; part++) {
             nf    = nface_for_zone_out[part];
             nnode_ea_face    = nnode_for_face_out[part];
             nodelist_ea_face = nodelist_for_face_out[part];
 
             offset = 0;
             nf_new = 0;
             for (f = 0; f < nf; f++) {
                 nn = nnode_ea_face[f];
                 nodes  = nodelist_ea_face + offset;
                 nn_new = 0;
                 for (i = 0; i < nn; i++) {
                     n = nodes[i];
                     if (include_node[n]) {
                         mynodes[nn_new] = old2new_node[n];
                         nn_new++;
                     }
                 }
                 memcpy(nodes, mynodes, (size_t)(nn_new * sizeof(int)));
                 mynnode[f] = nn_new;
                 if (nn_new < 3) {
                     include_face[f] = 0;
                 }
                 else {
                     include_face[f] = 1;
                     old2new_face[f] = nf_new;
                     new2old_face[nf_new] = f;
 
                     nf_new++;
                 }
                 offset += nn;
             }
//           consolidate nodes
 
             for (f = 0; f < nf; f++) {
                 nn = nnode_ea_face[f];
                 nn_new = mynnode[f];
                 offset     = 0;
                 offset_new = 0;
                 nodes_old = nodelist_ea_face + offset;
                 nodes_new = nodelist_ea_face + offset_new;
                 for (i = 0; i < nn_new; i++) {
                     nodes_new[i] = nodes_old[i];
                 }
                 nnode_ea_face[f] = nn_new;
                 offset += nn;
                 offset_new += nn_new;
             }
//           consolidate faces
 
             if (nf_new != nf) {
                 if (nf_new > 3) {
                     nf_new = 0;
                     offset = 0;
                     offset_new = 0;
                     for (f = 0; f < nf; f++) {
                         nn = nnode_ea_face[f];
                         if (include_face[f]) {
                             nodes_old = nodelist_ea_face + offset;
                             nodes_new = nodelist_ea_face + offset_new;
                             memcpy(nodes_new, nodes_old,
                                   (size_t)(nn * sizeof(int)));
                             nnode_ea_face[nf_new] = nn;
                             offset_new += nn;
                             nf_new++;
                         }
                         offset += nn;
                     }
                 }
                 else {
                     nface_for_zone_out[part] = 0;
                 }
             }
         }
         free(old2new_node);
     }
     else {
         free(old2new_node);
     }
/////////////////////////////////////////////////////////////////////////
 
     if (norm_direction == 0) {
         for (i = 0; i < *nnode_out; i++) {
             ct = coords1d_out + (dim * i);
             memcpy(c3, ct, (size_t)szdim);
             ct[0] = dx_scale *c3[0] + plane_location;
             ct[1] = dx_scale *c3[1] + x0_scale[1];
             ct[2] = dx_scale *c3[2] + x0_scale[2];
         }
     }
     else if (norm_direction == 1) {
         for (i = 0; i < *nnode_out; i++) {
             ct = coords1d_out + (dim * i);
             memcpy(c3, ct, (size_t)szdim);
             ct[0] = dx_scale *c3[2] + x0_scale[0];
             ct[1] = dx_scale *c3[0] + plane_location;
             ct[2] = dx_scale *c3[1] + x0_scale[2];
         }
     }
     else if (norm_direction == 2) {
         for (i = 0; i < *nnode_out; i++) {
             ct = coords1d_out + (i * dim);
             memcpy(c3, ct, (size_t)szdim);
             ct[0] = dx_scale *c3[1] + x0_scale[0];
             ct[1] = dx_scale *c3[2] + x0_scale[1];
             ct[2] = dx_scale *c3[0] + plane_location;
         }
     }
     if ((nedge0 == 0) && (nface != 0)) {
         free(edgelist_for_face);
         free(nodelist_for_edge);
     }
     free(my_coords1d);
 
     return;
 }
 
 
void poly3d_cut_by_xplane(int dim, double xplane,
                 int nface, int nedge0, int nnode,  double *coords1d,
                 int *num_node_for_face, int *edgelist_for_face0,
                 int *nodelist_for_edge0, int *nodelist_for_face,
                 int *nnode_out, double *coords1d_out,
                 int *nface_for_zone_out, int **nnode_for_face_out,
                 int **nodelist_for_face_out,
                 int *nnode_intface, int *nodelist_intface)
{
//   The outputs of the intersection between a plane x = xplane and polyhedron
//   nface, nedge, nnode, coords,num_node_for_face, edgelist_for_face, and
//   nodelist_for_edge are input for the polyhedron.
//
//   nnode_out, coords_out, nface_for_zone_out, nnode_for_face_out, and
//   nnode_for_face_out are output for two resulting polyhedrons.
//
//   nnode_out includes nnode and the new nodes involved in the two resulting polys.
//
//   Each edge is assumed to be unique in the input.
 
     int   szdim2, szdim, szint, failed;
     int   myside, found, edge_on_interface_done, need_edge_oninterface;
     int   i, i1, j, k, k1, ii, jj, idx, nlast, e, ee, f, np, np_added, np0, ip;
     int   n0, n, n1, nn, m0, m1, ne, nf, ne_tot;
 
     int   sznlist, sznlist_mx, offset, offset_edge, offset_face;
 
     int   ne_added[2], nf2[2];
     int   *mynodes;
     int   nn2[2], ne2[2], ne2_tot[2], nface_ea_zone2[2], sides[2];
     int   *nodelist2[2];
     int   *node_included, *edge_included;
     int   *offplane_ea_node, *side_ea_edge;
     int   *nn_onplane_ea_edge;
     int   *intersectd_ea_edge;
     int   *edge2intfnode;
     int   *nodeid_added2np;  //  old2new_edgeid[2*mxnum][2];
     int   *old2new_edgeid;
 
     int   *nodelist_for_face2[2];
     int   *nodelist_this_edge2[2];
     int   *nnode_for_face2[2], *nedge_ea_face2[2];
     int   *edgelist_ea_face2[2], *edgelist2[2];
     int   *nodelist_ea_edge2[2];
     int   *offsets_edgelist_for_face;
     int   *nodelist, *edgelist, *facelist, *edgelist_t, *nodelist_t;
     int   *nodelist_this_edge;
 
     double  small, eta;
     double  *c, *ci;
 
     double  *c2[2];
     double  *coords_int2d, *coords_int2d_work;
 
     int nedge, *edgelist_for_face, *nodelist_for_edge;
 
//   determin whether each node is on the plane
//   offplane_ea_node[i] =  0  :  node i is on the plane
//                       = -1  :  x[i] < xplane
//                       =  1  :  x[i] > xplane
//
     small  = 1.0e-12;
 
     szint  = sizeof(int);
     szdim2 = sizeof(double);
     szdim  = 3 * szdim2;
     szdim2 = szdim2 + szdim2;
 
     sznlist = 0;
     for (i = 0; i < nface; i++) {
         sznlist += num_node_for_face[i];
     }
     sznlist_mx = sznlist + sznlist + 256;
////////////////////////////////////////////////////////////////////////////////////////////
 
//   build edges if neccessary
 
     if ((nedge0 == 0) && (nface != 0)) {
 
//       find nedge, num_node_for_face, edgelist_for_face, and nodelist_for_edge,
 
         get_edges(nface, nnode, num_node_for_face,
                   nodelist_for_face, &nedge,
                   &edgelist_for_face, &nodelist_for_edge);
     }
     else {
         nedge = nedge0;
         edgelist_for_face = edgelist_for_face0;
         nodelist_for_edge = nodelist_for_edge0;
     }
//////////////////////////////////////////////////////////////////////////////
     mynodes            = (int *) malloc(nnode * sizeof(int));
     node_included      = (int *) malloc(nnode * sizeof(int));
     edge_included      = (int *) malloc(nedge * sizeof(int));
     offplane_ea_node   = (int *) malloc(nnode * sizeof(int));
     side_ea_edge       = (int *) malloc(nedge * sizeof(int));
     intersectd_ea_edge = (int *) malloc(nedge * sizeof(int));
     edge2intfnode      = (int *) malloc(nedge * sizeof(int));
 
     nodeid_added2np    = (int *) malloc(4*nnode * sizeof(int));
     old2new_edgeid     = (int *) malloc(4*nedge * sizeof(int));
 
     offsets_edgelist_for_face = (int *) malloc(nface * sizeof(int));
     coords_int2d = (double *) malloc(4*nnode * sizeof(double));
     coords_int2d_work = (double *) malloc(4*nnode * sizeof(double));
 
     nn_onplane_ea_edge = (int *) malloc(sznlist_mx * sizeof(int));
 
     for (i = 0; i < 2; i++) {
         nodelist_for_face2[i] = (int *) malloc(sznlist_mx * sizeof(int));
         nodelist_this_edge2[i] = (int *) malloc(sznlist_mx * sizeof(int));
         edgelist_ea_face2[i]   = (int *) malloc(sznlist_mx * sizeof(int));
         nodelist_ea_edge2[i]   = (int *) malloc(sznlist_mx * sizeof(int));
 
         nnode_for_face2[i]     = (int *) malloc(nface * sizeof(int));
         nedge_ea_face2[i]      = (int *) malloc(nface * sizeof(int));
     }
//////////////////////////////////////////////////////////////////////////////
 
     offset = 0;
     for (i = 0; i < nface; i++) {
         offsets_edgelist_for_face[i] = offset;
         offset += num_node_for_face[i];
     }
     myside = 0;
     for (i = 0; i < nnode; i++) {
         c = coords1d + (dim * i);
         offplane_ea_node[i] = 0;
         if (fabs(c[0] - xplane) < small) {
             offplane_ea_node[i] = 0;
         }
         else if (c[0] < xplane) {
             offplane_ea_node[i] = -1;
             myside = -1;
         }
         else {
             offplane_ea_node[i] = 1;
             myside = 1;
         }
     }
     nlast = offset;
//   check whether all nodes are on the same side of the plane
 
     found = 0;
     for (i = 0; i < nnode; i++) {
         if (offplane_ea_node[i] * myside < 0) {
             found = 1;
             break;
         }
     }
     if (!found) {
         idx = (1 + myside)/2;
         *nnode_out = nnode;
         memcpy(coords1d_out, coords1d, (size_t)(nnode * szdim));
         nface_for_zone_out[idx] = nface;
         nface_for_zone_out[1-idx] = 0;
         memcpy(nnode_for_face_out[idx], num_node_for_face,    (size_t)(nface * szint));
         memcpy(nodelist_for_face_out[idx], nodelist_for_face, (size_t)(nlast * szint));
         *nnode_intface = 0;
 
         if ((nedge0 == 0) && (nface != 0)) {
             free(edgelist_for_face);
             free(nodelist_for_edge);
         }
         return;
     }
 
//   determine whether each edge is cross the plane, excluding the touch
 
//   side_ea_edge[i] =  2: both nodes of the edge, x > xplane
//                   = -2: both nodes of the edge, x < xplane
//                   =  1: one node of the edge x > xplane, the other on the plane
//                   = -1: one node of the edge x < xplane, the other on the plane
//                   =  0: both nodes on the plane, or one node on one side of the plane
 
//   nn_onplane_ea_edge[i] =  2: both nodes are on the plane.
//   nn_onplane_ea_edge[i] =  1: one node is on the plane and the other x > xplane
//   nn_onplane_ea_edge[i] = -1: one node is on the plane and the other x < xplane
//   nn_onplane_ea_edge[i] =  0: no node on the plane
 
     for (i = 0; i < nedge; i++) {
         side_ea_edge[i] = 0;
 
         nodelist = nodelist_for_edge + (i + i);
         nn_onplane_ea_edge[i] = 0;
 
         n0 = nodelist[0];
         n  = nodelist[1];
 
         if ((offplane_ea_node[n0] < 0) && (offplane_ea_node[n] < 0)) {
             side_ea_edge[i]       = -2;
             nn_onplane_ea_edge[i] = 0;
         }
         else if ((offplane_ea_node[n0] > 0) && (offplane_ea_node[n] > 0)) {
             side_ea_edge[i]       = 2;
             nn_onplane_ea_edge[i] = 0;
         }
         else if (offplane_ea_node[n0] * offplane_ea_node[n] < 0) {
             side_ea_edge[i]       = 0;
             nn_onplane_ea_edge[i] = 0;
         }
         else if ((offplane_ea_node[n0] == 0) && (offplane_ea_node[n] == 0)) {
             side_ea_edge[i]       = 0;
             nn_onplane_ea_edge[i] = 2;
         }
         else if (offplane_ea_node[n0] == 0) {
             if (offplane_ea_node[n] < 0) {
                 side_ea_edge[i] = -1;
                 nn_onplane_ea_edge[i] = -1;
             }
             else if (offplane_ea_node[n] > 0) {
                 side_ea_edge[i] = 1;
                 nn_onplane_ea_edge[i] = 1;
             }
         }
         else if (offplane_ea_node[n] == 0) {
             if (offplane_ea_node[n0] < 0) {
                 side_ea_edge[i] = -1;
                 nn_onplane_ea_edge[i] = -1;
             }
             else if (offplane_ea_node[n0] > 0) {
                 side_ea_edge[i] = 1;
                 nn_onplane_ea_edge[i] = 1;
             }
         }
     }
//   calculate the intersectios between each edge and the plane
 
     for (i = 0; i < nnode; i++) {
         node_included[i] = 0;
     }
     np_added = 0;
     np       = 0;
     for (i = 0; i < nedge; i++) {
         intersectd_ea_edge[i] =  0;
         edge2intfnode[i]      = -1;
 
         nodelist = nodelist_for_edge + (i + i);
         if (nn_onplane_ea_edge[i] == 2) {
             intersectd_ea_edge[i] = 1;
             for (k = 0; k < 2; k++) {
                 n = nodelist[k];
                 if (!node_included[n]) {
                     memcpy(coords_int2d + (np + np), coords1d + (dim * n + 1), (size_t)szdim2);
                     nodelist_intface[np] = n;
                     np++;
                     node_included[n] = 1;
                 }
             }
         }
         else if ((nn_onplane_ea_edge[i] ==  1) ||
                  (nn_onplane_ea_edge[i] == -1)) {
             intersectd_ea_edge[i] = 1;
             for (k = 0; k < 2; k++) {
                 n = nodelist[k];
                 if ((offplane_ea_node[n] == 0) && !node_included[n]) {
                     memcpy(coords_int2d + (np + np), coords1d + (dim * n + 1), (size_t)szdim2);
                     nodelist_intface[np] = n;
                     np++;
                     node_included[n] = 1;
                 }
             }
         }
         else if ((nn_onplane_ea_edge[i] == 0)&&(side_ea_edge[i] == 0)) {
 
//           calculate the intersection
 
             for (k = 0; k < 2; k++) {
                 n  = nodelist[k];
                 c2[k] = coords1d + (dim * n);
                 sides[k] = offplane_ea_node[n];
             }
             assert(sides[0] * sides[1] != 0);
             if (sides[0] * sides[1] < 0) {
 
                 eta = (xplane - c2[0][0])/(c2[1][0] - c2[0][0]);
                 eta = MAX(0.0, MIN(1.0, eta));
                 c = coords_int2d + (np + np);
                 c[0] = c2[0][1] + eta *(c2[1][1] - c2[0][1]);
                 c[1] = c2[0][2] + eta *(c2[1][2] - c2[0][2]);
                 intersectd_ea_edge[i] = 1;
                 nodelist_intface[np]  = nnode + np_added;
                 edge2intfnode[i]      = nnode + np_added;
                 nodeid_added2np[np_added] = np;
 
                 np++;
                 np_added++;
             }
         }
     }
     n = 0;
     assert((np == 0) || (np >= 3));
 
//   put np nodes into the correct order
 
     *nnode_intface = 0;
     if (np >= 3) {
         memcpy(coords_int2d_work, coords_int2d,(size_t)(np * szdim2));
         mychull_sort(coords_int2d_work, nodelist_intface, np, &failed);
         *nnode_intface = np;
     }
///////////////////////////////////////////////////////////////////////
//   debug -- view the interface
 
 
 
//////////////////////////////////////////////////////////////////////
//   copy nodes to results
 
     memcpy(coords1d_out, coords1d, (size_t)(nnode * szdim));
     *nnode_out = nnode + np_added;
 
     for (i = 0; i < np_added; i++) {
         ip = nodeid_added2np[i];
         ci = coords_int2d + (ip + ip);
         c  = coords1d_out + (dim * (nnode + i));
         c[0] = xplane;
         c++;
         for (j = 0; j < 2; j++) {
             c[j] = ci[j];
         }
     }
//   assemble the polyhedron on each side of the plane
 
     for (idx = 0; idx < 2; idx++) {
         nface_ea_zone2[idx] = 0;
         ne2_tot[idx]        = 0;
         edgelist2[idx] = edgelist_ea_face2[idx];
     }
     for (i = 0; i < nedge + nedge; i++) {
         old2new_edgeid[i] = -1;
     }
//     for (i = 0; i < nedge; i++) {
//         old2new_edgeid[i][0] = -1;
//         old2new_edgeid[i][1] = -1;
//     }
     for (f = 0; f < nface; f++) {
 
         need_edge_oninterface  = 0;
         edge_on_interface_done = 0;
 
         nf2[0] = nface_ea_zone2[0];
         nf2[1] = nface_ea_zone2[1];
         nedge_ea_face2[0][nf2[0]] = 0;
         nedge_ea_face2[1][nf2[1]] = 0;
 
         ne_added[0] = 0;
         ne_added[1] = 0;
 
         ne2[0] = 0;             // for face on x < 0 side
         ne2[1] = 0;             // for face on x > 0 side
         nodelist2[0] = nodelist_this_edge2[0];
         nodelist2[1] = nodelist_this_edge2[1];
 
         ne       = num_node_for_face[f];
         offset   = offsets_edgelist_for_face[f];
         edgelist = edgelist_for_face + offset;
 
         for (i = 0; i < ne; i++) {
             e  = edgelist[i];
             nodelist = nodelist_for_edge + (e + e);
 
             if (side_ea_edge[e] == -2) {
                 if (old2new_edgeid[e+e] < 0) {
                     for (j = 0; j < 2; j++) {
                         n = nodelist[j];
                         nodelist2[0][j] = n;
                     }
                     *(edgelist2[0]) = ne2_tot[0] + ne2[0];
                     old2new_edgeid[e+e] = *(edgelist2[0]);
 
                     nodelist2[0] += 2;
                     ne2[0]++;
                 }
                 else {
                     *(edgelist2[0]) = old2new_edgeid[e+e];
                 }
                 (edgelist2[0])++;
                 (nedge_ea_face2[0][nf2[0]])++;
                 ne_added[0]++;
             }
             else if (side_ea_edge[e] == 2) {
                 if (old2new_edgeid[e+e+1] < 0) {
                     for (j = 0; j < 2; j++) {
                         n = nodelist[j];
                         nodelist2[1][j] = n;
                     }
                     *(edgelist2[1])= ne2_tot[1] + ne2[1];
                     old2new_edgeid[e+e+1] = *(edgelist2[1]);
 
                     nodelist2[1] += 2;
                     ne2[1]++;
                 }
                 else {
                     *(edgelist2[1]) = old2new_edgeid[e+e+1];
                 }
                 (edgelist2[1])++;
                 (nedge_ea_face2[1][nf2[1]])++;
                 ne_added[1]++;
             }
             else if ((side_ea_edge[e] == -1)||(side_ea_edge[e] == 1))  { // one node on the plane
                 idx = -1;
                 for (j = 0; j < 2; j++) {
                     n = nodelist[j];
                     if (offplane_ea_node[n] < 0) {
                         idx = 0;
                         break;
                     }
                     else if (offplane_ea_node[n] > 0) {
                         idx = 1;
                         break;
                     }
                 }
                 assert(idx >= 0);
                 if (old2new_edgeid[e+e+idx] < 0) {
                     for (j = 0; j < 2; j++) {
                         n = nodelist[j];
                         nodelist2[idx][j] = n;
                     }
                     *(edgelist2[idx]) = ne2_tot[idx] + ne2[idx];
                     old2new_edgeid[e+e+idx] = *(edgelist2[idx]);
 
                     nodelist2[idx] += 2;
                     ne2[idx]++;
                 }
                 else {
                     *(edgelist2[idx]) = old2new_edgeid[e+e+idx];
                 }
                 (edgelist2[idx])++;
                 (nedge_ea_face2[idx][nf2[idx]])++;
                 ne_added[idx]++;
 
                 need_edge_oninterface = 1;
             }
             else if (nn_onplane_ea_edge[e] == 2) {  // both nodes on the plane
 
                 for (idx = 0; idx < 2; idx++) {
                     if (old2new_edgeid[e+e+idx] < 0) {
                         for (j = 0; j < 2; j++) {
                             n = nodelist[j];
                             nodelist2[idx][j] = n;
                         }
                         *(edgelist2[idx])= ne2_tot[idx] + ne2[idx];
                         old2new_edgeid[e+e+idx] = *(edgelist2[idx]);
 
                         nodelist2[idx] += 2;
                         ne2[idx]++;
                     }
                     else {
                         *(edgelist2[idx]) = old2new_edgeid[e+e+idx];
                     }
                     (edgelist2[idx])++;
                     (nedge_ea_face2[idx][nf2[idx]])++;
                     ne_added[idx]++;
                 }
             }
             else if ((nn_onplane_ea_edge[e]==0)&&(side_ea_edge[e]==0)) {   // the edge intersects with the plane
                 idx = -1;
                 ip  = edge2intfnode[e];
                 assert(ip >= nnode);
 
                 for (j = 0; j < 2; j++) {
                     n = nodelist[j];
                     if (offplane_ea_node[n] < 0) {
                         idx = 0;
                     }
                     else if (offplane_ea_node[n] > 0) {
                         idx = 1;
                     }
                     assert(idx >= 0);
                     if (old2new_edgeid[e+e+idx] < 0) {
                         nodelist2[idx][0] = n;
                         nodelist2[idx][1] = ip;
                         *(edgelist2[idx])= ne2_tot[idx] + ne2[idx];
                         old2new_edgeid[e+e+idx] = *(edgelist2[idx]);
 
                         nodelist2[idx] += 2;
                         ne2[idx]++;
                     }
                     else {
                         *(edgelist2[idx]) = old2new_edgeid[e+e+idx];
                     }
                     (edgelist2[idx])++;
                     (nedge_ea_face2[idx][nf2[idx]])++;
                     ne_added[idx]++;
                 }
                 // find the edge on the interface
 
                 if (!edge_on_interface_done) {
                     for (k1 = 1; k1 < ne; k1++) {
                         k = (ne + i - k1) % ne;
                         ee  = edgelist[k];
                         if (ee == e) continue;
 
                         nodelist_this_edge = nodelist_for_edge + (ee + ee);
                         if (side_ea_edge[ee] == -2) continue;
                         if (side_ea_edge[ee] ==  2) continue;
 
                         if ((side_ea_edge[ee] == 1) || (side_ea_edge[ee] == -1)) {
                             // one node on the plane
                             for (j = 0; j < 2; j++) {
                                 n  = nodelist_this_edge[j];
                                 n0 = nodelist_this_edge[1-j];
                                 idx = -1;
                                 if (offplane_ea_node[n] < 0) {
                                     idx = 0;
                                 }
                                 else if (offplane_ea_node[n] > 0) {
                                     idx = 1;
                                 }
                                 if (idx >= 0) {
                                     nodelist2[idx][0] = ip;
                                     nodelist2[idx][1] = n0;
                                     *(edgelist2[idx])= ne2_tot[idx] + ne2[idx];
                                     nodelist2[idx] += 2;
                                     ne2[idx]++;
                                     (edgelist2[idx])++;
                                     (nedge_ea_face2[idx][nf2[idx]])++;
                                     ne_added[idx]++;
 
                                     idx = 1 - idx;
                                     nodelist2[idx][0] = ip;
                                     nodelist2[idx][1] = n0;
                                     *(edgelist2[idx])= ne2_tot[idx] + ne2[idx];
                                     nodelist2[idx] += 2;
                                     ne2[idx]++;
                                     (edgelist2[idx])++;
                                     (nedge_ea_face2[idx][nf2[idx]])++;
                                     ne_added[idx]++;
 
                                     edge_on_interface_done = 1;
                                     break;
 
                                 }
                             }
                             if (edge_on_interface_done) {
                                 break;
                             }
//                           edge_on_interface_done = 1;
//                           break; for the possible edge of the face of the other idx
                         }
                         else if ((side_ea_edge[ee] == 0)&&(nn_onplane_ea_edge[ee] == 0)) {
                             n0 = edge2intfnode[ee];
 
                             assert(n0 >= nnode);
 
                             for (idx = 0; idx < 2; idx++) {
                                 nodelist2[idx][0] = ip;
                                 nodelist2[idx][1] = n0;
                                 *(edgelist2[idx])= ne2_tot[idx] + ne2[idx];
                                 nodelist2[idx] += 2;
                                 ne2[idx]++;
                                 (edgelist2[idx])++;
                                 (nedge_ea_face2[idx][nf2[idx]])++;
                                 ne_added[idx]++;
                             }
                             edge_on_interface_done = 1;
                             break;
                         }
                     }
                 }
             }
         }
         if (!edge_on_interface_done) {
             for (i = 0; i < ne; i++) {
                 e  = edgelist[i];
                 if (nn_onplane_ea_edge[e] == 2) {  // both nodes on the plane
                     need_edge_oninterface = 0;
                     break;
                 }
             }
             if (need_edge_oninterface) {
//               find two nodes on the face that are on the interface
 
                 nodelist = nodelist_for_face + offsets_edgelist_for_face[f];
                 nn = 0;
                 for (i = 0; i < ne; i++) {
                     n  = nodelist[i];
                     if (offplane_ea_node[n] == 0) {
                         mynodes[nn] = n;
                         nn++;
                     }
                 }
                 if (nn == 2) {
                     for (idx = 0; idx < 2; idx++) {
                         nodelist2[idx][0] = mynodes[0];
                         nodelist2[idx][1] = mynodes[1];
                         *(edgelist2[idx])= ne2_tot[idx] + ne2[idx];
                         nodelist2[idx] += 2;
                         ne2[idx]++;
                         (edgelist2[idx])++;
                         (nedge_ea_face2[idx][nf2[idx]])++;
                         ne_added[idx]++;
                     }
                 }
                 else {
                    assert(nn < 2);
                 }
             }
         }
         for (idx = 0; idx < 2; idx++) {
             if (ne2[idx] > 0) {
                 ne     = ne2[idx];
                 ne_tot = ne2_tot[idx];
                 nf     = nface_ea_zone2[idx];
 
                 nodelist_t = nodelist_ea_edge2[idx] + (ne_tot + ne_tot);
                 memcpy(nodelist_t, nodelist_this_edge2[idx], (size_t)((ne + ne) * szint));
                 ne2_tot[idx] += ne;
             }
             if (ne_added[idx] >= 3) {
                 nface_ea_zone2[idx]++;
             }
             else {
                 edgelist2[idx] -= ne_added[idx];
             }
         }
     }             // f
 
//   remove the edgelist in the result
 
     for (idx = 0; idx < 2; idx++) {
         nf   = nface_ea_zone2[idx];
         if (nf < 2) continue;
 
         nodelist_t = nodelist_for_face2[idx];
 
         offset = 0;
         for (f = 0; f < nf; f++) {
             ne = nedge_ea_face2[idx][f];
             for (k = 0; k < ne; k++) {
                 edge_included[k] = 0;
             }
             edgelist = edgelist_ea_face2[idx] + offset;
             e = edgelist[0];
             nodelist_this_edge = nodelist_ea_edge2[idx] + (e + e);
             n0    = nodelist_this_edge[0];
             nlast = nodelist_this_edge[1];
 
             nodelist_t[0] = n0;
             nodelist_t[1] = nlast;
 
             nodelist_t += 2;
             nn = 2;
             for (i = 1; i < ne-1; i++) {
                 found = 0;
                 for (k = 1; k < ne; k++) {
                     if (edge_included[k]) continue;
                     e = edgelist[k];
                     nodelist_this_edge = nodelist_ea_edge2[idx] + (e + e);
                     for (j = 0; j < 2; j++) {
                         n  = nodelist_this_edge[j];
                         if (n == nlast) {
                             found = 1;
                             nlast = nodelist_this_edge[1-j];
                             nodelist_t[0] = nlast;
                             nodelist_t++;
                             nn++;
                             edge_included[k] = 1;
                             break;
                         }
                     }
                     if (found) break;
                 }
                 assert(found);
             }
             assert(nn > 2);
             nnode_for_face2[idx][f] = nn;
 
             offset += ne;
         }
     }
 
//   Put the interface into the two zones
 
     for (idx = 0; idx < 2; idx++) {
         nf   = nface_ea_zone2[idx];
         nface_for_zone_out[idx] = 0;
         if (nf < 2) continue;
 
         nodelist_t    = nodelist_for_face_out[idx];
 
         nface_for_zone_out[idx] = 1 + nf;
         nnode_for_face_out[idx][0] = np;
         memcpy(nodelist_t,  nodelist_intface, (size_t)(np * szint));
         nodelist_t += np;
 
         nodelist = nodelist_for_face2[idx];
         for (k = 0; k < nf; k++) {
             f = k + 1;
             nn = nedge_ea_face2[idx][k];
             nnode_for_face_out[idx][f] = nn;
             memcpy(nodelist_t, nodelist, (size_t)(nn * szint));
             nodelist_t += nn;
             nodelist   += nn;
         }
     }
 
     free(mynodes);
     free(node_included);
     free(edge_included);
     free(offplane_ea_node);
     free(side_ea_edge);
     free(intersectd_ea_edge);
     free(edge2intfnode);
 
     free(nodeid_added2np);
     free(old2new_edgeid);
     free(offsets_edgelist_for_face);
 
     free(coords_int2d);
     free(coords_int2d_work);
 
     free(nn_onplane_ea_edge);
     for (i = 0; i < 2; i++) {
         free(nodelist_for_face2[i]);
         free(nodelist_this_edge2[i]);
         free(edgelist_ea_face2[i]);
         free(nodelist_ea_edge2[i]);
         free(nnode_for_face2[i]);
         free(nedge_ea_face2[i]);
     }
     if ((nedge0 == 0) && (nface != 0)) {
         free(edgelist_for_face);
         free(nodelist_for_edge);
     }
 
     return;
 }
 
 
void grevolve_obj(geom_obj_type obj,
                  int ifinquiry, int *mixed,
                  double *pt_axis, int axis_rotate,
                  double angle_theta, double angle_phi,
                  double angle_begin, double angle_end,
            double *ctr, double radius,                   // for obj = 2, sphere
            int if_helix, double h0,   double dhdphi,     // for helix
            int nnode, double *xys, double *xys_scratch,  // for obj = 1, poly
            double *x2min, double *x2max,                 // for obj = 1, poly
            double *xl, double *dx, double *vol)
{
//   if axis_rotate = 0:  rotate around x-axis
//                  = 1:  y-axis
//                  = 2:  z-axis
//                  other : rotate around (angle_theta, angle_phi).
//
//   orig(0:1) is the reference point of vertices of the polygon, coords2d.
 
     int    list_default[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
     int    i, k, dim, nn, ne, nf, split, mymixed;
     double myxl[3], mydx[3], myctr[3];
     double factor, vol0, theta, phi, sinp, cosp, sint, cost;
     double xr[3], c8[8][3], cp8[8][3], *c, *cp;
     double mypt_axis[3];
 
     int    nnode8, nedge8, nface8;
     int    nface, nploy, nnode_out8[8], nedge_out8[8], nface_out8[8];
     int    nnode_for_face[6], nodelist_for_face[24];
     int    nodelist_for_edge[24], edgelist_for_face[24];
     int    **nnode_for_face_out8, **nodelist_for_face_out8;
     int    **edgelist_for_face_out8, **nodelist_for_edge_out8;
     int    *mynnode_for_face, *mynodelist_for_face;
     int    *myedgelist_for_face, *mynodelist_for_edge;
     double val, dvol, dvol0, **coords_out8;
 
     if (axis_rotate == 2) {
         grevolve_obj_z(obj,
                         ifinquiry, mixed, angle_begin, angle_end, pt_axis,
                         ctr, radius, if_helix, h0, dhdphi,
                         nnode, xys, xys_scratch, x2min, x2max,
                         xl, dx, vol);
     }
     else if (axis_rotate == 1) {
         myxl[0] = xl[2];
         myxl[1] = xl[0];
         myxl[2] = xl[1];
 
         mydx[0] = dx[2];
         mydx[1] = dx[0];
         mydx[2] = dx[1];
 
         mypt_axis[0] = pt_axis[2];
         mypt_axis[1] = pt_axis[0];
         mypt_axis[2] = pt_axis[1];
 
         grevolve_obj_z(obj,
                         ifinquiry, mixed, angle_begin, angle_end, mypt_axis,
                         ctr, radius, if_helix, h0, dhdphi,
                         nnode, xys, xys_scratch, x2min, x2max,
                         myxl, mydx, vol);
     }
     else if (axis_rotate == 0) {
         myxl[0] = xl[1];
         myxl[1] = xl[2];
         myxl[2] = xl[0];
 
         mydx[0] = dx[1];
         mydx[1] = dx[2];
         mydx[2] = dx[0];
 
         mypt_axis[0] = pt_axis[1];
         mypt_axis[1] = pt_axis[2];
         mypt_axis[2] = pt_axis[0];
 
         grevolve_obj_z(obj,
                        ifinquiry, mixed, angle_begin, angle_end, mypt_axis,
                        ctr, radius, if_helix, h0, dhdphi,
                        nnode, xys, xys_scratch, x2min, x2max,
                        myxl, mydx, vol);
     }
     else {
         factor = PI / 180.0;
         theta  = angle_theta * factor;
         phi    = angle_phi   * factor;
         sinp   = sin(phi);
         cosp   = sqrt(MAX(1.0 - sinp * sinp, 0.0));
         sint   = sin(theta);
         cost   = sqrt(MAX(1.0 - sint * sint, 0.0));
 
         sinp   = - sinp;
         sint   = - sint;
 
         dim = 3;
         vol0 = 1.0;
         for (i = 0; i < dim; i++) {
             xr[i] = xl[i] + dx[i];
             vol0 *= dx[i];
         }
         c8[0][0] = xl[0];
         c8[0][1] = xl[1];
         c8[0][2] = xl[2];
 
         c8[1][0] = xr[0];
         c8[1][1] = xl[1];
         c8[1][2] = xl[2];
 
         c8[2][0] = xr[0];
         c8[2][1] = xr[1];
         c8[2][2] = xl[2];
 
         c8[3][0] = xl[0];
         c8[3][1] = xr[1];
         c8[3][2] = xl[2];
 
         c8[4][0] = xl[0];
         c8[4][1] = xl[1];
         c8[4][2] = xr[2];
 
         c8[5][0] = xr[0];
         c8[5][1] = xl[1];
         c8[5][2] = xr[2];
 
         c8[6][0] = xr[0];
         c8[6][1] = xr[1];
         c8[6][2] = xr[2];
 
         c8[7][0] = xl[0];
         c8[7][1] = xr[1];
         c8[7][2] = xr[2];
 
         nnode8 = 8;
         for (i = 0; i < nnode8; i++) {
             c  = c8[i];
             cp = cp8[i];
             cp[0] =  cost * cosp * c[0] + cost * sinp * c[1] - sint * c[2];
             cp[1] = -       sinp * c[0] +        cosp * c[1];
             cp[2] =  sint * cosp * c[0] + sint * sinp * c[1] + cost * c[2];
         }
//       This polyhedron should be cut planes of axes, as done in grevolve_obj_z, and
//       then call grevolve_obj_p. This is to be done.
 
         nface8 = 6;
         for (i = 0; i < nface8; i++) {
             nnode_for_face[i] = 4;
         }
         nodelist_for_face[0] =  0;
         nodelist_for_face[1] =  3;
         nodelist_for_face[2] =  7;
         nodelist_for_face[3] =  4;
 
         nodelist_for_face[4] =  1;
         nodelist_for_face[5] =  5;
         nodelist_for_face[6] =  6;
         nodelist_for_face[7] =  2;
 
         nodelist_for_face[8]  =  0;
         nodelist_for_face[9]  =  4;
         nodelist_for_face[10] =  5;
         nodelist_for_face[11] =  1;
 
         nodelist_for_face[12] =  2;
         nodelist_for_face[13] =  6;
         nodelist_for_face[14] =  7;
         nodelist_for_face[15] =  3;
 
         nodelist_for_face[16] =  0;
         nodelist_for_face[17] =  1;
         nodelist_for_face[18] =  2;
         nodelist_for_face[19] =  3;
 
         nodelist_for_face[20] =  4;
         nodelist_for_face[21] =  7;
         nodelist_for_face[22] =  6;
         nodelist_for_face[23] =  5;
 
         nedge8 = 12;
         nodelist_for_edge[0]  = 0;
         nodelist_for_edge[1]  = 1;
         nodelist_for_edge[2]  = 1;
         nodelist_for_edge[3]  = 2;
         nodelist_for_edge[4]  = 2;
         nodelist_for_edge[5]  = 3;
         nodelist_for_edge[6]  = 3;
         nodelist_for_edge[7]  = 0;
 
         nodelist_for_edge[8]  = 4;
         nodelist_for_edge[9]  = 5;
         nodelist_for_edge[10] = 5;
         nodelist_for_edge[11] = 6;
         nodelist_for_edge[12] = 6;
         nodelist_for_edge[13] = 7;
         nodelist_for_edge[14] = 7;
         nodelist_for_edge[15] = 4;
 
         nodelist_for_edge[16] = 0;
         nodelist_for_edge[17] = 4;
         nodelist_for_edge[18] = 1;
         nodelist_for_edge[19] = 5;
         nodelist_for_edge[20] = 2;
         nodelist_for_edge[21] = 6;
         nodelist_for_edge[22] = 3;
         nodelist_for_edge[23] = 7;
 
         edgelist_for_face[0]  = 3;
         edgelist_for_face[1]  = 11;
         edgelist_for_face[2]  = 7;
         edgelist_for_face[3]  = 8;
 
         edgelist_for_face[4]  = 1;
         edgelist_for_face[5]  = 9;
         edgelist_for_face[6]  = 5;
         edgelist_for_face[7]  = 10;
 
         edgelist_for_face[8]  = 0;
         edgelist_for_face[9]  = 8;
         edgelist_for_face[10] = 4;
         edgelist_for_face[11] = 9;
 
         edgelist_for_face[12] = 2;
         edgelist_for_face[13] = 10;
         edgelist_for_face[14] = 6;
         edgelist_for_face[15] = 11;
 
         edgelist_for_face[16] = 0;
         edgelist_for_face[17] = 1;
         edgelist_for_face[18] = 2;
         edgelist_for_face[19] = 3;
 
         edgelist_for_face[20] = 4;
         edgelist_for_face[21] = 7;
         edgelist_for_face[22] = 6;
         edgelist_for_face[23] = 5;
 
//       The runs with and without split seem getting the same results for
//       oblique rotation, and both incorrect.
 
         split = 0;
         for (i = 0; i < dim; i++) {
             val = cp8[0][i];
             for (k = 0; k < 8; k++) {
                 cp = cp8[k];
                 if (val * cp[i] < 0.0) {
                     split = 1;
                     break;
                 }
             }
             if (split) break;
         }
         if (!split) {
             grevolve_obj_p(obj,
                         ifinquiry, mixed, angle_begin, angle_end,
                         ctr, radius, if_helix, h0, dhdphi,
                         nnode, xys, x2min, x2max,
                         cp8[0], nnode8, nedge8, nface8,
                         nnode_for_face, nodelist_for_face,
                         nodelist_for_edge, edgelist_for_face,
                         dx[0], vol0, vol);
         }
         else {
             coords_out8 = NULL;
             nnode_for_face_out8 = NULL;
             nodelist_for_face_out8 = NULL;
             nodelist_for_edge_out8 = NULL;
             edgelist_for_face_out8 = NULL;
 
             axis_decompose(cp8[0], nnode8, nedge8, nface8,
                   nnode_for_face, nodelist_for_face,
                   nodelist_for_edge, edgelist_for_face,
                   dx[0],
                   &nploy, nnode_out8, nedge_out8, nface_out8,
                   &coords_out8, &nnode_for_face_out8,
                   &nodelist_for_face_out8,
                   &nodelist_for_edge_out8, &edgelist_for_face_out8);
 
             if (nploy == 1) {
                 grevolve_obj_p(obj,
                         ifinquiry, mixed, angle_begin, angle_end,
                         ctr, radius, if_helix, h0, dhdphi,
                         nnode, xys, x2min, x2max,
                         cp8[0], nnode8, nedge8, nface8,
                         nnode_for_face, nodelist_for_face,
                         nodelist_for_edge, edgelist_for_face,
                         dx[0], vol0, vol);
             }
             else {
                 *mixed = 0;
                 *vol = 0.0;
                 for (i = 0; i < nploy; i++) {
                     cp = coords_out8[i];
                     nn = nnode_out8[i];
                     ne = nedge_out8[i];
                     nf = nface_out8[i];
                     mynnode_for_face    = nnode_for_face_out8[i];
                     mynodelist_for_face = nodelist_for_face_out8[i];
                     mynodelist_for_edge = nodelist_for_edge_out8[i];
                     myedgelist_for_face = edgelist_for_face_out8[i];
 
                     cal_ctr0(nf, nn, list_default, mynnode_for_face,
                              mynodelist_for_face, cp, myctr, &dvol0);
 
                     grevolve_obj_p(obj,
                             ifinquiry, &mymixed, angle_begin, angle_end,
                             ctr, radius, if_helix, h0, dhdphi,
                             nnode, xys, x2min, x2max,
                             cp, nn, ne, nf,
                             mynnode_for_face, mynodelist_for_face,
                             mynodelist_for_edge, myedgelist_for_face,
                             dx[0], dvol0, &dvol);
 
                     *vol += dvol;
                     if (mymixed) {
                         *mixed = 1;
                     }
                     free(mynnode_for_face);
                     free(mynodelist_for_face);
                     free(mynodelist_for_edge);
                     free(myedgelist_for_face);
                 }
                 if (!(*mixed) && (*vol > 0.0)) {
                     if (fabs(*vol - vol0) > small) {
                         *mixed = 1;
                     }
                 }
                 free(nnode_for_face_out8);
                 free(nodelist_for_face_out8);
                 free(nodelist_for_edge_out8);
                 free(edgelist_for_face_out8);
             }
         }
     }
     return;
 }
 
void axis_decompose(double *coords, int nnode, int nedge, int nface,
               int *nnode_for_face, int *nodelist_for_face,
               int *nodelist_for_edge, int *edgelist_for_face,
               double dx_scale, int *npoly_out,
               int *nnode_out, int *nedge_out, int *nface_out,
               double ***coords_out, int ***nnode_for_face_out,
               int ***nodelist_for_face_out,
               int ***nodelist_for_edge_out, int ***edgelist_for_face_out)
{
//   npoly_out[8], nnode_out[8], nedge_out[8], and nface_out[8]
//   are assumed to be allocated before calling this function.
//
//   If *npoly_out = 1 in output,
//       nnode_out[0] = nnode
//       nface_out[0] = nface
//       coords_out, nface_for_poly, nnode_for_face_out, and
//       nodelist_for_face_out are not allocated in the routine.
//   If *npoly_out > 1 in output, free allocateion after use
//       for (i = 0; i < npoly_out; i++) {
//           free(coords_out[i]);
//           free(nnode_for_face_out[i]);
//           free(nodelist_for_face_out[i]);
//           free(nodelist_for_edge_out[i]);
//           free(edgelist_for_face_out[i]);
//        }
//        free(coords_out);
//        free(nnode_for_face_out);
//        free(nodelist_for_face_out);
//        free(nodelist_for_edge_out);
//        free(edgelist_for_face_out);
//
     const int mxnum  = 32;
 
     double zeros[] = {0.0, 0.0, 0.0};
 
     int  dim, idim, i, i1, j, k, f, lh, idx, szint, szdim, offset;
     int  split, split3[3];
     int  r, w, sznl, p, ip, np, nf, ne, nn, mynf, myne, mynn;
     int  e, n, n0, n1, m0, m1;
     int  *nn_for_face, *nlist_for_face;
     int  *nodelist_s, *nodelist_t, *edgelist_t, *nodelist_exist;
     int  *elist_for_face, *nlist_for_edge;
     int  included[mxnum], new2old[mxnum], old2new[mxnum];
     int  mynodelist_for_face[mxnum*mxnum], mynodelist_for_edge[mxnum+mxnum];
     int  myedgelist_for_face[mxnum*mxnum];
 
     int     npoly_p[2], nnode_p[2][8], nedge_p[2][8], nface_p[2][8];
     int     nnode_for_face_p[2][8][mxnum];             // [rw][p][f]
     int     nodelist_for_face_p[2][8][mxnum*mxnum];    // [rw][p][f]
     int     edgelist_for_face_p[2][8][mxnum*mxnum];    // [rw][p][f]
     int     nodelist_for_edge_p[2][8][mxnum*mxnum];    // [rw][p][f]
     double  coords1d_p[2][8][3*mxnum];                 // [rw][p][n]
 
     int     nnode_out2;
     int     nface_for_zone_out2[2];
     double  coords1d_out2[3*(mxnum+mxnum)], coords1d_out[3*mxnum];
 
     int     nodelist_for_face_out1d[2*mxnum*mxnum], *nodelist_for_face_out2[2];
     int     nnode_for_face_out1d[mxnum+mxnum],      *nnode_for_face_out2[2];
 
     int     nnode_intface, nodelist_intface[mxnum];
 
     double  x, *cs, *p0;
 
     dim   = 3;
     szdim = dim * sizeof(double);
     szint = sizeof(int);
 
//   check whether each of axis-plane cut the polygons
 
     for (k = 0; k < dim; k++) {
         split3[k] = 0;
         x = coords[k];
         for (i = 0; i < nnode; i++) {
             p0 = coords + (dim * i);
             if (p0[k] * x < 0.0) {
                 split3[k] = 1;
                 break;
             }
         }
     }
     split = split3[0];
     for (k = 1; k < dim; k++) {
         split += split3[k];
     }
     if (!split) {
         *npoly_out = 1;
         nnode_out[0] = nnode;
         nedge_out[0] = nedge;
         nface_out[0] = nface;
 
         *coords_out            = NULL;
         *nnode_for_face_out    = NULL;
         *nodelist_for_face_out = NULL;
         *nodelist_for_edge_out = NULL;
         *edgelist_for_face_out = NULL;
 
         return;
     }
     for (i = 0; i < 2; i++) {
         nodelist_for_face_out2[i] = nodelist_for_face_out1d + (i*mxnum*mxnum);
         nnode_for_face_out2[i]    = nnode_for_face_out1d    + (i*mxnum);
     }
     sznl = 0;
     for (i = 0; i < nface; i++) {
         sznl += nnode_for_face[i];
     }
     r = 0;
     memcpy(nnode_for_face_p[r][0],    nnode_for_face,    (size_t)(nface *szint));
     memcpy(nodelist_for_face_p[r][0], nodelist_for_face, (size_t)(sznl  *szint));
     memcpy(nodelist_for_edge_p[r][0], nodelist_for_edge, (size_t)((nedge+nedge)*szint));
     memcpy(edgelist_for_face_p[r][0], edgelist_for_face, (size_t)(sznl  *szint));
 
     memcpy(coords1d_p[r][0],          coords,            (size_t)(nnode *szdim));
 
     npoly_p[r]    = 1;
     nnode_p[r][0] = nnode;
     nedge_p[r][0] = nedge;
     nface_p[r][0] = nface;
 
     x = 0.0;
     for (idim = 0; idim < dim; idim++) {
         if (!split3[idim]) continue;
         w = 1 - r;
         npoly_p[w] = 0;
 
         np = npoly_p[r];
 
         for (p = 0; p < np; p++) {
             nf = nface_p[r][p];
             ne = nedge_p[r][p];
             nn = nnode_p[r][p];
             cs = coords1d_p[r][p];
             nn_for_face    = nnode_for_face_p[   r][p];
             nlist_for_face = nodelist_for_face_p[r][p];
             elist_for_face = edgelist_for_face_p[r][p];
             nlist_for_edge = nodelist_for_edge_p[r][p];
 
             poly3d_cut_by_plane( dim, idim, x,
                             zeros, dx_scale, nf, ne, nn, cs,
                             nn_for_face, nlist_for_face,
                             elist_for_face, nlist_for_edge,
                             &nnode_out2, coords1d_out2,
                             nface_for_zone_out2, nnode_for_face_out2,
                             nodelist_for_face_out2,
                             &nnode_intface, nodelist_intface);
 
             assert(nnode_out2 < mxnum);
 
             for (lh = 0; lh < 2; lh++) {
                 if (!nface_for_zone_out2[lh]) continue;
 
                 mynn           = nnode_out2;
                 cs             = coords1d_out2;
                 mynf           = nface_for_zone_out2[lh];
                 nn_for_face    = nnode_for_face_out2[lh];
                 nlist_for_face = nodelist_for_face_out2[lh];
                 sznl = 0;
                 for (i = 0; i < mynf; i++) {
                     sznl += nn_for_face[i];
                 }
                 for (i = 0; i < nnode_out2; i++) {
                     included[i] = 0;
                 }
                 offset = 0;
                 for (f = 0; f < mynf; f++) {
                     nn = nn_for_face[f];
                     nodelist_s = nlist_for_face + offset;
                     for (i = 0; i < nn; i++) {
                         n = nodelist_s[i];
                         included[n] = 1;
                     }
                     offset += nn;
                 }
                 mynn = 0;
                 for (i = 0; i < nnode_out2; i++) {
                     if (included[i]) {
                         new2old[mynn] = i;
                         old2new[i]    = mynn;
                         memcpy(coords1d_out + (dim * mynn), coords1d_out2 + (dim * i), (size_t)szdim);
                         mynn++;
                     }
                     else {
                         old2new[i] = -1;
                     }
                 }
                 offset = 0;
                 for (f = 0; f < mynf; f++) {
                     nn = nn_for_face[f];
                     nodelist_s = nlist_for_face + offset;
                     nodelist_t = mynodelist_for_face + offset;
                     for (i = 0; i < nn; i++) {
                         n = nodelist_s[i];
                         idx = old2new[n];
                         nodelist_t[i] = idx;
                     }
                     offset += nn;
                 }
//               find edgelist_for_face and nodelist_for_edge
 
 
                 nodelist_s = mynodelist_for_face;
                 nodelist_t = mynodelist_for_edge;
                 edgelist_t = myedgelist_for_face;
 
                 myne = 0;
                 for (f = 0; f < mynf; f++) {
                     nn = nn_for_face[f];
                     for (i = 0; i < nn; i++) {
                         i1 = (i + 1) % nn;
                         n0 = nodelist_s[i];
                         n1 = nodelist_s[i1];
 
        //               check whether this edge is already exist
 
                         e = -1;
                         for (k = 0; k < myne; k++) {
                             nodelist_exist = mynodelist_for_edge + (k + k);
                             m0 = nodelist_exist[0];
                             m1 = nodelist_exist[1];
                             if (((m0 == n0)&&(m1 == n1)) || ((m0 == n1)&&(m1 == n0))) {
                                 e = k;
                                 break;
                             }
                         }
                         if (e >= 0) {
                             edgelist_t[i] = e;
                         }
                         else {
                             edgelist_t[i] = myne;
                             nodelist_t[0] = n0;
                             nodelist_t[1] = n1;
                             nodelist_t   += 2;
                             myne++;
                        }
                     }
                     nodelist_s += nn;
                     edgelist_t += nn;
                 }
 
                 ip = npoly_p[w];
                 memcpy(coords1d_p[w][ip],          coords1d_out, (size_t)(mynn * szdim));
                 memcpy(nnode_for_face_p[w][ip],    nn_for_face,  (size_t)(mynf * szint));
                 memcpy(nodelist_for_face_p[w][ip], mynodelist_for_face, (size_t)(sznl * szint));
                 memcpy(edgelist_for_face_p[w][ip], myedgelist_for_face, (size_t)(sznl * szint));
                 memcpy(nodelist_for_edge_p[w][ip], mynodelist_for_edge, (size_t)((myne+myne)*szint));
 
                 nnode_p[w][ip] = mynn;
                 nedge_p[w][ip] = myne;
                 nface_p[w][ip] = mynf;
                 npoly_p[w]++;
             }
         }
         r = w;
     }
     *npoly_out = npoly_p[r];
     if (*npoly_out == 1) {
         nnode_out[0] = nnode;
         nedge_out[0] = nedge;
         nface_out[0] = nface;
 
         *coords_out            = NULL;
         *nnode_for_face_out    = NULL;
         *nodelist_for_face_out = NULL;
         *nodelist_for_edge_out = NULL;
         *edgelist_for_face_out = NULL;
 
         return;
     }
     *coords_out            = (double **) malloc(*npoly_out * sizeof(double *));
     *nnode_for_face_out    = (int    **) malloc(*npoly_out * sizeof(int *));
     *nodelist_for_face_out = (int    **) malloc(*npoly_out * sizeof(int *));
 
     *nodelist_for_edge_out = (int    **) malloc(*npoly_out * sizeof(int *));
     *edgelist_for_face_out = (int    **) malloc(*npoly_out * sizeof(int *));
 
     for (p = 0; p < *npoly_out; p++) {
         nf = nface_p[r][p];
         ne = nedge_p[r][p];
         nn = nnode_p[r][p];
 
         cs = coords1d_p[r][p];
         nn_for_face    = nnode_for_face_p[   r][p];
         nlist_for_face = nodelist_for_face_p[r][p];
         elist_for_face = edgelist_for_face_p[r][p];
         nlist_for_edge = nodelist_for_edge_p[r][p];
 
         sznl = 0;
         for (i = 0; i < nf; i++) {
             sznl += nn_for_face[i];
         }
         (*coords_out)[p]            = (double *) malloc(nn   * szdim);
         (*nnode_for_face_out)[p]    = (int    *) malloc(nf   * szint);
         (*nodelist_for_face_out)[p] = (int    *) malloc(sznl * szint);
         (*edgelist_for_face_out)[p] = (int    *) malloc(sznl * szint);
         (*nodelist_for_edge_out)[p] = (int    *) malloc((ne+ne)*szint);
 
         memcpy((*coords_out)[p],            cs,             (size_t)(nn   * szdim));
         memcpy((*nnode_for_face_out)[p],    nn_for_face,    (size_t)(nf   * szint));
         memcpy((*nodelist_for_face_out)[p], nlist_for_face, (size_t)(sznl * szint));
 
         memcpy((*edgelist_for_face_out)[p], elist_for_face, (size_t)(sznl * szint));
         memcpy((*nodelist_for_edge_out)[p], nlist_for_edge, (size_t)((ne+ne)*szint));
 
         nface_out[p] = nf;
         nedge_out[p] = ne;
         nnode_out[p] = nn;
     }
     return;
 }
 
void grevolve_obj_z(geom_obj_type  obj,
              int ifinquiry, int *mixed,
              double angle_begin, double angle_end, double *pt_axis,
              double *ctr, double radius,
              int if_helix, double h0, double dhdphi,
              int nnode, double *xys, double *xys_scratch,
              double *x2min, double *x2max,
              double *xl, double *dx, double *vol)
{
     int    dim, szdim, i, j, k, n, idim;
     int    ncell_dim[3];
     double vol0, myctr[3], xr[3], myxl[3], myxr[3], xc[3], c8[8][3];
     double xxl[2][3], xxr[2][3];
     double my_x2min[2], my_x2max[2]; 
     double *p, tmp, a_begin, a_end, dvol, myh0, h0_mod, dhdphi_mod;
 
     int    nnode8, nedge8, nface8;
     int    nnode_for_face[6], nodelist_for_face[24];
     int    nodelist_for_edge[24], edgelist_for_face[24];
 
     dim = 3;
     szdim = dim * sizeof(double);
 
     vol0 = 1.0;
     for (i = 0; i < dim; i++) {
         xr[i] = xl[i] + dx[i];
         vol0 *= dx[i];
     }
     *vol = 0.0;
     memcpy(xxl[0], xl, (size_t)szdim);
     memcpy(xxr[0], xr, (size_t)szdim);
 
     if (if_helix) {
         myh0 = h0 - pt_axis[2];
     }
     for (i = 0; i < dim; i++) {
         if ((xl[i] < pt_axis[i]) && (pt_axis[i] < xr[i])) {
             memcpy(xxl[1], xxl[0], (size_t)szdim);
             memcpy(xxr[1], xxr[0], (size_t)szdim);
             xxr[0][i] = pt_axis[i];
             xxl[1][i] = pt_axis[i];
             ncell_dim[i] = 2;
         }
         else {
             ncell_dim[i] = 1;
         }
     }
     for (k = 0; k < ncell_dim[2]; k++) {
         myxl[2] = xxl[k][2] - pt_axis[2];
         myxr[2] = xxr[k][2] - pt_axis[2];
         for (j = 0; j < ncell_dim[1]; j++) {
             myxl[1] = xxl[j][1] - pt_axis[1];
             myxr[1] = xxr[j][1] - pt_axis[1];
             for (i = 0; i < ncell_dim[0]; i++) {
 
                 myxl[0] = xxl[i][0] - pt_axis[0];
                 myxr[0] = xxr[i][0] - pt_axis[0];
 
                 for (idim = 0; idim < dim; idim++) {
                     xc[idim]    = 0.5 * (myxl[idim] + myxr[idim]);
                     myctr[idim] = 0.0;
                 }
                 if (obj == polygon) {  // polygon
                     memcpy(xys_scratch, xys, (size_t)(2*(nnode+1)*sizeof(double)));
                     memcpy(my_x2min, x2min, (size_t)(2 * sizeof(double)));
                     memcpy(my_x2max, x2max, (size_t)(2 * sizeof(double)));
                 }
                 else if (obj == sphere) {  // sphere
                     for (idim = 0; idim < dim; idim++) {
                         myctr[idim] = ctr[idim] - pt_axis[idim];
                     }
                 }
                 if ((xc[0] > 0.0) && (xc[1] > 0.0) && (xc[2] > 0.0)) {
                     a_begin = angle_begin;
                     a_end   = angle_end;
 
                     if (if_helix) {
                         h0_mod     = myh0;
                         dhdphi_mod = dhdphi;
                     }
                 }
                 else if ((xc[0] < 0.0) && (xc[1] > 0.0) && (xc[2] > 0.0)) {
                     tmp     =   myxl[0];
                     myxl[0] = - myxr[0];
                     myxr[0] = - tmp;
 
//                   exchange x and y
 
//                     tmp     = myxl[0];
//                     myxl[0] = myxl[1];
//                     myxl[1] = tmp;
//                     tmp     = myxr[0];
//                     myxr[0] = myxr[1];
//                     myxr[1] = tmp;
 
//                   phi_old = 180 - phi_new
 
                     a_begin = 180.0 - angle_end;
                     a_end   = 180.0 - angle_begin;
 
                     if (if_helix) {
                         h0_mod     = myh0 + 180.0 * dhdphi;
                         dhdphi_mod = - dhdphi;
                     }
                 }
                 else if ((xc[0] > 0.0) && (xc[1] < 0.0) && (xc[2] > 0.0)) {
                     tmp     =   myxl[1];
                     myxl[1] = - myxr[1];
                     myxr[1] = - tmp;
 
//                   exchange x and y
 
//                     tmp     = myxl[0];
//                     myxl[0] = myxl[1];
//                     myxl[1] = tmp;
//                     tmp     = myxr[0];
//                     myxr[0] = myxr[1];
//                     myxr[1] = tmp;
 
//                   phi_old = 360 - phi_old
 
                     a_begin = 360.0 - angle_end;
                     a_end   = 360.0 - angle_begin;
 
                     if (if_helix) {
                         h0_mod     =    myh0 + 360.0 * dhdphi;
                         dhdphi_mod =  - dhdphi;
                     }
                 }
                 else if ((xc[0] < 0.0) && (xc[1] < 0.0) && (xc[2] > 0.0)) {
                     tmp     =   myxl[0];
                     myxl[0] = - myxr[0];
                     myxr[0] = - tmp;
                     tmp     =   myxl[1];
                     myxl[1] = - myxr[1];
                     myxr[1] = - tmp;
 
//                   phi_old  = phi_new - 180
 
                     a_begin = angle_begin + 180.0;
                     a_end   = angle_end   + 180.0;
 
                     if (if_helix) {
                         h0_mod     = myh0 - 180.0 * dhdphi;
                         dhdphi_mod = dhdphi;
                     }
                     if (a_end >= 360.0) {
                         a_begin = a_begin - 360.0;
                         a_end   = a_end   - 360.0;
 
                         h0_mod = h0_mod + 360.0 * dhdphi_mod;
                     }
                 }
                 else if ((xc[0] > 0.0) && (xc[1] > 0.0) && (xc[2] < 0.0)) {
                     tmp     =   myxl[2];
                     myxl[2] = - myxr[2];
                     myxr[2] = - tmp;
 
//                   exchange x and y
 
//                     tmp     = myxl[0];
//                     myxl[0] = myxl[1];
//                     myxl[1] = tmp;
//                     tmp     = myxr[0];
//                     myxr[0] = myxr[1];
//                     myxr[1] = tmp;
 
//                   phi_new = phi_old
 
                     a_begin = angle_begin;
                     a_end   = angle_end;
 
                     if (if_helix) {
                         h0_mod     = - myh0;
                         dhdphi_mod = - dhdphi;
                     }
                 }
                 else if ((xc[0] < 0.0) && (xc[1] > 0.0) && (xc[2] < 0.0)) {
                     tmp     =   myxl[0];
                     myxl[0] = - myxr[0];
                     myxr[0] = - tmp;
                     tmp     =   myxl[2];
                     myxl[2] = - myxr[2];
                     myxr[2] = - tmp;
 
                     a_begin = 180.0 - angle_end;
                     a_end   = 180.0 - angle_begin;
 
                     if (if_helix) {
                         h0_mod     = -(myh0 + dhdphi * 180.0);
                         dhdphi_mod = dhdphi;
                     }
                 }
                 else if ((xc[0] > 0.0) && (xc[1] < 0.0) && (xc[2] < 0.0)) {
                     tmp     =   myxl[1];
                     myxl[1] = - myxr[1];
                     myxr[1] = - tmp;
                     tmp     =   myxl[2];
                     myxl[2] = - myxr[2];
                     myxr[2] = - tmp;
 
                     a_begin = 360.0 - angle_end;
                     a_end   = 360.0 - angle_begin;
 
                     if (if_helix) {  // z_new = -z_old, phinew = 360 - phi_old
                         h0_mod     = -(myh0 + 360.0 * dhdphi);
                         dhdphi_mod =  dhdphi;
                     }
                 }
                 else if ((xc[0] < 0.0) && (xc[1] < 0.0) && (xc[2] < 0.0)) {
 
                     tmp     =   myxl[2];
                     myxl[2] = - myxr[2];
                     myxr[2] = - tmp;
 
                     tmp     =   myxl[0];
                     myxl[0] = - myxr[0];
                     myxr[0] = - tmp;
                     tmp     =   myxl[1];
                     myxl[1] = - myxr[1];
                     myxr[1] = - tmp;
 
//                   exchange x and y
 
//                     tmp     = myxl[0];
//                     myxl[0] = myxl[1];
//                     myxl[1] = tmp;
//                     tmp     = myxr[0];
//                     myxr[0] = myxr[1];
//                     myxr[1] = tmp;
 
                     a_begin = angle_begin + 180.0;
                     a_end   = angle_end   + 180.0;
 
                     if (if_helix) {
                         h0_mod     = -myh0 + 180.0 * dhdphi;
                         dhdphi_mod = -dhdphi;
                     }
                     if (a_end >= 360.0) {
                         a_begin = a_begin - 360.0;
                         a_end   = a_end   - 360.0;
 
                         h0_mod = h0_mod + 360.0 * dhdphi_mod;
                     }
                 }
//               if (a_end - a_begin + small >= 360.0) {
//                   a_begin = 0.0;
//                   a_end   = 360.0;
//               }
                 if (xc[2] < 0.0) {
                     if (obj == polygon) {    // for polygon
                         for (n = 0; n <= nnode; n++) {
                             xys_scratch[n+n+1] = - xys[n+n+1];

                         }
                         my_x2min[1] = - x2max[1];
                         my_x2max[1] = - x2min[1];
                     }
                     else if (obj == sphere) {   // for sphere
                         myctr[1] = - myctr[1];
                     }
                 }
                 c8[0][0] = myxl[0];
                 c8[0][1] = myxl[1];
                 c8[0][2] = myxl[2];
 
                 c8[1][0] = myxr[0];
                 c8[1][1] = myxl[1];
                 c8[1][2] = myxl[2];
 
                 c8[2][0] = myxr[0];
                 c8[2][1] = myxr[1];
                 c8[2][2] = myxl[2];
 
                 c8[3][0] = myxl[0];
                 c8[3][1] = myxr[1];
                 c8[3][2] = myxl[2];
 
                 c8[4][0] = myxl[0];
                 c8[4][1] = myxl[1];
                 c8[4][2] = myxr[2];
 
                 c8[5][0] = myxr[0];
                 c8[5][1] = myxl[1];
                 c8[5][2] = myxr[2];
 
                 c8[6][0] = myxr[0];
                 c8[6][1] = myxr[1];
                 c8[6][2] = myxr[2];
 
                 c8[7][0] = myxl[0];
                 c8[7][1] = myxr[1];
                 c8[7][2] = myxr[2];
 
 
                 nnode8 = 8;
                 nface8 = 6;
                 for (i = 0; i < nface8; i++) {
                     nnode_for_face[i] = 4;
                 }
                 nodelist_for_face[0] =  0;
                 nodelist_for_face[1] =  3;
                 nodelist_for_face[2] =  7;
                 nodelist_for_face[3] =  4;
 
                 nodelist_for_face[4] =  1;
                 nodelist_for_face[5] =  5;
                 nodelist_for_face[6] =  6;
                 nodelist_for_face[7] =  2;
 
                 nodelist_for_face[8]  =  0;
                 nodelist_for_face[9]  =  4;
                 nodelist_for_face[10] =  5;
                 nodelist_for_face[11] =  1;
 
                 nodelist_for_face[12] =  2;
                 nodelist_for_face[13] =  6;
                 nodelist_for_face[14] =  7;
                 nodelist_for_face[15] =  3;
 
                 nodelist_for_face[16] =  0;
                 nodelist_for_face[17] =  1;
                 nodelist_for_face[18] =  2;
                 nodelist_for_face[19] =  3;
 
                 nodelist_for_face[20] =  4;
                 nodelist_for_face[21] =  7;
                 nodelist_for_face[22] =  6;
                 nodelist_for_face[23] =  5;
 
                 nedge8 = 12;
                 nodelist_for_edge[0]  = 0;
                 nodelist_for_edge[1]  = 1;
                 nodelist_for_edge[2]  = 1;
                 nodelist_for_edge[3]  = 2;
                 nodelist_for_edge[4]  = 2;
                 nodelist_for_edge[5]  = 3;
                 nodelist_for_edge[6]  = 3;
                 nodelist_for_edge[7]  = 0;
 
                 nodelist_for_edge[8]  = 4;
                 nodelist_for_edge[9]  = 5;
                 nodelist_for_edge[10] = 5;
                 nodelist_for_edge[11] = 6;
                 nodelist_for_edge[12] = 6;
                 nodelist_for_edge[13] = 7;
                 nodelist_for_edge[14] = 7;
                 nodelist_for_edge[15] = 4;
 
                 nodelist_for_edge[16] = 0;
                 nodelist_for_edge[17] = 4;
                 nodelist_for_edge[18] = 1;
                 nodelist_for_edge[19] = 5;
                 nodelist_for_edge[20] = 2;
                 nodelist_for_edge[21] = 6;
                 nodelist_for_edge[22] = 3;
                 nodelist_for_edge[23] = 7;
 
                 edgelist_for_face[0]  = 3;
                 edgelist_for_face[1]  = 11;
                 edgelist_for_face[2]  = 7;
                 edgelist_for_face[3]  = 8;
 
                 edgelist_for_face[4]  = 1;
                 edgelist_for_face[5]  = 9;
                 edgelist_for_face[6]  = 5;
                 edgelist_for_face[7]  = 10;
 
                 edgelist_for_face[8]  = 0;
                 edgelist_for_face[9]  = 8;
                 edgelist_for_face[10] = 4;
                 edgelist_for_face[11] = 9;
 
                 edgelist_for_face[12] = 2;
                 edgelist_for_face[13] = 10;
                 edgelist_for_face[14] = 6;
                 edgelist_for_face[15] = 11;
 
                 edgelist_for_face[16] = 0;
                 edgelist_for_face[17] = 1;
                 edgelist_for_face[18] = 2;
                 edgelist_for_face[19] = 3;
 
                 edgelist_for_face[20] = 4;
                 edgelist_for_face[21] = 7;
                 edgelist_for_face[22] = 6;
                 edgelist_for_face[23] = 5;
 
 
                 if (a_begin < 0.0)   a_begin = 0.0;
                 if (a_end   >= 90.0) a_end   = 360.0;
 
                 dvol = 0.0;
                 grevolve_obj_p(obj,
                                ifinquiry, mixed, a_begin, a_end,
                                myctr, radius, if_helix, h0_mod, dhdphi_mod,
                                nnode, xys_scratch, my_x2min, my_x2max,
                                c8[0], nnode8, nedge8, nface8,
                                nnode_for_face, nodelist_for_face,
                                nodelist_for_edge, edgelist_for_face,
                                dx[0], vol0, &dvol);
             }
             *vol += dvol;
         }
     }
     return;
 }
 
void grevolve_obj_p(geom_obj_type  obj,
                    int ifinquiry, int *mixed,
                    double angle_begin, double angle_end,
                    double *ctr, double radius,
                    int if_helix, double h0, double dhdphi,
                    int nnode, double *coords2d, double *x2min, double *x2max,
                    double *coords, int nnode8, int nedge8, int nface8,
                    int *nnode_for_face, int *nodelist_for_face,
                    int *nodelist_for_edge, int *edgelist_for_face,
                    double scale, double vcell, double *vol)
{
//   The 8 vertices, coords, are in the assumed order.
 
     const int mxnum = 32;
     int err, whole_ring, anyinside, alloutside, allinside;
     int szdim2, isinside, intersected;
     int nnode_int, nnode_intface, v, i, k, k1, n, n2;
     int inside_v[mxnum], vofphi[mxnum];
 
     double factor, tmp, dz, rxy, r2, sinp, cosp, phi, xplane, area, ctroid[2];
     double cr2[mxnum][2], cp_rotated[mxnum][3], coords2d_intface[2*mxnum];
     double phi8[mxnum], phi_for_sec[mxnum], dc[2], dx2[2], *c, *ct;
     double cmin2[2], cmax2[2], my_x2min[2], my_x2max[2];
     double cosa, pt[3];
     double radius2, phimin, phimax, rmax, rmin, zmax, zmin, hmin, hmax;
     double myxyint[256], *xyint, *xyint_allocated, *my_coords2d, *pcoords2d;
     point_position pos, pos_old;
     double normal[3], point[3];
 
     assert(nnode8 < mxnum);
 
     xyint_allocated = NULL;
     szdim2 = 2 * sizeof(double);
 
     rmax = 0.0;
     rmin = 1.0e+30;
     zmin = 1.0e+30;
     zmax = - zmin;
 
     radius2 = radius * radius;
 
     if ((angle_begin <= 0.0) && (angle_end >= 360.0)) {
         whole_ring = 1;
     }
     else {
         whole_ring = 0;
     }
     factor = 180.0 / PI;
 
     for (v = 0; v < nnode8; v++) {
         c     = coords + (3 * v);
         ct    = cr2[v];
         rxy   = sqrt(c[0] * c[0] + c[1] * c[1]);
         ct[0] = rxy;
         ct[1] = c[2];
     }
     for (v = 0; v < nnode8; v++) {
         phi8[v] = -360.0;
         c    = coords + (3 * v);
         ct   = cr2[v];
         rxy  = ct[0];
         if (rxy != 0.0) {
             sinp = fabs(c[1]) / rxy;
             phi  = asin(sinp);
             phi *= factor;
             if ((c[0] < 0.0) && (c[1] >= 0.0)) {
                 phi = (180.0 - phi);
             }
             else if ((c[0] < 0.0) && (c[1] < 0.0)) {
                 phi += 180.0;
             }
             else if ((c[0] >= 0.0) && (c[1] < 0.0)) {
                 phi = 360.0 - phi;
             }
             phi8[v] = phi;
         }
         else {
             phi8[v] = 0.0;
         }
     }
//   choose between 0 and 360
 
     for (v = 0; v < nnode8; v++) {
         if (phi8[v] >= 0.0) {
             phimin = phi8[v];
             break;
         }
     }
     if (phimin > 90.0) {
         for (v = 0; v < nnode8; v++) {
             if (phi8[v] == 0.0) {
                 phi8[v] = 360.0;
             }
         }
     }
//   sort phi
 
     for (i = 0; i < nnode8; i++) {
         vofphi[i] = i;
     }
     for (v = 0; v < nnode8; v++) {
         phimin = phi8[v];
         for (i = v+1; i < nnode8; i++) {
             if (phimin  > phi8[i]) {
                 tmp     = phi8[i];
                 phi8[i] = phimin;
                 phi8[v] = tmp;
                 phimin  = tmp;
                 k       = vofphi[i];
                 vofphi[i] = vofphi[v];
                 vofphi[v] = k;
             }
         }
     }
     phimin = phi8[0];
     phimax = phi8[nnode8-1];
 
     if ((phimin >= angle_end) || phimax <= angle_begin) {
         *vol = 0,0;
         *mixed = 0;
         return;
     }
     else {
         for (i = 0; i < nnode8; i++) {
             v    = vofphi[i];
             ct   = cr2[v];
             if (ct[0] == 0.0) {
                 inside_v[v] = 1;
                 continue;
             }
             if ((phi8[i] < angle_begin) || (phi8[i] > angle_end)) {
                 inside_v[v] = 0;
             }
             else if ((angle_begin < phi8[i]) && (phi8[i] < angle_end)) {
                 inside_v[v] = 1;
             }
             else if ((phi8[i] == angle_begin) || (phi8[i] == angle_end)) {
                 inside_v[v] = -1;
             }
         }
     }
     my_coords2d = NULL; 
     if (obj == polygon)  {
         
         alloutside = 0;
         if (!if_helix) { 
             for (k = 0; k < 2; k++) {
                 cmin2[k] = cr2[0][k];
                 cmax2[k] = cr2[0][k];
             }
             for (v = 1; v < nnode8; v++) {
                 ct = cr2[v];
                 for (k = 0; k < 2; k++) {
                     if (cmin2[k] > ct[k]) {
                         cmin2[k] = ct[k];
                     }
                     else if (cmax2[k] < ct[k]) {
                         cmax2[k] = ct[k];
                     }
                 }
             }
             for (i = 0; i < 2; i++) {
                 if ((cmin2[i] >= x2max[i]) || (cmax2[i] <= x2min[i])) {
                     alloutside = 1;
                     break;
                 }
             }
             if (alloutside) {
                 *vol = 0.0;
                 *mixed = 0;
                 if (my_coords2d) free(my_coords2d);
                 return;
             }
             for (i = 0; i < nnode8; i++) {
                 v  = vofphi[i];
                 if (inside_v[v] == 0) continue;
    
                 c = cr2[v];
                 pos = point_in_polygon(coords2d, nnode, x2min, x2max, c, &k);
//               pos_old = point_in_polygon_old(coords2d, nnode, c, &k);
//               assert(pos == pos_old);
     
                 if (pos == inside) {
//                   keep original value, 1 or -1
                 }
                 else if (pos == outside) {
                     inside_v[v] = 0;
                 }
                 else if (pos == touching) {
                     inside_v[v] = -1;
                 }
             }  
         } 
         else {   // for helix  
             my_coords2d = (double *) malloc(nnode *szdim2);
             memcpy(my_coords2d, coords2d, (size_t)(nnode * szdim2));

             my_x2min[0] = x2min[0];
             my_x2max[0] = x2max[0]; 

             alloutside = 1;
             for (i = 0; i < nnode8; i++) {
                 v  = vofphi[i];
                 phi = phi8[i];
                 ct  = cr2[v]; 
                 
                 dz  = h0 + dhdphi * phi; 
                 my_x2min[1] = x2min[1] + dz;
                 my_x2max[1] = x2max[1] + dz; 

                 if ((   x2min[0] < ct[0]) && (ct[0] <    x2max[0]) && 
                     (my_x2min[1] < ct[1]) && (ct[1] < my_x2max[1])) { 
                     alloutside = 0;
                     break;
                 }
             }
             if (alloutside) {
                 *vol = 0.0;
                 *mixed = 0;
                 if (my_coords2d) free(my_coords2d);
                 return;
             }
             for (i = 0; i < nnode8; i++) {
                 v  = vofphi[i];
                 phi = phi8[i];

                 dz  = h0 + dhdphi * phi;
                 for (n = 0; n < nnode; n++) {
                     n2 = n + n;
                     my_coords2d[n2+1] = coords2d[n2+1] + dz;
                 }
                 c = cr2[v];
                 pos = point_in_polygon(my_coords2d, nnode, my_x2min, my_x2max, c, &k);
//               pos_old = point_in_polygon_old(coords2d, nnode, c, &k); 
//               assert(pos == pos_old);  

                 if (pos == inside) {
//                   keep original value, 1 or -1
                 }
                 else if (pos == outside) {
                     inside_v[v] = 0;
                 }
                 else if (pos == touching) {
                     inside_v[v] = -1;
                 }
             } 
         }         
     }
     else if (obj == sphere) {   // sphere
 
         for (i = 0; i < nnode8; i++) {
             v  = vofphi[i];
             phi = phi8[i];
 
             ct  = cr2[v];
 
    //       on the plane of rotated x and z
 
             dc[0] = ct[0] - ctr[0];  // rotated x
             if (if_helix) {
                 dc[1] = ct[1] - (h0 + dhdphi * phi + ctr[1]);
             }
             else {
                 dc[1] = ct[1] - ctr[1];    // z
             }
             r2 = dc[0] * dc[0] + dc[1] * dc[1];
 
             inside_v[v]  = 0;
 
             if (rmax < ct[0]) rmax = ct[0];
             if (rmin > ct[0]) rmin = ct[0];
             if (zmax < ct[1]) zmax = ct[1];
             if (zmin > ct[1]) zmin = ct[1];
 
             if (r2 < radius2) {
//               keep the original value, 1 or -1
//               inside_v[v] = 1;  // inside
             }
             else if (r2 > radius2) {
                 inside_v[v] = 0; // outside
             }
             else {  // r2 == radius2)
                 inside_v[v] = -1;  // touch
             }
         }
     }
     anyinside  = 0;
     alloutside = 1;
     allinside  = 1;
 
     for (i = 0; i < nnode8; i++) {
         if (inside_v[i] == 1) {
             anyinside = 1;
             break;
         }
     }
     for (i = 0; i < nnode8; i++) {
         if (inside_v[i]) {
             alloutside = 0;
             break;
         }
     }
     for (i = 0; i < nnode8; i++) {
         if (inside_v[i] == 0) {
             allinside = 0;
             break;
         }
     }
     if (allinside) {
//       Since the cricle or polygon is comvex, this is always r=true.
 
         *vol   = vcell;
         *mixed = 0;
         if (my_coords2d) free(my_coords2d); 
         return;
     }
     else if (anyinside) {
         if (ifinquiry) {
             *mixed = 1;
             *vol   = small;
             if (my_coords2d) free(my_coords2d); 
             return;
         }
     }
     if (alloutside && (obj == sphere)) {
 
         isinside = 0;
         tmp = 1.0/factor;
         for (i = 0; i < nnode8; i++) {
             v  = vofphi[i];
             phi = phi8[i];
 
//           check whether the center of the circle at phi is inside of the cube
 
             cosa  = cos(tmp * phi);
             pt[0] = ctr[0] * cosa;
             pt[1] = ctr[0] * sqrt(MAX(1.0 - cosa * cosa, 0.0));
             if (if_helix) {
                 pt[2] = h0 + dhdphi * phi + ctr[1];
             }
             else {
                 pt[2] = ctr[1];
             }
             isinside = is_node_within_3d_ucell(0, pt, nnode_for_face,
                               nface8, nodelist_for_face, coords, nnode8);
             if (isinside) {
                 *mixed = 1;
                 *vol   = small;
                 if (ifinquiry) {
                     if (my_coords2d) free(my_coords2d); 
                     return;
                 }
                 else {
                     break;
                 }
             }
         }
     }
     else if (alloutside && (obj == polygon)) {
 
         intersected = 0;
 
         for (i = 0; i < nnode8; i++) {
             phi  = phi8[i];
             sinp = sin(phi);
             cosp = cos(phi);
 
//           This is the new approach, and it is supposed to be more robust, accurate, and faster.
 
             normal[0] = sinp;
             normal[1] = cosp;
             normal[2] = 0.0;
             point[0]  = 0.0;
             point[1]  = 0.0;
             point[2]  = 0.0;
             hex_plane(point, normal, coords, coords2d_intface, &nnode_intface);
 
/************************** This is the old approach *******************************
//           rotate -(pi/2 - phi) along z so that phi-plane is the y-plane
 
             for (v = 0; v < nnode8; v++) {
                 c  = coords + (3 * v);
                 ct = cp_rotated[v];
                 ct[0] = sinp * c[0] - cosp * c[1];
                 ct[1] = cosp * c[0] + sinp * c[1];
                 ct[2] = c[2];
             }
             for (v = 0; v < nnode8; v++) {
                 ct = cp_rotated[v];
                 ct[1] -= ctr[0];
                 ct[2] -= ctr[1];
             }
             xplane = 0.0;
 
             polygon_fr_poly3d_cut_by_xplane(xplane, cp_rotated[0],
                                   nnode8, nedge8, nodelist_for_edge,
                                   &nnode_intface, coords2d_intface);
 
*************************************************************************************/
             if (nnode_intface > 2) {
 
                 if (if_helix) {
                     pcoords2d = my_coords2d;
                     for (n = 0; n < nnode; n++) {
                         n2 = n + n;
                         my_coords2d[n2+1] = coords2d[n2+1] + tmp;
                     }
                 }
                 else { 
                     pcoords2d = coords2d;
                 }
//               tmp1d is never used inside ihull2
 
                 memcpy(coords2d_intface + (nnode_intface+nnode_intface),
                        coords2d_intface, (size_t)szdim2);
 
                
                 xyint = NULL;
                 err = hull_test(pcoords2d, nnode, coords2d_intface, nnode_intface,
                                 &xyint, &nnode_int);
 
//                     err = myihull(coords2d, nnode, x2min, x2max, coords2d_intface, nnode_intface,
//                            scale, xyint, &nnode_int);
//                     if (err == -1) {
//                         ihull2(coords2d, nnode, coords2d_intface, nnode_intface,
//                                scale, xyint, &nnode_int, NULL);
//                     }
                 if (xyint) free(xyint);

                 if (nnode_int > 2) {
                     intersected = 1;
                     break;
                 }
             }
         }
         if (!intersected) {
             *vol = 0.0;
             *mixed = 0;
//           if (xyint_allocated) free(xyint_allocated);

             if (my_coords2d) free(my_coords2d); 
             return;
         }
         else if (ifinquiry) {
             *vol = small;
             *mixed = 1;
//           if (xyint_allocated) free(xyint_allocated);

             if (my_coords2d) free(my_coords2d); 
             return;
         }
     }
 
     revolve_obj_p(obj, angle_begin, angle_end,
                   ctr,  radius, if_helix, h0, dhdphi,
                   nnode, coords2d, x2min, x2max,
                   coords, phi8, inside_v, vofphi,
                   nnode8, nedge8, nface8, nodelist_for_edge, edgelist_for_face,
                   nodelist_for_face, nnode_for_face,
                   scale, vcell, vol);
 
     assert(*vol <= vcell + 5.0e-10);
     *vol = MIN(*vol, vcell);
     if ((*vol == 0.0) || (*vol == vcell)) {
         *mixed = 0;
     }
     else {
         *mixed = 1;
     }
//     if (xyint_allocated) free(xyint_allocated);
 
     if (my_coords2d) free(my_coords2d);  

     return;
 }
 
void revolve_obj_p(geom_obj_type  obj, double angle_begin, double angle_end,
                   double *ctr, double radius,
                   int if_helix, double h0, double dhdphi,
                   int nn, double *coords2dp, double *x2min, double *x2max,
                   double *coords8, double *myphis,
                   int *inside_v, int *vofphi,
                   int nnode8, int nedge, int nface,
                   int *nodelist_for_edge, int *edgelist_for_face,
                   int *nodelist_for_face, int *nnode_for_face,
                   double scale,  double vcell, double *vol)
{
     int    allinside, previous, ifstart, ifend;
     int    is_new_nnodes_inside;
     int    v, i, k, s, nsec, offset, nvert, nv, npart;
     int    nvert_sec[8], listv[32];
     int    nvert_merge[8];
     int    *list;
     double phimin, tmp, factor, phi0, phi1, dphi, dvol, ddvol;
     double phis[8], phi_for_sec[8], phi_for_sec_merge[8];
 
     double dvol_larger, dhdphi_rad;
 
     factor = PI / 180.0;
 
     if (if_helix) {
         dhdphi_rad = dhdphi / factor;
     }
     else {
         dhdphi_rad = 0.0;
     }
     memcpy(phis, myphis, (size_t)(nnode8 * sizeof(double)));
 
     for (i = 0; i < nnode8; i++) {
         phis[i] = MAX(phis[i], angle_begin);
         phis[i] = MIN(phis[i], angle_end);
     }
//   phis[0]    += small;
//   phis[nsec] -= small;
 
//   myphis is assumed to be sorted
 
     nsec = 0;
     phi0 = phis[0];
 
     phi_for_sec[0] = phis[0];
     nvert_sec[0]   = 1;
     listv[0]       = vofphi[0];
     offset         = 1;
 
     for (v = 1; v < nnode8; v++) {
         if (phis[v] > phi0) {
             nvert_sec[nsec]++;
             listv[offset] = vofphi[v];
             offset++;
 
             for (i = v+1; i < nnode8; i++) {
                 if (phis[i] == phis[v]) {
                     nvert_sec[nsec]++;
                     listv[offset] = vofphi[i];
                     offset++;
                 }
                 else {
                     break;
                 }
             }
             nsec++;
             phi_for_sec[nsec] = phis[v];
             phi0 = phis[v];
             nvert_sec[nsec] = 1;
             listv[offset] = vofphi[v];
             offset++;
         }
         else {
             nvert_sec[nsec]++;
             listv[offset] = vofphi[v];
             offset++;
         }
     }
//   merge sections
 
     for (i = 0; i < nsec; i++) {
         nvert_merge[i] = 0;
     }
     list = listv;
     s    = 0;
     nv = 0;
     phi_for_sec_merge[0] = phi_for_sec[0];
 
     for (i = 0; i < nsec; i++) {
         phi0 = phi_for_sec[i];
         phi1 = phi_for_sec[i+1];
 
         nvert = nvert_sec[i];
 
         allinside = 1;
         for (k = 0; k < nvert; k++) {
             v  = list[k];
             if (inside_v[v] == 0) {  // outside
                 allinside = 0;
                 break;
             }
         }
         if (!i) {
             nv += nvert;
         }
         else if (previous * allinside) {   // both sections are inside
             nv += nvert;
         }
         else {
             nvert_merge[s] = nv;
             s++;
             phi_for_sec_merge[s] = phi0;
             nv = nvert;
         }
         previous = allinside;
         list += nvert;
     }
     nvert_merge[s] = nv;
     s++;
     phi_for_sec_merge[s] = phi1;
     nsec = s;
     memcpy(nvert_sec, nvert_merge, (size_t)(nsec * sizeof(int)));
     memcpy(phi_for_sec, phi_for_sec_merge, (size_t)((nsec+1) * sizeof(double)));
 
     npart  = 1;
     *vol   = 0.0;
     list   = listv;
 
     for (s = 0; s < nsec; s++) {
 
         phi0  = factor * phi_for_sec[s];
         phi1  = factor * phi_for_sec[s+1];
         nvert = nvert_sec[s];
 
         allinside = 1;
         for (k = 0; k < nvert; k++) {
             v = list[k];
             if (inside_v[v] == 0) {
                 allinside = 0;
                 break;
             }
         }
         if (allinside) {
             ifstart = 0;
             if (s == 0) {
                 ifstart = 1;
             }
             ifend = 0;
             if (s == nsec - 1) {
                 ifend = 1;
             }
             revolve_poly_poly(obj,
                               phi0, phi1, ifstart, ifend,
                               ctr, radius,
                               nn, coords2dp, x2min, x2max,
                               coords8, nnode8, nedge, nface,
                               inside_v,
                               nodelist_for_edge, edgelist_for_face,
                               nodelist_for_face, nnode_for_face,
                               &dvol, &is_new_nnodes_inside);
 
             if (!is_new_nnodes_inside) {
                 revolve_obj_slice(obj,
                                   phi0, phi1, ctr, radius,
                                   if_helix, h0, dhdphi_rad,
                                   nn, coords2dp, x2min, x2max,
                                   coords8, nnode8, nedge, nface,
                                   nodelist_for_edge, edgelist_for_face,
                                   nodelist_for_face, nnode_for_face,
                                   scale, vcell, &dvol);
             }
         }
         else {
             dphi = (phi1 - phi0)/(double)npart;
             phi0 -= dphi;
             dvol = 0.0;
             for (i = 0; i < npart; i++) {
                 phi0 += dphi;
                 phi1  = phi0 + dphi;
                 revolve_obj_slice(obj, phi0, phi1, ctr, radius,
                                    if_helix, h0, dhdphi_rad,
                                    nn, coords2dp, x2min, x2max,
                                    coords8, nnode8, nedge, nface,
                                    nodelist_for_edge, edgelist_for_face,
                                    nodelist_for_face, nnode_for_face,
                                    scale, vcell, &ddvol);
                 dvol += ddvol;
             }
/*****
// DEBUG for sphere
//////////////////
             revolve_poly_poly(obj,
                               phi0, phi1, ifstart, ifend,
                               ctr, radius,
                               nn, coords2dp, x2min, x2max,
                               coords8, nnode8, nedge, nface,
                               inside_v,
                               nodelist_for_edge, edgelist_for_face,
                               nodelist_for_face, nnode_for_face,
                               &dvol_larger, &is_new_nnodes_inside);
 
             assert(dvol_larger + 1.0e-10 > dvol);
//////////////////
****************/
         }
         *vol += dvol;
         list += nvert;
     }
     return;
 }
 
void revolve_poly_poly(geom_obj_type obj,
                      double angle_begin, double angle_end,
                      int ifstart, int ifend,
                      double *ctr, double radius,
                      int nn2d, double *coords2d,
                      double *x2min, double *x2max,
                      double *coords8, int nnode8, int nedge, int nface,
                      int    *inside_v,
                      int *nodelist_for_edge, int *edgelist_for_face,
                      int *nodelist_for_face, int *nnode_for_face,
                      double *vol, int *is_new_nnodes_inside)
{
     int list_default[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
     int checkd, found;
     int i, k, v, f, szdim2, nn, mynn, n, junk;
     int nnode2, nedge2, nface2, nnode3, nedge3, nface3;
     int nnode_out,  nnode_out2, nnode_intface, nnode_intface2;
 
     int nface_for_zone_out[ 2];
     int nface_for_zone_out2[2];
 
     int *nnode_for_face_out[2], *nnode_for_face_out2[2];
     int my_nnode_for_face_out[ 2*16];
     int my_nnode_for_face_out2[2*16];
 
     int *nodelist_for_face_out[2], *nodelist_for_face_out2[2];
     int my_nodelist_for_face_out[ 2 * 64];
     int my_nodelist_for_face_out2[2 * 64];
 
     int nodelist_intface[ 16];
     int nodelist_intface2[16];
     int nodelist_for_zone_out[64], *nodelist_s, *nodelist_t;
 
     int *nnode_for_face2,    *nodelist_for_face2;
     int *edgelist_for_face2, *nodelist_for_edge2;
     int *nnode_for_face3,    *nodelist_for_face3;
     int *edgelist_for_face3, *nodelist_for_edge3;
 
     int inside_v_added[32];
     double phi0, phi1, sinp, cosp, xplane;
 
     double coords_work[3*8];
     double tmp2[2], ctr3[3], coords2d_intface[2*32];
     double coords1d_out[ 3 * 16];
     double coords1d_out2[3 * 16];
 
     double *p0, *p1, *coords2, *coords3, c2[2];
     point_position pos;
 
     szdim2  = sizeof(double);
     szdim2 += szdim2;
 
     for (i = 0; i < 2; i++) {
         k  = i*16;
         nnode_for_face_out[ i]   = my_nnode_for_face_out + k;
         nnode_for_face_out2[i]   = my_nnode_for_face_out2 + k;
         k = i * 64;
         nodelist_for_face_out[ i] = my_nodelist_for_face_out + k;
         nodelist_for_face_out2[i] = my_nodelist_for_face_out2 + k;
     }
     memcpy(inside_v_added, inside_v, (size_t)(nnode8 * sizeof(int)));
 
     phi1 = angle_end;
     phi0 = angle_begin;
 
     *vol = 0.0;
     if (!ifstart) {
 
         sinp = sin(phi0);
         cosp = cos(phi0);
 
//       rotate -(pi/2 - phi) along z so that phi-plane is the y-plane
 
         for (v = 0; v < nnode8; v++) {
             p0 = coords8  + (3 * v);
             p1 = coords_work + (3 * v);
             p1[0] = sinp * p0[0] - cosp * p0[1];
             p1[1] = cosp * p0[0] + sinp * p0[1];
             p1[2] = p0[2];
         }
         xplane = 0.0;
 
         poly3d_cut_by_xplane(3, xplane, nface, nedge, nnode8,
                              coords_work,
                              nnode_for_face,    edgelist_for_face,
                              nodelist_for_edge, nodelist_for_face,
                              &nnode_out, coords1d_out,
                              nface_for_zone_out, nnode_for_face_out,
                              nodelist_for_face_out,
                              &nnode_intface, nodelist_intface);
 
         nnode2             = nnode_out;
         nedge2             = 0;
         nface2             = nface_for_zone_out[0];
         coords2            = coords1d_out;
         nnode_for_face2    = nnode_for_face_out[0];
         nodelist_for_face2 = nodelist_for_face_out[0];
         edgelist_for_face2 = NULL;
         nodelist_for_edge2 = NULL;
 
//       Check whether the nodes newly generated are inside of the polygon
//       (nn2d, coords2dp)
 
         for (i = nnode8; i < nnode8 + nnode_out; i++) {
             inside_v_added[i] = 0;
         }
         check_inside(obj,
                      nface2, nnode_for_face2, nodelist_for_face2,
                      coords2, nnode2, nnode8, inside_v_added,
                      ctr, radius,
                      coords2d, nn2d, x2min, x2max,
                      is_new_nnodes_inside);
 
         if (!(*is_new_nnodes_inside)) {
             return;
         }
//       move back
 
         for (v = 0; v < nnode_out; v++) {
             p0 = coords1d_out + (3 * v);
             memcpy(tmp2, p0, (size_t)szdim2);
             p0[0] = -sinp * tmp2[0] - cosp * tmp2[1];
             p0[1] =  cosp * tmp2[0] - sinp * tmp2[1];
         }
     }
     else {
         nnode2             = nnode8;
         nedge2             = nedge;
         nface2             = nface;
         coords2            = coords8;
         nnode_for_face2    = nnode_for_face;
         nodelist_for_face2 = nodelist_for_face;
         edgelist_for_face2 = edgelist_for_face;
         nodelist_for_edge2 = nodelist_for_edge;
     }
     if (!ifend) {
         sinp = sin(phi1);
         cosp = cos(phi1);
 
//       rotate -(pi/2 - phi) along z so that phi-plane is the y-plane
 
         for (v = 0; v < nnode2; v++) {
             p0 = coords2 + (3 * v);
             p1 = coords_work  + (3 * v);
             p1[0] = sinp * p0[0] - cosp * p0[1];
             p1[1] = cosp * p0[0] + sinp * p0[1];
             p1[2] = p0[2];
         }
         xplane = 0.0;
 
         poly3d_cut_by_xplane(3, xplane, nface2, nedge2, nnode2,
                              coords_work,
                              nnode_for_face2,    edgelist_for_face2,
                              nodelist_for_edge2, nodelist_for_face2,
                              &nnode_out2, coords1d_out2,
                              nface_for_zone_out2, nnode_for_face_out2,
                              nodelist_for_face_out2,
                              &nnode_intface2, nodelist_intface2);
 
         nnode3             = nnode_out2;
         nedge3             = 0;
         nface3             = nface_for_zone_out2[1];
         coords3            = coords1d_out2;
         nnode_for_face3    = nnode_for_face_out2[1];
         nodelist_for_face3 = nodelist_for_face_out2[1];
         edgelist_for_face3 = NULL;
         nodelist_for_edge3 = NULL;
 
//       Check whether the nodes newly generated are inside of the polygon
//       (nn2d, coords2d)
 
         for (i = nnode2; i < nnode2 + nnode_out2; i++) {
             inside_v_added[i] = 0;
         }
         check_inside(obj,
                      nface3, nnode_for_face3, nodelist_for_face3,
                      coords3, nnode3, nnode2, inside_v_added,
                      ctr, radius,
                      coords2d, nn2d, x2min, x2max,
                      is_new_nnodes_inside);
 
         if (!(*is_new_nnodes_inside)) {
             return;
         }
     }
     else {
         nnode3 = nnode2;
         nedge3 = nedge2;
         nface3 = nface2;
         coords3 = coords2;
         nnode_for_face3    = nnode_for_face2;
         nodelist_for_face3 = nodelist_for_face2;
         edgelist_for_face3 = edgelist_for_face2;
         nodelist_for_edge3 = nodelist_for_edge2;
     }
     if (nface3 > 3) {
         assert(nface3 <= 16);
         cal_ctr0(nface3, nnode3, list_default, nnode_for_face3,
                  nodelist_for_face3, coords3,  ctr3, vol);
     }
 
     return;
 }
 
void check_inside(geom_obj_type obj,
                  int nface, int *nnode_for_face, int *nodelist_for_face,
                  double *coords3d, int nnode_tot, int nnode_checkd, int *inside_v,
                  double *ctr, double radius,
                  double *coords2d, int nn2d,
                  double *x2min, double *x2max,
                  int *is_new_nnodes_inside)
{
     int found, junk;
     int f, i, k, n, nn, mynn;
     int nodelist_for_zone[32], *nodelist_t, *nodelist_s;
     double radius2, r2, c2[2], *p0;
     point_position pos, pos_old;
 
     mynn       = 0;
     nodelist_t = nodelist_for_zone;
     nodelist_s = nodelist_for_face;
 
     for (f = 0; f < nface; f++) {
         nn = nnode_for_face[f];
         for (i = 0; i < nn; i++) {
             n  = nodelist_s[i];
             found  = 0;
             for (k = 0; k < mynn; k++) {
                 if (n == nodelist_for_zone[k]) {
                     found = 1;
                     break;
                 }
             }
             if (!found) {
                 *nodelist_t = n;
                 nodelist_t++;
                 mynn++;
             }
         }
         nodelist_s += nn;
     }
     if (obj == polygon) {
 
//       Check whether the nodes newly generated are inside of the polygon
//       (nn2d, coords2d)
 
         *is_new_nnodes_inside = 1;
 
         for (i = 0; i < mynn; i++) {
             n  = nodelist_for_zone[i];
             if (n < nnode_checkd) {
                 if (inside_v[n] == 0) {
                     *is_new_nnodes_inside = 0;
                     break;
                 }
             }
             else {
                 p0 = coords3d + (3 * n);
//               p0[0] should be zero
 
                 c2[0] = p0[1];
                 c2[1] = p0[2];
                 for (k = 0; k < 2; k++) {
                     if (c2[k] <= x2min[k]) {
                         *is_new_nnodes_inside = 0;
                         break;
                     }
                     else if (c2[k] >= x2max[k]) {
                         *is_new_nnodes_inside = 0;
                         break;
                     }
                 }
                 if (!(*is_new_nnodes_inside)) {
                     break;
                 }
                 pos = point_in_polygon(coords2d, nn2d, x2min, x2max, c2, &junk);
//               pos_old = point_in_polygon_old(coords2d, nn2d, c2, &junk);
//               assert(pos_old ==pos);
 
                 if (pos == outside) {
                     *is_new_nnodes_inside = 0;
                     break;
                 }
                 else {
                     inside_v[n] = 1;
                 }
             }
         }
     }
     else if (obj == sphere) {
         radius2 = radius * radius;
         for (i = 0; i < mynn; i++) {
             n  = nodelist_for_zone[i];
             if (n < nnode_checkd) {
                 if (inside_v[n] == 0) {
                     *is_new_nnodes_inside = 0;
                     break;
                 }
             }
             else {
                 p0 = coords3d + (3 * n);
//               p0[0] should be zero
 
                 c2[0] = p0[1] - ctr[0];
                 c2[1] = p0[2] - ctr[1];
                 r2    = c2[0] * c2[0] + c2[1] * c2[1];
                 if (r2 > radius2) {
                     *is_new_nnodes_inside = 0;
                     break;
                 }
                 else {
                     inside_v[n] = 1;
                 }
             }
         }
     }
     return;
 }
 
void revolve_obj_slice(int obj,
                      double angle_begin, double angle_end,
                      double *ctr0, double radius,
                      int if_helix, double h0, double dhdphi,
                      int nn, double *coords2dp, double *x2min, double *x2max,
                      double *coords0, int nnode, int nedge, int nface,
                      int *nodelist_for_edge, int *edgelist_for_face,
                      int *nodelist_for_face, int *nnode_for_face,
                      double scale, double vcell, double *vol)
{
     int mixed, allocated, szdim2;
     int i, i2, k, k1, v, v2, nit, npart, npart2, done;
     int nnode_intface, nnode_int;
 
     double phi0, phi1, dphi, dphi6, phi, sinp, cosp, xplane, sixth;
     double third, tmp, factor, vol_old, dvol, area, area_p, darea, dz, err;
     double ctr[3], ctroid[2], coords2d_intface[2*32];
     double coords[3 * 32], tmp1d[64];
     double *array, *fend, *fm, *p0, *p1, *my_coords2d, *pcoords2d;
     double myxyint[256], *xyint, *xyint_allocated;
     double ctr_p[2], normal[3], point[3];
 
double mytol;
mytol = 1.0e-07;
 
     xyint_allocated = NULL;
 
     third = 0.33333333333333333333333333333;
     szdim2 = 2 * sizeof(double);
 
     *vol  = 0.0;
     sixth = 1.0/6.0;
 
     phi1  = angle_end;
     phi0  = angle_begin;
     npart = 8;
     fend  = NULL;
 
     my_coords2d = NULL;
     if (obj == polygon) {
         if (if_helix) { 
             my_coords2d = (double *) malloc((nn + nn + 2) * sizeof(double));
             memcpy(my_coords2d, coords2dp, (size_t)((nn + nn) * sizeof(double)));
             pcoords2d = my_coords2d; 
         } 
         else { 
             pcoords2d = coords2dp;
         }
     }
     else if (obj == sphere) {
         memcpy(ctr, ctr0, (size_t)szdim2);
     }
/***
//////////////////////////////////////////////////////////////////////////////
//   Simpson integral
 
     fend = (double *) malloc((npart + npart + 1) * sizeof(double));
     fm   = fend + (npart + 1);
 
     dphi = (phi1 - phi0)/(double)npart;
     phi  = phi0 - dphi;
 
     for (i = 0; i <= npart; i++) {
         phi += dphi;
 
         sinp = sin(phi);
         cosp = cos(phi);
 
//       rotate -(pi/2 - phi) along z so that phi-plane is the y-plane
 
         for (v = 0; v < nnode; v++) {
             p0 = coords0 + (3 * v);
             p1 = coords  + (3 * v);
             p1[0] = sinp * p0[0] - cosp * p0[1];
             p1[1] = cosp * p0[0] + sinp * p0[1];
             p1[2] = p0[2];
         }
         xplane = 0.0;
         polygon_fr_poly3d_cut_by_xplane(xplane, coords,
                           nnode, nedge, nodelist_for_edge,
                           &nnode_intface, coords2d_intface);
         area      = 0.0;
         ctroid[0] = 0.0;
         ctroid[1] = 0.0;
 
//       find the area of intersection between coords_intface and the circle
 
         if (nnode_intface > 2) {
 
             memcpy(coords2d_intface + (nnode_intface + nnode_intface),
                    coords2d_intface, (size_t)szdim2);
 
             xyint = NULL;
             err = hull_test(coords2dp, nn, coords2d_intface, nnode_intface,
                             &xyint, &nnode_int);
 
//             err = myihull(coords2dp, nn, x2min, x2max, coords2d_intface, nnode_intface,
//                    scale, xyint, &nnode_int);
//             if (err == -1) {
//                 ihull2(coords2dp, nn, coords2d_intface, nnode_intface,
//                        scale, xyint, &nnode_int, tmp1d);
//             }
//             assert(nnode_intface < nn + nn + nnode);
 
             if (nnode_int > 2) {
                 for (k = 0; k < nnode_int; k++) {
                     k1 = (k + 1) % nnode_int;
                     p0 = xyint + (k  + k);
                     p1 = xyint + (k1 + k1);
                     darea = p0[0] * p1[1] - p0[1] * p1[0];
                     area += darea;
                     ctroid[0] += (p0[0] + p1[0]) * darea;
                     ctroid[1] += (p0[1] + p1[1]) * darea;
                 }
                 area *= 0.5;
                 if (area != 0.0) {
                     factor = 1.0/(6.0 * area);
                     ctroid[0] *= factor;
                     ctroid[1] *= factor;
                 }
                 if (area < 0.0) area = -area;
             }
             if (xyint) free(xyint);
         }
         fend[i] = ctroid[0] * area;
     }
////////////////////////////////////////////////////////////////////////
***/
     vol_old = 1.0e+30;
     done = 0;
     nit = 0;
 
     while (!done) {
          dphi = (phi1 - phi0)/(double)npart;
          dphi6 = sixth * dphi;
          *vol = 0.0;
          phi  = phi0 - 0.5 * dphi;
 
          for (i = 0; i < npart; i++) {
              phi += dphi;
 
              k1   = i % 3;
              if (nit && (k1 == 1)) continue;
 
              sinp = sin(phi);
              cosp = cos(phi);
 
/**********************************************************************************
//            this is a new approach, which is supposed to be more robust,
//            accurate, and fast, but, this is not for arbitrary polyhedron.
 
              normal[0] = sinp;
              normal[1] = cosp;
              normal[2] = 0.0;
              point[0]  = 0.0;
              point[1]  = 0.0;
              point[2]  = 0.0;
              hex_plane(point, normal, coords0, coords2d_intface, &nnode_intface);
 
***********************************************************************************/
//            rotate -(pi/2 - phi) along z so that phi-plane is the y-plane
 
              for (v = 0; v < nnode; v++) {
                  p0 = coords0 + (3 * v);
                  p1 = coords  + (3 * v);
                  p1[0] = sinp * p0[0] - cosp * p0[1];
                  p1[1] = cosp * p0[0] + sinp * p0[1];
                  p1[2] = p0[2];
              }
              xplane = 0.0;
              polygon_fr_poly3d_cut_by_xplane(xplane, coords,
                                nnode, nedge, nodelist_for_edge,
                                &nnode_intface, coords2d_intface);
 
              ctroid[0] = 0.0;
              ctroid[1] = 0.0;
              area      = 0.0;
              if (nnode_intface > 2) {
 
                  if (obj == polygon) {
//                    ctr2d[0] should be the centroid of the overlap in sph_poly2d.
 
                      memcpy(coords2d_intface + (nnode_intface + nnode_intface),
                             coords2d_intface, (size_t)szdim2);
 
                      if (if_helix) { 
                          dz = h0 + dhdphi * phi; 
                          for (v = 0; v < nn; v++) {
                              v2 = v + v;
                              my_coords2d[v2+1] = coords2dp[v2+1] + dz;
                          } 
                      }
                      xyint = NULL;
                      err = hull_test(pcoords2d, nn, coords2d_intface, nnode_intface,
                                      &xyint, &nnode_int);
 
//                    err = myihull(coords2dp, nn, x2min, x2max, coords2d_intface, nnode_intface,
//                             scale, xyint, &nnode_int);
//                      if (err == -1) {
//                          ihull2(coords2dp, nn, coords2d_intface, nnode_intface,
//                                 scale, xyint, &nnode_int, NULL);
//                      }
 
                      if (nnode_int > 2) {
                          for (k = 0; k < nnode_int; k++) {
                              k1 = (k + 1) % nnode_int;
                              p0 = xyint + (k  + k);
                              p1 = xyint + (k1 + k1);
                              darea = p0[0] * p1[1] - p0[1] * p1[0];
                              area += darea;
                              ctroid[0] += (p0[0] + p1[0]) * darea;
                              ctroid[1] += (p0[1] + p1[1]) * darea;
                          }
                          area *= 0.5;
                          if (area != 0.0) {
                              factor = 1.0/(6.0 * area);
                              ctroid[0] *= factor;
                              ctroid[1] *= factor;
                          }
                          if (area < 0.0) area = -area;
                      }
                      if (xyint) free(xyint);
 
                  }
//                fm[i] = ctroid[0] * area;
//                dvol  = dphi6 *(fend[i] + 4.0 * fm[i] + fend[i+1]);
 
                  else if (obj == sphere) {
 
//                    ctr2d[0] should be the centroid of the overlap in sph_poly2d.
 
                      if (if_helix) {
                          ctr[1] = h0 + dhdphi * phi + ctr0[1];
                      }
                      sph_poly2d(0, 1, 2, ctr, radius, coords2d_intface,
                                 nnode_intface, &mixed, &area, ctroid);
 
 
////////////////////////////////
    //                only for debug
                      area_p = 0.0;
                      ctr_p[0] = 0.0;
                      ctr_p[1] = 0.0;
                      for (k = 0; k < nnode_intface; k++) {
                          k1 = (k + 1) % nnode_intface;
                          p0 = coords2d_intface + (k + k);
                          p1 = coords2d_intface + (k1 + k1);
                          darea = p0[0] * p1[1] - p1[0] * p0[1];
                          ctr_p[0] += ((p0[0] + p1[0]) * darea);
                          ctr_p[1] += ((p0[1] + p1[1]) * darea);
                          area_p += darea;
                      }
                      area_p *= 0.5;
                      if (area_p != 0.0) {
                          factor = 1.0/(6.0 * area_p);
                          ctr_p[0] *= factor;
                          ctr_p[1] *= factor;
                      }
                      if (area_p < 0.0) area_p = - area_p;
                      assert(area < area_p + small);
////////////////////////////////
 
                  }
              }
              dvol  = ctroid[0] * area * dphi;
 
              *vol += dvol;
          }
          if (nit) {
              *vol += third * vol_old;
          }
          err = fabs(*vol - vol_old)/vcell;
          if (err > mytol) {
              npart2 = npart + npart + npart;
/****
//////////////
              npart2 = npart + npart;
              array  = (double *) malloc((npart2 + npart2 + 1) * sizeof(double));
              assert(array);
 
              for (i = 0; i < npart; i++) {
                  i2 = i + i;
                  array[i2]   = fend[i];
                  array[i2+1] = fm[i];
              }
              array[npart2] = fend[npart];
              free(fend);
              fend = array;
              fm   = fend + (npart2 + 1);
//////////////
***/
              npart = npart2;
              nit++;
          }
          else {
              done = 1;
          }
          vol_old = *vol;
     }
//   if (xyint_allocated) free(xyint_allocated);
     if (fend) free(fend);
     if (my_coords2d) free(my_coords2d);

     return;
 }
 
 
 
 
void polygon_fr_poly3d_cut_by_xplane(double xplane, double *coords,
                                     int nnode, int nedge, int *nodelist_for_edge,
                                     int *nnode_intface, double *coords2d_intface)
{
     const int mxnum = 16;
     int   found, szdim2, szdim;
     int   i, k, n0, n, np, np_added, failed;
 
     int   myside, sides[2];
     int   node_included[2*mxnum];
     int   offplane_ea_node[mxnum], side_ea_edge[8 * mxnum];
     int   nn_onplane_ea_edge[8 * mxnum];
     int   intersectd_ea_edge[8 * mxnum];
     int   nodelist_intface[2 * mxnum];
     int   edge2intfnode[8 * mxnum];
     int   nodeid_added2np[mxnum], old2new_edgeid[2*mxnum][2];
 
     int   *nodelist;
 
     double  eta;
     double  *c, *ci;
 
     double  *c2[2];
     double  coords_int2d[mxnum][2], coords_int2d_work[mxnum+mxnum];
 
//   determin whether each node is on the plane
//   offplane_ea_node[i] =  0  :  node i is on the plane
//                       = -1  :  x[i] < xplane
//                       =  1  :  x[i] > xplane
//
 
 
//   determine whether each edge is cross the plane, excluding the touch
 
//   side_ea_edge[i] =  2: both nodes of the edge, x > xplane
//                   = -2: both nodes of the edge, x < xplane
//                   =  1: one node of the edge x > xplane, the other on the plane
//                   = -1: one node of the edge x < xplane, the other on the plane
//                   =  0: both nodes on the plane, or one node on one side of the plane
 
//   nn_onplane_ea_edge[i] =  2: both nodes are on the plane.
//   nn_onplane_ea_edge[i] =  1: one node is on the plane and the other x > xplane
//   nn_onplane_ea_edge[i] = -1: one node is on the plane and the other x < xplane
//   nn_onplane_ea_edge[i] =  0: no node on the plane
 
     szdim  = sizeof(double);
     szdim2 = szdim + szdim;
     szdim  += szdim2;
 
///////////////////////////////////////////////////////////////////////////////////
//   This part is exactly the same as one part in poly3d_cut_by_xplane.
//   The two could be consolidated.
 
 
     myside = 0;
     for (i = 0; i < nnode; i++) {
         c = coords + (3 * i);
         offplane_ea_node[i] = 0;
         if (fabs(c[0] - xplane) < small) {
             offplane_ea_node[i] = 0;
         }
         else if (c[0] < xplane) {
             offplane_ea_node[i] = -1;
             myside = -1;
         }
         else {
             offplane_ea_node[i] = 1;
             myside = 1;
         }
     }
     *nnode_intface = 0;
 
//   check whether all nodes are on the same side of the plane
 
     found = 0;
     for (i = 0; i < nnode; i++) {
         if (offplane_ea_node[i] * myside < 0) {
             found = 1;
             break;
         }
     }
     if (!found) {
         return;
     }
     for (i = 0; i < nedge; i++) {
         side_ea_edge[i] = 0;
 
         nodelist = nodelist_for_edge + (i + i);
         nn_onplane_ea_edge[i] = 0;
 
         n0 = nodelist[0];
         n  = nodelist[1];
 
         if ((offplane_ea_node[n0] < 0) && (offplane_ea_node[n] < 0)) {
             side_ea_edge[i]       = -2;
             nn_onplane_ea_edge[i] = 0;
         }
         else if ((offplane_ea_node[n0] > 0) && (offplane_ea_node[n] > 0)) {
             side_ea_edge[i]       = 2;
             nn_onplane_ea_edge[i] = 0;
         }
         else if (offplane_ea_node[n0] * offplane_ea_node[n] < 0) {
             side_ea_edge[i]       = 0;
             nn_onplane_ea_edge[i] = 0;
         }
         else if ((offplane_ea_node[n0] == 0) && (offplane_ea_node[n] == 0)) {
             side_ea_edge[i]       = 0;
             nn_onplane_ea_edge[i] = 2;
         }
         else if (offplane_ea_node[n0] == 0) {
             if (offplane_ea_node[n] < 0) {
                 side_ea_edge[i] = -1;
                 nn_onplane_ea_edge[i] = -1;
             }
             else if (offplane_ea_node[n] > 0) {
                 side_ea_edge[i] = 1;
                 nn_onplane_ea_edge[i] = 1;
             }
         }
         else if (offplane_ea_node[n] == 0) {
             if (offplane_ea_node[n0] < 0) {
                 side_ea_edge[i] = -1;
                 nn_onplane_ea_edge[i] = -1;
             }
             else if (offplane_ea_node[n0] > 0) {
                 side_ea_edge[i] = 1;
                 nn_onplane_ea_edge[i] = 1;
             }
         }
     }
//   calculate the intersectios between each edge and the plane
 
     for (i = 0; i < nnode; i++) {
         node_included[i] = 0;
     }
     np_added = 0;
     np       = 0;
     for (i = 0; i < nedge; i++) {
         intersectd_ea_edge[i] =  0;
         edge2intfnode[i]      = -1;
 
         nodelist = nodelist_for_edge + (i + i);
         if (nn_onplane_ea_edge[i] == 2) {
             intersectd_ea_edge[i] = 1;
             for (k = 0; k < 2; k++) {
                 n = nodelist[k];
                 if (!node_included[n]) {
                     ci = coords + (3 * n + 1);
                     memcpy(coords_int2d[np], ci, (size_t)szdim2);
                     node_included[n] = 1;
//                   nodelist_intface[np] = n;
                     np++;
                 }
             }
         }
         else if ((nn_onplane_ea_edge[i] ==  1) ||
                  (nn_onplane_ea_edge[i] == -1)) {
             intersectd_ea_edge[i] = 1;
             for (k = 0; k < 2; k++) {
                 n = nodelist[k];
                 if ((offplane_ea_node[n] == 0) && !node_included[n]) {
                     ci = coords + (3 * n + 1);
                     memcpy(coords_int2d[np], ci, (size_t)szdim2);
                     node_included[n] = 1;
//                   nodelist_intface[np] = n;
                     np++;
                 }
             }
         }
         else if ((nn_onplane_ea_edge[i] == 0)&&(side_ea_edge[i] == 0)) {
 
//           calculate the intersection
 
             for (k = 0; k < 2; k++) {
                 n  = nodelist[k];
                 c2[k] = coords + (3 * n);
                 sides[k] = offplane_ea_node[n];
             }
             assert(sides[0] * sides[1] != 0);
             if (sides[0] * sides[1] < 0) {
 
                 eta = (xplane - c2[0][0])/(c2[1][0] - c2[0][0]);
                 eta = MAX(0.0, MIN(1.0, eta));
                 c = coords_int2d[np];
                 c[0] = c2[0][1] + eta *(c2[1][1] - c2[0][1]);
                 c[1] = c2[0][2] + eta *(c2[1][2] - c2[0][2]);
                 intersectd_ea_edge[i] = 1;
 
//               nodelist_intface[np]  = nnode + np_added;
//               edge2intfnode[i]      = nnode + np_added;
//               nodeid_added2np[np_added] = np;
//               np_added++;
 
                 np++;
             }
         }
     }
     assert((np == 0) || (np >= 3));
 
//   put np nodes into the correct order
 
     *nnode_intface = 0;
     if (np > 3) {
         for (i = 0; i < np; i++) {
             nodelist_intface[i] = i;
         }
         memcpy(coords_int2d_work, coords_int2d[0],(size_t)(np * szdim2));
         mychull_sort(coords_int2d_work, nodelist_intface, np, &failed);
         *nnode_intface = np;
         memcpy(coords2d_intface, coords_int2d_work, (size_t)(np * szdim2));
     }
/////////////////////////////////////////////////////////////////////////////
 
     return;
 }
 
void mychull_sort(double *ptr, int *map, int npoint, int *failed)
{
//   The first two vertices are assimes to be adjacent.
 
     int i, i0, i1, i2, i3, k, itmp, szdim, found;
     double  small2;
     double factor, a, b, axb, sina, ab, normz, mynorm, t1, t2, t3;
     double x0[2], dx[2], ctr[2], dxy0[2], dxy1[2], tmp2[2];
     double myangle[256], myproj1[256], myproj0[256], myuvector[2*256];
     double *angle, *angle_allocated, *proj1, *proj0, *uvector;
     double *p0, *p1, *p2, *p3, *c1, *c2, *c3;
 
     angle_allocated = NULL;
 
     small2 = small * small;
     szdim = 2 * sizeof(double);
 
     if (npoint < 256) {
        angle = myangle;
        proj1    = myproj1;
        proj0    = myproj0;
        uvector  = myuvector;
     }
     else {
        k = npoint + 1;
        angle_allocated = (double *) malloc(5 * k * sizeof(double));
        angle = angle_allocated;
        proj0   = angle + k;
        proj1   = proj0 + k;
        uvector = proj1 + k;
     }
     dx[0] = 0.0;
     dx[1] = 0.0;
     x0[0] = ptr[0];
     x0[1] = ptr[1];
     for (i = 0; i < npoint; i++) {
         i1 = (i + 1) % npoint;
         p0 = ptr + (i + i);
         p1 = ptr + (i1 + i1);
         for (k = 0; k < 2; k++) {
             a = fabs(p1[k] - p0[k]);
             if (a > dx[k]) dx[k] = a;
         }
     }
     for (i = 0; i < npoint; i++) {
         p0 = ptr + (i + i);
         for (k = 0; k < 2; k++) {
             p0[k] = (p0[k] - x0[k])/dx[k];
         }
     }
//   center of these points
 
     ctr[0] = 0.0;
     ctr[1] = 0.0;
     for (i = 0; i < npoint; i++) {
         p0 = ptr + (i + i);
         ctr[0] += p0[0];
         ctr[1] += p0[1];
     }
     factor = 1.0/(double) npoint;
     ctr[0] *= factor;
     ctr[1] *= factor;
 
     for (i = 0; i < npoint; i++) {
         p0 = ptr + (i + i);
         p1 = uvector + (i + i);
         p1[0] = p0[0] - ctr[0];
         p1[1] = p0[1] - ctr[1];
     }
//   find the unit vectors
     for (i = 0; i < npoint; i++) {
         p0 = uvector + (i + i);
         factor = 1.0/sqrt(p0[0] * p0[0] + p0[1] * p0[1]);
         p0[0] *= factor;
         p0[1] *= factor;
     }
     for (i = 0; i < npoint; i++) {
         proj1[i] = -2.0;
         proj0[i] = -2.0;
     }
     i0 = 0;
     p0 = uvector + (i0 + i0);
 
//   pick closest node as i1
 
     for (i = 1; i < npoint; i++) {
         p1 = uvector + (i + i);
         proj0[i] = p0[0] * p1[0] + p0[1] * p1[1];
     }
     i1 = 1;
     a = proj0[i1];
     for (i = 2; i < npoint; i++) {
         if (proj0[i] > a) {
             i1 = i;
             a  = proj0[i];
         }
     }
     if (i1 != 1) {
         p0 = ptr + 2;
         p1 = ptr + (i1 + i1);
         memcpy(tmp2, p0,   (size_t)szdim);
         memcpy(p0,   p1,   (size_t)szdim);
         memcpy(p1,   tmp2, (size_t)szdim);
 
         p0 = uvector + 2;
         p1 = uvector + (i1 + i1);
         memcpy(tmp2, p0,   (size_t)szdim);
         memcpy(p0,   p1,   (size_t)szdim);
         memcpy(p1,   tmp2, (size_t)szdim);
 
         itmp    = map[1];
         map[1]  = map[i1];
         map[i1] = itmp;
 
         i1 = 1;
     }
     for (i = 0; i < npoint - 2; i++) {
//       determine the vertex i2
         i1 = i + 1;
         p0 = ptr + (i  + i );
         p1 = ptr + (i1 + i1);
         a = 0.0;
         for (k = 0; k < 2; k++) {
             dxy0[k] = p1[k] - p0[k];
             a += (dxy0[k] * dxy0[k]);
         }
         a = sqrt(a);
         factor = 1.0/a;
         for (k = 0; k < 2; k++) {
             dxy0[k] *= factor;
         }
         for (i2 = 0; i2 <= i1; i2++) {
             angle[i2] = 1.0e+30;;
         }
         for (i2 = i1 + 1; i2 < npoint; i2++) {
             p2 = ptr + (i2 + i2);
             b = 0.0;
             for (k = 0; k < 2; k++) {
                 dxy1[k] = p2[k] - p1[k];
                 b += (dxy1[k] * dxy1[k]);
             }
             b   = sqrt(b);
             factor = 1.0/b;
             for (k = 0; k < 2; k++) {
                 dxy1[k] *= factor;
             }
             axb = dxy0[0] * dxy1[1] - dxy0[1] * dxy1[0];
             sina = fabs(axb);
             angle[i2] = sina; // asin(sina); I need only relative values
             ab = dxy0[0] * dxy1[0] + dxy0[1] * dxy1[1];
             if (ab < 0.0) {
                 angle[i2] = PI - angle[i2];
             }
         }
//       find the smallest angle
         i2 = i1 + 1;
         for (k = i1 + 1; k < npoint; k++) {
             if (angle[k] < angle[i2]) {
                 i2 = k;
             }
         }
//       whether there is another point with the angle
 
         a = angle[i2];
         found = 0;
         for (k = i1 + 1; k < npoint; k++) {
             if (k == i2) continue;
 
             if (a == angle[k]) {
                 found = 1;
                 i3 = k;
                 break;
             }
         }
         if (found) {  // i1, i2, and i3 are in a straight line
             c1 = ptr + (i1 + i1);
             c2 = ptr + (i2 + i2);
             c3 = ptr + (i3 + i3);
 
             if (fabs(c2[0] - c1[0]) > small) {
                 t1 = (c1[0] - c2[0])/(c3[0] - c2[0]);
                 t2 = (c2[0] - c1[0])/(c3[0] - c1[0]);
                 t3 = (c3[0] - c1[0])/(c2[0] - c1[0]);
             }
             else if (fabs(c2[1] - c1[1]) > small) {
                 t1 = (c1[1] - c2[1])/(c3[1] - c2[1]);
                 t2 = (c2[1] - c1[1])/(c3[1] - c1[1]);
                 t3 = (c3[1] - c1[1])/(c2[1] - c1[1]);
             }
             else {
                 printf("ERROR: i1, i2, and i3 are in a straight line\n");
             }
             if (((0.0 <= t3) && (t3 <= 1.0) && (i2 < i3)) ||
                 ((0.0 <= t2) && (t2 <= 1.0) && (i3 < i2)) ) { // swap i2 and i3
 
                 p2 = ptr + (i2 + i2);
                 p3 = ptr + (i3 + i3);
                 memcpy(tmp2, p2,   (size_t)szdim);
                 memcpy(p2,   p3,   (size_t)szdim);
                 memcpy(p3,   tmp2, (size_t)szdim);
 
                 p2 = uvector + (i2 + i2);
                 p3 = uvector + (i3 + i3);
                 memcpy(tmp2, p2,   (size_t)szdim);
                 memcpy(p2,   p3,   (size_t)szdim);
                 memcpy(p3,   tmp2, (size_t)szdim);
 
                 a = angle[i2];
                 angle[i2] = angle[i3];
                 angle[i3] = a;
             }
         }
         i0 = i2;
 
         i2 = i1 + 1;
         if (i0 != i2) {
             p0 = ptr + (i0 + i0);
             p2 = ptr + (i2 + i2);
             memcpy(tmp2, p0,   (size_t)szdim);
             memcpy(p0,   p2,   (size_t)szdim);
             memcpy(p2,   tmp2, (size_t)szdim);
 
             p0 = uvector + (i0 + i0);
             p2 = uvector + (i2 + i2);
             memcpy(tmp2, p0,   (size_t)szdim);
             memcpy(p0,   p2,   (size_t)szdim);
             memcpy(p2,   tmp2, (size_t)szdim);
 
             a = angle[i0];
             angle[i0] = angle[i2];
             angle[i2] = a;
 
             itmp    = map[i0];
             map[i0] = map[i2];
             map[i2] = itmp;
         }
     }
//   check the order of vertices
 
//   find a non-zero cross product
 
     for (i = 0; i < npoint; i++) {
         i1 = (i  + 1) % npoint;
         i2 = (i1 + 1) % npoint;
         p0 = ptr + (i  + i);
         p1 = ptr + (i1 + i1);
         p2 = ptr + (i2 + i2);
         for (k = 0; k < 2; k++) {
             dxy0[k] = p1[k] - p0[k];
             dxy1[k] = p2[k] - p1[k];
         }
         normz = dxy0[0] * dxy1[1] - dxy0[1] * dxy1[0];
         if (fabs(normz) > 1.0e-03) break;
     }
     *failed = 0;
     for (i = 0; i < npoint; i++) {
         i1 = (i  + 1) % npoint;
         i2 = (i1 + 1) % npoint;
         p0 = ptr + (i  + i );
         p1 = ptr + (i1 + i1);
         p2 = ptr + (i2 + i2);
         for (k = 0; k < 2; k++) {
             dxy0[k] = p1[k] - p0[k];
             dxy1[k] = p2[k] - p1[k];
         }
         mynorm = dxy0[0] * dxy1[1] - dxy0[1] * dxy1[0];
         if (fabs(mynorm) < 1.0e-03) mynorm = 0.0;
 
         if (normz * mynorm  < 0.0) {
             *failed = 1;
         }
//       assert(normz * mynorm >= 0.0);
//       due to rounding-off error, it is possible to have additional points
//       generated during the formation the npoint points.
     }
 
//   scale back
 
     for (i = 0; i < npoint; i++) {
         p0 = ptr + (i + i);
         for (k = 0; k < 2; k++) {
             p0[k] = p0[k] * dx[k] + x0[k];
         }
     }
 
     if (angle_allocated) free(angle_allocated);
 
     return;
 }
 
 
void is_point_inside_polygon(double *p, double *coords, int nnode, int *isinside)
{
     int i, i1, nn, nleft, nright;
     double *p0, *p1, *xint, t, x;
 
     xint = (double *) malloc(nnode * sizeof(double));
     nn = 0;
     for (i = 0; i < nnode; i++) {
         i1 = (i + 1) % nnode;
         p0 = coords + (i  + i);
         p1 = coords + (i1 + i1);
         if (p0[0] == p1[0]) {
             if (((p0[1] <= p[1])&&(p[1] <= p1[1])) ||
                 ((p0[1] >= p[1])&&(p[1] >= p1[1]))) {
                 xint[nn] = p0[0];
                 nn++;
             }
         }
         else if (p0[1] == p1[1]) {
         }
         else {
             t = (p[1] - p0[1])/(p1[1] - p0[1]);
             x = p0[0] + t * (p1[0] - p0[0]);
             if ((t >= 0.0) && (t <= 1.0)) {
                 xint[nn] = x;
                 nn++;
             }
         }
     }
     x = p[0];
     nleft  = 0;
     nright = 0;
 
     for (i = 0; i < nn; i++) {
         if (xint[i] < x) {
             nleft++;
         }
         else if (xint[i] > x) {
             nright++;
         }
     }
     if ((nleft % 2) && (nright % 2)) {
         *isinside = 1;
     }
     else {
         *isinside = 0;
     }
     free(xint);
 
     return;
 }
 
/////////////////////////////////////////////////////////////////////////////////
int hull_test(double *myr1, int nn1, double *myr2, int nn2, double **pr, int *n)
{
     int szdim, i, n1, n2, h, myn;
     double *p0, *p1, *p2, *r1, *r2, dp1[2], dp2[2], tmp[2], cross1, cross2;
     poly_t *poly1, *poly2, *polyi;
 
     *n = 0;
     *pr = NULL;
 
     szdim = 2 * sizeof(double);
 
     r1 = (double *) malloc((nn1 + nn2 + 2) * szdim);
     r2 = r1 + (nn1 + nn1 + 2);
 
//   remove redundant nodes
 
     memcpy(r1, myr1, (size_t)szdim);
     n1  = 1;
     p0  = r1;
     myn = nn1;
     for (i = nn1-1; i > 0; i--) {
         p1 = myr1 + (i + i);
         if ((p0[0] != p1[0]) || (p0[1] != p1[1])) {
             break;
         }
         else {
             myn--;
         }
     }
     for (i = 1; i < myn; i++) {
         p0 = r1 + (n1 + n1 - 2);
         p1 = myr1 + (i + i);
         if ((p0[0] != p1[0]) || (p0[1] != p1[1])) {
             memcpy(r1 + (n1 + n1), p1, (size_t)szdim);
             n1++;
         }
     }
     memcpy(r2, myr2, (size_t)szdim);
     n2  = 1;
     p0  = r2;
     myn = nn2;
     for (i = nn2-1; i > 0; i--) {
         p1 = myr2 + (i + i);
         if ((p0[0] != p1[0]) || (p0[1] != p1[1])) {
             break;
         }
         else {
             myn--;
         }
     }
     for (i = 1; i < myn; i++) {
         p0 = r2 + (n2 + n2 - 2);
         p1 = myr2 + (i + i);
         if ((p0[0] != p1[0]) || (p0[1] != p1[1])) {
             memcpy(r2 + (n2 + n2), p1, (size_t)szdim);
             n2++;
         }
     }
     if ((n1 < 3) || (n2 < 3)) {
         free(r1);
         *n = 0;
         return 0;
     }
//   check the order of nodes
     
     p0 = r1;
     p1 = p0 + 2;
     p2 = p1 + 2;
     for (i = 0; i < 2; i++) {
         dp1[i] = p1[i] - p0[i];
         dp2[i] = p2[i] - p0[i];
     }
     cross1 = dp1[0] * dp2[1] - dp2[0] * dp1[1];

     p0 = r2;
     p1 = p0 + 2;
     p2 = p1 + 2;
     for (i = 0; i < 2; i++) {
         dp1[i] = p1[i] - p0[i];
         dp2[i] = p2[i] - p0[i];
     }
     cross2 = dp1[0] * dp2[1] - dp2[0] * dp1[1];

     if (cross1 * cross2 < 0.0) {
         h = n1/2;
         for (i = 0; i < h; i++) {
             memcpy(tmp, r1 + (i+i), (size_t)szdim);
             memcpy(r1 + (i+i), r1 + (n1+n1-2-i-i), (size_t)szdim);
             memcpy(r1 +(n1+n1-2-i-i), tmp, (size_t)szdim);
         }
     }
     poly1 = (poly_t *) calloc((size_t)1, sizeof(poly_t));
     poly2 = (poly_t *) calloc((size_t)1, sizeof(poly_t));
 
     poly1->alloc = n1;
     poly1->v = (vec_t *)r1;
     poly1->len = n1;
     poly2->alloc = n2;
     poly2->v = (vec_t *)r2;
     poly2->len = n2;
 
     polyi = poly_clip(poly1, poly2);
 
     *n = polyi->len;
     *pr = (double *) malloc((*n + *n + 2) * sizeof(double));
     memcpy(*pr, polyi->v, (size_t)((*n) * szdim));

//   check redundant 


     memcpy((*pr) + (*n + *n), *pr, (size_t)szdim);
 
     free(polyi->v);
     free(polyi);
     free(poly1);
     free(poly2);
 
     free(r1);
 
     return 0;
 }
 
 
void hex_plane(double *point, double *normal, double *vertices, double *coord_res, int *nnode)
{
//   coord_res is at most (3*6) long.
//   The array, vertices, is in my normal order.
 
     int sets[3][4]      = {{0,1,5,6},{0,3,2,6},{0,4,7,6}};
     int dash_sets[3][2] = {{1,2},{3,7},{4,5}};
 
     int dim, nset, ne, set, i, e, n0, n1, n, same;
     int *nodelist, *dash;
     double tiny_for_same_coord, dx, dl, dlmax, top, bot, t;
     double norm[3], *p0, *p1, *p;
 
     tiny_for_same_coord = 1.0e-12;
 
     *nnode = 0;
     nset = 3;
     ne   = 3;
     dim  = 3;
 
     p  = vertices + (dim * 7);
     dlmax = 0.0;
     for (i = 0; i < dim; i++) {
         dx = p[i] - vertices[i];
         dlmax += (dx * normal[i]);
     }
     if (dlmax > 0.0) {
         for (i = 0; i < dim; i++) {
             norm[i] = normal[i];
         }
     }
     else {
         for (i = 0; i < dim; i++) {
             norm[i] = -normal[i];
         }
         dlmax = -dlmax;
     }
     dl = 0.0;
     for (i = 0; i < dim; i++) {
         dx = point[i] - vertices[i];
         dl += (dx * norm[i]);
     }
     if ((dl <= 0.0) || (dl >= dlmax)) {
         return;
     }
     for (set = 0; set < nset; set++) {
         nodelist = sets[set];
         dash     = dash_sets[set];
 
         for (e = 0; e < ne; e++) {
             n0 = nodelist[e];
             n1 = nodelist[e+1];
             p0 = vertices + (dim * n0);
             p1 = vertices + (dim * n1);
             top = 0.0;
             bot = 0.0;
             for (i = 0; i < dim; i++) {
                 top += (norm[i] * p0[i]);
                 bot += (norm[i] *(p1[i] - p0[i]));
             }
             top = dl - top;
             if (bot == 0.0) continue;
 
             t = top/bot;
             if ((t < 0.0) || (t > 1.0)) continue;
 
             p = coord_res + ((*nnode) * dim);
             for (i = 0; i < dim; i++) {
                 p[i] = p0[i] + t *(p1[i] - p0[i]);
             }
             same = 0;
             for (n = 0; n < *nnode; n++) {
                 p0 = coord_res + (dim * n);
                 same = 1;
                 for (i = 0; i < dim; i++) {
                     if (fabs(p0[i] - p[i]) > tiny_for_same_coord) {
                         same = 0;
                         break;
                     }
                 }
             }
             if (same) continue;
 
             (*nnode)++;
 
             n0 = dash[0];
             n1 = dash[1];
             p0 = vertices + (dim * n0);
             p1 = vertices + (dim * n1);
 
             top = 0.0;
             bot = 0.0;
             for (i = 0; i < dim; i++) {
                 top += (norm[i] * p0[i]);
                 bot += (norm[i] *(p1[i] - p0[i]));
             }
             top = dl - top;
             if (bot != 0.0) {
                 t = top/bot;
                 if ((t >= 0.0) && (t <= 1.0)) {
                     p = coord_res + ((*nnode) * dim);
                     for (i = 0; i < dim; i++) {
                         p[i] = p0[i] + t *(p1[i] - p0[i]);
                     }
                     same = 0;
                     for (n = 0; n < *nnode; n++) {
                         p0 = coord_res + (dim * n);
                         same = 1;
                         for (i = 0; i < dim; i++) {
                             if (fabs(p0[i] - p[i]) > tiny_for_same_coord) {
                                 same = 0;
                                 break;
                             }
                         }
                     }
                     if (!same) {
                         (*nnode)++;
                     }
                 }
             }
             break;
         }
     }
     return;
  }

point_position point_in_polygon(double *myr, int n, double *myxl, double *myxr,
                                double *p, int *i0_touch)
{ 
     int   szdim, i, i1, myn, nn;
     int   *new2old;
     double t, curl0, curl, x0, y0, x1, y1;
     double *r, *p0, *p1;

     szdim = 2 * sizeof(double);
     new2old = NULL;

     *i0_touch = -1;

//   check p is any vertex

     for (i = 0; i < n; i++) {
         p0 = myr + (i + i);
         if ((p0[0] == p[0]) && (p0[1] == p[1])) { 
             *i0_touch = i; 
             return touching; 
         }
     }

//   remove the redundant point

     r = (double *) malloc((n + 1) * 2 * sizeof(double));
     new2old = (int *) malloc(n * sizeof(int));

     memcpy(r, myr, (size_t)szdim);
     new2old[0] = 0;
     nn = 1;
     
     p0  = r; 
     myn = n;
     for (i = n-1; i > 0; i--) { 
         p1 = myr + (i + i);
         if ((p0[0] != p1[0]) || (p0[1] != p1[1])) { 
             break;
         }
         else { 
             myn--;
         }
     } 
     for (i = 1; i < myn; i++) { 
         p0 = r + (nn + nn - 2);
         p1 = myr + (i + i);
         if ((p0[0] != p1[0]) || (p0[1] != p1[1])) { 
             memcpy(r + (nn + nn), p1, (size_t)szdim);
             new2old[nn] = i;
             nn++;
         }
     }
     if (nn < 3) { 
         free(r);
         if (new2old) free(new2old);
         return outside; 
     } 

//   check whether p is along any edge 

     for (i = 0; i < nn; i++) { 
         i1 = (i + 1) % nn;
         p0 = r + (i  + i );
         p1 = r + (i1 + i1);
         x0 = p[0] - p0[0];
         y0 = p[1] - p0[1];
         x1 = p[0] - p1[0];
         y1 = p[1] - p1[1];
         curl = x0 * y1 - x1 * y0;
         if (fabs(curl) < tiny) {  
             if (p1[0] != p0[0]) { 
                 t = (p[0] - p0[0])/(p1[0] != p0[0]); 
             }
             else if (p1[1] != p0[1]) {
                 t = (p[1] - p0[1])/(p1[1] != p0[1]);
             }
             else { 
                 printf("ERROR: same nodes in point_in_polygon\n");
                 printf("i = %d of nn = %d\n", i, nn); 
                 printf("p0 = (%e, %e)\n", p0[0], p0[1]);
                 printf("p1 = (%e, %e)\n", p1[0], p1[1]);
             }
             free(r); 
             free(new2old);
             if ((t >= 0.0) && (t <= 1.0)) { 
                 *i0_touch = new2old[i]; 
                 return touching;
             } 
             else { 
                 return outside;
             }  
         } 
     }  
//   check inside or outside 

     p0 = r;
     p1 = r + 2;
     x0 = p[0] - p0[0];
     y0 = p[1] - p0[1];
     x1 = p[0] - p1[0];
     y1 = p[1] - p1[1];    
     curl0 = x0 * y1 - x1 * y0;
 
     for (i = 1; i < nn; i++) { 
         i1 = (i + 1) % nn;
         p0 = r + (i  + i );
         p1 = r + (i1 + i1);

         x0 = p[0] - p0[0];
         y0 = p[1] - p0[1];
         x1 = p[0] - p1[0];
         y1 = p[1] - p1[1];
         curl = x0 * y1 - x1 * y0;
         if (curl * curl0 < 0.0) { 
             free(r);
             if (new2old) free(new2old);
             return outside; 
         }
     }
     free(r);  
     if (new2old) free(new2old);

     return inside; 
  } 

int line_intersect2(double *linea, double *lineb, double *r, int *touch_only)
{
/***
    returned value: 
    0:      two lines do not intersect.
    1:      two lines intersect at one point r.
    2:      two lines have a common seqment, (r, r2).
***/

    int n, szdim;
    double dxa, dya, dxb, dyb, dxba, dyba, delt0, delta, deltb;
    double ta1, ta2, tb1, tb2, ta, tb;
    double *a1, *a2, *b1, *b2; 
    double amin[2], amax[2], bmin[2], bmax[2]; 
   
    szdim = 2 * sizeof(double);
    *touch_only = 0;

    amin[0] = linea[0];
    amin[1] = linea[1];
    memcpy(amax, amin, (size_t)szdim);
    bmin[0] = lineb[0];
    bmin[1] = lineb[1];
    memcpy(bmax, bmin, (size_t)szdim);
    a1 = linea + 2;
    b1 = lineb + 2;
    for (n = 0; n < 2; n++) { 
        if (amin[n] > a1[n]) { 
            amin[n] = a1[n];
        }
        else if (amax[n] < a1[n]) { 
            amax[n] =  a1[n];
        }
        if (bmin[n] > b1[n]) {
            bmin[n] = b1[n];
        }
        else if (bmax[n] < b1[n]) {
            bmax[n] =  b1[n];
        }
    }
    for (n = 0; n < 2; n++) { 
        if ((amin[n] > bmax[n]) || (amax[n] < bmin[n])) { 
            return 0;
        }
    }  
    a1 = linea;  
    a2 = a1 + 2;
    b1 = lineb;
    b2 = b1 + 2;

    n = 0;

    dxa = a2[0] - a1[0];
    dya = a2[1] - a1[1];
    dxb = b2[0] - b1[0];
    dyb = b2[1] - b1[1];
    dxba = b1[0] - a1[0];
    dyba = b1[1] - a1[1];

    delt0 = dxb * dya - dxa * dyb;
    if (delt0 == 0.0) { 
        if ((a2[0] != a1[0]) && (a2[1] != a1[1])) { 
            ta1 = (b1[0] - a1[0])/(a2[0] - a1[0]);
            ta2 = (b1[1] - a1[1])/(a2[1] - a1[1]);
            if ((ta1 == ta2) && (ta1 >= 0.0) && ta1 <= 1.0) { 
                // b1 is the common pint
                memcpy(r + (n+n), b1, (size_t)szdim);
                n++;
            }
            ta1 = (b2[0] - a1[0])/(a2[0] - a1[0]);
            ta2 = (b2[1] - a1[1])/(a2[1] - a1[1]); 
            if ((ta1 == ta2) && (ta1 >= 0.0) && ta1 <= 1.0) {
                memcpy(r + (n+n), b2, (size_t)szdim);
                n++;
            } 
        }
        if ((n < 2) && (b2[0] != b1[0]) && (b2[1] != b1[1])) {
            tb1 = (a1[0] - b1[0])/(b2[0] - b1[0]);
            tb2 = (a1[1] - b1[1])/(b2[1] - b1[1]);
            if ((tb1 == tb2) && (tb1 >= 0.0) && tb1 <= 1.0) {
                memcpy(r + (n+n), a1, (size_t)szdim);
                n++;
            }
            tb1 = (a2[0] - b1[0])/(b2[0] - b1[0]);
            tb2 = (a2[1] - b1[1])/(b2[1] - b1[1]);
            if ((tb1 == tb2) && (tb1 >= 0.0) && tb1 <= 1.0) {
                memcpy(r + (n+n), a2, (size_t)szdim);
                n++;
            }
        }
        if (n == 1) { 
            *touch_only = 1;
        } 
    }
    else { 
        delta = dxb * dyba - dxba * dyb;
        deltb = dxa * dyba - dxba * dya; 
        ta = delta/delt0;
        tb = deltb/delt0;
        if ((ta < 0.0) || (ta > 1.0) || (tb < 0.0) || (tb > 1.0)) { 
            n = 0;
        }    
        else { 
            r[0] = a1[0] + (a2[0] - a1[0]) * ta;
            r[1] = a1[1] + (a2[1] - a1[1]) * ta;
            n = 1;
            if ((ta == 0.0) || (ta == 1.0) || (tb == 0.0) || (tb == 1.0)) { 
                *touch_only = 1;
            }    
        }  
    } 
    return n;
 } 

int is_node_within_3d_ucell(int if_edge_included, double *pt, int *num_nodes_for_face,
                            int nface, int *nodelist_for_face,
                            double *coords, int nnode)
{
    int dim, i, i1, i2, f, d, n, nn, n0, n1, n2, same, iaminside;
    int *nodelist;
    double *c0, *c1, *c2, cross[3], dc1[3], dc2[3],  norm2, projection;

//    coords of the umesh cell is assumed to be in form
//    (x,y,z)_1, (x,y,z)_2, ..., (x,y,z)_nnode for 3D
//
//    nodelist of each face is either all inward, or all outward.
// 
//     for each node P of c8,  we do the following
//     for each face of the unstructured cell, f,
//         find the perpendicular line that pass through P, and the line intersect with the face at P_0
//         find the vector from P_0 to P,   p0_P
//         find the normal, normal,  of the face through the nodelist of the face for the order
//         find the dot product, P0_P dot normal, of the two vectors.
//
//     If this dot product has the same sign for all the faces, then the node P is inside of the
//     unstructured cell.
//     or
//     If we assume that the normal made from nodelist of each face is always inward, the node P
//     will be outside if this dot product is negative or zero.
//
//     The plane equation of face:   nx(x - x0) + ny(y - y0) + nz(z - z0) = 0;
//     Assume the vertical line passing through P(a,b,c) intersects with the plane at (x, y, z),
//     we have
//         (nx, ny, nz) = q[(a - x), (b - y), (c - z)], here q is the proportional constant,
//     and therefore,
//          x = a - q * nx,   y = b - q * ny,  z = c - q * nz.
//     Substituting these into the plane Eq, we get
//          q = [(nx(a - x0) + ny(b - y0) + nz(c - z0)]/(nx^2 + ny^2 + nz^2).
//     The dot product is
//         (nx, ny, nz) . [(a - x), (b - y), (c - z)] = q * (nx^2 + ny^2 + nz^2)
//         = nx (a - x0) + ny (b - y0) + nz * (c - z0).
//

     dim = 3;

     iaminside = 1;
    
// check whether pt is one of the nodes. If it is, it is NOT within the cell 

     iaminside = 1;
     c0 = coords;
     for (i = 0; i < nnode; i++) {
         same = 1;
         for (d = 0; d < dim; d++) {
             if (pt[d] != c0[d]) {
                 same = 0;
                 break;
             }
         }
         if (same) {
             iaminside = 0;
             return iaminside;
         }
         else { 
             c0 += dim;
         }
     }
     nodelist = nodelist_for_face;
     for (f = 0; f < nface; f++) {

//       find the normal of the face

         norm2 = 0.0;
         nn = num_nodes_for_face[f];
         for (i = 0; i < nn; i++) {
             i1 = (i  + 1) % nn;
             i2 = (i1 + 1) % nn;
             n0 = nodelist[i];
             n1 = nodelist[i1];
             n2 = nodelist[i2];
             c0 = coords + (n0 * dim);
             c1 = coords + (n1 * dim);
             c2 = coords + (n2 * dim);
             for (d = 0; d < dim; d++) {
                 dc1[d] = c1[d] - c0[d];
                 dc2[d] = c2[d] - c1[d];
             }
             cross[0] = dc1[1] * dc2[2] - dc2[1] * dc1[2];
             cross[1] = dc1[2] * dc2[0] - dc2[2] * dc1[0];
             cross[2] = dc1[0] * dc2[1] - dc2[0] * dc1[1];
             norm2 = 0.0;
             for (d = 0; d < dim; d++) {
                 norm2 += (cross[d] * cross[d]);
             }
             if (norm2 > 0.0) break;
         }
         assert(norm2);
         projection = 0.0;
         for (d = 0; d < dim; d++) {
             projection += (cross[d] *(pt[d] - c0[d]));
         }
         if (projection <= 0.0) {
             iaminside = 0;
             break;
         }
         nodelist += nn;
     }
    return iaminside;
 }
///////////////////////////////////////////////////////////////////
 
 
