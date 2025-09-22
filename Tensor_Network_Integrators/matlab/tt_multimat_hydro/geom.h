#ifndef _GEOM_
#define _GEOM_

#ifdef __cplusplus
extern "C" {
#endif

enum point_position {
     outside  = 0,
     inside   = 1,
     touching = 2,
};
typedef enum point_position point_position;
 
void gsph_rec(int ifinquiry, int *mixed,
              int geop, int dim, double *ctr, double r,
              double *xl, double *dx, double *vol);

void poly2d_rec(int ifinquiry, int *mixed, int geop, double vcell, int dim,
                int nnode, double *xys, double *x2min, double *x2max,
                double *xl, double *dx, double *vol);

void rec_rec(int ifinquiry, int *ifmixed, int geop,
             double vcell, int dim, double *xxl, double *xxr,
             double *xl, double *dx, double *vol);

void gconic_rec(int ifinquiry, int *mixed, int ifcyl,
             int dim, double r0, double r1, double *ctr0, double *ctr1,
             double *xl, double *dx, double *vol);

void poly3d_cut_by_plane(int dim, int norm_direction, double plane_location,
                 double *x0_scale, double dx_scale,
                 int nface, int nedge, int nnode,  double *coords1d,
                 int *num_node_for_face, int *nodelist_for_face,
                 int *edgelist_for_face, int *nodelist_for_edge,
                 int *nnode_out, double *coords1d_out,
                 int *nface_for_zone_out, int **nnode_for_face_out,
                 int **nodelist_for_face_out,
                 int *nnode_intface, int *nodelist_intface);
/****
void gpoly_rec(int mype,
               int nreg, int ifinquiry, int *ifmixed, int geop, double vcell,
               int ifcyl, int dim, double *coords, int nnode,
               double *x2min, double *x2max,
               int nedge, int nface, int *nodelist_for_edge,
               int *edgelist_for_face, int *nedge_ea_face,
               int *nodelist_for_face, int *num_nodes_for_face,
               double *ctr, double *norm_of_face,
               double *xl, double *dx, double *vol);

void sph_rec_int(int dim, int geop, double *ctr, double radius,
                 double vcell, double *xl, double *dx, double *vol);

void rec_rec(int ifinquiry, int *ifmixed, int geop,
             double vcell, int dim, double *xxl, double *xxr,
             double *xl, double *dx, double *vol);

****/

point_position point_in_polygon(double *myr, int n, double *myxl, double *myxr,
                                double *p, int *i0_touch); 

int line_intersect2(double *linea, double *lineb, double *r, int *touch_only);

int hull_test(double *myr1, int nn1, double *myr2, int nn2, double **pr, int *n);

int is_node_within_3d_ucell(int if_edge_included, double *pt, int *num_nodes_for_face,
                            int nface, int *nodelist_for_face,
                            double *coords, int nnode); 


#ifdef __cplusplus
}
#endif

#endif
