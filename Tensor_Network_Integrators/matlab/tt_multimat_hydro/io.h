#ifndef _IO_
#define _IO_

#ifdef __cplusplus
extern "C" {
#endif

#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
// #include "mio.h" 
#include "minip.h"
  void write_2Dllarray_to_ascii(const char *filename, long long **int_2dcell, int *ncell_ext);
  void write_2Dmat_variables(const char *suffix, double ***pres_2dmat, double ***ei_2dmat, double ***rho_2dmat, double ***vf_2dmat, int ny, int nx, int nz, int ncycle);
  void write_2Dintarray_to_ascii(const char *filename, int **int_2dcell, int *ncell_ext);

  void write_2Dcell_variables(const char *suffix, double **pres_2dcell, double **ei_2dcell, double **rho_2dcell, int ny, int nx, int nz, int ncycle);
void create_file(char *filename, double t, int ncycle, int *fileid);
void close_file(int *fileid); 

  void write_3Darray_as_matrix(const char *filename, double ***array, int n, int m, int l);

  void write_3Darray_to_ascii(const char *filename, double ***array, int rows, int cols, int depth);

  void write_2Darray_to_ascii(const char *filename, double **rho_2dcell, int *ncell_ext);

  void write_array_to_bin(const char *filename, void *array, size_t element_size, int size);

  /////////
  void create_file(char *filename, double t, int ncycle, int *fileid);

  void close_file(int *fileid);

void viz_dump(const char *probname, int *fileid,
                int dim, double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
                int nmat,
                double ***vf_2dmat, double ***rho_2dmat, double ***ei_2dmat, double ***pres_2dmat,
                double **rho_2dcell, double **ei_2dcell, double **pres_2dcell, double **divu_2dcell, double ***vel_2dnode,
                double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat, double ****pres_3dmat,
                double ***rho_3dcell, double ***ei_3dcell, double ***pres_3dcell, double ***divu_3dcell, double ****vel_3dnode);

void write_mesh_mat(int fileid, char *meshname,
                      int dim, double *xl_prob, double *xr_prob, int *ncell, int nbdry,
                      int nmat_prob,
                      double ***vf_2dmat, double ***rho_2dmat, double ***ei_2dmat, double ***pres_2dmat,
                      double **rho_2dcell, double **ei_2dcell, double **pres_2dcell, double **divu_2dcell, double ***vel_2dnode,
                      double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat, double ****pres_3dmat,
                      double ***rho_3dcell, double ***ei_3dcell, double ***pres_3dcell, double ***divu_3dcell, double ****vel_3dnode);

// void write_scalar(char *varname, int ghost_cell_included, mio_Data_Type datatype, mio_Mesh_Var_Type vartype,
//                  int fileid, int meshid, int dim, int *ncell, int nbdry,
//                  void *var1d, int *varid);

// void write_vecter(char *varname, int ghost_cell_included, mio_Data_Type datatype, mio_Mesh_Var_Type vartype,
//                  int fileid, int meshid, int dim, int *ncell, int nbdry,
//                  double ***var_for_2d, double ****var_for_3d, int *varid);

void write_mpoly(int fileid, int dim, int pass,
                   int nmat,
                   int nmixcell_mpoly, int *nmat_for_mixcell, int **matids_for_mixcell,
                   int *nnode_in_mixcell, double **coords_in_mixcell,
                   int **nnode_for_minterface_in_mixcell, int ***nodes_for_minterface_in_mixcell,
                   int **nnode_for_mpoly_in_mixcell, int ***nodes_for_mpoly_in_mixcell,         // for 2D
                   int **nface_for_mpoly_in_mixcell, int ***nnode_for_face_ea_mpoly_in_mixcell, // for 3D
                   int ***nodelist_for_face_ea_mpoly_in_mixcell);

// void write_a_polygon(int fid, char *name, int dim, int nnode_tot, double *coord,
//                     int nnode, int *nodelist_for_face);

#ifdef __cplusplus
}
#endif
#endif
