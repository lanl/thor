#ifdef MPI
#include <mpi.h>
#endif

#include "minip.h"
#include "mesh.h"
#include "util.h"
#include "mat.h"
#include "vof2d.h"
#include "io.h"
#include "eos.h"
#include "hdf5.h"
#include "xdmf.h"
#include "h5aux.h"
#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> // for access()


int      dim_prob     = 0;        // the dimensionality
int      nbdry_prob   = 0;        // the number of layers, typical 1
int      ncell_prob[] = {0,0,0};  // the number of cells three dimension
double   xl_prob[] = {0.0, 0.0, 0.0};  // the left ends of simulation domain, excluding the ghost cells
double   xr_prob[] = {0.0, 0.0, 0.0};  // the right ends of simulation domain, excluding the ghost cells

int      nmat_mesh    = 0;

// mesh variables
 
double **rho_2dcell  = NULL;
double **ei_2dcell   = NULL;  // internal energy density, per volume
double **pres_2dcell = NULL;
double **divu_2dcell = NULL; 
double **qvis_2dcell = NULL;
double **force_2dnode = NULL;
double ***vel_2dnode = NULL;
double ***vav_2dnode = NULL;  // time-averaged velocity for advection

double ***rho_3dcell  = NULL;
double ***ei_3dcell   = NULL;  // internal energy density, per volume
double ***pres_3dcell = NULL;
double ***divu_3dcell = NULL; 
double ***qvis_3dcell = NULL;
double ***force_3dnode = NULL;  
double ****vel_3dnode = NULL;
double ****vav_3dnode = NULL; // time-averaged velocity for advection

double ***vf_2dmat   = NULL; 
double ***rho_2dmat  = NULL;
double ***ei_2dmat   = NULL;  // internal energy density, per volume
double ***pres_2dmat = NULL;

double ****vf_3dmat   = NULL;
double ****rho_3dmat  = NULL;
double ****ei_3dmat   = NULL;  // internal energy density, per volume
double ****pres_3dmat = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////
void mesh_pass_mat(double **xl, double **xr, int **ncell, int *nbdry)
{
     *xl = xl_prob;
     *xr = xr_prob;
     *ncell = ncell_prob;
     *nbdry = nbdry_prob;

     return;
 }

void set_mesh(int dim, const double *xl, const double *xr, const int *ncell, const int nbdry,
          const int nmat)
{
    int i, j, k;
    int  ncell_ext[3], nnode_ext[3], sizes[4];
    long long lsize, offset;
    double *vf, *rho, *ei, *pres; 
    double **vf1d, **rho1d, **ei1d, **pres1d;     // for value per mat  
    double ***vf2d, ***rho2d, ***ei2d, ***pres2d; // for value per mat 
    double ***v3d, **v2d, *v1d, **v2dav, *v1dav;
    double **rc2d, **ec2d, **pc2d, **dc2d, *rc1d, *ec1d, *pc1d, *dc1d, **qc2d, *qc1d; // for value per cell

    assert((dim > 1) && (dim <= 3));

    dim_prob = dim;
    for (i = 0; i < dim; i++) {
        xl_prob[i] = xl[i];
        xr_prob[i] = xr[i];
        ncell_prob[i] = ncell[i];

        ncell_ext[i] = ncell[i] + nbdry + nbdry;
        nnode_ext[i] = ncell_ext[i] + 1;
    }
    nbdry_prob = nbdry;

    nmat_mesh = nmat;

    if (dim == 2) {
        vf_2dmat   = (double ***) malloc(ncell_ext[1] * sizeof(double **));
        rho_2dmat  = (double ***) malloc(ncell_ext[1] * sizeof(double **));
        ei_2dmat   = (double ***) malloc(ncell_ext[1] * sizeof(double **));
        pres_2dmat = (double ***) malloc(ncell_ext[1] * sizeof(double **));

        lsize = ncell_ext[1] * ncell_ext[0];
        vf_2dmat[0]   = (double **) malloc(lsize * sizeof(double *));
        rho_2dmat[0]  = (double **) malloc(lsize * sizeof(double *));
        ei_2dmat[0]   = (double **) malloc(lsize * sizeof(double *));
        pres_2dmat[0] = (double **) malloc(lsize * sizeof(double *));

        for (j = 1; j < ncell_ext[1]; j++) {
            vf_2dmat[j]   = vf_2dmat[j-1]  + ncell_ext[0];
            rho_2dmat[j]  = rho_2dmat[j-1] + ncell_ext[0];
            ei_2dmat[j]   = ei_2dmat[j-1]  + ncell_ext[0];
            pres_2dmat[j] = pres_2dmat[j-1] + ncell_ext[0];
        }
        lsize *= nmat; 
        vf   = (double *) malloc(lsize * sizeof(double));
        rho  = (double *) malloc(lsize * sizeof(double));
        ei   = (double *) malloc(lsize * sizeof(double));
        pres = (double *) malloc(lsize * sizeof(double));

        for (j = 0; j < ncell_ext[1]; j++) { 
            for (i = 0; i < ncell_ext[0]; i++) { 
                vf_2dmat[j][i]  = vf; 
                rho_2dmat[j][i] = rho;
                ei_2dmat[j][i]  = ei; 
                pres_2dmat[j][i] = pres; 

                vf  += nmat;  
                rho += nmat;
                ei  += nmat;
                pres += nmat; 
            }
        }  
        rho_2dcell   = (double **) malloc(ncell_ext[1] * sizeof(double *));
        ei_2dcell    = (double **) malloc(ncell_ext[1] * sizeof(double *));
        pres_2dcell  = (double **) malloc(ncell_ext[1] * sizeof(double *));
        divu_2dcell  = (double **) malloc(ncell_ext[1] * sizeof(double *));
        qvis_2dcell  = (double **) malloc(ncell_ext[1] * sizeof(double *));

        lsize = ncell_ext[1] * ncell_ext[0];       
        rho_2dcell[0]  = (double *) malloc(lsize * sizeof(double));
        ei_2dcell[0]   = (double *) malloc(lsize * sizeof(double));
        pres_2dcell[0] = (double *) malloc(lsize * sizeof(double));
        divu_2dcell[0] = (double *) malloc(lsize * sizeof(double));
        qvis_2dcell[0] = (double *) malloc(lsize * sizeof(double));

        for (j = 1; j < ncell_ext[1]; j++) {
            rho_2dcell[j]  = rho_2dcell[j-1] + ncell_ext[0];
            ei_2dcell[j]   = ei_2dcell[j-1]  + ncell_ext[0];
            pres_2dcell[j] = pres_2dcell[j-1] + ncell_ext[0];
            divu_2dcell[j] = divu_2dcell[j-1] + ncell_ext[0];
            qvis_2dcell[j] = qvis_2dcell[j-1] + ncell_ext[0]; 
        }
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {
                divu_2dcell[j][i] = 0.0;
            }
        }
        lsize = nnode_ext[1] * nnode_ext[0];
        force_2dnode    = (double **) malloc(nnode_ext[1] * sizeof(double *));
        force_2dnode[0] = (double  *) malloc(lsize * sizeof(double)); 
        for (j = 1; j < nnode_ext[1]; j++) { 
            force_2dnode[j] = force_2dnode[j-1] + nnode_ext[0]; 
        } 
        vel_2dnode = (double ***) malloc(nnode_ext[1] * sizeof(double **));
        vav_2dnode = (double ***) malloc(nnode_ext[1] * sizeof(double **));
       
        v2d   = (double **) malloc(lsize * sizeof(double *));
        v2dav = (double **) malloc(lsize * sizeof(double *));
        v1d   = (double  *) malloc(lsize * dim * sizeof(double));
        v1dav = (double  *) malloc(lsize * dim * sizeof(double));

        offset  = 0;
        for (j = 0; j < nnode_ext[1]; j++) {
            v2d[0] = v1d + offset;
            v2dav[0] = v1dav + offset;
            for (i = 1; i < nnode_ext[0]; i++) {
                v2d[i] = v2d[i-1] + dim;
                v2dav[i] = v2dav[i-1] + dim;
            }
            vel_2dnode[j] = v2d;
            vav_2dnode[j] = v2dav;
            offset += (nnode_ext[0] * dim);
            v2d    +=  nnode_ext[0];
            v2dav  += nnode_ext[0];
        }
    }
    else if (dim == 3) {
        vf_3dmat   = (double ****) malloc(ncell_ext[2] * sizeof(double ***));
        rho_3dmat  = (double ****) malloc(ncell_ext[2] * sizeof(double ***));
        ei_3dmat   = (double ****) malloc(ncell_ext[2] * sizeof(double ***));
        pres_3dmat = (double ****) malloc(ncell_ext[2] * sizeof(double ***));

        rho_3dcell  = (double ***) malloc(ncell_ext[2] * sizeof(double **));
        ei_3dcell   = (double ***) malloc(ncell_ext[2] * sizeof(double **));
        pres_3dcell = (double ***) malloc(ncell_ext[2] * sizeof(double **));
        divu_3dcell = (double ***) malloc(ncell_ext[2] * sizeof(double **));
        qvis_3dcell = (double ***) malloc(ncell_ext[2] * sizeof(double **));

        lsize  = ncell_ext[1] * ncell_ext[2];

        vf2d   = (double ***) malloc(lsize * sizeof(double **));
        rho2d  = (double ***) malloc(lsize * sizeof(double **));
        ei2d   = (double ***) malloc(lsize * sizeof(double **));
        pres2d = (double ***) malloc(lsize * sizeof(double **));

        rc2d  = (double **) malloc(lsize * sizeof(double **));
        ec2d  = (double **) malloc(lsize * sizeof(double **));
        pc2d  = (double **) malloc(lsize * sizeof(double **));
        dc2d  = (double **) malloc(lsize * sizeof(double **));
        qc2d  = (double **) malloc(lsize * sizeof(double **));

        lsize *= ncell_ext[0];
        vf1d   = (double **) malloc(lsize * sizeof(double *));
        rho1d  = (double **) malloc(lsize * sizeof(double *));
        ei1d   = (double **) malloc(lsize * sizeof(double *));
        pres1d = (double **) malloc(lsize * sizeof(double *));

        rc1d = (double *) malloc(lsize * sizeof(double));
        ec1d = (double *) malloc(lsize * sizeof(double));
        pc1d = (double *) malloc(lsize * sizeof(double));
        dc1d = (double *) malloc(lsize * sizeof(double));
        qc1d = (double *) malloc(lsize * sizeof(double));

        offset = 0;
        for (k = 0; k < ncell_ext[2]; k++) {
            vf2d[0]   = vf1d   + offset;
            rho2d[0]  = rho1d  + offset;
            ei2d[0]   = ei1d   + offset;
            pres2d[0] = pres1d + offset;

            rc2d[0] = rc1d + offset;
            ec2d[0] = ec1d + offset;
            pc2d[0] = pc1d + offset;
            dc2d[0] = dc1d + offset;
            qc2d[0] = qc1d + offset; 

            for (j = 1; j < ncell_ext[1]; j++) {
                vf2d[j]   =   vf2d[j-1] + ncell_ext[0]; 
                rho2d[j]  =  rho2d[j-1] + ncell_ext[0];
                ei2d[j]   =   ei2d[j-1] + ncell_ext[0];
                pres2d[j] = pres2d[j-1] + ncell_ext[0];

                rc2d[j] = rc2d[j-1] + ncell_ext[0];
                ec2d[j] = ec2d[j-1] + ncell_ext[0];
                pc2d[j] = pc2d[j-1] + ncell_ext[0];
                dc2d[j] = dc2d[j-1] + ncell_ext[0];
                qc2d[j] = qc2d[j-1] + ncell_ext[0]; 
            }
            vf_3dmat[k]   = vf2d;
            rho_3dmat[k]  = rho2d;
            ei_3dmat[k]   = ei2d;
            pres_3dmat[k] = pres2d;

            rho_3dcell[k]  = rc2d;
            ei_3dcell[k]   = ec2d;
            pres_3dcell[k] = pc2d;
            divu_3dcell[k] = dc2d; 
            qvis_3dcell[k] = qc2d;

            offset += (ncell_ext[0] * ncell_ext[1]);
            vf2d   += ncell_ext[1];
            rho2d  += ncell_ext[1];
            ei2d   += ncell_ext[1];
            pres2d += ncell_ext[1];

            rc2d += ncell_ext[1];
            ec2d += ncell_ext[1];
            pc2d += ncell_ext[1];
            dc2d += ncell_ext[1];
            qc2d += ncell_ext[1]; 
        }
        for (k = 0; k < ncell_ext[2]; k++) { 
            for (j = 0; j < ncell_ext[1]; j++) { 
                for (i = 0; i < ncell_ext[0]; i++) { 
                    divu_3dcell[k][j][i] = 0.0;
                }
            }
        } 
        lsize *= nmat; 
        vf   = (double *) malloc(lsize * sizeof(double));
        rho  = (double *) malloc(lsize * sizeof(double));
        ei   = (double *) malloc(lsize * sizeof(double));
        pres = (double *) malloc(lsize * sizeof(double));

        for (k = 0; k < ncell_ext[2]; k++) { 
            for (j = 0; j < ncell_ext[1]; j++) {
                for (i = 0; i < ncell_ext[0]; i++) {  
                    vf_3dmat[k][j][i]  = vf; 
                    rho_3dmat[k][j][i] = rho; 
                    ei_3dmat[k][j][i]  = ei;
                    pres_3dmat[k][j][i] = pres;

                    vf   += nmat;
                    rho  += nmat;
                    ei   += nmat;
                    pres += nmat; 
                }
            }
        }
        sizes[0] = dim;
        for (i = 0; i < dim; i++) {
            sizes[1+i] = nnode_ext[i];
        }
        vel_3dnode = (double ****) malloc(nnode_ext[2] * sizeof(double ***));
        v3d = NULL;
        v2d = NULL;
        v1d = NULL;
        ASSIGN_4D_FORM(double, sizes, vel_3dnode, v3d, v2d, v1d);

        vav_3dnode = (double ****) malloc(nnode_ext[2] * sizeof(double ***));
        v3d = NULL;
        v2d = NULL;
        v1d = NULL;
        ASSIGN_4D_FORM(double, sizes, vav_3dnode, v3d, v2d, v1d);

        force_3dnode = (double ***) malloc(nnode_ext[2] * sizeof(double **));
        lsize = nnode_ext[2] * nnode_ext[1];
        v2d = (double **) malloc(lsize * sizeof(double *)); 
        lsize *= nnode_ext[0];
        v1d = (double  *) malloc(lsize * sizeof(double));  

        offset = 0;
        for (k = 0; k < nnode_ext[2]; k++) {
            v2d[0] = v1d + offset;
            for (j = 1; j < nnode_ext[1]; j++) {
                v2d[j] = v2d[j-1] + nnode_ext[0];
            }
            force_3dnode[k] = v2d;
            offset += (nnode_ext[0] * nnode_ext[1]);
            v2d    +=  nnode_ext[1];
        }
    }

    return;
}

void set_mesh_mat(int dim, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                  int nmat, int *solid_num_ea_mat, double *gamma_ea_mat,
                  int nreg, int *reg2matids,
                  Region_Shape *reg_shape,
                  double *rho_ea_reg, double *pres_ea_reg, double *ei_ea_reg, double **v_ea_reg)
{
     int ncycle;
     double t;

     int fileid; 
     int i, j, k, m;
     int ncell_ext[3], sizes[4];
     long long lsize;
     double dx[3];
     double ***vel_2dcell, ****vel_3dcell, *v1d, **v2d, ***v3d;

     ncycle = 0;

     lsize = 1;
     for (i = 0; i < dim; i++) {
         dx[i] = (xr_prob[i] - xl_prob[i])/(double)ncell_prob[i];
         ncell_ext[i] = ncell_prob[i] + nbdry_prob + nbdry_prob;
         lsize *= ncell_ext[i];
     }
     vel_2dcell = NULL;
     vel_3dcell = NULL;
     if (dim == 2) {
         sizes[0] = dim;
         sizes[1] = ncell_ext[0];
         sizes[2] = ncell_ext[1];
         vel_2dcell = (double ***) malloc(ncell_ext[1] * sizeof(double **));
         v2d            = (double  **) malloc(lsize * sizeof(double *));
         v1d            = (double   *) malloc(lsize * dim * sizeof(double));
         ASSIGN_3D_FORM(double, vel_2dcell, v2d, v1d, sizes);
     }
     else if (dim == 3) {
         sizes[0] = dim;
         sizes[1] = ncell_ext[0];
         sizes[2] = ncell_ext[1];
         sizes[3] = ncell_ext[2];
         vel_3dcell  = (double ****) malloc(ncell_ext[2] * sizeof(double ***));
         v3d             = (double  ***) malloc(ncell_ext[2] * ncell_ext[1] * sizeof(double **));
         v2d             = (double   **) malloc(lsize * sizeof(double *));
         v1d             = (double    *) malloc(lsize * dim * sizeof(double));
         ASSIGN_4D_FORM(double, sizes, vel_3dcell, v3d, v2d, v1d);
     }
     set_mat(dim, ncell_prob, xl_prob, dx, nbdry_prob, 
             btype_lower, btype_upper,
             nmat, solid_num_ea_mat, gamma_ea_mat,
             nreg, reg2matids, reg_shape,
             rho_ea_reg, pres_ea_reg, ei_ea_reg, v_ea_reg,
             vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat, 
             rho_2dcell, ei_2dcell, pres_2dcell, vel_2dcell,
             vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat, 
             rho_3dcell, ei_3dcell, pres_3dcell, vel_3dcell);

     if (dim == 2) { 
//       get node velocity

         get_node_vel_2d(ncell_prob, nbdry_prob, rho_2dcell, vel_2dcell, vel_2dnode);
         bdry_node_2d(ncell_prob, nbdry_prob, btype_lower, btype_upper, vel_2dnode);
     }
     else if (dim == 3) { 
//       get node velocity

         get_node_vel_3d(ncell_prob, nbdry_prob, rho_3dcell, vel_3dcell, vel_3dnode);
         bdry_node_3d(ncell_prob, nbdry_prob, btype_lower, btype_upper, vel_3dnode);
     }
#ifdef MESHIO
//   write
     t = 0.0;
     ncycle = 0;

     create_file("file_initial", t, ncycle, &fileid);
     write_mesh_mat(fileid, "mesh", dim, xl_prob, xr_prob, ncell_prob, nbdry_prob,
                nmat_mesh,
                vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat, rho_2dcell, ei_2dcell, pres_2dcell, NULL, vel_2dnode,
                vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat, rho_3dcell, ei_3dcell, pres_3dcell, NULL, vel_3dnode);

     close_file(&fileid);
#endif 

     if (dim == 2) {
         free(vel_2dcell[0][0]);
         free(vel_2dcell[0]);
         free(vel_2dcell);
     }
     else if (dim == 3) {
         free(vel_3dcell[0][0][0]);
         free(vel_3dcell[0][0]);
         free(vel_3dcell[0]);
         free(vel_3dcell);
     }
     return;
 }


void cal_cell_zgrad3d(double *dx, double var[3][3][3], double *gradv)
{
     int i, j, k;
     double v0, dfsum;

     v0 = var[1][1][1];
     dfsum = 0.0;
     for (k = 0; k < 3; k++) {
         for (j = 0; j < 3; j++) {
             dfsum += (var[k][j][2] - v0);
             dfsum -= (var[k][j][0] - v0);
         }
     }
     gradv[0] = dfsum/(18.0 * dx[0]);

     dfsum = 0.0;
     for (k = 0; k < 3; k++) {
         for (i = 0; i < 3; i++) {
             dfsum += (var[k][2][i] - v0);
             dfsum -= (var[k][0][i] - v0);
         }
     }
     gradv[1] = dfsum/(18.0 * dx[1]);

     dfsum = 0.0;
     for (j = 0; j < 3; j++) {
         for (i = 0; i < 3; i++) {
             dfsum += (var[2][j][i] - v0);
             dfsum -= (var[0][j][i] - v0);
         }
     }
     gradv[2] = dfsum/(18.0 * dx[2]);

     return;
}


void cal_cell_zgrad2d(double *dx, double var[3][3], double *gradv)
{
    int i, j;
    double v0, dfsum;

    v0 = var[1][1];

    dfsum = 0.0;
    for (j = 0; j < 3; j++) {
        dfsum += (var[j][2] - v0);
        dfsum -= (var[j][0] - v0);
    }
    gradv[0] = dfsum/(6.0 * dx[0]);
    dfsum = 0.0;
    for (i = 0; i < 3; i++) {
        dfsum += (var[2][i] - v0);
        dfsum -= (var[0][i] - v0);
    }
    gradv[1] = dfsum/(6.0 * dx[1]);

    return;
}


void mesh_sound_speed(int dim, int *ncell, int nbdry, int nmat, 
                      int *solid_num_ea_mat, double *gamma_ea_mat,
                      double **cs_2dcell,  double ***cs_3dcell)
{
    int i, j, k, m;
    int ncell_ext[3];
    double cs; 

    for (i = 0; i < dim; i++) {
        ncell_ext[i] = ncell_prob[i] + nbdry_prob + nbdry_prob;
    }
    if (dim == 2) {
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[0]; i++) {
                cs_2dcell[j][i] = 0.0; 

                for (m = 0; m <  nmat; m++) { 
                    if (vf_2dmat[j][i][m] > 0.0) { 
                        if (!solid_num_ea_mat[m]) { 
                            cs = sqrt(gamma_ea_mat[m] * pres_2dmat[j][i][m]/
                                     (rho_2dmat[j][i][m] + tiny)); 
                        }
                        else { 
                            sspd_solid(solid_num_ea_mat[m], rho_2dmat[j][i][m], &cs); 
                        }
                        cs_2dcell[j][i] += (vf_2dmat[j][i][m] * cs);
                    }
                }  
            }
        }
    }
    else if (dim == 3) {
        for (k = 0; k < ncell_ext[2]; k++) {
            for (j = 0; j < ncell_ext[1]; j++) {
                for (i = 0; i < ncell_ext[0]; i++) {
                    cs_3dcell[k][j][i] = 0.0;

                    for (m = 0; m <  nmat; m++) {
                        if (vf_3dmat[k][j][i][m] > 0.0) {
                            if (!solid_num_ea_mat[m]) {
                                cs = sqrt(gamma_ea_mat[m] * pres_3dmat[k][j][i][m]/
                                     (rho_3dmat[k][j][i][m] + tiny));
                            }
                            else {
                                sspd_solid(solid_num_ea_mat[m], rho_3dmat[k][j][i][m], &cs);
                            }
                            cs_3dcell[k][j][i] += (vf_3dmat[k][j][i][m] * cs);
                        } 
                    }
                }
            }
        }
    }
    return;
}


void get_node_vel_2d(int *ncell, int nbdry, double **rho_for_cell, double ***vel_for_cell, double ***vel_for_node)
{

     int dim, i, j, ic, jc, dir;
     int nnode_last[2];
     double mass, massinv;

     dim = 2;
     for (dir = 0; dir < dim; dir++) {
         nnode_last[dir] = ncell[dir] + nbdry;
     }
     for (j = nbdry; j <= nnode_last[1]; j++) {
         for (i = nbdry; i <= nnode_last[0]; i++)  {
             mass = 0.0;
             for (dir = 0; dir < dim; dir++) {
                 vel_for_node[j][i][dir] = 0.0;
             }
             for (jc = j-1; jc <= j; jc++) {
                 for (ic = i - 1; ic <= i; ic++) {
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_node[j][i][dir] += (rho_for_cell[jc][ic] * vel_for_cell[jc][ic][dir]);
                     }
                     mass += rho_for_cell[jc][ic];
                 }
             }
             massinv = 1.0/(mass + tiny);
             for (dir = 0; dir < dim; dir++) {
                 vel_for_node[j][i][dir] *= massinv;
             }
         }
     }
     return;
  }


void get_node_vel_3d(int *ncell, int nbdry, double ***rho_for_cell, double ****vel_for_cell, double ****vel_for_node)
{
     int dim, i, j, k, ic, jc, kc, dir;
     int nnode_last[3];
     double mass, massinv;

     dim = 3;
     for (dir = 0; dir < dim; dir++) {
         nnode_last[dir] = ncell[dir] + nbdry;
     }
     for (k = nbdry; k <= nnode_last[2]; k++) {
         for (j = nbdry; j <= nnode_last[1]; j++) {
             for (i = nbdry; i <= nnode_last[0]; i++)  {
                 mass = 0.0;
                 for (dir = 0; dir < dim; dir++) {
                     vel_for_node[k][j][i][dir] = 0.0;
                 }
                 for (kc = k-1; kc <= k; kc++) {
                     for (jc = j-1; jc <= j; jc++) {
                         for (ic = i-1; ic <= i; ic++) {
                             for (dir = 0; dir < dim; dir++) {
                                 vel_for_node[k][j][i][dir] += (rho_for_cell[k][jc][ic] * vel_for_cell[kc][jc][ic][dir]);
                             }
                             mass += rho_for_cell[kc][jc][ic];
                         }
                     }
                 }
                 massinv  = 1.0/(mass + tiny);
                 for (dir = 0; dir < dim; dir++) {
                     vel_for_node[k][j][i][dir] *= massinv;
                 }
             }
         }
     }
     return;
 } //get_node_vel_3d


////////////////////////////////////////////////////////////////////////////////////////////////
void mesh_pass_mesh_data(double ***rho_ea_2dcell,  double ***ei_ea_2dcell,   double ***pres_ea_2dcell, 
                         double ***divu_ea_2dcell, double ***qvis_ea_2dcell,
                         double ***force_ea_2dnode, double ****vel_ea_2dnode, double ****vav_ea_2dnode,
                         double ****vf_ea_2dmat, double ****rho_ea_2dmat, 
                         double ****ei_ea_2dmat, double ****pres_ea_2dmat, 
                         double ****rho_ea_3dcell, double ****ei_ea_3dcell,   double ****pres_ea_3dcell,
                         double ****divu_ea_3dcell, double ****qvis_ea_3dcell, 
                         double ****force_ea_3dnode, double *****vel_ea_3dnode, double *****vav_ea_3dnode,
                         double *****vf_ea_3dmat, double *****rho_ea_3dmat,  
                         double *****ei_ea_3dmat, double *****pres_ea_3dmat) 
{
     *rho_ea_2dcell  = rho_2dcell;
     *ei_ea_2dcell   = ei_2dcell;
     *pres_ea_2dcell = pres_2dcell;
     *divu_ea_2dcell = divu_2dcell; 
     *qvis_ea_2dcell = qvis_2dcell; 
     *force_ea_2dnode = force_2dnode; 
     *vel_ea_2dnode = vel_2dnode;
     *vav_ea_2dnode = vav_2dnode;

     *vf_ea_2dmat    = vf_2dmat;
     *rho_ea_2dmat   = rho_2dmat;
     *ei_ea_2dmat    = ei_2dmat;
     *pres_ea_2dmat  = pres_2dmat; 
    
     *rho_ea_3dcell  = rho_3dcell;
      *ei_ea_3dcell   = ei_3dcell;
     *pres_ea_3dcell = pres_3dcell;
     *divu_ea_3dcell = divu_3dcell; 
     *qvis_ea_3dcell = qvis_3dcell;  
     *force_ea_3dnode = force_3dnode; 
     *vel_ea_3dnode  = vel_3dnode;
     *vav_ea_3dnode  = vav_3dnode;

     *vf_ea_3dmat    = vf_3dmat;
     *rho_ea_3dmat   = rho_3dmat;
     *ei_ea_3dmat    = ei_3dmat;
     *pres_ea_3dmat  = pres_3dmat;

     return;
 }

void xdmf_dump(const int ncycle, const double timestamp)
{
    int i, j, k, Nlines;
    const char *fname_h5 = "cs_output.h5";
    const char *fname_xdmf = "cs_output.xdmf";
    char timestep_group_name[20];
    char dataset_name[20];
    hid_t fid_h5, group_id;
    FILE *fid_xdmf;
    char **lines;

    int dim, ncell[3], nbdry;
    double xl[3], xr[3];
    int ncell_ext[3], nnode_ext[3];

    ///////////////////////////////
    dim = dim_prob;
    nbdry = nbdry_prob;
    for (i = 0; i < dim; ++i)
    {
        ncell[i] = ncell_prob[i];
        xl[i] = xl_prob[i];
        xr[i] = xr_prob[i];
        ncell_ext[i] = ncell[i] + 2 * nbdry;
        nnode_ext[i] = ncell_ext[i] + 1;
    }

    if (ncycle == 0)
    {
        // create a new HDF5 file; delete if the file already exists
        delete_file_if_exists(fname_h5);
        fid_h5 = h5_new_file(fname_h5);
        H5Fclose(fid_h5);

        // create a new XDMF file
        xdmf_new_file(fname_xdmf);
        printf("Files %s and %s created successfully\n", fname_h5, fname_xdmf);
    }

    // open an HDF5 group to write the current cycle
    sprintf(timestep_group_name, "Step#%d", ncycle);
    fid_h5 = h5_open_existing_rdwr(fname_h5);
    group_id = h5_open_group(fid_h5, timestep_group_name);

    // open XDMF and write mesh description
    fid_xdmf = xdmf_open_file_append(fname_xdmf);
    xdmf_write_mesh(fid_xdmf, dim, xl, xr, ncell, nbdry);

    if (dim == 2)
    {
        h5_write_2d(group_id, "rho", rho_2dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "rho", CELL_VARIABLE);

        h5_write_2d(group_id, "ei", ei_2dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "ei", CELL_VARIABLE);

        h5_write_2d(group_id, "pres", pres_2dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "pres", CELL_VARIABLE);

        // copy_velocity_component_2d(vel_2dnode, &var_for_2dnode, 0, nnode_ext);
        // h5_write_2d(group_id, "vx", var_for_2dnode, nnode_ext);
        // xdmf_write_dataset(fid_xdmf, dim, nnode_ext, xl, xr, fname_h5,
        //                    timestep_group_name, "vx", NODE_VARIABLE);

        // copy_velocity_component_2d(vel_2dnode, &var_for_2dnode, 1, nnode_ext);
        // h5_write_2d(group_id, "vy", var_for_2dnode, nnode_ext);
        // xdmf_write_dataset(fid_xdmf, dim, nnode_ext, xl, xr, fname_h5,
        //                    timestep_group_name, "vy", NODE_VARIABLE);

        /***
               // pass the data from mat
               int    nmixcell, nmixcell_int, nmixcell_mpoly;
               int    **ijk_in_mixcell, *nmat_in_mixcell, **matid_in_mixcell;
               double **vf_in_mixcell, **rho_in_mixcell, **pres_in_mixcell, **ei_in_mixcell;
               mat_pass_mix_data(&nmixcell, &nmixcell_int, &nmixcell_mpoly,
                                 &nmat_in_mixcell, &ijk_in_mixcell, &matid_in_mixcell,
                                 &vf_in_mixcell, &rho_in_mixcell,
                                 &pres_in_mixcell, &ei_in_mixcell);

               int lsize = 0;
               for (int mx = 0; mx < nmixcell; mx++)
                   lsize += nmat_in_mixcell[mx];
               h5_write_1d_int(group_id, "nmixcell", &nmixcell, 1);
               h5_write_1d_int(group_id, "nmixcell_int", &nmixcell_int, 1);
               h5_write_1d_int(group_id, "nmixcell_mpoly", &nmixcell_mpoly, 1);
               h5_write_1d_int(group_id, "lsize", &lsize, 1);
               if (nmixcell > 0) {
                   h5_write_1d_int(group_id, "nmat_in_mixcell", nmat_in_mixcell, nmixcell);
                   h5_write_1d_int(group_id, "matid_in_mixcell", &(matid_in_mixcell[0][0]), lsize);
                   h5_write_1d_int(group_id, "ijk_in_mixcell", &(ijk_in_mixcell[0][0]), dim*nmixcell);
                   h5_write_1d(group_id, "vf_in_mixcell", &(vf_in_mixcell[0][0]), lsize);
                   h5_write_1d(group_id, "rho_in_mixcell", &(rho_in_mixcell[0][0]), lsize);
                   h5_write_1d(group_id, "ei_in_mixcell", &(ei_in_mixcell[0][0]), lsize);
                   h5_write_1d(group_id, "pres_in_mixcell", &(pres_in_mixcell[0][0]), lsize);
               }
        ***/
    }
    else if (dim == 3)
    { // TODO: var_for_3dnode is not properly allocated
        h5_write_3d(group_id, "rho", rho_3dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "rho", CELL_VARIABLE);

        h5_write_3d(group_id, "ei", ei_3dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "ei", CELL_VARIABLE);

        h5_write_3d(group_id, "pres", pres_3dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "pres", CELL_VARIABLE);

        // copy_velocity_component_3d(vel_3dnode, &var_for_3dnode, 0, nnode_ext);
        // h5_write_3d(group_id, "vx", var_for_3dnode, nnode_ext);
        // xdmf_write_dataset(fid_xdmf, dim, nnode_ext, xl, xr, fname_h5,
        //                    timestep_group_name, "vx", NODE_VARIABLE);

        // copy_velocity_component_3d(vel_3dnode, &var_for_3dnode, 1, nnode_ext);
        // h5_write_3d(group_id, "vy", var_for_3dnode, nnode_ext);
        // xdmf_write_dataset(fid_xdmf, dim, nnode_ext, xl, xr, fname_h5,
        //                    timestep_group_name, "vy", NODE_VARIABLE);

        // copy_velocity_component_3d(vel_3dnode, &var_for_3dnode, 2, nnode_ext);
        // h5_write_3d(group_id, "vz", var_for_3dnode, nnode_ext);
        // xdmf_write_dataset(fid_xdmf, dim, nnode_ext, xl, xr, fname_h5,
        //                    timestep_group_name, "vz", NODE_VARIABLE);
    }
    H5Gclose(group_id);
    H5Fclose(fid_h5);

    xdmf_write_timestamp(fid_xdmf, timestamp);
    xdmf_close_group(fid_xdmf);
    xdmf_close_file(fid_xdmf);

    ///////////////////////////////
}

void xdmf_dump_basic(const char *basename, const int ncycle, const double timestamp)
{
    int i;
    char fname_h5[128], fname_xdmf[128], timestep_group_name[64];
    hid_t fid_h5, group_id;
    FILE *fid_xdmf;

    int dim = dim_prob;
    int ncell[3], ncell_ext[3], nnode_ext[3];
    double xl[3], xr[3];

    // Construct file names
    sprintf(fname_h5, "%s.h5", basename);
    sprintf(fname_xdmf, "%s.xdmf", basename);

    for (i = 0; i < dim; ++i)
    {
        ncell[i] = ncell_prob[i];
        xl[i] = xl_prob[i];
        xr[i] = xr_prob[i];
        ncell_ext[i] = ncell[i] + 2 * nbdry_prob;
        nnode_ext[i] = ncell_ext[i] + 1;
    }

    if (ncycle == 0)
    {
        delete_file_if_exists(fname_h5);
        fid_h5 = h5_new_file(fname_h5);
        H5Fclose(fid_h5);

        xdmf_new_file(fname_xdmf);
        printf("Created files: %s and %s\n", fname_h5, fname_xdmf);
    }

    sprintf(timestep_group_name, "Step#%d", ncycle);
    fid_h5 = h5_open_existing_rdwr(fname_h5);
    group_id = h5_open_group(fid_h5, timestep_group_name);

    fid_xdmf = xdmf_open_file_append(fname_xdmf);
    xdmf_write_mesh(fid_xdmf, dim, xl, xr, ncell, nbdry_prob);

    if (dim == 2)
    {
        h5_write_2d(group_id, "rho", rho_2dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "rho", CELL_VARIABLE);

        h5_write_2d(group_id, "ei", ei_2dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "ei", CELL_VARIABLE);

        h5_write_2d(group_id, "pres", pres_2dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "pres", CELL_VARIABLE);
    }
    else if (dim == 3)
    {
        h5_write_3d(group_id, "rho", rho_3dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "rho", CELL_VARIABLE);

        h5_write_3d(group_id, "ei", ei_3dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "ei", CELL_VARIABLE);

        h5_write_3d(group_id, "pres", pres_3dcell, ncell_ext);
        xdmf_write_dataset(fid_xdmf, dim, ncell_ext, xl, xr, fname_h5,
                           timestep_group_name, "pres", CELL_VARIABLE);
    }

    H5Gclose(group_id);
    H5Fclose(fid_h5);

    xdmf_write_timestamp(fid_xdmf, timestamp);
    xdmf_close_group(fid_xdmf);
    xdmf_close_file(fid_xdmf);
}

