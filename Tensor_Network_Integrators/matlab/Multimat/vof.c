#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>
#include <time.h>
#include "mesh.h" 
#include "mat.h" 
#include "vof3d.h"
#include "vof2d.h"
#include "util.h" 
#include "io.h"

static double vfmin = 1.0e-06;

void advected_vol2d(double *xl, double *dx, int direction, double advected_distance,
                    int nnode_mat, double *coords_mat, int nmpoly, int *matids,
                    int *nnode_for_face_mat, int **nodes_for_face_mat,
                    int nmat, double *vols_mat);

void advected_vol3d(double *xl, double *dx, int direction, double advected_distance,
                    int nnode_mat, double *coords_mat, int nmpoly, int *matids,
                    int *nface_for_mpoly_mat, int **nnode_for_face_mat,
                    int **nodelist_for_face_mat,
                    int nmat, double *vols_mat); 

void vof_init(int dim, double *xl_prob, double *xr_prob, int *ncell_prob, int nbdry_prob, int nmat, 
              double t, int ncycle, double *vf1d, double *rho1d, double *pres1d, double *ei1d)
{
    char filename[32]; 
    int i, j, k, s, nmix, m, nm, lsize, offset, fileid, meshid; 
    int ncell_ext[3], sizes[4]; 
    long long llsize;
    
    int *matids; 
    int **nmat_for_2dcell, **matid_for_2dcell, **mixcell_for_2dcell;
    int ***nmat_for_3dcell, ***matid_for_3dcell, ***mixcell_for_3dcell;  
    int *nmat_for_mixcell, **ijk_for_mixcell, **matids_for_mixcell;
    int **ibuffer2d, *ibuffer1d; 
    double dx[3]; 
    double ***vf2d, ****vf3d, ***rho2d, ****rho3d, ***pres2d, ****pres3d, ***ei2d, ****ei3d; 
    double **vf_for_mixcell, **rho_for_mixcell, **pres_for_mixcell, **ei_for_mixcell, *vfs; 
    double ***buffer3d, **buffer2d, *buffer1d; 

    int *ncell, nbdry, ifmesh_initialized_for_mat;
    double *xl, *xr; 
    
    mesh_pass_mat(&xl, &xr, &ncell, &nbdry,
                 &nmat_for_2dcell, &matid_for_2dcell, &mixcell_for_2dcell,
                 &nmat_for_3dcell, &matid_for_3dcell, &mixcell_for_3dcell);

    ifmesh_initialized_for_mat = 1; // assume nmat_for_2dcell etc are allocated. 

    llsize = 1; 
    for (i = 0; i < dim; i++) { 
        ncell_ext[i] = ncell_prob[i] + nbdry_prob + nbdry_prob;
        sizes[i+1] = ncell_ext[i]; 
        llsize *= ncell_ext[i]; 
    }
    matids = (int *) malloc(nmat * sizeof(int));
    for (i = 0; i < nmat; i++) { 
        matids[i] = i;
    } 
    sizes[0] = nmat;
    if (dim == 2) { 
        vf2d     = (double ***) malloc(sizes[2] * sizeof(double **));
        buffer2d = (double  **) malloc(sizes[2] * sizes[1] * sizeof(double *));
        ASSIGN_3D_FORM(double,  vf2d, buffer2d, vf1d, sizes);

	rho2d    = (double ***) malloc(sizes[2] * sizeof(double **)); 
	buffer2d = (double  **) malloc(sizes[2] * sizes[1] * sizeof(double *)); 
	ASSIGN_3D_FORM(double, rho2d, buffer2d, rho1d, sizes);

	pres2d    = (double ***) malloc(sizes[2] * sizeof(double **));
        buffer2d = (double  **) malloc(sizes[2] * sizes[1] * sizeof(double *));
        ASSIGN_3D_FORM(double, pres2d, buffer2d, pres1d, sizes);

	ei2d    = (double ***) malloc(sizes[2] * sizeof(double **));
        buffer2d = (double  **) malloc(sizes[2] * sizes[1] * sizeof(double *));
        ASSIGN_3D_FORM(double, ei2d, buffer2d, ei1d, sizes);

        if (!nmat_for_2dcell) {  
            ifmesh_initialized_for_mat = 0; 

            nmat_for_2dcell  = (int **) malloc(ncell_ext[1] * sizeof(int *));
            matid_for_2dcell = (int **) malloc(ncell_ext[1] * sizeof(int *));
            mixcell_for_2dcell = (int **) malloc(ncell_ext[1] * sizeof(int *));
    
            nmat_for_2dcell[0]  = (int *) malloc(llsize * sizeof(int));
            matid_for_2dcell[0] = (int *) malloc(llsize * sizeof(int));
            mixcell_for_2dcell[0] = (int *) malloc(llsize * sizeof(int));
    
            for (j = 1; j < ncell_ext[1]; j++) {
                nmat_for_2dcell[j]  =  nmat_for_2dcell[j-1] + ncell_ext[0];
                matid_for_2dcell[j] = matid_for_2dcell[j-1] + ncell_ext[0];
                mixcell_for_2dcell[j] = mixcell_for_2dcell[j-1] + ncell_ext[0];
            }
        } 
        for (j = 0; j < ncell_ext[1]; j++) {
            for (i = 0; i < ncell_ext[1]; i++) {
                nmat_for_2dcell[j][i] = 1;
                matid_for_2dcell[j][i] = 0;
                mixcell_for_2dcell[j][i] = -1; // no mixed cell  
            }
        }
    }
    else if (dim == 3) { 
        vf3d     = (double ****) malloc(sizes[3] * sizeof(double ***));
        buffer3d = (double  ***) malloc(sizes[3] * sizes[2] * sizeof(double **));
        buffer2d = (double   **) malloc(sizes[3] * sizes[2] * sizes[1] * sizeof(double *));
        ASSIGN_4D_FORM(double, sizes, vf3d, buffer3d, buffer2d, vf1d);

	rho3d     = (double ****) malloc(sizes[3] * sizeof(double ***));
        buffer3d = (double  ***) malloc(sizes[3] * sizes[2] * sizeof(double **));
        buffer2d = (double   **) malloc(sizes[3] * sizes[2] * sizes[1] * sizeof(double *));
        ASSIGN_4D_FORM(double, sizes, rho3d, buffer3d, buffer2d, rho1d);

	pres3d     = (double ****) malloc(sizes[3] * sizeof(double ***));
        buffer3d = (double  ***) malloc(sizes[3] * sizes[2] * sizeof(double **));
        buffer2d = (double   **) malloc(sizes[3] * sizes[2] * sizes[1] * sizeof(double *));
        ASSIGN_4D_FORM(double, sizes, pres3d, buffer3d, buffer2d, pres1d);

	ei3d     = (double ****) malloc(sizes[3] * sizeof(double ***));
        buffer3d = (double  ***) malloc(sizes[3] * sizes[2] * sizeof(double **));
        buffer2d = (double   **) malloc(sizes[3] * sizes[2] * sizes[1] * sizeof(double *));
        ASSIGN_4D_FORM(double, sizes, ei3d, buffer3d, buffer2d, ei1d);

        if (!nmat_for_3dcell) { 
            nmat_for_3dcell    = (int ***) malloc(ncell_ext[2] * sizeof(int **));
            matid_for_3dcell   = (int ***) malloc(ncell_ext[2] * sizeof(int **));
            mixcell_for_3dcell = (int ***) malloc(ncell_ext[2] * sizeof(int **));
    
            ibuffer2d = (int **) malloc(ncell_ext[2] * ncell_ext[1] * sizeof(int *));
            ibuffer1d = (int  *) malloc(llsize * sizeof(int));
            ASSIGN_3D_FORM(int, nmat_for_3dcell, ibuffer2d, ibuffer1d, ncell_ext); 
    
            ibuffer2d = (int **) malloc(ncell_ext[2] * ncell_ext[1] * sizeof(int *));
            ibuffer1d = (int  *) malloc(llsize * sizeof(int));
            ASSIGN_3D_FORM(int, matid_for_3dcell, ibuffer2d, ibuffer1d, ncell_ext);  
           
            ibuffer2d = (int **) malloc(ncell_ext[2] * ncell_ext[1] * sizeof(int *));
            ibuffer1d = (int  *) malloc(llsize * sizeof(int));
            ASSIGN_3D_FORM(int, mixcell_for_3dcell, ibuffer2d, ibuffer1d, ncell_ext);
        }
        for (k = 0; k < ncell_ext[2]; k++) {
            for (j = 0; j < ncell_ext[1]; j++) {
                for (i = 0; i < ncell_ext[1]; i++) {
                    nmat_for_3dcell[k][j][i] = 1;
                    matid_for_3dcell[k][j][i] = 0;
                    mixcell_for_3dcell[k][j][i] = -1;
                }
            }
        }
    }  
    nmix = 0;
    lsize = 0;
    if (dim == 2) { 
        for (j = 0; j < ncell_ext[1]; j++) { 
            for (i = 0; i < ncell_ext[0]; i++) { 
                nm = 0;
                vfs = vf2d[j][i];
                for (m = 0; m < nmat; m++) { 
                    if (vfs[m] >= vfmin) { 
                        nm++;
                        matid_for_2dcell[j][i] = m;
                    }    
                }
                nmat_for_2dcell[j][i] = nm;
                if (nm > 1) {  
                    mixcell_for_2dcell[j][i] = nmix;
                    nmix++;
                    lsize += nm; 
                }
            }
        }
    }    
    else if (dim == 3) { 
        for (k = 0; k < ncell_ext[2]; k++) {  
            for (j = 0; j < ncell_ext[1]; j++) { 
                for (i = 0; i < ncell_ext[0]; i++) { 
                    nm = 0;
                    vfs = vf3d[k][j][i];
                    for (m = 0; m < nmat; m++) { 
                        if (vfs[m] >= vfmin) { 
                            nm++;
                            matid_for_3dcell[k][j][i] = m;
                        }
                    }
                    nmat_for_3dcell[k][j][i] = nm;
                    if (nm > 1)  { 
                        mixcell_for_3dcell[k][j][i] = nmix;
                        nmix++;
                        lsize += nm; 
                    }
                }
            }
        }
    } 
    ijk_for_mixcell    = (int **) malloc(nmix * sizeof(int *));
    ijk_for_mixcell[0] = (int  *) malloc(nmix * dim * sizeof(int));
    for (s = 1; s < nmix; s++) {
        ijk_for_mixcell[s] = ijk_for_mixcell[s-1] + dim;
    }
    nmat_for_mixcell = (int *) malloc(nmix * sizeof(int));
  
    vf_for_mixcell      = (double **) malloc(nmix  * sizeof(double *));
    vf_for_mixcell[0]   = (double  *) malloc(lsize * sizeof(double));
    rho_for_mixcell     = (double **) malloc(nmix  * sizeof(double *));
    rho_for_mixcell[0]  = (double  *) malloc(lsize * sizeof(double)); 
    pres_for_mixcell    = (double **) malloc(nmix  * sizeof(double *)); 
    pres_for_mixcell[0] = (double  *) malloc(lsize * sizeof(double)); 
    ei_for_mixcell      = (double **) malloc(nmix  * sizeof(double *));
    ei_for_mixcell[0]   = (double  *) malloc(lsize * sizeof(double));

    matids_for_mixcell    = (int    **) malloc(nmix * sizeof(int *));
    matids_for_mixcell[0] = (int    *) malloc(lsize * sizeof(int));

    offset = 0;
    nmix = 0; 
    if (dim == 2) {    
        for (j = 0; j < ncell_ext[1]; j++) { 
            for (i = 0; i < ncell_ext[0]; i++) { 
                if (nmat_for_2dcell[j][i] < 2) continue;
                ijk_for_mixcell[nmix][0] = i;
                ijk_for_mixcell[nmix][1] = j; 
                nm = 0;
                vfs = vf2d[j][i];
                matids_for_mixcell[nmix] = matids_for_mixcell[0] + offset; 
                vf_for_mixcell[nmix]     = vf_for_mixcell[0]     + offset; 
		rho_for_mixcell[nmix]    = rho_for_mixcell[0]    + offset; 
		pres_for_mixcell[nmix]   = pres_for_mixcell[0]   + offset; 
		ei_for_mixcell[nmix]     = ei_for_mixcell[0]     + offset; 

                for (m = 0; m < nmat; m++) { 
                    if (vfs[m] >= vfmin) {
                        matids_for_mixcell[nmix][nm] = m;          
                        vf_for_mixcell[nmix][nm]   = vf2d[j][i][m]; 
			rho_for_mixcell[nmix][nm]  = rho2d[j][i][m]; 
			pres_for_mixcell[nmix][nm] = pres2d[j][i][m];   
			ei_for_mixcell[nmix][nm]   = ei2d[j][i][m]; 

                        nm++;
                    }
                }
                nmat_for_mixcell[nmix] = nm; 
                offset += nm; 
                nmix++;
            }
        }
    }
    else if (dim == 3) {
        for (k = 0; k < ncell_ext[2]; k++) {
            for (j = 0; j < ncell_ext[1]; j++) {
                for (i = 0; i < ncell_ext[0]; i++) {
                    if (nmat_for_3dcell[k][j][i] < 2) continue; 
                    ijk_for_mixcell[nmix][0] = i;
                    ijk_for_mixcell[nmix][1] = j;
                    ijk_for_mixcell[nmix][2] = k; 

                    nm = 0;
                    vfs = vf3d[k][j][i];
                    matids_for_mixcell[nmix] = matids_for_mixcell[0] + offset;
                    vf_for_mixcell[nmix]     = vf_for_mixcell[0]     + offset;
		    rho_for_mixcell[nmix]    = rho_for_mixcell[0]    + offset;
                    pres_for_mixcell[nmix]   = pres_for_mixcell[0]   + offset;
                    ei_for_mixcell[nmix]     = ei_for_mixcell[0]     + offset;

                    for (m = 0; m < nmat; m++) {
                        if (vfs[m] >= vfmin) {
                            matids_for_mixcell[nmix][nm] = m; 
                            vf_for_mixcell[nmix][nm]   = vf3d[k][j][i][m];
			    rho_for_mixcell[nmix][nm]  = rho3d[k][j][i][m]; 
			    pres_for_mixcell[nmix][nm] = pres3d[k][j][i][m]; 
			    ei_for_mixcell[nmix][nm]   = ei3d[k][j][i][m]; 
                            nm++;
                        }
                    }
                    nmat_for_mixcell[nmix] = nm;
                    offset += nm;
                    nmix++;
                }
            }
        }
    }
    if (!ifmesh_initialized_for_mat) {  
        mesh_get_mat(dim, xl_prob, xr_prob, ncell_prob, nbdry_prob, 
                 nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell, 
                 nmat_for_3dcell, matid_for_3dcell, mixcell_for_3dcell); 
    }   
    mat_get_mix(dim, nmat,  nmix, ijk_for_mixcell, nmat_for_mixcell, matids_for_mixcell,
                vf_for_mixcell, rho_for_mixcell, pres_for_mixcell, ei_for_mixcell);     

//  test material polyhedrons

    for (i = 0; i < dim; i++) { 
        dx[i] = (xr_prob[i] - xl_prob[i])/(double)ncell_prob[i];
    } 
    get_mpoly(dim, ncell_prob, nbdry_prob, xl_prob, dx, nmix,
              nmat_for_2dcell, matid_for_2dcell, mixcell_for_2dcell,
              nmat_for_3dcell, matid_for_3dcell, mixcell_for_3dcell); 

    if (dim == 2) { 
        free(vf2d[0]);
        free(vf2d);
	free(rho2d[0]);
        free(rho2d);
	free(pres2d[0]);
        free(pres2d);
	free(ei2d[0]);
        free(ei2d);
    } 
    else if (dim == 3) { 
        free(vf3d[0][0]);
        free(vf3d[0]);
        free(vf3d);
	free(rho3d[0][0]);
        free(rho3d[0]);
        free(rho3d);
	free(pres3d[0][0]);
        free(pres3d[0]);
        free(pres3d);
	free(ei3d[0][0]);
        free(ei3d[0]);
        free(ei3d);
    } 
//  write 
    sprintf(filename, "file_%08d", ncycle); 
    create_file(filename, t, 0, &fileid);
    write_mesh_mat(fileid, "mesh", dim, xl_prob, xr_prob, ncell_prob, nbdry_prob,
		    nmat, matids, 
		    nmat_for_2dcell, matid_for_2dcell,
                    nmat_for_3dcell, matid_for_3dcell, &meshid); 

    write_mat(fileid, dim); 
    
    close_file(fileid);

    free(matids);

    return;
} 

void advected_vol(int dim, int direction, 
                  int edge_index, int *cell_index, double adv_distance, 
                  int nmat, double *vols_advected)
{ 
     int  m, i, j, k, nm, matid, mix, nnode; 
     int  *matids, *nnode_for_face, **nodes_for_face;
     int  *nface_for_zone, **nnode_for_face_ea_mpoly, **nodelist_for_face_ea_mpoly;
   
     double small, delta, dvol, dx[3], xl[3], ijk[3];
     double *coords; 

     int    *ncell_prob, nbdry_prob; 
     int    **nmat_for_2dcell, **matid_for_2dcell, **mixcell_for_2dcell; 
     int    ***nmat_for_3dcell, ***matid_for_3dcell, ***mixcell_for_3dcell; 
     double *xl_prob, *xr_prob; 

     int   nmixcell;
     int   **ijk_for_mixcell, *nmat_for_mixcell, **matids_for_mixcell;
     int   *nnode_for_mixcell;
     int   **nnode_for_mpoly_for_mixcell, ***nodes_for_mpoly_for_mixcell;
     int   **nface_for_mpoly_for_mixcell, ***nnode_for_face_ea_mpoly_for_mixcell;
     int   ***nodelist_for_face_ea_mpoly_for_mixcell;  
     double **coords_for_mixcell;

     small = 1.0e-06; 

     mesh_pass_mat(&xl_prob, &xr_prob, &ncell_prob, &nbdry_prob, 
              &nmat_for_2dcell, &matid_for_2dcell, &mixcell_for_2dcell,
              &nmat_for_3dcell, &matid_for_3dcell, &mixcell_for_3dcell);

     for (m = 0; m < nmat; m++) {
         vols_advected[m] = 0.0;
     }  
     for (i = 0; i < dim; i++) { 
         dx[i] = (xr_prob[i] - xl_prob[i])/(double)ncell_prob[i];
     }
     dvol = 1.0;
     for (i = 0; i < dim; i++) {
         if (i == direction - 1) continue;
         dvol *= dx[i];
     }
     delta = fabs(adv_distance)/dx[direction-1]; 
     if (delta < small) return;
           
     mat_pass_mix(&nmixcell, &ijk_for_mixcell, &nmat_for_mixcell, 
                  &matids_for_mixcell, &nnode_for_mixcell, &coords_for_mixcell,
                  &nnode_for_mpoly_for_mixcell, &nodes_for_mpoly_for_mixcell,
                  &nface_for_mpoly_for_mixcell, &nnode_for_face_ea_mpoly_for_mixcell,
                  &nodelist_for_face_ea_mpoly_for_mixcell);  

//   detremine which cell 

     for (i = 0; i < dim; i++) { 
         ijk[i] = nbdry_prob + cell_index[i] - 1;  // from 1-based to 0-based 
     }  
     if (adv_distance > 0.0) { 
         ijk[direction-1] = nbdry_prob + edge_index - 2;   // cell_index is 1-based 
     } 
     else { 
         ijk[direction-1] = nbdry_prob + edge_index - 1; 
     } 
     for (i = 0; i < dim; i++) { 
         xl[i] = xl_prob[i] + (ijk[i]  - nbdry_prob) * dx[i];
     } 
     if (dim == 2) { 
         i     = ijk[0];
         j     = ijk[1];
         nm    =    nmat_for_2dcell[j][i];   
         matid =   matid_for_2dcell[j][i];  
         mix   = mixcell_for_2dcell[j][i];
     }
     else if (dim == 3) { 
         i     = ijk[0];
         j     = ijk[1]; 
         k     = ijk[2]; 
         nm    =    nmat_for_3dcell[k][j][i];
         matid =   matid_for_3dcell[k][j][i];
         mix   = mixcell_for_3dcell[k][j][i];
     } 
     if (nm == 1) { 
         vols_advected[matid] = fabs(adv_distance) * dvol;
         return;
     }
     nnode  = nnode_for_mixcell[mix];
     coords = coords_for_mixcell[mix];
     matids = matids_for_mixcell[mix];
 
     if (dim == 2) { 
         nnode_for_face = nnode_for_mpoly_for_mixcell[mix];
         nodes_for_face = nodes_for_mpoly_for_mixcell[mix];
     }
     else if (dim == 3) { 
         nface_for_zone             = nface_for_mpoly_for_mixcell[mix];
         nnode_for_face_ea_mpoly    = nnode_for_face_ea_mpoly_for_mixcell[mix];
         nodelist_for_face_ea_mpoly = nodelist_for_face_ea_mpoly_for_mixcell[mix]; 
     }
     if (dim == 2) { 
         advected_vol2d(xl, dx, direction-1, adv_distance, nnode, coords, 
                        nm, matids, nnode_for_face, nodes_for_face, nmat, vols_advected);  
     }
     else if (dim == 3) { 

         advected_vol3d(xl, dx, direction-1, adv_distance, nnode, coords, nm, matids, 
                        nface_for_zone, nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly,
                        nmat, vols_advected);

     } 

     return;
 } 
                
void advected_vol2d(double *xl, double *dx, int direction, double advected_distance,  
                    int nnode_mat, double *coords_mat, int nmpoly, int *matids,  
                    int *nnode_for_face_mat, int **nodes_for_face_mat,
                    int nmat, double *vols_mat)  
{
//   advected_distance = u_averaged * dt 

     int dim, i, k, n, m, dir_other, matid, nnode_new, nnode_upper, nnode_lower;
     int nnode_face, *nodes_face; 
     int *nodelist_default, *node_loc, *nodelist_upper, *nodelist_lower;
     int nnode_interface, nodes_interface[2]; 
     double small, advected_ds, distance, vol_previous, vol, vol_tot, vol_sum; 
     double myxl[2], myxr[2], mydx[2], norm_inward[2], coords_plane[4], *c0, *c1; 
     double *ds_ea_node, *coords_scaled, *coords_mpoly, *coords_new, *coords_sol;

     small = 1.0e-06;
     dim = 2; 

     node_loc = (int *) malloc((nnode_mat + nnode_mat) * sizeof(int));
     nodelist_default = node_loc + nnode_mat; 
     for (i = 0; i < nnode_mat; i++) { 
          nodelist_default[i] = i;
     }
//   scale 

     coords_scaled = (double *) malloc((nnode_mat+4*nnode_mat*dim + 1) * sizeof(double));
     coords_mpoly  = coords_scaled + (nnode_mat * dim);
     coords_new    = coords_mpoly  + (nnode_mat * dim);
     coords_sol    = coords_new    + (nnode_mat * dim);
     ds_ea_node    = coords_sol    + (nnode_mat * dim) + 1;

     for (n = 0; n < nnode_mat; n++) {
         c0 = coords_mat + (n * dim);
         c1 = coords_scaled + (n * dim);
         for (i = 0; i < dim; i++) { 
             c1[i] = (c0[i] - xl[i]) / dx[i];
         }
     }
     advected_ds = advected_distance/dx[direction];

     for (i = 0; i < dim; i++) { 
         myxl[i] = 0.0;
         myxr[i] = 1.0;
         mydx[i] = 1.0;
     } 
     dir_other = 1 - direction;
     for (i = 0; i < dim; i++) { 
         norm_inward[i] = 0.0;
     } 
     c0 = coords_plane;
     c1 = coords_plane + dim;
    
     c0[dir_other] = myxl[dir_other] - mydx[dir_other];
     c1[dir_other] = myxr[dir_other] + mydx[dir_other];

     if (advected_distance > 0.0) {
         c0[direction] = myxr[direction] - advected_distance;
         norm_inward[direction] = 1.0;
     } 
     else { 
         c0[direction] = myxl[direction] - advected_distance;  // + fabs(advected_distance) 
         norm_inward[direction] = -1.0;
     }
     c1[direction] = c0[direction];  

     distance = 0.5 *(norm_inward[0] *(c0[0] + c1[0]) + norm_inward[1] * (c0[1] + c1[1]));
     
     for (m = 0; m < nmpoly-1; m++) { 
         matid = matids[m];
         nnode_face = nnode_for_face_mat[m];
         nodes_face = nodes_for_face_mat[m];
         for (k = 0; k < nnode_face; k++) { 
             ds_ea_node[k] = 0.0;
             n  = nodes_face[k];
             c0 = coords_scaled + (n * dim); 
             for (i = 0; i < dim; i++) { 
                 ds_ea_node[k] += (norm_inward[i] * c0[i]);
             }
         } 
         for (k = 0; k < nnode_face; k++) {
            if (fabs(ds_ea_node[k] - distance) <= small) {
                node_loc[k] = 0;  // on the plane
            }
            else if (ds_ea_node[k] > distance)  {   // above the plane
                node_loc[k] = 1;
            }
            else {
                node_loc[k] = -1;
            }
        }
        nodelist_lower = NULL;
        nodelist_upper = NULL;

        for (k = 0; k < nnode_face; k++) { 
            n = nodes_face[k];
            c0 = coords_scaled + (n * dim);
            c1 = coords_mpoly + (k * dim);
            memcpy(c1, c0, (size_t)(dim * sizeof(double)));
        }  
        find_interface2d(nnode_face, coords_mpoly, nodelist_default, norm_inward,
                         distance, node_loc,
                         &nnode_new, coords_new,
                         &nnode_interface, nodes_interface,
                         &nnode_lower, &nodelist_lower,
                         &nnode_upper, &nodelist_upper);

        if (nnode_upper >= 0) { 
            for (k = 0; k < nnode_upper; k++) { 
                n = nodelist_upper[k];
                if (n < nnode_face) {
                    c0 = coords_mpoly + (n * dim);
                }
                else {
                    n -= nnode_face;
                    c0 = coords_new + (n * dim);
                }
                c1 = coords_sol + (k * dim);
                memcpy(c1, c0, (size_t)(dim * sizeof(double))); 
            } 
//          scale back 
            for (n = 0; n < nnode_new; n++) { 
                c0 = coords_sol + (n * dim);
                for (i = 0; i < dim; i++) { 
                    c0[i] = (c0[i] * dx[i] + xl[i]); 
                }
            }
            cal_poly_area(nnode_upper, coords_sol, nnode_upper, NULL, &vol);
            vols_mat[matid] = vol - vol_previous;  
            vol_previous    = vol; 
        } 
        if (nodelist_lower) free(nodelist_lower);
        if (nodelist_upper) free(nodelist_upper); 
     }  
     vol_tot = 1.0;
     for (i = 0; i < dim; i++) { 
         if (i == direction) continue;
         vol_tot *= dx[i];
     }
     vol_tot *= fabs(advected_distance);
     vol_sum = 0.0;
     for (m = 0; m < nmpoly-1; m++) { 
         vol_sum += vols_mat[m];
     } 
     matid = matids[nmpoly-1];
     vols_mat[matid] = vol_tot - vol_sum;
     
     free(node_loc);
     free(coords_scaled); 
     
     return;
  } 


void advected_vol3d(double *xl, double *dx, int direction, double advected_distance,  
                    int nnode_mat, double *coords_mat, int nmpoly, int *matids,  
                    int *nface_for_mpoly_mat, int **nnode_for_face_mat,
                    int **nodelist_for_face_mat, 
                    int nmat, double *vols_mat)  
{
//   advected_distance = u_averaged * dt 

     int dim, i, k, n, nn_ea_face, m, matid;
     int nface, nface_lower, nnode_lower, nface_upper, nnode_upper;
     int *nnode_for_face, *nodelist_for_face;
     int *nnode_for_face_lower, *nodelist_for_face_lower;
     int *nnode_for_face_upper, *nodelist_for_face_upper;  
     double small, advected_ds, distance, vol_previous, vol, vol_tot, vol_sum; 
     double myxl[3], myxr[3], mydx[3], norm_inward[3], ctr[3], c4[4][3], *c0, *c1; 
     double *coords_lower, *coords_upper, *coords_scaled, *coords_sol;

     dim = 3; 
     nn_ea_face = 4; 

//   scale 

     coords_scaled = (double *) malloc(2*nnode_mat*dim * sizeof(double));
     coords_sol    = coords_scaled + (nnode_mat * dim);

     for (n = 0; n < nnode_mat; n++) {
         c0 = coords_mat + (n * dim);
         c1 = coords_scaled + (n * dim);
         for (i = 0; i < dim; i++) { 
             c1[i] = (c0[i] - xl[i]) / dx[i];
         }
     }
     advected_ds = advected_distance/dx[direction];

     for (i = 0; i < dim; i++) { 
         myxl[i] = 0.0;
         myxr[i] = 1.0;
         mydx[i] = 1.0;
     } 
     for (i = 0; i < dim; i++) { 
         norm_inward[i] = 0.0;
     } 
     if (direction == 0) { 
      
//     2 ^   ----------  
//       |   !        ! 
//       |   !        ! 
//           ----------
//           .....>  1 

         c4[0][1] = myxl[1] - mydx[1]; 
         c4[0][2] = myxl[2] - mydx[2]; 
         c4[1][1] = myxr[1] + mydx[1];
         c4[1][2] = c4[0][2]; 
         c4[2][1] = c4[1][1];
         c4[2][2] = myxr[2] + mydx[2];
         c4[3][1] = c4[0][1];
         c4[3][2] = c4[2][2];
     }  
     else if (direction == 1) { 
//       
//     0 ^   ----------  
//       |   !        ! 
//       |   !        ! 
//           ----------
//           .....>  0 
// 
         c4[0][2] = myxl[2] - mydx[2]; 
         c4[0][0] = myxl[0] - mydx[0]; 
         c4[1][2] = myxr[2] + mydx[2];
         c4[1][0] = c4[0][0];
         c4[2][2] = c4[1][2];
         c4[2][0] = myxr[0] + mydx[0];
         c4[3][2] = c4[0][2];
         c4[3][0] = c4[2][0];
     }
     else if (direction == 2) { 
//       
//     2 ^   ----------  
//       |   !        ! 
//       |   !        ! 
//           ----------
//           .....>  1 
// 
         c4[0][0] = myxl[0] - mydx[0];
         c4[0][1] = myxl[1] - mydx[1];
         c4[1][0] = myxr[0] + mydx[0];
         c4[1][1] = c4[0][1];
         c4[2][0] = c4[1][0]; 
         c4[2][1] = myxr[1] + mydx[1];
         c4[3][0] = c4[0][0];
         c4[3][1] = c4[2][1];
     } 
     if (advected_distance > 0.0) {
         for (n = 0; n < nn_ea_face; n++) { 
             c4[n][direction] = myxr[direction] - advected_ds; 
         }
         norm_inward[direction] = 1.0;
     } 
     else { 
         for (n = 0; n < nn_ea_face; n++) { 
             c4[n][direction] = myxr[direction] + advected_ds;
         }
         norm_inward[direction] = -1.0;
     }
     distance = 0.0;
     for (n = 0; n < nn_ea_face; n++) {
         for (i = 0; i < dim; i++) { 
             distance += (c4[n][i] * norm_inward[i]);
         }
     }  
     distance /= ((double)nn_ea_face); 

     vol_previous = 0.0; 
     for (m = 0; m < nmpoly-1; m++) { 
         matid = matids[m];

         nface = nface_for_mpoly_mat[m];
         nnode_for_face = nnode_for_face_mat[m];
         nodelist_for_face = nodelist_for_face_mat[m];
                        
         nface_lower = 0;
         nnode_lower = 0;
         coords_lower = NULL;
         nnode_for_face_lower    = NULL;
         nodelist_for_face_lower = NULL;

         nface_upper = 0;
         nnode_upper = 0;
         coords_upper = NULL;
         nnode_for_face_upper    = NULL;
         nodelist_for_face_upper = NULL; 

         polyhedron_plane(nface, nnode_mat, coords_scaled,
                          nnode_for_face, nodelist_for_face,
                          norm_inward, distance,
                          &nface_lower, &nnode_lower, &coords_lower,
                          &nnode_for_face_lower, &nodelist_for_face_lower,
                          &nface_upper, &nnode_upper, &coords_upper,
                          &nnode_for_face_upper, &nodelist_for_face_upper);

        if (nface_upper > 3) { 
//          scale back 
            for (n = 0; n < nnode_upper; n++) { 
                c0 = coords_upper + (n * dim);
                c1 = coords_sol   + (n * dim);
                for (i = 0; i < dim; i++) { 
                    c1[i] = (c0[i] * dx[i] + xl[i]); 
                }
            }
            cal_vol(nface_upper, nnode_upper, coords_upper, 
                    nnode_for_face_upper, nodelist_for_face_upper,
                    &vol); 

            vols_mat[matid] = vol - vol_previous;  
            vol_previous    = vol; 
        } 
        if (coords_lower)            free(coords_lower);
        if (nnode_for_face_lower)    free(nnode_for_face_lower);
        if (nodelist_for_face_lower) free(nodelist_for_face_lower);

        if (coords_upper)            free(coords_upper);
        if (nnode_for_face_upper)    free(nnode_for_face_upper);
        if (nodelist_for_face_upper) free(nodelist_for_face_upper);
     }  
     vol_tot = 1.0;
     for (i = 0; i < dim; i++) { 
         if (i == direction) continue;
         vol_tot *= dx[i];
     }
     vol_tot *= fabs(advected_distance);
     vol_sum = 0.0;
     for (m = 0; m < nmpoly-1; m++) { 
         vol_sum += vols_mat[m];
     } 
     matid = matids[nmpoly-1];
     vols_mat[matid] = vol_tot - vol_sum;
     
     free(coords_scaled); 
     
     return;

  } 


