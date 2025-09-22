///////////////////////////////////////////////////////////////////////
//   This file contains fuctions to calculate the distance needed 
//   in the definition of the reconstructed material interface. 
///////////////////////////////////////////////////////////////////////
#include "globals.h"
#include "mesh.h" 
#include "vof2d.h"
#include "util.h"
#include "io.h" 

static double accuracy = 1.0e-10;
static int    niter_mx = 50; 

//////////////////////////////////////////////////////////////////////////////////////

void reconstruct2d_nmat(int geop,
                        double *xl, double *dx, int nmat, double *vf_ea_mat,
                        double *norm_ea_mat, double *priority_ea_mat,
                        int *nnode_final, double **coords_final,
                        int *nnode_for_interface, int **nodes_for_interface,
                        int *nnode_for_mpoly, int ***nodelist_for_mpoly);

void cal_distance2d(int geop,
                    double vf_to_match, double volume, double *norm,
                    int nnode, double *coords, int *nodelist,
                    int *nnode_new, double *coords_new,
                    double *distance,
                    int *nnode_interface, int *nodes_interface,
                    int *nnode_lower, int **nodelist_lower,
                    int *nnode_upper, int **nodelist_upper);

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

void remap2d_mat(int ifdump,
                 int nnode, double *coords, double *inward_norm_ea_face,
                 int nnode_m, double *coords_m,
                 int *nnodex, double **coordsx);

//////////////////////////////////////////////////////////////////////////////////////
void reconstruct2d_nmat_pagosa(int geop, 
                        double *xl, double *dx, 
                        int **nmat_mesh, int ***matid_mesh, double ***vf_mesh,  
                        int *nnode_final, double *coords_final, 
                        int *nnode_for_interface, int **nodes_for_interface, 
                        int *nnode_for_mpoly, int **nodelist_for_mpoly)
{ 
//    nmat_mesh[3][3], matid_mesh[3][3][nm], vf_mesh[3][3][nm]

     int nodelist_default[] = {0,1,2,3}; 
     int dim, i, j, nnode, nnode_new, nnode_previous, nnode_lower, nnode_upper;
     int nmat, nm_neighb, m, m_neighb, m_other, idx;
     int nnode_interface; 
     int n, nnew;   
     int *ids, *ids_neighb; 
     int *nodelist_lower, *nodelist_upper;
     int nodelist_interface[2], mydx[2];
     double dx_mx, cell_vol, vf, distance; 
     double *c, *coords,  coords_new[4]; 
     double priority_max, priority, vsq, factor;  
     double vfs[2][3][3], grad[2], norms[2][2], norm[2], *mynorm; 
     
     dim = 2;

     dx_mx = dx[0];
     if (dx_mx < dx[1]) dx_mx = dx[1];
 
     cell_vol = 1.0;
     for (i = 0; i < dim; i++) { 
         mydx[i] = dx[i]/dx_mx;
         cell_vol *= mydx[i];
     } 
     nnode  = 4;
     coords = (double *) malloc(nnode * dim * sizeof(double));
     c = coords;
     c[0] = 0.0; 
     c[1] = 0.0;
     c += dim;
     c[0] = mydx[0];
     c[1] = 0.0;
     c += dim;
     c[0] = mydx[0];
     c[1] = mydx[1];
     c += dim;
     c[0] = 0.0;
     c[1] = mydx[1];

     nmat = nmat_mesh[1][1];
     assert(nmat > 1); 

     memcpy(coords_final, coords, (size_t)(nnode * dim * sizeof(double)));
     *nnode_final = nnode;

     for (j = 0; j <= 2; j++) { 
         for (i = 0; i <= 2; i++) { 
             vfs[0][j][i] = 0.0;
         }
     }  
     for (m = 0; m < nmat-1; m++) { 
         
         for (j = 0; j <= 2; j++) {
             for (i = 0; i <= 2; i++) { 
                 if ((i == 1) && (j == 1)) continue;
                 vfs[1][j][i] = 0.0;  
             }
         } 
         vfs[0][1][1] += vf_mesh[1][1][m];   // vfs[inner:outer][i][j]
         vfs[1][1][1]  = 1.0 - vfs[0][1][1]; 
          
         ids = matid_mesh[1][1];
         for (j = 0; j <= 2; j++) { 
             for (i = 0; i <= 2; i++) {  
                 if ((i == 1) && (j == 1)) continue; 

                 nm_neighb  = nmat_mesh[j][i]; 
                 ids_neighb = matid_mesh[j][i];  
                 for (m_neighb = 0; m_neighb < nm_neighb; m_neighb++) { 
                     if (ids_neighb[m_neighb] == ids[m]) { 
                         vfs[0][j][i] += vf_mesh[j][i][m_neighb];     
                     } 
                     else { 
                         for (m_other = m + 1; m_other < nmat; m_other++) { 
                             if (ids_neighb[m_neighb] == ids[m_other]) { 
                                 vfs[1][j][i] += vf_mesh[j][i][m_neighb]; 
                                 break;
                             }
                         }
                     } 
                 }
             }
         } 
//       now, we have the distribution of vfs[0] and vfs[1] each of which is a 2x2 mesh 
         
         priority_max = 0.0;
         idx = 0;
         for (m_other = 0; m_other < 2; m_other++) { 
             mynorm = norms[m_other]; 
             cal_cell_zgrad2d(dx, vfs[m_other], grad);
             vsq = 0.0;
             for (i = 0; i < dim; i++) {
                 mynorm[i] = - grad[i];
                 vsq += (mynorm[i] * mynorm[i]);
             }
             factor = sqrt(vsq); 
             priority = factor * sqrt(vfs[m_other][1][1]);
             if (priority > priority_max) {  
                 idx = m_other;
             }
             factor = 1.0/(factor + tiny); 
             for (i = 0; i < dim; i++) {
                 mynorm[i] *= factor;
             }
         } 
         mynorm = norms[idx];
         if (idx == 1) {
             for (i = 0; i < dim; i++) { 
                 norm[i] = - mynorm[i];
             }
         }
         else { 
             for (i = 0; i < dim; i++) { 
                 norm[i] = mynorm[i];
             }
         } 
         nnode_new  = 0;
         nnode_interface = 0;
         nnode_lower = 0;
         nnode_upper = 0;
         nodelist_lower = NULL;
         nodelist_upper = NULL;

         cal_distance2d(geop, vfs[0][1][1], cell_vol, norm, nnode, coords, nodelist_default,
                        &nnode_new, coords_new, &distance, 
                        &nnode_interface, nodelist_interface,
                        &nnode_lower,     &nodelist_lower,
                        &nnode_upper,     &nodelist_upper); 

         memcpy(coords_final + (*nnode_final * dim), coords_new, (size_t)(nnode_new * dim * sizeof(double)));

         nnode_for_interface[m] = nnode_interface;
         if (nnode_interface == 2) { 
             for (i = 0; i < 2; i++) { 
                 n = nodelist_interface[i];
                 if (n < nnode) { 
                     nodes_for_interface[m][i] = n; 
                 }
                 else { 
                     nnew = n - nnode; 
                     nodes_for_interface[m][i] = *nnode_final + nnew;
                 } 
             }
         }  
         nnode_for_mpoly[m] = nnode_lower;
         
         for (i = 0; i < nnode_lower; i++) { 
             n = nodelist_lower[i];
             if (n < nnode) {                    
                 nodelist_for_mpoly[m][i] = n;
             }
             else { 
                 nnew = n - nnode;  
                 nodelist_for_mpoly[m][i] = *nnode_final + nnew;
             } 
         }
         if (nodelist_lower) { 
             free(nodelist_lower);
             nodelist_lower = NULL;
         }
         nodelist_for_mpoly[m+1] = nodelist_for_mpoly[m] + nnode_lower;

         if (m < nmat - 2) { 
             if (nodelist_upper) free(nodelist_upper);
         }
         nnode_previous = *nnode_final;    // for the last nodelist_upper 
         *nnode_final  += nnode_new;
     }   // m 
     nnode_for_mpoly[m] = nnode_upper; 
     for (i = 0; i < nnode_upper; i++) { 
         n = nodelist_upper[i];
         if (n < nnode) { 
             nodelist_for_mpoly[m][i] = n;
         }
         else { 
             nnew = n - nnode;
             nodelist_for_mpoly[m][i] = nnode_previous + nnew;
         } 
     }  
     free(coords);

//   scale coords back 

     for (n = 0; n < *nnode_final; n++) {
         c = coords_final + n * dim;
         for (i = 0; i < dim; i++) { 
             c[i] = xl[i] + dx_mx * c[i]; 
         }
     }
     return; 
 } 
 

void reconstruct2d_nmat(int geop, 
                        double *xl, double *dx, int nmat, double *vf_ea_mat, 
                        double *norm_ea_mat, double *priority_ea_mat,  
                        int *nnode_final, double **coords_final, 
                        int *nnode_for_interface, int **nodes_for_interface, 
                        int *nnode_for_mpoly, int ***nodelist_for_mpoly)
{ 
     int dim, i, m, nnode, nnode_new, nnode_tot, nnode_lower, nnode_upper;
     int nnode_interface; 
     int *nodelist_lower, *nodelist_upper;
     int nodelist_interface[2], mydx[2];
     double dx_mx, vol, vf, distance; 
     double *c, *norm, *coords,  coords_new[4]; 
    
     int lsize_mx, lsize, nnode_current, n, nnew, nold; 
     int *node_new2old, *node_new2old_tmp, *nodelist_current, *ip; 
     
     double vf_current, vol_current; 
     double *coords_current, *coords_tmp; 

     dim = 2;

     dx_mx = dx[0];
     if (dx_mx < dx[1]) dx_mx = dx[1];
 
     vol = 1.0;
     for (i = 0; i < dim; i++) { 
         mydx[i] = dx[i]/dx_mx;
         vol *= mydx[i];
     } 
     nnode  = 4;
     coords = (double *) malloc(nnode * dim * sizeof(double));
     c = coords;
     c[0] = 0.0; 
     c[1] = 0.0;
     c += dim;
     c[0] = mydx[0];
     c[1] = 0.0;
     c += dim;
     c[0] = mydx[0];
     c[1] = mydx[1];
     c += dim;
     c[0] = 0.0;
     c[1] = mydx[1];

     *coords_final = (double *) malloc((nnode + nmat + nmat) * dim * sizeof(double));
     memcpy(*coords_final, coords, (size_t)(nnode * dim * sizeof(double)));
     *nnode_final = nnode;

     lsize_mx = nmat * nnode;
     *nodelist_for_mpoly      = (int **) malloc(nmat * sizeof(double *));
     (*nodelist_for_mpoly)[0] = (int  *) malloc(lsize_mx * sizeof(double));

     node_new2old = (int *) malloc(nnode * nmat * sizeof(int));
     node_new2old_tmp    = (int *) malloc(nnode * nmat * sizeof(int));
     for (i = 0; i < nnode; i++) { 
         node_new2old[i] = i;
     } 
     nnode_current = nnode;
     coords_current = (double *) malloc(nnode *nmat * dim * sizeof(double));
     coords_tmp     = (double *) malloc(nnode *nmat * dim * sizeof(double));

     memcpy(coords_current, coords, (size_t)(nnode_current * dim * sizeof(double)));
     nodelist_current = (int *) malloc(nnode *nmat * sizeof(int));
     
     for (i = 0; i < nnode; i++) { 
         nodelist_current[i] = i;
     }
     vf_current  = vf_ea_mat[0];
     vol_current = vol;
     
     lsize = 0;
     for (m = 0; m < nmat-1; m++) { 
          
         norm = norm_ea_mat + dim * m;

         nnode_new  = 0;
         nnode_interface = 0;
         nnode_lower = 0;
         nnode_upper = 0;
         nodelist_lower = NULL;
         nodelist_upper = NULL;

         cal_distance2d(geop, vf_current, vol, norm, nnode_current, coords_current, nodelist_current,
                        &nnode_new, coords_new, &distance, 
                        &nnode_interface, nodelist_interface,
                        &nnode_lower,     &nodelist_lower,
                        &nnode_upper,     &nodelist_upper); 

         nnode_for_interface[m] = nnode_interface;
         if (nnode_interface == 2) { 
             for (i = 0; i < 2; i++) { 
                 n = nodelist_interface[i];
                 if (n < nnode_current) { 
                     nold = node_new2old[n];       
                     nodes_for_interface[m][i] = nold; 
                 }
                 else { 
                     nnew = n - nnode_current; 
                     nodes_for_interface[m][i] = *nnode_final + nnew;
                 } 
             }
         }  
         nnode_for_mpoly[m] = nnode_lower;
         
         for (i = 0; i < nnode_lower; i++) { 
             n = nodelist_lower[i];
             if (n < nnode_current) {                    
                 nold = node_new2old[n]; 
                 (*nodelist_for_mpoly)[m][i] = nold;
             }
             else { 
                 nnew = n - nnode_current;  
                 (*nodelist_for_mpoly)[m][i] = *nnode_final + nnew;
             } 
         }
         if (nodelist_lower) { 
             free(nodelist_lower);
             nodelist_lower = NULL;
         }
         (*nodelist_for_mpoly)[m+1] = (*nodelist_for_mpoly)[m] + nnode_lower; 

         for (i = 0; i < nnode_upper; i++) { 
             nodelist_current[i] = i; 
             n = nodelist_upper[i];
             if (n < nnode_current) { 
                 memcpy(coords_tmp + i*dim, coords_current+n*dim, (size_t)(dim*sizeof(double)));
                 node_new2old_tmp[i] = node_new2old[n]; 
             }
             else {  
                 nnew = n - nnode_current;  
                 memcpy(coords_tmp + i*dim, coords_new+nnew*dim, (size_t)(dim*sizeof(double)));
                 node_new2old_tmp[i] = *nnode_final + nnew;
             } 
         }
         free(nodelist_upper);
         nodelist_upper = NULL;

         c              = coords_tmp;
         coords_tmp     = coords_current;
         coords_current = c;
   
         ip           = node_new2old;
         node_new2old = node_new2old_tmp;
         node_new2old_tmp = ip;

         memcpy(*coords_final+(*nnode_final)*dim,coords_new,(size_t)(nnode_new*dim*sizeof(double)));
         *nnode_final += nnode_new;

         vol_current -= (vf_ea_mat[m] * vol); 
         vf = (vf_ea_mat[m+1] * vol)/vol_current; 
 
         lsize += nnode_lower; 
     }   // m 
     if (lsize + nnode_upper > lsize_mx) { 
         lsize_mx = lsize + nnode_upper; 
         (*nodelist_for_mpoly)[0] = (int *)realloc((*nodelist_for_mpoly)[0],(size_t)(lsize_mx*sizeof(int)));
         for (i = 1; i < nmat; i++) { 
             (*nodelist_for_mpoly)[i] = (*nodelist_for_mpoly)[i-1] + nnode_for_mpoly[i-1]; 
         } 
     }
     nnode_for_mpoly[m] = nnode_upper; 
     for (i = 0; i < nnode_upper; i++) { 
         (*nodelist_for_mpoly)[m][i] = node_new2old[i]; 
     }  

     free(nodelist_current); 
     free(coords_tmp);
     free(coords_current);
     free(node_new2old);
     free(node_new2old_tmp); 
  
     free(coords);

//   scale coords back 

     for (n = 0; n < *nnode_final; n++) {
         c = *coords_final + n * dim;
         for (i = 0; i < dim; i++) { 
             c[i] = xl[i] + dx_mx * c[i]; 
         }
     }
 

     return; 
 } 
 

void cal_distance2d(int geop,
                    double vf_to_match, double volume, double *norm, 
                    int nnode, double *coords, int *nodelist, 
                    int *nnode_new, double *coords_new, 
                    double *distance, 
                    int *nnode_interface, int *nodes_interface,
                    int *nnode_lower, int **nodelist_lower,
                    int *nnode_upper, int **nodelist_upper)
{ 
//   norm[0] * x + norm[1] * y = distance 

     int dim, szdim, i, n, niter, niter_mx, if_vol_matched, done; 
     int *node_order_for_ds, *node_loc; 
     double err, ds0, ds1, vf0, vf1, vf, vol; 
     double ds_lower, ds_upper, vf_lower, vf_upper; 
     double *c0, *c1, *ds_ea_node, *coords_work; 

     niter_mx = 100;

     dim = 2;
     szdim = dim * sizeof(double); 

     *nnode_interface = 0;
     *nnode_lower     = 0;
     *nnode_upper     = 0;
     *nodelist_lower  = NULL;
     *nodelist_upper  = NULL;

     coords_work = (double *) malloc((nnode + nnode) * szdim);

     node_loc          = (int    *) malloc(nnode * sizeof(int));
     node_order_for_ds = (int    *) malloc(nnode * sizeof(int));
     ds_ea_node        = (double *) malloc(nnode * sizeof(double));
     order_nodes_along_norm(dim, norm, nnode, coords, node_order_for_ds, ds_ea_node);

     
     bounds_2d(geop, vf_to_match, volume, norm, 
               nnode, coords, nodelist,
               node_order_for_ds, ds_ea_node, 
               &ds_lower, &ds_upper, &vf_lower, &vf_upper, 
               nnode_new, coords_new, 
               &if_vol_matched, nnode_interface, nodes_interface, 
               nnode_lower, nodelist_lower, nnode_upper, nodelist_upper);

     if (if_vol_matched) { 
         *distance = ds_upper;
     }
     else { 
    
         err = 1.0;
         ds1 = ds_upper;
         ds0 = ds_lower;
         vf0 = vf_lower;
         vf1 = vf_upper; 
         vf  = vf0; 

         *distance = 0.5 *(ds0 + ds1); 
         niter = 0;
         done  = 0;
         while (!done && (niter < niter_mx)) { 

               if (*nodelist_lower) { 
                   free(*nodelist_lower);
                   *nodelist_lower = NULL;
               }  
               if (*nodelist_upper) { 
                   free(*nodelist_upper);  
                   *nodelist_upper = NULL;
               }
               for (i = 0; i < nnode; i++) { 
                   if (fabs(ds_ea_node[i] - *distance) <= accuracy) { 
                       node_loc[i] = 0; 
                   }
                   else if (ds_ea_node[i] > *distance)  { 
                       node_loc[i] = 1;
                   }
                   else { 
                       node_loc[i] = -1;
                   }
               } 
               find_interface2d(nnode, coords, nodelist, norm,
                                *distance, node_loc, 
                                nnode_new, coords_new, 
                                nnode_interface, nodes_interface, 
                                nnode_lower, nodelist_lower, 
                                nnode_upper, nodelist_upper);

               if (*nnode_interface < 2) {  
                   *nnode_interface = 0;
                   *nnode_new = 0; 
     
                   if (*nodelist_lower) free(*nodelist_lower);
                   *nodelist_lower = NULL;
                   if (*nodelist_upper) free(*nodelist_upper); 
                   *nodelist_upper = NULL;

                   if (vf_to_match > 0.5) { 
                       *nnode_lower = 0; 
                       
                       *nnode_upper = nnode;
                       *nodelist_upper = (int *) malloc(nnode * sizeof(int));
                       memcpy(*nodelist_upper, nodelist, (size_t)(nnode * sizeof(int)));

                       n = node_order_for_ds[nnode-1];
                       *distance = ds_ea_node[n];
                   }
                   else { 
                       *nnode_upper = 0; 

                       *nnode_lower = nnode;
                       *nodelist_lower = (int *) malloc(nnode * sizeof(int));
                       memcpy(*nodelist_lower, nodelist, (size_t)(nnode * sizeof(int)));

                       n = node_order_for_ds[0];
                       *distance = ds_ea_node[n];
                   }
                   done = 1;
               }
               else { 
                   for (i = 0; i < *nnode_lower; i++) { 
                       n = (*nodelist_lower)[i]; 
                       if (n < nnode) { 
                           c0 = coords + (n + n);
                       }
                       else { 
                           n -= nnode;
                           c0 = coords_new + (n + n);
                       }
                       c1 = coords_work + (i + i);
                       memcpy(c1, c0, (size_t)szdim);
                   }
                   if (geop == 1) {
                       cal_poly_area(*nnode_lower, coords_work, *nnode_lower, NULL, &vol);
                   }
                   else if (geop == 2) {
//                     comment out for TT effort 
//                     rz_area(*nnode_lower, coords_work, &vol);
                       printf("ERROR: activate the line above\n");
                   }
                   vf = vol/volume;
                   err = fabs(vf - vf_to_match);

                  if ((err <= accuracy) || (ds1 - ds0 <= accuracy)) {
                      done = 1;
                  }
                  else {
                      if (vf > vf_to_match) {
                          vf1 = vf;
                          ds1 = *distance;
                      }
                      else {
                          vf0 = vf;
                          ds0 = *distance;
                      }
                      *distance = 0.5 *(ds0 + ds1);
                      niter++;
                  }
              }   // else   
         }        // while  
//       assert(niter < niter_mx);
         if (niter >= niter_mx) {
	     printf("WARNING: %d iteration reached\n", niter_mx);
	 }
     }            // else  

     free(node_loc);
     free(node_order_for_ds);
     free(ds_ea_node); 
     free(coords_work); 

     return;
 } 

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
               int *nnode_upper, int **nodelist_upper)
{ 
//   coords_new:         output, but is assumed to have been allocated long enough  
//   nodes_interface:    output, but is assumed to have been allocated as 2 long  
//   nodelist_lower:     output
//   nodelist_upper:     output

     int szdim;
     int dim, i, n, n0, n1, nn, idx, idx_next, offset, ndistance, found; 
     int *node_loc, *nnode_ea_ds, *nodelist_ea_ds;

     double ds, vol, vf;
     double *c0, *c1, *ds_ordered, *coords_work;  

     dim = 2;
     szdim = dim * sizeof(double); 

     ds_ordered     = (double *) malloc(nnode * sizeof(double));
     node_loc       = (int    *) malloc((3 * nnode) * sizeof(int));
     nnode_ea_ds    = node_loc    + nnode;
     nodelist_ea_ds = nnode_ea_ds + nnode;

     coords_work = (double *) malloc((nnode + nnode) * szdim);

     n  = node_order_for_ds[0];
     ds = ds_ea_node[n];
     ds_ordered[0] = ds;
     nnode_ea_ds[0] = 1;
     nodelist_ea_ds[0] = n;
     offset    = 1;
     ndistance = 1;
     idx_next  = 1;
     found     = 1;

     *nnode_lower = 0;
     *nnode_upper = 0;
     *nodelist_lower = NULL;
     *nodelist_upper = NULL; 

     while (found) {
         found = 0;
         for (i = idx_next; i < nnode; i++) {
             n = node_order_for_ds[i];
             if (ds_ea_node[n] - ds > accuracy) { 
                 ds = ds_ea_node[n];
                 ds_ordered[ndistance]  = ds;
                 nnode_ea_ds[ndistance] = 1;
                 nodelist_ea_ds[offset] = n;

                 ndistance++;
                 offset++;
                 idx_next = i + 1;
                 found    = 1;
                 break;
             }
             else {
                 nnode_ea_ds[ndistance-1]++;
                 nodelist_ea_ds[offset] = n;
                 offset++;
             }
         }
     }
     idx          = 1;
     *ds_lower  = ds_ordered[0];
     *vf_lower  = 0.0;
     *vol_matched = 0;
 
     while (idx < ndistance) { 

         if (*nodelist_lower) { 
             free(*nodelist_lower);
             *nodelist_lower = NULL;
         }
         if (*nodelist_upper) { 
             free(*nodelist_upper);
             *nodelist_upper = NULL;
         }
         ds = ds_ordered[idx];

         for (i = 0; i < nnode; i++) {
             node_loc[i] = 1;          // above the interface    
         }
         offset = 0;
         for (i = 0; i < idx; i++) {
             offset += nnode_ea_ds[i];
         }
         for (i = 0; i < offset; i++) {
             n = nodelist_ea_ds[i];

assert(n < nnode);

             node_loc[n] = -1;         // below the interface 
         }
         nn = nnode_ea_ds[idx];
         for (i = 0; i < nn; i++) {
             n  = nodelist_ea_ds[offset+i];
             nodes_interface[i] = n;

assert(n < nnode);

             node_loc[n] = 0;          // on the interface  
         }
         if (idx == ndistance - 1) { 
             *nnode_new = 0;
       
             *nnode_interface = nn;
             *nnode_upper = 0;
             *nnode_lower = nnode;
             *nodelist_lower = (int *) malloc(nnode * sizeof(int));
             memcpy(*nodelist_lower, nodelist, (size_t)(nnode * sizeof(int))); 
             vf = 1.0; 
         } 
         else { 
             find_interface2d(nnode, coords, nodelist, norm, ds, node_loc,
                           nnode_new, coords_new, 
                           nnode_interface, nodes_interface,
                           nnode_lower, nodelist_lower, 
                           nnode_upper, nodelist_upper);  
    
    //       calculate the volume of the lower polygon
    
             vol = 0.0;
             for (i = 0; i < *nnode_lower; i++) { 
                 n0 = (*nodelist_lower)[i];
                 if (n0 < nnode) { 
                     c0 = coords + (n0 + n0);
                 }
                 else { 
                     n1 = n0 - nnode;
                     c0 = coords_new + (n1 + n1);
                 }
                 c1 = coords_work + (i + i);
                 memcpy(c1, c0, (size_t)szdim); 
             }
             if (geop == 1) { 
                 cal_poly_area(*nnode_lower, coords_work, *nnode_lower, NULL, &vol);
             }
             else if (geop == 2) { 
//               comment out for TT effort  
//                 rz_area(*nnode_lower, coords_work, &vol);
                 printf("ERROR: activate the line above\n");
             } 
             vf = vol/volume;
         }                       // idx  != ndistance - 1 
 
         if (fabs(vf - vf_to_match) < accuracy) {
             *vol_matched = 1;
             *ds_lower    = ds;
             *ds_upper    = ds;

             *vf_lower = vf;
             *vf_upper = vf;
             break;
         }
         else if (idx == ndistance-1) {
             *ds_upper = ds;
             *vf_upper = vf;
             break;
         }
         else if (vf > vf_to_match) {
             *ds_upper = ds;
             *vf_upper = vf;
             break;
         }
         else if (vf < vf_to_match) {
             *vf_lower = vf;
             *ds_lower = ds;
             idx++;
         }
     }                   
     free(ds_ordered); 
     free(node_loc);  
     free(coords_work);
     
     return;
   } 

void find_interface2d(int nnode, double *coords, int *nodelist, 
                      double *norm, double distance, 
                      int *node_loc, 
                      int *nnode_new, double *coords_new, 
                      int *nnode_interface, int *nodes_interface,
                      int *nnode_lower, int **nodelist_lower, 
                      int *nnode_upper, int **nodelist_upper)
{ 
//   node_loc[0:nnode) : input, node_loc[n] =  0 if the node n is on the interface 
//                              node_loc[n] = -1 if the node n is belower the interface
//                              node_loc[n] =  1 if the node n is above the interface  
//
//   nodes_interface:    output, but is assumed to have been allocated as 2 long  
//
//   if node (= nodelist_lower[i] or nodelist_upper[i]) < nnode, node is one of 
//   the originel nnode nodes. If node >= nnode, it is (node - nnode)th new node.

     int i0, i1, n0, n1, nn;
     int *nodelist_low, *nodelist_up; 
     double t, *c0, *c1, *ci; 

     *nnode_interface = 0; 
     *nnode_new       = 0;
     *nnode_lower     = 0;
     *nnode_upper     = 0;
     *nodelist_lower  = NULL;
     *nodelist_upper  = NULL;
     
//   nodelist_low and nodelist_low are working arrays, each is assumed to be
//   at least (nnode + nnode) long. 

     nn = nnode + nnode;
     nodelist_low = (int *) malloc((nn + nn) * sizeof(int));
     nodelist_up  = nodelist_low + nn;
     
     for (i0 = 0; i0 < nnode; i0++) { 
         i1 = (i0 + 1) % nnode;   
         n1 = nodelist[i1];
         n0 = nodelist[i0];

         if (!node_loc[n0] && !node_loc[n1]) { // both points are on the line 
             nodelist_low[*nnode_lower]     = n0;
             nodelist_low[(*nnode_lower)+1] = n1;
             *nnode_lower += 2;
             nodelist_up[*nnode_upper]     = n0;
             nodelist_up[(*nnode_upper)+1] = n1;
             *nnode_upper += 2;

             nodes_interface[0] = n0;
             nodes_interface[1] = n1; 
             *nnode_interface = 2;
         }
         else if (!node_loc[n0]) { 
             nodelist_low[*nnode_lower] = n0; 
             nodelist_up[ *nnode_upper] = n0; 
             (*nnode_lower)++;
             (*nnode_upper)++;
 
             if (*nnode_interface) { 
                 if (nodes_interface[(*nnode_interface)-1] != n0) { 
                     nodes_interface[*nnode_interface]      = n0;
                     (*nnode_interface)++; 
                 }
             } 
             else { 
                 nodes_interface[*nnode_interface] = n0;
                 (*nnode_interface)++;
             } 
             if (node_loc[n1] < 0) { 
                 nodelist_low[*nnode_lower] = n1; 
                 (*nnode_lower)++;
             } 
             else if (node_loc[n1] > 0) { 
                 nodelist_up[*nnode_upper] = n1;
                 (*nnode_upper)++;
             } 
         }
         else if (!node_loc[n1]) {  
              
             if (node_loc[n0] < 0) {
                 nodelist_low[*nnode_lower] = n0;
                 (*nnode_lower)++; 
             } 
             else if (node_loc[n0] > 0) {
                 nodelist_up[*nnode_upper] = n0;
                 (*nnode_upper)++;
             }
             nodelist_low[*nnode_lower] = n1;
             nodelist_up[ *nnode_upper] = n1;  
             (*nnode_lower)++;
             (*nnode_upper)++; 

             if (*nnode_interface < 2) { 
                 if (*nnode_interface) { 
                     if (nodes_interface[(*nnode_interface)-1] != n1) { 
                         nodes_interface[*nnode_interface]      = n1;
                         (*nnode_interface)++;
                     } 
                 }
                 else { 
                     nodes_interface[*nnode_interface] = n1;
                     (*nnode_interface)++;
                 } 
             } 
         } 
         else if ((node_loc[n0] > 0) && (node_loc[n1] > 0)) {  
             nodelist_up[*nnode_upper]     = n0;
             nodelist_up[(*nnode_upper)+1] = n1;
             *nnode_upper += 2; 
         }
         else if ((node_loc[n0] < 0) && (node_loc[n1] < 0)) { 
             nodelist_low[*nnode_lower]     = n0;
             nodelist_low[(*nnode_lower)+1] = n1;
             *nnode_lower += 2;
         }
         else {   // find the intersection 
             c0 = coords + (n0 + n0);
             c1 = coords + (n1 + n1);
             t = (distance - (norm[0] * c0[0] + norm[1] * c0[1]))/
                 (norm[0] * (c1[0] - c0[0]) + norm[1] * (c1[1] - c0[1]));
             t = MAX(0.0, MIN(1.0, t)); 
             ci = coords_new + (*nnode_new + *nnode_new); 
             ci[0] = c0[0] + (c1[0] - c0[0]) * t;
             ci[1] = c0[1] + (c1[1] - c0[1]) * t;  

             if (node_loc[n0] < 0) { 
                 nodelist_low[*nnode_lower]     = n0; 
                 nodelist_low[(*nnode_lower)+1] = nnode + *nnode_new;
                 *nnode_lower += 2; 

                 nodelist_up[*nnode_upper]     = nnode + *nnode_new;
                 nodelist_up[(*nnode_upper)+1] = n1;
                 *nnode_upper += 2;
             }
             else if (node_loc[n1] < 0) { 
                 nodelist_low[*nnode_lower]     = nnode + *nnode_new;
                 nodelist_low[(*nnode_lower)+1] = n1;
                 *nnode_lower += 2; 
             
                 nodelist_up[*nnode_upper]     = n0;
                 nodelist_up[(*nnode_upper)+1] = nnode + *nnode_new;
                 *nnode_upper += 2;
             }
             nodes_interface[*nnode_interface] = nnode + *nnode_new;
             (*nnode_interface)++;

             (*nnode_new)++; 
         }  
     } 
     if (nodelist_low[0] == nodelist_low[(*nnode_lower)-1]) { 
         (*nnode_lower)--;
     }
     if (nodelist_up[0] == nodelist_up[(*nnode_upper)-1]) { 
         (*nnode_upper)--;
     } 
//       remove the reduntant codes

     nn = 1;
     for (i0 = 1; i0 < *nnode_lower; i0++) { 
         if (nodelist_low[i0] != nodelist_low[nn-1]) { 
             nodelist_low[nn]  = nodelist_low[i0];
             nn++; 
         }
     }  
     *nnode_lower = nn; 
     if (*nnode_lower < 3) { 
         *nnode_lower = 0;
     }
     else { 
         *nodelist_lower = (int *) malloc((*nnode_lower) * sizeof(int));
         memcpy(*nodelist_lower, nodelist_low, (size_t)((*nnode_lower) * sizeof(int)));
     }  
     nn = 1; 
     for (i0 = 1; i0 < *nnode_upper; i0++) {
         if (nodelist_up[i0] != nodelist_up[nn-1]) {                    
             nodelist_up[nn]  = nodelist_up[i0];
             nn++;
         }
     }
     *nnode_upper = nn;
     if (*nnode_upper < 3) { 
         *nnode_upper = 0;
     }
     else { 
         *nodelist_upper = (int *) malloc((*nnode_upper) * sizeof(int));
         memcpy(*nodelist_upper, nodelist_up, (size_t)((*nnode_upper) * sizeof(int))); 
     }
     free(nodelist_low);

     return;
  } 



void remap2d_scaled(int ifdump,
                    int nnode, double *coords, double *inward_norm_ea_face,
                    int nnode_m, double *coords_m, 
                    int *nnodex, double **coordsx)
{ 
     int dim, n, i;
     double dx_mx, tmp; 
     double *coords_scaled, *coords_m_scaled, *c0, *c, *ct; 

     dim = 2;

     coords_scaled   = (double *) malloc((nnode + nnode_m) * dim * sizeof(double));
     coords_m_scaled = coords_scaled + (nnode * dim);  

     c0 = coords; 
     dx_mx = 0;
     for (n = 1; n < nnode; n++) { 
         c = coords + (n * dim);
         for (i = 0; i < dim; i++) { 
             dx_mx = MAX(dx_mx, fabs(c[i] - c0[i]));
         }
     }

//   shift and scale the coordinates

     tmp = 1.0/dx_mx; 
     for (n = 0; n < nnode; n++) { 
         c  = coords + (n * dim);
         ct = coords_scaled + (n * dim);
         for (i = 0; i < dim; i++) { 
             ct[i] = tmp * (c[i] - c0[i]);
         }
     }
     for (n = 0; n < nnode_m; n++) {
         c  = coords_m + (n * dim);
         ct = coords_m_scaled + (n * dim);
         for (i = 0; i < dim; i++) {
             ct[i] = tmp * (c[i] - c0[i]);
         }
     }
     remap2d_mat(ifdump, nnode, coords_scaled, inward_norm_ea_face, 
                 nnode_m, coords_m_scaled, 
                 nnodex, coordsx);

//   shift and scale back 

     for (n = 0; n < *nnodex; n++) { 
         c = *coordsx + (n * dim);
         for (i = 0; i < dim; i++) { 
             c[i] = c0[i] + dx_mx * c[i];
         }
     }
     free(coords_scaled); 

     return;
 }  
    
                       


void remap2d_mat(int ifdump, 
                 int nnode, double *coords, double *inward_norm_ea_face, 
                 int nnode_m, double *coords_m,   
                 int *nnodex, double **coordsx)  
{ 
    char name[64];

    int dim, szint, szdim, n, n1, n0, nnode_mx, i, k, fileid;
    int my_nnode, nnode_new, nnode_lower, nnode_upper, nnode_interface;
    int nodes_interface[2]; 
    int *nodelist_upper, *nodelist_lower, *nodelist_default, *node_loc;
    
    double distance; 
    double *ds_ea_node, *coords_work, *coords_out, *coords_new, *my_coords;
    double *norm, *c0, *c1, *c;  

    dim = 2;
    szdim = dim * sizeof(double);
    szint = sizeof(int);

    *nnodex  = 0;
    *coordsx = NULL; 

    nnode_mx = MAX(nnode, nnode_m);

//  nodelist_lower and nodelist_upper need to be at least 2 * nnode_mx long.
//  set 3 * nnode_mx since thr resulting polygon could have more than nnode_mx nodes. 

    n = nnode_mx + nnode_mx + nnode_mx;
    nodelist_default = (int *) malloc((n + n) * sizeof(int));
    node_loc         = nodelist_default + n;

    for (i = 0; i < n; i++) { 
        nodelist_default[i] = i;
    } 
    ds_ea_node = (double *) malloc(n *(sizeof(double) + 3 * szdim));
    coords_work = ds_ea_node  + n; 
    coords_out  = coords_work + (dim * n); 
    coords_new  = coords_out + (dim * n); 

    my_nnode    = nnode_m;
    my_coords   = coords_m; 

if (ifdump) { 
// sprintf(name,  "%s", "file_2dmap");
// mio_open_file(name, mio_file_create, &fileid);
// sprintf(name,  "%s", "poly0");
// write_a_polygon(fileid, name, 2, nnode, coords, nnode, nodelist_default);
// sprintf(name,  "%s", "polym");
// write_a_polygon(fileid, name, 2, nnode_m, coords_m, nnode_m, nodelist_default);
} 


    nodelist_lower = NULL;
    nodelist_upper = NULL;

    for (n = 0; n < nnode; n++) { 
        n1 = (n + 1) % nnode;
        norm = inward_norm_ea_face + (n * dim);

        c  = coords + (dim * n);
        c1 = coords + (dim * n1);
        distance = 0.5 *(norm[0] *(c[0] + c1[0]) + norm[1] * (c[1] + c1[1]));
        
        for (i = 0; i < my_nnode; i++) { 
            ds_ea_node[i] = 0.0;
            c = my_coords + (i * dim);
            for (k = 0; k < dim; k++) {  
                ds_ea_node[i] += (norm[k] * c[k]);
            }
        }
        for (i = 0; i < my_nnode; i++) {   
            if (fabs(ds_ea_node[i] - distance) <= accuracy) {
                node_loc[i] = 0;  // on the plane
            }
            else if (ds_ea_node[i] > distance)  {   // above the plane 
                node_loc[i] = 1;
            }
            else {
                node_loc[i] = -1;
            }
        } 
        find_interface2d(my_nnode, my_coords, nodelist_default, norm,
                         distance, node_loc,
                         &nnode_new, coords_new,
                         &nnode_interface, nodes_interface,
                         &nnode_lower, &nodelist_lower,
                         &nnode_upper, &nodelist_upper);

        for (i = 0; i < nnode_upper; i++) { 
            n0 = nodelist_upper[i]; 
            if (n0 < my_nnode) { 
                c0 = my_coords + (n0 + n0);
            }
            else { 
                n0 -= my_nnode;
                c0 = coords_new + (n0 + n0); 
            }
            c1 = coords_out + (i + i);
            memcpy(c1, c0, (size_t)szdim); 
        } 
        if (nnode_upper < 3) { 
            my_nnode = 0;
            break;
        }
        else { 
            my_nnode = nnode_upper;
            memcpy(coords_work, coords_out, (size_t)(nnode_upper * szdim));
            my_coords = coords_work;
            

if (ifdump) { 
// sprintf(name,  "poly_cut_f_%d", n);
// write_a_polygon(fileid, name, 2, my_nnode, my_coords, my_nnode, nodelist_default);
} 



        }
        if (nodelist_lower) free(nodelist_lower); 
        if (nodelist_upper) free(nodelist_upper);
        nodelist_lower = NULL;
        nodelist_upper = NULL;

    } 

if (ifdump) {
// mio_close_file(fileid);
fileid = -1;
} 





    *nnodex = my_nnode;
    if (my_nnode >= 3) { 
        *coordsx = (double *) malloc(my_nnode * szdim); 
        memcpy(*coordsx, my_coords, (size_t)(my_nnode * szdim));
    }
    free(ds_ea_node);
    free(nodelist_default); 

    return;
 }  
