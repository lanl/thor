
#ifdef MPI
#include "mpi.h"
#endif

#include "minip.h"
#include "bdry.h" 



void bdry_cell_2d(int nmat, int *ncell, int nbdry, 
                    Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                    double ***vf_2dmat, double ***rho_2dmat, 
                    double ***ei_2dmat, double ***pres_2dmat,
                    double **rho_2dcell, double **ei_2dcell, double **pres_2dcell) 
{

     int dim, i, i0, i1, j, j0, j1, dir, n0, n1, m;
     int ncell_ext[2];

     dim = 2;
     for (i = 0; i < dim; i++) { 
         ncell_ext[i] = ncell[i] + nbdry + nbdry;
     }
//   xl edge without corners 
     if (btype_lower[0] == bdry_transmitted) {
         i0 = nbdry; 
         for (j = nbdry; j < ncell[1]+nbdry; j++) { 
             for (i = 0; i < nbdry; i++) { 
                 rho_2dcell[j][i] = 0.0;
                 ei_2dcell[j][i] = 0.0; 
                 pres_2dcell[j][i] = 0.0;                  

                 for (m = 0; m < nmat; m++) { 
                     vf_2dmat[j][i][m]   = vf_2dmat[j][i0][m];     
                     rho_2dmat[j][i][m]  = rho_2dmat[j][i0][m]; 
                     ei_2dmat[j][i][m]   = ei_2dmat[j][i0][m]; 
                     pres_2dmat[j][i][m] = pres_2dmat[j][i0][m]; 

                     rho_2dcell[j][i] += (vf_2dmat[j][i][m] * rho_2dmat[j][i][m]);
                     ei_2dcell[j][i]  += (vf_2dmat[j][i][m] * ei_2dmat[j][i][m]);
                     pres_2dcell[j][i]+= (vf_2dmat[j][i][m] * pres_2dmat[j][i][m]); 
                 } 
             }
         }
     } 
     else { 
     
     }
//   xr edge without corners 
     if (btype_upper[0] == bdry_transmitted) { 
         i0 = ncell[0] + nbdry - 1; 
         for (j = nbdry; j < ncell[1]+nbdry; j++) {
             for (i = 0; i < nbdry; i++) { 
                 i1 = i0 + 1 + i;
                 rho_2dcell[j][i1]  = 0.0;
                 ei_2dcell[j][i1]   = 0.0;
                 pres_2dcell[j][i1] = 0.0;

                 for (m = 0; m < nmat; m++) { 
                     vf_2dmat[j][i1][m]   = vf_2dmat[j][i0][m];
                     rho_2dmat[j][i1][m]  = rho_2dmat[j][i0][m]; 
                     ei_2dmat[j][i1][m]   = ei_2dmat[j][i0][m];
                     pres_2dmat[j][i1][m] = pres_2dmat[j][i0][m];

                     rho_2dcell[j][i1] += (vf_2dmat[j][i1][m] * rho_2dmat[j][i1][m]);
                     ei_2dcell[j][i1]  += (vf_2dmat[j][i1][m] * ei_2dmat[j][i1][m]);
                     pres_2dcell[j][i1]+= (vf_2dmat[j][i1][m] * pres_2dmat[j][i1][m]); 
                 } 
             }
         }
     }
     else { 

     }
//   yl edge including corners 
     if (btype_lower[1] == bdry_transmitted) {  
         j0 = nbdry; 
         for (j = 0; j < nbdry; j++) { 
             for (i = 0; i < ncell_ext[0]; i++) { 
                 rho_2dcell[j][i] = 0.0;
                 ei_2dcell[j][i] = 0.0;
                 pres_2dcell[j][i] = 0.0;
 
                 for (m = 0; m < nmat; m++) {
                     vf_2dmat[j][i][m]  = vf_2dmat[j0][i][m];
                     rho_2dmat[j][i][m] =  rho_2dmat[j0][i][m];
                      ei_2dmat[j][i][m] =   ei_2dmat[j0][i][m];
                    pres_2dmat[j][i][m] = pres_2dmat[j0][i][m];

                    rho_2dcell[j][i] += (vf_2dmat[j][i][m] * rho_2dmat[j][i][m]);
                     ei_2dcell[j][i]  += (vf_2dmat[j][i][m] * ei_2dmat[j][i][m]);
                     pres_2dcell[j][i]+= (vf_2dmat[j][i][m] * pres_2dmat[j][i][m]);
                 }
             } 
         } 
     }
     else {

     }
//   yr edge including corners  
     if (btype_upper[1] == bdry_transmitted) {
         j0 = ncell[1] + nbdry - 1; 
         for (j = 0; j < nbdry; j++) {  
             j1 = j0 + 1 + j;
             for (i = 0; i < ncell_ext[0]; i++) {  
                 rho_2dcell[j1][i]  = 0.0;
                 ei_2dcell[j1][i]   = 0.0;
                 pres_2dcell[j1][i] = 0.0;

                 for (m = 0; m < nmat; m++) {  
                     vf_2dmat[j1][i][m]   = vf_2dmat[j0][i][m];
                      rho_2dmat[j1][i][m] = rho_2dmat[j0][i][m];
                       ei_2dmat[j1][i][m] = ei_2dmat[j0][i][m];
                     pres_2dmat[j1][i][m] = pres_2dmat[j0][i][m];

                     rho_2dcell[j1][i] += (vf_2dmat[j1][i][m] * rho_2dmat[j1][i][m]);
                     ei_2dcell[j1][i]  += (vf_2dmat[j1][i][m] * ei_2dmat[j1][i][m]);
                     pres_2dcell[j1][i]+= (vf_2dmat[j1][i][m] * pres_2dmat[j1][i][m]);
                 }
             } 
         }
     }
     else {

     }

     return;
 }

void bdry_cell_3d(int nmat, int *ncell, int nbdry, 
                Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                double ****vf_3dmat, double ****rho_3dmat, 
                double ****ei_3dmat, double ****pres_3dmat,
                double ***rho_3dcell, double ***ei_3dcell, double ***pres_3dcell)
{
     int dim, i, i0, i1, j, j0, j1, k, k0, k1, dir, n0, n1, m;
     int ncell_ext[3];

     dim = 3; 
     for (i = 0; i < dim; i++) { 
         ncell_ext[i] = ncell[i] + nbdry + nbdry;
     }
//   xl face without edeges 
     if (btype_lower[0] == bdry_transmitted) {
         i0 = nbdry; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) { 
             for (j = nbdry; j < ncell[1]+nbdry; j++) { 
                 for (i = 0; i < nbdry; i++) { 
                      rho_3dcell[k][j][i]   = 0.0;
                      ei_3dcell[k][j][i]    = 0.0;
                      pres_3dcell[k][j][i]  = 0.0; 
                      for (m = 0; m < nmat; m++) { 
                          vf_3dmat[k][j][i][m]   = vf_3dmat[k][j][i0][m];     
                          rho_3dmat[k][j][i][m]  = rho_3dmat[k][j][i0][m]; 
                          ei_3dmat[k][j][i][m]   = ei_3dmat[k][j][i0][m]; 
                          pres_3dmat[k][j][i][m] = pres_3dmat[k][j][i0][m]; 

                          rho_3dcell[k][j][i]   += (vf_3dmat[k][j][i][m] * rho_3dmat[k][j][i][m]);
                          ei_3dcell[k][j][i]    += (vf_3dmat[k][j][i][m] *  ei_3dmat[k][j][i][m]);
                          pres_3dcell[k][j][i]  += (vf_3dmat[k][j][i][m] * pres_3dmat[k][j][i][m]);
                      }
                  }
              }
         }
     } 
     else { 
     
     }
//   xr face without edeges 
     if (btype_upper[0] == bdry_transmitted) { 
         i0 = ncell[0] + nbdry - 1; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) { 
             for (j = nbdry; j < ncell[1]+nbdry; j++) {
                 for (i = 0; i < nbdry; i++) { 
                     i1 = i0 + 1 + i;

                     rho_3dcell[k][j][i1]   = 0.0;
                     ei_3dcell[k][j][i1]    = 0.0;
                     pres_3dcell[k][j][i1]  = 0.0;
       
                     for (m = 0; m < nmat; m++) {  
                         vf_3dmat[k][j][i1][m]   = vf_3dmat[k][j][i0][m];
                         rho_3dmat[k][j][i1][m]  = rho_3dmat[k][j][i0][m]; 
                         ei_3dmat[k][j][i1][m]   = ei_3dmat[k][j][i0][m];
                         pres_3dmat[k][j][i1][m] = pres_3dmat[k][j][i0][m];

                         rho_3dcell[k][j][i1]   += (vf_3dmat[k][j][i1][m] * rho_3dmat[k][j][i1][m]);
                          ei_3dcell[k][j][i1]    += (vf_3dmat[k][j][i1][m] *  ei_3dmat[k][j][i1][m]);
                          pres_3dcell[k][j][i1]  += (vf_3dmat[k][j][i1][m] * pres_3dmat[k][j][i1][m]);
                     }
                 }
             } 
         }
     }
     else { 

     }
//   yl face with edges  
     if (btype_lower[1] == bdry_transmitted) {  
         j0 = nbdry; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) {
             for (j = 0; j < nbdry; j++) { 
                 for (i = 0; i < ncell_ext[0]; i++) { 
                     rho_3dcell[k][j][i]   = 0.0;
                     ei_3dcell[k][j][i]    = 0.0;
                     pres_3dcell[k][j][i]  = 0.0;

                     for (m = 0; m < nmat; m++) { 
                         vf_3dmat[k][j][i][m]   = vf_3dmat[k][j0][i][m];
                         rho_3dmat[k][j][i][m]  = rho_3dmat[k][j0][i][m];
                         ei_3dmat[k][j][i][m]   = ei_3dmat[k][j0][i][m];
                         pres_3dmat[k][j][i][m] = pres_3dmat[k][j0][i][m];

                         rho_3dcell[k][j][i]   += (vf_3dmat[k][j][i][m] * rho_3dmat[k][j][i][m]);
                         ei_3dcell[k][j][i]    += (vf_3dmat[k][j][i][m] *  ei_3dmat[k][j][i][m]);
                         pres_3dcell[k][j][i]  += (vf_3dmat[k][j][i][m] * pres_3dmat[k][j][i][m]);
                     }
                 } 
             } 
         } 
     }
     else {

     }
//   yr face with edge 
     if (btype_upper[1] == bdry_transmitted) {
         j0 = ncell[1] + nbdry - 1; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) {
             for (j = 0; j < nbdry; j++) {  
                 j1 = j0 + 1 + j;
                 for (i = 0; i < ncell_ext[0]; i++) {  
                     rho_3dcell[k][j1][i]   = 0.0;
                     ei_3dcell[k][j1][i]    = 0.0;
                     pres_3dcell[k][j1][i]  = 0.0;

                     for (m = 0; m < nmat; m++) { 
                         vf_3dmat[k][j1][i][m]   = vf_3dmat[k][j0][i][m];
                         rho_3dmat[k][j1][i][m]  = rho_3dmat[k][j0][i][m];
                         ei_3dmat[k][j1][i][m]   = ei_3dmat[k][j0][i][m];
                         pres_3dmat[k][j1][i][m] = pres_3dmat[k][j0][i][m];

                         rho_3dcell[k][j1][i]   += (vf_3dmat[k][j1][i][m] * rho_3dmat[k][j1][i][m]);
                         ei_3dcell[k][j1][i]    += (vf_3dmat[k][j1][i][m] *  ei_3dmat[k][j1][i][m]);
                         pres_3dcell[k][j1][i]  += (vf_3dmat[k][j1][i][m] * pres_3dmat[k][j1][i][m]);
                     }
                 }
             }
         }
     }
     else {

     }
//   zl face including edges and corners
     if (btype_lower[2] == bdry_transmitted) {
         k0 = nbdry;
         for (k = 0; k < nbdry; k++) {
             for (j = 0; j < ncell_ext[1]; j++) {
                 for (i = 0; i < ncell_ext[0]; i++) {
                     rho_3dcell[k][j][i]   = 0.0;
                     ei_3dcell[k][j][i]    = 0.0;
                     pres_3dcell[k][j][i]  = 0.0;

                     for (m = 0; m < nmat; m++) { 
                         vf_3dmat[k][j][i][m]   = vf_3dmat[k0][j][i][m];
                         rho_3dmat[k][j][i][m]  = rho_3dmat[k0][j][i][m];
                         ei_3dmat[k][j][i][m]   = ei_3dmat[k0][j][i][m];
                         pres_3dmat[k][j][i][m] = pres_3dmat[k0][j][i][m];

                         rho_3dcell[k][j][i]   += (vf_3dmat[k][j][i][m] * rho_3dmat[k][j][i][m]);
                         ei_3dcell[k][j][i]    += (vf_3dmat[k][j][i][m] *  ei_3dmat[k][j][i][m]);
                         pres_3dcell[k][j][i]  += (vf_3dmat[k][j][i][m] * pres_3dmat[k][j][i][m]);
                     }
                 }
             } 
         }
     }
     else {

     }
//   zr face including edges and corners 
     if (btype_upper[2] == bdry_transmitted) {
         k0 = ncell[2] + nbdry - 1;
         for (k = 0; k < nbdry; k++) { 
             k1 = k0 + 1 + k;
             for (j = 0; j < ncell_ext[1]; j++) {
                 for (i = 0; i < ncell_ext[0]; i++) {
                     rho_3dcell[k1][j][i]   = 0.0;
                     ei_3dcell[k1][j][i]    = 0.0;
                     pres_3dcell[k1][j][i]  = 0.0;

                     for (m = 0; m < nmat; m++) { 
                         vf_3dmat[k1][j][i][m]   = vf_3dmat[k0][j][i][m];
                         rho_3dmat[k1][j][i][m]  = rho_3dmat[k0][j][i][m];
                         ei_3dmat[k1][j][i][m]   = ei_3dmat[k0][j][i][m];
                         pres_3dmat[k1][j][i][m] = pres_3dmat[k0][j][i][m];

                         rho_3dcell[k1][j][i]   += (vf_3dmat[k1][j][i][m] * rho_3dmat[k1][j][i][m]);
                         ei_3dcell[k1][j][i]    += (vf_3dmat[k1][j][i][m] *  ei_3dmat[k1][j][i][m]);
                         pres_3dcell[k1][j][i]  += (vf_3dmat[k1][j][i][m] * pres_3dmat[k1][j][i][m]);
                     }
                 } 
             }
         }
     }
     else {

     }

     return;
 } 

void bdry_node_2d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                  double ***vel_for_2dnode)
{
     int dim, i, i0, i1, j, j0, j1, dir, n0, n1;
     int ncell_ext[2];

     dim = 2;
     for (i = 0; i < dim; i++) { 
         ncell_ext[i] = ncell[i] + nbdry + nbdry;
     }
//   xl edge without corners 
     if (btype_lower[0] == bdry_transmitted) {
         n0 = nbdry; 
         for (j = nbdry; j <= ncell[1]+nbdry; j++) {
             for (i = 0; i < nbdry; i++) { 
                 for (dir = 0; dir < dim; dir++) { 
                     vel_for_2dnode[j][i][dir] = vel_for_2dnode[j][n0][dir]; 
                 } 
             }
         }          
     } 
     else { 
     
     }
//   xr edge without corners 
     if (btype_upper[0] == bdry_transmitted) { 
         n0 = ncell[0] + nbdry; 
         for (j = nbdry; j <= ncell[1]+nbdry; j++) {
             for (i = 0; i < nbdry; i++) {       
                 n1 = n0 + 1 + i; 
                 for (dir = 0; dir < dim; dir++) {
                     vel_for_2dnode[j][n1][dir] = vel_for_2dnode[j][n0][dir];
                 }
             }
         }
     }
     else { 

     }
//   yl edge including corners 
     if (btype_lower[1] == bdry_transmitted) {  
         n0 = nbdry; 
         for (j = 0; j < nbdry; j++) {  
             for (i = 0; i <= ncell_ext[0]; i++) {   
                 for (dir = 0; dir < dim; dir++) {
                     vel_for_2dnode[j][i][dir] = vel_for_2dnode[n0][i][dir]; 
                 }
             } 
         } 
     }
     else {

     }
//   yr edge including corners  
     if (btype_upper[1] == bdry_transmitted) {
         n0 = ncell[1] + nbdry;
         for (j = 0; j < nbdry; j++) {
             n1 = n0 + 1 + j; 
             for (i = 0; i <= ncell_ext[0]; i++) {  
                 for (dir = 0; dir < dim; dir++) { 
                     vel_for_2dnode[n1][i][dir] = vel_for_2dnode[n0][i][dir]; 
                 } 
             }
         }
     }
     else {

     }

     return;
 } 

void bdry_node_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                  double ****vel_for_3dnode)
{
     int dim, i, i0, i1, j, j0, j1, k, k0, k1, dir, n0, n1;
     int ncell_ext[3];

     dim = 3; 
     for (i = 0; i < dim; i++) { 
         ncell_ext[i] = ncell[i] + nbdry + nbdry;
     }
//   xl face without edeges 
     if (btype_lower[0] == bdry_transmitted) {
         n0 = nbdry; 
         for (k = nbdry; k <= ncell[2]+nbdry; k++) { 
             for (j = nbdry; j <= ncell[1]+nbdry; j++) {
                 for (i = 0; i < nbdry; i++) { 
                     for (dir = 0; dir < dim; dir++) { 
                         vel_for_3dnode[k][j][i][dir] = vel_for_3dnode[k][j][n0][dir]; 
                     }
                 } 
             }
         }          
     } 
     else { 
     
     }
//   xr face without edeges 
     if (btype_upper[0] == bdry_transmitted) { 
         n0 = ncell[0] + nbdry;
         for (k = nbdry; k <= ncell[2]+nbdry; k++) { 
             for (j = nbdry; j <= ncell[1]+nbdry; j++) {
                 for (i = 0; i < nbdry; i++) {       
                     n1 = n0 + 1 + i; 
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dnode[k][j][n1][dir] = vel_for_3dnode[k][j][n0][dir];
                     }
                 }
             }
         }
     }
     else { 

     }
//   yl face with edges  
     if (btype_lower[1] == bdry_transmitted) {  
         n0 = nbdry; 
         for (k = nbdry; k <= ncell[2]+nbdry; k++) {
             for (j = 0; j < nbdry; j++) {  
                 for (i = 0; i <= ncell_ext[0]; i++) {   
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dnode[k][j][i][dir] = vel_for_3dnode[k][n0][i][dir]; 
                     }
                 }
             } 
         } 
     }
     else {

     }
//   yr face with edge 
     if (btype_upper[1] == bdry_transmitted) {
         n0 = ncell[1] + nbdry;
         for (k = nbdry; k <= ncell[2]+nbdry; k++) { 
             for (j = 0; j < nbdry; j++) {
                 n1 = n0 + 1 + j; 
                 for (i = 0; i <= ncell_ext[0]; i++) {  
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dnode[k][n1][i][dir] = vel_for_3dnode[k][n0][i][dir]; 
                     }
                 }
             }
         }
     }
     else {

     }
//   zl face including edges and corners
     if (btype_lower[2] == bdry_transmitted) {
         n0 = nbdry;
         for (k = 0; k < nbdry; k++) {
             for (j = 0; j <= ncell_ext[1]; j++) {
                 for (i = 0; i <= ncell_ext[0]; i++) { 
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dnode[k][j][i][dir] = vel_for_3dnode[n0][j][i][dir];
                     }
                 }
             }
         }
     }
     else {

     }
//   zr face including edges and corners 
     if (btype_upper[2] == bdry_transmitted) {
         n0 = ncell[2] + nbdry;
         for (k = 0; k < nbdry; k++) {
             n1 = n0 + 1 + k; 
             for (j = 0; j <= ncell_ext[1]; j++) {
                 for (i = 0; i <= ncell_ext[0]; i++) {
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dnode[n1][j][i][dir] = vel_for_3dnode[n0][j][i][dir];
                     }
                 }
             }
         }
     }
     else {

     }

     return;
 } 

void bdry_cell_vel_2d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                      double ***vel_for_2dcell)
{
     int dim, i, i0, i1, j, j0, j1, dir, n0, n1;
     int ncell_ext[2];

     dim = 2;
     for (i = 0; i < dim; i++) { 
         ncell_ext[i] = ncell[i] + nbdry + nbdry;
     }
//   xl edge without corners 
     if (btype_lower[0] == bdry_transmitted) {
         i0 = nbdry; 
         for (j = nbdry; j < ncell[1]+nbdry; j++) { 
             for (i = 0; i < nbdry; i++) { 
                 for (dir = 0; dir < dim; dir++) { 
                     vel_for_2dcell[j][i][dir] = vel_for_2dcell[j][i0][dir];     
                 } 
             }
         }
     } 
     else { 
     
     }
//   xr edge without corners 
     if (btype_upper[0] == bdry_transmitted) { 
         i0 = ncell[0] + nbdry - 1; 
         for (j = nbdry; j < ncell[1]+nbdry; j++) {
             for (i = 0; i < nbdry; i++) { 
                 i1 = i0 + 1 + i;
                 for (dir = 0; dir < dim; dir++) { 
                     vel_for_2dcell[j][i1][dir] = vel_for_2dcell[j][i0][dir];
                 } 
             }
         }
     }
     else { 

     }
//   yl edge including corners 
     if (btype_lower[1] == bdry_transmitted) {  
         j0 = nbdry; 
         for (j = 0; j < nbdry; j++) { 
             for (i = 0; i < ncell_ext[0]; i++) { 
                 for (dir = 0; dir < dim; dir++) {
                     vel_for_2dcell[j][i][dir] = vel_for_2dcell[j0][i][dir];
                 }
             } 
         } 
     }
     else {

     }
//   yr edge including corners  
     if (btype_upper[1] == bdry_transmitted) {
         j0 = ncell[1] + nbdry - 1; 
         for (j = 0; j < nbdry; j++) {  
             j1 = j0 + 1 + j;
             for (i = 0; i < ncell_ext[0]; i++) {  
                 for (dir = 0; dir < dim; dir++) {
                     vel_for_2dcell[j1][i][dir] = vel_for_2dcell[j0][i][dir];
                 }
             }
         }
     }
     else {

     }
     return;
 } 

void bdry_cell_vel_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
                      double ****vel_for_3dcell) 
{
     int dim, i, i0, i1, j, j0, j1, k, k0, k1, dir, n0, n1;
     int ncell_ext[3];

     dim = 3; 
     for (i = 0; i < dim; i++) { 
         ncell_ext[i] = ncell[i] + nbdry + nbdry;
     }
//   xl face without edeges 
     if (btype_lower[0] == bdry_transmitted) {
         i0 = nbdry; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) { 
             for (j = nbdry; j < ncell[1]+nbdry; j++) { 
                 for (i = 0; i < nbdry; i++) { 
                     for (dir = 0; dir < dim; dir++) { 
                         vel_for_3dcell[k][j][i][dir] = vel_for_3dcell[k][j][i0][dir];     
                     } 
                 }
             }
         }
     } 
     else { 
     
     }
//   xr face without edeges 
     if (btype_upper[0] == bdry_transmitted) { 
         i0 = ncell[0] + nbdry - 1; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) { 
             for (j = nbdry; j < ncell[1]+nbdry; j++) {
                 for (i = 0; i < nbdry; i++) { 
                     i1 = i0 + 1 + i;
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dcell[k][j][i1][dir] = vel_for_3dcell[k][j][i0][dir];
                     }
                 }
             }
         }
     }
     else { 

     }
//   yl face with edges  
     if (btype_lower[1] == bdry_transmitted) {  
         j0 = nbdry; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) {
             for (j = 0; j < nbdry; j++) { 
                 for (i = 0; i < ncell_ext[0]; i++) { 
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dcell[k][j][i][dir] = vel_for_3dcell[k][j0][i][dir];
                     }
                 }
             } 
         } 
     }
     else {

     }
//   yr face with edge 
     if (btype_upper[1] == bdry_transmitted) {
         j0 = ncell[1] + nbdry - 1; 
         for (k = nbdry; k < ncell[2]+nbdry; k++) {
             for (j = 0; j < nbdry; j++) {  
                 j1 = j0 + 1 + j;
                 for (i = 0; i < ncell_ext[0]; i++) {  
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dcell[k][j1][i][dir] = vel_for_3dcell[k][j0][i][dir];
                     }
                 }
             }
         }
     }
     else {

     }
//   zl face including edges and corners
     if (btype_lower[2] == bdry_transmitted) {
         k0 = nbdry;
         for (k = 0; k < nbdry; k++) {
             for (j = 0; j < ncell_ext[1]; j++) {
                 for (i = 0; i < ncell_ext[0]; i++) {
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dcell[k][j][i][dir] = vel_for_3dcell[k0][j][i][dir];
                     }
                 }
             }
         }
     }
     else {

     }
//   zr face including edges and corners 
     if (btype_upper[2] == bdry_transmitted) {
         k0 = ncell[2] + nbdry - 1;
         for (k = 0; k < nbdry; k++) { 
             k1 = k0 + 1 + k;
             for (j = 0; j < ncell_ext[1]; j++) {
                 for (i = 0; i < ncell_ext[0]; i++) {
                     for (dir = 0; dir < dim; dir++) {
                         vel_for_3dcell[k1][j][i][dir] = vel_for_3dcell[k0][j][i][dir];
                     }
                 }
             }
         }
     }
     else {

     }


     return;
 }  

