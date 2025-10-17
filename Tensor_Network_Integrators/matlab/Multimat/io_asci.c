
#include <assert.h> 
#include <stdio.h>
#include "io_asci.h" 

void dump2d_asci(const char *probname, double t,  int ncycle,  
		 double *xl_prob, double *xr_prob,
 		 int *ncell, int nbdry, int nmat_prob, 
                 double **rho_for_2dcell,  double **ei_for_2dcell, double **pres_for_2dcell, 
                 int **nmat_for_2dcell, int **mixcell_for_2dcell, int **matid_for_2dcell, 
                 double ***vel_for_2dnode,
                 int nmixcell, int **ijk_in_mixcell, int *nmat_in_mixcell, 
                 int **matids_in_mixcell, double **vf_in_mixcell, double **rho_in_mixcell, 
                 double **ei_in_mixcell, double **pres_in_mixcell)
{ 
     char filename[128]; 
     int dim, i, i0, j, j0, nm, mx, m; 
     FILE *fp;

     dim = 2;
   
     sprintf(filename, "%s_%08d_asci", probname, ncycle);

     fp = fopen(filename, "w");
     
     fprintf(fp, "t = %e\n", t);
     fprintf(fp, "ncycle = %d\n", ncycle);

     fprintf(fp, "xl_prob =");
     for (i = 0; i < dim; i++) { 
	 fprintf(fp, " %e", xl_prob[i]);
     }
     fprintf(fp,"\n");
     fprintf(fp, "xr_prob =");
     for (i = 0; i < dim; i++) {
         fprintf(fp, " %e", xr_prob[i]);
     }
     fprintf(fp,"\n");

     fprintf(fp, "ncell = %d, %d\n", ncell[0], ncell[1]);
     fprintf(fp, "nmat_prob = %d\n", nmat_prob);

     fprintf(fp, "     j      i       vx            vy           rho            ei          pres     nmat id1    vf1          rho1          ei1          pres1    id2      vf2          rho2          ei2          pres2   ....  \n");

     for (j0 = 0; j0 < ncell[1]; j0++) {
	 j = j0 + nbdry;
	 for (i0 = 0; i0 < ncell[0]; i0++) { 
             i = i0 + nbdry;
             fprintf(fp, "%6d %6d %13e %13e %13e %13e %13e %2d ", j0, i0,
                     vel_for_2dnode[j][i][0], 
		     vel_for_2dnode[j][i][1], 
	             rho_for_2dcell[j][i],
		      ei_for_2dcell[j][i], 
                    pres_for_2dcell[j][i], 
		    nmat_for_2dcell[j][i]);
	     if (nmat_for_2dcell[j][i] == 1) { 
		 fprintf(fp, "%2d\n", matid_for_2dcell[j][i]);
             }
	     else { 
                 nm = nmat_for_2dcell[j][i];
		 mx = mixcell_for_2dcell[j][i];
		 assert(mx >= 0);
                 assert(nm == nmat_in_mixcell[mx]);
                 for (m = 0; m < nm; m++) { 
		     fprintf(fp, "%2d %13e %13e %13e %13e", 
		             matids_in_mixcell[mx][m], 
			     vf_in_mixcell[mx][m],
			     rho_in_mixcell[mx][m],
			     ei_in_mixcell[mx][m], 
			     pres_in_mixcell[mx][m]); 
		 }
		 fprintf(fp, "\n");
             } 
         }
//       wrte the node velocity at i0 = ncell[0], j0. 
         i0 = ncell[0];
         i  = i0 + nbdry;
         fprintf(fp, "%6d %6d %13e %13e\n", j0, i0, 
                 vel_for_2dnode[j][i][0],
                 vel_for_2dnode[j][i][1]);
     }
//   write the velocity at j0 = ncell[1]
     j0 = ncell[1];
     j = j0 + nbdry;
     for (i0 = 0; i0 <= ncell[0]; i0++) {
         i = i0 + nbdry;
         fprintf(fp, "%6d %6d %13e %13e\n", j0, i0,
                     vel_for_2dnode[j][i][0],  
                     vel_for_2dnode[j][i][1]);
     }
     fclose(fp); 

     return;
 }

void dump3d_asci(const char *probname, double t, int ncycle, 
		 double *xl_prob, double *xr_prob,
		 int *ncell, int nbdry, int nmat_prob,
                 double ***rho_for_3dcell,  double ***ei_for_3dcell, double ***pres_for_3dcell,
                 int ***nmat_for_3dcell, int ***mixcell_for_3dcell, int ***matid_for_3dcell, 
                 double ****vel_for_3dnode,
                 int nmixcell, int **ijk_in_mixcell, int *nmat_in_mixcell, 
                 int **matids_in_mixcell, double **vf_in_mixcell, double **rho_in_mixcell, 
		 double **ei_in_mixcell, double **pres_in_mixcell)
{ 
     char filename[128];
     int dim, i, i0, j, j0, k, k0, nm, mx, m;
     FILE *fp;
     
     dim = 2; 

     sprintf(filename, "%s_%08d_asci", probname, ncycle);
     fp = fopen(filename, "w");
     
     fprintf(fp, "t = %e\n", t); 
     fprintf(fp, "ncycle = %d\n", ncycle); 

     fprintf(fp, "xl_prob =");
     for (i = 0; i < dim; i++) {
         fprintf(fp, " %e", xl_prob[i]);
     }
     fprintf(fp,"\n");
     fprintf(fp, "xr_prob =");
     for (i = 0; i < dim; i++) {
         fprintf(fp, " %e", xr_prob[i]);
     }
     fprintf(fp,"\n");

     fprintf(fp, "ncell = %d %d %d\n", ncell[0], ncell[1], ncell[2]);
     fprintf(fp, "nmat_prob = %d\n", nmat_prob);

     fprintf(fp, "     k      j      i        vx            vy            vz           rho            ei           pres     nmat id1    vf1          rho1          ei1          pres1    id2      vf2          rho2          ei2          pres2   ....  \n");

     for (k0 = 0; k0 < ncell[2]; k0++) { 
	 k = k0 + nbdry; 
         for (j0 = 0; j0 < ncell[1]; j0++) {
             j = j0 + nbdry;
             for (i0 = 0; i0 < ncell[0]; i0++) {
                 i = i0 + nbdry;
                 fprintf(fp, "%6d %6d %6d   %e %13e %13e %13e %13e %13e   %d ", k0, j0, i0,
                         vel_for_3dnode[k][j][i][0],  
                         vel_for_3dnode[k][j][i][1], 
			 vel_for_3dnode[k][j][i][2],
                         rho_for_3dcell[k][j][i],
                          ei_for_3dcell[k][j][i],  
                        pres_for_3dcell[k][j][i],  
                        nmat_for_3dcell[k][j][i]);

                 if (nmat_for_3dcell[k][j][i] == 1) {  
                     fprintf(fp, "%2d\n", matid_for_3dcell[k][j][i]);
                 }
                 else {
                     nm = nmat_for_3dcell[k][j][i];
                     mx = mixcell_for_3dcell[k][j][i];
		     assert(mx >= 0); 
                     assert(nm == nmat_in_mixcell[mx]);
                     for (m = 0; m < nm; m++) {
                         fprintf(fp, "%2d %13e %13e %13e %13e",
                                 matids_in_mixcell[mx][m], 
				 vf_in_mixcell[mx][m],
                                 rho_in_mixcell[mx][m],
                                 ei_in_mixcell[mx][m],
                                 pres_in_mixcell[mx][m]);
                     }
                     fprintf(fp, "\n");
                 }
             }
//           write the node velocity at i0 = ncell[0], j0, and k0 
             i0 = ncell[0];
             i = i0 + nbdry;
             fprintf(fp, "%6d %6d %6d   %e %13e %13e\n", k0, j0, i0, 
                     vel_for_3dnode[k][j][i][0],
                     vel_for_3dnode[k][j][i][1],
                     vel_for_3dnode[k][j][i][2]);
         }
//       write the node velocity at j = ncell[1] and k0 
         j0 = ncell[1];
         j  = j0 + nbdry;
         for (i0 = 0; i0 <= ncell[0]; i0++) {
             i = i0 + nbdry;
             fprintf(fp, "%6d %6d %6d   %e %13e %13e\n", k0, j0, i0,
                     vel_for_3dnode[k][j][i][0],
                     vel_for_3dnode[k][j][i][1],
                     vel_for_3dnode[k][j][i][2]);
         }
     } 
//   write the node velocity at k0 = ncell[2] 
     k0 = ncell[2];
     k  = k0 + nbdry; 
     for (j0 = 0; j0 <= ncell[1]; j0++) {
         j = j0 + nbdry;
         for (i0 = 0; i0 <= ncell[0]; i0++) {
             i = i0 + nbdry;
             fprintf(fp, "%6d %6d %6d  %13e %13e %13e\n", k0, j0, i0,
                         vel_for_3dnode[k][j][i][0],
                         vel_for_3dnode[k][j][i][1],
	                 vel_for_3dnode[k][j][i][2]); 
         }
     } 
     fclose(fp);

     return;
  }

 
 
