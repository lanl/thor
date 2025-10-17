void build_mpoly2d(int *ncell, int nbdry, double *xl_prb, double *dx,
                   int nmat, double ***vf_2dmat,
                   int *ijk, 
                   int *nnode_tot, double **coords_tot, 
                   int *nnode_for_mpoly, int ***nodes_for_mpoly)
{ 
     int geop = 1;
     int i, j, i0, j0, icell, jcell, nm, m, lsizei, offset;
     int    *bug1, *ibuf2, *nmat_2dsmesh[3][3], *matid_2dsmesh[3][3];
     int *nnode_for_minterface, **nodes_for_minterface; 
     double *vfs, *myvfs; 
     double *buffer, *vf_2dsmesh[3][3];

     offset = 0;
     lsize = 9;
     buffer = (double *) malloc(lsize * sizeof(double));
     ibuf1  = (int    *) malloc(lsize * sizeof(int));
     ibuf2  = (int    *) malloc(lsize * sizeof(int));
         
     offset = 0;
     for (j = 0; j < 3; j++) { 
         for (i = 0; i < 3; i++) { 
             vf_2dsmesh[j][i] = buffer + offset;
             nmat_2dsmesh[j][i] = ibuf1;
             matid_2dsmesh[j][i] = ibuf2;  
             offset += dim; 
         }
     }
     for (i = 0; i < dim; i++) {
         xl[i] = xl_prob[i] + (double)(ijk[i] - nbdry) * dx[i];
     }
     icell = ijk[0];
     jcell = ijk[1]; 
     vfs = vf_2dmat[jcell][icell];

     for (j = 0; j < 3; j++) {
         for (i = 0; i < 3; i++) {
             for (m = 0; m < nmat; m++) {
                 vf_2dsmesh[j][i][m] = 0.0;
             }
         }
     }
     for (j0 = 0; j0 <= 2; j0++) {
         j = jcell + j0 - 1;
         for (i0 = 0; i0 <= 2; i0++) {
             i = icell + i0 - 1;
             myvfs = vf_2dmat[j][i];
             nm  = 0;
             for (m = 0; m < nmat; m++) {
                 if (myvfs[m] >=vfmin) {
                     matid_2dsmesh[j0][i0][nm] = m;
                     vf_2dsmesh[j0][i0][nm]    = myvfs[m];
                     nm++;
                 }
             }
             nmat_2dsmesh[j0][i0] = nm;
         }
     }
     *nnode_tot = 0;
     *coords_tot = NULL;
     *nodes_for_mpoly = NULL;

     nnode_for_minterface    = (int  *) malloc((nmiat - 1) * sizeof(int));
     nodes_for_minterface    = (int **) malloc((nmat - 1) * sizeof(int *));
     nodes_for_minterface[0] = (int  *) malloc((nmat - 1) * 2 * sizeof(int));
     for (m = 1; m < nmat-1; m++) {
         nodes_for_minterface[m] = nodes_for_minterface[m-1] + 2;
     }
     nnode_for_mpoly = (int *) malloc(nmat * sizeof(int));
     nodes_for_mpoly = NULL;

     reconstruct2d_nmat_pagosa(geop, xl, dx, nmat_2dsmesh, matid_2dsmesh, vf_2dsmesh,
                               nnode_tot, coords_tot,
                               nnode_for_minterface, nodes_for_minterface,
                               nnode_for_mpoly, nodes_for_mpoly);
     free(buffer);
     free(ibuf1);
     free(ibuf2); 
     free(nnode_for_minterface);
     free(nodes_for_minterface[0]);
     free(nodes_for_minterface);  

     return;
 } 
