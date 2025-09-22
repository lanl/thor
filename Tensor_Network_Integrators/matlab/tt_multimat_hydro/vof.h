/* ======================= copyright begin ========================  */
/* Copyright (C) 2010 Los Alamos National Laboratory.                */
/* All rights Reserved.                                              */
/* Export Controlled Information                                     */
/* ======================== copyright end =========================  */

#ifndef _VOF_
#define _VOF_

#ifdef __cplusplus
extern "C" {
#endif

void vof_init(int dim, double *xl_prob, double *xr_prob, int *ncell, int nbdry, int nmat,
              double t, int ncycle, double *vf1d, double *rho1d, double *pres1d, double *ei1d);

//  This function build material interfaces im mixed cells.

//  dim:     input dimensionality of the probem.
//  xl_prob: input, the left bounds of the simulation domain, excluding the ghost cells. 
//           There are dim values, xl_prob[0:dim).  
//  xr_prob: input, the right bounds of the simulation domain, excluding the ghost cells.
//           There are dim values in it, xr_prob[0:dim).
//  ncell:   the number of cells in each dimension, excluding the ghost cells.
//           There are dim values in it, ncell[0:dim). 
//  nbdry:   input, the number of the layer of ghost cells.
//  nmat:    input, the number of materials in the problem.
//  t:       input, simulation time.
//  ncycle:  input, the current time step.
//  vf1d:    input, the starting address of valume factions of cells. Here is memory layout
//           in C/C++: for 2D vf1d = vf[0][0], here 
//                            vf[0:NY+NBDRY+NBDRY)[0:NX+NBDRY+NBDRY)[0:NMAT) for 2D 
//                     for 3D vf1d = vf[0][0][0], here
//                            vf[0:NZ+NBDRY+NBDRY)[0:NY+NBDRY+NBDRY)[0:NX+NBDRY+NBDRY)[0:NMAT) for 3D   
//           in Fortran: vf1d = vf, here
//                       vf(NMAT, NX+NBDRY+NBDRY, NY+NBDRY+NBDRY) for 2D 
//                       vf(NMAT, NX+NBDRY+NBDRY, NY+NBDRY+NBDRY, NZ+NBDRY+NBDRY) 
// 
//                       or
//                       vf(NMAT, -NBDRY+1:NX+NBDRY, -NBDRY+1:NY+NBDRY, -NBDRY+1:NZ+NBDRY) 

void advected_vol(int dim, int dir,
                  int edge_index, int *cell_index, double adv_distance,
                  int nmat, double *vols_advected);

//  This function calculates advected volume of each material. 

// dim: input, the dimensionalilty of the problem.
// dir: input, 1 or 2 or 3 for the pass in dimensional pass. 
// edge_index: input, 1-based edge number in the "dir" direction, excluding ghost cell. 
// cell_index: input, dim values, cell_index[0:dim), to indicate the cells to be considere, excluding
//             ghost cells. For example, if dir = 1, edge_index = 2, cell_index = (0, 3, 4),
//             we are considering the advection between the cell (1,3,4) and the cell (2,3,4) 
// adv_distance: input, the distance advected. It could be negative. 
//               For the example,  dir = 1, edge_index = 2, cell_index = (0, 3, 4), a positive
//               adv_distance indicates the case that materials move from cell (1,3,4) to the cell (2,3,4);
//               a negative value of adv_distance is for the case materials move from the cell (2,3,4)
//               to the cell (1,3,4).
// nmat: input, the number of the materials in the problem.
// vols_advected: output of nmat values. The volumes advected of each material.      
//  
#ifdef __cplusplus
}
#endif
#endif

