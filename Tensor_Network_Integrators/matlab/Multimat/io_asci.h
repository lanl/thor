#ifndef _IO_ASCI_
#define _IO_ASCI_

#ifdef __cplusplus
extern "C" {
#endif

void dump2d_asci(const char *probname, double t,  int ncycle,
		 double *xl_prob, double *xr_prob,
		 int *ncell, int nbdry, int nmat_prob,
                 double **rho_for_2dcell,  double **ei_for_2dcell, double **pres_for_2dcell,
                 int **nmat_for_2dcell, int **mixcell_for_2dcell, int **matid_for_2dcell,
                 double ***vel_for_2dnode,
                 int nmixcell, int **ijk_in_mixcell, int *nmat_in_mixcell,
                 int **matids_in_mixcell, double **vf_in_mixcell, double **rho_in_mixcell, 
		 double **ei_in_mixcell, double **pres_in_mixcell);

void dump3d_asci(const char *probname, double t, int ncycle,  
		 double *xl_prob, double *xr_prob,
		 int *ncell, int nbdry, int nmat_prob,
                 double ***rho_for_3dcell, double ***ei_for_3dcell, double ***pres_for_3dcell,
                 int ***nmat_for_3dcell, int ***mixcell_for_3dcell, int ***matid_for_3dcell, 
                 double ****vel_for_3dnode,
                 int nmixcell, int **ijk_in_mixcell, int *nmat_in_mixcell, 
                 int **matids_in_mixcell, double **vf_in_mixcell, double **rho_in_mixcell, 
                 double **ei_in_mixcell, double **pres_in_mixcell);

#ifdef __cplusplus
}
#endif
#endif
