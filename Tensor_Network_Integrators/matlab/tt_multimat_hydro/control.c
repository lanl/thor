
#ifdef MPI
#include "mpi.h"
#endif

#include "minip.h"
#include "control.h"
#include "advection.h"
#include "mesh.h"
#include "mat.h"
#include "update.h"
// #include "mio.h"
#include "io.h"
#include "bdry.h"
#include "time.h"

double **cs_2dcell = NULL;
double ***cs_3dcell= NULL;

int fix_dt = 1;
int CSDUMP = 1;
int debug_control = 1;
// #define H5DUMP_DEBUG
#ifdef H5DUMP_DEBUG
  #include "hdf5.h"
  #include "xdmf.h"
  #include "h5aux.h"
  hid_t debug_gid;
#endif

void courant_from_cs(int dim, int *ncell, int nbdry, double *dx, double **cs_2dcell, double ***cs_for_3dcell,
				double dt, double *courant_cs);

void courant_from_vel(int dim, int *ncell, int nbdry, double *dx,
					  double ***vel_2dnode, double ****vel_3dnode,
					  double dt, double *courant);

void control(const char *probname,
			int dim, int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
			double *xl_prob, double *xr_prob,
			int nmat, int *matids, int *solid_num_ea_mat, double *gamma_ea_mat,
			int *if_fixed_state_ea_mat,
			double *rho_fixed_state_ea_mat, double *ei_fixed_state_ea_mat, double *pres_fixed_state_ea_mat,
			double dt_initial, double tmax, int ncycle_max, double courant,
			int ncycle_viz_freq, double dt_viz_freq, int *ncycle_final, int mesh_scaling)
{
	char filename[32], tempfile[256];
	int  fileid; 
	int  ncycle, ncycle_to_dump, i, j, k, pass_start, pass, dir, nm_tot, mx, m;
	int  ncell_ext[3], nnode_ext[3], ncell_nmat[4], dims[4];
	long long lsize, lsize2, offset, l, lsize_node;

	double t, dt, ddt,  qvis_coef, dx[3]; 
	double cournmx, courant_cs, courant_div, courant_adv, courant_ener;
	double t_to_dump;

	double *cs1d, **cs2d, *rho1d, **rho2d, *es1d, **es2d;

	double **rho_2dcell,  ***rho_3dcell;
	double **pres_2dcell, ***pres_3dcell;
	double **ei_2dcell,   ***ei_3dcell;
	double **divu_2dcell, ***divu_3dcell; 
	double **qvis_2dcell, ***qvis_3dcell; 
	
	double ***vf_2dmat,   ****vf_3dmat; 
	double ***rho_2dmat,  ****rho_3dmat;
	double ***ei_2dmat,   ****ei_3dmat;
	double ***pres_2dmat, ****pres_3dmat;

	double **force_2dnode, ***force_3dnode; 
	double ***vel_2dnode, ***velav_2dnode;
	double ****vel_3dnode, ****velav_3dnode;

	int    nmixcell_mpoly;
	int    **ijk_in_mixcell;

	double **var_for_2dnode_debug = NULL;  // auxiliary variable for writing velocity components
	double ***var_for_3dnode_debug = NULL; // auxiliary variable for writing velocity components
	double **v2d = NULL, *v1d = NULL;      // auxiliary variables for allocation

	qvis_coef = 0.1;

	fileid = -1;

	lsize = 1.0;
	lsize_node = 1.0;
	for (i = 0; i < dim; i++) {
		dx[i] = (xr_prob[i] - xl_prob[i])/(double)ncell[i];

		ncell_ext[i] = ncell[i] + nbdry + nbdry;
		nnode_ext[i] = ncell[i] + nbdry + nbdry + 1;
		lsize *= ncell_ext[i];
		lsize_node *= nnode_ext[i];
	}
	if (dim == 2) {
		if (!cs_2dcell) {
			cs_2dcell     = (double **) malloc(ncell_ext[1] * sizeof(double *));
			cs_2dcell[0]  = (double  *) malloc(lsize * sizeof(double));
			for (j = 1; j < ncell_ext[1]; j++) {
				cs_2dcell[j]  = cs_2dcell[j-1] + ncell_ext[0];
			}
		}
	}
	else if (dim == 3) {
		lsize2 = ncell_ext[1] * ncell_ext[2];
		if (!cs_3dcell) {
			cs_3dcell  = (double ***) malloc(ncell_ext[2] * sizeof(double **));
			cs2d  = (double **) malloc(lsize2 * sizeof(double *));
			cs1d  = (double  *) malloc(lsize * sizeof(double));

			offset = 0;
			for (k = 0; k < ncell_ext[2]; k++) {
				cs2d[0]  = cs1d  + offset;
				for (j = 1; j < ncell_ext[1]; j++) {
					cs2d[j] = cs2d[j-1] + ncell_ext[0];
				}
				cs_3dcell[k] = cs2d;
				offset += (ncell_ext[0] * ncell_ext[1]);
				cs2d   += ncell_ext[1];
			}
		} 
	}

#ifdef H5DUMP_DEBUG
	const char *fname_h5 = "debug.h5";
	char timestep_group_name[20];
	hid_t fid_h5;
	char dataset_name[64];

	int ncell_ext_inv[3], nnode_ext_inv[3], dims_vel[4];

	if (dim == 3)
	{
		// Allocate 3D debug array var_for_3dnode_debug[z][y][x]
		var_for_3dnode_debug = (double ***)malloc(nnode_ext[2] * sizeof(double **));								 // z
		double **v2d_alloc = (double **)malloc(nnode_ext[1] * nnode_ext[2] * sizeof(double *));			 // y
		double *v1d = (double *)malloc(nnode_ext[0] * nnode_ext[1] * nnode_ext[2] * sizeof(double)); // x

		if (!var_for_3dnode_debug || !v2d_alloc || !v1d)
		{
			fprintf(stderr, "Memory allocation failed for debug velocity arrays.\n");
			exit(1);
		}

		int offset = 0;
		for (int k = 0; k < nnode_ext[2]; k++) // z
		{
			double **v2d = v2d_alloc + k * nnode_ext[1];
			v2d[0] = v1d + offset;
			for (int j = 1; j < nnode_ext[1]; j++) // y
			{
				v2d[j] = v2d[j - 1] + nnode_ext[0]; // x
			}

			var_for_3dnode_debug[k] = v2d;
			offset += nnode_ext[0] * nnode_ext[1];
		}

		// Save inverse dimensions for HDF5 writes
		ncell_ext_inv[0] = ncell_ext[2];
		ncell_ext_inv[1] = ncell_ext[1];
		ncell_ext_inv[2] = ncell_ext[0];

		nnode_ext_inv[0] = nnode_ext[2];
		nnode_ext_inv[1] = nnode_ext[1];
		nnode_ext_inv[2] = nnode_ext[0];

		// dims_vel = [Z, Y, X, 3]
		dims_vel[0] = nnode_ext[2];
		dims_vel[1] = nnode_ext[1];
		dims_vel[2] = nnode_ext[0];
		dims_vel[3] = 3;

		// For writing [Z, Y, X, nmat] data later
		ncell_nmat[0] = ncell_ext[2];
		ncell_nmat[1] = ncell_ext[1];
		ncell_nmat[2] = ncell_ext[0];
		ncell_nmat[3] = nmat;
	}
#endif

	//   pass the data from mesh
	mesh_pass_mesh_data(&rho_2dcell, &ei_2dcell, &pres_2dcell, &divu_2dcell, &qvis_2dcell, 
						&force_2dnode, &vel_2dnode, &velav_2dnode,
						&vf_2dmat,   &rho_2dmat, &ei_2dmat,    &pres_2dmat,   
						&rho_3dcell, &ei_3dcell, &pres_3dcell, &divu_3dcell, &qvis_3dcell, 
						&force_3dnode, &vel_3dnode, &velav_3dnode,
						&vf_3dmat,   &rho_3dmat, &ei_3dmat,    &pres_3dmat);

	t = 0.0;
	ncycle = 0;
	ncycle_to_dump = ncycle_viz_freq;
	t_to_dump      = dt_viz_freq;

	//? Write/Dump data
	#ifdef H5DUMP_DEBUG
		delete_file_if_exists(fname_h5);
		fid_h5 = h5_new_file(fname_h5);
		H5Fclose(fid_h5);
	#endif

	if (CSDUMP)
	{
		xdmf_dump(ncycle, dt);
	}

	char savefilename[256];
	sprintf(savefilename, "%s_mesh_%d", probname, mesh_scaling);
	xdmf_dump_basic(savefilename, ncycle, t);


	//*
	mesh_sound_speed(dim, ncell, nbdry, nmat, solid_num_ea_mat, gamma_ea_mat, cs_2dcell, cs_3dcell);

	courant_from_cs(dim, ncell, nbdry, dx, cs_2dcell, cs_3dcell, dt_initial, &courant_cs);
	
	courant_from_vel(dim, ncell, nbdry, dx, vel_2dnode, vel_3dnode, dt, &courant_adv);
	courant_cs = MAX(courant_cs, courant_adv);

	if (courant_cs > courant)
	{
		dt = courant * dt_initial / courant_cs;
		printf("WARNING: dt_initial too big, used dt = %e\n", dt);
	}
	else {
		dt = dt_initial;
	}
	pass_start = 0;
	pass = pass_start;

//     tmax = 5.0 * dt;
	// printf("tmax = %f\n", tmax);
	double dt0;
	dt0 = dt;

	clock_t loop_start = clock();

	while ((t < tmax) && (ncycle < ncycle_max))
	{
#ifdef H5DUMP_DEBUG
		// open an HDF5 group to write the current cycle
		sprintf(timestep_group_name, "Step#%d", ncycle);
		fid_h5 = h5_open_existing_rdwr(fname_h5);
		debug_gid = h5_open_group(fid_h5, timestep_group_name);
		if (ncycle == 0)
		{
			h5_write_4d(debug_gid, "init:vel_3dnode", vel_3dnode, dims_vel);
			h5_write_4d(debug_gid, "init:val_3dnode", velav_3dnode, dims_vel);
		}
#endif
		cournmx = 0.0;
		if (fix_dt)
		{
			dt = dt0; //! set fixed time step
			printf("dt = %f -- ", dt);
		}

		//         Lagrangian phase

		compute_divu(dim, ncell, nbdry, dx, vel_2dnode, vel_3dnode, divu_2dcell, divu_3dcell);
#ifdef H5DUMP_DEBUG
		h5_write_3d(debug_gid, "divu", divu_3dcell, ncell_ext_inv);
#endif
		compute_qvis(dim, ncell, nbdry, dx[0],
			rho_2dcell, cs_2dcell, divu_2dcell, vel_2dnode,
			rho_3dcell, cs_3dcell, divu_3dcell, vel_3dnode,
			qvis_2dcell, qvis_3dcell);

#ifdef H5DUMP_DEBUG
		h5_write_3d(debug_gid, "qvis", qvis_3dcell, ncell_ext_inv);
#endif

		for (dir = 0; dir < dim; dir++)
		{
			compute_force(dim, ncell, nbdry, dx, dir,
						  pres_2dcell, pres_3dcell,
						  rho_2dcell, rho_3dcell,
						  qvis_2dcell, qvis_3dcell,
						  force_2dnode, force_3dnode);

			update_vel_comp(dim, ncell, nbdry, dt, dir,
							force_2dnode, force_3dnode,
							vel_2dnode, velav_2dnode, vel_3dnode, velav_3dnode);
		}

		update_energy(nmat, dim, ncell, nbdry, dt,
					  rho_2dcell, divu_2dcell, vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat,
					  rho_3dcell, divu_3dcell, vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat,
					  &courant_ener);

		cournmx = MAX(cournmx, courant_ener);

		update_density(nmat, dim, ncell, nbdry, dt,
					   if_fixed_state_ea_mat, rho_fixed_state_ea_mat,
					   ei_fixed_state_ea_mat, pres_fixed_state_ea_mat,
					   rho_2dcell, divu_2dcell, vf_2dmat, rho_2dmat,
					   rho_3dcell, divu_3dcell, vf_3dmat, rho_3dmat,
					   &courant_div);

		cournmx = MAX(cournmx, courant_div);

		//     change ei_2dmat or ei_3dmat from specific internal energy density
		//     to internal energy density,
		//     Also, update pressure
#ifdef H5DUMP_DEBUG
	h5_write_4d(debug_gid, "before_update_pressure:pres_3dmat", pres_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "before_update_pressure:ei_3dmat", ei_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "before_update_pressure:rho_3dmat", rho_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "before_update_pressure:vf_3dmat", vf_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "before_update_pressure:vel_3dnode", vel_3dnode, dims_vel);
	h5_write_4d(debug_gid, "before_update_pressure:velav_3dnode", velav_3dnode, dims_vel);
#endif
		update_pressure(nmat, solid_num_ea_mat, gamma_ea_mat,
						dim, ncell, nbdry,
						if_fixed_state_ea_mat, rho_fixed_state_ea_mat,
						ei_fixed_state_ea_mat, pres_fixed_state_ea_mat,
						ei_2dcell, pres_2dcell,
						vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat,
						ei_3dcell, pres_3dcell,
						vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat);

#ifdef H5DUMP_DEBUG
		h5_write_4d(debug_gid, "after_update_pressure:pres_3dmat", pres_3dmat, ncell_nmat);
#endif

		if (dim == 2)
		{
			bdry_cell_2d(nmat, ncell, nbdry, btype_lower, btype_upper,
						 vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat,
						 rho_2dcell, ei_2dcell, pres_2dcell);

			bdry_node_2d(ncell, nbdry, btype_lower, btype_upper, vel_2dnode);
			bdry_node_2d(ncell, nbdry, btype_lower, btype_upper, velav_2dnode);
		}
		else if (dim == 3)
		{
			bdry_cell_3d(nmat, ncell, nbdry, btype_lower, btype_upper,
						 vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat,
						 rho_3dcell, ei_3dcell, pres_3dcell);

			bdry_node_3d(ncell, nbdry, btype_lower, btype_upper, vel_3dnode);
			bdry_node_3d(ncell, nbdry, btype_lower, btype_upper, velav_3dnode);
		}

#ifdef H5DUMP_DEBUG
	h5_write_4d(debug_gid, "Lagrange:pres_3dmat", pres_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "Lagrange:ei_3dmat", ei_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "Lagrange:rho_3dmat", rho_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "Lagrange:vf_3dmat", vf_3dmat, ncell_nmat);
	h5_write_4d(debug_gid, "Lagrange:vel_3dnode", vel_3dnode, dims_vel);
	h5_write_4d(debug_gid, "Lagrange:velav_3dnode", velav_3dnode, dims_vel);
#endif

		//!         advection phase

		for (dir = 0; dir < dim; dir++)
		{

			advection(fileid, dim, ncell, nbdry, xl_prob, dx,
					  nmat, solid_num_ea_mat, gamma_ea_mat,
					  if_fixed_state_ea_mat,
					  rho_fixed_state_ea_mat, ei_fixed_state_ea_mat, pres_fixed_state_ea_mat,
					  ncycle, pass, dt,
					  btype_lower, btype_upper,
					  vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat,
					  rho_2dcell, ei_2dcell, pres_2dcell, vel_2dnode, velav_2dnode,
					  vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat,
					  rho_3dcell, ei_3dcell, pres_3dcell, vel_3dnode, velav_3dnode,
					  &courant_adv);

			cournmx = MAX(cournmx, courant_adv);

			pass = (pass + 1) % dim;
		}
		pass_start = (pass_start + 1) % dim;
		pass = pass_start;

		t += dt;
		ncycle++;

		// determine the next dt

		mesh_sound_speed(dim, ncell, nbdry, nmat, solid_num_ea_mat, gamma_ea_mat, cs_2dcell, cs_3dcell);

		courant_from_cs(dim, ncell, nbdry, dx, cs_2dcell, cs_3dcell, dt, &courant_cs);
		cournmx = MAX(cournmx, courant_cs);
		courant_from_vel(dim, ncell, nbdry, dx, vel_2dnode, vel_3dnode, dt, &courant_adv);
		cournmx = MAX(cournmx, courant_adv);

		printf(" ncycle = %d,  t = %e, dt = %e, courant = %e\n", ncycle, t, dt, cournmx);

		ncycle_to_dump--;
		t_to_dump -= dt;
		if ((ncycle_to_dump <= 0) || (t_to_dump <= 0) || (t >= tmax))
		{
			xdmf_dump(ncycle, t);
			ncycle_to_dump = ncycle_viz_freq;
			t_to_dump = dt_viz_freq;
		}
		//         determine the next dt

		if (cournmx > courant)
		{
			dt *= (courant / cournmx);
		}
		else if (cournmx < courant)
		{
			ddt = (courant / cournmx - 1.0) * dt;
			dt += (0.1 * ddt);
		}

		//! round dt
		// dt = ceil(dt * 100000) / 100000;
		dt = MIN(dt, tmax - t);
#ifdef H5DUMP_DEBUG
		H5Gclose(debug_gid);
		H5Fclose(fid_h5);
#endif
		if (ncycle_final != NULL)
			*ncycle_final = ncycle;
	}

	clock_t loop_end = clock();
  double loop_elapsed = (double)(loop_end - loop_start) / CLOCKS_PER_SEC;

  printf("=== Timing Summary ===\n");
  printf("Elapsed time in loop = %.6f seconds\n", loop_elapsed);
  printf("Total cycles = %d\n", ncycle);
  if (ncycle > 0)
  {
	printf("Average time per cycle = %.6f seconds\n", loop_elapsed / ncycle);
  }
  //? write timing data out
  char csv_filename[128];
  sprintf(csv_filename, "%s_loop_timing.csv", probname);
  FILE *f = fopen(csv_filename, "a");
  if (f != NULL)
  {
	// write header only if file is empty
	static int header_written = 0;
	if (!header_written)
	{
	  fprintf(f, "mesh_scaling,cycles,elapsed_sec,sec_per_cycle\n");
	  header_written = 1;
	}

	fprintf(f, "%d,%d,%.6f,%.6f\n", mesh_scaling, ncycle, loop_elapsed, loop_elapsed / ncycle);
	fclose(f);
  }
  //? export the data
  // char savefilename[256];
  // sprintf(savefilename, "%s_mesh%d", probname, mesh_scaling);
  xdmf_dump_basic(savefilename, ncycle, t);

  if (dim == 2)
  {
	free(cs_2dcell[0]);
	free(cs_2dcell);
  }
  else if (dim == 3)
  {
	free(cs_3dcell[0][0]);
	free(cs_3dcell[0]);
	free(cs_3dcell);
  }

return;
}

void courant_from_cs(int dim, int *ncell, int nbdry, double *dx, double **cs_2dcell, double ***cs_3dcell,
					 double dt, double *courant_cs)
{

	int  i, j, k;
	long long lsize, c;
	double dxmin, cournt, *cs1d;

	*courant_cs = 0.0;

	lsize = 1;
	dxmin = 1.0e+10;
	for (i = 0; i < dim; i++) {
		lsize *= (ncell[i] + nbdry + nbdry); 
		if (dxmin > dx[i]) dxmin = dx[i];
	}
	if (dim == 2) {
		cs1d = cs_2dcell[0];
	}
	else if (dim == 3) {
		cs1d = cs_3dcell[0][0];
	}
	for (c = 0; c < lsize; c++) {
		cournt = dt * cs1d[c]/dxmin;
		*courant_cs = MAX(*courant_cs, cournt);
	}
	return;
}

void courant_from_vel(int dim, int *ncell, int nbdry, double *dx,   
					  double ***vel_2dnode, double ****vel_3dnode,
					  double dt, double *courant)
{     
	int  i, j, k, dir;
	int  ncell_ext[3];  
	double c; 

	for (i = 0; i < dim; i++) { 
		ncell_ext[i] = ncell[i] + nbdry + nbdry;
	}
	*courant = 0.0; 
	if (dim == 2) {
		for (j = 0; j <= ncell_ext[1]; j++) {
			for (i = 0; i <= ncell_ext[0]; i++) {
				for (dir = 0; dir < dim; dir++) {
					c = fabs(vel_2dnode[j][i][dir]) * dt / dx[dir];
					*courant = MAX(c, *courant);
				}
			}
		}
	}
	else if (dim == 3) {
		for (k = 0; k <= ncell_ext[2]; k++) {
			for (j = 0; j <= ncell_ext[1]; j++) {
				for (i = 0; i <= ncell_ext[0]; i++) {
					for (dir = 0; dir < dim; dir++) {
						c = fabs(vel_3dnode[k][j][i][dir]) * dt / dx[dir];
						*courant = MAX(c, *courant);
					}
				}
			}
		}
	}
	return;
} 


