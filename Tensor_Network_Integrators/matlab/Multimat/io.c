#ifdef MPI
#include <mpi.h>
#endif

#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
// #include "mio.h"
#include "io.h"
#include <string.h>

static double vfmin = 1.0e-08;

#include <stdio.h>
#include <stdlib.h>

void write_2Dllarray_to_ascii(const char *filename, long long **int_2dcell, int *ncell_ext)
{
    if (int_2dcell == NULL) return;

	char folder[256] = "../output/"; // Ensure folder array has enough space
	strcat(folder, filename);				 // Concatenate folder path and filename

	// Open the file for writing
	FILE *file = fopen(folder, "w");
	if (file == NULL)
	{
		perror("Error opening file");
		exit(EXIT_FAILURE);
	}

	// Write the data to the file
	for (int i = 0; i < ncell_ext[1]; i++)
	{
		for (int j = 0; j < ncell_ext[0]; j++)
		{
			fprintf(file, "%lld ", int_2dcell[i][j]); // Write each element as an integer
		}
		fprintf(file, "\n"); // Newline after each row
	}

	// Close the file
	fclose(file);
	printf("Successfully written data to %s\n", folder);
}

// Function to write pres_2dmat, ei_2dmat, rho_2dmat, and vf_2dmat
void write_2Dmat_variables(const char *suffix, double ***pres_2dmat, double ***ei_2dmat, double ***rho_2dmat, double ***vf_2dmat, int ny, int nx, int nz, int ncycle)
{
	char tempfile[256];
    if (pres_2dmat == NULL) return;

	// Write pres_2dmat
	sprintf(tempfile, "pres_2dmat_%d_%s.txt", ncycle, suffix);
	write_3Darray_as_matrix(tempfile, pres_2dmat, ny, nx, nz);

	// Write ei_2dmat
	sprintf(tempfile, "ei_2dmat_%d_%s.txt", ncycle, suffix);
	write_3Darray_as_matrix(tempfile, ei_2dmat, ny, nx, nz);

	// Write rho_2dmat
	sprintf(tempfile, "rho_2dmat_%d_%s.txt", ncycle, suffix);
	write_3Darray_as_matrix(tempfile, rho_2dmat, ny, nx, nz);

	// Write vf_2dmat
	sprintf(tempfile, "vf_2dmat_%d_%s.txt", ncycle, suffix);
	write_3Darray_as_matrix(tempfile, vf_2dmat, ny, nx, nz);
}

// Function to write pres_2dmat, ei_2dmat, rho_2dmat, and vf_2dmat
void write_2Dcell_variables(const char *suffix, double **pres_2dcell, double **ei_2dcell, double **rho_2dcell, int ny, int nx, int nz, int ncycle)
{
    if (pres_2dcell == NULL) return;
	char tempfile[256];
	int tempsize[2];
	tempsize[0] = nx; tempsize[1] = ny;
	// Write pres_2dmat
	sprintf(tempfile, "pres_2dcell_%d_%s.txt", ncycle, suffix);
	write_2Darray_to_ascii(tempfile, pres_2dcell, tempsize);

	// Write ei_2dmat
	sprintf(tempfile, "ei_2dcell_%d_%s.txt", ncycle, suffix);
	write_2Darray_to_ascii(tempfile, ei_2dcell, tempsize);

	// Write rho_2dmat
	sprintf(tempfile, "rho_2dcell_%d_%s.txt", ncycle, suffix);
	write_2Darray_to_ascii(tempfile, rho_2dcell, tempsize);
}

void write_3Darray_as_matrix(const char *filename, double ***array, int n, int m, int l)
{
	char folder[256] = "../output/"; // Ensure folder array has enough space
	strcat(folder, filename);				 // Concatenate folder path and filename
	FILE *file = fopen(folder, "w");
	if (file == NULL)
	{
		perror("Error opening file");
		exit(EXIT_FAILURE);
		}

	// Iterate over the first dimension
	for (int i = 0; i < n; i++)
	{
			for (int j = 0; j < m; j++)
			{
					for (int k = 0; k < l; k++)
					{
						fprintf(file, "%.15f ", array[i][j][k]);
					}
			}
			fprintf(file, "\n"); // New row after each i
	}

	fclose(file);
	printf("Successfully written data to %s\n", folder);
}
void write_3Darray_to_ascii(const char *filename, double ***array, int rows, int cols, int depth)
{
	char folder[256] = "../output/"; // Ensure folder array has enough space
	strcat(folder, filename);        // Concatenate folder path and filename
	// Open the file for writing
	FILE *file = fopen(folder, "w");
	if (file == NULL)
	{
			perror("Error opening file");
			exit(EXIT_FAILURE);
	}

	// Write the 3D array to the file
	for (int k = 0; k < depth; k++)
	{
		fprintf(file, "Slice %d:\n", k + 1); // Indicate the current slice
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
					fprintf(file, "%.15f ", array[k][i][j]); // Write each element
			}
			fprintf(file, "\n"); // Newline after each row
		}
		fprintf(file, "\n"); // Separate slices with a blank line
	}

	// Close the file
	fclose(file);
	printf("Successfully written data to %s\n", folder);
}

void write_2Darray_to_ascii(const char *filename, double **rho_2dcell, int *ncell_ext)
{
        if (rho_2dcell == NULL) return;
		char folder[256] = "../output/"; // Ensure folder array has enough space
		strcat(folder, filename);        // Concatenate folder path and filename

		// Open the file for writing
		FILE *file = fopen(folder, "w");
		if (file == NULL)
		{
				perror("Error opening file");
				exit(EXIT_FAILURE);
		}

		// Write the data to the file
		for (int i = 0; i < ncell_ext[1]; i++)
		{
				for (int j = 0; j < ncell_ext[0]; j++)
				{
						fprintf(file, "%.15f ", rho_2dcell[i][j]); // Write each element with 15 decimal places
				}
				fprintf(file, "\n"); // Newline after each row
		}

		// Close the file
		fclose(file);
		printf("Successfully written data to %s\n", folder);
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void write_2Dintarray_to_ascii(const char *filename, int **int_2dcell, int *ncell_ext)
{
    if (int_2dcell == NULL) return;
	char folder[256] = "../output/"; // Ensure folder array has enough space
	strcat(folder, filename);				 // Concatenate folder path and filename

	// Open the file for writing
	FILE *file = fopen(folder, "w");
	if (file == NULL)
	{
		perror("Error opening file");
		exit(EXIT_FAILURE);
	}

	// Write the data to the file
	for (int i = 0; i < ncell_ext[1]; i++)
	{
		for (int j = 0; j < ncell_ext[0]; j++)
		{
			fprintf(file, "%d ", int_2dcell[i][j]); // Write each element as an integer
		}
		fprintf(file, "\n"); // Newline after each row
	}

	// Close the file
	fclose(file);
	printf("Successfully written data to %s\n", folder);
}

void write_array_to_bin(const char *filename, void *array, size_t element_size, int size)
{
		FILE *file = fopen(filename, "wb");
		if (file == NULL)
		{
				perror("Error opening file");
				return;
		}

		// // Write the dimensions first
		// fwrite(&rows, sizeof(int), 1, file);
		// fwrite(&cols, sizeof(int), 1, file);

		// Write the 2D array data
		fwrite(array, element_size, size, file);

		fclose(file);
		printf("Successfully written data to %s\n", filename);
}

void write_mats_2d(int fileid, int nmat, int *ncell, int nbdry,
									 long long nelem_prob, long long nnode_prob, double *xy_prob[2],
									 int *nodelist_for_elem_prob,
									 double ***vf_2dmat, double ***rho_2dmat, double ***ei_2dmat, double ***pres_2dmat); 

void write_mats_3d(int fileid, int nmat, int *ncell, int nbdry,
									 long long nelem_prob, long long nnode_prob, double *xyz_prob[3],
									 int *nodelist_for_elem_prob,
									 double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat, double ****pres_3dmat); 


void write_line_cell_2d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												int nmat, double ***vf_2dmat,
												double **rho_2dcell, double **ei_2dcell,
												double **pres_2dcell, double **divu_2dcell);

void write_line_cell_3d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												int nmat, double ****vf_3dmat,
												double ***rho_3dcell, double ***ei_3dcell,
												double ***pres_3dcell, double ***divu_3dcell);

void write_line_node_2d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												double ***vel_2dnode);

void write_line_node_3d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												double ****vel_3dnode);

void create_file(char *filename, double t, int ncycle, int *fileid)
{ 
/****
		 double sim_time;
		 mio_Attr attr;
		 char asim_time[] = "simulation_time";
		
		 sim_time = t;

		 mio_open_file(filename, mio_file_create, fileid);
		 attr.name = asim_time;
		 attr.datatype = mio_double;
		 attr.values = &sim_time;
		 attr.count = 1;
		 mio_write(*fileid, mio_attr, *fileid, &attr);
****/
		 return;
} 

void close_file(int *fileid)
{

/****
		if (*fileid >= 0) { 
				mio_close_file(*fileid);
				*fileid = -1;
		} 
o****/
		return;
 } 

void viz_dump(const char *probname, int *fileid,
						 int dim, double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle, 
							int nmat, 
							double ***vf_2dmat, double ***rho_2dmat, double ***ei_2dmat, double ***pres_2dmat,
							double **rho_2dcell, double **ei_2dcell, double **pres_2dcell, double **divu_2dcell, double ***vel_2dnode,
							double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat, double ****pres_3dmat,
							double ***rho_3dcell, double ***ei_3dcell, double ***pres_3dcell, double ***divu_3dcell, double ****vel_3dnode)
{
/******
		 int meshid, varid; 
		 char filename[32];
 
		 mio_Mesh_Var_Type vartype;
		 mio_Mesh_Var var;

		 if (dim == 2) {
				 vartype = mio_face;
		 }
		 else if (dim == 3) {
				 vartype = mio_zone;
		 }
		 sprintf(filename, "%s_%08d", probname, ncycle);
		 printf("create file = %s\n", filename);

		 create_file(filename, t, ncycle, fileid);
		
		 write_mesh_mat(*fileid, "mesh", 
										dim, xl_prob, xr_prob, ncell, nbdry,
										nmat, 
										vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat, rho_2dcell, ei_2dcell, pres_2dcell, divu_2dcell, vel_2dnode,
										vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat, rho_3dcell, ei_3dcell, pres_3dcell, divu_3dcell, vel_3dnode);

//     close_file(fileid);

		if (dim == 2) { 
				write_line_cell_2d(probname, xl_prob, xr_prob, ncell, nbdry, t, ncycle,
													 nmat, vf_2dmat, rho_2dcell, ei_2dcell, pres_2dcell, divu_2dcell);
				 
				write_line_node_2d(probname, xl_prob, xr_prob, ncell, nbdry, t, ncycle,
														vel_2dnode);

		}
		else if (dim == 3) { 
				write_line_cell_3d(probname, xl_prob, xr_prob, ncell, nbdry, t, ncycle,
													 nmat, vf_3dmat, rho_3dcell, ei_3dcell, pres_3dcell, divu_3dcell);
				 
				write_line_node_3d(probname, xl_prob, xr_prob, ncell, nbdry, t, ncycle,
														vel_3dnode);

		} 

		return;
 } 

void write_smesh(int fileid, char *meshname,
								int dim, double *xl_prob, double *xr_prob, int *ncell, int *meshid)
{
		 int i;
		 mio_Structured_Mesh mesh;

		 mio_init(fileid, mio_smesh, -1, &mesh);
		 mesh.name = meshname;
		 mesh.dims = dim;
		 mesh.datatype = mio_double;

		 mesh.element_centered = 0;

		 for (i = 0; i < dim; i++) {
				 mesh.sizes[i] = ncell[i];
				 mesh.dcoord[i] = (xr_prob[i] - xl_prob[i])/(double)ncell[i];
				 mesh.coordmin[i] = xl_prob[i];
		 }
		 mio_write(fileid, mio_smesh, fileid, &mesh);
		 *meshid = mesh.id;

*****/
		 return;
 }


void write_mesh_mat(int fileid, char *meshname, 
										int dim, double *xl_prob, double *xr_prob, int *ncell, int nbdry, 
				int nmat,  
										double ***vf_2dmat, double ***rho_2dmat, double ***ei_2dmat, double ***pres_2dmat,
										double **rho_2dcell, double **ei_2dcell, double **pres_2dcell, double **divu_2dcell, double ***vel_2dnode,
										double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat, double ****pres_3dmat,
										double ***rho_3dcell, double ***ei_3dcell, double ***pres_3dcell, double ***divu_3dcell, double ****vel_3dnode)
{ 

/******
		 int ghost_included = 1;
		 int i, j, k, m, meshid, varid;
		 int ncell_ext[3]; 
		 int nn, nn_ea_elem, nn_ea_k, n0, n00;  
		 int *nodelist_for_elem, *nodelist; 
		 int **nmat_2dcell, ***nmat_3dcell, **nmat2d, *nmat1d;  
		 long long nelem, nnode, lsize, lsize_node, offset; 
		 double dx[3], x, y, z; 
		 double *xyz[3]; 

		 mio_Unstructured_Mesh mesh;
		 mio_Coord *coord;
		 mio_Unstructured_Mesh_Type mesh_type; 
		 mio_Mesh_Var_Type vartype;

		 lsize = 1;
		 lsize_node = 1;
		 for (i = 0; i < dim; i++) { 
				 ncell_ext[i] = ncell[i] + nbdry + nbdry;
				 lsize *= ncell_ext[i];
				 lsize_node *+ (ncell_ext[i] + 1); 
		 } 
		 nelem = 1; 
		 nnode = 1;
		 for (i = 0; i < dim; i++) { 
				 nnode *= (ncell_ext[i] + 1);
				 nelem *= ncell_ext[i];
				 dx[i] = (xr_prob[i] - xl_prob[i])/(double)ncell[i];
		 } 
		 xyz[0] = (double *) malloc(nnode * dim * sizeof(double));
		 for (i = 1; i < dim; i++) {
				 xyz[i] = xyz[i-1] + nnode;
		 }
		 nn = 0;
		 if (dim == 2) { 
				 y = xl_prob[1] - (double)nbdry * dx[1];
				 for (j = 0; j <= ncell_ext[1]; j++) { 
						 x = xl_prob[0] - (double)nbdry * dx[0];
						 for (i = 0; i <= ncell_ext[0]; i++) { 
								 xyz[0][nn] = x;
								 xyz[1][nn] = y;
								 nn++;
								 x += dx[0];
						 }
						 y += dx[1];
				 }
				 mesh_type = mio_quad; 
				 nn_ea_elem = 4;
		 }  
		 else if (dim == 3) {
				 z = xl_prob[2] - (double)nbdry * dx[2];
				 for (k = 0; k <= ncell_ext[2]; k++) { 
						 y = xl_prob[1] - (double)nbdry * dx[1];
						 for (j = 0; j <= ncell_ext[1]; j++) {
								 x = xl_prob[0] - (double)nbdry * dx[0];
								 for (i = 0; i <= ncell_ext[0]; i++) {
										 xyz[0][nn] = x;
										 xyz[1][nn] = y;
										 xyz[2][nn] = z; 
										 nn++;
										 x += dx[0];
								 }
								 y += dx[1];
						 }  
						 z += dx[2]; 
				 }
				 mesh_type = mio_hex;
				 nn_ea_elem = 8; 
		 }
		 nodelist_for_elem = (int *) malloc(nelem * nn_ea_elem * sizeof(int));
		 nodelist = nodelist_for_elem;
		 if (dim == 2) { 
				 for (j = 0; j < ncell_ext[1]; j++) {
						 n0 = j * (ncell_ext[0] + 1); 
						 for (i = 0; i < ncell_ext[0]; i++) { 
								 nodelist[0] = n0 + i; 
								 nodelist[1] = nodelist[0] + 1;
								 nodelist[2] = nodelist[1] + (ncell_ext[0] + 1); 
								 nodelist[3] = nodelist[2] - 1;
								 
								 nodelist += nn_ea_elem;   
						 }
				 }
		 }
		 else if (dim == 3) {
				 nn_ea_k = (ncell_ext[1] + 1) * (ncell_ext[0] + 1);
				 for (k = 0; k < ncell_ext[2]; k++) {
						 n00 = k * nn_ea_k;
						 for (j = 0; j < ncell_ext[1]; j++) {
								 n0 = n00 + j * (ncell_ext[0] + 1);
								 for (i = 0; i < ncell_ext[0]; i++) {
										 nodelist[0] = n0 + i; 
										 nodelist[1] = nodelist[0] + 1;
										 nodelist[2] = nodelist[1] + (ncell_ext[0] + 1);
										 nodelist[3] = nodelist[2] - 1;

										 nodelist[4] = nodelist[0] + nn_ea_k;
										 nodelist[5] = nodelist[1] + nn_ea_k;
										 nodelist[6] = nodelist[2] + nn_ea_k;
										 nodelist[7] = nodelist[3] + nn_ea_k; 
										 
										 nodelist += nn_ea_elem;
								 } 
						 }
				 }
		 }
		 mio_init(fileid, mio_umesh, -1, &mesh);
		 coord = &(mesh.coord);
		 coord->datatype = mio_double;
		 for (i = 0; i < dim; i++) { 
				 coord->coord[i] = xyz[i];
		 }
		 mesh.name = meshname;
		 mesh.dims = dim;
		 mesh.idmin = 0;
		 mesh.datatype = mio_int;
		 mesh.type = mesh_type;
		 mesh.order_for_nodelist = 1;   

		 mesh.sizes[3] = nnode;
		 if (dim == 2) { 
				 mesh.sizes[1] = nelem;
				 mesh.nodelist_for_face = nodelist_for_elem;
		 }
		 else if (dim == 3) { 
				 mesh.sizes[0] = nelem;
				 mesh.nodelist_for_zone = nodelist_for_elem; 
		 }      
		 mio_write(fileid, mio_umesh, fileid, &mesh);
		 meshid = mesh.id;

		 if (dim == 2) { 
				 nmat_2dcell = (int **) malloc(ncell_ext[1] * sizeof(int *));
				 nmat_2dcell[0] = (int *) malloc(lsize * sizeof(int));
				 for (j = 1; j < ncell_ext[1]; j++) { 
						 nmat_2dcell[j] = nmat_2dcell[j-1] + ncell_ext[0];
				 }
				 for (j = 0; j < ncell_ext[1]; j++) { 
						 for (i = 0; i < ncell_ext[0]; i++) { 
								 nmat_2dcell[j][i] = 0;
								 for (m = 0; m < nmat; m++) { 
										 if (vf_2dmat[j][i][m] >= vfmin) { 
													nmat_2dcell[j][i]++;
										 }
								 }
						 }
				 }  
		 }
		 else if (dim == 3) { 
				 nmat_3dcell = (int ***) malloc(ncell_ext[2] * sizeof(int **));
				 nmat2d = (int **) malloc(ncell_ext[2] * ncell_ext[1] * sizeof(int *));
				 nmat1d = (int  *) malloc(lsize * sizeof(int));
				 
				 offset = 0;
				 for (k = 0; k < ncell_ext[2]; k++) { 
						 nmat2d[0] = nmat1d + offset;
						 for (j = 1; j < ncell_ext[1]; j++) {
								 nmat2d[j] = nmat2d[j-1] + ncell_ext[0]; 
						 }
						 nmat_3dcell[k] = nmat2d;
						 offset += (ncell_ext[0] * ncell_ext[1]);
						 nmat2d += ncell_ext[1]; 
				 }  
				 for (k = 0; k < ncell_ext[2]; k++) { 
						 for (j = 0; j < ncell_ext[1]; j++) { 
								 for (i = 0; i < ncell_ext[0]; i++) { 
										 nmat_3dcell[k][j][i] = 0;
										 for (m = 0; m < nmat; m++) { 
												 if (vf_3dmat[k][j][i][m] >= vfmin) { 
														 nmat_3dcell[k][j][i]++;
												 }
										 }
								 }
						 }
				 }
		 } 
		 if (dim == 2) {
				 vartype = mio_face;
		 }
		 else if (dim == 3) {
				 vartype = mio_zone;
		 }
		 if (dim == 2) { 
				 write_scalar("nmat",ghost_included, mio_int,    vartype, fileid, meshid, dim, ncell, nbdry, (void *)(nmat_2dcell[0]), &varid); 
				 write_scalar("rho", ghost_included, mio_double, vartype, fileid, meshid, dim, ncell, nbdry, (void *)(rho_2dcell[0]), &varid);
				 write_scalar("ei",  ghost_included, mio_double, vartype, fileid, meshid, dim, ncell, nbdry, (void *)(ei_2dcell[0]),  &varid);
				 write_scalar("pres",ghost_included, mio_double, vartype, fileid, meshid, dim, ncell, nbdry, (void *)(pres_2dcell[0]),&varid);

				 free(nmat_2dcell[0]);
				 free(nmat_2dcell); 
		 }
		 else if (dim == 3) { 
				 write_scalar("nmat",ghost_included, mio_int,    vartype, fileid, mesh.id, dim, ncell, nbdry, (void *)(nmat_3dcell[0][0]), &varid);
				 write_scalar("rho", ghost_included, mio_double, vartype, fileid, meshid, dim, ncell, nbdry, (void *)(rho_3dcell[0][0]), &varid);
				 write_scalar("ei",  ghost_included, mio_double, vartype, fileid, meshid, dim, ncell, nbdry, (void *)(ei_3dcell[0][0]),  &varid);
				 write_scalar("pres",ghost_included, mio_double, vartype, fileid, meshid, dim, ncell, nbdry, (void *)(pres_3dcell[0][0]),&varid);

				 free(nmat_3dcell[0][0]);
				 free(nmat_3dcell[0]);
				 free(nmat_3dcell); 
		 }
		 write_vecter("vel", ghost_included, mio_double, mio_node, fileid, meshid, dim, ncell, nbdry, vel_2dnode, vel_3dnode, &varid);

		 if (dim == 2) { 
	 write_mats_2d(fileid, nmat, ncell, nbdry,
								 nelem, nnode, xyz, nodelist_for_elem, 
											 vf_2dmat, rho_2dmat, ei_2dmat, pres_2dmat);
		 }
		 else if (dim == 3) { 
				 write_mats_3d(fileid, nmat, ncell, nbdry,
												 nelem, nnode, xyz, nodelist_for_elem,
												 vf_3dmat, rho_3dmat, ei_3dmat, pres_3dmat);
		 } 
		 free(xyz[0]);
		 free(nodelist_for_elem);

*****/

		 return;
 } 
/*****

void write_scalar(char *varname, int ghost_cell_included, mio_Data_Type datatype, mio_Mesh_Var_Type vartype,
									int fileid, int meshid, int dim, int *ncell, int nbdry,
									void *var1d, int *varid)
{   
		int i0, j0, k0, i, j, k, nbyte;
		int ncell_ext[3]; 
		long long lsize; 
		char *buffer, *buf, *v;
		mio_Mesh_Var var;
		
		for (i = 0; i < dim; i++) {
				ncell_ext[i] = ncell[i] + nbdry + nbdry;
		}
		if (datatype == mio_int) {
				nbyte = sizeof(int);
		}
		else if (datatype == mio_double) {
				nbyte = sizeof(double);
		}  
		if (!ghost_cell_included) {
				lsize = 1;  
				for (i = 0; i < dim; i++) { 
						lsize *= ncell[i];
				} 
				buffer = (char *) malloc((size_t)(lsize * nbyte)); 
				buf = buffer; 
				
				if (dim == 2) {  
						v = (char *)var1d + (long long)(ncell_ext[0] * nbdry * nbyte); 
						for (j0 = 0; j0 < ncell[1]; j0++) {
								v += (nbdry * nbyte); 
								for (i0 = 0; i0 < ncell[0]; i0++) {
										memcpy(buf, v, (size_t)nbyte);
										buf += nbyte;
										v   += nbyte; 
								}
						}
				}  
				else if (dim == 3) {  
						v = (char *)var1d + (long long)(ncell_ext[0] * ncell_ext[1] * nbyte);
						for (k0 = 0; k0 < ncell[2]; k0++) {
								v += (nbdry * ncell_ext[0] * nbyte); 
								for (j0 = 0; j0 < ncell[1]; j0++) {
										v += (nbdry * nbyte); 
										for (i0 = 0; i0 < ncell[0]; i0++) {
												memcpy(buf, v, (size_t)nbyte);
												buf += nbyte;
												v   += nbyte;
										}
								}
						}
				}
		 }
		 else {
				buffer = (char *)var1d;
		 } 
		 mio_init(fileid, mio_mesh_var, -1, &var);
		 var.name = varname;
		 var.mesh_ids[0] = meshid;
		 var.num_meshes = 1;
		 var.type = vartype;
		 var.datatype = datatype;
		 var.rank = 0;
		 var.comps[0].buffer = (void *)buffer;
		 mio_write(fileid, mio_mesh_var, fileid, &var);
		 *varid = var.id;
		 
		 if (!ghost_cell_included) { 
				 free(buffer);
		 } 
		 return;
}


void write_vecter(char *varname, int ghost_cell_included, mio_Data_Type datatype, mio_Mesh_Var_Type vartype,
									int fileid, int meshid, int dim, int *ncell, int nbdry, 
									double ***var_for_2d, double ****var_for_3d, int *varid)
{  
		int i0, j0, k0, i, j, k, c, nbyte; 
		int nnode_ext[3]; 
		long long lsize, offset; 
		char *buffer, *mybuffer[3], *buf[3];  
		mio_Mesh_Var var;
	
		if (datatype == mio_int) {
				nbyte = sizeof(int);
		}
		else if (datatype == mio_double) {
				nbyte = sizeof(double);
		} 
		lsize = 1; 
		if (!ghost_cell_included) {
				for (i = 0; i < dim; i++) {
						lsize *= (ncell[i] + 1);
				}
		}
		else { 
				for (i = 0; i < dim; i++) { 
						nnode_ext[i] = ncell[i] + nbdry + nbdry + 1;
						lsize       *= nnode_ext[i];
				}
		}
		buffer = (char *)malloc((size_t)(dim * lsize * nbyte)); 

		offset = 0; 
		for (i = 0; i < dim; i++) {
				mybuffer[i] = buffer + offset;
				buf[i]      = mybuffer[i];
				offset += (lsize * nbyte); 
		}
		if (!ghost_cell_included) { 
				if (dim == 2) { 
						for (j0 = 0; j0 < ncell[1]; j0++) {
								j = j0 + nbdry;
								for (i0 = 0; i0 < ncell[0]; i0++) {
										i = i0 + nbdry;
										for (c = 0; c < dim; c++) { 
												memcpy(buf[c], &(var_for_2d[j][i][c]), (size_t)nbyte); 
												buf[c] += nbyte;
										}
								}
						}  
				}  
				else if (dim == 3) { 
						for (k0 = 0; k0 < ncell[2]; k0++) { 
								k = k0 + nbdry; 
								for (j0 = 0; j0 < ncell[1]; j0++) {
										j = j0 + nbdry;
										for (i0 = 0; i0 < ncell[0]; i0++) {
												i = i0 + nbdry;
												for (c = 0; c < dim; c++) { 
														memcpy(buf[c], &(var_for_3d[k][j][i][c]), (size_t)nbyte);
														buf[c] += nbyte;
												} 
										}
								}
						}
				}
		}
		else {
			 if (dim == 2) {
					 for (j = 0; j < nnode_ext[1]; j++) {
							 for (i = 0; i < nnode_ext[0]; i++) {
									 for (c = 0; c < dim; c++) {
											 memcpy(buf[c], &(var_for_2d[j][i][c]), (size_t)nbyte);
											 buf[c] += nbyte;
									 }
							}
					 }
				}
				else if (dim == 3) {
						for (k = 0; k < nnode_ext[2]; k++) {
								for (j = 0; j < nnode_ext[1]; j++) {
										for (i = 0; i < nnode_ext[0]; i++) {
												for (c = 0; c < dim; c++) { 
														memcpy(buf[c], &(var_for_3d[k][j][i][c]), (size_t)nbyte);
														buf[c] += nbyte;
												} 
										}
								}
						}
				}
		} 
		mio_init(fileid, mio_mesh_var, -1, &var);
		var.name = varname;
		var.mesh_ids[0] = meshid;
		var.num_meshes = 1;
		var.type = vartype;
		var.datatype = datatype;
		var.rank = 1;
		var.comp_sizes[0] = dim;
		for (i = 0; i < dim; i++) {
				var.comps[i].which_comp[0] = i;
				var.comps[i].buffer = mybuffer[i];
		}
		mio_write(fileid, mio_mesh_var, fileid, &var);
		*varid = var.id;

		free(buffer);

		return;
}

void write_mpoly(int fileid, int dim, int pass,
									 int nmat, 
									 int nmixcell_mpoly, int *nmat_in_mixcell, int **matids_in_mixcell, 
									 int *nnode_in_mixcell, double **coords_in_mixcell,
									 int **nnode_for_minterface_in_mixcell, int ***nodes_for_minterface_in_mixcell,
									 int **nnode_for_mpoly_in_mixcell, int ***nodes_for_mpoly_in_mixcell,  // for 2D        
									 int **nface_for_mpoly_in_mixcell, int ***nnode_for_face_ea_mpoly_in_mixcell, // for 3D
									 int ***nodelist_for_face_ea_mpoly_in_mixcell)
{ 
		 char xpass[] = "xpass";
		 char ypass[] = "ypass";
		 char zpass[] = "zpass";
		 char *apass, name[16]; 
		 int i, j, k, n, nn, nn_max, nnode, nedge, nf, nface, nf_max, nface_max, nzone, offset, offset_node;
		 int m, nm, nm_max, mxpoly, mix;
		 int lsize, lsize_max, lsize_nodelist, lsize_facelist; 
		 
		 int *nface_for_zone, *facelist_for_zone, **nodes_for_face;
		 int *nodelist_for_edge, *nnode_for_face, *nodelist_for_face;
		 int *nodelist_s, *nodelist, *facelist, *nnode_face;
		 int *matids; 

		 double *xyz[3], *coords, *coords_s;

		 mio_Unstructured_Mesh mesh;
		 mio_Coord *coord_mio;
		 mio_Mesh_Var_Type vartype;
		 mio_Mesh_Var var;

		 if (fileid < 0) return;

		 if (nmat < 2)      return;
		 if (nmixcell_mpoly < 1) return;

		 if ((pass < 0) || (pass > 2)) { 
				 apass = NULL;
		 }
		 else if (pass == 0) { 
				 apass = xpass;
		 }
		 else if (pass == 1) { 
				 apass = ypass;
		 }
		 else if (pass == 2) { 
				 apass = zpass;
		 } 
		 if (dim == 2) {
				 vartype = mio_face;
		 }
		 else if (dim == 3) {
				 vartype = mio_zone;
		 }
		 nnode = 0;
		 for (mix = 0; mix < nmixcell_mpoly; mix++) {
				 nnode += nnode_in_mixcell[mix];
		 }
		 xyz[0] = (double *) malloc(nnode * dim * sizeof(double));
		 for (i = 1; i < dim; i++) {
				 xyz[i] = xyz[i-1] + nnode;
		 }
		 offset_node = 0;
		 for (mix = 0; mix < nmixcell_mpoly; mix++) {
				 nn = nnode_in_mixcell[mix];
				 coords_s = coords_in_mixcell[mix];
				 for (k = 0; k < nn; k++) {
						 coords = coords_s + (k * dim);
						 n = offset_node + k;
						 for (i = 0; i < dim; i++) {
									xyz[i][n] = coords[i];
						 }
				 }
				 offset_node += nn;
		 }
		 if (dim == 2) { 

//       write material interafces    

//         nodelist_for_edge = (int *) malloc((nmixcell_mpoly + nmixcell_mpoly) * sizeof(int));
//         for (m = 0; m < nmat - 1; m++) {
// 
//             offset_node = 0;
//             nedge = 0;
//             nodelist = nodelist_for_edge;
//             for (mxpoly = 0; mxpoly < nmixcell_mpoly; mxpoly++) {
//                 mix = mix_mpoly_to_mix[mxpoly];
//                 matids = matids_in_mixcell[mix];
//                 nmat   = nmat_n_mixcell[mix];
//                 for (k = 0; k < nmat-1; k++) { 
//                     if (matids[k] == matid) { 
//                         nn = nnode_for_minterface_in_mixcell[mxpoly][k];
//                         if (nn == 2) { 
//                             for (i = 0; i < nn; i++) { 
//                                 nodelist[i] = offset_node + nodes_for_minterface_in_mixcell[mxpoly][k][i];
//                             }
//                             nodelist += nn;
//                             nedge++;
//                         }
//                         break;
//                     }
//                 }
//                 offset_node += nnode_in_mixcell[mxpoly];
//             } 
//             mio_init(fileid, mio_umesh, -1, &mesh);
//             coord_mio = &(mesh.coord);
//             for (i = 0; i < dim; i++) {
//                 coord_mio->coord[i] = (void *) xyz[i];
//             }
//             coord_mio->datatype = mio_double;
//             if (apass) { 
//                 sprintf(name, "interface_mat_%d_%s", m, apass);
//             }
//             else { 
//                 sprintf(name, "interface_mat_%d", m);
//             } 
//             mesh.name = name;
//             mesh.dims = dim;
//             mesh.datatype = mio_int;
//             mesh.idmin    = 0; // 0-based
//             mesh.type = mio_bar;
//             mesh.sizes[3] = nnode;
//             mesh.sizes[2] = nedge;
//             mesh.nodelist_for_edge = nodelist_for_edge;
//             mio_write(fileid, mio_umesh, fileid, &mesh);
//         }
//         free(nodelist_for_edge); 
//
//

//       write material polygons

				 nm_max = 0;
				 nn_max = 0;
				 for (mxpoly = 0; mxpoly < nmixcell_mpoly; mxpoly++) { 
						 nm = nmat_in_mixcell[mxpoly];
						 if (nm_max < nm) nm_max = nm; 
						 for (k = 0; k < nm; k++) {
								 if (nn_max < nnode_for_mpoly_in_mixcell[mxpoly][k]) nn_max = nnode_for_mpoly_in_mixcell[mxpoly][k]; 
						 }
				 }
				 nnode_for_face    = (int *) malloc(nmixcell_mpoly * sizeof(int));
				 nodelist_for_face = (int *) malloc(nmixcell_mpoly * nm_max * nn_max * sizeof(int));
				 
				 for (m = 0; m < nmat; m++) {
						 nface = 0;
						 nodelist = nodelist_for_face;

						 offset_node = 0;
						 for (mxpoly = 0; mxpoly < nmixcell_mpoly; mxpoly++) {
								 nm = nmat_in_mixcell[mxpoly]; 
								 matids = matids_in_mixcell[mxpoly];
									
								 for (k = 0; k < nm; k++) { 
										 if (matids[k] == m) { 
												 nn = nnode_for_mpoly_in_mixcell[mxpoly][k];
												 nnode_for_face[nface] = nn;
												 for (i = 0; i < nn; i++) {
														 nodelist[i] = offset_node + nodes_for_mpoly_in_mixcell[mxpoly][k][i];
												 }
												 nodelist += nn; 
												 nface++;
												 break;
										 } 
								 }
								offset_node += nnode_in_mixcell[mxpoly];
						}  
						mio_init(fileid, mio_umesh, -1, &mesh);
						coord_mio = &(mesh.coord);
						for (i = 0; i < dim; i++) {
								coord_mio->coord[i] = (void *) xyz[i];
						}
						coord_mio->datatype = mio_double;
						if (apass) { 
								sprintf(name, "mat_%d_%s", m, apass);
						}
						else { 
								sprintf(name, "mat_%d", m);
						}
						mesh.name = name;
						mesh.dims = dim;
						mesh.datatype = mio_int;
						mesh.idmin    = 0; // 0-based
						mesh.type = mio_general_mesh;
						mesh.sizes[3] = nnode;
						mesh.sizes[1] = nface;
						mesh.num_nodes_for_face = nnode_for_face;
						mesh.nodelist_for_face  = nodelist_for_face;

						mio_write(fileid, mio_umesh, fileid, &mesh);
				}
				free(nnode_for_face);
				free(nodelist_for_face);
		 }
		 else if (dim == 3) { 

//      write material interfaces 



//        lsize_max = 0;  
//        nface_max = 0; 
//        for (m = 0; m < nmat; m++) {
//            lsize = 0;  
//            nface = 0; 
//            matid = matids_prob[m];
//            for (mxpoly = 0; mxpoly < nmixcell_mpoly; mxpoly++) {
//
//                nm = nmat_in_mixcell[mxpoly]; 
//                matids = matids_in_mixcell[mxpoly];
//                for (k = 0; k < nmat-1; k++) {
//                    if (matids[k] == m) { 
//                        lsize += nnode_for_minterface_in_mixcell[mxpoly][k]; 
//                        nface++; 
//			break;
//                    } 
//                }
//            }
//            if (nface > nface_max) nface_max = nface; 
//            if (lsize > lsize_max) lsize_max = lsize;
//        } 
//        nnode_for_face    = (int *) malloc(nface_max * sizeof(int)); 
//        nodelist_for_face = (int *) malloc(lsize_max * sizeof(int));
//
//        for (m = 0; m < nmat; m++) {
//            matid = matids_prob[m]; 
//            nodelist = nodelist_for_face; 
//            offset_node = 0;
//            nface = 0;
//
//            for (mxpoly = 0; mxpoly < nmixcell_mpoly; mxpoly++) {
//                mix = mix_mpoly_to_mix[mxpoly];
//                matids = matids_in_mixcell[mix];
//                nmat   = nmat_in_mixcell[mix];
//
//                for (k = 0; k < nmat-1; k++) {
//                    if (matids[k] == matid) { 
//                        nn = nnode_for_minterface_in_mixcell[mxpoly][k];
//                        for (i = 0; i < nn; i++) { 
//                            nodelist[i] = offset_node + nodes_for_minterface_in_mixcell[mxpoly][k][i]; 
//                        }
//                        nnode_for_face[nface] = nn; 
//                        nodelist += nn;
//                        nface++;
//                        break;
//                    }
//                }
//                offset_node += nnode_in_mixcell[mxpoly]; 
//            } 
//	    if (!nface) continue;
//
//            mio_init(fileid, mio_umesh, -1, &mesh);
//            coord_mio = &(mesh.coord);
//            for (i = 0; i < dim; i++) {
//                coord_mio->coord[i] = (void *) xyz[i];
//            }
//            coord_mio->datatype = mio_double;
//            if (apass) { 
//                sprintf(name, "interface_mat_%d_%s", m, apass);
//            }
//            else { 
//                sprintf(name, "interface_mat_%d", m);
//            }
//            mesh.name = name;
//            mesh.dims = dim;
//            mesh.datatype = mio_int;
//            mesh.idmin    = 0; // 0-based
//            mesh.type = mio_general_mesh;
//            mesh.sizes[3] = nnode;
//            mesh.sizes[1] = nface;
//            mesh.num_nodes_for_face = nnode_for_face; 
//            mesh.nodelist_for_face  = nodelist_for_face;
//            mio_write(fileid, mio_umesh, fileid, &mesh);
//        }
//        free(nnode_for_face); 






 
//      write material polyhedrons 
 
				lsize_max = 0; 
				nface_max = 0;
	nf_max = 0; 
				nm_max = 0;

				for (m = 0; m < nmat; m++) {
						nface = 0;
						lsize = 0; 
						for (mxpoly = 0; mxpoly < nmixcell_mpoly; mxpoly++) {
								nm = nmat_in_mixcell[mxpoly];
								matids = matids_in_mixcell[mxpoly];
								for (k = 0; k < nm; k++) {
										if (matids[k] == m) {
												nf = nface_for_mpoly_in_mixcell[mxpoly][k]; 
			if (nf > nf_max) nf_max = nf; 
												for (i = 0; i < nf; i++) { 
														nn = nnode_for_face_ea_mpoly_in_mixcell[mxpoly][k][i];           
														lsize += nn;
												}
												nface += nf;
												break;
										}
								}
						}
						if (nface > nface_max) nface_max = nface;
						if (lsize > lsize_max) lsize_max = lsize; 
				}
				nface_for_zone    = (int *) malloc(nmixcell_mpoly * sizeof(int));
				facelist_for_zone = (int *) malloc(nmixcell_mpoly * nf_max * sizeof(int)); 
				nnode_for_face    = (int *) malloc(nface_max * sizeof(int)); 
				nodelist_for_face = (int *) malloc(lsize_max * sizeof(int));

				nodes_for_face    = (int **) malloc(nf_max * sizeof(int *));

				for (m = 0; m < nmat; m++) {
						nzone = 0;
						nface = 0; 
						facelist   = facelist_for_zone;
						nnode_face = nnode_for_face; 
						nodelist   = nodelist_for_face;

			lsize_facelist = 0;  // for debug 
			lsize_nodelist = 0;  // for debug  
						
						offset_node = 0;
						for (mxpoly = 0; mxpoly < nmixcell_mpoly; mxpoly++) {
								nmat = nmat_in_mixcell[mxpoly];
								matids = matids_in_mixcell[mxpoly];

								for (k = 0; k < nmat; k++) {
										if (matids[k] == m) {
												nodelist_s = nodelist_for_face_ea_mpoly_in_mixcell[mxpoly][k]; 
												nf = nface_for_mpoly_in_mixcell[mxpoly][k];
												nface_for_zone[nzone] = nf; 
												offset = 0;
												for (j = 0; j < nf; j++) { 
														nodes_for_face[j] = nodelist_s + offset;
														offset += nnode_for_face_ea_mpoly_in_mixcell[mxpoly][k][j]; 
												} 
												for (j = 0; j < nf; j++) {
														nn = nnode_for_face_ea_mpoly_in_mixcell[mxpoly][k][j];
														nnode_face[j] = nn;
														facelist[j]   = nface + j;
														for (i = 0; i < nn; i++) { 
																nodelist[i] = offset_node + nodes_for_face[j][i];  
														}
														nodelist += nn;
					lsize_nodelist += nn; 
												}
												facelist   += nf; 
												nnode_face += nf; 
												nface      += nf;

			lsize_facelist += nf;

												nzone++;
												break;
										} 
								}    // k 
								offset_node += nnode_in_mixcell[mxpoly]; 
						}   // mix
						mio_init(fileid, mio_umesh, -1, &mesh);
						coord_mio = &(mesh.coord);
						for (i = 0; i < dim; i++) {
								coord_mio->coord[i] = (void *) xyz[i];
						}
						coord_mio->datatype = mio_double;
						if (apass) { 
								sprintf(name, "mat_%d_%s", m, apass);
						}
						else { 
								sprintf(name, "mat_%d", m);
						} 
						mesh.name = name;
						mesh.dims = dim;
						mesh.datatype = mio_int;
						mesh.idmin    = 0; // 0-based
						mesh.type = mio_general_mesh;
						mesh.sizes[3] = nnode;
						mesh.sizes[1] = nface;
						mesh.sizes[0] = nzone; 
						mesh.num_faces_for_zone = nface_for_zone;
						mesh.facelist_for_zone  = facelist_for_zone;
						mesh.num_nodes_for_face = nnode_for_face;
						mesh.nodelist_for_face  = nodelist_for_face;

						mio_write(fileid, mio_umesh, fileid, &mesh);

				}       // m 
				free(nface_for_zone);
				free(facelist_for_zone);
				free(nnode_for_face);
				free(nodelist_for_face);
				free(nodes_for_face); 
		 }  
		 free(xyz[0]);

		 return;
	}

void write_mats_2d(int fileid, int nmat, int *ncell, int nbdry, 
			 long long nelem_prob, long long nnode_prob, double *xy_prob[2],
			 int *nodelist_for_elem_prob,  
									 double ***vf_2dmat, double ***rho_2dmat, double ***ei_2dmat, double ***pres_2dmat)
{
		char meshname[32]; 
		int dim, i, j, m, e, n, nn, meshid, varid; 
		int nn_ea_elem;
		int ncell_ext[2];
		int *included_elem, *included_node, *node_prob2mat;
		int *nodes, *nodes_s; 
		int *nodelist_for_elem; 
		long long nelem, nnode, lsize; 
		double *vf, *rho, *ei, *pres, *xy[2];

		mio_Unstructured_Mesh mesh;
		mio_Coord *coord;

		int ghost_included = 1;

		dim = 2;

		lsize = 1;
		for (i = 0; i < dim; i++) {
				ncell_ext[i] = ncell[i] + nbdry + nbdry;
				lsize *= ncell_ext[i];
		} 
		nn_ea_elem = 4; 

		included_elem = (int *) malloc(nelem_prob * sizeof(int));  
		included_node = (int *) malloc(nnode_prob * sizeof(int)); 
		node_prob2mat = (int *) malloc(nnode_prob * sizeof(int));  
		
		for (m = 0; m < nmat; m++) { 
				for (e = 0; e < nelem_prob; e++) { 
						included_elem[e] = 0;
				} 
	for (n = 0; n < nnode_prob; n++) {
						included_node[n] = 0;
				} 
	nelem = 0; 
				for (j = 0; j < ncell_ext[1]; j++) { 
						for (i = 0; i < ncell_ext[0]; i++) { 
								if (vf_2dmat[j][i][m] >= vfmin) { 
				included_elem[nelem] = 1;
										nelem++;
								}
						}
				} 
				vf = (double *) malloc(nelem * 4 * sizeof(double));
				rho  = vf  + nelem;
				ei   = rho + nelem;
				pres = ei  + nelem;

	nelem = 0;
	for (e = 0; e < nelem_prob; e++) { 
			if (included_elem[e]) { 
					nodes = nodelist_for_elem_prob + (e * nn_ea_elem);
								for (i = 0; i < nn_ea_elem; i++) { 
							n = nodes[i];
										included_node[n] = 1;
								}
								nelem++;
						}
				}
	nnode = 0;
	for (n = 0; n < nnode_prob; n++) { 
						if (included_node[n]) { 
		node_prob2mat[n] = nnode; 
		nnode++;
			}
				}
	if (!nelem) continue;

	xy[0] = (double *) malloc(nnode * dim * sizeof(double));
	xy[1] = xy[0] + nnode; 
	for (i = 0; i < dim; i++) { 
						nn = 0;
			for (n = 0; n < nnode_prob; n++) {
								if (included_node[n]) { 
							xy[i][nn] = xy_prob[i][n];
				nn++; 
								} 
						} 
				} 
	nodelist_for_elem = (int *) malloc(nelem * nn_ea_elem * sizeof(int));

				nelem = 0; 	
				for (e = 0; e < nelem_prob; e++) {
						if (included_elem[e]) { 
								nodes_s = nodelist_for_elem_prob + (nn_ea_elem * e);
					nodes   = nodelist_for_elem      + (nn_ea_elem * nelem); 	
								for (i = 0; i < nn_ea_elem; i++) { 
							n = nodes_s[i];
				nodes[i] = node_prob2mat[n];
		}
		nelem++;
						} 
				} 
				mio_init(fileid, mio_umesh, -1, &mesh);
				coord = &(mesh.coord);
				coord->datatype = mio_double;
				for (i = 0; i < dim; i++) {
						coord->coord[i] = xy[i];
				}
	sprintf(meshname, "cleanCell_mat%d", m); 
				mesh.name = meshname;
				mesh.dims = dim;
				mesh.idmin = 0;
				mesh.datatype = mio_int;
				mesh.type = mio_quad;
				mesh.order_for_nodelist = 1;   

				mesh.sizes[3] = nnode;
				mesh.sizes[1] = nelem;
				mesh.nodelist_for_face  = nodelist_for_elem;
//	mesh.num_nodes_for_face = nnode_for_face;
				mio_write(fileid, mio_umesh, fileid, &mesh);
				meshid = mesh.id; 

				free(xy[0]);
				free(nodelist_for_elem);

				nelem = 0;
				for (j = 0; j < ncell_ext[1]; j++) {
						for (i = 0; i < ncell_ext[0]; i++) {
								if (vf_2dmat[j][i][m] >= vfmin) {
										vf[nelem]  = vf_2dmat[j][i][m]; 
										rho[nelem] = rho_2dmat[j][i][m];
										ei[nelem] = ei_2dmat[j][i][m];
										pres[nelem] = pres_2dmat[j][i][m]; 
										
										nelem++;
								}
						}
				}
				write_scalar("vf",  ghost_included, mio_double, mio_face, fileid, meshid, dim, ncell, nbdry, (void *)vf,  &varid);
				write_scalar("rho", ghost_included, mio_double, mio_face, fileid, meshid, dim, ncell, nbdry, (void *)rho,  &varid);
				write_scalar("ei",  ghost_included, mio_double, mio_face, fileid, meshid, dim, ncell, nbdry, (void *)ei,  &varid);
				write_scalar("pres",ghost_included, mio_double, mio_face, fileid, meshid, dim, ncell, nbdry, (void *)pres,&varid); 

				free(vf);
		} 
		free(included_elem);
		free(included_node);
		free(node_prob2mat); 

		return;

	} 

void write_mats_3d(int fileid, int nmat, int *ncell, int nbdry, 
			 long long nelem_prob, long long nnode_prob, double *xyz_prob[3],
			 int *nodelist_for_elem_prob,  
									 double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat, double ****pres_3dmat)
{
		char meshname[32]; 
		int dim, i, j, k, m, e, n, nn, meshid, varid; 
		int nn_ea_elem;
		int ncell_ext[3]; 
		int *included_elem, *included_node, *node_prob2mat;
		int *nodes, *nodes_s; 
		int *nodelist_for_elem; 
		long long nelem, nnode;
		double *vf, *rho, *ei, *pres, *xyz[3];

		mio_Unstructured_Mesh mesh;
		mio_Coord *coord;

		int ghost_included = 1; 

		dim = 3;
		for (i = 0; i < dim; i++) {
				ncell_ext[i] = ncell[i] + nbdry + nbdry;
		}
		nn_ea_elem = 8;

		included_elem = (int *) malloc(nelem_prob * sizeof(int));  
		included_node = (int *) malloc(nnode_prob * sizeof(int)); 

		node_prob2mat = (int *) malloc(nnode_prob * sizeof(int));  
		
		for (m = 0; m < nmat; m++) { 
				for (e = 0; e < nelem_prob; e++) { 
						included_elem[e] = 0;
				} 
	for (n = 0; n < nnode_prob; n++) {
						included_node[n] = 0;
				} 
	nelem = 0; 
	for (k = 0; k < ncell_ext[2]; k++) {  
						for (j = 0; j < ncell_ext[1]; j++) { 
								for (i = 0; i < ncell_ext[0]; i++) { 
										if (vf_3dmat[k][j][i][m] >= vfmin) {                 
						included_elem[nelem] = 1;
						nelem++; 
										} 
								} 
						}
				} 
				vf = (double *) malloc(nelem * 4 * sizeof(double));
				rho  = vf  + nelem;
				ei   = rho + nelem;
				pres = ei  + nelem;
 
	nelem = 0;
				for (e = 0; e < nelem_prob; e++) {
						if (included_elem[e]) {
								nodes = nodelist_for_elem_prob + (e * nn_ea_elem);
								for (i = 0; i < nn_ea_elem; i++) {
										n = nodes[i];
										included_node[n] = 1;
								}
								nelem++;
						}
				}
	nnode = 0;
	for (n = 0; n < nnode_prob; n++) { 
						if (included_node[n]) { 
		node_prob2mat[n] = nnode; 
		nnode++;
			}
				}
	if (!nelem) continue;

	xyz[0] = (double *) malloc(nnode * dim * sizeof(double));
	for (i = 1; i < dim; i++) { 
						xyz[i] = xyz[i-1] + nnode;
				}
	for (i = 0; i < dim; i++) { 
						nn = 0;
			for (n = 0; n < nnode_prob; n++) {
								if (included_node[n]) { 
							xyz[i][nn] = xyz_prob[i][n];
				nn++; 
								} 
						} 
				} 
	nodelist_for_elem = (int *) malloc(nelem * nn_ea_elem * sizeof(int));

				nelem = 0; 	
				for (e = 0; e < nelem_prob; e++) {
						if (included_elem[e]) { 
		nodes_s = nodelist_for_elem_prob + (nn_ea_elem * e); 
		nodes   = nodelist_for_elem      + (nn_ea_elem * nelem);
				 	for (i = 0; i < nn_ea_elem; i++) { 
										n = nodes_s[i]; 
							nodes[i] = node_prob2mat[n]; 
								}
		nelem++;
						} 
				} 
				mio_init(fileid, mio_umesh, -1, &mesh);
				coord = &(mesh.coord);
				coord->datatype = mio_double;
				for (i = 0; i < dim; i++) {
						coord->coord[i] = xyz[i];
				}
	sprintf(meshname, "cleanCell_mat%d", m); 
				mesh.name = meshname;
				mesh.dims = dim;
				mesh.idmin = 0;
				mesh.datatype = mio_int;
				mesh.type = mio_hex;
				mesh.order_for_nodelist = 1;   

				mesh.sizes[3] = nnode;
				mesh.sizes[0] = nelem;
				mesh.nodelist_for_zone = nodelist_for_elem;

				mio_write(fileid, mio_umesh, fileid, &mesh);
				meshid = mesh.id;

				free(xyz[0]);
				free(nodelist_for_elem);

				nelem = 0;
				for (k = 0; k < ncell_ext[2]; k++) { 
						for (j = 0; j < ncell_ext[1]; j++) {
								for (i = 0; i < ncell_ext[0]; i++) {
										if (vf_3dmat[k][j][i][m] >= vfmin) {
												vf[nelem]  = vf_3dmat[k][j][i][m];
												rho[nelem] = rho_3dmat[k][j][i][m];
												ei[nelem] = ei_3dmat[k][j][i][m];
												pres[nelem] = pres_3dmat[k][j][i][m];
		
												nelem++;
										} 
								}
						}
				}
				write_scalar("vf",  ghost_included, mio_double, mio_zone, fileid, meshid, dim, ncell, nbdry, (void *)vf,  &varid);
				write_scalar("rho", ghost_included, mio_double, mio_zone, fileid, meshid, dim, ncell, nbdry, (void *)rho, &varid);
				write_scalar("ei",  ghost_included, mio_double, mio_zone, fileid, meshid, dim, ncell, nbdry, (void *)ei,  &varid);
				write_scalar("pres",ghost_included, mio_double, mio_zone, fileid, meshid, dim, ncell, nbdry, (void *)pres,&varid);

				free(vf);

		} 
		free(included_elem);
		free(included_node);
		free(node_prob2mat); 

		return;

	} 


void write_a_polygon(int fid, char *name, int dim, int nnode_tot, double *coord,
										 int nnode, int *nodelist_for_face)
{ 
		 char name_nodeid[] = "nodeid";
		 int  nn, i, i1, n, n1, nedge, fileid;
		 int *nodelist_for_edge, *list_default;
		 double *x, *y, *z;
		 mio_Coord *mcoord;
		 mio_Unstructured_Mesh face;
		 mio_Mesh_Var var;

		 if (fid < 0) {
				 mio_open_file(name,mio_file_create, &fileid);
		 }
		 else {
				 fileid = fid;
		 }
		 x = (double *) malloc(dim * nnode_tot * sizeof(double)); 
		 y = x + nnode_tot;

		 if (dim == 2) { 
				 for (i = 0; i < nnode_tot; i++) {
						 x[i] = coord[i*dim];
						 y[i] = coord[i*dim+1];
				 }
		 }
		 else if (dim == 3) {
				 z = y + nnode_tot;
				 for (i = 0; i < nnode_tot; i++) {
						 x[i] = coord[i*dim];
						 y[i] = coord[i*dim+1];
						 z[i] = coord[i*dim+2];
				 }
		 }
		 mio_init(fileid, mio_umesh, -1, &face);
		 mcoord = &(face.coord);
		 mcoord->coord[0] = x;
		 mcoord->coord[1] = y;
		 if (dim == 3) {
				 mcoord->coord[2] = z;
		 }
		 face.name       = name;
		 face.dims       = dim;
		 face.idmin      = 0;
		 face.datatype   = mio_int;
		 mcoord->datatype = mio_double;

		 nedge = nnode;
		 nodelist_for_edge = (int *) malloc((nedge + nedge) * sizeof(int));
		 for (i = 0; i < nedge; i++) {
				 i1 = (i + 1) % nnode;
				 n  = nodelist_for_face[i];
				 n1 = nodelist_for_face[i1];
				 nodelist_for_edge[i+i]   = n;
				 nodelist_for_edge[i+i+1] = n1;
		 }
		 face.type     = mio_bar;
		 face.sizes[2] = nedge;
		 face.sizes[3] = nnode_tot;
		 face.nodelist_for_edge = nodelist_for_edge;

		 mio_write(fileid, mio_umesh, fileid, &face);

		 list_default = (int *) malloc(nnode_tot * sizeof(int));
		 for (i = 0; i < nnode_tot; i++) { 
				 list_default[i] = i;
		 }
		 mio_init(fileid, mio_mesh_var, -1, &var);
		 var.name = name_nodeid;
		 var.mesh_ids[0] = face.id;
		 var.num_meshes = 1;
		 var.type = mio_node;
		 var.rank = 0;
		 var.datatype = mio_int;
		 (var.comps[0]).buffer = list_default;
		 mio_write(fileid, mio_mesh_var, fileid, &var);

		 free(x);
		 free(nodelist_for_edge);

		 if (fid < 0) {
				 mio_close_file(fileid);
		 }
		 free(list_default);

		 return;
 }  

void write_line_cell_2d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												int nmat, double ***vf_2dmat, 
												double **rho_2dcell, double **ei_2dcell, 
												double **pres_2dcell, double **divu_2dcell)
{  

		 char filename[128];
		 int ncell_ext[2], dim, i, j, m;
		 int **nmat_2dcell;
		 double dx, x; 
		 FILE *fp;

		 dim = 2;

		 for (i = 0; i < dim; i++) {
				 ncell_ext[i] = ncell[i] + nbdry + nbdry;
		 }
		 nmat_2dcell    = (int **) malloc(ncell_ext[1] * sizeof(int *));
		 nmat_2dcell[0] = (int  *) malloc(ncell_ext[1] * ncell_ext[0] * sizeof(int));
		 for (j = 1; j < ncell_ext[1]; j++) { 
				 nmat_2dcell[j] = nmat_2dcell[j-1] + ncell_ext[0];
		 }
		 for (j = 0; j < ncell_ext[1]; j++) {
				 for (i = 0; i < ncell_ext[0]; i++) {
						 nmat_2dcell[j][i] = 0;
						 for (m = 0; m < nmat; m++) {
								 if (vf_2dmat[j][i][m] >= vfmin) {
										 nmat_2dcell[j][i]++;
								 }
						 }
				 }
		 }
		 sprintf(filename, "%s_%08d_cell.txt", probname, ncycle);

		 fp = fopen(filename, "w");

		 fprintf(fp, "# t = %e\n", t);
		 fprintf(fp, "# ncycle = %d\n", ncycle);
		 fprintf(fp, "# ncell = %d\n", ncell[0]);

		 fprintf(fp, "    i       x         rho         ei        pres    divu    nmat \n");

		 dx = (xr_prob[0] - xl_prob[0])/ncell[0];
		 j  = nbdry + ncell[1]/2;
		 x  = - (double) nbdry * dx + 0.5 * dx;
		 for (i = 0; i < ncell_ext[0]; i++) {
				 fprintf(fp, "%6d %13e %13e %13e %13e %13e %2d \n", i, x,
										 rho_2dcell[j][i],
											ei_2dcell[j][i],
										pres_2dcell[j][i],
				divu_2dcell[j][i],
										nmat_2dcell[j][i]);
				 x += dx; 
		 } 
		 fclose(fp);

		 free(nmat_2dcell[0]);
		 free(nmat_2dcell);

		 return;
	} 

void write_line_cell_3d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												int nmat, double ****vf_3dmat, 
												double ***rho_3dcell, double ***ei_3dcell,
												double ***pres_3dcell, double ***divu_3dcell)
{

		 char filename[128];
		 int ncell_ext[3], dim, i, j, k, m;
		 int ***nmat_3dcell, **nmat2d, *nmat1d; 
		 long long offset, lsize; 
		 double dx, x;
		 FILE *fp;

		 dim = 3;

		 lsize = 1;
		 for (i = 0; i < dim; i++) {
				 ncell_ext[i] = ncell[i] + nbdry + nbdry;
				 lsize *= ncell_ext[i]; 
		 }
		 nmat_3dcell = (int ***) malloc(ncell_ext[2] * sizeof(int **));
		 nmat2d = (int **) malloc(ncell_ext[2] * ncell_ext[1] * sizeof(int *));
		 nmat1d = (int  *) malloc(lsize * sizeof(int));

		 offset = 0;
		 for (k = 0; k < ncell_ext[2]; k++) {
				 nmat2d[0] = nmat1d + offset;
				 for (j = 1; j < ncell_ext[1]; j++) {
						 nmat2d[j] = nmat2d[j-1] + ncell_ext[0];
				 }
				 nmat_3dcell[k] = nmat2d;
				 offset += (ncell_ext[0] * ncell_ext[1]);
				 nmat2d += ncell_ext[1];
		 }
		 for (k = 0; k < ncell_ext[2]; k++) {
				 for (j = 0; j < ncell_ext[1]; j++) {
						 for (i = 0; i < ncell_ext[0]; i++) {
								 nmat_3dcell[k][j][i] = 0;
								 for (m = 0; m < nmat; m++) {
										 if (vf_3dmat[k][j][i][m] >= vfmin) {
												 nmat_3dcell[k][j][i]++;
										 }
								 }
						 }
				 }
		 }
		 sprintf(filename, "%s_%08d_cell.txt", probname, ncycle);

		 fp = fopen(filename, "w");

		 fprintf(fp, "# t = %e\n", t);
		 fprintf(fp, "# ncycle = %d\n", ncycle);
		 fprintf(fp, "# ncell = %d\n", ncell[0]);

		 fprintf(fp, "  ic     xc        rho        ei        p     divu   nmat \n");

		 dx = (xr_prob[0] - xl_prob[0])/ncell[0];
		 j  = nbdry + ncell[1]/2;
		 k  = nbdry + ncell[2]/2; 
		 x  = - (double) nbdry * dx + 0.5 * dx;
		 for (i = 0; i < ncell_ext[0]; i++) {
				 fprintf(fp, "%6d %13e %13e %13e %13e %13e %2d \n", i, x,
										 rho_3dcell[k][j][i],
											ei_3dcell[k][j][i],
										pres_3dcell[k][j][i],
				divu_3dcell[k][j][i], 
										nmat_3dcell[k][j][i]);
				 x += dx;
		 }
		 fclose(fp);

		 free(nmat_3dcell[0][0]);
		 free(nmat_3dcell[0]);
		 free(nmat_3dcell); 

		 return;
	}

void write_line_node_2d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												double ***vel_2dnode)
{
		 char filename[128];
		 int nnode_ext[2], dim, i, j;
		 double dx, x;
		 FILE *fp;

		 dim = 2;

		 for (i = 0; i < dim; i++) {
				 nnode_ext[i] = ncell[i] + nbdry + nbdry + 1;
		 }
		 sprintf(filename, "%s_%08d_node.txt", probname, ncycle);

		 fp = fopen(filename, "w");

		 fprintf(fp, "# t = %e\n", t);
		 fprintf(fp, "# ncycle = %d\n", ncycle);
		 fprintf(fp, "# ncell = %d\n", ncell[0]);

		 fprintf(fp, "    in       xl         ux          uy\n");

		 dx = (xr_prob[0] - xl_prob[0])/ncell[0];
		 j  = nbdry + ncell[1]/2;
		 x  = - (double) nbdry * dx;
		 for (i = 0; i < nnode_ext[0]; i++) {
				 fprintf(fp, "%6d %13e %13e %13e\n", i, x,
										vel_2dnode[j][i][0], vel_2dnode[j][i][1]);
				 x += dx;
		 }
		 fclose(fp);

		 return;
	}


												
void write_line_node_3d(const char *probname,
												double *xl_prob, double *xr_prob, int *ncell, int nbdry, double t, int ncycle,
												double ****vel_3dnode)
{  
		 char filename[128];
		 int nnode_ext[3], dim, i, j, k;
		 double dx, x; 
		 FILE *fp;

		 dim = 3;

		 for (i = 0; i < dim; i++) { 
				 nnode_ext[i] = ncell[i] + nbdry + nbdry + 1;
		 } 
		 sprintf(filename, "%s_%08d_node.txt", probname, ncycle);
		 
		 fp = fopen(filename, "w");
		 
		 fprintf(fp, "# t = %e\n", t);
		 fprintf(fp, "# ncycle = %d\n", ncycle);
		 fprintf(fp, "# ncell = %d\n", ncell[0]);

		 fprintf(fp, "    i       x         ux          uy        uz\n");

		 dx = (xr_prob[0] - xl_prob[0])/ncell[0];
		 j  = nbdry + ncell[1]/2;
		 k  = nbdry + ncell[2]/2; 
		 x  = - (double) nbdry * dx;
		 for (i = 0; i < nnode_ext[0]; i++) {
				 fprintf(fp, "%6d %13e %13e %13e %13e\n", i, x,
										vel_3dnode[k][j][i][0], vel_3dnode[k][j][i][1], vel_3dnode[k][j][i][2]);
				 x += dx; 
		 }
		 fclose(fp);

		 return;
	}


******/ 
