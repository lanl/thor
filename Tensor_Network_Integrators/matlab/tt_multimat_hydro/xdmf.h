int delete_file_if_exists(const char *filename);

void xdmf_new_file(const char *filename);


int count_lines_in_file(const char *filename);


int read_file_into_lines(const char *filename, char ***lines);


FILE* xdmf_open_file_append(const char *filename);


int xdmf_close_file(FILE* file);

void xdmf_write_mesh(FILE *file_id, const int dim, const double xl[], 
        const double xr[], const int ncell[], const int nbdry);

void xdmf_write_dataset(FILE *file_id, 
        const int dim, const int *ncell,
        const double *xl, const double *xr,
        const char *fname_h5, const char *timestep_group, 
        const char *varname, const int var_type);

void xdmf_write_dataset_int(FILE *file_id, 
        const int dim, const int *ncell,
        const double *xl, const double *xr,
        const char *fname_h5, const char *timestep_group, 
        const char *varname, const int var_type);

void xdmf_write_timestamp(FILE *file_id, const double timestamp);

void xdmf_close_group(FILE *file_id);

void write_mockup_mesh_data();

