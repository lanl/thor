#define CELL_VARIABLE 101
#define NODE_VARIABLE 102

#define ASSIGN_2D_FORM(type, rho2d, sizes)     \
{   int j;                                     \
    long long offset;                          \
    rho2d = (type **) malloc(sizes[1] * sizeof(type *));   \
    rho2d[0] = (type *) malloc(sizes[0] * sizes[1] * sizeof(type)); \
    offset = 0;                                \
    for (j = 1; j < sizes[1]; j++) {           \
        rho2d[j] = rho2d[j-1] + sizes[0];      \
    }                                          \
}

hid_t h5_new_file(const char *filename);
hid_t h5_open_existing_rdwr(const char *filename);
hid_t h5_open_existing_rdonly(const char *filename);

int h5_count_groups(const hid_t fogr_id);
hid_t h5_open_group(hid_t file_id, const char *group_name);

herr_t h5_write_1d(hid_t fogr_id, const char *dataset_name, double *data, const int len);
herr_t h5_write_2d(hid_t fogr_id, const char *dataset_name, double **data, const int ncell[]);
herr_t h5_write_3d(hid_t fogr_id, const char *dataset_name, double ***data, const int ncell[]);
herr_t h5_write_4d(hid_t fogr_id, const char *dataset_name, double ****data, const int ncell[]);
herr_t h5_write_1d_int(hid_t fogr_id, const char *dataset_name, int *data, const int len);
herr_t h5_write_2d_int(hid_t fogr_id, const char *dataset_name, int **data, const int ncell[]);
herr_t h5_write_3d_int(hid_t fogr_id, const char *dataset_name, int ***data, const int ncell[]);
herr_t h5_write_2d_llong(hid_t fogr_id, const char *dataset_name, long long **data, const int ncell[]);
herr_t h5_write_3d_llong(hid_t fogr_id, const char *dataset_name, long long ***data, const int ncell[]);

double** h5_read_2d(hid_t group_id, const char *dataset_name, int *dims_exp);
int** h5_read_2d_int(hid_t group_id, const char *dataset_name, int *dims_exp);

herr_t h5_write_attr(hid_t fogr_id, const char *attr_name, double val);
herr_t h5_write_attr_int(hid_t fogr_id, const char *attr_name, int ival);
herr_t h5_write_attr_1d(hid_t fogr_id, const char *attr_name, double *data, const int len);
herr_t h5_write_attr_1d_int(hid_t fogr_id, const char *attr_name, int *data, const int len);

int h5_read_attr_int(hid_t group_id, const char *attr_name);
void h5_read_attr_1d(hid_t group_id, const char *attr_name, double *attr_data, const int dim);
void h5_read_attr_1d_int(hid_t group_id, const char *attr_name, int *attr_data, const int dim);

void copy_velocity_component_2d (double ***vel_for_2dnode,
        double ***var_for_2dnode, const int c, const int *dims);
void copy_velocity_component_3d (double ****vel_for_3dnode,
        double ****var_for_3dnode, const int c, const int *dims);

