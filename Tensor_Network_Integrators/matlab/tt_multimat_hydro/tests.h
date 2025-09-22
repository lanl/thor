#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _b : _a; })

double random_int(const int imin, const int imax);
double random_double(const double xmin, const double xmax);
void create_random_mesh(int dim, int *ncell, int *nbdry,
                        int *ncell_ext, int *nnode, int *nnode_ext, 
                        double *xl, double *xr, double *dx, const int verbose);
void random_field_3dmesh(const int *dims, double ****field, 
                         const double xmin, const double xmax);
void random_field_3dmesh_int(const int *dims, int ****field, 
                             const int imin, const int imax);
void allocate_velocity_3dmesh(const int *dims, double *****vv);
