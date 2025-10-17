#ifndef _UTIL_
#define _UTIL_

#ifdef __cplusplus
extern "C" {
#endif

//          | x1  y1  z1 |
//    v =   | x2  y2  z2 | 
//          | x3  y3  z3 |

#define DET3(v, x1, y1, z1, x2, y2, z2, x3, y3, z3)   \
{   v = x1 * y2 * z3 + x3 * y1 * z2 + x2 * y3 * z1 -  \
        x3 * y2 * z1 - x1 * y3 * z2 - x2 * y1 * z3;   \
} 

#define ASSIGN_3D_FORM(type, rho3d, rho2d, rho, sizes)         \
{   int j, k;                                                  \
    long long offset;                                          \
    if (!(rho3d)) {                                            \
        rho3d = (type ***)malloc(sizes[2] * sizeof(type **));  \
    }                                                          \
    if (!rho2d) {                                              \
        rho2d = (type **) malloc(sizes[1] * sizes[2] * sizeof(type *));   \
    }                                              \
    if (!rho) {                                    \
        rho = (type *) malloc(sizes[0] * sizes[1] * sizes[2] * sizeof(type)); \
    }                                              \
    offset = 0;                                    \
    for (k = 0; k < sizes[2]; k++) {               \
        rho2d[0] = rho + offset;                   \
        for (j = 1; j < sizes[1]; j++) {           \
            rho2d[j] = rho2d[j-1] + sizes[0];      \
        }                                          \
        (rho3d)[k] = rho2d;                          \
        offset += (sizes[0] * sizes[1]);           \
        rho2d  += sizes[1];                        \
    }                                              \
}

#define ASSIGN_4D_FORM(type, sizes, rho4d, tmp3d, tmp2d, tmp1d_to_clean)    \
{   int i, j, k;                                                            \
    long long offset;                                              \
    if (!(rho4d)) {                                                  \
        rho4d = (type ****) malloc(sizes[3] * sizeof(type ***));   \
    }                                                              \
    if (!tmp3d) {                                                                  \
        tmp3d = (type ***) malloc(sizes[3] * sizes[2] * sizeof(type **));          \
    }                                                                              \
    if (!tmp2d) {                                                                  \
        tmp2d = (type **) malloc(sizes[3] * sizes[2] * sizes[1] * sizeof(type *)); \
    }                                                                              \
    if (!tmp1d_to_clean) {                               \
        tmp1d_to_clean = (type *) malloc(sizes[3] * sizes[2] * sizes[1] * sizes[0] * sizeof(type)); \
    }                                                    \
    offset = 0;                                          \
    for (k = 0; k < sizes[3]; k++) {                     \
        tmp3d[0] = tmp2d + offset;                       \
        for (j = 1; j < sizes[2]; j++) {                 \
            tmp3d[j] = tmp3d[j-1] + sizes[1];            \
        }                                                \
        (rho4d)[k] = tmp3d;                              \
        offset += (sizes[1] * sizes[2]);                 \
        tmp3d  += sizes[2];                              \
    }                                                    \
    if (!tmp1d_to_clean) {                               \
        tmp1d_to_clean = (type *) malloc(sizes[3] * sizes[2] * sizes[1] * sizes[0] * sizeof(type)); \
    }                                                    \
    offset = 0;                                          \
    for (k = 0; k < sizes[3]; k++) {                     \
        for (j = 0; j < sizes[2]; j++) {                 \
            for (i = 0; i < sizes[1]; i++) {             \
                (rho4d)[k][j][i] = tmp1d_to_clean + offset; \
                offset += sizes[0];                      \
            }                                            \
        }                                                \
    }                                                    \
} 

struct ipair_t {
       int val;
       int index;
};
typedef struct ipair_t ipair_t;

struct rpair_t {
       double val;
       int index;
};
typedef struct rpair_t rpair_t;

typedef struct { double x, y; } vec_t;
typedef struct { int len, alloc; vec_t *v; } poly_t;

/////////////////////////////////////////////////////////////////////
double dot(vec_t *a, vec_t *b);
double cross(vec_t *a, vec_t *b);
vec_t *vsub(vec_t *a, vec_t *b, vec_t *res);

int left_of(vec_t *a, vec_t *b, vec_t *c);

int line_sect(vec_t *x0, vec_t *x1, vec_t *y0, vec_t *y1, vec_t *res);

poly_t *poly_new();

void poly_free(poly_t *p);

void poly_append(poly_t *p, vec_t *v);
int poly_winding(poly_t *p);

void poly_edge_clip(poly_t *sub, vec_t *x0, vec_t *x1, int left, poly_t *res);

poly_t *poly_clip(poly_t *sub, poly_t *clip);

void order_nodes_along_norm(int dim, double *norm,
                            int nnode, double *coords,
                            int *node_order, double *ds_ea_node);

void cal_poly_area(int nnode, double *coords, int nnode_poly, int *nodelist, double *area);
int rz_area(int nn, double *rz, double *vol);

void cal_ctr0(int nface, int nnode, int *facelist_for_the_zone,
             int *num_nodes_for_face, int *nodelist_for_face,
             double *coord, double *ctr, double *vol); 

int ipaircmp(const void *pointa, const void *pointb);
int rpaircmp(const void *pointa, const void *pointb);


void binarySearch(int *loc, int *found, int Item, int *Array, int right, int left);

int binaryInsert(int *list, int *size, int value, int *found);
int insertItem(int loc, int value, int *list, int *size);
int llicomp(const void *pa, const void *pb);

void r8sort(double *psi, int *order_tar_to_src, int n);

void find_min_max(int dim, int nnode, double *coords, double *xmin, double *xmax); 

#ifdef __cplusplus
}
#endif
#endif

