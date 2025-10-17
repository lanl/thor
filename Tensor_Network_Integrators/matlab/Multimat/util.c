#include "globals.h"
#include "util.h"

static double tiny   = 1.0e-30; 

void r8sort(double *psi, int *order_tar_to_src, int n)
{
    rpair_t *p, *list;
    size_t size;
    int  i;

    list = (rpair_t *) malloc(n * sizeof(rpair_t));

    for (i = 0; i < n; i++) {
        p = list + i;
        p->val = psi[i];
        p->index = i;
    }
    size = sizeof(rpair_t);
    qsort(list, (size_t)n, size, rpaircmp);
    for (i = 0; i < n; i++) {
        order_tar_to_src[i] = list[i].index;
    }
    free(list);

    return;
  }

void binarySearch(int *loc, int *found, int Item, int *Array, int right, int left)
{
/****  found is 1 if the item exist otherwise 0
    * loc is where the item is if found
    * Item is the value you are looking for
    * right is the size -1 at first when called
    * left is 0 when called
    * Array is the sorted array of elements.
****/
    int middle = (right + left) /2;
    *found = 0;
    if (right < left) { 
        (*loc) = 0; 
        return;
    } 
    (*loc) = right;
    if (*loc < 0) {
       (*loc) = 0;
    }
    if ((middle >= 0) && (Array[middle] == Item)) {
        (*loc) = middle;
        (*found) = 1;
    }
    else if (right > left) {
        if (Array[middle] < Item) {
            binarySearch(loc, found, Item, Array, right, middle+1);
        }
        else {
            binarySearch(loc, found, Item, Array, middle-1, left);
        }
    }
    else if (Array[middle] < Item) {
       (*loc) = right + 1;
    }
}

int binaryInsert(int *list, int *size, int value, int *found)
{
/*  If value is found, do nothing. Otherwise insert value */

    int loc, right, left;

    *found = 0;
    right = *size - 1;
    left  = 0;
    binarySearch(&loc, found, value, list, right, left);

    if (!(*found)) {
        insertItem(loc, value, list, size);
    }
    return 0;
}

int insertItem(int loc, int value, int *list, int *size)
{
    int k;

    for (k = *size - 1; k >= loc; k--) {
        list[k+1] = list[k];
    }
    list[loc] = value;
    (*size)++;
  
    return 0;
}

int llicomp(const void *pa, const void *pb)
{
    int a, b, va, vb;
    a  = *((int *)pa);
    b  = *((int *)pb);

    if (a < b) {
        return -1;
    }
    else if (a == b) {
        return 0;
    }
    else if (a > b) {
       return 1;
    }
    return 0;
 }

void order_nodes_along_norm(int dim, double *norm,
                            int nnode, double *coords, 
                            int *node_order, double *ds_ea_node)
{ 
     int     i, k, n; 
     double  *c; 
     rpair_t *p, *list; 
   
     list = (rpair_t *) malloc(nnode * sizeof(rpair_t));

     for (i = 0; i < nnode; i++) { 
         c = coords + (i * dim);
         ds_ea_node[i] = 0.0;  
         for (k = 0; k < dim; k++) { 
             ds_ea_node[i] += (norm[k] * c[k]);
         }
     }
     for (i = 0; i < nnode; i++) {
         p = list + i;
         p->val = ds_ea_node[i];
         p->index = i;
     }
     n = (int) sizeof(rpair_t);
     qsort(list, nnode, n, rpaircmp);

     for (i = 0; i < nnode; i++) {
         node_order[i] = list[i].index;
     }
     free(list);

     return;
  }  

int ipaircmp(const void *pointa, const void *pointb)
{
     ipair_t *a, *b;

     a = (ipair_t *)pointa;
     b = (ipair_t *)pointb;

     if (a->val > b->val) {
         return 1;
     }
     else if (a->val < b->val) {
         return -1;
     }
     return 0;
 }

int rpaircmp(const void *pointa, const void *pointb)
{
     rpair_t *a, *b;

     a = (rpair_t *)pointa;
     b = (rpair_t *)pointb;

     if (a->val > b->val) {
         return 1;
     }
     else if (a->val < b->val) {
         return -1;
     }
     return 0;
 }


// TODO: ask Will why nnode is not used here
void cal_poly_area(int nnode, double *coords, int nnode_poly, int *nodelist, double *area)
{
    int k, k1, n, n1;
    double darea, *p0, *p1;

    *area  = 0.0;

    if (nodelist) { 
        for (k = 0; k < nnode_poly; k++) {
            k1 = (k + 1) % nnode_poly;
            n  = nodelist[k];
            n1 = nodelist[k1]; 
            p0 = coords + (n  + n);
            p1 = coords + (n1 + n1);
            darea = p0[0] * p1[1] - p1[0] * p0[1];
            *area += darea;
        }
    }
    else { 
        for (k = 0; k < nnode_poly; k++) {
            k1 = (k + 1) % nnode_poly;
            p0 = coords + (k  + k);
            p1 = coords + (k1 + k1);
            darea = p0[0] * p1[1] - p1[0] * p0[1];
            *area += darea;
        }
    }
    *area *= 0.5;
    if (*area < 0.0) *area = - (*area);

    return;
}
          

int rz_area(int nn, double *rz, double *vol)
{
    int k, k1;
    double *pts[20], *p0, *p1;
    double r0, z0, r1, z1, dz, r02, r12, r0r1, dv;

    for (k = 0; k < nn; k++) {
        pts[k] = rz + (k + k);
    }
    *vol = 0.0;
    for (k = 0; k < nn; k++) {
        k1 = (k + 1) % nn;

        p0 = pts[k];
        p1 = pts[k1];
        r0 = p0[0];
        z0 = p0[1];
        r1 = p1[0];
        z1 = p1[1];

//      dv = 2 * integral (r dr dz)

        dz   = z1 - z0;
        r02  = r0 * r0;
        r12  = r1 * r1;
        r0r1 = r0 * r1;


        dv = dz *(r12 + r0r1 + r02);
        *vol += dv;
     }
     *vol /= 3.0;
     if (*vol < 0.0) *vol = -(*vol);

     return 0;
 }

void cal_ctr0(int nface, int nnode, int *facelist_for_the_zone,  
             int *num_nodes_for_face, int *nodelist_for_face,  
             double *coord, double *ctr, double *vol)
{ 
//  This function calculate centroid nad volume of the polyhedron poly.

    int  sz3int, found;
    int  fid, dim, nnode_tet, nnode_tri;
    int  d, i, k, n, n0, nn, ntri, itri, ntri_this_face, ntet, offset; 
    int  *nodelist_for_tri;
    int  nodelist_for_tet[4];
    int  *nodelist, *nodelist_s;
    double factor, v;     
    double *c, *c0, *r1, *r2, *r3, *r[4];
    double *ctr_tet, *vol_tet;

    sz3int = 3 * sizeof(int);
    nnode_tet = 4;
    nnode_tri = 3;
    dim       = 3;

    ntri = 0;
    for (i = 0; i < nface; i++) { 
        nn = num_nodes_for_face[i];
        ntri += (nn - 2);
    } 
    nodelist_for_tri = (int *) malloc(nnode_tri * ntri * sizeof(int));

    n0 = nodelist_for_face[0];
    itri = 0;
    nodelist   = nodelist_for_tri;
    nodelist_s = nodelist_for_face;

    offset = 0;
    for (i = 0; i < nface; i++) { 
        nn = num_nodes_for_face[i];
        found = 0;
        for (k = 0; k < nn; k++) { 
            if (nodelist_s[k] == n0) { 
                found = 1;
                break;
            }
        }
        if (!found) { 
            ntri_this_face = nn - 2;
            for (k = 0; k < ntri_this_face; k++) {
                n = nodelist_s[0]; 
                nodelist[0] = n; 
                n = nodelist_s[k+1];
                nodelist[1] = n; 
                n = nodelist_s[k+2];
                nodelist[2] = n;
                nodelist += nnode_tri;
                offset   += nnode_tri; 
            }
            itri += ntri_this_face;
        }
        nodelist_s += nn;
    }  
    assert(offset <= nnode_tri * ntri); 

    assert(itri < ntri);

    ntet = itri;
    vol_tet = (double *) malloc(ntet * sizeof(double));

    ctr_tet = (double *) malloc(ntet * sizeof(double) * dim);
 
    for (k = 0; k < ntet; k++) { 
        nodelist_s = nodelist_for_tri + (nnode_tri * k); 
        memcpy(nodelist_for_tet, nodelist_s, (size_t)sz3int);
        nodelist_for_tet[3] = n0;
        for (i = 0; i < nnode_tet; i++) {
            n = nodelist_for_tet[i];
            r[i] = coord + (dim * n);
        }
        vol_tet[k] = 0.0;
        r1 = r[0];
        r2 = r[1];
        r3 = r[2]; 
        DET3(v, r1[0],r1[1],r1[2], r2[0],r2[1],r2[2], r3[0],r3[1],r3[2])
        vol_tet[k] += v;
        r1 = r[0];
        r2 = r[1];
        r3 = r[3];
        DET3(v, r1[0],r1[1],r1[2], r2[0],r2[1],r2[2], r3[0],r3[1],r3[2])
        vol_tet[k] -= v;
        r1 = r[0];
        r2 = r[2];
        r3 = r[3];
        DET3(v, r1[0],r1[1],r1[2], r2[0],r2[1],r2[2], r3[0],r3[1],r3[2])
        vol_tet[k] += v;
        r1 = r[1];
        r2 = r[2];
        r3 = r[3];
        DET3(v, r1[0],r1[1],r1[2], r2[0],r2[1],r2[2], r3[0],r3[1],r3[2])
        vol_tet[k] -= v;

        if (vol_tet[k] < 0.0) { 
            vol_tet[k] = -vol_tet[k];
        }
        vol_tet[k] /= 6.0;
        c0 = ctr_tet + (dim * k); 
        for (d = 0; d < dim; d++) { 
            c0[d] = 0.0;
        }
        for (i = 0; i < nnode_tet; i++) { 
            n = nodelist_for_tet[i]; 
            c = coord + (dim * n);
            for (d = 0; d < dim; d++) { 
                c0[d] += c[d];
            }
        }
        factor = 1.0/(double)nnode_tet;
        for (d = 0; d < dim; d++) { 
            c0[d] *= factor;
        } 
    }  
    *vol = 0.0;
    for (d = 0; d < dim; d++) { 
        ctr[d] = 0.0;
    }
    for (i = 0; i < ntet; i++) {
        *vol += vol_tet[i];
        c = ctr_tet + (dim * i);
        for (d = 0; d < dim; d++) {  
            ctr[d] += (vol_tet[i] * c[d]); 
        }
    }
    factor = 1.0/(*vol);
    for (d = 0; d < dim; d++) { 
        ctr[d] *= factor;
    }
    free(ctr_tet);
    free(vol_tet);
    free(nodelist_for_tri);

    return;
 } 

int cal_polygon_area(int dim, double **r, int nnode, int d_excluded, double *area)
{
     int    i, i1;

     *area = 0.0;

     if ((dim == 2) || ((dim == 3) && (d_excluded == 2))) {
         for (i = 0; i < nnode; i++) {
             i1 = (i + 1) % nnode;
             *area += (r[i][0] * r[i1][1] - r[i1][0] * r[i][1]);
         }
     }
     else if (dim == 3) {
         if (d_excluded == 0) {
             for (i = 0; i < nnode; i++) {
                 i1 = (i + 1) % nnode;
                 *area += (r[i][1] * r[i1][2] - r[i1][1] * r[i][2]);
             }
         }
         else if (d_excluded == 1) {
            for (i = 0; i < nnode; i++) {
                i1 = (i + 1) % nnode;
                *area += (r[i][2] * r[i1][0] - r[i1][2] * r[i][0]);
            }
         }
     }
     *area *= 0.5;
     if (*area < 0.0) *area = -(*area);

     return 0;
 }

void cal_ctr_vol_2d(int geop, int nnode, int *nodelist, 
                    int nnode_totol, double *coords, double *ctr, double *area)
{ 
//   geop: input. 1 for Cartisan coordinate
//                1 for cylincrical coornidante. 

     int    dim, i, i1, n, n1;
     double third, sixth, twlvth;
     double r0, z0, r1, z1, dz, r02, r12, r0r1, da, rc, r01, zc, darea, factor;
     double *c0, *c1; 

     dim = 2;

     third = 1.0/3.0;
     sixth = 0.5 * third;
     twlvth = 0.5 * sixth;

     ctr[0] = 0.0;
     ctr[1] = 0.0;
     *area  = 0;

     if (geop == 1) {  

         for (i = 0; i < nnode; i++) { 
             i1 = (i + 1) % nnode; 
             n  = nodelist[i]; 
             n1 = nodelist[i1];   
             c0 = coords + (n  * dim);
             c1 = coords + (n1 * dim);
             darea   = c0[0] * c1[1] - c1[0] * c0[1];
             *area  += darea; 
             ctr[0] += ((c0[0] + c1[0]) * darea);
             ctr[1] += ((c0[1] + c1[1]) * darea);
         }
         (*area) *= 0.5;
         if (*area != 0.0) { 
             factor = 1.0/(6.0 * (*area)); 
             ctr[0] *= factor;
             ctr[1] *= factor;
         }
         if (*area < 0.0) { 
             *area = - (*area);
         }  
     }
     else if (geop == 2) {  // cylindrical 
         for (i = 0; i < nnode; i++) { 
             i1 = (i + 1) % nnode;
             n  = nodelist[i];
             n1 = nodelist[i1];
             c0 = coords + (n  * dim);
             c1 = coords + (n1 * dim); 
     
             r0 = c0[0];
             z0 = c0[1];
             r1 = c1[0];
             z1 = c1[1];
     
             dz   = z1 - z0;
             r02  = r0 * r0;
             r12  = r1 * r1; 
             r0r1 = r0 * r1;
     
     //      da = 2 * integral (r dr dz)
     
             darea = dz *(r12 + r0r1 + r02);
             *area += darea;
     
     //      zc = 2 integral (r r dr dz)
     
             rc = dz *(r02 * r0 + r02 * r1 + r0 * r12 + r12 * r1);
             ctr[0] += rc;
        
     //      zc = 2 integral (z r dr dz)
     
             r01 = r0 + r1;
             zc = dz *(r01 * r01 *(z0 + z1) + 2.0 *(r02 * z0 + r12 * z1)); 
             ctr[1] += zc;
         }
         *area  *= third;
         ctr[0] *= sixth;
         ctr[1] *= twlvth;
         if (*area != 0.0) { 
             factor = 1.0/(*area);
             ctr[0] *= factor;
             ctr[1] *= factor;
         } 
         if (*area < 0.0) { 
            *area = -(*area);
         }
     } 

     return;
  } 

int ifnan(double val)
{ 
    if (val == val) { 
        return 0;
    }
    else { 
        return 1;
    }
}


 
double dot(vec_t *a, vec_t *b)
{
	return a->x * b->x + a->y * b->y;
}
 
double cross(vec_t *a, vec_t *b)
{
	return a->x * b->y - a->y * b->x;
}
 
vec_t *vsub(vec_t *a, vec_t *b, vec_t *res)
{
	res->x = a->x - b->x;
	res->y = a->y - b->y;
	return res;
}
 
/* tells if vecT c lies on the left side of directed edge a->b
 * 1 if left, -1 if right, 0 if colinear
 */
int left_of(vec_t *a, vec_t *b, vec_t *c)
{
	vec_t tmp1, tmp2;
	double x;
	vsub(b, a, &tmp1);
	vsub(c, b, &tmp2);
	x = cross(&tmp1, &tmp2);
	return x < 0 ? -1 : x > 0;
}
 
int line_sect(vec_t *x0, vec_t *x1, vec_t *y0, vec_t *y1, vec_t *res)
{
	vec_t dx, dy, d;
	vsub(x1, x0, &dx);
	vsub(y1, y0, &dy);
	vsub(x0, y0, &d);
	/* x0 + a dx = y0 + b dy ->
	   x0 X dx = y0 X dx + b dy X dx ->
	   b = (x0 - y0) X dx / (dy X dx) */
	double dyx = cross(&dy, &dx);
	if (!dyx) return 0;
	dyx = cross(&d, &dx) / dyx;
	if (dyx <= 0 || dyx >= 1) return 0;
 
	res->x = y0->x + dyx * dy.x;
	res->y = y0->y + dyx * dy.y;
	return 1;
}
 
poly_t *poly_new()
{
	return (poly_t *)calloc(1, sizeof(poly_t));
}
 
void poly_free(poly_t *p)
{
	free(p->v);
	free(p);
}
 
void poly_append(poly_t *p, vec_t *v)
{
	if (p->len >= p->alloc) {
		p->alloc *= 2;
		if (!p->alloc) p->alloc = 4;
		p->v = (vec_t *)realloc(p->v, sizeof(vec_t) * p->alloc);
	}
	p->v[p->len++] = *v;
}
 
/* this works only if all of the following are true:
 *   1. poly has no colinear edges;
 *   2. poly has no duplicate vertices;
 *   3. poly has at least three vertices;
 *   4. poly is convex (implying 3).
*/
int poly_winding(poly_t *p)
{
	return left_of(p->v, p->v + 1, p->v + 2);
}
 
void poly_edge_clip(poly_t *sub, vec_t *x0, vec_t *x1, int left, poly_t *res)
{
	int i, side0, side1;
	vec_t tmp;
	vec_t *v0 = sub->v + sub->len - 1, *v1;
	res->len = 0;
 
	side0 = left_of(x0, x1, v0);
	if (side0 != -left) poly_append(res, v0);
 
	for (i = 0; i < sub->len; i++) {
		v1 = sub->v + i;
		side1 = left_of(x0, x1, v1);
		if (side0 + side1 == 0 && side0)
			/* last point and current straddle the edge */
			if (line_sect(x0, x1, v0, v1, &tmp))
				poly_append(res, &tmp);
		if (i == sub->len - 1) break;
		if (side1 != -left) poly_append(res, v1);
		v0 = v1;
		side0 = side1;
	}
}
 
poly_t *poly_clip(poly_t *sub, poly_t *clip)
{
	int i;
	poly_t *p1 = poly_new(), *p2 = poly_new(), *tmp;
 
	int dir = poly_winding(clip);
	poly_edge_clip(sub, clip->v + clip->len - 1, clip->v, dir, p2);
	for (i = 0; i < clip->len - 1; i++) {
		tmp = p2; p2 = p1; p1 = tmp;
		if(p1->len == 0) {
			p2->len = 0;
			break;
		}
		poly_edge_clip(p1, clip->v + i, clip->v + i + 1, dir, p2);
	}
 
	poly_free(p1);
	return p2;
}
 
/////////////////////////////////////////////////////////////////////////
void find_min_max(int dim, int nnode, double *coords, double *xmin, double *xmax)
{
    int n, i; 
    double *c;

    for (i = 0; i < dim; i++) { 
	xmin[i] = coords[i];
	xmax[i] = coords[i]; 
    }
    for (n = 1; n < nnode; n++) { 
	c = coords + (dim * n);
	for (i = 0; i < dim; i++) { 
            if (c[i] < xmin[i]) { 
		xmin[i] = c[i];
            }
	    else if (c[i] > xmax[i]) { 
		xmax[i] = c[i];
            }
        }
    }
    return;
 } 
