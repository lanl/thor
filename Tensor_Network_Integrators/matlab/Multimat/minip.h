#ifndef _MINIP_
#define _MINIP_

#include "globals.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NBDRY 2

enum Bdry_Type {
         bdry_reflected   = 1,
         bdry_transmitted = 2,
	 bdry_periodic    = 3
      };
typedef enum Bdry_Type Bdry_Type;

enum Shape_Type {
     shape_universe     = 0,
     shape_sphere       = 1,
     shape_cylinder     = 2,
     shape_rectangular  = 3, 
     shape_quad         = 4,
     shape_triangle     = 5,
     shape_hex          = 6,
     shape_poly         = 7    // ?  
};
typedef enum Shape_Type Shape_Type;

struct Region_Shape { 
       Shape_Type type;      // one of shapes
       double *parameters;   // parameters describing the shape 
};
typedef struct Region_Shape Region_Shape; 

// for shape_universe, parameters = NULL;
// for shape_sphere,   parameters[0] = radius, parameters[1:dim] = center of radius 
// for shape_cylinder, parameters[0:dim] = radius and center of the one circle 
//                     parameters[dim+1:dim+dim] = radius and center of another circle
// for shape_quad,     parameters[0:7] = (x0,y0; x1,y1; x2,y2; x3,y3) of vertices
// for shape_triangle, parameters[0:5] = (x0,y0; x1,y1; x2,y2) of three vertices
// for shape_hex,      parameters[0:15] = (x0,y0; x1,y1; ...; x15,y15) in the following order
// 
//             7__________6  
//            /.         /|   
//           / .        / |
//          4----------5  |  
//          |  3.......|..2   
//          | .        |  /       
//          |.         | /     
//          0----------1/       
//                      
//     z         
//     |   y  
//     |  /  
//     | / 
//      -------> x   

static const double tiny  = 1.0e-30; 
static const double small = 1.0e-08; 

#ifdef __cplusplus
}
#endif

#endif

