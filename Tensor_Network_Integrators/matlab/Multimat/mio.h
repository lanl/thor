/* ======================= copyright begin ========================  */
/* Copyright (C) 2010 Los Alamos National Security, LLC.             */
/* All rights Reserved.  See Copyright Notice File.                  */
/* Export Controlled Information                                     */
/* ======================== copyright end =========================  */

#ifndef MIO_H
#define MIO_H


#ifdef __cplusplus
extern "C" {
#endif

#define MAX_COMP_RANK  3
#define MAX_MATS       100
#define MAX_MESHES_IN_UMESH 16

/**< max rank for a tensor, it may be
     changed to a larger number for the maximum
     value of rank in a tensor.     
*/

/**
\mainpage MeshIO Interface User's Guide
\author William W. Dai, HPC-4, Los Alamos National Laboratory
\version V-2007-July
\section Copyright Copyright

  Copyright (c) 2006, The Regents of the University of California.
  All rights reserved.

  Copyright (2006). The Regents of the University of California.
  This software was produced under U.S. Government contract
  W-7405-ENG-36 for Los Alamos National Laboratory (LANL),
  which is operated by the University of California for the U.S.
  Department of Energy. The U.S. Government has rights to use,
  reproduce, and distribute this software. NEITHER THE GOVERNMENT
  NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
  ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

\section acknowledgements Acknowledgements
 This research and development are funded by the Department of Energy's ASCI program.
 Los Alamos National Laboratory is operated by the University of California
 for the National Nuclear Security Administration of the United States
 Department of Energy under contract W-7405-ENG-36.

\section abstract Abstract
 This MIO Interface User's Guide provides the information needed to write
 and read data on parallel or serial computer platforms, and to query files
 through MIO library interface functions. Example codes are included in this
 User's Guide for the usage.

\section overview Overview

 The MIO library is an IO library built directly on MPI and MPI-IO.
 The goal of the MIO library is to provide simulation codes
 an efficient IO capability on parallel or serial computer platforms
 for various data structures defined in applications. The
 data structures include one- and multi-dimensional arrays, structured meshes
 and variables, unstructured mesh and variables, and adaptive mesh refinements. 
 The files created through the MIO library may contain hierarchical structures and
 descriptions of arrays, meshes, variables, and parallel computer environments.
 Since the files are self-describing,
 users may read any part of an array, or mesh, or variable,
 without knowledge of the file.
 The goal of the MIO Interface User's Guide is to introduce the main functionality
 of the library through a number of C/C++ callable functions.

 Although the MIO library doesn't depend on the hierarchical structure,
 a file written through the MIO library may have the structure.
 The hierarchical structure is similar to a Unix file system.
 In a Unix file system, directories can
 contain other files. The top directory is referred to as the
 "root" directory and has the "/" designation. Each file has
 certain attributes, e.g., read/write permissions, etc.
 In the file written through the MIO library, the object corresponding to a
 directory is a group, the object corresponding to a Unix file is an array.
 or mesh, or variable. A group is a container for other objects.
 In this users guide, groups, arrays, meshes, and variables will be called objects.

 Attributes are associated with objects. Unlike groups, an attribute can only
 be written once the object the attribute will be associated with has been created.
 Likewise, an attribute are read by referring to the location of the object to which it is
 associated. Attributes allow additional information
 to be associated with an object. With attributes, users may add more descriptions
 to any object. 

 Example codes are in the directory ../examples.

 For all the C/C++ interface functions in this users' guide,
 an id of an object is unique,
 and it is non-negative. A function will return 0
 if no errors were encountered; otherwise, -1 will be returned.
 All constants shown in the User's Guide are included in the file

 - mio.h
*/  

/**
\defgroup file Functions for Opening and Closing Files

 This section covers the functions to open and
 close files.

*/

/**
\ingroup file Functions for Opening and Closing Files

 mio_IO_Mode contains the values for possible IO mechanisms.

*/

enum mio_IO_Mode{
/** - MPI_File_write_all */
 mio_collective = 1,
/** -  MPI_File_write_at */
 mio_independent = 2,
/** - fopen, without MPI-IO  */
 mio_use_fopen = 3,
/** - open, without MPI-IO  */
 mio_use_open = 4
};
typedef enum mio_IO_Mode mio_IO_Mode;

/**
\ingroup file Functions for Opening and Closing Files

 mio_File_Mode contains the values of open modes for files.

*/

enum mio_File_Mode{
/** - file create */
 mio_file_create = 1,
/** - file read only */
 mio_file_read_only = 2,
/** - file read and write, not implemented yet */
 mio_file_read_write = 3
};
typedef enum mio_File_Mode mio_File_Mode;

/**
\ingroup file Functions for Opening and Closing Files

 Given a relative path or full path for a file,
 mio_create_file will create a new file or open an existing file. See example_open.c.
 This function
 also automatically creates or opens the root group. The output, file_id,
 is also the id of the root group.

 - file_name: an input, the relative path or full path for a file,
              such as restart/file1.

 - file_mode: an input, one of the three values of mio_File_Mode.

 - io_mode: an input, one of the three values of mio_IO_Mode.

 - file_id: an output, the id of the file just created or opened.
*/
int mio_create_file(char *file_name, mio_File_Mode file_mode, mio_IO_Mode io_mode, int *file_id);

/**
\ingroup file Functions for Opening and Closing Files

 mio_open_file is the old version of mio_create_file with io_mode = mio_independent.
*/
int mio_open_file(char *file_name, mio_File_Mode file_mode, int *file_id);

/**
\ingroup file Functions for Opening and Closing Files

 mio_close_file will close the file associated with a given id, and return
 the size of the file just written in bytes. See example_open.c.

 - file_id: an input, the id of the file which was created or opened before.

*/
long long mio_close_file(int file_id);

/**
\ingroup file Functions for Opening and Closing Files

 mio_correct_file returns 1 if an existing file was written through mio. Otherwise,
 it returns 0. 
       
 - filename: an input, a relative path or full path of an existing file.
*/
int mio_correct_file(char *filename);

/**
\ingroup file Functions for Opening and Closing Files

 mio_resilient is only for writing data. By default, mio library doesn't
 do any thing for resilience. But if this function is called with a non-zero
 argument, the following files to be created and data to be written will be
 resilient in the sense that all the data previously written could be read
 should io crashes before a file is closed.

 - ifresilient: an input, 0 or 1.
*/

int mio_resilient(int ifresilient);

/**
\ingroup file Functions for Opening and Closing Files


 mio_backup_file is used to setup the filename of backup file only when
 mio library is used for resilience. If this function is not called, the
 backup file will be the filename of the file to be created with the
 extension "backup".

 - backup_file: the filename of the backup file.
*/

 int mio_backup_file(char *backup_file);

/**
 \ingroup file Functions for Opening and Closing Files
   
   By default, mio reads all the buffered data through buffering, i.e, 
   bio reads the whole buffer and keeps the buffer in memory for any
   array. mio_buffer_mode resets the mode to read any arrays in buffers.
   This function has no effects in writing. 
 
   - flg_for_buffered_read : an input
                             0          for non-buffered read
                             otherwise, for buffered read.
*/

int mio_buffer_mode(int flg_for_buffered_read);

int mio_init_compression(int fileid, char *name, int collective, long long nbytes);

int mio_finalize_compression(int fileid);

/**
\defgroup io Functions for Writing Data

  This section covers the functions to open an group,
  to write and read arrays, structured meshes,
  unstructured meshes, AMR meshes, variables, and attributes.

*/

/**
\ingroup io Functions for Writing Data

 mio_Object_Type contains the values of mio_Object_Type.

*/

enum mio_Object_Type {
      mio_grp    = 1,
      mio_array  = 2,
      mio_smesh  = 3,
      mio_umesh  = 4,
      mio_smesh_cell_amr  = 5,
      mio_umesh_cell_amr  = 6,
      mio_smesh_block_amr = 7,

      mio_mesh_var = 8,
      mio_coord    = 9,
      mio_array_struct = 10,
      mio_attr         = 11,

      mio_mat_group    = 12,
      mio_matvar_group = 13,
      mio_matvar       = 14,

      mio_obj_type_invalid = 128 
};
typedef enum mio_Object_Type mio_Object_Type;

/**
\ingroup io Functions for Writing Data

 mio_Group defines the object, group.

*/
struct mio_Group {
       int  id; 
       char *name;  
      };
typedef struct mio_Group mio_Group;

/**
\ingroup io Functions for Writing Data

 mio_Data_Type contains the values of data types.

*/

enum mio_Data_Type { 
     /** char */ 
     mio_char             = 1,     
     /** double */ 
     mio_double           = 2,     
     /** float  */ 
     mio_float            = 3,     
     /** int    */
     mio_int              = 4,     
     /** long   */
     mio_long             = 5,     
     /** long double */
     mio_long_double      = 6,     
     /** long long   */
     mio_long_long        = 7,     
     /* invalid      */ 
     mio_datatype_invalid = 16
};
typedef enum mio_Data_Type mio_Data_Type;

/**
\ingroup io Functions for Writing Data

 mio_Attr defines the structure for attributes.
*/

struct mio_Attr {
       char *name;
       void *values;
       int  count;
       mio_Data_Type datatype;
      };
typedef struct mio_Attr mio_Attr;

/**
\ingroup io Functions for Writing Data

 mio_Index_Domain defines the structure for a domain of a structured mesh. 
*/

struct mio_Index_Domain {
       long long offsets[3];
       long long sizes[3];
       long long nbdyl[3];
       long long nbdyr[3];
      };  
typedef struct mio_Index_Domain mio_Index_Domain;

/**
\ingroup io Functions for Writing Data

 mio_Elem_Domain defines a part of an unstructured mesh in terms of global element ids.
*/

struct mio_Elem_Domain {
       long long elem_offset;
       long long elem_size;
     }; 
typedef struct mio_Elem_Domain mio_Elem_Domain;

/**
\ingroup io Functions for Writing Data

 mio_Space_Relation defines a two possible relations of two space domains,  
 which will be used in mio_query. 
*/

enum mio_Space_Relation {
     mio_touch = 1,
     mio_cover = 2
};    
typedef enum mio_Space_Relation mio_Space_Relation;

/**
\ingroup io Functions for Writing Data

 mio_Coord_Domain defines a part of an unstructured mesh in terms of coordinates. 
*/
struct mio_Coord_Domain {
       int dims; 
       double coordmin[3];
       double coordmax[3];
       mio_Space_Relation criterion;
     };
typedef struct mio_Coord_Domain mio_Coord_Domain;

/**
\ingroup io Functions for Writing Data

 mio_Mat_Domain defines a part of mesh in terms of material properties.
 Not implemented yet. 
*/
struct mio_Mat_Domain {
      char *mat_name;
      double vmin;
      double vmax;
     };
typedef struct mio_Mat_Domain mio_Mat_Domain;

/**
\ingroup io Functions for Writing Data

 mio_PE_Domain defines a part of mesh in terms of the rank of a PE.
*/
struct mio_PE_Domain {
     int rank;
    };
typedef struct mio_PE_Domain mio_PE_Domain; 

/**
\ingroup io Functions for Writing Data

 mio_Partition_Method defines a few possible methods for domain partition.
 Not implemented yet. 
*/
enum mio_Partition_Method {
     mio_pe_partition     = 600,
     mio_elem_partition   = 601,
     mio_coord0_partition = 602,
     mio_coord1_partition = 603,
     mio_coord2_partition = 604,
     mio_coord_partition  = 605,
     mio_mat_partition    = 606,

     mio_partition_method_invalid = 699
};
typedef enum mio_Partition_Method mio_Partition_Method;

/**
\ingroup io Functions for Writing Data

 mio_Partition_Domain defines a part of an unstructured mesh in terms of domain partition.
 Not implemented yet.
*/
struct mio_Partition_Domain {
     int mycpu;
     int ncpus;
     mio_Partition_Method method;
};
typedef struct mio_Partition_Domain mio_Partition_Domain;

/**
\ingroup io Functions for Writing Data

 mio_Domain_Type and mio_Domain define a few possible methods to define
 a part of an unstructured mesh. 
*/

enum mio_Domain_Type {
     mio_index_domain  = 301,
     mio_elem_domain   = 302,
     mio_coord_domain  = 303,
     mio_pe_domain     = 304,
     mio_mat_domain    = 305,

     mio_partition_domain = 305,

     mio_domain_invalid = 400 
};
typedef enum mio_Domain_Type mio_Domain_Type;


struct mio_Domain {
       mio_Domain_Type type; 
       void *obj;
     };
typedef struct mio_Domain mio_Domain;

/**
\ingroup io Functions for Writing Data

 mio_Array_Struct defines the structure of one- or multi-dimensional arrays.
*/

struct mio_Array_Struct {
           int id;
           int dims;           /**< The dimension 0 is the one changing most
                                    slowly, and the dimension (dims-1) is
                                    the one changing most fast.           */
 
           long long *sizes;   /**< the sizes of an array on the current
                                    PE in dims dimensions, which exclude
                                    any possible ghost points.            */
   
           long long *nbdyl;   /**< an array with dims elements for
                                    the sizes of ghost (fake) points at
                                    the lower ends of dims dimensions.    */

           long long *nbdyr;   /**< an array with dims elements for
                                    the sizes of ghost (fake) points at
                                    the higher ends of dims dimensions.   */

           long long *offsets; /**< offsets of the part of the array on the
                                    current PE in dims dimensions.        */

           long long *gsizes;  /**< the total sizes of an array in dims
                                    dimensions, excluding ghost points.   */
   
           long long *psizes;  /**< an array with dims elements, which is
                                    for the sizes of PE configuration in 
                                    dims dimensions.
                                                                          */
           long long *plist;   /**< an array with its size
                                    psizes[0] * psizes[1] * ... * psizes[dims-1]
                                    for the order of PEs in the
                                    PE configuration. plist sweeps the
                                    dimension dims-1 first, and sweeps the
                                    dimension 0 last. By default, plist[i] = i,
                                    for i = 0, 1, ..., npes-1.
\verbatim 
           ----------------------------------------
    dim 1  |            |            |            |
           | plist(1)   | plist(3)   | plist(5)   |
      ^    |            |            |            |
      |    ----------------------------------------
      |    |            |            |            |
      |    | plist(0)   | plist(2)   | plist(4)   |
           |            |            |            |
           ----------------------------------------
                        ----> dim 0
\endverbatim
                                                                          */
           int npes;               /**< The number of PEs which wrote arrays 
                                        with this array_struct. This member
                                        is ignored in mio_write.     
                                                                    */
           int multiple;           /**< the number of values at each point.
                                    By default, it is set to 1.     */
};   
typedef struct mio_Array_Struct mio_Array_Struct;

/**
\ingroup io Functions for Writing Data
 mio_Array defines the structure of arrays.
*/

struct mio_Array {
       int              id;
       char             *name; 
       mio_Array_Struct array_struct; 
       mio_Data_Type    datatype;
       void             *buffer;
      };
typedef struct mio_Array mio_Array;

/**
\ingroup io Functions for Writing Data
 mio_Structured_Mesh defines the structure of structured meshes. 
*/

struct mio_Structured_Mesh {
       int  id;     
       char *name;  
       int  dims;           /**< the dimensionality of a structured mesh */
       long long sizes[3];  /**< sizes of mesh elements in dims dimension (NOT
                                 sizes of grid points). Dimension 0 is the
                                 dimension changing most slowly. It does
                                 not include any possible ghost mesh elements.*/

       long long offsets[3];/**< offsets of mesh elements in dims dimension for
                                 the current PE. It doesn't include
                                 any possible ghost mesh elements.            */

       long long gsizes[3]; /**< total sizes of mesh elements in dims dimensions. 
                                 It doesn't include any possible ghost mesh elements.  */

       long long psizes[3]; /**< sizes of PE-configuration in dims dimensions.
                                 It will be used when offsets are not provided
                                 in writing a mesh.  By default,
                                 psizes[0] = npes, psizes[i] = 1 for i != 0.
                                                                         */
       long long *plist;    /**< This is an array to indicate the order of PEs
                                 when offsets are not
                                 provided for writing a mesh.
\verbatim
           ----------------------------------------
    dim 1  |            |            |            |
           | plist(1)   | plist(3)   | plist(5)   |
      ^    |            |            |            |
      |    ----------------------------------------
      |    |            |            |            |
      |    | plist(0)   | plist(2)   | plist(4)   |
           |            |            |            |
           ----------------------------------------
                        ----> dim 0
\endverbatim
                                                                             */

       int npes;               /**< The number of PEs which wrote the structured
                                    mesh. This member is ignored in mio_write.  
                                                                    */

       long long nbdyl[3];  /**< the number of ghost mesh elements at the lower
                                 ends of dims dimensions.                    */
       long long nbdyr[3];  /**< the number of ghost mesh elements at the
                                 higher ends of dims dimensions.        */

/*  for coordinate points being element centered or not  */

       int element_centered;/**< It should be 1 if a coordinate point,
                                 for example in a 2D mesh,
                                 (coord[0][k],coord[1][j]), is
                                 a center of a mesh element, and it should be
                                 0 if (coord[0][k],coord[1][j])
                                 is a grid point.               */

/* for a set of uniform grids at any dimension  */

       double dcoord[3];    /**< If dcoord[i] > 0 (i = 0, 1, 2), dimension
                                 i will be considered uniform in dimension i
                                 in the current PE.                      */

       double coordmin[3];  /**< coordmin[i] (i = 0, 1, 2) is the minimum
                                 value of coordinates of real
                                 mesh elements in the current PE in dimension i.
                                 coordmin[i] will be used for the following
                                 two situations: (a) for a uniform mesh, i.e.,
                                 dcoord[i] > 0, (b) element_centered = 1.
                                                                          */
/* for a set of nonuniform grids at any dimension  */

       void *coord[3];      /**< coord[i] (i = 0, 1, 2) should be valid
                                 if dcoord[i] < 0. If element_centered = 0,
                                 the array coord[i] is
                                 (sizes[i] + nbdyl[i] + nbdyr[i] + 1)
                                 long, otherwise it is
                                 (sizes[i] + nbdyl[i] + nbdyr[i]) long.  */

       mio_Data_Type datatype; /**<  datatype of coord[i], i = 0, 1, 2.
                                  Currently, only mio_float, mio_double,
                                  and mio_long_double are supported.     */
     };
typedef struct mio_Structured_Mesh mio_Structured_Mesh;

/**
\ingroup io Functions for Writing Data

 mio_Structured_Cell_AMR defines the structure of structured meshes with
 adaptive mesh refinements.
*/

struct mio_Structured_Cell_AMR {
       int   id; 
       char  *name; 

       int  dims;  

       long long offset; 
       long long gsize; 

       long long size; 

       long long fsize;

       long long bsize;

       long long *plist; 
 
       int npes;               /**< The number of PEs which wrote the mesh.  
                                    This member is ignored in mio_write.
                                                                    */
       mio_Data_Type datatype_ids;  
       int  idmin;  
       void *fakeelems; 
       void *masterelems; 
       void *bdryelems; 

       void *gid;
       void *nelems_connected_to_elem; /**< an array whose length is number of
                                       mesh elements for the number of mesh elements
                                       each mesh element is connected to.
                                       In zone-meshes, zone is the element,
                                       in face-meshes, face is the element,
                                       and in edge-meshes, edge is the element.
                                       If this array is not
                                       provided, don't set this array and
                                       the following array.       */

       void *elems_connected_to_elem;  /**< an array whose length is the sum
                                        of nelems_connected_to_elem for the
                                        list of mesh elements connected by each
                                        element.                   */

       int element_centered;    /**< It should be 1 if a coordinate point,
                                 for example in a 2D mesh,
                                 (coord[0][k],coord[1][j]) is
                                 a center of a mesh element, and it should be
                                 0 if (coord[0][k],coord[1][j])
                                 is a grid point.               */
     
       mio_Data_Type datatype_coord;  
       void *coord[3];
       void *dcoord[3];

       long long size_gid;
       long long size_masterelems;
       long long size_elems_connected_to_elem;
      };
typedef struct mio_Structured_Cell_AMR mio_Structured_Cell_AMR;

/**
\ingroup io Functions for Writing Data

 mio_Unstructured_Mesh_Type defines several types of unstructured meshes. 
*/

enum mio_Unstructured_Mesh_Type {
                    /** - points               */ 
                    mio_point    = 401,        
                    /** - bars with two nodes  */
                    mio_bar      = 402,
                    /** - triangles with three nodes */
                    mio_triangle = 403,
                    /** - quadrangles with four nodes */
                    mio_quad     = 404,
                    /** - tetrahedron with four nodes */
                    mio_tet      = 405,
                    /** - pyramids with five nodes       */
                    mio_pyramid  = 406,
                    /** - hexahedron with eight nodes */
                    mio_hex      = 407,
                    /** - wedges with six nodes       */
                    mio_wedge    = 408,
                    /** - mixture of above, not implemented yet */
                    /** - pentagon with five nodes    */
                    mio_pentagon = 411,
                    /** - pentagon prism with 10 nodes   */
                    mio_pentagon_prism = 412,
                    /** - mixture of above, not implemented yet */
                    mio_mixed_elements = 409,
                    /** - general polyhedrons  */
                    mio_general_mesh = 410,
                    /** - invalid unstructured mesh type */
                    mio_meshtype_invalid = 500 
                   };
typedef enum mio_Unstructured_Mesh_Type mio_Unstructured_Mesh_Type;

/**
\ingroup io Functions for Writing Data

 mio_Unstructured_Mesh_Size defines a structure to hold the lengths of the arrays
 in a unstructured meshes.
*/

struct mio_Unstructured_Mesh_Size {

  long long size_nodelist_for_zone; /**< the size of the array nodelist_for_zone
                                         for a part of mesh or whole mesh.        */
  long long totalsize_nodelist_for_zone;  
  long long size_edgelist_for_zone; /**< the size of the array edgelist_for_zone
                                         for a part of mesh or whole mesh.        */
  long long totalsize_edgelist_for_zone;
  long long size_facelist_for_zone; /**< the size of the array facelist_for_zone
                                         for a part of mesh or whole mesh.        */
  long long totalsize_facelist_for_zone;
  long long size_nodelist_for_face; /**< the size of the array nodelist_for_face
                                         for a part of mesh or whole mesh.        */
  long long totalsize_nodelist_for_face; 
  long long size_edgelist_for_face; /**< the size of the array edgelist_for_face
                                         for a part of mesh or whole mesh.        */
  long long totalsize_edgelist_for_face;
  long long size_elems_connected_to_elem;
                                    /**< The size of the array
                                         elems_connected_to_elem for a part of mesh 
                                         or whole mesh. If it is zero, there are no
                                         arrays, nelems_connected_to_elem and
                                         elems_connected_to_elem, defined in the
                                         mesh or the part of the mesh.        */
  long long totalsize_elems_connected_to_elem;
  long long size_elems_connected_to_node;
                                    /**< The size of the array
                                         elems_connected_to_node for a part of mesh
                                         or whole mesh. If it is zero, there are no
                                         arrays, nelems_connected_to_node and
                                         elems_connected_to_node, defined in the
                                         mesh or the part of the mesh.        */
  long long totalsize_elems_connected_to_node;
  long long size_masterzones;      /**< The size of the array masterzones for a
                                        part of mesh.                           */
  long long totalsize_masterzones; 
  long long size_masterfaces;      /**< The size of the array masterfaces for a 
                                        part of mesh.                           */
  long long totalsize_masterfaces; 
  long long size_masteredges;      /**< The size of the array masteredges for a
                                        part of mesh.                           */
  long long totalsize_masteredges; 
  long long size_masternodes;      /**< The size of the array masternodes for a 
                                        part of mesh.                           */
  long long totalsize_masternodes;

  int gids_found[4];                /**< if gids_found[i] != 0, gids[i] is found
                                         in the mesh, i = 0, 1, 2, 3. Currently,
                                         only gids[3], i.e., gids of nodes, is
                                         supported.                            */
 };
typedef struct mio_Unstructured_Mesh_Size mio_Unstructured_Mesh_Size;

/**
\ingroup io Functions for Writing Data
 mio_Coord defines the structure of coordinates used in unstructured meshes.
*/

struct mio_Coord {
            int id;
            int mesh_ids[MAX_MESHES_IN_UMESH];
                                 /**< a list of ids of a set of unstructured meshes
                                      the set of coordinates associated with. */
            int num_meshes;      /**< the size of the array mesh_ids.         */
            void *coord[3];      /**< coord[i] (i = 0, 1, 2) should be valid
                                    if dcoord[i] < 0. If element_centered = 0,
                                    the array coord[i] is
                                    (sizes[i] + nbdyl[i] + nbdyr[i] + 1)
                                    long, otherwise it is
                                    (sizes[i] + nbdyl[i] + nbdyr[i]) long.  */

            mio_Data_Type datatype; /**<  datatype of coord[i], i = 0, 1, 2.
                                     Currently, only mio_float, mio_double,
                                     and mio_long_double are supported.     */
           };
typedef struct mio_Coord mio_Coord;

/**
\ingroup io Functions for Writing Data
 mio_Unstructured_Mesh defines a structure for unstructured meshes. 
*/

struct mio_Unstructured_Mesh {
       int   id; 
       char  *name; 
       mio_Unstructured_Mesh_Type type; /**< type of unstructured mesh */

       int dims;  /**< Dimensionality of the mesh. It should be 3
                       if coord[0], coord[1], and coord[2] (below) are used to
                       specify the coordinates of nodes, and it should
                       be 2 if only coord[0] and coord[1] are used, and it
                       should be 1 if only coord[0] is used.        */

       long long offsets[4]; /**< offsets[i] with i from 0 to 3 are
                                  offsets of zones, faces, edges and
                                  nodes in the current PE.
                                  Currently, only the offset for mesh
                                  elements has been implemented, i.e.,
                                  offsets[0] for a zone-mesh,
                                  offsets[1] for a face-mesh, offsets[2]
                                  for an edge_mesh, and offsets[3] for
                                  a node mesh.               */

       long long gsizes[4]; /**< gsizes[i] with i from 0 to 3 are
                                 the total number of zones, faces,
                                 edges and nodes. Currently, only the
                                 total number of mesh elements are
                                 implemented, i.e., gsizes[0] for a zone-mesh,
                                 offsets[1] for a face-mesh, offsets[2]
                                 for an edge-mesh, and gsizes[3] for
                                 a node-mesh.                 */
                                                                
       long long sizes[4]; /**< sizes[i] with i from 0 to 3 are
                                the numbers of zones, faces, edges
                                and nodes in the current PE.       */

       long long fsizes[4]; /**< not tested. fsizes[i] with i from 0 to 3 are
                                 the numbers of ghost zones, faces,
                                 edges and nodes in the current PE.
                                 If fsizes[0] is set to non-zero,
                                 the array, fakezones, must be fsizes[0]
                                 long. If fsizes[1] is non-zero,
                                 the array, fakefaces, must be fsizes[1] long,
                                 and the same logic holds for
                                 (fsizes[2], fakeedges)
                                 and (fsizes[3], fakenodes).        */

       long long ssizes[4]; /**< not tested. ssizes[i] with i from 1 to 3 are
                                 the numbers of slip faces, edges, and
                                 nodes in the current PE. If ssizes[1]
                                 is set to non-zero, the array, slipfaces, 
                                 must be ssizes[1] long. If ssizes[2]
                                 non-zero, the array, slipedges, must be ssizes[2]
                                 long, and if ssizes[3] is non-zero,
                                 the array, slipnodes, must be ssizes[3] long.
                                 ssizes[0] is ignored.              */

       long long bsizes[4]; /**< not tested. bsizes[i] with i from 1 to 3 are
                                 the numbers of faces, edges, and nodes
                                 on a boundary. If bsizes[1] is set to
                                 non-zero, the array, bdryfaces, must be
                                 bsizes[1] long. If bsizes[2] is non-zero,
                                 the array, bdryedges, must be bsizes[2]
                                 long, and if bsizes[3] is non-zero,
                                 the array, bdrynodes, must be bsizes[3] long.
                                 bsizes[0] is ignored.              */

       mio_Data_Type datatype;  /**< datatype of all void integer
                                         arrays used in the structure.
                                                                    */
       long long *plist; /**< plist is used for the order of elements
                              when a mesh is written and when
                              the offsets of elements on PEs are not given.
                              For zone-meshes, the elements are zones; for
                              face-meshes, the elements are faces; and for
                              edge-meshes, the elements are edges.
                                                                    */

       int npes;         /**< The number of PEs which wrote the mesh.
                              This member is ignored in mio_write.
                                                                    */
       int idmin; /**< 1 for one-based ids, and 0 for zero-based ids.
                               This should be the same for all PEs.
                               For example, in a mesh with 0-based ids, a
                               node with id 1 has coordinates
                               (coord[0][1], coord[1][1]). But, a mesh with
                               1-based ids, a node with id 1 has coordinates
                               (coord[0][0], coord[1][0]).
                               If users don't explicitly specify idmin in
                               a call of mio_mesh_write, idmin is the
                               one in the most recent call of mio_mesh_write. */

       int order_for_nodelist; /**< This field will be used zone-elements are
                                    directly made from nodes  
                                    to indicate the order of nodes for
                                    each mesh element. Only three values
                                    are allowed: 1 for the order of Ensight
                                    (right hand rule), -1 for the
                                    order used in GMV (left hand rule),
                                    0 for others or unspecified. By default,
                                    -1 will be used. If users don't explicitly
                                    specify it in a call of mio_mesh_write,
                                    the one in the most recent call of
                                    mio_mesh_write will be used.        */ 

       void *num_faces_for_zone; /**< an array of sizes[0] long
                                  for the number of faces of each zone
                                  for general polyhedrons. This field could be
                                  used only when the type of mesh is
                                  mio_general_mesh or mio_mixed_elements.  */

       void *num_edges_for_zone; /**< an array of sizes[0] long
                                  for the number of edges for each zone
                                  for general polyhedrons. This field could be used
                                  only when the type of mesh is mio_general_mesh
                                  or mio_mixed_elements. */

       void *num_nodes_for_zone; /**< an array of sizes[0] long
                                  for the number of nodes for each zone
                                  for general polyhedrons. This field could be used
                                  only when the type of mesh is mio_general_mesh
                                  or mio_mixed_elements.        */

       void *num_edges_for_face; /**< an array of sizes[1] long
                                      for the number of edges
                                      for each face.                 */

       void *num_nodes_for_face; /**< an array of sizes[1] long
                                      for the number of nodes
                                      for each face.                 */

       void *facelist_for_zone;  /**< an array of the face list for
                                  zones. The length of the array is
                                  the sum of the numbers of faces of each zone.
                                  This field is ignored if no faces are involved
                                  in a mesh.                 */

       void *edgelist_for_zone;  /**< an array of the edge list for
                                  zones. The length of the array is the sum
                                  of the numbers of edges of each zone.
                                  This field is ignored if no edges are involved
                                  in a mesh.                     */

       void *nodelist_for_zone;  /**< an array of the node list for
                                  zones. The length of the array is the sum of
                                  the numbers of nodes of each zone.   */

       void *edgelist_for_face;  /**< an array of the edge list for
                                  faces. The length of the array is the sum of the
                                  numbers of edges of each face.       */

       void *nodelist_for_face;  /**< an array of the node list for
                                  faces. The length of the array is the sum of the
                                  numbers of nodes of each face.       */

       void *nodelist_for_edge;  /**< an array of the node list for
                                  edges. The length of the array is the sum of
                                  the numbers of edges.                */

       void *gid_zone; /**< an array of length sizes[0] for global ids for zones.
                            not implemented yet */

       void *gid_face; /**< an array of length sizes[1] for global ids for faces.
                            not implemented yet */
       void *gid_edge; /**< an array of length sizes[2] for global ids of edges,
                            not implemented yet */
       void *gid_node; /**< an array of length sizes[3] for global ids of nodes */

       void *fakezones; /**< not tested yet. an array of length fsizes[0] for the
                             list of ghost zones.                     */
       void *masterzones; /**< not tested yet. 
                          an array of length 2 * fsizes[0]. The first fsizes[0]
                          values are the ids of master zones associated with the
                          ghost zones, and the next fsizes[0] values are the PE 
                          ranks associated with the master zones.      */

       void *fakefaces; /**< not tested yet.
                             an array of length fsizes[1] for the list of ghost faces.            
                         */
       void *masterfaces; /**< not tested yet.
                           an array of length 2 * fsizes[1]. The first fsizes[1]
                           values are the ids of master faces associated with the
                           ghost faces, and the next fsizes[1] values are the PE
                           ranks associated with the master faces.      */

       void *fakeedges; /**< an array of length fsizes[2] for the
                             list of ghost edges.                     */
       void *masteredges; /**< not tested yet.
                           an array of length 2 * fsizes[2]. The first fsizes[2]
                           values are the ids of master edges associated with the
                           ghost edges, and the next fsizes[2] values are the PE 
                           ranks associated with the master edges.      */

       void *fakenodes; /**< an array of length fsizes[3] for the
                             list of ghost nodes.                     */
       void *masternodes; /**< not tested yet.
                           an array of length 2 * fsizes[3]. The first fsizes[3]
                           values are the ids of master nodes associated with the
                           ghost nodes, and the next fsizes[3] values are the PE 
                           ranks associated with the master nodes.      */

       void *slipzones; /**< not tested yet.
                             an array of length ssizes[0] for a list of
                             zones specified.                */
       void *slipfaces; /**< not tested yet.
                             an array of length ssizes[1] for the
                             list of slip faces.             */
       void *slipedges; /**< not tested yet.
                             an array of length ssizes[2] for the
                             list of slip edges.             */
       void *slipnodes; /**< not tested yet.
                             an array of length ssizes[3] for the
                             list of slip nodes.             */

       void *bdryzones; /**< not tested yet.
                             an array of length bsizes[0] for a list
                             of zones specified. */ 
       void *bdryfaces; /**< not tested yet.
                             an array of length bsizes[1] for the list
                             of faces on the boundary.                  */
       void *bdryedges; /**< not tested yet.
                             an array of length bsizes[2] for the list
                             of edges on the boundary.                  */
       void *bdrynodes; /**< not tested yet.
                             an array of length bsizes[3] for the list
                             of nodes on the boundary.                  */

       void *nelems_connected_to_elem; /**< not tested yet.
                                       an array, whose length is number of
                                       mesh elements, for the number of mesh elements
                                       each mesh element connected to.
                                       In zone-meshes, zone is the element,
                                       in face-meshes, face is the element,
                                       and 
                                       in edge-meshes, edge is the element.
                                       If this array is not
                                       provided, don't set this array and
                                       the following array.       */

       void *elems_connected_to_elem;  /**< not tested yet.
                                        an array whose length is the sum
                                        of nelems_connected_to_elem for the
                                        list of mesh elements connected by each
                                        element.                   */

       void *nelems_connected_to_node; /**< not tested yet.
                                        an array of length sizes[3]
                                        for the number of mesh elements
                                        each node connected to.
                                        If this array is not
                                        provided, don't set this array and
                                        the following array.       */

       void *elems_connected_to_node;  /**< not tested yet.
                                        an array whose length is the sum
                                        of nelems_connected_to_node for the
                                        list of elements connected by each
                                        node.                   */

       mio_Unstructured_Mesh_Type *types; /**< an array whose length is the number
                                          of mesh elements.
                                          In zone-meshes, zone is the element,
                                          in face-meshes, face is the element,
                                          and in edge-meshes, edge is the element,
                                          This array is used if
                                          and only if mesh type is
                                          mio_mixed_elements.
                                          The values of the array are not
                                          allowed to be mio_mixed_elements
                                          and mio_general_mesh. The values
                                          are either all zone-types (such as
                                          mio_tet, mio_hex, mio_wedge,
                                          mio_pentagon_prism) or all
                                          face-types (such as mio_triangle,
                                          mio_quad, mio_pentagon).
                                          This field has not been implemented yet.*/
       mio_Unstructured_Mesh_Size msize; 
       mio_Coord coord;
      };
typedef struct mio_Unstructured_Mesh mio_Unstructured_Mesh;

/**
\ingroup io Functions for Writing Data
 mio_Mesh_Var_Type defines types of variables defined on a mesh. 
*/

enum mio_Mesh_Var_Type {
     /** - variables defined on 3D zones */
     mio_zone = 501,
     /** - variables defined on faces, including faces in a 2D or 3D space */
     mio_face = 502,
     /** - variables defined on edges */
     mio_edge = 503,
     /** - variables defined on nodes    */
     mio_node = 504,
     /** - invalid mesh variable type    */
     mio_meshvartype_invalid  = 600 
};
typedef enum mio_Mesh_Var_Type mio_Mesh_Var_Type;

/**
\ingroup io Functions for Writing Data
 mio_Var_Comps defines a scalar, or each component of a vector or tensor.
*/

struct mio_Var_Comps {
       int which_comp[MAX_COMP_RANK]; /**< This field is used to indicate a specific
                                     component of a vector or tensor.
                                     The values of the array is
                                     zero-based.
                                     For a scalar this field is ignored.
                                     For a vector, which_comp[0] is used
                                     to specify a component. For tensors,
                                     for example,
                                     (which_comp[0], which_comp[1]) =
                                     (0,1) refers to the component
                                     T_01 of a tensor T_ij.          */

       void *buffer; /**< It is the starting address of the values for a 
                          scalar variable, or a component of a variable.
                          The size of buffer should be consistent
                          with the associated mesh and the type of
                          the variable. For example, 
                          for an element variable on a structured mesh,
                          the expected size
                          should be product of 
                          (sizes[i] + nbdyl[i] + nbdyr[i])
                          (i = 0, 1, 2) of the associated 3D mesh. For a node
                          variable, the expected size should be the product of
                          (sizes[i] + nbdyl[i] + nbdyr[i] + 1) 
                          (i = 0, 1, 2) of the associated 3D mesh if
                          element_centered of the mesh is 0, and should
                          be the product of (sizes[i] + nbdyl[i] + nbdyr[i])
                          of the associated 3D mesh if element_centered
                          of the mesh is 1. */
      };  
typedef struct mio_Var_Comps mio_Var_Comps; 

/**
\ingroup io Functions for Writing Data
 mio_Mesh_Var defines a variable associated with meshes.
*/

struct mio_Mesh_Var {
       int  id;
       char *name;  
       int  mesh_ids[MAX_MESHES_IN_UMESH]; /**< It is possible that a node-variable
                                                is associated with a set of unstructured
                                                meshes if these meshes share one set of
                                                coordinates.
                                            */ 
       int  num_meshes;        /**< default is one  */
       mio_Object_Type mesh_type; /**< mio_smesh or mio_umesh, the type of the mesh the 
                               variable is defined on. This member is used only in the
                               output of the function, mio_query. */         

       mio_Mesh_Var_Type type; /**< zone-, face-, edge- or node-variable */
       int rank;               /**< 0 for a scalar, 
                                    1 for a vector;
                                    > 1 for a tensor. */
       mio_Data_Type datatype; /**< datatype comp[i].buffer           */

       int  compressed; /**< not implemented yet, intent to cover the
                             case with buffer[0] being all values on
                             the domain of this PE.           */

       int comp_sizes[MAX_COMP_RANK]; /**< For a scalar, this field is ignored.
                                For a vector, comp_sizes[0] is the number of
                                components of the vector.
                                For a tensor, comp_sizes[0], comp_sizes[1], ...,
                                comp_sizes[rank-1] are the numbers of
                                components for each index. For example,
                                for a tensor Tij (i = 0, 1, 2; 
                                j = 0, 1, 2), rank = 2, comp_sizes[0] = 3;
                                comp_sizes[1] = 3.                        */
       
       mio_Var_Comps comps[MAX_COMP_RANK*MAX_COMP_RANK*MAX_COMP_RANK];
      };
typedef struct mio_Mesh_Var mio_Mesh_Var;

/**
\ingroup io Functions for Writing Data
 mio_Mat_Group defines a set of materials with a mesh.
*/

struct mio_Mat_Group {
       int id;                        /**< The id of this group of materials */ 
       char *name;
       int meshid;
       int nmat;
       char *names[MAX_MATS];         /**< name of the material                      */
       int  matids[MAX_MATS];         /**< material id assigned by mio function calls  */

/***   The following line should be changed to nelems[MAX_MATS] even for structured meshes. 
       AND, mio.c and test codes have to be modified. Sepetember 9, 2008. 
 **/
       long long nelems[MAX_MATS][3]; /**< the number of elements where the material is
                                            possiblly present for each dimension. For unstructured
                                           meshes, only nelems[MAX_MATS][0] will be used.  */
       long long offsets[MAX_MATS];   /**< output, the offset of the elements on this pe   */
       long long gsizes[MAX_MATS];    /**< output, the total number of elements */
       void *indices[MAX_MATS][3]; 
                                     /**< For structured meshes,
                                          indices[0] is an array for indices of elements
                                          in the dimension 0, indices[1] is for dimension 1,
                                          and indices[2] is for the dimension 3.
                                          For unstructured meshes, only indices[0] is
                                          used which is the ids of elements.                 */
       mio_Data_Type datatype;       /**< data type of indices   */
      };
typedef struct mio_Mat_Group mio_Mat_Group;

/**
\ingroup io Functions for Writing Data
 mio_Matvar defines a variable defined on a material region.
*/

struct mio_Matvar_Group {
       char *name;             /**< name of the variable */
       int  id; 
       int  matgrpid;      
       mio_Data_Type datatype;
    };
typedef struct mio_Matvar_Group mio_Matvar_Group;

struct mio_Matvar {
       int  matvargrpid;       /**< the id of a mio_matvar_group */ 
       char *name;             /**< name of the variable */
       int  matid;              /**< input   */ 
       int  rank;               /**< 0 for a scalar, 
                                    1 for a vector;
                                    > 1 for a tensor. */
       mio_Data_Type datatype;     
       int comp_sizes[MAX_COMP_RANK]; /**< For a scalar, this field is ignored.
                                For a vector, comp_sizes[0] is the number of
                                components of the vector.
                                For a tensor, comp_sizes[0], comp_sizes[1], ...,
                                comp_sizes[rank-1] are the numbers of
                                components for each index. For example,
                                for a tensor Tij (i = 0, 1, 2; 
                                j = 0, 1, 2), rank = 2, comp_sizes[0] = 3;
                                comp_sizes[1] = 3.                        */

       mio_Var_Comps comps[MAX_COMP_RANK*MAX_COMP_RANK*MAX_COMP_RANK];
     };
typedef struct mio_Matvar mio_Matvar;

/**
\ingroup io Functions for Writing Data 

  There are two usages of mio_init. The first one is to initialize a new object
  before the object is created in a file. The object may be a group, or array,
  or mesh, or variable. The second usage is to get the information of an existing object
  through a given name.  
   
 - type: an input, which must be one of the following:
         mio_group, mio_array, mio_smesh, mio_smesh_cell_amr, mio_umesh, mio_mesh_var.     

 - group_id: an input, the id of an existing group, which
         may be the id of a file. A valid group_id is
         required only for obtaining the information of an existing 
         object. To initialize an object, set group_id to a 
         negative integer.

 - obj: an input and output. It must be an object of the
        following: mio_group, mio_array, mio_smesh, mio_smesh_cell_amr, mio_umesh, 
        mio_mesh_var. The object must be allocated before the call.
        If group_id >= 0, obj->name must be a valid name.

 - To initialize a new object, set group_id to -1, and this function will set
   the members of the object to their default values, and assign an id to obj->id.

 - If group_id is a non-negative, and obj->name != NULL,
   this function will try to find the id of the object through group_id
   and obj->name. If no object is found with
   obj->name under group_id, obj->id will be set to -1 in the output, and
   -1 will be returned. 
    
-  To find an existing array, set a valid group_id and obj->name
   before the call. If the array is found, id and datatype of the array will be given,
   and the members of obj->array_struct will be allocated which include
   sizes, offsets, nbdyl, nbdyr, psizes, and plist. For this case,  sizes, 
   offsets, nbdyl, and nbdyr have different meanings depending on the number
   of PEs calling this function. If the number of PEs calling this function
   is the same as the original number of PEs through which the array was written,
   the function sets their values on the original PE. If the number of PEs is different,
   this function sets the members to those on the first original PE on which
   the array has non-zero sizes. After the call, the allocated members of array_struct
   may be freed through mio_clean,
\verbatim
   mio_clean(mio_array_struct, 1, &(array->array_struct))
\endverbatim
or 
\verbatim
   if (sizes)   free(sizes);
   if (offsets) free(offsets);
   if (nbdyl)   free(nbdyl);
   if (nbdyr)   free(nbdyr);
   if (psizes)  free(psizes);
   if (plist)   free(plist);
\endverbatim

 - If type is mio_smesh, and an existing mesh is found, this function will assign
   the values for the following members of the mesh: 
   id, dims, npes, datatype, element_centered, gsizes, sizes, offsets, nbdyl and nbdyr
   (if there are ghost elements in the mesh), dcoord and coordmin (if the mesh is
    uniform). The members of the mesh in output, sizes, offsets, nbdyl and nbdyr, 
   dcoord and coordmin, have different meanings depending on the number of PEs
   calling this function. If the number of PEs calling this function is the same
   as the original number of PEs which wrote the mesh, the function sets the
   members to the values on the original PE. If the number of PEs is different,
   this function sets the members to those on the first original PE on which 
   there are mesh elements.

 - If type is mio_umesh, and an existing mesh is found, this function will
   assign the following members of the mesh:  id, dims, npes, idmin,
   order_for_nodelist, type, datatype, offsets, sizes, gsizes, ssizes, bsizes,
   msize, and coord. In gsizes[i] and offsets[i], i = 0, 1, 2, 3, only total size
   of elements, and offset of elements are meaningful, and they are gsizes[0] and
   offsets[0] if gsizes[0] != 0, they are gsizes[1] and offsets[1] if
   gsizes[0] = 0 and gsizes[1] != 0, and so on. The meaning of the members of
   the mesh, offsets, sizes, fsizes, ssizes, bsizes, and msizes, have different 
   meanings depending on the number of PEs calling this function.
   If the number of PEs calling this function is the same as the original number
   of PEs which wrote the mesh, the function sets the members to the values on
   the original PE. If the number of PEs is different, this function sets the
   members to those on the first original PE on which there are mesh elements.

 - If type is mio_smesh_cell_amr, and an existing mesh is found, this function
   will give the following members, datatype_coord, datatype_ids, npes, gsize, dims,
   idmin, element_centered, offset, size, fsize, size_masterelems, bsize, size_gid,
   and size_elems_connected_to_elm. The meaning of the members of
   the mesh, offset, size, fsize, ssizes, size_gid, and size_elems_connected_to_elm,
   depends on the number of PEs calling this function.
   If the number of PEs calling this function is the same as the original number
   of PEs which wrote the mesh, the function sets the members to the values on
   the original PE. If the number of PEs is different, this function sets the
   members to those on the first original PE on which there are mesh elements.

 - If type is mio_mesh_var, and an existing variable is found, this function
   will give the following members of the variable, type, rank, comp_sizes,
   datatype, num_meshes, mesh_ids, and mesh_type.  

 */  
int mio_init(int fileid, mio_Object_Type type, int group_id, void *obj);

/**
\ingroup io Functions for Writing Data

 THIS IS A COLLECTIVE CALL, i.e.,
 all the PEs in the communicator must call this function.

 mio_write writes one of the following objects:
 mio_Group, mio_Attr, mio_Array, mio_Structured_Mesh, mio_Unstructured_Mesh, mio_Coord,
 and mio_Mesh_Var. If the object is one of mio_Array, mio_Structured_Mesh,
 mio_Unstructured_Mesh, mio_Mesh_Var, 
 it must be initialized through mio_init before the call of this function.
 To write a mio_Group, or
 mio_Attr, or mio_Coord, the initialization is not needed. 

 - type: an input. It must be one of the following: mio_grp, mio_attr, mio_array,
   mio_smesh, mio_smesh_cell_amr, mio_umesh and mio_mesh_var.

 - obj_id: an input, the id of an object. If type is not mio_attr, id must be for
       a group or file, and a new object will be written under the group or file.
       If type is mio_attr,
       id may be the one for a group, array, structure mesh,
       unstructured mesh, variable, and the attribute will be attached to the object. 

 - Except for an attribute,
   the member of the object, name, may contain sub-groups, for example, 
   obj->name = "grp1/grp2/object1", obj->name = "/grp0/grp1/grp2/object1".
   A leading "/" in obj->name stands for an absolute path starting
   from the root group, i.e., the file. 
   If any group or sub-group contained in obj->name
   doesn't exist in the file, this function will create the
   group or sub-group. 

 - For the case with type = mio_grp and grp = (mio_Group *) obj, 
   if grp->name is given before the
   call, eg, grp->name = "/grp1/grp2/mygrp", or 
   grp->name = "grp1/grp2/mygrp", this function will set
   the id of the group, mygrp, in grp->id. See example_open.c.

 - If type = mio_attr, this function will attach the attribute, 
   obj, to another object specified by id. See example_attr.c.

 - To write an array, the object, obj, has to be initialized through mio_init befoe
   the call. After the initialization, the following members of mio_Array and
   members of mio_Array_Struct must be set before the call: name, datatype, buffer
   of array, dims and sizes of array_struct. If one PE gives offsets of
   array_struc, all the PEs must give offsets, and if one PE does not provide
   offsets, all PEs must set offsets to NULL (or each PE leaves offsets alone
   after mio_init). This rule is also true for gsizes,
   nbdyl, nbdyr, psizes, and plist. The member, npes, of array_struct desn't
   have to be given. The member of an array, buffer, is the starting address of the
   part of array in the current PE. The dimension 0 is the one along which the address
   changes most slowly, and dimension (dims-1) is the one the address changes fastest.
   see example_write_array.c.
   
 - If one PE contains ghost points, all PEs must set nbdyl or nbdyr, and values
   of nbdyl[i] or nbdyr[i] may be zero. Even if there are ghost points, sizes, offsets,
   and gsizes of array_struct do not include ghost points. The starting address,
   buffer of array_struct, includes ghost points. 

 - If type is mio_array, and offsets is not provided, the function will
   calculate offsets according to psizes and plist. 
   psizes is for the sizes of PEs in a PE-configuration.
   For example, for the PE-configuration shown in the table below,
   psizes[0] = 3, and psizes[1] = 2. If psizes is not
   provided, it will be assumed that psizes[0] = npes, and psizes[i] = 1
   for i != 0. The member, plist, of array_struct is a list of PE ranks in the
   order of dims dimension as shown the table below.
   The size of plist is psizes[0] * psizes[1] * ... * [dims-1].
   If plist is not provided, a default plist will be assumed which is
   the layout in the order of dimensions (dims-1), (dims-2),...,0.
\verbatim
           ----------------------------------------
    dim 1  |            |            |            |
           | plist(1)   | plist(3)   | plist(5)   |
      ^    |            |            |            |
      |    ----------------------------------------
      |    |            |            |            |
      |    | plist(0)   | plist(2)   | plist(4)   |
           |            |            |            |
           ----------------------------------------
                        ----> dim 0
\endverbatim

 - For the case with type = mio_array and array = (mio_Array *) obj,
   it is recommended that all the arrays with the same PE configuration use
   the same array_struct, which is shown in example ??. In the case, the first array
   is written as specified, and mio_write will set array_struct.id of
   the first array. For the other arrays with the same PE configuration as the first
   array, users have to set only array_struct.id to the one in the first array.
   See example_write_array2.c.

 - If type is mio_smesh,  on each PE, a mesh may be either uniform or nonuniform in any
   dimension. But, coord and/or (dcoord, coordmin) on all PEs should
   constitute a structured mesh.
   If one PE specifies offsets, or nbdyl, or nbdyr, or psizes, or
   plist of mio_Structured_Mesh, all the PEs have to specify the member.
   If some PE doesn't set
   the member, all the PEs leave the member alone after mio_init.
   If the member, offsets, is set, psizes and plist are ignored.
   If the member, offsets, is not given, psizes and plist will be used to
   calculate offsets. psizes is for the sizes of PEs in a PE-configuration.
   For example, for the PE-configuration shown in the table below,
   psizes[0] = 3, and psizes[1] = 2.
   plist is a list of PE ranks in the order of dims dimensions as shown in
   the table below.
   If psizes and/or plist are not given, default psizes and/or plist will be used. The
   default psizes is psizes[0] = the number of PEs and
   psizes[i] = 1; i = 1, 2. The default plist is
   plist[i] = i, i = 0, 1, ....
   See example_write_uniform_smesh.c and example_write_non_uniform_smesh.c.

\verbatim
           ----------------------------------------
    dim 1  |            |            |            |
           | plist(1)   | plist(3)   | plist(5)   |
      ^    |            |            |            |
      |    ----------------------------------------
      |    |            |            |            |
      |    | plist(0)   | plist(2)   | plist(4)   |
           |            |            |            |
           ----------------------------------------
                        ----> dim 0
\endverbatim

- For the case with type = mio_umesh, obj needs to be initialized through mio_init, but
   coord (= obj->coord) doesn't have to be initialized.
   There are two ways users may setup and write
   the coordinates of a mesh.
   a) After a call of mio_init, users set coord->datatype and coord.coord[i].
   In this case, coord.id, coord.num_meshes, and coord.mesh_ids[0] will be
   set in coord after the call of this function.
   b) After a call of mio_init, users set a valid coord.id
   which was output in a previous call of mio_write. In this case,
   the same coordinate id will be used
   for the current mesh, and num_meshes and mesh_ids[0] of the current mesh
   will be set in the function.
   If users write a mesh without coordinates specified,
   don't touch coord after mio_init.
   See example_write_coord.c.

 - If type is mio_umesh, and offset of mesh elements is not provided for some PE,
   plist of the mesh object will be used to calculate the element offset of each PE. 
   In this case the element order is determined through plist[0], plist[1], ..., plist[npes-1].
   If plist is not specified,
   a default plist will be used, which is plist[i] = i, i = 0, 1,..., npes-1.
   For a zone-mesh, the offset of mesh elements is offsets[0], for a face-mesh,
   the offset is offsets[1], for an edge-mesh, the offset is offsets[2], and for node-mesh,
   the offset is offsets[3]. The function will give the offset of mesh elements in the current PE, 
   and the total number of the elements.

 - If users want to write a same set of coordinates to a number of unstructured meshes which
   were previously written without coordinates specified, these meshes must have the same number of
   nodes. In this case, before the call, users must set datatype, coord[0], ..., coord[dims-1],
   num_meshes, and mesh_ids of obj.        
   See example_write_coord.c.

 - If type is mio_var or mio_matvar, it is suggested that obj be initialized through mio_init
   before the calling. The initialization is required mio_write writes only a part of components
   of a variable. 

 - To write a variable associated with a mesh, see example_write_smesh_variable.c and
   example_write_umesh_variable.c.

 */ 

int mio_write(int fileid, mio_Object_Type type, int obj_id, void *obj);

/** 
\ingroup io Functions for Writing Data

 mio_init_append prepares an operation which writes an object through more than
 one call of mio_write. Currently, the object can only be mio_umesh, mio_var, mio_array.
 Other objects have not been tested yet.

- type: an input. It must be one of the following: mio_array, mio_smesh, mio_umesh. 

- group_id: an input, the id of an existing group, which may be the id of a file.

- obj: an output. This function sets obj->id, which will be used as the id of the object.

- If type is mio_umesh, users don't have to set any member of obj.
  This function sets obj->id as the id of the object. 
  The procedure to write a mesh through
  mio_init_append should be something like the following
  See example_append_umesh.c and example_append_umeshvar.c.

\verbatim               
       mio_Unstructured_Mesh mesh;
       mio_init_append(mio_umesh, &mesh);
       // set members of mesh 
       ....
       mio_write(mio_umesh, &mesh);
       // set members of mesh
       .... 
       mio_write(mio_umesh, &mesh);
       ...
       mio_finalize_append(mio_umesh, mesh.id);
\endverbatim

The call of mio_init_append sets mesh.id, and the first call of mio_write sets
mesh.coord.id. The coordinate id, mesh.coord.id, should not be changed before
any subsequent call of mio_write. 
 
*/ 
int mio_init_append(int fileid, mio_Object_Type type, void *obj); 

/**
\ingroup io Functions for Writing Data
 
  mio_finalize_append finalizes an object initiated through mio_init_append. 
  All objects initiated through mio_init_append must be finalized through
  this function before the file is closed.
  See example_append_umesh.c and example_append_umeshvar.c.
*/  
int mio_finalize_append(int fileid, mio_Object_Type type, int id);

int mio_finalize(int fileid, mio_Object_Type type, int id);

/**
\defgroup query Functions for Querying Relationships

  This section covers the functions to query files, objects, and relationships
  between objects. 
*/

/**
\ingroup query Functions for Querying Relationships 

 mio_query gives the information of all the objects related to a particular object.
 Its functionality includes the following:
 1) Find the number of the objects that meet the given criteria. 
 2) If the particular object is a group, this function gives all the information of
 its immediate children (non-recursive), or all the generations under the group
 (recursive).
 3) If the particular object is a mesh, this function gives the
 information of all the variables associated with the mesh. 
 4) If the particular object
 is a variable, this function gives the information about the mesh on which the
 variable is defined. 
 5) If the particular object is a coordinate of an unstructured
 mesh, this function gives the information of all the meshes which use the 
 coordinate. 
 6) If the particular object is an array, or mesh, or variable,
 this function gives all attributes attached to the object.
 In the cases 1 to 5, the search may be limited to a specific name of
 object to be found. 

 - type: an input, the object type of the objects the function will list.
   It must be one of the following: mio_grp, mio_array, 
   mio_smesh, mio_umesh, mio_mesh_var, mio_coord, and mio_attr. 

 - obj_id: an input, the id of an object. 

 - path: an input, the path of an object relative to obj_id. It
   may be one of the following forms, for examples, 
   1) path = NULL; 
   2) path = *; 
   2) path = grp/sgrp/obj_name;
   3) path = grp/sgrp/[*].  
   Also, if there is a leading '/' in path, such as /grp/sgrp, the path is
   an absolute path, relative to the root group.  In this case, obj_id may be
   the id of any object in this file.

 - filter: an input, the name to be matched with.
       A non-null filter limits the scope of a search in which
       names of objects should match with the filter.

 - nobjs: an input or output, the number of objects which meet
   the specifications of the search. If nobjs = -1 in the input,
   through nobjs, this function will set nobjs
   to the number of objects which satisfy the conditions.
   In this case, obj is ignored, and this function only finds the
   number of relavant objects.
   If nobjs is not -1 in the input, nobjs is an input if objs is allocated before
   the call, and nobjs is an output if objs is set to NULL before the call.

 - objs: an output, an array of objects of the given type.  

 - If nobjs is not set to -1, and objs is set to NULL before the call,
   this function will calculate
   nobjs and allocate the memory for objs, and users are expected
   to cleanup through mio_clean after the use.
\verbatim
      mio_clean(type, *nobjs, objs);
\endverbatim

 - If the particular object defined
   by obj_id and path is a group, for example,
   path = grp/sgrp, and type is mio_grp, or mio_array,
   or mio_umesh, or mio_smesh, or mio_smesh_cell_amr, or mio_mesh_var, this function will give
   all the immediate children of sub-groups specified by type. If path is something like
   path = grp/sgrp/[*], this function will give all the children under sgrp specified
   by type (recursive), and the member, name,
   of each object in the output, objs, is relative to the particular object, sgrp. 
   See example_query.c.

 - If the particular object defined by obj_id
   and path is a variable and type is
   a mesh, mio_umesh, or mio_smesh, or mio_smesh_cell_amr, and type is mio_mesh_var,
   this function will give all the variables specified by type, and 
   defined on the particular object (mesh).
   In this case, the member, name,
   of any object in the output, objs, is an absolute path relative to the root.
   See example_query_association.c.  

 - If the particular object defined by obj_id and path
   is a group, and type is mio_smesh,
   the function sets the following members of mio_Structured_Mesh,
   id, name,  datatype, element_centered,
   dims, gsizes, sizes, offsets, nbdyl and nbdyr (if there is ghost elements).
   This function also gives
   smesh->dcoord[j] and smesh->coordmin[j] if the mesh is globally uniform in dimension j, otherwise
   smesh->dcoord[j] < 0 in output. The meanings of the members, sizes, offsets, nbdyl and nbdyr, 
   depend on the number of PEs calling this function.
   If the number of PEs calling this function is the same as the original number of PEs
   which wrote the mesh, the function sets the members to the values on the original PE. If the number
   of PEs is different, this function sets the members to those on the first original PE on which
   the mesh has elements. The members, psizes and plist, are set
   to NULL in this function.

 - If the particular object defined by obj_id and path 
   is a group, and type is mio_umesh, the function sets the following members of
   mio_Unstructured_Mesh, id, name, type, dims, datatype, npes, idmin,
   order_for_nodelist, msize, coord.datatype, gsizes, offsets, sizes, fsizes, 
   ssizes, and bsizes. The meanings of the members,   
   offsets, sizes, fsizes, ssizes, and bsizes, 
   depend on the number of PEs calling this function.
   If the number of PEs calling this function is the same as the original number of PEs
   which wrote the mesh, the function sets the members to the values on the original PE. If the number
   of PE is different, this function sets the members to those on the first original PE on which
   the mesh has elements. The member, plist, is set
   to NULL in this function.

 - If the particular object defined by obj_id and path is a variable, and type is
   mio_umesh, or mio_smesh, or mio_smesh_cell_amr, this function gives all the
   meshes on which the particular variable
   is defined. In this case, the member, name,
   of an object in the output, objs, is an absolute path relative to the root. 

 - If the particular object defined by obj_id and path is a group, and type is mio_coord,
   this function will give all the coordinates of unstructured meshes.
   In this case, filter is ignored.

 - If the particular object defined by obj_id is a coordinate, and type is mio_umesh,
   this function will give all the unstructured meshes specified by filter, which use the coordinate.

 - If the particular object defined by obj_id and path is a group, and type is mio_array, 
   this function will allocate the members of array->array_struct,  
   gsizes, sizes, offsets, nbdyl and nbdyr (if there are ghost
   elements). In the output, dims and multiple of array_struct are set in output. The meanings
   of the members of array_struct, sizes, offsets, nbdyl and nbdyr, depend 
   on the number of PEs calling this function. If the number of PEs calling
   this function is the same as the original number of PEs
   which wrote the array, the function sets the members to the values on the original
   PE. If the number of PEs is different, this function sets the members to those on
   the first original PE on which the array has non-zero sizes. The members, psizes
   plist, of array_struct, are set to NULL in this function. This function also sets
   array->datatype. If objs is set to NULL before the call, users should clean objs
   through mio_clean(type, nobjs, objs) after the use. If objs is allocated before the call,
   users should cleanup array->array_struct of each array through mio_clean.
   See example_query_array.c.
\verbatim
      mio_Array *array;
      int n = 1;
      mio_Array *arrays = (mio_Array *) objs;
      for (i = 0; i < *nobjs; i++) {  
          array = arrays + i;  
          mio_clean(mio_array_struct, n, &(array->array_struct));
      }
\endverbatim
 
 - If type is mio_attr, this function will try to find attributes attached to the object
   specified through obj_id and path. In this case,
   filter may be used to specify the name of a particular attribute.
   If objs is set to NULL before the call, this function will allocate objs for a list
   of attributes, and name and values of each attribute. If attr->datatype of an attribute 
   is mio_char, attr->count of the attribute is the number of characters in
   attr->values, and attr->count
   doesn't include the terminating character '\\0', although attr->values contains the terminating
   character. 
   See example_query_attr.c. 

 - If type is mio_obj_type_invalid, obj->id is given before the call,
   this function will give the name and general information of the object.

*/

int mio_query(int fileid, mio_Object_Type type, int obj_id, char *path, 
              char *filter, int *nobjs, void **objs);

/**
\ingroup query Functions for Querying Relationships
 
  mio_clean cleans the memory allocated in mio_query.
  Function mio_clean deallocates all the memory pointed through
  the pointers in obj. 
*/

int mio_clean(mio_Object_Type type, int nobjs, void *objs);

/** 
\ingroup query Functions for Querying Relationships

  mio_get_size gets the sizes of a part of an object. The object must be one of mio_Array, 
  or mio_Structured_Mesh, or mio_Unstructured_Mesh, or mio_Structured_Mesh_Cell_AMR.
  This function will first try to find the object
  through obj->id. If the function failed to find the object, it will try to find the object
  through group_id together with obj->name. The part of an object may be specified through a rank of
  the PE which previously wrote the object, or through offset and size of elements, or through
  mio_coord_domain. The possible ranks of PEs may be obtained through
  mio_query. This function also gives the sizes of connectivity arrays for
  the part of an unstructured mesh. 

  - type: an input, which must be one of the following, mio_array, mio_smesh,
          mio_umesh, mio_smesh_cell_amr.

  - group_id: an input, the id of a group. If obj->id is a valid id for a specified object,
    group_id is ignored. This input is used only when obj->id is not a valid id of an existing
    object.

 - obj: an input and output. It must be allocated before the call.

 - domain: an input. It specifies a part of the object. 

 - The object, obj, must be allocated before the call, but doesn't have to be initialized.
   There are two ways to specify an object: a) obj->id, or b) group_id together with obj->name.
   for the first case, obj->name is ignored, and for the second case, set obj->id to -1 before
   the call. This function will first try to find the object through the first approach,
   and if it failed, it will try the second approach, and a valid id will be assigned to obj->id
   if it is successful.

 - If domain->type is mio_pe_domain, domain->obj 
   is expected to be a pointer pointing to an object of mio_PE_Domain.
   The possible values for PE rank
   may be obtained through the function, mio_query. If the object was written through
   one call of mio_write in parallel or serial, rank may be the one of any PE. If the object
   was written through mio_init_append, rank has a similar meaning. 

 - If type = mio_array and array = (mio_Array *) obj,
   domain->type can only be mio_pe_domain. For this case, this function will
   allocate the members of array->array_struct, gsizes,
   sizes, offsets, nbdyl and nbdyr (if there are ghost elements),
   set array_struct->id to -1, and  
   set obj->datatype, array_struct->dims, array_struct->multiple. The members, psizes and plist,
   of array_struct are set to NULL in this function. The allocated members
   of array_struct may be cleaned through mio_clean(mio_array_struct, n, &(array->array_struct)).
   See example_get_size_array.c.
 
 - If type = mio_smesh, domain->type can be only mio_pe_domain.
   For this case, this function will set the following elements of mesh, id, dims,
   element_centered, datatype, gsizes, offsets, sizes, nbdyl and nbdyr (if there are ghost
   elements), dcoord, coordmin. A negative (or positive) dcoord[i] indicates that
   the part of mesh is non-uniform (or uniform with width dcoord[i]) in
   dimension i.
   See example_get_size_smesh.c.

 - If type = mio_umesh, domain->type may be mio_pe_domain, or mio_elem_domain, or mio_coord_domain.
   All the pointers in obj are reset to NULL, and previously allocated pointers
   before the call will be lost. After the call, mesh->datatype and mesh->coord.datatype are set
   to the original data types. 
   See example_get_size_umesh.c. 
   
 - For the case with type = mio_umesh and domain->type = mio_coord_domain, this function
   has not been tested yet. 
   If domain->criterion = mio_touch,
   bsizes[0:2] will be set in output for the numbers of zones, faces, and edges, which are intersected
   by the domain, and bsizes[3] will be set for the number of nodes which are out of the domain, but are
   associated with the included elements. 
   If domain->criterion = mio_cover, normal bsizes[1:3] will be set in output for
   the numbers of boundary faces, edges and nodes specified in the original data.  

 - If type = mio_grp, and domain->type is mio_PE_Domain, this function will find the extents of
   all the meshes under the group, obj. The group may be specified through obj->id, or
   obj->name. 
 */

int mio_get_size(int fileid, mio_Object_Type type, mio_Domain *domain, int group_id, void *obj);

/**
\defgroup read Functions for Reading Data

  This section covers the functions to read attributes, arrays, structured meshes,
  structured meshes with AMR, unstructured meshes, and variables defined on meshes.
*/
 
/**
\ingroup read Functions for Reading Data 

   mio_read reads a part or a whole object for one or more of following objects:
   mio_Array, mio_Structured_Mesh,
   mio_Structured_Mesh_Cell_AMR, mio_Unstructured_Mesh, mio_Mesh_Var. The part of the object
   may be specified through each object itself, obj, or the input, domain.
   Except for arrays, all the necessary buffers in the object must be allocated before
   the call, and sizes of the buffers may be obtained through mio_get_size. 
   The datatype, obj->datatype,
   doesn't have to be the same as that stored in the file.  

 - type: an input, which must be one of the following:
   mio_attr, mio_array, mio_smesh, mio_umesh, mio_smesh_cell_amr.

 - domain: an input, which specifies the part of objects. 

 - objs: an input and output, a list of objects,
         which must be one of following: mio_attr,
         mio_Array, mio_Structured_Mesh, mio_Structured_Mesh_Cell_AMR, mio_Unstructured_Mesh, 
         mio_Mesh_Var, and it must be allocated before the call.
 
 - nobjs: an input, the number of objects in the list specified by objs. 

 - It is recommended that nobjs is set to 1 except for the case
   with domain->type = mio_coord_domain for unstructured meshes.
 
 - The each of object, obj,  in the list, objs, must be allocated before the call,
   but doesn't have to be initialized. 
   This function will first try to use obj->id to locate the object. If the object is found,
   and obj->name
   is null before the call, the function will allocate obj->name and set obj->name. In this
   case, users are expected to free the allocation through free( ).
   The obj->id may be obtained through mio_query, or mio_init, or mio_get_size. 
   If obj->id failed to locate the object, for example,
   if obj->id is negative, group_id together 
   with obj->name will be used to locate the object, and a valid id will be assigned
   to obj->id.

 - If type = mio_array, and array is an object in objs,
   the members, sizes and offsets, of
   array->array_struct must be allocated, and the values of them must be valid.
   The member, nbdyl or nbdyr, of array->array_struct may be set to either NULL
   or valid. The input, domain, will be ignored,
   and therefore, domain may be set to NULL for arrays. The part of
   the array may be specified through the
   members, offsets, sizes, nbdyl and nbdyr, of array_struct. The buffer, array->buffer, 
   may or may not be allocated before the call. If array->buffer is allocated before the call, the size of buffer
   must be consistent with sizes, nbdyl and nbdyr of array_struct. If array->buffer is not allocated before
   the call, array->buffer must be set to NULL before the call, and this function will allocate the memory
   for array->buffer, and users are expected to free the memory through free(array->buffer). 

 - If type = mio_smesh, and the structured mesh is an object in the list,
   currently only mio_index_domain is allowed for domain->type.
   If domain is set to NULL before the call, the part of mesh queried is specified
   through the members, offsets, sizes, nbdyl, and nbdyr, of smesh. 

 - If type = mio_mesh_var, var is an object in the list,
   there are two ways to identify the variable, var->id and var->name
   together with var->mesh_ids[0]. 
   To read a variable, the memberis, datatype and comps, of var have to
   be specified before the call.   
   This function may read all the components of a vector/tensor if buffers of the
   components are allocated before the call. If the buffer of a component is explicitly set to NULL,
   this function will not read the component. This function
   will read any component with non-null (var->comps[i]).buffer. 
   This function will not use num_meshes, mesh_type, and comp_sizes of var. 

 - If type = mio_mesh_var and the variable is defined on a structured mesh, the part of mesh on which the variable
   is defined is determined through the input domain. Currently, domain->type can be only mio_index_domain
   for this case. 

 - If type = mio_unesh and domain->type = mio_coord_domain, this function will read all the mesh elements
   which satisfy the criterion, domain->criterion. If domain->criterion = mio_cover, members of mesh, 
   bdrynodes, bdryedges, bdryfaces, bdryzones, fakenodes, fakeedges, fakefaces, fakezones, slipnodes,
   slipedges, slipfaces, slipzones in output are the information in the original file. 
   If domain->criterion is mio_touch, some zones, faces, edges which are read from the file are
   intersected by the domain, domain->obj. In this case,
   the members of mesh, bdryzones, bdryfaces, and bdryedges, 
   are used to identify these zones, faces and edges respectively, the sizes of these array are stored
   in mesh->bsizes[0,1,2], mesh->bdrynodes and mesh->bsizes[3] are for the nodes which are
   outside the domain, but are involved in the mesh elements within the domain.

 - If type = mio_coord, domain->type should be mio_pe_domain, and objs should be of
   mio_Unstructured_Mesh.

 - If type = mio_partition_domain, domain->obj should be an obj of mio_Partition_Domain (not implemented yet).

 - If type = mio_attr, and the data type of attr is set to mio_datatype_invalid, the function will
   allocate the necessary memory for attr->values, and read the attribute in the original datatype.
   Attributes can be read only in their original data types.

*/

int mio_read(int fileid, mio_Object_Type type, mio_Domain *domain, int group_id, void *objs, int nobjs);

/**
 This function has not been implemented yet.
 mio_init_buffer creates a buffer for the subsequent writing of small datasets, or opens
 an existing buffer. Currently, this function supports only arrays.

 - type: an input. Currently, it can be only mio_array.

 - group_id: an input, an id of a group, or a file. A collection of nobjs small datasets will be written
   under the group.  

 - obj: an input and output, a pointer pointing to an object. 
   Currently, it can be only pointer of array.

 - nobjs: an input or output, the number of datasets this buffer. For creating a buffer,
   nobj is an input, but for opening an existing buffer, nobjs is an output.

 - total_size: an input or output. For creating a buffer, total_size is an input,
   and it is the number of values of all the datasets on the current PE.
   For opening an existing buffer (nobjs must be set to -1 before the call),
   total_size is an output, and it is the total
   number of values of all datasets in the buffer on all the PEs.

 - For writing arrays, the following members of obj and array_stryct (== obj->array_struct)
   must be given before the call, obj->name, obj->datatype, array_struct->dims.
   Also, array_struct->psizes and array_struct->plist must be given before the call. All
   the datasets in the buffer are assumed to have the same dims, psizes (if specified),
   plist (if specified), and datatype. The member of obj, name,
   will be used as the path and name of the buffer relative to group_id, and obj->id will be an
   output, which will be used for the id of the buffer in the subsequent calls of mio_write.

- To open an existing buffer, users must give group_id, obj->name, and set nobjs 
   to -1. The correct nobjs will be returned, obj->id will be assigned for the id of buffer
   in the output of obj,  
   and obj->id will be used for the id of the buffer in the subsequent calls of mio_read.
   If the buffer can not be found, -1 will be returned.

*/ 

/****
int mio_init_buffer(mio_Object_Type type, int group_id, 
                    void *obj, int *nobjs, long long *total_size);

int mio_close_buffer(mio_Object_Type type, void *obj);

***/
/***
  The function, mio_get_domain, has not been tested yet.
  This function has been implemented only for the case with type being mio_umesh and
  domain->type being mio_coord_domain.
 ***/

int mio_get_domain(int fileid, mio_Object_Type type, int group_id, void *obj, mio_Domain *domain);

/**
\defgroup mpi Functions Related to MPI

 This section covers the functions related to MPI.
*/

/**
\ingroup mpi Functions Related to MPI

mio_set_info sets the keyword-value pair for MPI_INFO members.

 - keyword: an input.
 
 - value: an input.
*/
int mio_set_info(char *keyword, char *value);

/**
\ingroup mpi Functions Related to MPI

 mio_delete_info deletes the mpi info associated with keyword.
 
 - keyword: an input.
*/
int mio_delete_info(char *keyword);

/**
\ingroup mpi Functions Related to MPI

 mio_free_info frees MPI_info.

*/
int mio_free_info(void); 

/**
\ingroup mpi Functions Related to MPI

 mio_comm_create creates a MPI commucator.

 - n: an input, the number of PEs in the communicator to be created, or 0.

 - ranks: an input, the list of PE ranks.

 - If n is 0, or the total number of processors in MPI_COMM_WORLD, ranks
   will be ignored. In this case, MPI_COMM_WORLD, will be duplicated through 
   MPI_Comm_dup.
*/

int mio_comm_create(int n, int *ranks);

/**
\ingroup mpi Functions Related to MPI

 mio_comm_free free the communicator created before.
*/
int mio_comm_free();

/**
\page appdx Example Codes

 Example codes are in the directory ../examples.

*/

#ifdef __cplusplus
}
#endif

#endif
