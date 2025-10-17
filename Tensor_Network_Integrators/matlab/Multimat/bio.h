#ifndef bio_H
#define bio_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MPI
#include <mpi.h>
#endif

/**
\mainpage BIO Interface Users' Guide 
\author William W. Dai, HPC-4, Los Alamos National Laboratory
\version V-2006-001
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

\section acknowledgements Acknowledgments
 This research and development are funded by the Department of Energy's ASCI program.
 Los Alamos National Laboratory is operated by the University of California
 for the National Nuclear Security Administration of the United States
 Department of Energy under contract W-7405-ENG-36.

\section abstract Abstract
 This BIO users' Guide provides the information needed to write and
 read data on parallel or serial computer platforms, and to query files which
 are created through BIO library interface functions.
 This Guide also contains the typical usages
 of BIO interface functions and example odes for the usages. 
 
\section overview Overview

 BIO library is an IO library built directly on MPI and MPI-IO. 
 The goal of the BIO library is to provide simulation codes
 an efficient IO capability on parallel or serial computer platforms.
 The files created through BIO may contains hierarchical structures and
 descriptions of each array or group. Since the files are self-describing,
 users may read any part of a file without any knowledge of the file.
 The goal of the BIO Interface is to introduce the main functionality
 of BIO through a number of C/C++ callable functions.

 A file written through BIO library may have a hierarchical structure,
 which is similar to a Unix file system. In a Unix file system, directories can
 contain other sub-directories and files. The top directory is referred to as the
 "root" directory and has the "/" designation. Each file has
 certain attributes, e.g., read/write permissions, etc.
 In the file written through BIO library, the object corresponding to a
 directory is a group, the object corresponding to a Unix file is an array.
 A group is a container for other groups and arrays. 

 Attributes are associated with a group or array. Unlike groups
 and arrays, attributes can only written once the group or array
 they will be associated with has been created. Likewise, Attributes
 are read by referring to the location of the object to which they are
 associated. Attributes allow additional information
 to be associated with a group or array. With groups, arrays, and attributes
 users may create structured and self-describing files for their own data structures.

 For all the C/C++ interface functions, an id of an object is unique,
 and it is non-negative. A function will return 0
 if no errors were encountered; otherwise, -1 will be returned.
 All constants shown in the Users' Guide are included in the file

 - bio.h

*/

/**
\defgroup file Functions for Opening and Closing Files
 
 This section covers the functions to open and close files. 
 
*/

/**
\ingroup file Functions for Opening and Closing Files

 bio_IO_Mode contains the values for possible IO mechanisms.
 
*/

enum bio_IO_Mode {
     /** MPI_File_write_all */
     bio_collective   = 1,
     /** - MPI_File_Write_at */
     bio_independent = 2, 
     /** use fopen without MPI_IO */
     bio_use_fopen = 3,
     /** use open, without MPI_IO */
     bio_use_open = 4
}; 
typedef enum bio_IO_Mode bio_IO_Mode;

/**
\ingroup file Functions for Opening and Closing Files

 bio_File_Mode contains the values of open modes for files.
 
*/

enum bio_File_Mode {
     /** - create */
     bio_file_create   = 1,
     /** - read only */
     bio_file_read_only = 2, 
     /** - read and write, not implemented yet */
     bio_file_read_write = 3
}; 
typedef enum bio_File_Mode bio_File_Mode;

/**
\ingroup file Functions for Opening and Closing Files

 Given a relative path or full path of an existing file,
 bio_correct_file return 1 if the file is wriiten through bio.
 Otherwise, the function returns 0.

 - file_name: an input, the relative path or full path of an existing file.
*/

int bio_correct_file(char *filename);

/**
\ingroup file Functions for Opening and Closing Files

 bio_resilient is only for writing data. By default, bio library doesn't
 do any thing for resilience. But if this function is called with a non-zero
 argument, the following files to be created and data to be written will be
 resilient in the sense that all the data previously written could be read
 should io crashes before a file is closed.

 - ifresilient: an input, 0 or 1, and others will be ignored.
*/

int bio_resilient(int ifresilient);

/**
\ingroup file Functions for Opening and Closing Files

 bio_backup_file is used to setup the filename of backup file only when
 bio library is used for resilience. If this function is not called, the
 backup file will be the filename of the file to be created with the
 extension "backup". 

 - backup_file: the filename of the backup file.
*/

int bio_backup_file(char *backup_file);

/**
 * \ingroup file Functions for Opening and Closing Files
 *
 *  By default, bio reads all the buffered data through buffering, i.e, 
 *  bio reads the whole buffer and keeps the buffer in memory for any
 *  array. bio_buffer_mode resets the mode to read any arrays in buffers.
 *  This function has no effects in writing. 
 *
 *     - flg_for_buffered_read : an input
 *                               0          for non-buffered read
 *                               otherwise, for buffered read.
 *     */

int bio_buffer_mode(int flg_for_buffered_read);

/**
\ingroup file Functions for Opening and Closing Files

 Given a relative path or full path for a file,
 mio_file_open will create a new file or open an existing file.
 This function
 also automatically creates or opens the root group. The output, file_id,
 is the id of the file, and it is also the id of the root group.

 - file_name: an input, the relative path or full path for a file,
              such as restart/file1.
 
 - mode: an input, one of the three values of bio_File_Mode.
 
 - file_id: an output, the id of the file just created or opened.

*/

int bio_file_open(char *name, bio_File_Mode mode, bio_IO_Mode io, int* file_id);

/* This function is for the use in Sage. If buf_size < 0, no buffering. Otherwise,
   buf_size will be used for the number of bytes of the buffer. */

int bio_file_open_for_sage(char *name, bio_File_Mode mode, bio_IO_Mode io, 
                           long long bufsize, int *file_id);

/**
\ingroup file Functions for Opening and Closing Files
 
 bio_group_open will open a sub-group under a group or the root group.
 
 - parent_id: an input, the id of the parent group.  
   To open a group under a root, file_id should be used as parent_id.
 
 - name: an input, the name of a group to be opened/created.
 
 - group_id: an output, the id of the group just created or opened.

 - If there already exists a group with the name, name, under parent_id in the file,
   no new groups will be created, but a group_id will be assigned. Otherwise
   the function will create a group with the name.

*/
int bio_group_open(int fileid, int parent_id, char* name, int *group_id);
/**
\ingroup file Functions for Opening and Closing Files

 bio_file_close will close the file associated with a given id, and reurns
 the size of the file in bytes if the file was newly created or rewritten
 (or 0 if the file was open for read).

 - file_id: an input, the id of a file which was created or opened before.

*/
long long  bio_file_close(int file_id);

/**
\ingroup file Functions for Opening and Closing Files

 bio_file_flush will flush the file associated with a given id. This function is
 similar to bio_file_close, except for that a user may continue to write data to
 the file after calling bio_file_flush. This function returns the size of the file
 just written in bytes.

*/

long long bio_file_flush(int file_id);

/**
\defgroup array Functions for Writing and Reading

 This section covers the functions to write and read arrays.
*/

/**
\ingroup array Functions for Writing and Reading

 bio_Data_Type contains the allowed data types.

*/

enum bio_Data_Type {
     /** - char */
     bio_char             = 1,
     /** double */
     bio_double           = 2, 
     /** float  */
     bio_float            = 3,     
     /** int    */
     bio_int              = 4,     
     /** long   */
     bio_long             = 5,     
     /** long double */
     bio_long_double      = 6,     
     /** long long   */
     bio_long_long        = 7,     
     /* invalid      */ 
     bio_datatype_invalid = 16
}; 
typedef enum bio_Data_Type bio_Data_Type;

/**
\ingroup array Functions for Writing and Reading

 bio_Object_Type defines pre-defined descriptions of groups and arrays.

*/

enum bio_Object_Type {
     /** - group  */
     bio_group  = 1,
     /** - array */
     bio_array  = 2,
     /** - structured mesh */
     bio_smesh  = 3,   
     /** - unstructured mesh */
     bio_umesh  = 4, 
     /** - structured mesh cell AMR       */
     bio_smesh_cell_amr = 5,
     /** - unstructured_mesh_cell_amr */
     bio_umesh_cell_amr = 6,  
     /** - variable on structured mesh */
     bio_smeshvar = 7,
     /** - variable on unstructured mesh */
     bio_umeshvar = 8,
     /** - variable on structured mesh cell amr */
     bio_smesh_cell_amr_var = 9,
     /** - variable on unstructured mesh cell amr */
     bio_umesh_cell_amr_var = 10,

     bio_mat     = 11,
     bio_mat_var = 12, 

     bio_coord   = 13, 

     /** - collected data */
     bio_buffer = 64, 
     /** - not specified, don't care */
     bio_object_type_undefined = 128
};  
typedef enum bio_Object_Type bio_Object_Type;

/**
\ingroup array Functions for Writing and Reading

   bio_write is to write an array into a file, or copy an array into a buffer.
   This function is a collective call unless it is used together with bio_buffer_init.

   - grp_id: an input, the id of a group or buffer under which the array will be written.

   - type:  an input, any value of bio_Object_Type.
  
   - datatype: an input, the data type of the array.

   - offset: an input and output, the offset of the part of the array on this PE.

   - gsize: an input and output, the total size of the array. 

   - buffer: an input, the starting address of the part of the array.

   - array_id: an output, the id of the array.

   - For writing an array into a file, grp_id, name, type, data type and size must be given,
   and grp_id, name, type, datatype, and gsize, must be the same for all the PEs. 
   If one PE gives valid offset and gsize, all other PEs must give valid offset and gsize too.
   For this case, the function will use these values for offset and gsize. Users also may set
   offset and gsize to -1 on all the PEs so that the function will calculate offset and gsize
   according to PE ranks, and set offset and gsize in output. 

   - If the file has already an array with the name under the group,
   gsize must be given before the call, and
   gsize, type, and datatype, must be the same as those of the existing array. Otherwise.
   For this case, the function will rewrite the array.

   - If offset of own PE is known before the call, each PE should set its own valid offset. If
   any one of PEs doesn't know its own offset, all PEs should set offset to -1 before call. If
   all the PEs know gsize before the call, set the correct gsise. If any one of PEs doesn't
   know gsize, all the PEs set gsize to -1 before the call.
 
   - If grp_id is the id of a buffer, the array will be written to the buffer. For this case,
   gsize are ignored, and array_id is set to the id of the buffer.

   - To copy a array to a buffer, users must first initialize the buffer through bio_buffer_init.
   For this case, grp_id should be the id a buffer. If the argument collective in bio_buffer
   is set to 1, the subsequent calls of bio_write for coping arrays to the buffer are collective
   calls. For this case, name must be the same for all PEs. 

   - If the argument collective in bio_buffer_init is set to 0, the
   subsequent calls of bio_write for coping arrays to the buffer is non-collective calls. For
   this case, name in bio_write must be different on different PEs, ie, data on each PE is
   considered as a complete array, and offset and gsize are ignored, 
   and offset of the each array on each PE is considered zero.    
*/ 
int bio_write(int fileid, int grp_id, char *name,
             bio_Object_Type type, bio_Data_Type datatype,
             long long *offset, long long size, long long *gsize,
             void *buffer, int *array_id);

/**
\ingroup array Functions for Writing and Reading

   bio_buffer_init initializes buffers for subsequent calls of bio_write. bio_buffer_init,
   together with bio_write and bio_buffer_finalize, is designed for writting many small arrays
   for the best possible IO performance. 

   - grp_id: an input, the id of a group under which the buffer will be written.
  
   - name: an input, the name of the buffer to be written.

   - collective: an input. If collective is 1, the subsequent calls of bio_write are assumed
     to be collective calls. If collective is 0, the calls of bio_write are non-collective
     calls. 

   - datatype: an input. The data type of all the arrays in the subsequent collective calls.
     If is assumed that data types are the same for all the arrays to be written in one buffer.

   - buffer_size: an input, the size of the part of the buffer on this PE. buffer_size
     should be at least as large as the size of the longest array in this PE in all
     the subsequent calls of bio_write. The subsequent calls of bio_write will fill the buffer
     with buffer_size before they are actually written into a file. After the buffer is full,
     the arrays in the buffer are written to a file, and the buffer becomes empty to hold
     more arrays.

   - buffer_id: an output, the id of the buffer, which is needed in subsequent calls of
     bio_write and bio_buffer_finalize. 

   - If collective is 1, each subsequent call of bio_write(buffer_id, array_name, ...) is a
   collective call, and array_name must be the same for all PEs. If collective is 0, each
   subsequent call of bio_write is an independent call, and array_name on different PEs must
   be different.   

*/

int bio_buffer_init(int fileid, int grp_id, char *name, int collective,
                    long long buffer_size, int *buffer_id);

/**
\ingroup array Functions for Writing and Reading

   bio_buffer_finalize closeis a buffer which is opened through bio_buffer_init. This
   function is a collective call. Any buffer opened should be finalized before a file
   is closed.
*/

int bio_buffer_finalize(int fileid, int buffer_id);

/**
\ingroup array Functions for Writing and Reading

 - mio_attr_write writes an attribute to an object which may be a group or array.
 To write an attribute,
 all PEs are supposed to have the same name and values of the attribute.
 This function is a collective call.

 - id: an input, the id of a group or array.

 - name: an input, the name of an attribute.

 - datatype: an input, the data type of the attribute values.

 - size: an input, the number of elements in buffer for attribute values.

 - buffer: an input, the starting address of the values of the attribute.

*/
int bio_attr_write(int fileid,
         int id, char *name, bio_Data_Type datatype, int size, void *buffer);

/**
\ingroup array Functions for Writing and Reading

   bio_read reads a part, or whole, of an array under a group.

 - group_id: an input, the id of a group.

 - name: an input, the name of the array to be read.
   
 - datatype: an input, the requested data type, which may be different from 
             the one stored in the file.
 
 - offset: an input, the offset of the part of the array requested.

 - size: an input, the size of the part of the array requested.
   
 - buffer: an input or output, the starting address of the part of the array.
   
 - array_id, an output, the id of the array. 
   
 - If buffer is set to NULL before the call, this function will allocate the
   necessary memory for buffer.
\verbatim
       int *buffer = NULL;
       bio_read(..., (void **)&buffer);
       if (buffer) free(buffer);
\endverbatim

 - If buffer is not explicitly set to NULL, buffer should be allocated before
   the call, for example, 
\verbatim
       int *buffer = (int *) malloc(100 * sizeof(int));.
       bio_read(..., (void **)&buffer).
\endverbatim

 - buffer may be statically allocated before the call, for example,
\verbatim
       int p[100]; 
       int *buffer = p;
       bio_read(..., (void **)&buffer).
\endverbatim

*/ 

int bio_read(int fileid, int group_id, char *name, bio_Data_Type datatype, 
             long long offset, long long size, void **buffer, int *array_id);

/**
\ingroup array Functions for Writing and Reading

 mio_attr_read reads a given attribute of a given object.
 
 - id: an input, the id of an object.
   
 - name: an input, the name of an attribute.
   
 - datatype: an input or output, the data type of attribute values.
 
 - size: an input or output, the number of elements in attribute values.

 - buffer: output, the starting address of attribute values.
 
 - If buffer is set NULL before the call, this function will allocate the
   necessary memory for buffer, and the arguments, datatype and size, are
   output, for example, 
\verbatim
       int size;
       bio_Data_Type datatype;
       void *buffer = NULL;
       bio_attr_read(id, name, &datatype, &size, (void **)&buffer).
\endverbatim

 - If buffer is not explicitly set to NULL, buffer should be allocated before
   the call, and data type and size are input. buffer may be dynamically
   allocated before the call, for example, 
\verbatim
       int *buffer = (int *) malloc(100 * sizeof(int));.
       bio_attr_read(..., (void **)&buffer).
\endverbatim

 - buffer may be statically allocated before the call, for example,
\verbatim
       int p[100]; 
       int *buffer = p;
       bio_attr_read(..., (void **)&buffer).
\endverbatim
*/

int bio_attr_read(int fileid,
            int id, char *name, bio_Data_Type *datatype, int *size, void **buffer);

/**
\defgroup query Functions for Query

  This section covers the functions to query a file.

*/

/**
\ingroup query Functions for Query

 bio_List_Struct defines possible outputs from function
 bio_list.

*/
struct bio_List_Struct {
       char *name;               /**< the name of an object */
       int  id;                  /**< id of an object, invalid for an attribute */
       bio_Object_Type type;      /**< type of an object, invalid for an attribute */
       bio_Data_Type   datatype;  /**< data type if an object is not a group */
       long long       size;      /**< total size of the array if the object is not a group */
       void *value;              /**< values of an array or attribute  */

       int narrays;              /**< the number of arrays in a buffer */
       char **names;             /**< names of the arrays in a buffer */
       bio_Data_Type *datatypes; /**< the datatypes of these arrays */
       long long *sizes;         /**< the sizes of the arrays in a buffer */
                                 /**< names[0], names, and sizes need to be freed
                                      if they are not null */
     };
typedef struct bio_List_Struct bio_List_Struct;

/**
\ingroup query Functions for Query

  bio_count gives the number of groups and arrays under a given group.

 - group_id: an input, the id of a group.

 - filter: an input, a name to be matched with.

 - n: an output, the number of objects found.

 - If filter is set to NULL, this function gives the number of arrays
   and groups in a group. If filter is not NULL,
   the number of arrays and groups includes
   only those with their names matched with filter.

 */

int bio_count(int fileid, int grp_id, char *filter, int *n);

/**
\ingroup query Functions for Query

 bio_list gives the information of the groups, arrays, and buffers under a given group.

 - group_id: an input, the id of a group.
 
 - filter: an input, a name to be matched with. This input will not be used to match
   any array contained in a buffer. Instead, filter may be used to match the name of
   buffer.

 - to_read_value: an input. If it is not zero, the values of any array matching
   with the criterion will be read, otherwise, the values will not be read.
   This input doesn't have effect on buffers.

 - n: an input and output, the number of objects in the list, list.

 - list: an output, the list of objects found. 

 - If filter is set to NULL, this function gives the information of all the arrays
   and groups in a group. If filter is not NULL,
   the information of only those groups and arrays which matches with filter
   are given. 

 - If list is allocated before the call, list[i].name will be allocated in this function. 
   Users are expected to free the memory through free().  
   Any previous allocated memory for list[i].name will be lost. 

 - If to_read_value is set to 1, list[i].value will be allocated and be read in this function.
   Users are expected to free the memory through free().

 - Users may set *list explicitly to NULL. In this case, this function will
   allocate necessary memory for *list if there are objects found, and n is
   an output in this case. Users are expected to free allocated memory, for example,
\verbatim
       bio_List_Struct *list, *a;
       int i, n;
       list = NULL;
       bio_list(..., &n, &list);
       for (i = 0; i < n; i++) {
           a = list + i;
           if (a->name) free(a->name);
           if (a->value) free(a->value);
       }
       if (list) free(list);
\endverbatim

- The outputs, n and list, contain the information of buffers, but not any array written
   through bio_buffer_init. The arrays written through bio_buffer_init are shown through
   the members, names and sizes, of bio_List_Struct. For buffers, the member, value, of
   bio_List_Struct is ignored.

*/ 

int bio_list(int fileid,
             int group_id, char *filter, int to_read_value, int *n, bio_List_Struct **list);


/**
\ingroup query Functions for Query

 bio_free_list deallocates the spaces allocated in list through the call bio_list.
 
 n:     an input, the number of objects in the second argument, list.
 list:  an input. 
*/ 

int bio_free_list(int n, bio_List_Struct *list);

/**
\ingroup query Functions for Query

 bio_attr_list lists all the attributes associated with a given object. It also gives
 the name, data type, size, and values of each attribute.

 - id: an input, the id of a group or array.

 - filter: an input, a name to be matched with.

 - n: an input and output. 

 - list: an output.        

 - If filter is NULL, this function gives the information of all the attributes
   associated with the object. If filter is not null,
   the information of only those attributes which matches with filter
   are given. 
 
 - Users may set *list explicitly to NULL. In this case, this function will
   allocate necessary memory for *list if there are attributes found, and n is 
   an output. Users are expected to free the memory, for example,
\verbatim
       bio_List_Struct *list, *a;
       int i, n;
       list = NULL;
       bio_attr_list(..., &n, &list);
       for (i = 0; i < n; i++) { 
           a = list + i;
           if (a->name) free(a->name);
           if (a->value) free(a->value);
       }
       if (list) free(list);
\endverbatim
*/

int bio_attr_list(int fileid, int id, char *filter, int *n, bio_List_Struct **list); 

/**
\defgroup mpi Functions Related to MPI
   
 This section covers the functions which are directly related to MPI.
 
*/ 
   
/**
\ingroup mpi Functions Related to MPI 
   
 bio_comm_create creates a user-defined MPI communication world, which
 may be different from the default one, MPI_COMM_WORLD.
   
 - n: an input, the number of PEs in a user-defined MPI communication world.
   
 - ranks: an input, a list of ranks with the length of n.
   
 - If the user-defined communication world is the MPI default,
   i.e, MPI_COMM_WORLD, there is no need to call this function. If
   the user-defined communication world is a sub-set of MPI_COMM_WORLD,
   this function should be called before bio_open_file. 
   
 - The files which are all open must have the same
   user-defined communication world.
   
 - If mio_comm_create is called. mio_comm_free should be called
   after a file is closed.
*/ 
int bio_comm_create(int n, int *ranks);

/**
\ingroup mpi Functions Related to MPI
 
  bio_comm_free frees the communication world created through bio_comm_create.
*/

int bio_comm_free();

/**
\ingroup mpi Functions Related to MPI
 
  bio_get_comm get the communication world created through bio_comm_create.

  - comm: output, it must allocated before the call.
*/

#ifdef MPI
int bio_get_comm(MPI_Comm *comm);
#endif

/**
\ingroup mpi Functions Related to MPI

 bio_info_set sets keyword-value pairs for MPI_INFO members.

 - keyword: an input, keyword.

 - value: an input, value.
*/

int bio_info_set(char *keyword, char *value);

/**
\ingroup mpi Functions Related to MPI

 bio_info_delete delete a keyword from MPI_INFO members.

 - kw: an input, keyword.
*/
int bio_info_delete(char *keyword);

/**
\ingroup mpi Functions Related to MPI

  bio_info_free calls MPI_Info_free.

*/
int bio_info_free();

/**
\defgroup mis Functions for Miscellaneous Use

 This section covers the functions for some miscellaneous uses.

*/

/**
\ingroup mis Functions for Miscellaneous Use
 
  Given the id of an object, bio_get_name gives the information of the
  object.

  - id: an input, the id of an object.
  
  - obj: an output. 
*/

int bio_get_name(int fileid, int id, bio_List_Struct *obj);

/**
\ingroup mis Functions for Miscellaneous Use
 
  bio_get_next_id returns the next available id.
*/

int bio_get_next_id(void);

/**
\ingroup mis Functions for Miscellaneous Use


 bio_group_open_i is the same as bio_group_open
 except for that group_id in bio_group_open_i is given
 before the call. 

 - file_id: an input.
 - name:    an input.
 - grp_id:  an input. 
 - objtype: an input,  for buffering of mio  
 - datatype: an input, for buffering of mio
*/

int bio_group_open_i(int fileid, int parentid, char* name, int *grp_id,
                     bio_Object_Type objtype, bio_Data_Type datatype);

/**
\ingroup mis Functions for Miscellaneous Use

 bio_write_i is the same as bio_write
 except for that array_id in bio_write_i is given
 before the call.

 - grp_id: an input.
 - name:   an input.
 - type:   an input.
 - datatype: an input.
 - offset: an input or output.
 - size:   an input.
 - gsize:  an input or output.
 - buffer: an input.
 - array_id: an input.
*/

int bio_write_i(int fileid, int grp_id, char *name,
                bio_Object_Type type, bio_Data_Type datatype,
             long long *offset, long long size, long long *gsize,
             void *buffer, int *array_id);

int bio_debug(int fileid, int id);

#ifdef __cplusplus
}
#endif

#endif
