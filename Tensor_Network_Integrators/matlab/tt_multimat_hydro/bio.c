/********************************************************************
  Copyright (c) 2004, The Regents of the University of California. 
  All rights reserved.                                             
                                                           
             **************************************************
             **************************************************
             *****                                        *****
             *****          THIS SOFTWARE WAS BIO         *****
             *****       THIS SOFTWARE WAS WRITTEN BY     *****
             *****                                        *****
             *****            William W. Dai              *****
             *****                                        *****
             *****             November 2004              *****
             *****          Revised, August  2006         *****
             *****          Revised, October 2007         *****
             *****          Revised, October 2012         *****
             *****                                        *****
             *****         COPYRIGHT,  2004-2012          *****
             *****          ALL RIGHTS RESERVED           *****
             *****                                        *****
             **************************************************
             **************************************************

********************************************************************/

#include <stdio.h>   
#include <string.h>  
#include <stdlib.h>  
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h> 

#ifdef USE_MPI_IO
#define MPI
#endif

#ifdef MPI
#define __MPI__
#define USE_MPI_IO
#endif

// #ifdef __MPI__
// #define MPI
// #define USE_MPI_IO
// #endif

#ifdef MPI
#define USE_MPI_IO 
#endif

/****
#if defined(aix)
#define GLOBAL_FILESYSTEM 1
#elif defined(osf1)
#define GLOBAL_FILESYSTEM 1
#elif defined(linux)
#define GLOBAL_FILESYSTEM 1
#else
#define GLOBAL_FILESYSTEM 0 
#endif 
****/

#define PANFS

#ifdef PANFS
#define O_CONCURRENT_WRITE 020000000000
#else
#define O_CONCURRENT_WRITE 0
#endif

#ifdef MPI
#include <mpi.h> 
// #include <mpio.h>
#endif

#include "bio.h"

#define MAXLNAME 2048 

/**    #include "io_int.h"    **/

struct io_Child;
struct io_Tree;

struct io_Tree { int  id;
                 int  idx;
                 char *p;
                 struct io_Tree  *parent;
                 struct io_Child *child;
                };
typedef struct io_Tree io_Tree;

struct io_Child { io_Tree  *tree;
                  struct io_Child *next;
                 };
typedef struct io_Child io_Child;

struct io_Attr {
/*
           size,           szllong 
           id,             szllong 
           num,            szllong 
           size of name,   szllong 
           name,           slen*szchar
           datatype,       szllong 
           size of value,  szllong 
           values,         n*sizeof(datatype)
*/ 
       long long max_size;   /**< max size of the value of an object in terms of char */
       void *value;
       long long *size;   /**< current location in value  */
       long long *id;     /**< the id of the attrs attached to  */
       long long *num;    /**< the number of attrs  */
};
typedef struct io_Attr io_Attr;

struct io_Obj {
/* 
          size_meta,          szllong 
          fid,                szllong 
          num,                szllong 
          size_data,          szllong
          size of name,       szllong 
          name,               slen*szchar
          id,                 szllong 
          parent_id,          szllong  
          obj_type,           szllong 
          datatype,           szllong  
          offset,             szllong
          size,               szllong
*/ 
       long long  max_num;     /**< the max number of objects    */
       long long  max_size;    /**< the maxinum size of value of an object */
       void *value;
       long long *size_meta; /**< the current location in value  */
       long long *fileid;
       long long *num;       /**< the number of objects        */
       long long *size_data;  /**< the current location of data   */

       io_Attr **attr;       /**< a pointer array   */
      };
typedef struct io_Obj io_Obj;

struct io_File;
struct io_Buffer;

struct io_File {

#ifdef MPI 
       MPI_File  mfile;      /**< mpi file handler     */
       double    t_open;     /**< a specific time before opening a file   */  
       double    t_close;    /**< a specific time after closing a file    */  
#else 
       time_t t_open;     
       time_t t_close;  
#endif
       double mbytes;   /**< the size of datafile       */

       FILE         *fp;
       int           fd;        /**< file descriptor      */ 
       bio_File_Mode mode;      /**< file modes           */
       bio_IO_Mode   io;        /**< io modes             */
       int           id;        /**< the ids of files     */
       int           id_min;    /**< smallest id used in this file */
       io_Obj       *object;    /**< a pointer array      */
       char          *name;     /**< names of files       */

       struct io_Tree *root;

       struct io_Buffer *buffer;  /**< first buffer of the file if there are buffers. */

       struct io_File *next;   /**< the next file        */
     };
typedef struct io_File io_File;

struct io_Buffer {
       char *name;
       int  id;              /* used in helper array */
       int  grp_id;          /* id of the paremt  */ 
       io_File *file; 
       int  collective;      /* 1 for all pes to call each array */
                             /* 0 for not   */
    
       int  narrays;         /* for writing: total num in all buffers of this pe */ 
                             /* for reading: total number of arrays of this buffer.
                                If collective = 1, narrays is the same as the number
                                of arrays on any pe. But, if collective = 0, narrays
                                is the sum of the numbers of arrays of all pes. 
                             */ 
       int  max_narrays;     /* max numbers of arrays, used in datatypes */

       bio_Data_Type *datatypes;     /* data types of each array */ 

       int  max_nbufs; 
       int  *narrays_buf;    /*  narrays_buf is an array for 
                                 narrays of each nbufs_written */
       int  nbufs_written; 
       long long max_size_in_name;
       long long size_in_name;
       char *names;                      /* obj_type follows each name with szint */ 
       long long max_size_in_helper;
       long long *helper; /* (offset,size), (offset, size) ... for collective = 1 */
                          /* size, size, size, ....            for collective = 0 */ 
                          /* in terms of the datatype of each array  */
       long long max_size_in_data;   /* in the unit of bytes */
       long long size_in_data;       /* in the unit of bytes */ 
       long long tot_size_in_data; /* for collective = 1: all buffers, all pes, in bytes   */ 
                                   /* for collective = 0: all buffers of mype, in bytes    */
       void *data; 

       long long *helper_for_read; 

       /* for reading */

       int npes;
       long long **narrays_pe_buf;   /* to be freed, but not narrays_pe_buf[0] */ 
       long long **sizes_pe;     /* In the datatype of each array,
                                    for collective sizes_pe[0] is an array for pe 0,
                                    which contains (offset, size), (offset, size) ...
                                    for noncollective, sizes_pe[0] is ian array,
                                    which contains size, size, size, ....
                                  */ 
                                   /* sizes_pe is to be freed, but not sizes_pe[0] */ 
       char ***data_pe_buf;       /* data_pe_buf, data_pe_buf[k], data_pe_buf[k][b]
                                     are to be freed if they are not null. */ 

       bio_Data_Type ***datatypes_pe_buf; /* datatypes_pe_buf, datatypes_pe_buf[k], 
                                             to be freed */

       struct io_Buffer *next;
      };
typedef struct io_Buffer io_Buffer;  

static int io_nshared_in_buf = 16; /* npes,collective, narrays,nbufs,
                                   tot_size_in_data, rsved */

static int io_nshared_in_ibuf = 16;
 
static int io_2swap    = 0;
static int io_next_id  = 0;
static long long io_initial_size = 100000; /* size used in io_attr_init, io_object_init */ 

static char io_buffer_pre[] = "bio_buffer$";
static const int io_lbuffer_pre = 11; 
static char io_bufnames_pre[] = "bio_buf_names$";
static const int io_lbufnames_pre = 14; 
static char io_bufhelper_pre[] = "bio_buf_helper$";
static const int io_lbufhelper_pre = 15;

static char io_ibuffer_pre[] = "bio_ibuffer$";
static const int io_libuffer_pre = 12; 
static char io_ibufnames_pre[] = "bio_ibuf_names$";
static const int io_libufnames_pre = 15; 
static char io_ibufhelper_pre[] = "bio_ibuf_helper$";
static const int io_libufhelper_pre = 16;

static int  flag_for_resilience   = 0;
static char *io_filename_backup   = NULL; 
int  flag_for_buffered_read = 0; 


int io_read_hdrtail(int ifbackup, io_File *file, bio_IO_Mode io, io_Obj *obj, int *successful);
int io_read(io_File *file, io_Buffer *buffer,
            int grp_id, char *name, bio_Data_Type datatype,
            long long offset, long long size, void **array,
            int *array_id); 

int io_attr_init(int id, io_Attr *attr, long long max_size);

int io_clean_attr(io_Attr *a);

int io_object_init(io_Obj *s, int fileid, 
                    int max_num, long long max_size);

int io_clean_obj(io_Obj *s);

int io_init();

int io_attr_allocate(io_Attr *s, long long size);

int io_object_num_allocate(io_Obj *s);

int io_object_value_allocate(io_Obj *s, long long size);

int io_sizeof(bio_Data_Type d);

int io_add_object(int to_write, io_Tree **my_parent,
                  io_File *file, char *name, bio_Object_Type type,
                  bio_Data_Type datatype, int *id,  
                  long long offset, long long size, long long gsize,
                  void *buffer);

int io_add_grp(io_Tree **parent, io_File *file, char *name, int *grp_id,
               bio_Object_Type type, bio_Data_Type datatype);

int io_array_copy(void *src, void *target, long long size, 
                  bio_Data_Type type_s, bio_Data_Type type_t);

int io_get_array(io_File *file, io_Tree *ptree, char *name, bio_Data_Type datatype,
                 long long offset, long long size, 
                 void **buffer, int *array_id);

int io_add_file(char *name, bio_File_Mode mode, bio_IO_Mode io, io_File *file);

int io_add_attr(io_Attr *a, char *name,
                bio_Data_Type datatype, int size, void *buffer);

int io_get_attr(io_Attr *a, char *name, bio_Data_Type *datatype, int *num, void **buffer); 

int io_get_list(int count_only, io_File *file, io_Tree *ptree, int read_value_too,
                char *filter, int *n, bio_List_Struct **list);

int io_get_list_from_buffer(io_File *file, int count_only, io_Buffer *buffer, int read_value_too,
                          char *filter, int *n, bio_List_Struct **list);

char *mystrchar(char *s, char separator); 

int io_set_this_obj(io_Child *this_child, char *fullname, io_Tree *ptree,
                    io_File *file, int read_value_too, int *num, bio_List_Struct *obj); 

int io_set_this_buffer(io_File *file, int parentid, char *purename, int id, 
                         int collective, bio_List_Struct *list);

int io_get_attr_list(io_Attr *a, char *filter, int *num, bio_List_Struct **list);

int io_get_object(io_File *file, int id, int *obj_idx,
                  int *parent_id, bio_Object_Type *type, bio_Data_Type *datatype,
                  long long *offset, long long *gsize, char **name);

void io_list_null(int n, bio_List_Struct *list);
int io_file_null(io_File *file);
int io_buffer_null(io_Buffer *buffer);
int io_buffer_clean(io_Buffer *buffer);

int io_buffer_copy(io_File *file, io_Buffer *buffer, char *name, bio_Object_Type obj_type,
                   bio_Data_Type datatype, long long offset, long long size, void *array);

int io_buffer_commit(io_File *file, io_Buffer *buffer, int *commited);
int io_buffer_commit_nonc(io_File *file, io_Buffer *buffer, int *commited);

int io_buffer_finalize(io_File *file, io_Buffer *buffer);

int io_buffer_finalize_nonc(io_File *file, io_Buffer *buffer); 

int io_get_buffer_info(long long *helper,
                    long long ***my_narrays_pe_buf, long long ***my_offsets_pe,
                    char ****my_data_pe_buf, bio_Data_Type ****my_datatypes_pe_buf);

int io_get_buffer_info_nonc(long long *helper,
                    long long ***my_narrays_pe_buf, long long ***my_sizes_pe,
                    char ****my_data_pe_buf, bio_Data_Type ****my_datatypes_pe_buf);

int io_buffer_getname(char *name, char *purename, int *pe);

int io_offset_gsize_get(long long size, long long *offset, long long *gsize);
int io_calc_datatype_stord2new(void);

int io_convert_metadata(int fid_stord, io_Obj *obj, int num_obj, int id_adjust,
                        void *meta_data, long long size_meta_stord);

int io_construct_tree(io_File *file, io_Obj *obj);

int io_file_open(char *name, bio_File_Mode mode, bio_IO_Mode io, io_File **myfile); 
long long io_file_flush(io_File *file); 

io_Tree *io_find_file(int fileid, io_File **file);
io_Tree *io_find_tree(io_File *file, int obj_id);
io_Tree *io_find_tree_helper(int obj_id, io_Tree *tree); 

int io_clean_tree(io_Tree *ptree);

int io_swap(void *buffer, int unitsize, long long size);
int io_swap_metadata(void *metadata);

/**** end of io_int.h  *****/  


static const int io_version_num   = 1;
static const int io_hdr_size      = 864; /* in terms of char */ 
static const int io_nchars_in_hdr = 64;  /* the number of chars at the beginning of io_hdr */  
static char *io_hdr = NULL; 

static int io_my_version_num;            /* the version number in a file to be read */
static int io_my_hdr_size;        
static int io_case_num = 0;    

static int io_szchar;
static int io_szshort; 
static int io_szint;
static int io_szlong;
static int io_szllong;
static int io_szfloat;
static int io_szdouble;
static int io_szldouble;

static int io_szint_stord;
static int io_szlong_stord;
static int io_szllong_stord;
static int io_szfloat_stord;
static int io_szdouble_stord;
static int io_szldouble_stord;

static io_Buffer *io_buffers = NULL;
static io_File *io_files = NULL; 
static int io_num_file  = 0;
static int mype = 0;
static int npes = 1; 

/* static size_t io_sizeof[1 + (int)bio_datatype_invalid - (int) bio_char]; */ 
static bio_Data_Type io_datatype_stord2new[1 + (int)bio_datatype_invalid - (int) bio_char];

#ifdef MPI 
static MPI_Comm bio_comm = MPI_COMM_WORLD;
static MPI_Status  *mpi_status = NULL;
static MPI_Request *mpi_reqs = NULL;
#endif 
#ifdef MPI
static MPI_Info bio_info = MPI_INFO_NULL;
#endif

/****************************************************************************************/
int bio_comm_create(int n, int *ranks)
{
#ifdef MPI   
    int i, k, ncpus, myrank;
    MPI_Group world_comm_group,comm_group;

    MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
    if ((n == 0) || (n == ncpus)) {
        MPI_Comm_dup(MPI_COMM_WORLD, &bio_comm);
    } 
    else { 
        if (!ranks)
           { printf("ERROR: null ranks in bio_comm_create\n");
             return -1;
            }
        if (n < 0)
           { printf("ERROR: n < = 0 in bio_comm_create\n");
             return -1;
            }
        if (n > ncpus)
           { printf("ERROR: n larger than the number of pes allocated in bio_comm_create\n");
             return -1;
            }
        for (i = 0; i < n; i++)
          { if (ranks[i] >= ncpus)
               { printf("ERROR: ranks[%d] = %d too large in bio_comm_create\n", i, ranks[i]);
                 return -1;
                }
           }
        for (i = 0; i < n; i++)
          { myrank = ranks[i];
            for (k = i+1; k < n; k++)
              { if (ranks[k] == myrank)
                   { printf("ERROR: ranks[%d] = ranks[%d] in bio_comm_create\n", k, i);
                     i = n+1;
                     break;
                    }
               }
           }
        if (bio_comm != MPI_COMM_WORLD) { 
            MPI_Comm_free(&bio_comm);
            bio_comm = MPI_COMM_WORLD;
        }
        MPI_Comm_group(MPI_COMM_WORLD, &world_comm_group);
        MPI_Group_incl(world_comm_group, n, ranks, &comm_group);
        MPI_Comm_create(MPI_COMM_WORLD, comm_group, &bio_comm);
    
        MPI_Group_free(&comm_group);
        MPI_Group_free(&world_comm_group);
    } 
#endif 
    return 0;
} 
/****************************************************************************************/
int bio_get_next_id()
{  
    int id;
    id = io_next_id;
    io_next_id++;
    return id;
 }
/****************************************************************************************/
#ifdef MPI
int bio_get_comm(MPI_Comm *comm)
{ 
    *comm = bio_comm;
    return 0;
 } 
#endif
/****************************************************************************************/
int bio_comm_free()
{ 
#ifdef MPI 
    if (bio_comm != MPI_COMM_WORLD)
       { MPI_Comm_free(&bio_comm);
         bio_comm = MPI_COMM_WORLD;
        }
#endif 
    return 0;
 } 
/****************************************************************************************/
int bio_info_set(char *keyword, char *value)
{ 
#ifdef MPI 
    if (!keyword)
       { printf("ERROR: null keyword in bio_info_set\n");
         return -1;
        }
    if (!value)
       { printf("ERROR: null value in bio_info_set\n");
         return -1;
        }
    if (bio_info == MPI_INFO_NULL)
       { if (MPI_Info_create(&bio_info) != MPI_SUCCESS)
            { printf("ERROR: MPI_Info_create failed\n");
              return -1;
             }
        }
    if (MPI_Info_set(bio_info, keyword, value) != MPI_SUCCESS)
       { printf("ERROR: MPI_Info_set failed\n");
         return -1;
        }
#endif 
    return 0;
   }
/****************************************************************************************/
int bio_info_delete(char *keyword)
{ 
#ifdef MPI 
    if (!keyword)
       { printf("ERROR: null keyword in bio_info_delete\n");
         return -1;
        }
    if (MPI_Info_delete(bio_info, keyword) != MPI_SUCCESS)
       { printf("ERROR: MPI_Info_delete failed\n");
         return -1;
        }
#endif
    return 0;
   }
/****************************************************************************************/
int bio_info_free()
{   
#ifdef MPI 
    if (bio_info != MPI_INFO_NULL)
       { if (MPI_Info_free(&bio_info) != MPI_SUCCESS)
            { printf("ERROR: MPI_Info_free failed\n");
              return -1;
             }
         bio_info = MPI_INFO_NULL; 
        }
#endif
    return 0;
   }
/****************************************************************************************/
int io_attr_init(int id, io_Attr *attr, long long max_size)
{
   attr->max_size = max_size;
   attr->value = malloc((size_t)(max_size + 3 * io_szllong));
   assert(attr->value); 
/*
   size,            szllong 
   id,              szllong 
   num,             szllong  
   size of name,    szllong 
   name,            slen*szchar
   datatype,        szllong 
   size of value,   szllong 
   values,          n*io_sizeof(datatype) 
*/
   attr->size = (long long *)attr->value;
   attr->id   = attr->size + 1;
   attr->num  = attr->id + 1;
   
   *(attr->size) =  3 * io_szllong; /* size, id, num, size_name, size_value */ 
   *(attr->id)   = id;
   *(attr->num)  =  0;

   return 0;
  }
/****************************************************************************************/
int io_object_init(io_Obj *s, int fileid, int max_num, long long max_size)
{ 
    int n, i;
    long long a4[4];
    void *p; 

    s->max_num = max_num;
    s->max_size = max_size; 

    s->value = malloc((size_t)(max_size  + 4 * io_szllong)); 
    assert(s->value);
/*
    size_meta,          szllong 
    fid,                szllong   
    num,                szllong 
    size_data,          szllong
    size of name,       szllong
    name,               slen*szchar
    id,                 szllong 
    parent_id,          szllong  
    obj_type,           szllong 
    datatype,           szllong  
    offset,             szllong
    size,               szllong
*/   

    s->size_meta = (long long *) s->value;
    s->fileid   = s->size_meta + 1;
    s->num       = s->fileid   + 1;
    s->size_data = s->num + 1;

/****
    *(s->size_meta) = 4 * io_szllong;  
    *(s->fileid)   = fileid;
    *(s->num)       = 0;
    *(s->size_data) = 0;
****/

    a4[0] = 4 * io_szllong;
    a4[1] = fileid;
    a4[2] = 0;
    a4[3] = 0;
    memcpy(s->value, a4, (size_t)(4 * io_szllong));

    s->attr  = (io_Attr **) malloc((size_t)max_num * sizeof(io_Attr*));
    for (i = 0; i < max_num; i++)
      { s->attr[i] = NULL;
       }
    return 0;
   }  
/****************************************************************************************/
int io_init()
{   int i, n;
    long long *p;

    io_szchar = sizeof(char);
    io_szshort = sizeof(short int); 
    io_szint  = sizeof(int);
    io_szlong = sizeof(long);
    io_szllong   = sizeof(long long);
    io_szfloat   = sizeof(float);
    io_szdouble  = sizeof(double);
    io_szldouble = sizeof(long double);

/***
    i = (int)bio_char;
    io_sizeof[0] = 1;
    n = (int)bio_double - i;
    io_sizeof[n] = io_szdouble;
    n = (int)bio_float - i;
    io_sizeof[n] = io_szfloat;
    n = (int)bio_int - i;
    io_sizeof[n] = io_szint;
    n = (int)bio_long - i;
    io_sizeof[n] = io_szlong;
    n = (int)bio_long_double - i;
    io_sizeof[n] = io_szldouble;
    n = (int)bio_long_long - i;
    io_sizeof[n] = io_szllong;
***/

    for (i = 0; i < io_nchars_in_hdr; i++) 
      { io_hdr[i] = 'a';
       }
    p = (long long *)(io_hdr + io_nchars_in_hdr);
    n = (io_hdr_size - io_nchars_in_hdr)/io_szllong;
    for (i = 0; i < n; i++) 
      { p[i] = 0;
       }  
    io_my_version_num = io_version_num;

    return 0;
  }  
/****************************************************************************************/
int io_attr_allocate(io_Attr *s, long long size)
{ 
    long long max_size;
    void *v;
 
    if (size > s->max_size) { 
        max_size = size + size;
    }
    else { 
        max_size = s->max_size + s->max_size;
    }
    v = malloc((size_t)(max_size + 3 * io_szllong));

    memcpy(v, s->value, (size_t)(*(s->size)));
    
    s->max_size = max_size;
    s->size = (long long *) v;
    s->id   = s->size + 1;
    s->num  = s->id   + 1;

    free(s->value);
    s->value = v;

    return 0;
   } 
/****************************************************************************************/
int io_object_num_allocate(io_Obj *s)
{   
    long long i, max_num, num;
    io_Attr **attr;

    max_num = s->max_num + s->max_num;
    num = *(s->num); 

    attr = (io_Attr **) malloc((size_t)max_num * sizeof(io_Attr *));
    for (i = 0; i < num; i++)
      { attr[i] = s->attr[i];
       }
    for (i = num; i < max_num; i++)
      { attr[i] = NULL;
       } 
    free(s->attr);
    s->attr = attr;
    s->max_num = max_num;

    return 0;
   }
/****************************************************************************************/
int io_object_value_allocate(io_Obj *s, long long size)
{ 
    long long max_size;
    void *v;

    if (size > s->max_size) { 
        max_size = size + size;
    }
    else { 
        max_size = s->max_size + s->max_size;
    }
    v = malloc((size_t)(max_size + 4 * io_szllong));

    memcpy(v, s->value, (size_t)(*(s->size_meta)));
    free(s->value);
    s->value = v;
    s->max_size = max_size; 

/*  re-points to the new locations */ 

    s->size_meta = (long long *) s->value;
    s->fileid   = s->size_meta + 1;
    s->num       = s->fileid + 1;
    s->size_data = s->num + 1;

    return 0;
   }
/****************************************************************************************/
int io_sizeof(bio_Data_Type d)
{   if (d == bio_char) 
       return io_szchar;
    else if (d == bio_int)
       return io_szint;
    else if (d == bio_long)
       return io_szlong;
    else if (d == bio_long_long)
       return  io_szllong;
    else if (d == bio_float) 
       return io_szfloat;
    else if (d == bio_double)
       return io_szdouble;
    else if (d == bio_long_double)
       return  io_szldouble;
    else 
        return 1LL;
    
    return 1LL;
   } 
/****************************************************************************************/
int io_add_object(int to_write, io_Tree **my_parent,
                  io_File *file, char *name, bio_Object_Type type,
                  bio_Data_Type datatype, int *id, 
                  long long offset, long long size, long long gsize,
                  void *buffer)
{ /** 
    If to_write = 0,  this routine only updates it database of this pe without actual
                      writing, one or more other pes may actually write data to disk.
    If to_write != 0, all pes writes.

    If *id >= 0, *id is used the id of the dataset. 
    If *id <  0, the next available id will be used. 
  **/   

    int  i, mytype, myslen, slen, num, s, slen_stord;
    int  sztype, id_stord, parentid;
    long long size_meta_added, gsize_stord, mysize, myoffset, size_written;
    long long a6[6];
    char *pc, *myname, *c, *ptr;
    long long *pll, *po, *pd, *pid, *poffset;
    bio_Data_Type datatype_stord;
    bio_Object_Type  type_stord;
    io_Attr *a; 
    io_Obj *obj; 
    io_Tree *parent, *tree;
    io_Child *child, *next;

#ifdef MPI 
    MPI_Datatype mpidatatype;
    MPI_Status status;
#endif 

    if (to_write && gsize && (type != bio_group)) { 
        file->mbytes += (0.000001 * (double)(gsize * io_sizeof(datatype)));
    } 
    child = NULL;
    slen = strlen(name);
    obj = file->object; 
    num = *(obj->num);

    if (my_parent) { 
        parent = *my_parent;
        if (parent) { 
            parentid = parent->id;
        }
        else { 
            parentid = -1;
        } 
    }
    else { 
        parent = NULL;
        parentid = -1;
    }
    if (parent && to_write) {

/*      check whether or not this array already exists under the group */

        child = parent->child;
        while (child) { 
           tree = child->tree;
           if (!tree) {
               printf("ERROR: null tree in io_add_object\n");
               return -1;
           }
           c = (char *)tree->p;
           memcpy(a6, c, (size_t)io_szllong);
           myslen = a6[0];
           myname = c + io_szllong;
           c = myname + myslen;
           if (slen == myslen && !strncmp(name, myname, (size_t)myslen)) {
               memcpy(a6, c, (size_t)(6 * io_szllong));
               id_stord = (int) a6[0];
               type_stord = (bio_Object_Type) ((int)a6[2]);

               if (type_stord != type) {
                   printf("ERROR: type incorrect in io_add_object\n");
                   return -1;
               }
               datatype_stord = (bio_Data_Type) ((int)a6[3]);
               mysize = a6[5];
               sztype = io_sizeof(datatype_stord);
               gsize_stord = mysize/sztype;

               *id = id_stord;
               if (file->id_min > *id) { 
                   file->id_min = *id;
               }
               break;
            }
            next = child->next;
            child = next;
        }
        if (child && type != bio_group) {   
            if (datatype_stord != datatype ||  gsize_stord != gsize) {  
                printf("ERROR: rewrite the same array = %s, but different type/gsize\n", name);
                return -1;
            }
            myoffset = (long long)io_my_hdr_size + a6[4] + offset * sztype;
#ifdef MPI 
            if (datatype_stord == bio_char) { 
                mpidatatype = MPI_CHAR;
            }
            else if (datatype_stord == bio_int) {
                mpidatatype = MPI_INT;
            }
            else if (datatype_stord == bio_long) {
                mpidatatype = MPI_LONG;
            }
            else if (datatype_stord == bio_long_long) {
                mpidatatype = MPI_LONG_LONG_INT;
            }
            else if (datatype_stord == bio_float) {
                mpidatatype = MPI_FLOAT;
            }
            else if (datatype_stord == bio_double) {
                mpidatatype = MPI_DOUBLE;
            }
            else if (datatype_stord == bio_long_double) {
                mpidatatype = MPI_LONG_DOUBLE;
            }
            else { 
                printf("ERROR: datatype_stord incorrect in io_add_object\n");
                return -1;
            } 
            if (file->io == bio_collective) {  
               /****
               MPI_File_seek(file->mfile, (MPI_Offset)myoffset, MPI_SEEK_SET);
               MPI_File_write_all(file->mfile, buffer, (int)size, mpidatatype, &status);
               ***/
               MPI_File_write_at_all(file->mfile, (MPI_Offset)myoffset,buffer, (int)size, mpidatatype, &status);
            }
            else if (file->io == bio_independent) {   
               MPI_File_write_at(file->mfile,(MPI_Offset)myoffset,buffer,(int)size,mpidatatype,&status); 
            }
            else if (file->io == bio_use_fopen) {
               fseek(file->fp, (long)myoffset, SEEK_SET);
               size_written = (long long)fwrite(buffer, (size_t)sztype, (size_t)size, file->fp);
               if (size_written != size) { 
                   printf("ERROR:(size_written != szie in io_add_object\n");
                   return -1;
               }
            } 
            else if (file->io == bio_use_open) {
               size_written = (long long) pwrite(file->fd, buffer, (size_t)(size*sztype),(off_t)myoffset);
               if (size_written != size * (long long)sztype) { 
                   printf("ERROR:(size_written != szie in io_add_object\n");
                   return -1;
               }
            }             
#else
            if ((file->io == bio_use_fopen) || (file->io == bio_collective) || (file->io == bio_independent)) {
                fseek(file->fp, (long)myoffset, SEEK_SET);
                size_written = (long long) fwrite(buffer, (size_t)sztype, (size_t)size, file->fp);
                if (size_written != size) { 
                    printf("ERROR:(size_written != szie in io_add_object\n");
                   return -1;
                }
            }
            else if (file->io == bio_use_open) {
               size_written = (long long) pwrite(file->fd, buffer, (size_t)(size * sztype), (off_t)myoffset);
               if (size_written != size * (long long) sztype) { 
                   printf("ERROR:(size_written != szie in io_add_object\n");
                   return -1;
               }
            }             
#endif 
            return 0;
        } 
    }        
    if (!child) {  
        if (*id < 0) {
           *id = io_next_id;
            io_next_id++;
        } 
        mytype = (int) type; 
        size_meta_added = slen + 7 * io_szllong;  
        mysize = size_meta_added + *(obj->size_meta);
        if (mysize > obj->max_size) { 
            io_object_value_allocate(obj, mysize);
            io_construct_tree(file, obj); 
        /*  since obj->root has been re-created. Original parent was deallocated. */

            if (parentid >= 0) { 
                parent = io_find_tree(file, parentid);
                *my_parent = parent;
            }
            else { 
                printf("ERROR: parentid < 0 in io_add_object\n");
                return 0;
            }
        } 
        if ((long long) num == obj->max_num) {  
            io_object_num_allocate(obj);
        }
/*      In the order of
        size of name,       szllong 
        name,               slen*szchar
        id,                 szllong 
        parent_id,          szllong 
        obj_type,           szllong 
        datatype,           szllong 
        offset,             szllong
        size,               szllong
*/ 

/****
        pll  = (long long *)((char *)obj->value + *(obj->size_meta));
        *pll = (long long) slen;
        pc   = (char *)(pll + 1);
        memcpy(pc, name, (size_t)slen);
        pll = (long long *)(pc + slen);
        pll[0] = *id;
        pll[1] = parentid;
        po  = pll + 2;
        *po = (long long) mytype;
        pd  = po + 1;
        *pd = (long long) datatype;
        sztype = io_sizeof(datatype);
        pll = pd + 1;
        pll[0] = *(obj->size_data);
        pll[1] = gsize * sztype;
***/
        sztype = io_sizeof(datatype);
        mysize = (long long) slen;
        pc = (char *)obj->value + *(obj->size_meta);
        memcpy(pc, &mysize, (size_t)io_szllong);

        ptr = pc;

        pc += io_szllong;
        memcpy(pc, name, (size_t)slen);
        pc += slen;
        a6[0] = *id;
        a6[1] = parentid;
        a6[2] = (long long) mytype;
        a6[3] = (long long) datatype;
        a6[4] = *(obj->size_data);
        if (type == bio_group) { 
            a6[5] = 0;
        }
        else { 
            a6[5] = gsize * sztype; 
        }
        memcpy(pc, a6, (size_t)(6 * io_szllong));
    
        if ((type != bio_group) && to_write) {  

            myoffset = (long long)io_hdr_size + *(obj->size_data) + offset * sztype;
#ifdef MPI
            if (datatype == bio_char) {
                mpidatatype = MPI_CHAR;
            }
            else if (datatype == bio_int) {
                mpidatatype = MPI_INT;
            }
            else if (datatype == bio_long) {
                mpidatatype = MPI_LONG;
            }
            else if (datatype == bio_long_long) {
                mpidatatype = MPI_LONG_LONG_INT;
            }
            else if (datatype == bio_float) {
                mpidatatype = MPI_FLOAT;
            }
            else if (datatype == bio_double) {
                mpidatatype = MPI_DOUBLE;
            }
            else if (datatype == bio_long_double) {
                mpidatatype = MPI_LONG_DOUBLE;
            }
            else {
                printf("ERROR: datatype_stord incorrect in io_add_object\n");
                return -1;
            }
            if (file->io == bio_collective) {  
               /****
               MPI_File_seek(file->mfile, (MPI_Offset)myoffset, MPI_SEEK_SET);
               MPI_File_write_all(file->mfile, buffer, (int)size, mpidatatype, &status);
               ***/
               MPI_File_write_at_all(file->mfile, (MPI_Offset)myoffset, buffer, (int)size, mpidatatype, &status);
            } 
            else if (file->io == bio_independent) {  
               MPI_File_write_at(file->mfile,(MPI_Offset)myoffset,buffer,(int)size,mpidatatype,&status);
            }
            else if (file->io == bio_use_fopen) { 
               fseek(file->fp, (long) myoffset, SEEK_SET);
               size_written = (long long) fwrite(buffer, (size_t)sztype, (size_t)size, file->fp);
               if (size_written != size) { 
                   printf("ERROR: size_written != size in io_add_object\n");
                   return -1;
               }
            } 
            else if (file->io == bio_use_open) { 
               size_written = (long long) pwrite(file->fd, buffer, (size_t)(size * sztype), (off_t)myoffset);
               if (size_written != size * (long long) sztype) { 
                   printf("ERROR:(size_written != szie in io_add_object\n");
                   return -1;
               }
            } 
#else
            if ((file->io == bio_use_fopen) || (file->io == bio_collective) || (file->io == bio_independent)) { 
                fseek(file->fp, (long) myoffset, SEEK_SET);
                size_written = fwrite(buffer, (size_t)sztype, (size_t)size, file->fp);
                if (size_written != size) { 
                   printf("ERROR: size_written != size in io_add_object\n");
                   return -1;
                }
            }
            else if (file->io == bio_use_open) { 
               size_written = (long long) pwrite(file->fd, buffer, (size_t)(size * sztype), (off_t)myoffset);
               if (size_written != size * (long long)sztype) { 
                   printf("ERROR:(size_written != szie in io_add_object\n");
                   return -1;
               }
            } 
#endif 
        }
        a = (io_Attr *) malloc(sizeof(io_Attr));
        io_attr_init(*id, a, io_initial_size);
        obj->attr[num] = a;
    
        *(obj->size_meta) += size_meta_added;
        *(obj->num) += 1; 
        if (type != bio_group) {  
            *(obj->size_data) += a6[5];
        } 
        tree = (io_Tree *) malloc(sizeof(io_Tree));
        tree->id = *id;
        tree->idx = num;
        tree->p = ptr;
        tree->parent = parent;
        tree->child = NULL;
        if (parent) {
            child = (io_Child *) malloc(sizeof(io_Child));
            child->tree = tree;
            child->next = parent->child;
            parent->child = child;
        }
        else { 
            file->root = tree;
        } 
        if (file->id_min > *id) { 
            file->id_min = *id;
        }
    } 
    return 0;
  } 
/****************************************************************************************/
#define IO_ARRAY_COPY(type_t, type_s) \
{                            \
type_s *s;                   \
type_t *t;                   \
s = (type_s *) src;       \
t = (type_t *) target;       \
for (i = 0; i < size; i++) {   \
    t[i] = (type_t) s[i];    \
}                            \
}
int io_array_copy(void *src, void *target, long long size, bio_Data_Type type_s, bio_Data_Type type_t)
{
   // This is only for different types   
 
   long long i;
   char *tc, *sc;
   int  *ti, *si;
   long *tl, *sl;
   long long *tll, *sll;
   float *tf, *sf;
   double *td, *sd;
   long double *tld, *sld;
   
    if (type_s == bio_datatype_invalid || type_t == bio_datatype_invalid)  { 
        printf("ERROR: type_s or type_t incorrect in io_array_copy\n");
        return -1;
    } 
    if (type_t == bio_char)  {  
       if (type_s == bio_int) {  
           IO_ARRAY_COPY(char, int);
       }
       else if (type_s == bio_long)  { 
           IO_ARRAY_COPY(char, long);
       }
       else if (type_s == bio_long_long) { 
           IO_ARRAY_COPY(char, long long);
       }
    }
    else if (type_t == bio_int) { 
       if (type_s == bio_char) { 
           IO_ARRAY_COPY(int, char);
       }
       else if (type_s == bio_long) { 
           IO_ARRAY_COPY(int, long);
       } 
       else if (type_s == bio_long_long) {   
           IO_ARRAY_COPY(int, long long);
       }
       else if (type_s == bio_float) {  
           IO_ARRAY_COPY(int, float);
       }
       else if (type_s == bio_double) {  
           IO_ARRAY_COPY(int, double);
       }
       else if (type_s == bio_long_double) {  
           IO_ARRAY_COPY(int, long double);
       }   
    } 
    else if (type_t == bio_long)  {  
       if (type_s == bio_char)  { 
           IO_ARRAY_COPY(long, char);
       }
       else if (type_s == bio_int) { 
           IO_ARRAY_COPY(long, int);
       }
       else if (type_s == bio_long_long)  {  
           IO_ARRAY_COPY(long, long long);
       }
       else if (type_s == bio_float) {  
           IO_ARRAY_COPY(long, float);
       }
       else if (type_s == bio_double) {  
           IO_ARRAY_COPY(long, double);
       }
       else if (type_s == bio_long_double) {  
           IO_ARRAY_COPY(long, long double);
       } 
    } 
    else if (type_t == bio_long_long) {  
       if (type_s == bio_char) { 
          IO_ARRAY_COPY(long long, char);
       }
       else if (type_s == bio_int) {  
          IO_ARRAY_COPY(long long, int);
       }
       else if (type_s == bio_long) { 
          IO_ARRAY_COPY(long long, long);
       }
       else if (type_s == bio_float)  {  
          IO_ARRAY_COPY(long long, float);
       }
       else if (type_s == bio_double) { 
          IO_ARRAY_COPY(long long, double);
       }
       else if (type_s == bio_long_double) {  
          IO_ARRAY_COPY(long long, long double);
       }
   } 
   else if (type_t == bio_float) {  
       if (type_s == bio_int) { 
           IO_ARRAY_COPY(float, int);
       }
       else if (type_s == bio_long) { 
           IO_ARRAY_COPY(float, long);
       }
       else if (type_s == bio_long_long)  {  
           IO_ARRAY_COPY(float, long long);
       }
       else if (type_s == bio_double) {  
           IO_ARRAY_COPY(float, double);
       }
       else if (type_s == bio_long_double)  {  
           IO_ARRAY_COPY(float, long double);
       }
   }
   else if (type_t == bio_double)  { 
       if (type_s == bio_int)  {  
           IO_ARRAY_COPY(double, int);
       }
       else if (type_s == bio_long) { 
           IO_ARRAY_COPY(double, long);
       }
       else if (type_s == bio_long_long) {  
           IO_ARRAY_COPY(double, long long);
       }
       else if (type_s == bio_float)  {  
           IO_ARRAY_COPY(double, float);
       }
       else if (type_s == bio_long_double)  {  
           IO_ARRAY_COPY(double, long double);
       }
   }
   else if (type_t == bio_long_double) {  
       if (type_s == bio_int)  { 
           IO_ARRAY_COPY(long double, int);
       } 
       else if (type_s == bio_long)  {  
           IO_ARRAY_COPY(long double, long);
       } 
       else if (type_s == bio_long_long)  { 
           IO_ARRAY_COPY(long double, long long);
       }
       else if (type_s == bio_float)  {  
           IO_ARRAY_COPY(long double, float);
       }
       else if (type_s == bio_double) {  
           IO_ARRAY_COPY(long double, double);
       }
   }
   return 0;
  } 
/****************************************************************************************/
int io_get_array(io_File *file, io_Tree *ptree, char *name, bio_Data_Type datatype,    
                 long long offset, long long size, void **buffer, int *array_id)
{ 
    int  m, num, i, slen, slen_stord, sztype_stord, parentid_stord, found;
    long long llsize, gsize, offset_stord, myoffset, size_read; 
    long long a6[6];
    bio_Data_Type datatype_stord; 
    bio_Object_Type type_stord;
    io_Obj *obj;
    io_Tree *tree;
    io_Child *child, *next;
    char *c, *myname;
    void *v;

#ifdef MPI
    MPI_Datatype mpidatatype;
    MPI_Status status;
#endif 
    
    m = 7 * io_szllong; 
         /* sz_name, id, parent_id, obj_type, datatype, offset, size */ 
    obj = file->object;
    num = *(obj->num);
    
    slen = strlen(name);

    if (!ptree) { 
        printf("NULL ptree in io_get_array\n");
        return -1;
    } 
    child = ptree->child;
    while (child) { 
       tree = child->tree;
       c = (char *) tree->p;
       memcpy(a6, c, (size_t)io_szllong);
       slen_stord = a6[0];
       myname = c + io_szllong;
       c = myname + slen_stord;
       memcpy(a6, c, (size_t)(6 * io_szllong));
       type_stord = (bio_Object_Type) ((int)a6[2]);

       if ((type_stord != bio_group) && 
           (slen == slen_stord)      &&
           !strncmp(name, myname, (size_t)slen_stord)) {

           *array_id = (int) a6[0];
           datatype_stord = (bio_Data_Type) ((int)a6[3]);
           offset_stord = a6[4];
           sztype_stord = io_sizeof(datatype_stord);
           gsize = a6[5]/sztype_stord;
        
           break;
       }
       next = child->next;
       child = next;
    } 
    if (!child) { 
        printf("Warning: no %s found in io_get_array\n", name);
        return -2;
    } 
    if (offset + size > gsize)
       { printf("ERROR: offset + size > total size of the array in io_get_array\n");
         return -1;
        }
    if (datatype == datatype_stord)
       { if (*buffer)
            { v = *buffer;
             }
         else 
            { llsize = size * sztype_stord;
              v = malloc((size_t)llsize);
              assert(v);
              *buffer = v;
             }
         }
     else 
        { llsize = size * sztype_stord;
          v = malloc((size_t)llsize);
          assert(v);
          if (!(*buffer))
             { llsize = size * io_sizeof(datatype);
               *buffer = malloc((size_t)llsize);
               assert(*buffer);
              }
         } 
    /* mysize = size * sztype_stord;
    */
    myoffset = (long long)io_my_hdr_size + offset_stord + offset * sztype_stord;
#ifdef MPI 
    if ((file->io == bio_collective) || (file->io == bio_independent)) {
        if (datatype_stord == bio_char) {
            mpidatatype = MPI_CHAR;
        }
        else if (datatype_stord == bio_int) {
            mpidatatype = MPI_INT;
        }
        else if (datatype_stord == bio_long) {
            mpidatatype = MPI_LONG;
        }
        else if (datatype_stord == bio_long_long) {
            mpidatatype = MPI_LONG_LONG_INT;
        }
        else if (datatype_stord == bio_float) {
            mpidatatype = MPI_FLOAT;
        }
        else if (datatype_stord == bio_double) {
            mpidatatype = MPI_DOUBLE;
        }
        else if (datatype_stord == bio_long_double) {
            mpidatatype = MPI_LONG_DOUBLE;
        }
        else {  
            printf("ERROR: datatype_stord incorrect in io_get_array\n");
            return -1;
        } 
        MPI_File_read_at(file->mfile,(MPI_Offset)myoffset,v,(int)size,mpidatatype,&status);
    }
    else if (file->io == bio_use_fopen) { 
        fseek(file->fp, (long)myoffset, SEEK_SET);
        size_read = (long long) fread(v, (size_t)sztype_stord, (size_t)size, file->fp); 
        if (size_read != size) { 
            printf("ERROR: size_read != size in io_get_array at A, size_read = %lld, size_requested = %lld\n",
                           size_read, size);
            return -1;
        } 
        else if (file->io == bio_use_open) { 
            size_read = (long long) pread(file->fd, v, (size_t)(size * sztype_stord), (off_t) myoffset);
            if (size_read != size * (long long)sztype_stord) { 
                printf("ERROR: size_read != size in io_get_array at B, size_read = %lld, size_requested = %lld\n",   
                         size_read, size);
                return -1;
            }
        } 
    }
#else 
    if ((file->io == bio_use_fopen) || (file->io == bio_collective) || (file->io == bio_independent)) { 
        fseek(file->fp, (long)myoffset, SEEK_SET);
        size_read = (long long) fread(v, (size_t)sztype_stord, (size_t)size, file->fp); 
        if (size_read != size) { 
            printf("ERROR: size_read != size in io_get_array at C, size_read = %lld, size_requested = %lld\n",   
                         size_read, size);
            return -1;
        } 
    }
    else if (file->io == bio_use_open) { 
        size_read = (long long) pread(file->fd, v, (size_t)(size * sztype_stord), (off_t) myoffset);
        if (size_read != size * (long long)sztype_stord) { 
            printf("ERROR: size_read != size in io_get_array at D, size_read = %lld, size_requested = %lld\n",   
                         size_read, size);
            return -1;
        }
    } 
#endif 
    if (io_2swap && datatype_stord != bio_char) { 
        io_swap(v, (int)sztype_stord, (long long) size);
    } 
    if (datatype != datatype_stord)
       { io_array_copy(v, *buffer, size, datatype_stord, datatype);
         free(v);
        } 
    return 0;
   } 
/****************************************************************************************/
int io_add_file(char *name, bio_File_Mode mode, bio_IO_Mode io, io_File *file)
{  
    long long llversion_num;
    long       lversion_num;
    int        iversion_num;
    short      sversion_num;

    char pathname[512], command[512], *myname;
    char achar;
    int  successful; 
    int i, k, n, id, parent_id, slen, num_obj, num_attr, fid_stord, fd;
    int my_id_min, my_id_max, id_adjust, m;
    int n1, n2, n3, n4, n5, n6;
    char *pc, *c;
    short int a5s[5];
    int a5i[5];
    long a5l[5];
    long long a5ll[5];
    void *pv;
    int  *pi;
    long *pl;
    short int *ps;
    long long *pll;
    bio_Data_Type datatype_i, datatype_ll;
    void *meta_data;
    long long offset, myoffset, size, gsize;
    long long size_data, size_meta_block, size_meta_stord, size_attr; 
    io_Attr *a; 
    io_Obj  *obj;

#ifdef MPI 
    MPI_Status status;
    file->t_open = MPI_Wtime();  
#else 
    time(&(file->t_open));
#endif 

    if (mode == bio_file_create) {
        strcpy(pathname, name);
        slen = (int) strlen(pathname);
        n = 0;
        for (i = 0; i < slen; i++) { 
            if (pathname[i] == '/') { 
                pathname[i] = '\0';
                n++;
            }
        }
        if (n) { 
            slen = (int) strlen(pathname);
            myname = pathname;
            for (i = 0; i < n; i++) { 
                 if (slen == 0) { 
                     sprintf(command, "cd /");
                     system(command);
                     myname++;
                     continue;
                 } 
                 if (slen == 1) { 
                     if (myname[0] == '.') { 
                         myname += (slen + 1);
                         slen = (int) strlen(myname);
                         continue;
                     }
                 } 
                 else if (slen == 2) { 
                     if (myname[0] == '.' && myname[1] == '.') { 
                         sprintf(command, "cd ..");
                         system(command);
                         myname += (slen + 1);
                         slen = (int) strlen(myname);
                         continue;
                     }
                 } 
                 sprintf(command, "/bin/mkdir -p %s", myname);
                 system(command);
                 sprintf(command, "cd %s", myname);
                 myname += (slen + 1);
                 slen = (int) strlen(myname);
            }
        } 
    } 
    if (mode == bio_file_create) {  
#ifdef MPI 
        if ((io == bio_collective) || (io == bio_independent)) { 
            MPI_File_open(bio_comm, name, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                      bio_info, &(file->mfile));
        }
        else if (io == bio_use_fopen) { 
            if (mype == 0) { 
                file->fp = fopen(name, "wb");
            }
            if (npes > 1) MPI_Barrier(bio_comm);
            if (mype) { 
                file->fp = fopen(name, "r+b");
            }
            if (npes > 1) MPI_Barrier(bio_comm);
        } 
        else if (io == bio_use_open) { 
            if (mype == 0) { 
                file->fd = open(name, O_WRONLY|O_CREAT|O_TRUNC|O_CONCURRENT_WRITE, 0644);
            } 
            if (npes > 1) MPI_Barrier(bio_comm);

            if (mype != 0) { 
                file->fd = open(name, O_WRONLY|O_CONCURRENT_WRITE,0644);
            }  
            if (npes > 1) MPI_Barrier(bio_comm);
            if (file->fd == -1) { 
                printf("ERROR: open file at 1 %s\n", name);
                return -1;
            } 
        } 
#else 
        if (npes == 1) { 
            if ((io == bio_use_fopen) || (io == bio_collective) || (io == bio_independent)) { 
                file->fp = fopen(name, "wb");
            }
            else if (io == bio_use_open) { 
                file->fd = open(name, O_WRONLY|O_CREAT|O_TRUNC|O_CONCURRENT_WRITE, 0644);
                if (file->fd == -1) { 
                    printf("ERROR: open file at 2 %s\n", name);
                    return -1;
                } 
            }
        }
        else { 
            if (mype == 0) { 
                printf("ERROR: MPI not defined in the compile line\n");
            }
            return -1;
        }
#endif 
    }
    else if (mode == bio_file_read_only) {  
#ifdef MPI 
        if ((io == bio_collective) || (io == bio_independent)) {
            MPI_File_open(bio_comm,name,MPI_MODE_RDONLY,bio_info,&(file->mfile));
        }
        else if (io == bio_use_fopen) { 
            file->fp = fopen(name, "rb"); 
            if (!(file->fp)) { 
                printf("ERROR: open %s\n", name);
                return -1;
            }
        }
        else if (io == bio_use_open) { 
            file->fd = open(name, O_RDONLY|O_CONCURRENT_WRITE, 0);
            if (file->fd == -1) { 
                printf("ERROR: open file at 3 %s\n", name);
                return -1;
            }
        }
#else 
        if ((io == bio_use_fopen) || (io == bio_collective) || (io == bio_independent) ) { 
            file->fp = fopen(name, "rb"); 
            if (!(file->fp)) { 
                printf("ERROR: open file at 4 %s\n", name);
                return -1;
            }
        } 
        else if (io == bio_use_open) { 
/**
            file->fd = open(name, O_RDONLY|O_CONCURRENT_WRITE, 0);
**/
            file->fd = open(name, O_RDONLY);

            if (file->fd == -1) { 
                printf("ERROR: open file at 5 %s\n", name);
                return -1;
            }
        }
#endif 
    }
    else if (mode == bio_file_read_write) {   
#ifdef MPI 
        if ((io == bio_collective) || (io == bio_independent)) {
            MPI_File_open(bio_comm, name, MPI_MODE_RDWR, bio_info, &(file->mfile));
        }
        else if (io == bio_use_fopen) { 
            file->fp = fopen(name, "r+b");
            if (!(file->fp)) { 
                printf("ERROR: open file at 6 %s\n", name);
                return -1; 
            }
        }
        else if (io == bio_use_open) { 
/***
            file->fd = open(name, O_RDWR|O_CONCURRENT_WRITE, 0);
***/
            file->fd = open(name, O_RDWR, 0644);

            if (file->fd == -1) { 
                printf("ERROR: open file at 7 %s \n", name);
                return -1;
            }
        } 
#else 
        if ((io == bio_use_fopen) || (io == bio_collective) || (io == bio_independent)) { 
            file->fp = fopen(name, "r+b"); 
            if (!(file->fp)) { 
                printf("ERROR: open file at 8 %s\n", name);
                return -1; 
            }
        }
        else if (io == bio_use_open) { 
            file->fd = open(name, O_RDWR|O_CONCURRENT_WRITE, 0);
            if (file->fd == -1) { 
                printf("ERROR: open file at 9 %s \n", name);
                return -1;
            }
        }
#endif 
    }
    slen = strlen(name) + 1;
    file->name = (char *) malloc((size_t)slen);
    strcpy(file->name, name);
    file->mode = mode;
    file->io   = io;
    obj = (io_Obj *) malloc(sizeof(io_Obj));
    file->object = obj;

    if (mode == bio_file_create) {  
      
        io_my_hdr_size = io_hdr_size;
        file->id  = io_next_id;
        file->id_min =  file->id;  /* future data sets may have smaller ids,
                                      if they are initialized before the file
                                      is created */ 
        io_next_id++;

        num_obj = 128;
        io_object_init(obj, file->id, num_obj, io_initial_size);
        parent_id = -1;
        offset = 0;
        size   = 0; 
        gsize  = 0;

/**     file->id doesn't have to reset in the following io_add_object **/
         
        io_add_object(1, NULL,
                      file, name, bio_group, bio_datatype_invalid, &(file->id), 
                      offset, size, gsize, NULL);
    } 
    else if (mode == bio_file_read_only || mode == bio_file_read_write)  {  

         io_read_hdrtail(0, file, io, obj, &successful);

         if (!successful) { 
             io_read_hdrtail(1, file, io, obj, &successful);
         }
         if (!successful) { 
             if (!mype) { 
                 printf("Warning: failed to open file %s\n", file->name);
             }
             return -1;
         } 
    }   
    return 0;
   } 
/****************************************************************************************/
int io_read_hdrtail(int ifbackup, io_File *file, bio_IO_Mode io, io_Obj *obj, int *successful)
{    
    char filename_backup[128];
    long long llversion_num;
    long       lversion_num;
    int        iversion_num;
    short      sversion_num;

    char achar;
    char *data_backup;
    int err, i, k, n, id, parent_id, slen, size, num_obj, num_attr, fid_stord, fd;
    int size_backup, version_num, my_id_min, my_id_max, id_adjust, m;
    int i2[2], n1, n2, n3, n4, n5, n6;
    char *pc, *c;
    short int a5s[5];
    int a5i[5];
    long a5l[5];
    long long a5ll[5];
    void *pv;
    int  *pi;
    long *pl;
    short int *ps;
    long long *pll;
    bio_Data_Type datatype_i, datatype_ll;
    void *meta_data;
    long long offset, myoffset, gsize;
    long long size_data, size_meta_block, size_meta_stord, size_attr;
    io_Attr *a;
    FILE *fp;
    struct stat buffer;

    *successful = 0;
    data_backup = NULL;
    meta_data   = NULL;

#ifdef MPI 
    MPI_Status status;
    file->t_open = MPI_Wtime();
#else
    time(&(file->t_open));
#endif

    if (!ifbackup) { 
  
        size = 0;
#ifdef MPI
        if ((io == bio_collective) || (io == bio_independent)) {
            myoffset = 0;
            size = MPI_File_read_at(file->mfile,(MPI_Offset)myoffset,io_hdr,(int)io_hdr_size,MPI_CHAR,&status);
            if (!size) size = io_hdr_size;
        }
        else if (io == bio_use_fopen) {
            fseek(file->fp, 0L, SEEK_SET);
            size = fread(io_hdr, sizeof(char), (size_t)io_hdr_size, file->fp);
        }
        else if (io == bio_use_open) {
            size = pread(file->fd, io_hdr, (size_t)io_hdr_size,(off_t)0);
        }
#else
        if ((io == bio_use_fopen) || (io == bio_collective) || (io == bio_independent)) {
            fseek(file->fp, 0L, SEEK_SET);
            size = fread(io_hdr, sizeof(char), (size_t)io_hdr_size, file->fp);
        }
        else if (io == bio_use_open) {
            lseek(file->fd, (off_t)0, 0);
            size = pread(file->fd, io_hdr, (size_t)io_hdr_size,(off_t)0);
        }
#endif
        if (size != io_hdr_size) {
            /* if (mype == 0) printf("not bio file\n"); */
            *successful = 0;
            return  0;
        }
    }
    else { 
        data_backup = NULL;
        if (mype == 0) { 
            if (io_filename_backup) {   
                strcpy(filename_backup, io_filename_backup);
            }
            else { 
                sprintf(filename_backup, "%s.backup", file->name);
            } 
            fp  = fopen(filename_backup, "rb");
            err = stat(filename_backup, &buffer);

            if (!err) {
                size_backup = buffer.st_size;
                data_backup = (char *) malloc(size_backup);
                size = fread(data_backup, sizeof(char), (size_t)size_backup, fp);
                assert(size == size_backup);
                assert(size >= io_hdr_size);
                memcpy(io_hdr, data_backup, (size_t)io_hdr_size); 
                achar = io_hdr[0];
                if (achar == 'a') {
                    version_num = (int)(io_hdr[1] - achar);
                    if (version_num == io_version_num) {
                        *successful = 1;
                    }
                }
            } 
            else { 
                *successful = 0;
                if (fp) fclose(fp);
                return -1;
            }
            fclose(fp);
            i2[0] = *successful;
            i2[1] = size_backup;
        }
        if (npes > 1) { 
#ifdef MPI
            MPI_Bcast(i2, 2, MPI_INT, 0, bio_comm);
            *successful = i2[0];
            size_backup = i2[1];
            if (*successful) {
                if (!data_backup) { 
                    data_backup = (char *) malloc(size_backup);
                }
                MPI_Bcast(data_backup, size_backup, MPI_CHAR, 0, bio_comm);
                if (!mype) {
                    memcpy(io_hdr, data_backup, (size_t)io_hdr_size);
                }
            } 
#endif
        }
        if (*successful) { 
            memcpy(io_hdr, data_backup, (size_t)io_hdr_size); 
            meta_data = data_backup + io_hdr_size;
        }
        else { 
            if (data_backup) free(data_backup);
            return 0;
        }
    }
    achar              = io_hdr[0];
    io_my_version_num  = (int)(io_hdr[1] - achar);
    io_szint_stord     = (int)(io_hdr[2] - achar);
    io_szlong_stord    = (int)(io_hdr[3] - achar);
    io_szllong_stord   = (int)(io_hdr[4] - achar);
    io_szfloat_stord   = (int)(io_hdr[5] - achar);
    io_szdouble_stord  = (int)(io_hdr[6] - achar);
    io_szldouble_stord = (int)(io_hdr[7] - achar);

    c = io_hdr + 8;
    if (io_szllong_stord == io_szllong) {
        memcpy(&llversion_num, c, (size_t)io_szllong_stord);
        if (llversion_num != io_my_version_num) {
            io_2swap = 1;
        }
    }
    else if (io_szllong_stord == io_szint) {
        memcpy(&iversion_num, c, (size_t)io_szllong_stord);
        if (iversion_num != (int) io_my_version_num) {
            io_2swap = 1;
        }
    }
    else if (io_szllong_stord == io_szlong) {
        memcpy(&lversion_num, c, (size_t)io_szllong_stord);
    }
    else if (io_szllong_stord == io_szshort) {
        memcpy(&sversion_num, c, (size_t)io_szllong_stord);
        if (sversion_num != (short) io_my_version_num) {
            io_2swap = 1;
        }
    }
    else {
        *successful = 0;
        if (data_backup) free(data_backup);
        return 0;
    }
    io_calc_datatype_stord2new();
    pc = io_hdr + io_nchars_in_hdr;
    if (io_2swap) {
        io_swap(pc, (int) io_szllong_stord, 5ll);
    }
    if (io_szllong_stord == io_szllong) {
        /***
        pll = (long long *)(io_hdr + io_nchars_in_hdr);
        io_my_hdr_size  = (int) pll[0];
        my_id_min       = (int) pll[1];
        my_id_max       = (int) pll[2];
        size_data       =       pll[3];
        size_meta_block =       pll[4];
        ***/

        memcpy(a5ll, pc, (size_t)(5 * io_szllong_stord));
        io_my_hdr_size  = (int) a5ll[0];
        my_id_min       = (int) a5ll[1];
        my_id_max       = (int) a5ll[2];
        size_data       =       a5ll[3];
        size_meta_block =       a5ll[4];
    }
    else if (io_szllong_stord == io_szint) {
        /*****
         pi = (int *)(io_hdr + io_nchars_in_hdr);
         io_my_hdr_size  = pi[0];
         my_id_min       = pi[1];
         my_id_max       = pi[2];
         size_data       = (long long) pi[3];
         size_meta_block = pi[4];
         ****/
         memcpy(a5i, pc, (size_t)(5 * io_szllong_stord));
         io_my_hdr_size  = a5i[0];
         my_id_min       = a5i[1];
         my_id_max       = a5i[2];
         size_data       = (long long) a5i[3];
         size_meta_block = a5i[4];
    }
    else if (io_szllong_stord == io_szlong) {
         /***
         pl = (long *)(io_hdr + io_nchars_in_hdr);
         io_my_hdr_size  = (int)pl[0];
         my_id_min       = (int)pl[1];
         my_id_max       = (int)pl[2];
         size_data       = (long long) pl[3];
         size_meta_block =  pl[4];
         ****/
         memcpy(a5l, pc, (size_t)(5 * io_szllong_stord));
         io_my_hdr_size  = (int) a5l[0];
         my_id_min       = (int) a5l[1];
         my_id_max       = (int) a5l[2];
         size_data       = (long long) a5l[3];
         size_meta_block = a5l[4];
   }
    else if (io_szllong_stord == io_szshort) {
         /*****
         ps = (short int *)(io_hdr + io_nchars_in_hdr);
         io_my_hdr_size  = (int)ps[0];
         my_id_min       = (int)ps[1];
         my_id_max       = (int)ps[2];
         size_data       = (long long) ps[3];
         size_meta_block = ps[4];
         ****/
         memcpy(a5s, pc, (size_t)(5 * io_szllong_stord));
         io_my_hdr_size  = (int) a5s[0];
         my_id_min       = (int) a5s[1];
         my_id_max       = (int) a5s[2];
         size_data       = (long long) a5s[3];
         size_meta_block = a5s[4];
     }
     else {  
         *successful = 0;
         if (data_backup) free(data_backup); 
         return 0;
     } 
/*   all the ids in this file will be subtracted by id_adjust */

     id_adjust    = io_next_id - my_id_min;
     file->id_min = io_next_id;

     io_next_id += (my_id_max - my_id_min + 1);

     /* read all the meta data */

     if (!ifbackup) { 
         meta_data = malloc((size_t)size_meta_block);
         myoffset  = (long long)io_my_hdr_size + size_data;
         size      = 0;
#ifdef MPI 
         if ((io == bio_collective) || (io == bio_independent)) {
             size = MPI_File_read_at(file->mfile, (MPI_Offset)myoffset, meta_data, (int)size_meta_block, MPI_CHAR, &status);
             if (!size) size = size_meta_block;
         }
         else if (io == bio_use_fopen) {
             fseek(file->fp, (long)myoffset, SEEK_SET);
             size = (long long) fread(meta_data, sizeof(char), (size_t) size_meta_block, file->fp);
         }
         else if (file->io == bio_use_open) {
             size = (long long) pread(file->fd, meta_data, (size_t)size_meta_block, (off_t)myoffset);
         }
#else
         if ((io == bio_use_fopen) ||(io == bio_collective) || (io == bio_independent) ) {
             fseek(file->fp, (long)myoffset, SEEK_SET);
             size = (long long) fread(meta_data, sizeof(char), (size_t) size_meta_block, file->fp);
         }
         else if (file->io == bio_use_open) {
             size = (long long) pread(file->fd, meta_data, (size_t)size_meta_block,(off_t)myoffset);
         }
#endif
         if (size != size_meta_block) {
             free(meta_data);
             *successful = 0;
             return 0;
         }
     }
     if (io_2swap) {
         io_swap_metadata(meta_data);
     }
     datatype_ll = io_datatype_stord2new[(int)bio_long_long - (int)bio_char];

     *successful = 1;

     if (datatype_ll == bio_long_long) {
        /***
         pll = (long long *) meta_data;
         size_meta_stord = pll[0];
         fid_stord = pll[1];
         num_obj   = pll[2];
         size_data = pll[3];
         pv = (void *) (pll + 4);
        ****/
  
         memcpy(a5ll, meta_data, (size_t)(4 * io_szllong));
         size_meta_stord = a5ll[0];
         fid_stord = (int) a5ll[1];
         num_obj   = (int) a5ll[2];
         size_data = (long long) a5ll[3];
     }
     else if (datatype_ll == bio_long) {
         /***
         pl = (long *)  meta_data;
         size_meta_stord = pl[0];
         fid_stord = pl[1];
         num_obj   = pl[2];
         size_data = pl[3];
         pv = (void *) (pl + 4);
         ****/
  
         memcpy(a5l, meta_data, (size_t)(4 * io_szlong));
         size_meta_stord = a5l[0];
         fid_stord = (int) a5l[1];
         num_obj   = (int) a5l[2];
         size_data = (long long) a5l[3];
     }
     else if (datatype_ll == bio_int) {
         /****
         pi = (int *)  meta_data;
         size_meta_stord = pi[0];
         fid_stord = pi[1];
         num_obj   = pi[2];
         size_data = pi[3];
         pv = (void *) (pi + 4);
         ****/
  
         memcpy(a5i, meta_data, (size_t)(4 * io_szint));
         size_meta_stord = a5i[0];
         fid_stord = (int) a5i[1];
         num_obj   = (int) a5i[2];
         size_data = (long long) a5i[3];
     }
     else if (io_szint_stord == io_szshort) {
         /***
         ps = (short int *) meta_data;
         size_meta_stord = ps[0];
         fid_stord = ps[1];
         num_obj   = ps[2];
         size_data = ps[3];
         pv = (void *) (ps + 4);
         ****/
  
         memcpy(a5s, meta_data, (size_t)(4 * io_szshort));
         size_meta_stord = a5s[0];
         fid_stord = (int) a5s[1];
         num_obj   = (int) a5s[2];
         size_data = (long long) a5s[3];
     }
     else {
         *successful = 0;
         if (!ifbackup) { 
             if (meta_data) free(meta_data);
         }
         if (data_backup) free(data_backup);
         return 0;
     }
     file->id = fid_stord + id_adjust;
     io_convert_metadata(fid_stord, obj, num_obj, id_adjust, meta_data, size_meta_stord);

     io_szint_stord   = io_szint;
     io_szlong_stord  = io_szlong;
     io_szllong_stord = io_szllong;
     io_szfloat_stord = io_szfloat;
     io_szdouble_stord  = io_szdouble;
     io_szldouble_stord = io_szldouble;

     file->root = NULL;
     file->buffer = NULL;
     io_construct_tree(file, obj);

     if (!ifbackup) {
         if (meta_data) free(meta_data);
     }
     if (data_backup) free(data_backup);

     return 0;
  } 
/****************************************************************************************/
int io_convert_metadata(int fid_stord, io_Obj *obj, int num_obj, int id_adjust,
                        void *meta_data, long long size_meta_stord)
{

    int i, k, idx, slen, num_attr, nc_skip, same_datatypes, id_stord;
    char *t, *s;
    io_Attr *a;

    long long max_size_meta, size_meta;
    long long size_attr_stord,  max_size_attr, size_attr, size_value;
    long long a6[6];
    int a6i[6];
    int a6l[6];

    bio_Data_Type datatype;
    
    if (io_szllong_stord == io_szllong &&
        io_szlong_stord  == io_szlong  && 
        io_szint_stord   == io_szint   &&
        io_szfloat_stord == io_szfloat &&
        io_szdouble_stord == io_szdouble &&
        io_szldouble_stord == io_szldouble) {
      
        same_datatypes = 1;
    }
    else {
        same_datatypes = 0;
    }
/*
    size_meta,          szllong
    fid,                szllong
    num,                szllong 
    size_data,          szllong
    size of name,       szllong
    name,               slen*szchar
    id,                 szllong 
    parent_id,          szllong  
    obj_type,           szllong
    datatype,           szllong
    offset,             szllong
    size,               szllong
*/
    nc_skip = 7 * io_szllong;

    size_meta = 4 * io_szllong;
    s = (char *) meta_data;

    if (io_szllong_stord < io_szllong) {
        max_size_meta = 1 + size_meta_stord / io_szllong_stord;
        max_size_meta *= io_szllong;
    }
    else {  
        max_size_meta = size_meta_stord;
    }   
    if (io_szllong_stord == io_szllong && !id_adjust) {

        io_object_init(obj, fid_stord, num_obj, size_meta_stord);
        memcpy(obj->value, meta_data, (size_t)size_meta_stord);
    }
    else if (io_szllong_stord == io_szllong) {
        
        io_object_init(obj, fid_stord, num_obj, max_size_meta); /* fid_stord will be corrected */ 
        t = (char *) obj->value;

        memcpy(a6, s, (size_t)(4 * io_szllong_stord));
        s += (4 * io_szllong_stord);
        a6[1] += id_adjust;
        memcpy(t, a6, (size_t)(4 * io_szllong));
        t += (4 * io_szllong);

        for (i = 0; i < num_obj; i++) {
            memcpy(a6, s, (size_t) io_szllong_stord);
            s += io_szllong_stord;
            slen = a6[0];
            memcpy(t, a6, (size_t) io_szllong);
            t += io_szllong;
            memcpy(t, s, (size_t)slen);
            t += slen;
            s += slen;
            memcpy(a6, s, (size_t)(6 * io_szllong_stord));
            s += (6 * io_szllong_stord);
            a6[0] += id_adjust;
            a6[1] += id_adjust;
            idx = a6[3] - (int) bio_char;
            datatype = io_datatype_stord2new[idx];
            a6[3] = (long long) datatype;
            memcpy(t, a6, (size_t)(nc_skip - io_szllong));
            t += (nc_skip - io_szllong);
            size_meta += (slen + nc_skip);
        }
        memcpy(obj->value, &size_meta, (size_t) io_szllong);
    } 
    else if (io_szllong_stord == io_szint) {

        io_object_init(obj, fid_stord, num_obj, max_size_meta); /* fid_stord will be corrected */
        t = (char *) obj->value;

        memcpy(a6i, s, (size_t)(4 * io_szllong_stord));
        s += (4 * io_szllong_stord);
        a6[0] = a6i[0];
        a6[1] = a6i[1] + id_adjust;
        a6[2] = a6i[2];
        a6[3] = a6i[3];
        memcpy(t, a6, (size_t)(4 * io_szllong));
        t += (4 * io_szllong);
       
        for (i = 0; i < num_obj; i++) { 
            memcpy(a6i, s, (size_t) io_szllong_stord);
            s += io_szllong_stord;
            slen = a6i[0];
            a6[0] = a6i[0];
            memcpy(t, a6, (size_t) io_szllong);
            t += io_szllong;
            memcpy(t, s, (size_t)slen);
            t += slen;
            s += slen;
            memcpy(a6i, s, (size_t)(6 * io_szllong_stord));
            s += (6 * io_szllong_stord);
            a6[0] = a6i[0] + id_adjust;
            a6[1] = a6i[1] + id_adjust;
            a6[2] = a6i[2];
            idx = a6i[3] - (int) bio_char;
            datatype = io_datatype_stord2new[idx]; 
            a6[3] = (long long) datatype;  
            a6[4] = a6i[4];
            a6[5] = a6i[5];
            memcpy(t, a6, (size_t)(nc_skip - io_szllong));
            t += (nc_skip - io_szllong);  
            size_meta += (slen + nc_skip); 
        }
        memcpy(obj->value, &size_meta, (size_t) io_szllong);  
    }  
    else if (io_szllong_stord == io_szlong) {


    }
    else { 
       if (!mype) printf("ERROR: io_szllong_stord not covered in io_convert_metadata\n");
       return -1;
    } 

/****
    for attrs

    size,            szllong
    id,              szllong
    num,             szllong
    size of name,    szllong
    name,            slen*szchar
    datatype,        szllong
    size of value,   szllong
    values,          n*io_sizeof(datatype)
***/
    s = (char *) meta_data + size_meta_stord;

    if (same_datatypes && !id_adjust) {

        for (i = 0; i < num_obj; i++) {
            memcpy(a6, s, (size_t)(io_szllong_stord + io_szllong_stord));
            size_attr_stord = a6[0];
            id_stord  = a6[1];
            a = (io_Attr *) malloc(sizeof(io_Attr));
            io_attr_init(id_stord, a, size_attr_stord);
            memcpy(a->value, s, (size_t)size_attr_stord);
            obj->attr[i] = a;
            s += size_attr_stord;
        }
    }
    else if (io_szllong_stord == io_szllong) { 

        for (i = 0; i < num_obj; i++) {
            memcpy(a6, s, (size_t)(3 * io_szllong_stord));
            s += (3 * io_szllong_stord);
            size_attr_stord = a6[0];
            id_stord = a6[1] + id_adjust;
            a6[1] = id_stord;

            num_attr = a6[2];
            a = (io_Attr *) malloc(sizeof(io_Attr));
            io_attr_init(id_stord, a, size_attr_stord);
            obj->attr[i] = a;

            t = (char *) a->value;
            memcpy(t, a6, (size_t)(3 * io_szllong));
            t += (3 * io_szllong);

            for (k = 0; k < num_attr; k++) {
                memcpy(a6, s, (size_t) io_szllong_stord);
                s += io_szllong_stord;
                slen = a6[0];
                memcpy(t, a6, (size_t)io_szllong);
                t += io_szllong;
                memcpy(t, s, (size_t) slen);
                t += slen;
                s += slen;
                memcpy(a6, s, (size_t) (io_szllong_stord + io_szllong_stord));
                s += (io_szllong_stord + io_szllong_stord);
                idx = a6[0] - (int) bio_char;
                datatype = io_datatype_stord2new[idx];
                a6[0] = (long long) datatype;
                size_value = a6[1];
                memcpy(t, a6, (size_t)(io_szllong + io_szllong));
                t += (io_szllong + io_szllong);
                memcpy(t, s, (size_t) size_value);
                s += size_value;
                t += size_value;
            }
        }
    }
    else if (io_szllong_stord == io_szint) {

        for (i = 0; i < num_obj; i++) {
            memcpy(a6i, s, (size_t)(3 * io_szllong_stord));
            s += (3 * io_szllong_stord);
            size_attr_stord = a6i[0];
            id_stord = a6i[1] + id_adjust;

            if (same_datatypes) { 
                max_size_attr = size_attr_stord;
            }
            else { 
                max_size_attr = 4 * size_attr_stord;
            }
            num_attr = a6i[2];
            a = (io_Attr *) malloc(sizeof(io_Attr));
            io_attr_init(id_stord, a, max_size_attr);
            obj->attr[i] = a;
            
            t = (char *) a->value;
            a6[0] = size_attr_stord;  /* this will be corrected later */  
            a6[1] = (long long) (a6i[1] + id_adjust);
            a6[2] = (long long) num_attr;
            memcpy(t, a6, (size_t)(3 * io_szllong));
            t += (3 * io_szllong);
 
            size_attr = 3 * io_szllong;

            for (k = 0; k < num_attr; k++) { 
                memcpy(a6i, s, (size_t) io_szllong_stord);
                s += io_szllong_stord;
                slen = a6i[0];
                a6[0] = slen;
                memcpy(t, a6, (size_t)io_szllong);
                t += io_szllong;
                memcpy(t, s, (size_t) slen);
                t += slen;
                s += slen;
                memcpy(a6i, s, (size_t) (io_szllong_stord + io_szllong_stord));
                s += (io_szllong_stord + io_szllong_stord);
                idx = a6i[0] - (int) bio_char;
                datatype = io_datatype_stord2new[idx];
                a6[0] = (long long) datatype;  
                a6[1] = a6i[1];
                size_value = a6[1];
                memcpy(t, a6, (size_t)(io_szllong + io_szllong));
                t += (io_szllong + io_szllong); 
                memcpy(t, s, (size_t) size_value);
                s += size_value;
                t += size_value;

                size_attr += (3 * io_szllong + slen + size_value); 
            } 
            memcpy(a->value, &size_attr, (size_t) io_szllong);
        }
    }
  
     return 0;
 }
/****************************************************************************************/
int io_construct_tree(io_File *file, io_Obj *obj)
{
    int i, id, num, sztype, slen, parentid;
    char *c, *name, *ptr;

    long long size_meta, size_data, size, gsize;
    long long a6[6];
    size_t sztree, szchild;

    bio_Object_Type type; 
    bio_Data_Type datatype;
    io_Tree *tree, *ptree;
    io_Child *child;
    
/*
    size_meta,          szllong
    fid,                szllong
    num,                szllong 
    size_data,          szllong
    size of name,       szllong
    name,               slen*szchar
    id,                 szllong 
    parent_id,          szllong  
    obj_type,           szllong
    datatype,           szllong
    offset,             szllong
    size,               szllong
*/
    
     if (file->root) {
         io_clean_tree(file->root);
     } 
     file->root = NULL;
     file->buffer = NULL;

     sztree = sizeof(io_Tree);
     szchild = sizeof(io_Child);

     c = (char *) obj->value;
     memcpy(a6, c, (size_t)(4 * io_szllong));
     size_meta = a6[0];
     id  = a6[1];
     num = a6[2];
     size_data = a6[3];

     c += (4 * io_szllong);
     for (i = 0; i < num; i++) { 
         ptr = c;
         memcpy(a6, c, (size_t)io_szllong);
         slen = a6[0];
         name = c + io_szllong;
         c = name + slen;
         memcpy(a6, c, (size_t)(6 * io_szllong));
         id = a6[0];
         parentid = a6[1];
         type = (bio_Object_Type) ((int)a6[2]);
         datatype = (bio_Data_Type) ((int)a6[3]);
         sztype = io_sizeof(datatype);
         size = a6[5];
         gsize = size/sztype;
         
         tree = (io_Tree *) malloc(sztree);
         tree->id = id;
         tree->idx = i;
         tree->p = ptr;

         ptree = io_find_tree(file, parentid);
         tree->parent = ptree;
         tree->child = NULL;

         if (ptree) { 
             child = (io_Child *) malloc(szchild);
             child->tree = tree;
             child->next = ptree->child;
             ptree->child = child;
         }
         else { 
             file->root = tree;
         } 
     
         c += (6 * io_szllong);
     }     
  
     return 0;
 }
/****************************************************************************************/
int io_calc_datatype_stord2new(void)
{  
    int i, n;

    n = (int) bio_datatype_invalid - (int) bio_char;
    for (i = 0; i < n; i++) {
        io_datatype_stord2new[i] = bio_datatype_invalid;
    }

    i = (int)bio_char;
    io_datatype_stord2new[0] = bio_char;     

    n = (int)bio_double - i;
    if (io_szdouble_stord == io_szdouble)
       io_datatype_stord2new[n] = bio_double;
    else if (io_szdouble_stord == io_szfloat)
       io_datatype_stord2new[n] = bio_float;
    else if (io_szdouble_stord == io_szldouble)
       io_datatype_stord2new[n] = bio_long_double;

    n = (int)bio_float - i;
    if (io_szfloat_stord == io_szfloat)
       io_datatype_stord2new[n] = bio_float;
    else if (io_szfloat_stord == io_szdouble)
       io_datatype_stord2new[n] = bio_double;
    else if (io_szfloat_stord == io_szldouble)
       io_datatype_stord2new[n] = bio_long_double;

    n = (int)bio_int - i;
    if (io_szint_stord == io_szint) 
       io_datatype_stord2new[n] = bio_int;
    else if (io_szint_stord == io_szlong) 
       io_datatype_stord2new[n] = bio_long; 
    else if (io_szint_stord == io_szllong) 
       io_datatype_stord2new[n] = bio_long_long;

    n = (int)bio_long - i;
    if (io_szlong_stord == io_szlong) 
       io_datatype_stord2new[n] = bio_long;
    else if (io_szlong_stord == io_szint) 
       io_datatype_stord2new[n] = bio_int;                            
    else if (io_szlong_stord == io_szllong) 
       io_datatype_stord2new[n] = bio_long_long;

    n = (int)bio_long_double - i;
    if (io_szldouble_stord == io_szldouble)
       io_datatype_stord2new[n] = bio_long_double;
    else if (io_szldouble_stord == io_szfloat)
       io_datatype_stord2new[n] = bio_float;
    else if (io_szldouble_stord == io_szdouble)
       io_datatype_stord2new[n] = bio_double;

    n = (int)bio_long_long - i;
    if (io_szllong_stord == io_szllong)
       io_datatype_stord2new[n] = bio_long_long;
    else if (io_szllong_stord == io_szint) 
       io_datatype_stord2new[n] = bio_int; 
    else if (io_szllong_stord == io_szlong)
       io_datatype_stord2new[n] = bio_long;

    return 0;
   } 
/****************************************************************************************/
int bio_correct_file(char *filename)
{
    char filename_backup[128];
    char achar, hdr[io_hdr_size], *data_backup;
    int  err, size_to_read, size_read, size_backup, answer, version_num;
    FILE *fp;
    struct stat buffer;

    answer = 0;

    if (mype == 0) {
        assert(filename);
        fp = fopen(filename, "rb");
        assert(fp);
        size_to_read = 2;
        size_read = fread(hdr, sizeof(char), (size_t) size_to_read, fp);
        if (size_read != size_to_read) { 
            fclose(fp);
            return answer;
        }
        achar = hdr[0];
        if (achar == 'a') {
            version_num = (int)(hdr[1] - achar);
            if (version_num == io_version_num) {
                answer = 1;
            }
        }
        fclose(fp);
        if (!answer) {
            if (io_filename_backup) {
                strcpy(filename_backup, io_filename_backup);
            }
            else {
                sprintf(filename_backup, "%s.backup", filename);
            }
            err = 1;
            fp  = fopen(filename_backup, "rb");
            if (fp) {
                err = stat(filename, &buffer);
            }
            if (!err) {
                size_backup = buffer.st_size;
                data_backup = (char *) malloc(size_backup);
                size_read   = fread(data_backup, sizeof(char), (size_t)size_backup, fp);
                assert(size_read == size_backup);
                if (achar == 'a') {
                    version_num = (int)(hdr[1] - achar);
                    if (version_num == io_version_num) {
                        answer = 1;
                    }
                }
                free(data_backup);
            }
            if (fp) fclose(fp);
        }
    }
    if (npes > 1) {
#ifdef MPI
        MPI_Bcast(&answer, 1, MPI_INT, 0, bio_comm);
#endif
    }
    return answer;
 }
/****************************************************************************************/
int bio_resilient(int ifresilient)
{
    if (ifresilient >= 0) { 
        flag_for_resilience = ifresilient;
    }
    else { 
        if (!mype) { 
            printf("Warning: bio_resilient ignored with argument = %d\n",
                    ifresilient); 
        }
    } 
    return 0;
}
/****************************************************************************************/
int bio_backup_file(char *backup_file)
{
    int slen;

    if (io_filename_backup) free(io_filename_backup);
    if (!backup_file) { 
        io_filename_backup = NULL;
    }
    else { 
        slen = strlen(backup_file) + 1;
        io_filename_backup = (char *) malloc(slen);
        strcpy(io_filename_backup, backup_file);
    }
    return 0;
}  
/****************************************************************************************/
int bio_buffer_mode(int flg_for_buffered_read)
{
    flag_for_buffered_read = flg_for_buffered_read;
    return 0;
} 
/****************************************************************************************/
int bio_file_open(char *name, bio_File_Mode mode, bio_IO_Mode io, int *fileid)
{   

    int err;
    io_File *file;

    if (!fileid) {  
        printf("ERROR: null fileid in bio_file_open\n");
        return -1;
    }
    err = io_file_open(name, mode, io, &file);
    if (!err) { 
        *fileid = file->id;
    }
    else { 
        *fileid = -1;
    }
    return err;
  }
/****************************************************************************************/
int io_file_open(char *name, bio_File_Mode mode, bio_IO_Mode io, io_File **myfile)
{   
#ifdef MPI
    MPI_Comm world_comm;
    MPI_Group world_comm_group,comm_group;
#endif 
    int err, k, n;
    int *ranks;
    io_File *file;

    if (!name) { 
       printf("ERROR: null name in bio_file_open\n");
       return -1;
    }
#ifdef MPI  
    MPI_Comm_rank(bio_comm, &mype);
    MPI_Comm_size(bio_comm, &npes);
#endif
    if (!io_hdr) {  
         io_hdr = (char *) malloc((size_t)io_hdr_size);
         io_init();
    }
    if (io_num_file == 0 && mode == bio_file_create)
       { io_szint_stord     = io_szint;
         io_szlong_stord    = io_szlong;
         io_szllong_stord   = io_szllong;
         io_szfloat_stord   = io_szfloat;
         io_szdouble_stord  = io_szdouble;
         io_szldouble_stord = io_szldouble;
        } 
    file = (io_File *) malloc(sizeof(io_File));
    io_file_null(file);
    *myfile = file;
    err  = io_add_file(name, mode, io, file);  
    
    if (err != 0)
       { free(file);
         return err; 
        }
    if (!io_files)
       { io_files = file;
         io_files->next = NULL;
        }
    else
       { file->next = io_files;
         io_files   = file;
        }
    io_num_file++; 

#ifdef MPI
    if (mpi_status == NULL) {  
        mpi_status = (MPI_Status *) malloc((size_t)npes * sizeof(MPI_Status));
    }
    if (mpi_reqs == NULL) {
        mpi_reqs = (MPI_Request *) malloc((size_t)npes * sizeof(MPI_Request));
    } 
#endif

    return err;
  }
/****************************************************************************************/
int io_file_null(io_File *file) 
{
#ifdef MPI 
       file->mfile = (MPI_File)0;     
       file->t_open = 0.0;
       file->t_close = 0.0;
#endif
       file->fp = NULL;
       file->fd = 0;
       file->mode = bio_file_create;
       file->io  = bio_use_fopen;
       file->id  = -1;
       file->id_min = -1;
       file->object = NULL;
       file->name = NULL;
       file->mbytes = 0.0;
       
       file->root = NULL;
       file->buffer = NULL;
       file->next = NULL;
  
       return 0;
 } 
       

/****************************************************************************************/
int bio_group_open(int fileid, int parent_id, char *name, int *grp_id)
{ 
    int obj_idx, gparent_id, id_stord, i, slen, myslen, err;
    bio_Object_Type type_stord;
    bio_Data_Type datatype;
    long long offset, gsize, a6[6];
    char *c, *myname;
    io_File *file;
    io_Tree *tree;
    
    err = 0;

    if (!name) {   
        printf("ERROR: null name in bio_group_open\n");
        return -1;
    }
    if (!grp_id) {  
        printf("ERROR: null grp_id in bio_group_open\n");
        return -1;
    }
    slen = strlen(name);
    for (i = 0; i < slen; i++) {   
        if (name[i] == '/') { 
            printf("ERROR: '/' in name not allowd in bio_group_open\n");
            return -1;
        }
    }
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != parent_id) {
        tree = io_find_tree(file, parent_id);
    }
    if (tree) {
        c = (char *)tree->p;
        memcpy(a6, c, (size_t)io_szllong);
        myslen = a6[0];
        myname = c + io_szllong;
        c = myname + myslen;
        memcpy(a6, c, (size_t)(6 * io_szllong));
        id_stord = (int) a6[0];

        if (id_stord != parent_id) { 
            if (!mype) printf("ERROR: parent_is != id_stord\n");
            return -1;
        } 
        gparent_id = (int)a6[1];
        type_stord = (bio_Object_Type) ((int)a6[2]);
        datatype = (bio_Data_Type) ((int)a6[3]);
        offset = a6[4];
        gsize  = a6[5];
    }
    else { 
        printf("ERROR: parent_id not found in bio_group_open at 1, name = %s\n", name);
        return -1;
    }   
    if (type_stord != bio_group) { 
        printf("ERROR: parent_id is not for a group in bio_group_open\n");
        return -1;
    }
    *grp_id = -1;

    err = io_add_grp(&tree, file, name, grp_id, bio_group, bio_datatype_invalid);

    return err;
 } 
/****************************************************************************************/
int bio_group_open_i(int fileid, int parent_id, char *name, int *grp_id, 
                     bio_Object_Type objtype, bio_Data_Type datatype) 
{ 
/** This function is the same as bio_group_open, except that this function
    keep any non-negative grp_id 
 **/
    int obj_idx, gparent_id, id_stord, i, slen, myslen;
    bio_Object_Type type_stord;
    bio_Data_Type pdatatype;
    long long a6[6], offset, gsize;
    char *c, *myname;
    io_File *file;
    io_Tree *tree;

    if (!name)
       { printf("ERROR: null name in bio_group_open\n");
         return -1;
        }
    if (!grp_id)
       { printf("ERROR: null grp_id in bio_group_open\n");
         return -1;
        }
    slen = strlen(name);
    for (i = 0; i < slen; i++)
      { if (name[i] == '/')
           { printf("ERROR: '/' in name not allowd in bio_group_open\n");
             return -1;
            }
       }
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != parent_id) {
        tree = io_find_tree(file, parent_id);
    }
    if (tree) {
        c = (char *)tree->p;
        memcpy(a6, c, (size_t)io_szllong);
        myslen = a6[0];
        myname = c + io_szllong;
        c = myname + myslen;
        memcpy(a6, c, (size_t)(6 * io_szllong));
        id_stord = (int) a6[0];

        if (id_stord != parent_id) {
            if (!mype) printf("ERROR: parent_is != id_stord\n");
            return -1;
        }
        gparent_id = (int)a6[1];
        type_stord = (bio_Object_Type) ((int)a6[2]);
        pdatatype = (bio_Data_Type) ((int)a6[3]);
        offset = a6[4];
        gsize  = a6[5];
    }
    else {
        printf("ERROR: parent_id not found in bio_group_open at 2, name = %s\n", name);
        return -1;
    }
    if (!file)
       { printf("ERROR: parent_id not found in bio_group_open at 3, name = %s\n", name);
         return -1;
        }
    if (type_stord != bio_group)
       { printf("ERROR: parent_id is not for a group in bio_group_open\n");
         return -1;
        }

    io_add_grp(&tree, file, name, grp_id, objtype, datatype);

    return 0;
 }
/****************************************************************************************/
int io_add_grp(io_Tree **parent, io_File *file, char *name, int *grp_id,
               bio_Object_Type type, bio_Data_Type datatype)
{
    int  i, mytype, slen, myslen, num, s, slen_stord;
    int  sztype, grpid_stord, parentid;
    long long a6[6];
    long long size_meta_added, gsize_stord, mysize;
    char *pc, *c, *myname, *ptr;
    long long *pll, *po, *pd, *pid, *poffset;
    bio_Data_Type datatype_stord;
    bio_Object_Type   type_stord;
    io_Attr *a; 
    io_Obj *obj; 
    io_Tree *tree, *ptree; 
    io_Child *child, *next;

    obj = file->object; 
    num = *(obj->num); 
    slen = strlen(name);

    /* check whether this group already existed */

    if (parent) { 
        ptree = *parent;
        parentid = ptree->id;
    }
    else { 
        ptree = NULL;
        parentid = -1;
    }
    child = ptree->child;
    while (child) {
       tree = child->tree;
       c = (char *) tree->p;
       memcpy(a6, c, (size_t)io_szllong);
       myslen = a6[0];
       myname = c + io_szllong;
       c = myname + myslen;
       if (slen == myslen && !strncmp(name, myname, (size_t)myslen)) {
           memcpy(a6, c, (size_t)(6 * io_szllong));
           grpid_stord = (int) a6[0];
           type_stord = (bio_Object_Type) ((int)a6[2]);
           if (type_stord != bio_group) {
               printf("ERROR: an object with this name already exits in io_add_group\n");
               return -1;
           }
           *grp_id = grpid_stord;
           break;
       }
       next = child->next;
       child = next;
    }
    if (child) return 0;

    if (file->mode == bio_file_read_only) return -1;

    s  = 7 * io_szllong;

        if (*grp_id < 0) { 
            *grp_id = io_next_id;
            io_next_id++;
        }
        size_meta_added = slen + s;
        mysize = size_meta_added + *(obj->size_meta);
        if (mysize > obj->max_size) {
            io_object_value_allocate(obj, mysize);
            io_construct_tree(file, obj);
            if (parentid >= 0) { 
                ptree = io_find_tree(file, parentid);
                *parent = ptree;
            } 
        }
        if ((long long) num == obj->max_num) {
            io_object_num_allocate(obj);
        }

/*      In the order of
        size of name,       szllong
        name,               slen*szchar
        id,                 szllong
        parent_id,          szllong
        obj_type,           szllong
        datatype,           szllong
        offset,             szllong
        size,               szllong
*/
/*******
        pll  = (long long *)((char *)obj->value + *(obj->size_meta));
        *pll = (long long) slen;
        pc   = (char *)(pll + 1);
        memcpy(pc, name, (size_t)slen);
        pll = (long long *)(pc + slen);
        pll[0] = *grp_id;
        pll[1] = parent->id;
        po  = pll + 2;
        *po = (long long) type;
        pd  = po + 1;
        *pd = (long long) datatype;
        pll = pd + 1;
        pll[0] = 0;
        pll[1] = 0;
*****/

        pc = (char *)obj->value + *(obj->size_meta);
        ptr = pc;

        mysize = slen;
        memcpy(pc, &mysize, (size_t)io_szllong);
        pc += io_szllong;
        memcpy(pc, name, (size_t)slen);
        pc += slen;
        a6[0] = *grp_id;
        a6[1] = ptree->id;
        a6[2] = (long long) type;
        a6[3] = (long long) datatype;
        a6[4] = 0;
        a6[5] = 0; 
        memcpy(pc, a6, (size_t)(6 * io_szllong));

        a = (io_Attr *) malloc(sizeof(io_Attr));

        io_attr_init(*grp_id, a, io_initial_size);
        obj->attr[num] = a;
   
        *(obj->size_meta) += size_meta_added;
        *(obj->num) += 1;
    /*  *(obj->size_data) += a6[5];  */  

        /* create an object */
         
        tree = (io_Tree *) malloc(sizeof(io_Tree));
        tree->id = *grp_id;
        tree->idx = num; 
        tree->p = ptr;
        tree->parent = ptree;
        tree->child = NULL;
        if (ptree) { 
            child = (io_Child *) malloc(sizeof(io_Child));
            child->tree = tree;
            child->next = ptree->child;
            ptree->child = child;
        }
        else {
           file->root = tree;
        }
    if (file->id_min > *grp_id) { 
        file->id_min = *grp_id;
    }

    return 0;
 }
/****************************************************************************************/
int io_add_attr(io_Attr *a, char *name,  
                bio_Data_Type datatype, int size, void *buffer)
{   
    int  m, attr_num, i, slen, slen_stord;
    long long size_to_store, size_stord, size_add, mysize;
    long long *pll, *pd;
    long long a2[2];
    char *pc, *myname, *c;
    bio_Data_Type datatype_stord;
/*
   size,            szllong
   id,              szllong
   num,             szllong
   size of name,    szllong
   name,            slen*szchar
   datatype,        szllong
   size of value,   szllong
   values,          n*io_sizeof(datatype)
*/
    /* to check whether or not this attr already exist associated with this object */
   
    slen = strlen(name);
    size_to_store = size * io_sizeof(datatype);

    m = 3 * io_szllong; /* size_of_name, size_of_value, szdatatype */ 
    attr_num = *(a->num);
    pc = (char *) a->value + (3 * io_szllong);

    for (i = 0; i < attr_num; i++) { 
        memcpy(a2, pc, (size_t)io_szllong);
        slen_stord = (int) a2[0];
        myname = pc + io_szllong;
        c = myname + slen_stord;
        memcpy(a2, c, (size_t)(io_szllong + io_szllong));  
        datatype_stord = (bio_Data_Type) ((int)a2[0]);
        size_stord = a2[1];
        
        pc += (slen_stord + (int)size_stord + m);
        if (slen_stord != slen) continue;
        if (!strncmp(myname, name, (size_t)slen_stord)) { 
            if (datatype_stord != datatype) {
                if (!mype) { 
                    printf("ERROR: attr with name %s but different data type already exist for the object, ", name);
                    printf(" please use different attribute name\n");
                }
                return -1;
            } 
            if (size_stord != size_to_store) { 
                if (!mype) { 
                    printf("ERROR: attr with name %s but different length of values already exist for the object, ", name);
                    printf(" please use different attribute name\n");
                }
                return -1;
            } 
            c += (io_szllong + io_szllong);
            memcpy(c, buffer, (size_t)size_stord);
            return 0;
        }
    }
    /* a new attribute */

    if (datatype == bio_datatype_invalid) {
       printf("ERROR: datatype incorrect in io_add_attr\n");
       return -1;
    }
    slen = strlen(name);
    size_add = slen + (3 * io_szllong + size * io_sizeof(datatype));
    mysize = size_add + *(a->size);
    if (mysize > a->max_size)
       { io_attr_allocate(a, mysize);
        }
/****
    pll  = (long long *)((char *)a->value + *(a->size));
    *pll = slen;
    pc   = (char *)(pll + 1);
    memcpy(pc, name, (size_t)slen);
    pd  = (long long *)(pc + slen);
    *pd = (long long) datatype;
    pll = pd + 1;
    pll[0] = size * io_sizeof(datatype);
    memcpy(pll+1, buffer, (size_t)(pll[0]));
****/

    pc = (char *)a->value + *(a->size);
    mysize = slen;
    memcpy(pc, &mysize, (size_t)io_szllong);            /* size of the name */ 
    pc += io_szllong;
    memcpy(pc, name, (size_t)slen);                     /* name */ 
    pc += slen;
    a2[0] = (long long) datatype;
    a2[1] = size * io_sizeof(datatype);
    memcpy(pc, a2, (size_t)(io_szllong + io_szllong));  /* name and size of value */ 
    pc += (io_szllong + io_szllong);
    memcpy(pc, buffer, (size_t)a2[1]);                  /* values */  

    *(a->size) += size_add; 
    *(a->num)  += 1;  

    return 0;
  }
/****************************************************************************************/
int io_get_attr(io_Attr *a, char *name, bio_Data_Type *datatype, int *num, void **buffer)
{
    int m, attr_num, i, slen, slen_stord, n;
    int szdatatype_stord;
    char *pc, *c, *myname;
    long long a2[2];
    long long *pll, *pd;
    bio_Data_Type datatype_stord;

/*
   size,            szllong
   id,              szllong
   num,             szllong
   size of name,    szllong
   name,            slen*szchar
   datatype,        szllong
   size of value,   szllong
   values,          n*io_sizeof(datatype)
*/

    m = 3 * io_szllong;
        /* size_of_name, size_of_value, szdatatype */ 

    slen = strlen(name);
    *num = 0; 
    attr_num = *(a->num);
    pc = (char *) a->value + (3 * io_szllong);

    for (i = 0; i < attr_num; i++)
      { 
        /*****
        pll = (long long *) pc;
        slen_stord = *pll;

        c = (char *) (pll + 1);
        pd = (long long *) (c + slen_stord);
        datatype_stord = (bio_Data_Type) ((int)*pd);
        pll = pd + 1;
        n = (int) *pll;
        ****/
        memcpy(a2, pc, (size_t)io_szllong);
        slen_stord = (int) a2[0];
        myname = pc + io_szllong;
        c = myname + slen_stord;
        memcpy(a2, c, (size_t)(io_szllong + io_szllong));  
        datatype_stord = (bio_Data_Type) ((int)a2[0]);
        n = (int) a2[1];
        
        pc += (slen_stord + n + m);
        if (slen_stord != slen) continue;
        if (!strncmp(myname, name, (size_t)slen_stord))
           { 
             if (!(*buffer))
                { *buffer = malloc((size_t)n);
                 }
             c += (io_szllong + io_szllong);
             memcpy(*buffer, c, (size_t)n);
             szdatatype_stord = io_sizeof(datatype_stord);
             *num = n /szdatatype_stord;
             *datatype = datatype_stord; 

             break;
            }
        }
    /****
    if (!(*num)) 
       {  printf("ERROR: no attr %s in io_get_attr\n", name);
          return -1;
        }
    ****/
    return 0;
   }
/****************************************************************************************/
int io_get_object(io_File *file, int id, int *obj_idx, 
                  int *parent_id, bio_Object_Type *type, bio_Data_Type *datatype, 
                  long long *offset, long long *gsize, char **name)
{
    /* (1) offset and gsize are referred to the data chunck and in the unit of char  
       (2) If *name == NULL, *name will be allocated and should be cleaned after use. 
     */

/* 
          size_meta,          szllong 
          fid,                szllong 
          num,                szllong 
          size_data,          szllong
          size of name,       szllong 
          name,               slen*szchar
          id,                 szllong 
          parent_id,          szllong  
          obj_type,           szllong 
          datatype,           szllong  
          offset,             szllong
          size,               szllong
*/ 

    int m, i, num, myslen, id_stord;
    char *pc, *c, *myname;
    long long a6[6]; 
    long long *pll, *po, *pd;
    io_Obj *obj;

    m = 7 * io_szllong;
/*      size_of_name, id, parent_id, szobjtype, szdatatype,offset, size */  

    assert(file); 
        obj = file->object;
        num = *(obj->num);
        pc = (char *) obj->value + (4 * io_szllong);

        for (i = 0; i < num; i++)
          {  /****** 
             pll = (long long *) pc;
             slen = *pll;
             c = (char *)(pll + 1); 
             pll = (long long *)(c + slen);
             myid = (int) *pll;
             if (myid == id)
                { *parent_id = (int) pll[1];
                  po = pll + 2;
                  *type = (bio_Object_Type) ((int)*po);
                  pd = po + 1;
                  *datatype = (bio_Data_Type) ((int)*pd);
                  pll = pd + 1;
                  *offset = pll[0];
                  *gsize  = pll[1];
                  if (*name == NULL)  
                     { *name = (char *) malloc((size_t)(slen + 1));
                       memcpy(*name, c, (size_t) slen);
                       (*name)[slen] = '\0';
                      }
                  *obj_idx = i;
                  break;
                 }
             *****/

             memcpy(a6, pc, (size_t) io_szllong);
             myslen = a6[0];
             myname = pc + io_szllong;
             c = myname + myslen;
             memcpy(a6, c, (size_t)(6 * io_szllong));
             id_stord = (int) a6[0];
             if (id_stord == id) 
                { *parent_id = (int) a6[1];
                  *type = (bio_Object_Type) ((int)a6[2]);
                  *datatype = (bio_Data_Type) ((int)a6[3]);
                  *offset = a6[4];
                  *gsize  = a6[5];    
                  if (*name == NULL)
                     { *name = (char *) malloc((size_t)(myslen + 1));
                       memcpy(*name, myname, (size_t) myslen);
                       (*name)[myslen] = '\0';
                      }
                  *obj_idx = i;
                  break;
                 }
             pc += (myslen + m);
           }
       
    return 0;
   }    

io_Tree *io_find_file(int fileid, io_File **file)
{
   io_File *nextfile;
   io_Tree *root;

   *file = io_files;
   while (*file)
     { root = (*file)->root;
       if (root->id == fileid) { 
           return root;
       }
       nextfile = (*file)->next;
       (*file) = nextfile;
      }
   return NULL;
  }

io_Tree *io_find_tree(io_File *file, int obj_id)
{
   io_Tree *me;
   io_File *nextfile;

   assert(file);
   me = io_find_tree_helper(obj_id, file->root);
   return me;
 }

io_Tree *io_find_tree_helper(int obj_id, io_Tree *tree)
{
   io_Tree *me;
   io_Child *child, *nextchild;

   if (!tree) return NULL;

   if (tree->id == obj_id)
      { return tree;
       }
   else
      { child = tree->child;
        while (child)
           { if (!(child->tree))
                { printf("ERROR: null child->tree in io_find_tree_helper\n");
                  return NULL;
                 }
             me = io_find_tree_helper(obj_id, child->tree);
             if (me) return me;
             nextchild = child->next;
             child = nextchild;
            }
        return NULL;
      }
    return NULL;
  }
/****************************************************************************************/
int io_clean_tree(io_Tree *ptree)
{
   io_Tree *tree;
   io_Child *child, *next;

   child = ptree->child;
   while (child) {  
      tree = child->tree;
      if (tree) { 
          io_clean_tree(tree);
      }
      next = child->next;
      free(child);
      child = next; 
   }
   free(ptree);
  
   return 0;
} 
/****************************************************************************************/
int bio_attr_write(int fileid, int id, char *name, bio_Data_Type datatype, int size, void *buffer)
{   
    int obj_idx, parentid;
    bio_Object_Type type;
    bio_Data_Type dtype;
    long long offset, gsize, size_attr; 
    int k,  num_onj, slen, num_obj; 
    io_Obj  *obj;
    io_File *file;
    io_Attr *a; 
    io_Tree *tree;

    if (!name)
       { printf("ERROR: null name in bio_attr_write\n");
         return -1;
        }
    if (datatype == bio_datatype_invalid)
       { printf("ERROR: datatype invalid in bio_attr_write\n");
         return -1;
        }
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != id) {
        tree = io_find_tree(file, id);
    }
    if (!tree || !file) { 
        if (!mype) printf("ERROR: id not found in bio_attr_write, name = %s\n", name);
        return -1; 
    }
    obj_idx = tree->idx;

/*    io_get_object(file, id, &obj_idx, &parentid, &type, &dtype, &offset, &gsize, &name);
**/

    obj = file->object;
    a = obj->attr[obj_idx];  

    io_add_attr(a, name, datatype, size, buffer);

    return 0;
   }
/****************************************************************************************/
int bio_attr_read(int fileid, int id, char *name, bio_Data_Type *datatype, int *size, void **buffer)
{
    int obj_idx, parentid, err;
    bio_Object_Type type;
    bio_Data_Type dtype;
    long long offset, gsize;
    io_File *file;
    io_Attr   *a;
    io_Tree *tree;

    err = 0;
    if (!name)
       { printf("ERROR: null name in bio_attr_read\n");
         return -1;
        }
    if (!datatype) 
       { printf("ERROR: null datatype in bio_attr_read\n");
         return -1;
        }
    if (!size)
       { printf("ERROR: null size in bio_attr_read\n");
         return -1;
        }
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != id) {
        tree = io_find_tree(file, id);
    }
    if (!tree || !file) { 
        if (!mype) printf("ERROR: grp_id not found in bio_attr_read, name = %s\n", name);
        return -1; 
    }
    obj_idx = tree->idx;
    a = (file->object)->attr[obj_idx];
    err = io_get_attr(a, name, datatype, size, buffer);

    if (!(*size)) { 
         err = -2;
/*       printf("Warning: no attribute %s\n", name);  */ 
    } 
    return err;
   }
/****************************************************************************************/
int bio_get_name(int fileid, int id, bio_List_Struct *obj)
{ 
    char *c, *name;
    int slen, id_stord, parentid;
    bio_Object_Type type;
    bio_Data_Type datatype;
    long long a6[6], offset, gsize;
    io_File *file;
    io_Tree *tree;

    if (id < 0)  
       { printf("ERROR: id < 0 in bio_get_name\n");
         return -1;
        }
    if (obj == NULL) 
       { printf("ERROR: null obj in bio_get_name\n");
         return -1;
        }
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != id) {
        tree = io_find_tree(file, id);
    }
    if (!tree || !file) { 
        if (!mype) printf("ERROR: id not found in bio_get_name\n");
        return -1; 
    }
/*
          size_meta,          szllong
          fid,                szllong
          num,                szllong
          size_data,          szllong
          size of name,       szllong
          name,               slen*szchar
          id,                 szllong
          parent_id,          szllong
          obj_type,           szllong
          datatype,           szllong
          offset,             szllong
          size,               szllong
*/
    c = (char *)tree->p;
    memcpy(a6, c, (size_t)io_szllong);
    slen = a6[0];
    name = c + io_szllong;
    c = name + slen;
    memcpy(a6, c, (size_t)(6 * io_szllong));
    id_stord  = (int)a6[0];
    parentid = (int)a6[1];
    type = (bio_Object_Type)((int)a6[2]);
    datatype  = (bio_Data_Type) ((int)a6[3]);
    offset = a6[4];
    gsize  = a6[5];

    obj->name = (char *) malloc((size_t)(slen + 1));
    strncpy(obj->name, name, (size_t)slen);
    obj->name[slen] = '\0';
    obj->id = id;  
    obj->type = type;
    obj->datatype = datatype; 
      
/***
       name = NULL;
       io_get_object(file, id, &obj_idx, &parentid, &type, &datatype, &offset, &gsize, &name);
***/
    
    if (datatype == bio_char) 
       { obj->size = gsize;
        }
    else if (datatype == bio_int) 
       { obj->size = gsize / io_szint;
        }
    else if (datatype == bio_long) 
       { obj->size = gsize / io_szlong;
        }
    else if (datatype == bio_long_long) 
       { obj->size = gsize / io_szllong;
        }
    else if (datatype == bio_float) 
       { obj->size = gsize / io_szfloat;
        }
    else if (datatype == bio_double) 
       { obj->size = gsize / io_szdouble;
        }
    else if (datatype == bio_long_double) 
       { obj->size = gsize / io_szldouble;
        }   
    else 
       { obj->size = 0;
        } 
    obj->value = NULL;
    obj->narrays = 0;
    obj->names = NULL;
    obj->sizes = NULL;
   
    return 0;   
 }
/****************************************************************************************/
int io_buffer_null(io_Buffer *buffer)
{ 
    buffer->id = -1;
    buffer->grp_id = -1;
    buffer->narrays = -1;
    buffer->max_narrays = -1;
    buffer->datatypes = NULL;

    buffer->collective = 1;
    buffer->max_nbufs = 0;
    buffer->nbufs_written = 0;
    buffer->max_size_in_name = 0;
    buffer->size_in_name = 0;
    buffer->max_size_in_helper = 0;
    buffer->max_size_in_data = 0;
    buffer->size_in_data = 0;
    buffer->tot_size_in_data = 0;

    buffer->name = NULL;
    buffer->file = NULL;
    buffer->narrays_buf = NULL;
    buffer->names  = NULL;
    buffer->helper = NULL;
    buffer->helper_for_read = NULL;
    buffer->data = NULL;
    
    buffer->npes  = 0;
    buffer->datatypes_pe_buf = NULL;
    buffer->narrays_pe_buf = NULL;
    buffer->sizes_pe = NULL;
    buffer->data_pe_buf = NULL;

    buffer->next = NULL;

    return 0;
 } 
/****************************************************************************************/
int io_buffer_clean(io_Buffer *buffer)
{
    int n, k, b, i, collective;
    int ncpu, nbufs; 
    long long *helper, *nbufs_pe;
    char ***data_pe_buf;

    ncpu = buffer->npes;
    nbufs = buffer->nbufs_written;
    collective = buffer->collective;
    helper = buffer->helper_for_read; 
    data_pe_buf = buffer->data_pe_buf;

    if (collective == 1) { 
        if (data_pe_buf) {
            for (k = 0; k < ncpu; k++) {
                for (i = 0; i < nbufs; i++) {
                    if (data_pe_buf[k][i]) free(data_pe_buf[k][i]);
                }
            }
            if (data_pe_buf[0]) free(data_pe_buf[0]);
            free(data_pe_buf);
            buffer->data_pe_buf = NULL;
        }
    }   
    else if (collective == 0) {  

/*      structure of helper:
           header (ie, npes, collective, tot_nbufs, tot_size_in_data, tot_size_in_name, tot_narrays, datatype)

           part1:
                 nbufs_pe0,        nbufs_pe1,        ...
                 narrays_pe0,      narrays_pe1,      ...
                 size_in_data_pe0, size_in_data_pe1, ,,,
                 size_in_name_pe0, size_in_name_pe1, ...

           pe0:    narrays_buf0, narrays_buf1, ....
           pe1:    narrays_buf0, narrays_buf1, ....
            ....
           npes-1: narrays_buf0, narrays_buf1,

           pe0: size_array0, size_array1, .....
           pe1: size_array1, size_array1, .....
            ....
           last_pe: size_array0, aize_array1, ...
*/
        if (helper) { 
/***
            nbufs_pe = helper + io_nshared_in_ibuf;
            for (k = 0; k < ncpu; k++) { 
                nbufs = nbufs_pe[k];
                for (b = 0; b < nbufs; b++) { 
                    if (data_pe_buf[k][b]) { 
                        free(data_pe_buf[k][b]);
                        data_pe_buf[k][b] = NULL;
                    }
                }
            }
***/
        }
    }
    if (buffer->narrays_pe_buf) { 
        free(buffer->narrays_pe_buf);
        /* each element points to some part of helper, no free(narrays_pe_buf[k]) */
        buffer->narrays_pe_buf = NULL;
    }
    if (buffer->sizes_pe) { 
        free(buffer->sizes_pe);
        /* each element points to some part of helper, no free(sizes_pe[k]) */
        buffer->sizes_pe = NULL;
    } 
    if (buffer->data_pe_buf) { 
        if (buffer->data_pe_buf[0]) free(buffer->data_pe_buf[0]);
        free(buffer->data_pe_buf);
        buffer->data_pe_buf = NULL;
    }
    if (buffer->name)            free(buffer->name);
    if (buffer->narrays_buf)     free(buffer->narrays_buf);
    if (buffer->names)           free(buffer->names);
    if (buffer->helper)          free(buffer->helper);
    if (buffer->helper_for_read) free(buffer->helper_for_read);
    if (buffer->data)            free(buffer->data);
    if (buffer->datatypes)       free(buffer->datatypes); 

    buffer->name            = NULL;
    buffer->narrays_buf     = NULL;
    buffer->names           = NULL;
    buffer->helper          = NULL;
    buffer->helper_for_read = NULL;
    buffer->data            = NULL;
    buffer->datatypes       = NULL;
    
    return 0;
 }
/****************************************************************************************/
int bio_buffer_init(int fileid, int grp_id, char *name, 
                    int collective, long long nbytes, int *buffer_id) 
{ 
    const int max_size_in_name = 3200;
    const int max_size_in_data = 0; 
    const int max_size_in_helper = 400; 
    const int max_narrays = 200;  
    const int max_nbufs = 50;

    long long nchar, tot_size;
    int n, obj_idx, parentid, i, k, slen, which_buffer;
//  int *narrays;
    int *narrays_buf;
    bio_Object_Type obj_type;
    bio_Data_Type dtype;
    long long offset, gsize;
    io_Obj *obj;
    io_File *file;
    io_Tree *tree; 
    io_Buffer *buffer;

    if (!name) { 
       printf("ERROR: null name in bio_buffer_init\n");
       return -1;
    }
    slen = strlen(name);
    for (i = 0; i < slen; i++)
      { if (name[i] == '/')
           { printf("ERROR: '/' in name not allowd in bio_buffer_init\n");
             return -1;
            }
       }
    if (!buffer_id)
       { printf("ERROR: null buffer_id in bio_buffer_init\n");
         return -1;
        }
    if (nbytes < max_size_in_data) { 
       tot_size = max_size_in_data;
    }
    else { 
       tot_size = nbytes;
    } 
    file = NULL;
    tree = io_find_file(fileid, &file);
    io_get_object(file, grp_id, &obj_idx, &parentid, &obj_type,
                   &dtype, &offset, &gsize, &name);
    if (obj_type != bio_group)
       { printf("ERROR: grp_id is not for a group in bio_buffer_init\n");
         return -1;
        }
    buffer = (io_Buffer *) malloc(sizeof(io_Buffer));
    io_buffer_null(buffer);

    buffer->name = (char *) malloc((size_t)(slen+1));
    strcpy(buffer->name, name);
    buffer->id = io_next_id;
    io_next_id++;
    buffer->grp_id = grp_id;

    buffer->max_nbufs = max_nbufs;
    n = max_nbufs * npes;
/****
    narrays = (int *) malloc((size_t)(n * io_szint));
    for (i = 0; i < n; i++) {  
      narrays[i] = 0;
    } 
*****/
    narrays_buf = (int *) malloc((size_t)(max_nbufs * io_szint));
    for (i = 0; i < max_nbufs; i++) {
      narrays_buf[i] = 0;
    }
    buffer->narrays_buf = narrays_buf;
    buffer->nbufs_written = 0;

    buffer->file     = file;  
    file->buffer     = buffer;

    buffer->narrays  = 0;
    buffer->max_narrays = (long long) max_narrays;
    buffer->datatypes = (bio_Data_Type *)malloc(max_narrays * sizeof(bio_Data_Type)); 

    buffer->collective = collective; 
    buffer->max_size_in_name = max_size_in_name;
    buffer->size_in_name = 0;
    buffer->names = (char *)malloc((size_t)max_size_in_name);
    buffer->max_size_in_helper = max_size_in_helper;
    buffer->helper = (long long *)malloc((size_t)(buffer->max_size_in_helper * io_szllong));

    buffer->max_size_in_data = tot_size;
    buffer->size_in_data = 0;
    buffer->tot_size_in_data = 0;
    buffer->data = malloc((size_t)(tot_size)); 
    assert(buffer->data); 

    buffer->next = io_buffers;
    io_buffers = buffer;
    
    *buffer_id = buffer->id;

    if (file->id_min > buffer->id) { 
        file->id_min = buffer->id;  
    }
//    bio_group_open_i(grp_id, name, buffer_id); 

    return 0;
 } 
/****************************************************************************************/
int bio_write(int fileid, int grp_id, char *name, 
             bio_Object_Type type, bio_Data_Type datatype,
             long long *offset, long long size, long long *gsize, 
             void *array, int *array_id)
{ 
/*  This function may be used to copy an array into a buffer. For this case,
    grp_id is an id of a buffer, and gsize is ignored which may be set to
    zero before the call. 

    This function will reset array_id for writing an array, 
    but ignore array_id for buffer. 
*/  
    int err, myid;

    assert((type != bio_group) && (type != bio_object_type_undefined));
    if (!array_id) {  
        printf("ERROR: null array_id in bio_write\n");
        return -1;
    }
    myid = -1;

    err = bio_write_i(fileid,
            grp_id, name, type, datatype, offset, size, gsize, array, &myid);

    if (myid > 0) { 
        *array_id = myid;
    }        
    return 0;
   }   
/****************************************************************************************/
int bio_write_i(int fileid, int grp_id, char *name, 
             bio_Object_Type type, bio_Data_Type datatype,
             long long *offset, long long size, long long *gsize, 
             void *array, int *array_id)
{ 
/*  This function is the same as bio_write except that this function will use array_id
    as an input if *array_id > 0. 
*/  
    int obj_idx, parentid, i, myslen, id_stord, gparent_id;
    char *c, *myname;
    bio_Object_Type obj_type;
    bio_Data_Type dtype;
    long long a6[6], myoffset, mygsize, filesize_written; 
    io_Obj *obj;
    io_File *file; 
    io_Buffer *buffer, *next;
    io_Tree *tree;
    
    assert((type != bio_group) && (type != bio_object_type_undefined));

    if (!name)
       { printf("ERROR: null name in bio_write_i\n");
         return -1;
        }
    if (!array_id)
       { printf("ERROR: null array_id in bio_write_i\n");
         return -1;
        }
    /*****
    myslen = strlen(name);
    for (i = 0; i < myslen; i++)
      { if (name[i] == '/')
           { printf("ERROR: '/' in name not allowd in bio_write_i\n");
             return -1;
            }
       } 
    *****/
    if (size > 0 && !array)
       { printf("ERROR: null array in bio_write_i\n");
         return -1;
        }
    if (datatype == bio_datatype_invalid)
       { printf("ERROR: invalid datatype in bio_write_i\n");
         return -1;
        }
    if (type == bio_group || type == bio_object_type_undefined)
       { printf("ERROR: type incorrect in bio_write_i\n");
         return 0;
        }
    file = NULL;
    tree = io_find_file(fileid, &file);

//  buffer = io_buffers;   This is more safe
    buffer = file->buffer;
    while (buffer)
       { if (buffer->id == grp_id) break; 
         next = buffer->next;
         buffer = next;
        }
    if (buffer)
       { if (!buffer->helper) {
             printf("ERROR: null buffer->helper in bio_write_i\n");
             return -1;
         }
         obj_type = (bio_Object_Type) ((int)type);
         if (buffer->collective) { 
             if (*offset < 0) { 
                 io_offset_gsize_get(size, offset, gsize);
             } 
             myoffset = *offset;
         }
         else { 
             myoffset = 0;
         } 
         io_buffer_copy(file, buffer, name, obj_type, datatype, myoffset, size, array);
        }
    else 
       { if (*offset < 0 || *gsize < 0) { 
            io_offset_gsize_get(size, offset, gsize);
         }   
         else if (*offset + size > *gsize)
            { printf("ERROR: gsize too small in bio_write_i\n");
              return -1;
             }
         file = NULL;
         tree = io_find_file(fileid, &file);
         if (fileid != grp_id) {
            tree = io_find_tree(file, grp_id);
         }
         if (!tree || !file) { 
             if (!mype) printf("ERROR: grp_id not found in bio_write_i\n");
             return -1;         } 
         c = (char *)tree->p;
         memcpy(a6, c, (size_t)io_szllong);
         myslen = a6[0];
         myname = c + io_szllong;
         c = myname + myslen;
         memcpy(a6, c, (size_t)(6 * io_szllong));
         id_stord = (int) a6[0];
         gparent_id = (int)a6[1];
         obj_type = (bio_Object_Type) ((int)a6[2]);

         if (obj_type != bio_group) {
             if (!mype) 
                 printf("ERROR: gparent_id is not for a group in bio_write_i, name = %s\n",
                 name);
             return -1;
         }  
         io_add_object(1, &tree,
                       file, name, type, datatype, array_id, 
                       *offset, size, *gsize, array);
        } 
    if (flag_for_resilience == 1) { 
        filesize_written = io_file_flush(file); 
    } 
    return 0;
   }   
/****************************************************************************************/
int bio_write_i_old(int fileid, int grp_id, char *name, 
             bio_Object_Type type, bio_Data_Type datatype,
             long long *offset, long long size, long long *gsize, 
             void *array, int *array_id)
{ 
/* This function is the same as bio_write,  but it will use array_id in
   in the input.  
*/  
    int obj_idx, parentid, i, myslen, id_stord, gparent_id;
    char *c, *myname;
    bio_Object_Type obj_type;
    bio_Data_Type dtype;
    long long a6[6], myoffset, mygsize; 
    io_Obj *obj;
    io_File *file; 
    io_Buffer *buffer, *next;
    io_Tree *tree;
    
    if (!name)
       { printf("ERROR: null name in bio_write\n");
         return -1;
        }
    if (!array_id)
       { printf("ERROR: null array_id in bio_write\n");
         return -1;
        }
    if (*array_id < 1) { 
         printf("ERROR: array_id < 1 in  bio_write_i\n");
         return -1;
    }
    myslen = strlen(name);
    for (i = 0; i < myslen; i++)
      { if (name[i] == '/')
           { printf("ERROR: '/' in name not allowd in bio_write\n");
             return -1;
            }
       } 
    if (size > 0 && !array)
       { printf("ERROR: null array in bio_write\n");
         return -1;
        }
    if (datatype == bio_datatype_invalid)
       { printf("ERROR: invalid datatype in bio_write\n");
         return -1;
        }
    if (type == bio_group || type == bio_object_type_undefined)
       { printf("ERROR: type incorrect in bio_write\n");
         return 0;
        }
    file = NULL;
    tree = io_find_file(fileid, &file);

//  buffer = io_buffers;   This is more safe 
    buffer = file->buffer;
    while (buffer)
       { if (buffer->id == grp_id) break; 
         next = buffer->next;
         buffer = next;
        }
    if (buffer)
       { if (!buffer->helper) {
             printf("ERROR: null buffer->helper in bio_write\n");
             return -1;
         }
         obj_type = (bio_Object_Type) ((int)type);
         io_buffer_copy(file, buffer, name, obj_type, datatype, *offset, size, array);
        }
    else 
       { if (*offset < 0 || *gsize < 0) { 
            io_offset_gsize_get(size, offset, gsize);
         }    
         else if (*offset + size > *gsize)
            { printf("ERROR: gsize too small in bio_write\n");
              return -1;
             }
         file = NULL;
         tree = io_find_file(fileid, &file);
         if (fileid != grp_id) {
             tree = io_find_tree(file, grp_id);
         }
         if (!tree || !file) { 
             if (!mype) printf("ERROR: grp_id not found in bio_write\n");
             return -1; 
         }
         c = (char *)tree->p;
         memcpy(a6, c, (size_t)io_szllong);
         myslen = a6[0];
         myname = c + io_szllong;
         c = myname + myslen;
         memcpy(a6, c, (size_t)(6 * io_szllong));
         id_stord = (int) a6[0];
         gparent_id = (int)a6[1];
         obj_type = (bio_Object_Type) ((int)a6[2]);

         if (obj_type != bio_group) {
             if (!mype) printf("ERROR: gparent_id is not for a group in bio_write\n");
             return -1;
         }
         io_add_object(1, &tree,
                       file, name, type, datatype, array_id, 
                       *offset, size, *gsize, array);
        } 
    return 0;
   }   
/****************************************************************************************/
int io_offset_gsize_get(long long size, long long *offset, long long *gsize)
{
    int k; 
    long long *sizes;
  
#ifdef MPI
    MPI_Status status;

     sizes = (long long *) malloc((size_t)(npes * io_szllong));
     MPI_Allgather(&size, (int)1, MPI_LONG_LONG_INT, 
                   sizes, (int)1, MPI_LONG_LONG_INT, bio_comm); 
     *offset = 0;
     for (k = 0; k < mype; k++) { 
         *offset += sizes[k];
     }
     *gsize = *offset;
     for (k = mype; k < npes; k++) { 
         *gsize += sizes[k];
     }
     free(sizes);
#else 
    *gsize = size;
    *offset = 0; 
#endif

    return 0;
} 
/****************************************************************************************/
int io_buffer_copy(io_File *file,
                   io_Buffer *buffer, char *name, bio_Object_Type obj_type, 
                   bio_Data_Type datatype, long long offset, long long size, void *array)
{
    int  mycommit, to_commit;
    int  i, signal, slen, k, max_nbuf, max_narrays, commited;
    long long n, nchar, max_size_in_helper;
    char *myname, *names, *pc;
    char extendedname[MAXLNAME];
    void *data, *p;
    long long *helper, *pll;
    bio_Data_Type *datatypes;
    
    nchar = size * io_sizeof(datatype);

    if (buffer->size_in_data + nchar > buffer->max_size_in_data) { 
        mycommit = 1;
    }
    else { 
        mycommit = 0;
    }
    if (npes > 1) { 
#ifdef MPI
//      MPI_Allreduce(&mycommit, &to_commit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&mycommit, &to_commit, 1, MPI_INT, MPI_MAX, bio_comm);
#endif 
    } 
    else {   
        to_commit = mycommit;
    } 
    if (to_commit) {    /* this is newly added to avoid the hang, 6/20/2012 */  
        if (buffer->collective) { 
            io_buffer_commit(file, buffer, &commited);
        }
        else { 
            io_buffer_commit_nonc(file, buffer, &commited);
        } 
    }
    if (buffer->size_in_data + nchar > buffer->max_size_in_data) { 
        buffer->max_size_in_data = buffer->size_in_data + nchar + 8;
        data = (void *) malloc(buffer->max_size_in_data * sizeof(char));
        memcpy(data, buffer->data, (size_t)buffer->size_in_data);
        free(buffer->data);
        buffer->data = data;
    } 
    p = (void *)((char *) buffer->data + buffer->size_in_data);

    memcpy(p, array, (size_t)nchar);
    buffer->size_in_data += nchar;

    if (buffer->collective) { 
        myname = name;
    }
    else {
        sprintf(extendedname,"%s$%d", name, mype);
        myname = extendedname;
    }
    slen = strlen(myname);
    n = buffer->size_in_name + (slen + 1 + io_szint);
    if (n >= buffer->max_size_in_name)
       { names = (char *) malloc((size_t)(n + n));
         memcpy(names, buffer->names, (size_t)buffer->size_in_name);
         buffer->max_size_in_name = n + n;
         if (buffer->names) free(buffer->names);
         buffer->names = names;
        }
    pc = buffer->names + buffer->size_in_name;
    strcpy(pc, myname);
    pc += (slen + 1);
    i = (int) obj_type;
    memcpy(pc, &i, (size_t)io_szint); 
    buffer->size_in_name += (slen + 1 + io_szint);
    
    if (buffer->collective) { 
        n = buffer->narrays + buffer->narrays;
        if (n + 2 > buffer->max_size_in_helper)
           { max_size_in_helper = n + n;
             helper = (long long *) malloc((size_t) (max_size_in_helper * io_szllong));
             memcpy(helper, buffer->helper, (size_t)(n * io_szllong));
             free(buffer->helper);
             buffer->helper = helper;
             buffer->max_size_in_helper = max_size_in_helper; 
            }
        buffer->helper[n] = offset;
        buffer->helper[n+1] = size;
    }
    else {   
        n = buffer->narrays;
        if (n >= buffer->max_size_in_helper)
           { max_size_in_helper = n + n;
             helper = (long long *) malloc((size_t) (max_size_in_helper * io_szllong));
             memcpy(helper, buffer->helper, (size_t)(n * io_szllong));
             free(buffer->helper);
             buffer->helper = helper;
             buffer->max_size_in_helper = max_size_in_helper; 
            }
        buffer->helper[n] = size;
    }
    if (buffer->narrays >= buffer->max_narrays) { 
        max_narrays = buffer->narrays + buffer->narrays;
        datatypes = (bio_Data_Type *)malloc(max_narrays * sizeof(bio_Data_Type));
        memcpy(datatypes, buffer->datatypes, (size_t)(buffer->narrays * sizeof(bio_Data_Type)));
        free(buffer->datatypes);
        buffer->datatypes = datatypes;
        buffer->max_narrays = max_narrays;
    }
    buffer->datatypes[buffer->narrays] = datatype;

    (buffer->narrays)++;

    return 0;
   }    
/****************************************************************************************/
int io_buffer_commit(io_File *file, io_Buffer *buffer, int *commited)
{  
    int  i, k, k3, n, narrays_commited, narrays, id;
    int  max_nbufs, max_nbufs_old;
    int  *narrays_buf;
    char buffername[MAXLNAME];
    long long offset, size, gsize;
    long long sendbuf[3], *recvbuf, *sizes_pe;
    char *ptr;
    io_Tree *tree_grp;
 
    recvbuf = (long long *) malloc((size_t)(4 * npes) * io_szllong);
    sizes_pe = recvbuf + 3 * npes; 

#ifdef MPI
    sendbuf[0] = buffer->size_in_data;
    sendbuf[1] = io_next_id;
    sendbuf[2] = buffer->nbufs_written;
    MPI_Gather(sendbuf, (int)3, MPI_LONG_LONG_INT, 
               recvbuf, (int)3, MPI_LONG_LONG_INT, 0, bio_comm); 
#else 
    recvbuf[0] = buffer->size_in_data;
    recvbuf[1] = io_next_id;
    recvbuf[2] = buffer->nbufs_written;
#endif 
    if (!mype) {  
        id = io_next_id;
        for (k = 0; k < npes; k++) { 
            k3 = k + k + k + 1;
            if ((int)recvbuf[k3] > id) id = recvbuf[k3];
        }
        for (k = 0; k < npes; k++) { 
            k3 = k + k + k;
            if (!recvbuf[k3]) { 
                recvbuf[k3+1] = -1;
            }
            else { 
                recvbuf[k3+1] = id;
                id++;
            } 
        }
    }  
#ifdef MPI
    MPI_Bcast(recvbuf, (int)(npes + npes + npes), MPI_LONG_LONG_INT, (int)0, bio_comm);
#endif
    id = io_next_id; 
    gsize = 0;
    for (k = 0; k < npes; k++) { 
        k3 = k + k + k;
        sizes_pe[k] = recvbuf[k3]; 
        gsize += recvbuf[k3];
        if (id < recvbuf[k3+1]) id = recvbuf[k3+1]; 
    } 
    if (!gsize) { 
        *commited = 0;
        free(recvbuf);
        return 0;
    }  
    *commited = 1;
    io_next_id = id + 1;  
 
    if (buffer->nbufs_written + 2 >= buffer->max_nbufs) { 
        max_nbufs_old = buffer->max_nbufs;
        max_nbufs = max_nbufs_old + max_nbufs_old;
        narrays_buf = (int *) malloc((size_t)(max_nbufs * io_szint));
 
        memcpy(narrays_buf, buffer->narrays_buf, (size_t)(io_szint * max_nbufs_old));
        for (i = max_nbufs_old; i < max_nbufs; i++) {
            narrays_buf[i] = 0;
        }
        free(buffer->narrays_buf);
        buffer->narrays_buf = narrays_buf;
        buffer->max_nbufs   = max_nbufs;
    }
    offset = 0;
    for (i = 0; i < mype; i++) {
        offset += sizes_pe[i];
    }
    size = sizes_pe[mype];
    buffer->tot_size_in_data += gsize; 

    narrays_buf = buffer->narrays_buf;
    n = buffer->nbufs_written; 
    narrays_commited = 0;
    for (i = 0; i < n; i++) { 
        narrays_commited += narrays_buf[i];
    }
    narrays = buffer->narrays - narrays_commited;
    narrays_buf[n] = narrays;

    if (buffer->nbufs_written == 0) {  
        sprintf(buffername,"%s%s", io_buffer_pre,buffer->name);
        id = buffer->id;
    }
    else {   
        sprintf(buffername,"%s%d$%s", 
                io_buffer_pre,buffer->nbufs_written,buffer->name);
    } 
/*  id doesn't need to be reset in the following io_add_object */

    tree_grp = io_find_tree(file, buffer->grp_id);

    io_add_object(1, &tree_grp,
                 buffer->file, buffername, bio_buffer, bio_char,
                 &id, offset, size, gsize, buffer->data);

    (buffer->nbufs_written)++;             
    buffer->size_in_data = 0;

    free(recvbuf);
 
    return 0;
 }
/****************************************************************************************/
int io_buffer_commit_nonc(io_File *file, io_Buffer *buffer, int *commited)
{
    int  n, to_write, i, k, k3, id, idmax;
    int  max_nbufs, max_nbufs_old, nbufs_written;
    int  *narrays_buf, *nrrays_commitedarrays_buf;
    char buffername[MAXLNAME];
    long long narrays, narrays_commited, size, gsize;
    long long sendbuf[3], *recvbuf;
    io_Tree *tree_grp;

    tree_grp = io_find_tree(file, buffer->grp_id);
 
        recvbuf = (long long *) malloc((size_t)(npes + npes + npes) * io_szllong);
#ifdef MPI
        sendbuf[0] = buffer->size_in_data;
        sendbuf[1] = io_next_id;
        sendbuf[2] = buffer->nbufs_written;
        MPI_Gather(sendbuf, (int)3, MPI_LONG_LONG_INT, 
                   recvbuf, (int)3, MPI_LONG_LONG_INT, 0, bio_comm); 
#else 
        recvbuf[0] = buffer->size_in_data;
        recvbuf[1] = io_next_id;
        recvbuf[2] = buffer->nbufs_written;
#endif 
        if (!mype) {  
            id = io_next_id;
            for (k = 0; k < npes; k++) { 
                k3 = k + k + k + 1;
                if ((int)recvbuf[k3] > id) id = recvbuf[k3];
            }
            for (k = 0; k < npes; k++) { 
                k3 = k + k + k;
                if (!recvbuf[k3]) { 
                    recvbuf[k3+1] = -1;
                }
                else { 
                    recvbuf[k3+1] = id;
                    id++;
                } 
            }
        }  
#ifdef MPI
        MPI_Bcast(recvbuf, (int)(npes + npes + npes), MPI_LONG_LONG_INT, (int)0, bio_comm);
#endif
    
    idmax = io_next_id; 
    gsize = 0;
    for (k = 0; k < npes; k++) { 
        k3 = k + k + k;
        gsize += recvbuf[k3];
        if (idmax < recvbuf[k3+1]) idmax = recvbuf[k3+1]; 
    } 
    if (!gsize) { 
        *commited = 0;
        if (recvbuf) free(recvbuf);
        return 0;
    }  
    *commited = 1;
    io_next_id = idmax + 1;  

    if (buffer->nbufs_written + 2 >= buffer->max_nbufs) { 
        max_nbufs_old = buffer->max_nbufs;
        max_nbufs = max_nbufs_old + max_nbufs_old;
        narrays_buf = (int *) malloc((size_t)(max_nbufs * io_szint));
 
        memcpy(narrays_buf, buffer->narrays_buf, (size_t)(io_szint * max_nbufs_old));
        for (i = max_nbufs_old; i < max_nbufs; i++) {
            narrays_buf[i] = 0;
        }
        free(buffer->narrays_buf);
        buffer->narrays_buf = narrays_buf;
        buffer->max_nbufs = max_nbufs;
    }
    narrays_buf = buffer->narrays_buf;
    n = buffer->nbufs_written;
    narrays_commited = 0;
    for (i = 0; i < n; i++)
      { narrays_commited += narrays_buf[i];
       }
    narrays = buffer->narrays - narrays_commited;
    narrays_buf[n] = narrays;

    for (k = 0; k < npes; k++) { 
        k3 = k + k + k;
        size = recvbuf[k3];
        if (!size) continue;

        /** io_ibuffer$10$name  or io_ibuffer_$10$77$name  */

        nbufs_written = recvbuf[k3 + 2];
        if (nbufs_written == 0) { 
            sprintf(buffername,"%s%d$%s", io_ibuffer_pre, k, buffer->name);
            id = buffer->id;
        }
        else {  
            sprintf(buffername,"%s%d$%d$%s",
                    io_ibuffer_pre, k, nbufs_written,buffer->name);
        }
        if (mype == k) { 
            to_write = 1;
        }
        else { 
            to_write = 0;
        }
        io_add_object(to_write, &tree_grp,
               buffer->file, buffername, bio_buffer, bio_char,
               &id, 0ll, size, size, buffer->data);
         
        if (k == mype) { 
           (buffer->nbufs_written)++;             
            buffer->size_in_data = 0;
            buffer->tot_size_in_data += size; 

        }
    }
    free(recvbuf);

    return 0;
} 
/****************************************************************************************/
int bio_buffer_finalize(int fileid, int buffer_id)
{ 
    int err;
    io_Buffer *buffer, *next, *previous;
    io_File *file;

    err = 0;
    previous = NULL;
    buffer = io_buffers;  
    while (buffer) {  
      if (buffer_id == buffer->id) break;   
      next = buffer->next;
      previous = buffer;
      buffer = next;
    }  
    if (!buffer) {  
       printf("ERROR: buffer_id not found in bio_buffer_finalize\n");
       return -1;
    } 
    file = buffer->file;
    if (file->mode == bio_file_create) { 
        if (buffer->collective) { 
            err = io_buffer_finalize(file, buffer);
        }
        else { 
            err = io_buffer_finalize_nonc(file, buffer);  
        }
    } 
    io_buffer_clean(buffer);
    if (previous) {
        previous->next = buffer->next;
    }
    else {
        io_buffers = buffer->next;
    }
    free(buffer); 
    file->buffer = NULL;
    buffer = io_buffers;
    while (buffer && (file->buffer == NULL)) {
         if ((buffer->file)->id == file->id) {
              file->buffer = buffer;
         }
         buffer = buffer->next;
    }
    return err;
} 
/****************************************************************************************/
int io_buffer_finalize_nonc(io_File *file, io_Buffer *buffer)
{ 
    int ndbg;
    int i, k, ioffset, id, narrays, commited;
    int *recvcount, *displs, *recvbuf, *narrays_buf;
    long long size, offset, gsize;
    long long nbufs, nbufs_max, dnbufs;
    long long tot_narrays, tot_nbufs, tot_size_in_data, tot_size_in_name;
    
    long long sendbuf[4];
    long long *llrecvbuf;
    long long *nbufs_pe, *narrays_pe, *size_in_data_pe, *size_in_name_pe;
    long long *t, *s;
    long long *helper;
    char *namelist, name[MAXLNAME];
    bio_Data_Type *datatypes;
    io_Tree *tree_grp;

    tree_grp = io_find_tree(file, buffer->grp_id);

    if (!tree_grp) { 
        printf("ERROR: null tree_grp in io_buffer_finalize_nonc\n");
        return -1;
    }
    ndbg = 0;
    commited = 1;
    while (commited) { 
        io_buffer_commit_nonc(file, buffer, &commited);
        ndbg++;
    }
    llrecvbuf = (long long *) malloc(8 * npes * io_szllong);
    nbufs_pe = llrecvbuf  + (4 * npes);
    narrays_pe = nbufs_pe + npes;
    size_in_data_pe = narrays_pe + npes;
    size_in_name_pe = size_in_data_pe + npes;
    
    narrays = buffer->narrays;
    datatypes = buffer->datatypes;
    nbufs   = buffer->nbufs_written;
  
#ifdef MPI
    sendbuf[0] = nbufs;                    /* total number of bufs of this pe */
    sendbuf[1] = narrays;                  /* total number of arrays on this pe */
    sendbuf[2] = buffer->tot_size_in_data; /* total size in this pe */
    sendbuf[3] = buffer->size_in_name; 
   
    MPI_Allgather(sendbuf, (int)4, MPI_LONG_LONG_INT, 
                  llrecvbuf, (int)4, MPI_LONG_LONG_INT, bio_comm); 
#else
    llrecvbuf[0] = nbufs;                    /* total number of bufs of this pe */
    llrecvbuf[1] = narrays;                  /* total number of arrays on this pe */
    llrecvbuf[2] = buffer->tot_size_in_data; /* total size in this pe */
    llrecvbuf[3] = buffer->size_in_name;
#endif   
    t = llrecvbuf;
    for (k = 0; k < npes; k++) { 
        nbufs_pe[k] = t[0];
        narrays_pe[k] = t[1];
        size_in_data_pe[k] = t[2];
        size_in_name_pe[k] = t[3];
        t += 4;
    }  
    tot_narrays = 0;
    tot_nbufs   = 0;
    tot_size_in_data = 0;
    tot_size_in_name = 0;
    for (k = 0; k < npes; k++) { 
        tot_nbufs += nbufs_pe[k];
        tot_narrays += narrays_pe[k];
        tot_size_in_data += size_in_data_pe[k];
        tot_size_in_name += size_in_name_pe[k];
    }        
    if (!tot_nbufs) { 
        free(llrecvbuf);
        return 0;
    } 
    nbufs     = nbufs_pe[mype];
    recvcount = (int *) malloc((npes + npes + tot_nbufs) * io_szint);
    displs    = recvcount + npes;
    recvbuf   = displs    + npes; 

    ioffset = 0; 
    for (k = 0; k < npes; k++) {
        recvcount[k] = nbufs_pe[k]; 
        displs[k] = ioffset;
        ioffset += nbufs_pe[k];
    }
#ifdef MPI
    MPI_Allgatherv(buffer->narrays_buf, (int)(nbufs), MPI_INT, recvbuf, recvcount, displs, MPI_INT, bio_comm); 
#else 
    memcpy(recvbuf, buffer->narrays_buf, (size_t)(nbufs * io_szint));
#endif
   /*  structure of helper:

       header (ie, npes, collective, tot_nbufs, tot_size_in_data, tot_size_in_name, tot_narrays)  
    
       part1: 
           nbufs_pe0,   nbufs_pe1,        ... 
           narrays_pe0, narrays_pe1       ...

           size_in_data_pe0, size_in_data_pe1, ,,, 
           size_in_name_pe0, size_in_name_pe1, ...

       pe0:    narrays_buf0, narrays_buf1, ....
       pe1:    narrays_buf0, narrays_buf1, ....
       ....
       npes-1: narrays_buf0, narrays_buf1, .... 

       pe0:  datatypes of narrays,  size_array0, size_array1, .....
       pe1:  datatypes of narrays,  size_array1, size_array1, .....
       ....
       last_pe: size_array0, aize_array1, ...
*/
    gsize  = (io_nshared_in_ibuf + 4 * npes + tot_nbufs + tot_narrays + tot_narrays);  
    helper = (long long *) malloc(gsize * io_szllong);

    helper[0] = npes;
    helper[1] = 0;
    helper[2] = tot_nbufs;
    helper[3] = tot_size_in_data;
    helper[4] = tot_size_in_name;
    helper[5] = tot_narrays;

    for (i = 6; i < io_nshared_in_ibuf; i++) {
        helper[i] = 0;
    }
    t = helper + io_nshared_in_ibuf;
    for (k = 0; k < npes; k++) { 
        t[k] = nbufs_pe[k];
    }
    t += npes;
    for (k = 0; k < npes; k++) { 
        t[k] += narrays_pe[k];
    }   
    t += npes;
    
    for (k = 0; k < npes; k++) { 
        t[k] = size_in_data_pe[k];
    }   
    t += npes;
    for (k = 0; k < npes; k++) { 
        t[k] = size_in_name_pe[k];
    }   
    t += npes;
    ioffset = 0; 
    for (k = 0; k < npes; k++) { 
        nbufs = nbufs_pe[k];
        narrays_buf = recvbuf + ioffset;
        for (i = 0; i < nbufs; i++) { 
            t[i] = narrays_buf[i];
        }
        t += nbufs;
        ioffset += nbufs; 
    } 
    if (npes == 1) {  
        for (k = 0; k < narrays; k++) { 
            t[k] = (long long)((int)datatypes[k]);
        }
        t += narrays;
        memcpy(t, buffer->helper, (size_t)(narrays * io_szllong));
    } 
    else { 
        s = (long long *)malloc((narrays + narrays + 1) * io_szllong);
        for (k = 0; k < narrays; k++) { 
            s[k] = (long long)((int)datatypes[k]);
        }
        memcpy(s + narrays, buffer->helper, (size_t)(narrays * io_szllong));
        ioffset = 0;
        for (k = 0; k < npes; k++) { 
            displs[k] = ioffset;
            i = narrays_pe[k] + narrays_pe[k];            
            ioffset  += i;
            recvcount[k] = i;
        } 
#ifdef MPI

        MPI_Allgatherv(s, (int)narrays, MPI_LONG_LONG_INT, 
                       t, recvcount, displs, MPI_LONG_LONG_INT, bio_comm); 
#endif
        free(s);
    }  
    if (gsize) { 
        sprintf(name,"%s%s", io_ibufhelper_pre, buffer->name);
        bio_attr_write(file->id, tree_grp->id, name, bio_long_long, (int) gsize, helper);
    } 
    if (helper) free(helper);

    /* write the name list  */

    sprintf(name,"%s%s", io_ibufnames_pre, buffer->name);

    if (npes == 1) {  
        bio_attr_write(file->id, tree_grp->id, name, bio_char, (int)tot_size_in_name, buffer->names);
    }
    else {   
        namelist = (char *) malloc(tot_size_in_name + 1);  
        ioffset = 0;
        for (k = 0; k < npes; k++) { 
            displs[k] = ioffset;
            ioffset  += size_in_name_pe[k];
            recvcount[k] = size_in_name_pe[k];
        } 
#ifdef MPI
        MPI_Allgatherv(buffer->names, (int)size_in_name_pe[mype], MPI_CHAR, 
                       namelist, recvcount, displs, MPI_CHAR, bio_comm); 
      /***
        MPI_Allgatherv(buffer->names, (int)size_in_name_pe[mype], MPI_CHAR,
                       namelist, recvcount, displs, MPI_CHAR, 0, bio_comm);
        ***/
#endif

        bio_attr_write(file->id, tree_grp->id, name, bio_char, (int)tot_size_in_name, namelist);
        free(namelist);
    }
    free(llrecvbuf);
    free(recvcount);

    return 0;
  }  

/****************************************************************************************/
/* This routine writes helper and names as datasets, but the one above writes them as attrs. */

int io_buffer_finalize_nonc_old(io_File *file, io_Buffer *buffer)
{ 
    int ndbg;
    int i, k, k3, ioffset, id, commited;
    int *recvcount, *displs, *recvbuf, *narrays_buf;
    long long size, offset, gsize;
    long long narrays, nbufs, nbufs_max, dnbufs;
    long long tot_narrays, tot_nbufs, tot_size_in_data, tot_size_in_name;
    
    long long sendbuf[4];
    long long *llrecvbuf;
    long long *nbufs_pe, *narrays_pe, *size_in_data_pe, *size_in_name_pe;
    long long *t, *s;
    long long *helper;
    char name[MAXLNAME];
    io_Tree *tree_grp;

    tree_grp = io_find_tree(file, buffer->grp_id);
  
    ndbg = 0;
    commited = 1;
    while (commited) { 
        io_buffer_commit_nonc(file, buffer, &commited);
        ndbg++;
    }
    llrecvbuf = (long long *) malloc(8 * npes * io_szllong);
    nbufs_pe = llrecvbuf  + (4 * npes);
    narrays_pe = nbufs_pe + npes;
    size_in_data_pe = narrays_pe + npes;
    size_in_name_pe = size_in_data_pe + npes;
    
    narrays = buffer->narrays;
    nbufs   = buffer->nbufs_written;
  
#ifdef MPI
    sendbuf[0] = nbufs;                    /* total number of bufs of this pe */
    sendbuf[1] = narrays;                  /* total number of arrays on this pe */
    sendbuf[2] = buffer->tot_size_in_data; /* total size in this pe */
    sendbuf[3] = buffer->size_in_name; 
   
    MPI_Allgather(sendbuf, (int)4, MPI_LONG_LONG_INT, 
                  llrecvbuf, (int)4, MPI_LONG_LONG_INT, bio_comm); 
#else
    llrecvbuf[0] = nbufs;                    /* total number of bufs of this pe */
    llrecvbuf[1] = narrays;                  /* total number of arrays on this pe */
    llrecvbuf[2] = buffer->tot_size_in_data; /* total size in this pe */
    llrecvbuf[3] = buffer->size_in_name;
#endif   
    t = llrecvbuf;
    for (k = 0; k < npes; k++) { 
        nbufs_pe[k] = t[0];
        narrays_pe[k] = t[1];
        size_in_data_pe[k] = t[2];
        size_in_name_pe[k] = t[3];
        t += 4;
    }  
    tot_narrays = 0;
    tot_nbufs   = 0;
    tot_size_in_data = 0;
    tot_size_in_name = 0;
    for (k = 0; k < npes; k++) { 
        tot_nbufs += nbufs_pe[k];
        tot_narrays += narrays_pe[k];
        tot_size_in_data += size_in_data_pe[k];
        tot_size_in_name += size_in_name_pe[k];
    }        
    nbufs     = nbufs_pe[mype];
    recvcount = (int *) malloc((npes + npes + tot_nbufs) * io_szint);
    displs    = recvcount + npes;
    recvbuf   = displs    + npes; 

    ioffset = 0; 
    for (k = 0; k < npes; k++) {
        recvcount[k] = nbufs_pe[k]; 
        displs[k] = ioffset;
        ioffset += nbufs_pe[k];
    }
#ifdef MPI
    MPI_Allgatherv(buffer->narrays_buf, (int)(nbufs), MPI_INT, recvbuf, recvcount, displs, MPI_INT, bio_comm); 
#else 
    memcpy(recvbuf, buffer->narrays_buf, (size_t)(nbufs * io_szint));
#endif

   /*  structure of helper:
       header (ie, npes, collective, tot_nbufs, tot_size_in_data, tot_size_in_name, tot_narrays, datatype)  
    
       part1: 
           nbufs_pe0,        nbufs_pe1,        ... 
           narrays_pe0,      narrays_pe1,      ...
           size_in_data_pe0, size_in_data_pe1, ,,, 
           size_in_name_pe0, size_in_name_pe1, ...

       pe0:    narrays_buf0, narrays_buf1, ....
       pe1:    narrays_buf0, narrays_buf1, ....
       ....
       npes-1: narrays_buf0, narrays_buf1, 

       pe0: size_array0, size_array1, .....
       pe1: size_array1, size_array1, .....
       ....
       last_pe: size_array0, aize_array1, ...
*/
    gsize = (io_nshared_in_ibuf + 4 * npes + tot_nbufs + tot_narrays);  
    
    if (mype == 0) { 
        offset = 0;
        size   = io_nshared_in_ibuf + 4 * npes + tot_nbufs + narrays_pe[mype];
        helper = (long long *) malloc(size * io_szllong);
        helper[0] = npes;
        helper[1] = 0;
        helper[2] = tot_nbufs;
        helper[3] = tot_size_in_data;
        helper[4] = tot_size_in_name;
        helper[5] = tot_narrays;
        for (i = 6; i < io_nshared_in_ibuf; i++) {
            helper[i] = 0;
        }
        t = helper + io_nshared_in_ibuf;
        for (k = 0; k < npes; k++) { 
            t[k] = nbufs_pe[k];
        }
        t += npes;
        for (k = 0; k < npes; k++) { 
            t[k] += narrays_pe[k];
        }   
        t += npes;
        for (k = 0; k < npes; k++) { 
            t[k] = size_in_data_pe[k];
        }   
        t += npes;
        for (k = 0; k < npes; k++) { 
            t[k] = size_in_name_pe[k];
        }   
        t += npes;
        ioffset = 0; 
        for (k = 0; k < npes; k++) { 
            nbufs = nbufs_pe[k];
            narrays_buf = recvbuf + ioffset;
            for (i = 0; i < nbufs; i++) { 
                t[i] = narrays_buf[i];
            }
            t += nbufs;
            ioffset += nbufs; 
        } 
        memcpy(t, buffer->helper, (size_t)(narrays_pe[mype] * io_szllong));
    } 
    else { 
        size   = narrays_pe[mype];
        offset = io_nshared_in_ibuf + 4 * npes + tot_nbufs;
        for (k = 0; k < mype; k++) { 
            offset += narrays_pe[k];
        }
        helper = buffer->helper;
    }  
    if (gsize) {  
       sprintf(name,"%s%s", io_ibufhelper_pre, buffer->name);

       id = buffer->id;   /* buffer helper is using the id of the buffer */   

/*     id doesn't have to be reset in the following io_add_object   */
                         
       io_add_object(1, &tree_grp,
                     buffer->file, name, bio_buffer, bio_long_long,
                     &id, offset, size, gsize, helper);
    }
    if ((mype == 0) && helper ) {
        free(helper);
        helper = NULL;
    }

   /*  write the name list  */
 

    offset = 0;
    for (k = 0; k < mype; k++) { 
        offset += size_in_name_pe[k];
    } 
    if (tot_size_in_name) { 
        sprintf(name,"%s%s", io_ibufnames_pre, buffer->name);

        id = io_next_id;  
        io_next_id++;

        /*  id doesn't have to be reset in the following io_add_object  */  

        io_add_object(1, &tree_grp,
                      buffer->file, name, bio_buffer, bio_char,
                      &id, offset, size_in_name_pe[mype],
                      tot_size_in_name, buffer->names);
    }
   
    free(llrecvbuf);
    free(recvcount);

    return 0;
  }  
/****************************************************************************************/
int io_buffer_finalize(io_File *file, io_Buffer *buffer)
{ 
    int ndbg, commited;
    int i, k, n, m, nbuf, npesnb, id, narrays, flag;
    int ioffset;
    int *narrays_buf, *narrays_buf_allpe;
    long long size, offset, gsize;
    long long *helper, *t;
    char name[MAXLNAME];
    char *ptr;
    bio_Data_Type *datatypes; 
    io_Tree *tree_grp;

    tree_grp = io_find_tree(file, buffer->grp_id);
  
    if (!tree_grp) { 
        printf("ERROR: null tree_grp in io_buffer_finalize\n");
        return -1;
    }
    ndbg = 0;
    commited = 1;
    while (commited) { 
        io_buffer_commit(file, buffer, &commited);
        ndbg++;
    }
    nbuf   = buffer->nbufs_written;
    npesnb = npes * nbuf; 
    if (!nbuf) return 0;

    narrays_buf = (int *) buffer->narrays_buf;
    narrays_buf[nbuf] = buffer->size_in_name;  /* the space was ensured in bio_buffer_commit*/

    narrays_buf_allpe = (int *) malloc((size_t)((npes * nbuf + npes) * io_szint));

    /* npes in narrays_buf_allpe is for size_in_names of all pes */

    /* 1 is for size_in_name  in case with collective = 0 
       narrays_buf_allpe is available only to pe 0 if buffer->collective != 0  */

    n = nbuf + 1;
#ifdef MPI
    MPI_Gather(narrays_buf,n,MPI_INT,narrays_buf_allpe,n,MPI_INT,0,bio_comm);
#else 
    offset = n * mype;
    memcpy(narrays_buf_allpe + offset, narrays_buf, (size_t)(n * io_szint)); 
#endif
    
/*  structure of helper:
       header:                                   \ 
       narrays numbers for datatypes              \ 
       npes numbers for size_in_name_each_pe       \ 
       pe0: narrays_buf0, narrays_buf1, ....        \  m as follows
       pe1: narrays_buf0, narrays_buf1, ....        /  
       ....                                        /   
       lastpe: narrays_buf0, narrays_buf1, ...    / 

       pe0:  (offset, size) of array0, (offset, size) of array1, ....
       pe1:  (offset, size) of array0, (offset, size) of array1, ....
       ... 
       pelast: (offset, size) of array0, (offset, size) of array1, ....
*/
    narrays = buffer->narrays;
    m = io_nshared_in_buf + narrays + npes + npes * nbuf; /* narrays for datatypes,
                                                             npes is for size_in_name */
    size = narrays + narrays; 
    gsize = (long long)(m + npes * size);
   
    helper = (long long *) malloc((size_t)(gsize * io_szllong));
    for (i = 0; i < io_nshared_in_buf; i++) {
        helper[i] = -1;
    }
    helper[0] = (long long) npes;
    helper[1] = (long long) buffer->collective; 
    helper[2] = (long long) narrays;
    helper[3] = (long long) buffer->nbufs_written;
    helper[4] = buffer->tot_size_in_data;

    t = helper + io_nshared_in_buf;
    datatypes = buffer->datatypes;
    for (k = 0; k < narrays; k++) { 
        t[k] = (long long) ((int) datatypes[k]);
    }
    t += narrays;
    ioffset = nbuf; 
    for (k = 0; k < npes; k++) {
        t[k] = (long long) narrays_buf_allpe[ioffset]; /* size_in_name */
        ioffset += (nbuf + 1); 
    } 
    t += npes; 
    ioffset = 0;
    for (k = 0; k < npes; k++) {
        narrays_buf = narrays_buf_allpe + ioffset;
        for (i = 0; i < nbuf; i++) {
            t[i] = narrays_buf[i];
        }
        t += nbuf; 
        ioffset += (nbuf + 1);
    }  
    t = helper + m;
    if (npes == 1) { 
        memcpy(t, buffer->helper, (size_t)(io_szllong * size));
    } 
    else { 

#ifdef MPI
        MPI_Allgather(buffer->helper, (int)size, MPI_LONG_LONG_INT, 
                      t, (int)size, MPI_LONG_LONG_INT, bio_comm); 
       /****
        MPI_Gather(buffer->helper, (int)size, MPI_LONG_LONG_INT,
                      t, (int)size, MPI_LONG_LONG_INT, 0, bio_comm);
        ***/
#endif
    }
    if (gsize) {  
        sprintf(name,"%s%s", io_bufhelper_pre, buffer->name);
        bio_attr_write(file->id, tree_grp->id, name, bio_long_long, (int)gsize, helper);

        /**  swrite the name list    **/  

        /** if wanted, it may be activated 
        if (mype == 0) {
            ioffset = nbuf;
            for (k = 0; k < npes; k++) {
                if (buffer->size_in_name != narrays_buf_allpe[ioffset]) {
                    printf("ERROR: collective, but names different in bio_buffer_finalize\n");
                    return -1;
                }
                ioffset += (nbuf + 1);
            }
        }
        ****/
        sprintf(name,"%s%s", io_bufnames_pre, buffer->name);
        bio_attr_write(file->id, tree_grp->id, name, bio_char, (int)buffer->size_in_name, buffer->names);
    } 
    free(helper);
    free(narrays_buf_allpe);
   
    return 0;
  }  
/****************************************************************************************/
/****
  This writes helper and names as datasets, but the one above writes them as atttrs 
***/
int io_buffer_finalize_old(io_File *file, io_Buffer *buffer)
{ 
#ifdef MPI
    MPI_Status status;
    MPI_Request *reqs_recv;
    int *signals_recvd; 
#endif
    int ndbg, commited;
    int i, k, n, m, nbuf, npesnb, id, narrays, flag;
    int ioffset;
    int *narrays_buf, *narrays_buf_allpe;
    long long size, offset, gsize;
    long long *helper, *llp;
    char name[MAXLNAME];
    char *ptr;
    io_Tree *tree_grp;

    tree_grp = io_find_tree(file, buffer->grp_id);
  
    ndbg = 0;
    commited = 1;
    while (commited) { 
        io_buffer_commit(file, buffer, &commited);
        ndbg++;
    }
    nbuf   = buffer->nbufs_written;
    npesnb = npes * nbuf; 

    narrays_buf = (int *) buffer->narrays_buf;
    narrays_buf[nbuf] = buffer->size_in_name;  /* the space was ensured in bio_buffer_commit*/

    narrays_buf_allpe = (int *) malloc((size_t)((npes * nbuf + npes) * io_szint));

    /* npes in narrays_buf_allpe is for size_in_names of all pes */

/*  1 is for size_in_name  in case with collective = 0 
    narrays_buf_allpe is available only to pe 0 if buffer->collective != 0*/

    n = nbuf + 1;
#ifdef MPI
    MPI_Gather(narrays_buf,n,MPI_INT,narrays_buf_allpe,n,MPI_INT,0,bio_comm);
#else 
    offset = n * mype;
    memcpy(narrays_buf_allpe + offset, narrays_buf, (size_t)(n * io_szint)); 
#endif
    
/*  structure of helper:
       header                                     \ 
       npes numbers for size_in_name_each_pe       \ 
       pe0: narrays_buf0, narrays_buf1, ....        \  m as follows
       pe1: narrays_buf0, narrays_buf1, ....        /  
       ....                                        /   
       lastpe: narrays_buf0, narrays_buf1, ...    / 

       pe0:  (offset, size) of array0, (offset, size) of array1, ....
       pe1:  (offset, size) of array0, (offset, size) of array1, ....
       ... 
       pelast: (offset, size) of array0, (offset, size) of array1, ....
*/
    m = io_nshared_in_buf + npes + npes * nbuf; /* npes is for size_in_name */

    narrays = buffer->narrays;
    offset = 0;
    if (mype) {
       offset = (long long)(m  + mype *(narrays + narrays));
    }
    gsize = (long long)(m + npes *(narrays + narrays));
   
    if (mype == 0) {  
       /* npes, collective, narrays, nbuffers, tot_size_in_data, 
          datatype, nbuf_each_pe, empty_spaces(-1), 
          size_in_name of each pe,
          (offset, size) of array1
          (offset, size) of array2
           ....
          (offset, size) of array_narrays 
       */
       size = (long long) (m + buffer->narrays + buffer->narrays);
       helper = (long long *) malloc((size_t)(gsize * io_szllong));
       for (i = 0; i < io_nshared_in_buf; i++) {
         helper[i] = -1;
       }
       helper[0] = (long long) npes;
       helper[1] = (long long) buffer->collective; 
       helper[2] = (long long) narrays;
       helper[3] = (long long) buffer->nbufs_written;
       helper[4] = buffer->tot_size_in_data;

       llp = helper + io_nshared_in_buf;
       ioffset = nbuf; 
       for (k = 0; k < npes; k++) {
           llp[k] = (long long) narrays_buf_allpe[ioffset]; /* size_in_name */
           ioffset += (nbuf + 1); 
       } 
       llp += npes; 
       ioffset = 0;
       for (k = 0; k < npes; k++) {
           narrays_buf = narrays_buf_allpe + ioffset;
           for (i = 0; i < nbuf; i++) {
               llp[i] = narrays_buf[i];
           }
           llp += nbuf; 
           ioffset += (nbuf + 1);
       }  
       memcpy(helper+m, buffer->helper, (size_t)(io_szllong *(buffer->narrays + buffer->narrays)));
    }
    else { 
       size = (long long)(buffer->narrays + buffer->narrays);
       helper = buffer->helper; 
    } 
/******* 
     

***/

    if (gsize) {  
       sprintf(name,"%s%s", io_bufhelper_pre, buffer->name);

       id = buffer->id;   /* buffer helper is using the id of the buffer */   

/*     id doesn't have to be reset in the following io_add_object   */
                         
       io_add_object(1, &tree_grp,
                     buffer->file, name, bio_buffer, bio_long_long,
                     &id, offset, size, gsize, helper);
    }
/*  write the name list:
    if collective != 0: only pe 0 writes names.
    Otherwise, all pes write buffer->names, and the number of names is all
    added together.
 */

    offset = 0;
    if (mype == 0) {
       size = (long long) buffer->size_in_name;
    }
    else {
       size = 0;
    }
    gsize = (long long) buffer->size_in_name;
    if (mype == 0) {
        ioffset = nbuf;
        for (k = 0; k < npes; k++) {
            if (buffer->size_in_name != narrays_buf_allpe[ioffset]) {
                printf("ERROR: collective, but names different in bio_buffer_finalize\n");
                return -1;
            }
            ioffset += (nbuf + 1);
        }
    }
    if (gsize) {
       sprintf(name,"%s%s", io_bufnames_pre, buffer->name);

       id = io_next_id;  
       io_next_id++;

/*     id doesn't have to be reset in the following io_add_object  */  

       io_add_object(1, &tree_grp,
                     buffer->file, name, bio_buffer, bio_char,
                     &id, offset, size, gsize, buffer->names);
    }
    if (mype == 0) {
       free(helper);
    }
    return 0;
  }  
/****************************************************************************************/
int io_get_buffer_info(long long *helper, 
                    long long ***my_narrays_pe_buf, long long ***my_sizes_pe,
                    char ****my_data_pe_buf, bio_Data_Type ****my_datatypes_pe_buf)
{
  /* for collective buffer 
` 
     structure of helper:
     header                                    \
     narrays numbers for datatypes              \ 
     npes numbers for size_in_name_each_pe       \
     pe0: narrays_buf0, narrays_buf1, ....        \  m as follows
     pe1: narrays_buf0, narrays_buf1, ....        /
     ....                                        /
     lastpe: narrays_buf0, narrays_buf1, ...    /

     pe0:  (offset, size) of array0, (offset, size) of array1, ....
     pe1:  (offset, size) of array0, (offset, size) of array1, ....
     ... 
     pelast: (offset, size) of array0, (offset, size) of array1, ....
 */ 
    int ncpu, k, b, collective, offset, m, narrays, narrays_tot, npesnb, nbufs;
    long long i, tot_size, size;
    long long *lldatatypes, *t, *s;
    bio_Data_Type *datatypes, ***datatypes_pe_buf;

    long long **narrays_pe_buf, **sizes_pe;
    char ***data_pe_buf;

    ncpu       = (int) helper[0];
    collective = (int) helper[1];
    narrays_tot  = helper[2];
    nbufs    = helper[3];
    tot_size = helper[4];
 
    narrays_pe_buf = (long long **) malloc((size_t) ncpu * sizeof(long long*));
    npesnb = (long long) ncpu * nbufs;
    narrays_pe_buf[0] = helper + (io_nshared_in_buf + narrays_tot + ncpu);
    for (k = 1; k < ncpu; k++) {
        narrays_pe_buf[k] = narrays_pe_buf[k-1] + nbufs;
    }
    *my_narrays_pe_buf = narrays_pe_buf;

    lldatatypes = helper + io_nshared_in_buf;
    datatypes = (bio_Data_Type *)malloc(narrays_tot * sizeof(bio_Data_Type));
    for (k = 0; k < narrays_tot; k++) {
        datatypes[k] = (bio_Data_Type)((int)lldatatypes[k]);
    }
    datatypes_pe_buf = (bio_Data_Type ***)malloc(ncpu * sizeof(bio_Data_Type **));

    for (k = 0; k < ncpu; k++) { 
        offset = 0;
        datatypes_pe_buf[k] = (bio_Data_Type **)malloc(nbufs * sizeof(bio_Data_Type *));
        for (b = 0; b < nbufs; b++) { 
            datatypes_pe_buf[k][b] = datatypes + offset;
            offset += narrays_pe_buf[k][b]; 
        }
    }  
    *my_datatypes_pe_buf = datatypes_pe_buf;
        
    sizes_pe = (long long **) malloc((size_t) ncpu * sizeof(long long*));

    /* sizes_pe[k] is the starting address of the meta data (offset, size)
       of the 1st array on k-th pe.
    */
    
    m = io_nshared_in_buf + narrays_tot + ncpu + ncpu * nbufs; /* narrays for datatypes,
                                                              npes is for size_in_name */
    
//    sizes_pe[0] = narrays_pe_buf[0] + ncpu * nbufs;
    sizes_pe[0] = helper + m;

    for (k = 1; k < ncpu; k++) { 
        narrays = 0;
        for (i = 0; i < nbufs; i++) {
            narrays += narrays_pe_buf[k-1][i];
        }
        sizes_pe[k] = sizes_pe[k-1] + (narrays + narrays);
    }
    *my_sizes_pe = sizes_pe;

    data_pe_buf = (char ***) malloc((size_t) ncpu * sizeof(char **));
    data_pe_buf[0] = (char **) malloc((size_t)npesnb * sizeof(char*));
    for (i = 0; i < npesnb; i++) { 
        data_pe_buf[0][i] = NULL;
    }
    for (k = 1; k < ncpu; k++) {
        data_pe_buf[k] = data_pe_buf[k-1] + nbufs;
    }
    *my_data_pe_buf = data_pe_buf;

    return 0;
} 
/****************************************************************************************/
int io_get_buffer_info_nonc(long long *helper, 
                    long long ***my_narrays_pe_buf, long long ***my_sizes_pe, 
                    char ****my_data_pe_buf, bio_Data_Type ****my_datatypes_pe_buf)
{
    /*  structure of helper:
        header (ie, npes, collective, tot_nbufs, tot_size_in_data,
                tot_size_in_name, tot_narrays, datatype)

    part1:
        nbufs_pe0,        nbufs_pe1,        ...
        narrays_pe0,      narrays_pe1,      ...
        size_in_data_pe0, size_in_data_pe1, ,,,
        size_in_name_pe0, size_in_name_pe1, ...

        pe0:    narrays_buf0, narrays_buf1, ....
        pe1:    narrays_buf0, narrays_buf1, ....
        ....
        npes-1: narrays_buf0, narrays_buf1,

        pe0:     datatypes of narrays,  size_array0, size_array1, .....
        pe1:     datatypes of narrays,  size_array1, size_array1, .....
        ....
        last_pe: datatypes of narrays,  size_array0, aize_array1, ...
    */

    int ncpu, k, b, collective;
    long long i, narrays, nbufs, tot_nbufs, tot_narrays;
    long long tot_size_in_data, tot_size_in_name, offset;
    long long *nbufs_pe, *narrays_pe, *t, *s;
    bio_Data_Type ***datatypes_pe_buf, *datatypes, datatype;

    long long *lldatatypes;
    long long **narrays_pe_buf, **sizes_pe;
    char ***data_pe_buf;

    ncpu = helper[0];
    tot_nbufs = helper[2];
    tot_size_in_data = helper[3];
    tot_size_in_name = helper[4];
    tot_narrays      = helper[5];

    nbufs_pe = helper + io_nshared_in_ibuf;
    narrays_pe = nbufs_pe + ncpu;

    narrays_pe_buf = (long long **)malloc((size_t)(ncpu) * sizeof(long long *));
    narrays_pe_buf[0] = helper + (io_nshared_in_ibuf + 4 * ncpu);
    for (k = 1; k < ncpu; k++) {
        nbufs = nbufs_pe[k-1];
        narrays_pe_buf[k] = narrays_pe_buf[k-1] + nbufs;
    }
    *my_narrays_pe_buf = narrays_pe_buf;

    lldatatypes = helper + ((long long)(io_nshared_in_ibuf + 4 * ncpu) + tot_nbufs);
    datatypes = (bio_Data_Type *) malloc(tot_narrays * sizeof(bio_Data_Type));
    for (k = 0; k < tot_narrays; k++) { 
        datatypes[k] = (bio_Data_Type)((int)lldatatypes[k]);
    }
    datatypes_pe_buf = (bio_Data_Type ***)malloc(ncpu * sizeof(bio_Data_Type **));
    offset = 0;
    for (k = 0; k < ncpu; k++) {
        nbufs = nbufs_pe[k];
        datatypes_pe_buf[k] = (bio_Data_Type **)malloc(nbufs * sizeof(bio_Data_Type *));
        for (b = 0; b < nbufs; b++) {
            datatypes_pe_buf[k][b] = datatypes + offset;
            offset += narrays_pe_buf[k][b];
        }
    }
    *my_datatypes_pe_buf = datatypes_pe_buf;
 
    /** datatypes_pe_buf[0][0] points to datatypes **/ 

    sizes_pe = (long long **)malloc(ncpu * sizeof(long long *));
    sizes_pe[0] = lldatatypes + narrays_pe[0];
     
    for (k = 1; k < ncpu; k++) { 
        sizes_pe[k] = sizes_pe[k-1] + (narrays_pe[k-1] + narrays_pe[k]);
    } 
    *my_sizes_pe = sizes_pe;                   /* in the unit of datatype, e.x, double */ 

    data_pe_buf = (char ***) malloc((size_t) ncpu * sizeof(char **));
    data_pe_buf[0] = (char **) malloc((size_t)tot_nbufs * sizeof(char*));
    for (i = 0; i < tot_nbufs; i++) { 
        data_pe_buf[0][i] = NULL;
    }
    for (k = 1; k < ncpu; k++) {
        nbufs = nbufs_pe[k-1];
        data_pe_buf[k] = data_pe_buf[k-1] + nbufs;
    }
    *my_data_pe_buf = data_pe_buf;

    return 0;
} 
/****************************************************************************************/
int io_buffer_getname(char *name, char *purename, int *pe)
{
    int k, s;
    char *c, *cpe;
 
    cpe = NULL;

    strcpy(purename, name);
    s = strlen(purename);
    for (k = s-1; k >= 0; k--) {
        if (purename[k] == '$') {
            purename[k] = '\0';
            cpe = &(purename[k+1]); 
            break;
        }
    }
    if (cpe) { 
       *pe = atoi(cpe);
    }
    else { 
       *pe = -1;
    }
    return 0;
}

/****************************************************************************************/
int bio_read(int fileid, int grp_id, char *name, bio_Data_Type datatype,
            long long offset, long long size, void **array,
            int *array_id)
{ 
    int err; 
    io_Buffer *buffer, *next;
    io_File *file; 
    io_Tree *tree;
    
    file = NULL;
    buffer = io_buffers;
    while (buffer) {
        if (buffer->id == grp_id) {
            break;
        }
        next   = buffer->next;
        buffer = next;
    }
    if (buffer) { 
        file = buffer->file;
    }
    else { 
        tree = io_find_file(fileid, &file);
    }
    err = io_read(file, buffer, grp_id, name, datatype,
            offset, size, array, array_id); 

    return err;
   }   
/****************************************************************************************/
int io_read(io_File *file, io_Buffer *buffer, 
            int grp_id, char *name, bio_Data_Type datatype,
            long long offset, long long size, void **array,
            int *array_id)
{ 
    int obj_idx, id_stord, parentid, cpu, err, failed, found;
    bio_Object_Type type_stord;
    long long a6[6], myoffset, mysize, mygsize; 
    int i, j, slen, myslen;
    char *c, *myname;
    io_Obj *obj;
    io_Tree *tree;
    
    char *namelist;
    char bufname[MAXLNAME], purename[MAXLNAME];
    void *data, *t, *s;
    long long start, remaining, size_to_get, offset_this_array, offset_narray; 
    long long offset_in_this_pe, size_this_pe_buf, offset_this_pe_buf, nchars;
    long long **narrays_pe_buf, **sizes_pe, *offsets_array;
    long long *offsets_array_other_pe;
    char ***data_pe_buf;
    int  iallocated, tot_nbufs, nbufs, mynbufs; 
    int  n, na, myna, na_offset, dn, pid, ncpu, narrays;
    int  collective, k, b, kk, bb, myid;
    long long *helper, *nbufs_pe, *narrays_buf, *narrays_pe, *sizes, *llp;
    bio_Data_Type datatype_stord, *datatypes, ***datatypes_pe_buf;

    err = 0;

    *array_id = -1;
    iallocated = 0;

    if (size == 0) return 0;
    else if (size < 0) 
       { printf("ERROR: size < 0 in bio_read\n");
         return -1;
        }
    if (!name)
       { printf("ERROR: null name in bio_write\n");
         return -1;
        }
    slen = strlen(name);
   
    if (offset < 0) 
       { printf("ERROR: offset < 0 in bio_read\n");
         return -1;
        } 
    if (datatype == bio_datatype_invalid)
       { printf("ERROR: invalid datatype in bio_read\n");
         return -1;
        }
    if (!(*array))  { 
        nchars = size * io_sizeof(datatype);
        *array = malloc(nchars);
        iallocated = 1;
    } 
    if (buffer) {
        namelist = buffer->names;
        narrays_pe_buf = buffer->narrays_pe_buf;
        sizes_pe = buffer->sizes_pe;

        /* for collective buffer,
              sizes_pe[k] is the starting address of the meta data (offset, size) ... 
           for non-collective buffer,
              sizes_pe[k] is the starting address of meta data size, size, ....
         */ 
        data_pe_buf = buffer->data_pe_buf;
        
        if (!namelist || !narrays_pe_buf || !sizes_pe) { 
            printf("ERROR: bio_list must be called before reading an array in buffer\n");
            if (iallocated) { 
                free(*array); 
                *array = NULL;
            }
            return -1;
        }
        pid  = buffer->grp_id;
        
        ncpu = buffer->npes;
        tot_nbufs = buffer->nbufs_written;
        narrays = buffer->narrays;
        collective = buffer->collective;
        datatypes_pe_buf = buffer->datatypes_pe_buf;
        helper = buffer->helper_for_read;

        myname = namelist;
        n = 0;
        while (n < narrays) { 
            slen = strlen(myname);
            if (collective == 1) { 
                if (!strcmp(myname, name)) {
                    break;
                }
            }
            else if (collective == 0) { 
                io_buffer_getname(myname, purename, &cpu);
                if (!strcmp(purename, name)) {
                    break;
                }
            }
            else { 
                printf("ERROR: collective in bio_read\n");
                if (iallocated) { 
                    free(*array);
                    *array = NULL;
                }
                return -1;
            } 
            myname += (slen + 1 + io_szint);  /* there is an int after each name,
                                                 which is obj_type of the array */
            n++;
        } 
        if (n >= narrays) { 
//          printf("ERROR: name not found in the buffer\n");
            if (iallocated) {
                free(*array);
                *array = NULL;
            }
            return -1;
        }
        if (collective == 1) { 
            start = offset;
            remaining = size;
       
            for (k = 0; k < ncpu; k++) { 
                offset_narray = 0;
                for (b = 0; b < tot_nbufs; b++) { 
                    na = narrays_pe_buf[k][b];
                    if ((n >= offset_narray + na) || (n < offset_narray)) { 
                       offset_narray += na;                      
                       continue;
                    } 
                    datatypes = datatypes_pe_buf[k][b];
                    size_this_pe_buf = 0;                                          /* in terms of bytes */
                    offsets_array = sizes_pe[k] + (offset_narray + offset_narray); /* for this b */
                    for (i = 0; i < na; i++) { 
                        size_this_pe_buf += (offsets_array[i+i+1] * io_sizeof(datatypes[i]));  /* offset, size */ 
                    }
                /** offset_narray <=  n < (offset_narray + na)       **/

                    dn = n - offset_narray;
                    datatype_stord = datatypes[dn];

                    myoffset = offsets_array[dn + dn];
                    mysize   = offsets_array[dn + dn + 1];
                    if ((start + remaining <= myoffset) || (myoffset + mysize <= start)) continue;

                    if (!data_pe_buf[k][b] || !flag_for_buffered_read) { 

                        tree = io_find_tree(file, pid);
                        if (!tree) {
                            printf("ERROR: pid not found in bio_write\n");
                            if (iallocated) {
                                free(*array);
                                *array = NULL;
                            }
                            return -1;
                        }
                        /* free other allocated spaces */

                        if (flag_for_buffered_read) { 
                            for (kk = 0; kk < ncpu; kk++) { 
                                 for (bb = 0; bb < tot_nbufs; bb++) { 
                                     if (bb == b) continue;    /* keep the same buffers on other PEs, in case the data
                                                                  requested across different PEs.
                                                                */       
                                     if (buffer->data_pe_buf[kk][bb]) { 
                                         free(buffer->data_pe_buf[kk][bb]);
                                         buffer->data_pe_buf[kk][bb] = NULL;
                                     } 
                                 }
                            }
                            data_pe_buf[k][b] = malloc(size_this_pe_buf);
    
                            if (!data_pe_buf[k][b]) { 
                                /* if malloc failed, free other previously allocated spaces */
    
                                failed = 1;
                                while (failed) {  
                                   for (kk = 0; kk < ncpu; kk++) { 
                                        for (bb = 0; bb < tot_nbufs; bb++) { 
                                            if (buffer->data_pe_buf[kk][bb]) { 
                                                free(buffer->data_pe_buf[kk][bb]);
                                                buffer->data_pe_buf[kk][bb] = NULL;
                                                data_pe_buf[k][b] = malloc(size_this_pe_buf);
                                                if (data_pe_buf[k][b]) { 
                                                    failed = 0;
                                                    break;
                                                }
                                            }
                                        } 
                                        if (!failed) break;
                                    }
                                }
                                if (failed) { 
                                    printf("ERROR: malloc failed in XXX\n");
                                    if (iallocated) {
                                        free(*array);
                                        *array = NULL;
                                    }
                                    return -1;
                                } 
                            }
                        }
                        else {
                            for (kk = 0; kk < ncpu; kk++) {  
                                for (bb = 0; bb < tot_nbufs; bb++) {  
                                    if (buffer->data_pe_buf[kk][bb]) { 
                                        free(buffer->data_pe_buf[kk][bb]);
                                        buffer->data_pe_buf[kk][bb] = NULL;
                                    } 
                                } 
                            } 
                        }
                        if (b == 0) {
                            sprintf(bufname,"%s%s", io_buffer_pre, buffer->name);
                        }
                        else { 
                            sprintf(bufname,"%s%d$%s", io_buffer_pre, b, buffer->name);
                        } 
                        offset_this_pe_buf = 0;     /* ncpu pes wrote the buffer b, and offset_this_pe_buffer
                                                       is the offset of the part in the b written by k-th pe */
                       
                        for (kk = 0; kk < k; kk++) { 
                            offsets_array_other_pe = sizes_pe[kk];
    
    /*                      sizes_pe[kk] is sum_bb(narrays_pe_buf[kk][bb]) long,
                            and includes all bbs in this buffer */
    
                            na_offset = 0;
                            for (bb = 0; bb < b; bb++) {
                                na_offset += narrays_pe_buf[kk][bb];
                            }
                            offsets_array_other_pe += (na_offset + na_offset);
                    
                            myna = narrays_pe_buf[kk][b];
                            for (j = 0; j < myna; j++) { 
                                offset_this_pe_buf += (offsets_array_other_pe[j+j+1] * 
                                                      io_sizeof(datatypes_pe_buf[kk][b][j]));  /* offset, size */ 
                            }
                        } 
                        if (flag_for_buffered_read) {
                            io_get_array(file, tree, bufname, bio_char,
                                         offset_this_pe_buf, size_this_pe_buf, (void **)&data_pe_buf[k][b], &myid);
                        }
                    }
                    offset_this_array = (start - myoffset);
                    size_to_get = mysize - offset_this_array;
                    if (size_to_get > remaining) {
                        size_to_get = remaining;
                    }
                    offset_this_array *= io_sizeof(datatype_stord);
                    nchars = size_to_get * io_sizeof(datatype_stord);

                    for (i = 0; i < dn; i++) { 
                        offset_this_array += (offsets_array[i+i+1] * io_sizeof(datatypes[i]));  
                    }
                    if (datatype == bio_char) { 
                        t = (void *)((char *) (*array) + (size - remaining));
                    } 
                    else if (datatype == bio_int) {
                        t = (void *)((int *) (*array) + (size - remaining));
                    } 
                    else if (datatype == bio_long) {
                        t = (void *)((long *) (*array) + (size - remaining));
                    }
                    else if (datatype == bio_long_long) {
                        t = (void *)((long long *) (*array) + (size - remaining));
                    }
                    else if (datatype == bio_float) {
                        t = (void *)((float *) (*array) + (size - remaining));
                    }
                    else if (datatype == bio_double) {
                        t = (void *)((double *) (*array) + (size - remaining));
                    }
                    else if (datatype == bio_long_double) {
                        t = (void *)((long double *) (*array) + (size - remaining));
                    }  
                    if (!flag_for_buffered_read) { 
                        offset_this_array += offset_this_pe_buf; 
                        if (datatype == datatype_stord) {
                            io_get_array(file, tree, bufname, bio_char, offset_this_array, nchars, (void **)&t, &myid);
                        }
                        else {
                            s = malloc(nchars);
                            io_get_array(file,tree,bufname,bio_char,offset_this_array,nchars,(void **)&s, &myid);
                            io_array_copy(s, t, size_to_get, datatype_stord, datatype);
                            free(s);
                        }
                    }  
                    else { 
                        s = (void *)((char *)data_pe_buf[k][b] + offset_this_array);
                        
                        if (datatype == datatype_stord) { 
                            nchars = size_to_get * io_sizeof(datatype_stord);
                            memcpy(t, s, (size_t) nchars);
                        } 
                        else { 
                            io_array_copy(s, t, size_to_get, datatype_stord, datatype);
                        } 
                    }
                    start += size_to_get; 
                    remaining -= size_to_get;
                    break;
                }
                if (remaining < 1) break;
            }
        }
        else if (collective == 0) { 

/*          structure of helper: 
               header (ie, npes, collective, tot_nbufs, tot_size_in_data, tot_size_in_name, tot_narrays, datatype)
            
               part1: 
                   nbufs_pe0,        nbufs_pe1,        ...
                   narrays_pe0,      narrays_pe1,      ...
                   size_in_data_pe0, size_in_data_pe1, ,,,
                   size_in_name_pe0, size_in_name_pe1, ...
        
               pe0:    narrays_buf0, narrays_buf1, ....
               pe1:    narrays_buf0, narrays_buf1, ....
               ....
               npes-1: narrays_buf0, narrays_buf1,
        
               pe0:     datatypes of narrays,  size_array0, size_array1, .....
               pe1:     datatypes of narrays,  size_array1, size_array1, .....
               ....
               last_pe: datatypes of narrays,  size_array0, aize_array1, ...
*/  
            nbufs_pe   = helper + io_nshared_in_ibuf;
            narrays_pe = nbufs_pe + ncpu;

            start = offset;
            offset_narray = 0;
 
            found = 0;
            for (k = 0; k < ncpu; k++) {
                narrays_buf = narrays_pe_buf[k];
                nbufs = nbufs_pe[k];
                sizes = sizes_pe[k];

                offset_in_this_pe = 0;

                /* pick the buffer */

                for (b = 0; b < nbufs; b++) { 
                    na = narrays_buf[b];
                    if ((n >= offset_narray + na) || (n < offset_narray)) {
                       offset_narray += na;
                       offset_in_this_pe += na;
                       continue;
                    } 
                    /*   offset_narray <= n < offset_narray + na            */
                    /* The array requested belongs to rank = k and buffer b */
                  
                    found = 1;
                    datatypes = datatypes_pe_buf[k][b];

                    size_this_pe_buf = 0;
                    llp = sizes_pe[k] + offset_in_this_pe; 
                    for (i = 0; i < na; i++) {   
                        size_this_pe_buf += (llp[i] * io_sizeof(datatypes[i]));
                    }
                    offset_this_pe_buf = n - offset_narray;   
                    datatype_stord = datatypes[offset_this_pe_buf];
 
                    if (!data_pe_buf[k][b]) {
                       
                        /* free other allocated memory  */
          
                        for (kk = 0; kk < ncpu; kk++) { 
                            mynbufs = nbufs_pe[kk];
                            for (bb = 0; bb < mynbufs; bb++) { 
                                if (buffer->data_pe_buf[kk][bb]) {  
                                    free(buffer->data_pe_buf[kk][bb]);
                                    buffer->data_pe_buf[kk][bb] = NULL;
                                }
                            }
                        } 
                        tree = io_find_tree(file, pid);
                        if (!tree) { 
                            printf("ERROR: pid not found in bio_write\n");
                            if (iallocated) {
                                free(*array);
                                *array = NULL;
                            }
                            return -1;
                        }
                        data_pe_buf[k][b] = malloc(size_this_pe_buf);

                        myoffset = 0;
                        if (b == 0) {
                            sprintf(bufname,"%s%d$%s", io_ibuffer_pre, k, buffer->name);
                        }
                        else { 
                            sprintf(bufname,"%s%d$%d$%s", io_ibuffer_pre, k, b, buffer->name);
                        } 
                        io_get_array(file, tree, bufname, bio_char, 
                                     myoffset, size_this_pe_buf, (void **)&data_pe_buf[k][b], &myid);

                    }
                    data = (void *) data_pe_buf[k][b];

                    offset_this_array = 0;
                    sizes = sizes_pe[k];
                    for (i = 0; i < offset_this_pe_buf; i++) { 
                        offset_this_array += (sizes[i] * io_sizeof(datatypes[i]));
                    }
                    offset_this_array += (offset * io_sizeof(datatype_stord));

                    if (offset + size > sizes[offset_in_this_pe]) { 
                        printf("ERROR: offset + size > size of the array\n");
                        if (iallocated) {
                            free(*array);
                            *array = NULL;
                        }
                        return -1;
                    }
                    else {
                        size_to_get = size;
                    }
                    s = (void *)((char *) data + offset_this_array);  

                    if (datatype == datatype_stord) {
                        nchars = size_to_get * io_sizeof(datatype_stord);
                        memcpy(*array, s, (size_t) nchars);
                    }
                    else {
                        io_array_copy(s, *array, size_to_get, datatype_stord, datatype);
                    }
                    break;
                }           /* for (b = 0; b < nbuf; b++) */         
                if (found) break;
            }               /* for (k = 0; k < ncpu; k++) */
        }   
        *array_id = buffer->id; 
    } 
    else { 
        tree = io_find_tree(file, grp_id);
        if (!tree) {  
            printf("ERROR: grp_id not found in bio_write\n");
            if (iallocated) {
                free(*array);
                *array = NULL;
            }
            return -1;
        }  
        c = (char *)tree->p;
        memcpy(a6, c, (size_t)io_szllong);
        myslen = a6[0];
        myname = c + io_szllong;
        c = myname + myslen;
        memcpy(a6, c, (size_t)(6 * io_szllong));
        id_stord = (int) a6[0];
        type_stord = (bio_Object_Type) ((int)a6[2]);
    
        if (type_stord != bio_group) {
            printf("ERROR: grp_id is not for a group in bio_read\n");
            if (iallocated) {
                free(*array);
                *array = NULL;
            }
            return -1;
        }
        err = io_get_array(file, tree, name, datatype, offset, size, array, array_id);
    } 
    return err;
   }   
/****************************************************************************************/
long long bio_file_close(int fileid)
{ 
    int err, k, num_obj;
    char *pc;
    long long filesize_written, *pll;
    double dt;
    io_File   *file, *next, *previous;
    io_Obj *obj;
    io_Attr   *a;  
    void      *meta_data;
    io_Buffer *pre_buffer, *buffer, *next_buffer; 
    io_Tree *tree;

#ifdef MPI
    MPI_Status status; 
#endif 

     previous = NULL;
     file = io_files;
     while (file) { 
         if (file->id == fileid) break;
         next = file->next;
         previous = file;
         file = next;
     }
     if (!file) { 
         printf("ERROR: file not exist in bio_file_close\n");
         return -1;
     }
     pre_buffer = NULL;
     buffer = io_buffers;
     while (buffer) {
          next_buffer = buffer->next; 
          if (buffer->file == file) {
             
              if ((buffer->helper) && ((file->mode == bio_file_create)||(file->mode == bio_file_read_write))) {  
                  if (buffer->collective == 1) { 
                      err = io_buffer_finalize(file, buffer);
                  }
                  else { 
                      err = io_buffer_finalize_nonc(file, buffer);  
                  }
              } 
              io_buffer_clean(buffer);
              free(buffer); 
          
              if (pre_buffer) {
                  pre_buffer->next = next_buffer;
              }
              else { 
                  io_buffers = next_buffer;
              } 
          }
          else {
              pre_buffer = buffer;
          }
          buffer = next_buffer;
     }  
     if ((file->mode == bio_file_create) || (file->mode == bio_file_read_write))  {  
         filesize_written = io_file_flush(file);
     }
     else { 
         filesize_written = 0;
     }
#ifdef MPI 
     if ((file->io == bio_collective) || (file->io == bio_independent)) { 
         MPI_File_close(&(file->mfile));
     }
     else if (file->io == bio_use_fopen) { 
        /*** 
         if (mype) { 
            fflush( file->fp);
         }
         else { 
            fclose(file->fp);
         }
         *****/
         fclose(file->fp);
     }          
     else if (file->io == bio_use_open) { 
         close(file->fd);
     } 
     file->t_close = MPI_Wtime();
     dt = file->t_close - file->t_open;
     if (!mype) { 
//         printf("%s: %12.4e (sec) from open to close, Mflop = %12.4e\n", file->name, dt, file->mbytes/dt);
     } 
#else 
     if ((file->io == bio_use_fopen) || (file->io == bio_collective) || (file->io == bio_independent)) { 
         fclose(file->fp);
     }
     else if (file->io == bio_use_open) { 
         close(file->fd);
     } 
     time(&(file->t_close));
     dt = difftime(file->t_close, file->t_open);
#endif

/*   clean up  */  

     obj = file->object;
     num_obj = *(obj->num);

     for (k = 0; k < num_obj; k++)  { 
         a = obj->attr[k];
         if (a) {  
             free(a->value);
             free(a);
         } 
     }
     free(file->name);
     if (obj) {
         free(obj->value);
         free(obj->attr);
         free(obj);
     } 
     tree = file->root;
     io_clean_tree(tree);
     
     if (!previous)
        { io_files = file->next;
         }
     else 
        { previous->next = file->next;
         }
     free(file);
     io_num_file--; 

     if (io_num_file == 0) {  
         io_next_id = 0;

/****
         if (io_hdr) free(io_hdr);
         io_hdr = NULL;
         io_next_id = 0;
#ifdef MPI 
         if (bio_info != MPI_INFO_NULL) { 
              MPI_Info_free(&bio_info);
              bio_info = MPI_INFO_NULL;
             } 
         if (mpi_status) { 
             free(mpi_status);
             mpi_status = NULL;
         }
         if (mpi_reqs) {
            free(mpi_reqs);
            mpi_reqs = NULL;
         }
#endif 
********/

     }
     return filesize_written;
 }   
/****************************************************************************************/
long long bio_file_flush(int fileid)
{ 
    io_File   *file, *next;
    long long filesize_written;

    file = io_files;
    while (file)
       { if (file->id == fileid) break;
         next = file->next;
         file = next;
        }
    if (!file)
       { if (!mype) printf("ERROR: file not exist in bio_file_flush\n");
         return -1;
        }
    if ((file->mode == bio_file_create) || (file->mode == bio_file_read_write))  {  
        filesize_written = io_file_flush(file);
    } 
    else { 
        filesize_written = 0;
    }
    return filesize_written;
 } 
/****************************************************************************************/
long long io_file_flush(io_File *file)
{ 
    int k, num_obj, size_backup;
    char *pc, filename_backup[128];
    long long filesize_written, size_written, myoffset, llversion_num, size_meta_block;
    long long a5[5];
    long long *pll;
    bio_File_Mode mode;
    io_Obj   *obj;
    io_Attr  *a;  
    char     *meta_data, *data_backup;
    FILE     *fp;
#ifdef MPI
    MPI_Status status; 
#endif 

    mode = file->mode;
    obj  = file->object;
    num_obj = *(obj->num);

    if ((mode != bio_file_create) &&  (mode != bio_file_read_write))  { 
         printf("ERROR: bio_open_mode not consistent with io_file_flush\n");
         return -1;
    } 
/*  size_meta_block includes size_meta, and size of all attributes  of each object */ 

    size_meta_block = *(obj->size_meta);
    for (k  = 0; k < num_obj; k++) {  
        a = obj->attr[k];
        size_meta_block += *(a->size);
    }
    io_hdr[0] = 'a';
    io_hdr[1] = 'a' + io_version_num;
    io_hdr[2] = 'a' + io_szint; 
    io_hdr[3] = 'a' + io_szlong; 
    io_hdr[4] = 'a' + io_szllong; 
    io_hdr[5] = 'a' + io_szfloat; 
    io_hdr[6] = 'a' + io_szdouble; 
    io_hdr[7] = 'a' + io_szldouble; 

    llversion_num = io_version_num;
    pc = io_hdr + 8;
    memcpy(pc, &llversion_num, (size_t) io_szllong);

    pc = io_hdr + io_nchars_in_hdr;

    /*****
    pll    = (long long *)(io_hdr + io_nchars_in_hdr); 
    pll[0] = (long long) io_my_hdr_size;
    pll[1] = (long long) file->id_min;        io_id_min 
    pll[2] = (long long) io_next_id - 1;      io_id_max poosible 
    pll[3] = *(obj->size_data);  
    pll[4] = size_meta_block;
    ******/

    a5[0] = (long long) io_my_hdr_size;
    a5[1] = (long long) file->id_min;   /* io_id_min */
    a5[2] = (long long) io_next_id - 1; /* io_id_max poosible */
    a5[3] = *(obj->size_data);
    a5[4] = size_meta_block;
    memcpy(pc, a5, (size_t)(5 * io_szllong));

    myoffset = 0;
    if (mype == 0) { 
        size_written = 0;
        file->mbytes += (0.000001 * (double)io_my_hdr_size); 
#ifdef MPI 
            
        if ((file->io == bio_collective) || (file->io == bio_independent)) { 
            /***
            MPI_File_seek(file->mfile, (MPI_Offset)myoffset, MPI_SEEK_SET);
            MPI_File_write(file->mfile, io_hdr, (int)size, MPI_CHAR, &status);
            ***/
            MPI_File_write_at(file->mfile, (MPI_Offset)myoffset, io_hdr, (int)io_my_hdr_size, MPI_CHAR, &status);
        }
        else if (file->io == bio_use_fopen) { 
            rewind(file->fp);
            size_written = (long long) fwrite(io_hdr, sizeof(char), (size_t)(io_my_hdr_size), file->fp); 
            if (size_written != (long long)io_my_hdr_size) { 
                printf("ERROR: size_written != io_my_hdr_size in io_flush\n");
                return -1;
            }
        }
        else if (file->io == bio_use_open) {  
            size_written = (long long) pwrite(file->fd, io_hdr, (size_t)io_my_hdr_size, (off_t)0);
            if (size_written != (long long)io_my_hdr_size) { 
                printf("ERROR: size_written != io_my_hdr_size in io_flush\n");
                return -1;
            }
        }
#else 
        if ((file->io == bio_use_fopen) || (file->io == bio_collective) || (file->io == bio_independent)) { 
            rewind(file->fp);
            size_written = (long long) fwrite(io_hdr, sizeof(char), (size_t)(io_my_hdr_size), file->fp); 
        }
        else if (file->io == bio_use_open) {  
            size_written = (long long) pwrite(file->fd, io_hdr, (size_t)io_my_hdr_size,(off_t)0);
        }
        if (size_written != (long long)io_my_hdr_size) { 
            printf("ERROR: size_written != io_my_hdr_size in io_flush\n");
            return -1;
        }
#endif
    } 
    myoffset = (long long)io_my_hdr_size + *(obj->size_data);
    filesize_written = myoffset + size_meta_block;
    if (mype == 0) {
        meta_data = malloc((size_t)size_meta_block);
        memcpy(meta_data, obj->value, (size_t)(*(obj->size_meta)));
        pc = (char *) meta_data + *(obj->size_meta);

        for (k = 0; k < num_obj; k++) {
            a = obj->attr[k];
            memcpy(pc, a->value, (size_t)(*(a->size)));
            pc += *(a->size);
        }
    }
    else { 
        meta_data = NULL;
    }
    if (mype == 0) { 
        file->mbytes += (0.000001 * (double)size_meta_block); 
        size_written = 0;
#ifdef MPI 
        if ((file->io == bio_collective) || (file->io == bio_independent)) { 
           /***
            MPI_File_seek(file->mfile, (MPI_Offset)myoffset, MPI_SEEK_SET);
            MPI_File_write(file->mfile, meta_data, (int)size_meta_block, MPI_CHAR, &status);
            ****/
            MPI_File_write_at(file->mfile, (MPI_Offset)myoffset, meta_data, (int)size_meta_block, MPI_CHAR, &status);
        }
        else if (file->io == bio_use_fopen) { 
           fseek(file->fp, (long)(io_my_hdr_size + *(obj->size_data)), SEEK_SET);
           size_written = (long long) fwrite(meta_data, sizeof(char), (size_t) size_meta_block, file->fp);  
           if (size_written != size_meta_block) {  
               printf("ERROR: size_written != size_meta_block in io_flush\n");
               return -1;
           }  
        }
        else if (file->io == bio_use_open) { 
            size_written = (long long) pwrite(file->fd, meta_data, (size_t)size_meta_block,
                                            (off_t)(io_my_hdr_size + *(obj->size_data)));
            if (size_written != size_meta_block) {  
               printf("ERROR: size_written != size_meta_block in io_flush\n");
               return -1;
           }  
        }
#else 
        if ((file->io == bio_use_fopen) || (file->io == bio_collective) || (file->io == bio_independent)) {  
            fseek(file->fp, (long)(io_my_hdr_size + *(obj->size_data)), SEEK_SET);
            size_written = (long long) fwrite(meta_data, sizeof(char), (size_t) size_meta_block, file->fp);  
        }
        else if (file->io == bio_use_open) { 
            size_written = (long long) pwrite(file->fd, meta_data, (size_t)size_meta_block,
                                            (off_t)(io_my_hdr_size + *(obj->size_data)));
        }
        if (size_written != size_meta_block) {  
            printf("ERROR: size_written != size_meta_block in io_flush\n");
            return -1;
        }  
#endif
        if (flag_for_resilience) { 
            if (io_filename_backup) { 
                strcpy(filename_backup, io_filename_backup);
            }
            else { 
                sprintf(filename_backup, "%s.backup", file->name); 
            } 
            fp = fopen(filename_backup, "w");
            size_backup = io_my_hdr_size + size_meta_block;
            data_backup = (char *) malloc(size_backup);
            memcpy(data_backup, io_hdr, (size_t) io_my_hdr_size);
            memcpy(data_backup + io_my_hdr_size, meta_data, (size_t) size_meta_block); 
            fwrite(data_backup, sizeof(char), (size_t)size_backup, fp);
            fclose(fp); 
            free(data_backup);
        } 
        free(meta_data);
    } 
    return filesize_written;
 }   
/****************************************************************************************/
int io_get_list(int count_only, io_File *file, io_Tree *ptree, int read_value_too,
                char *filter, int *n, bio_List_Struct **list)
{ 
    int narrays, num, id, s, i, k, slen, myslen, collective;
    int id_stord, grpid_stord;
    bio_Object_Type type_stord;
    bio_Data_Type datatype_stord;    
    long long a6[6], offset, nchars, size;

    long long  *llp;

    char name[MAXLNAME]; 
    char *pc, *c, *myname, *purename;
    void *v; 
    bio_Data_Type datatype; 
    long long tot_size, myoffset; 
    io_Obj *obj;
    bio_List_Struct this_obj, *array_list, *mylist;

    io_Tree  *tree;
    io_Child *child, *nextc;

/* 
          size_meta,          szllong 
          fid,                szllong 
          num,                szllong 
          size_data,          szllong
          size of name,       szllong 
          name,               slen*szchar
          id,                 szllong 
          parent_id,          szllong  
          obj_type,           szllong 
          datatype,           szllong  
          offset,             szllong
          size,               szllong
*/ 
    *n = 0;

   /*** for buffers  */

    if (filter) { 
       slen = strlen(filter);

       /*** for buffers */

        sprintf(name, "%s%s", io_buffer_pre, filter);

        io_set_this_obj(NULL, name, ptree, file, 0, &num, &this_obj);
        if (num == 1) { 
            collective = 1;
        }
        else { 
            sprintf(name, "%s%s", io_ibuffer_pre, filter);
            io_set_this_obj(NULL, name, ptree, file, 0, &num, &this_obj);
            if (num == 1) { 
                collective = 0;
            }
        } 
        if (num == 1) { 
            (*n)++;
            if (this_obj.name) free(this_obj.name);
            
            if (count_only) return 0;

            if (!(*list)) { 
                 *list = (bio_List_Struct *) malloc((size_t)num * sizeof(bio_List_Struct));
            }
            io_set_this_buffer(file, ptree->id, filter, this_obj.id, collective, *list);

            return 0;
        }
       /** for non-buffers **/

        num = 0;
        child = ptree->child;
        while (child) { 
            tree = child->tree;
            if (!tree) { 
                printf("ERROR: null tree in io_get_list\n");
                return -1;
            }
            c = (char *)tree->p;
            memcpy(a6, c, (size_t) io_szllong);
            myslen = a6[0];
            myname = c + io_szllong;
            c = myname + myslen;
            memcpy(a6, c, (size_t)(6 * io_szllong));
            id_stord = (int) a6[0];
            type_stord = (bio_Object_Type) ((int)a6[2]);
    
            if ((slen == myslen) && !strncmp(myname,filter,(size_t)myslen)) { 
                num++;
                break;
            } 
            nextc = child->next;
            child = nextc;
        } 
        if (num == 1) { 
            (*n)++;
            if (count_only) return 0;
              
            if (!(*list)) { 
                *list = (bio_List_Struct *) malloc((size_t)num * sizeof(bio_List_Struct));
            }
            io_set_this_obj(child, NULL, NULL, file, read_value_too, &num, *list);
            return 0;
        } 
    }
    else {  
       /*** for buffer **/

        num = 0;
        child = ptree->child;
        while (child) {
            collective = -1;
            tree = child->tree;
            if (!tree) {
                printf("ERROR: null tree in io_get_list\n");
                return -1;
            }
            c = (char *)tree->p;
            memcpy(a6, c, (size_t) io_szllong);
            myslen = a6[0];
            myname = c + io_szllong;
            c = myname + myslen;
            memcpy(a6, c, (size_t)(6 * io_szllong));
            id_stord = (int) a6[0];
            type_stord = (bio_Object_Type) ((int)a6[2]);
    
            if (type_stord == bio_buffer) { 
                memcpy(name, myname, (size_t)myslen);
                name[myslen] = '\0';
                c = mystrchar(name, '$');
                if (!c) { 
                    printf("ERROR: block name for buffers in datafile\n");
                    return -1;
                }
                pc = mystrchar(c+1, '$');
                if (!pc) { 
                    if (!strncmp(name, io_buffer_pre, (size_t)io_lbuffer_pre)) { 
                        collective = 1;
                        num++;
                    }
                }
                else { 
                   c = mystrchar(pc+1, '$');
                   if (!c && !strncmp(name, io_ibuffer_pre, (size_t)io_libuffer_pre)) { 
                        collective = 0;
                        num++;
                   }
                }                          
            } 
            nextc = child->next;
            child = nextc;
        } 
       /*** for non-buffers  ***/

        child = ptree->child;
        while (child) { 
            tree = child->tree;
            if (!tree) { 
                printf("ERROR: null tree in io_get_list\n");
                return -1;
            }
            c = (char *)tree->p;
            memcpy(a6, c, (size_t) io_szllong);
            myslen = a6[0];
            myname = c + io_szllong;
            c = myname + myslen;
            memcpy(a6, c, (size_t)(6 * io_szllong));
            id_stord = (int) a6[0];
            type_stord = (bio_Object_Type) ((int)a6[2]);
  
            if (type_stord != bio_buffer) { 
                num++;
            }  
            nextc = child->next;
            child = nextc;
        } 
        (*n) += num;
        if (!num || count_only) {  
            return 0;
        }
        if (!(*list)) {  
            *list = (bio_List_Struct *) malloc((size_t)num * sizeof(bio_List_Struct));
        }
        mylist = *list;

       /** get information for buffers **/

        child = ptree->child;
        while (child) {
            collective = -1;
            tree = child->tree;
            if (!tree) {
                printf("ERROR: null tree in io_get_list\n");
                return -1;
            }
            c = (char *)tree->p;
            memcpy(a6, c, (size_t) io_szllong);
            myslen = a6[0];
            myname = c + io_szllong;
            c = myname + myslen;
            memcpy(a6, c, (size_t)(6 * io_szllong));
            id_stord = (int) a6[0];
            type_stord = (bio_Object_Type) ((int)a6[2]);
            datatype_stord = (bio_Data_Type) ((int)a6[3]);
    
            if (type_stord == bio_buffer) { 
                memcpy(name, myname, (size_t)myslen);
                name[myslen] = '\0';
                c = mystrchar(name, '$');
                if (!c) { 
                    printf("ERROR: block name for buffers in datafile\n");
                    return -1;
                }
                pc = mystrchar(c+1, '$');
                if (!pc) { 
                    if (!strncmp(name, io_buffer_pre, (size_t)io_lbuffer_pre)) { 
                        purename = c + 1;
                        collective = 1;
                        io_set_this_buffer(file, ptree->id, purename, id_stord, collective, mylist);
                        mylist++;
                    }
                }
                else { 
                   c = mystrchar(pc+1, '$');
                   if (!c && !strncmp(name, io_ibuffer_pre, (size_t)io_libuffer_pre)) { 
                        purename = pc + 1;
                        collective = 0;
                        io_set_this_buffer(file, ptree->id, purename, id_stord, collective, mylist);
                        mylist++;
                   }
                }             
            } 
            nextc = child->next;
            child = nextc;
        } 
       /** get information for non-buffers */

        child = ptree->child;
        while (child) {
    
            tree = child->tree;
            c = (char *)tree->p;
            memcpy(a6, c, (size_t) io_szllong);
            myslen = a6[0];
            myname = c + io_szllong;
            c = myname + myslen;
            memcpy(a6, c, (size_t)(6 * io_szllong));
            id_stord = (int) a6[0];
            grpid_stord = (int) a6[1];
            type_stord  = (bio_Object_Type) ((int)a6[2]);
    
            if (type_stord != bio_buffer) { 
    
                io_set_this_obj(child, NULL, NULL, file, read_value_too, &num, mylist);
                mylist++;
            }                    
            nextc = child->next;
            child = nextc;
        }    
    }  

    return 0;
  } 

int io_get_list_from_buffer(io_File *file, int count_only, io_Buffer *buffer, int read_value_too, 
                          char *filter, int *n, bio_List_Struct **list)
{
    char purename[MAXLNAME];
    char *names, *name;
    int  i, i1, k, mytype, narrays, collective, slen, ncpu, cpu, id;
    long long       **sizes_pe, *sizes, size;
    bio_List_Struct *obj;
    bio_Data_Type   *datatypes;
    
    *n = 0;
    if (!buffer) return 0;

    narrays    = buffer->narrays;
    names      = buffer->names;
    collective = buffer->collective;
    datatypes  = buffer->datatypes;
    sizes_pe   = buffer->sizes_pe;
    ncpu       = buffer->npes; 
    
    if (filter) {
        name = names;
        i = 0;
        while (i < narrays) { 
            slen = strlen(name);
            if (collective == 1) {
                if (!strcmp(filter, name)) {
                    break;
                }
            }
            else if (collective == 0) {
                io_buffer_getname(name, purename, &cpu);
                if (!strcmp(purename, filter)) {
                    break;
                }
            }
            else {
                printf("ERROR: collective in io_get_list_from_buffer\n");
                return -1;
            }
            name += (slen + 1 + io_szint); /* there is an int after each name,
                                               which is obj_type of the array */
            i++;
        }
        if (i >= narrays) return 0;

        *n = 1;
        if (count_only) return 0;
 

        if (!(*list)) {
            *list = (bio_List_Struct *) malloc(*n * sizeof(bio_List_Struct));
        }
        io_list_null(*n, *list);
                
        slen = strlen(name) + 1;
        obj = *list;
        obj->name = (char *) malloc(slen);
        strcpy(obj->name, name);
//      obj->id = buffer->id;
        obj->id = bio_object_type_undefined;
        memcpy(&mytype, name + (slen + 1), (size_t)io_szint);
        obj->type = (bio_Object_Type) mytype;
        obj->datatype = datatypes[i];
        if (collective == 1) { 
            i1 = i + i + 1;
            size = 0;
            for (k = 0; k < ncpu; k++) { 
                sizes = sizes_pe[k];
                size += sizes[i1];
            }
            obj->size = size;
        }
        else { 
            sizes = sizes_pe[cpu];
            obj->size = sizes[i];
        }
        if (read_value_too) { 
            io_read(file, buffer,
            buffer->id, filter, datatypes[i], (long long)0, obj->size, obj->value, &id);
        }
    }
    else { 
        *n = narrays;
        if (count_only) return 0;

        if (!(*list)) {
            *list = (bio_List_Struct *) malloc(*n * sizeof(bio_List_Struct));
        }
        io_list_null(*n, *list);
        
        name = names;
        for (i = 0; i < *n; i++) {
            obj = *list + i;
            slen = strlen(name);
            obj->name = (char *) malloc(slen + 1);
            strcpy(obj->name, name);
//          obj->id = buffer->id;
            obj->id = bio_object_type_undefined;
            memcpy(&mytype, name + (slen + 1), (size_t)io_szint);
            obj->type = (bio_Object_Type) mytype;
            obj->datatype = datatypes[i];
            if (collective == 1) {
                i1 = i + i + 1;
                size = 0;
                for (k = 0; k < ncpu; k++) {
                    sizes = sizes_pe[k];
                    size += sizes[i1];
                }
                obj->size = size;
            }
            else { 
                obj->size = sizes_pe[cpu][i];
            }  
            if (read_value_too) { 
                io_read(file, buffer,
                buffer->id, obj->name, datatypes[i], (long long)0, obj->size, obj->value, &id);
            }
            name += (slen + 1 + io_szint); /* there is an int after each name,
                                               which is obj_type of the array */
        }
    } 
    return 0;
 } 
/****************************************************************************************/
int io_set_this_obj(io_Child *this_child, char *fullname, io_Tree *ptree,
                    io_File *file, int read_value_too,
                    int *num, bio_List_Struct *obj)
{  
    /* If this_child is not null, take its members, and ptree and fullname may be null.
       If this_child is null, take fullname and ptree, and search fullname.

       num :  output
              = 0 : nothing is set
              = 1 : obj has been set 
     */  

    int slen, myslen;
    int sztype, id_stord, grpid_stord;
    bio_Object_Type type_stord;
    bio_Data_Type datatype_stord;
    long long a6[6], offset, nchars, size, myoffset, size_read;
    
    long long *llp;
    
    char *pc, *c, *myname;

    io_Tree  *tree;
    io_Child *child, *nextc;

#ifdef MPI
    MPI_Datatype mpidatatype;
    MPI_Status status;
#endif

/*
          size_meta,          szllong
          fid,                szllong
          num,                szllong
          size_data,          szllong
          size of name,       szllong
          name,               slen*szchar
          id,                 szllong
          parent_id,          szllong
          obj_type,           szllong
          datatype,           szllong
          offset,             szllong
          size,               szllong
*/
    child = NULL;
    if (this_child) { 
        child = this_child;
 
        tree = child->tree;
        c = (char *)tree->p;
        memcpy(a6, c, (size_t) io_szllong);
        myslen = a6[0];
        myname = c + io_szllong;
        c = myname + myslen;
        memcpy(a6, c, (size_t)(6 * io_szllong));
        id_stord = (int) a6[0];
        grpid_stord = (int) a6[1];
        type_stord  = (bio_Object_Type) ((int)a6[2]);
        datatype_stord = (bio_Data_Type) ((int)a6[3]);
        offset = a6[4];
        nchars = a6[5];
    } 
    else { 

        if (!fullname) { 
            printf("ERROR: null fullname in io_get_this_obj\n");
            return -1;
        }  
        slen = strlen(fullname);

        child = ptree->child;
        while (child) {
    
            tree = child->tree;
            c = (char *)tree->p;
            memcpy(a6, c, (size_t) io_szllong);
            myslen = a6[0];
            myname = c + io_szllong;
            c = myname + myslen;
            memcpy(a6, c, (size_t)(6 * io_szllong));
            id_stord = (int) a6[0];
            grpid_stord = (int) a6[1];
            type_stord  = (bio_Object_Type) ((int)a6[2]);
            datatype_stord = (bio_Data_Type) ((int)a6[3]);
            offset = a6[4];
            nchars = a6[5];
            if ((slen == myslen) && !strncmp(myname,fullname,(size_t)myslen)) {
                break;
            }
            nextc = child->next;
            child = nextc;
        }
    } 
    if (!child) { 
        *num = 0;
        return 0;
    }
    *num = 1;
    obj->name = (char *)malloc(myslen + 1);
    memcpy(obj->name, myname, (size_t)myslen);
    obj->name[myslen] = '\0';
    obj->id   = id_stord;
    obj->type = type_stord;
    obj->datatype = datatype_stord;
    sztype = io_sizeof(datatype_stord);
    size   = nchars / sztype;
    obj->size = size;

    obj->value = NULL;
    obj->narrays = 0;
    obj->names = NULL;
    obj->sizes = NULL;
    obj->datatypes = NULL;

    if (read_value_too && (type_stord != bio_group)) {

/*      This is to read an actual dataset from a disk. */
  
        obj->value = malloc((size_t)(nchars + 1));
        myoffset = (long long)io_my_hdr_size + offset;
        size_read = 0;
#ifdef MPI
        if ((file->io == bio_collective) || (file->io == bio_independent)) {
            if (datatype_stord == bio_char) {
                mpidatatype = MPI_CHAR;
            }
            else if (datatype_stord == bio_int) {
                mpidatatype = MPI_INT; 
            }
            else if (datatype_stord == bio_long) {
                mpidatatype = MPI_LONG;
            }
            else if (datatype_stord == bio_long_long) {
                mpidatatype = MPI_LONG_LONG_INT;
            }
            else if (datatype_stord == bio_float) {
                mpidatatype = MPI_FLOAT;
            }
            else if (datatype_stord == bio_double) {
                mpidatatype = MPI_DOUBLE;
            }
            else if (datatype_stord == bio_long_double) {
                mpidatatype = MPI_LONG_DOUBLE;
            }
            else {
                printf("ERROR: datatype_stord incorrect in io_get_this_obj\n");
                return -1;
            }
            MPI_File_read_at(file->mfile, (MPI_Offset)myoffset, obj->value, (int)size,
                             mpidatatype, &status);
        }
        else if (file->io == bio_use_fopen) {
            fseek(file->fp, (long)myoffset, SEEK_SET);
            size_read = (long long) fread(obj->value, (size_t)sztype, (size_t)size, file->fp);
            if (size_read != size) {
                printf("ERROR: szie_read != size in io_get_this_obj\n");
                return -1;
            }
        }
        else if (file->io == bio_use_open) {
            size_read = (long long) pread(file->fd, obj->value, (size_t)nchars, (off_t) myoffset);
            if (size_read != nchars) {
                printf("ERROR: szie_read != nchars in io_get_this_obj\n");
                return -1;
            }
        }
#else
        if ((file->io == bio_use_fopen) || (file->io == bio_collective) || (file->io == bio_independent)) {
            fseek(file->fp, (long)myoffset, SEEK_SET);
            size_read = (long long) fread(obj->value, (size_t)sztype, (size_t)size, file->fp);   
            if (size_read != size) {
                printf("ERROR: szie_read != size in io_get_list\n");
                return -1;
            }
        }
        else if (file->io == bio_use_open) {
            lseek(file->fd, (off_t) myoffset, 0);
            size_read = (long long) pread(file->fd, obj->value, (size_t)nchars, (off_t) myoffset);
            if (size_read != nchars) {
                printf("ERROR: szie_read != nchars in io_get_list\n");
                return -1;
            }
        }
#endif
        if (io_2swap && datatype_stord != bio_char) {
            io_swap(obj->value, (int)sztype, (long long) size);
        }
    } 

    return 0;
} 
/****************************************************************************************/

int io_set_this_buffer(io_File *file, int parentid, char *purename, int id, 
                           int collective, bio_List_Struct *list)
{
   /* id is supposed to be the id of the first block of the buffer  */

    int nattr, num, i, k, b, obj_idx, slen, coll;
    int ncpu, narrays, nbuf, na;
    bio_Data_Type type;
    long long tot_nbufs, tot_size_in_data, tot_size_in_name, tot_narrays;
    long long tot_size, offset, myoffset;

    long long *helper, *size_in_name_pe, *sizes,  *offsets_array, *llp;
    long long *nbufs_pe, *narrays_pe, *lldatatypes;
    long long **narrays_pe_buf, **sizes_pe;
    char ***data_pe_buf;
    bio_Data_Type *datatypes, ***datatypes_pe_buf;

    char helpername[MAXLNAME], name[MAXLNAME];
    char *pc, *c, *namelist, *myname;
    char **names;
    bio_List_Struct *attr_list;
    io_Buffer *buffer, *next;
    io_Tree *tree;
    io_Obj *obj;
    io_Attr *a; 

    buffer = io_buffers;
    while (buffer) {
        if (buffer->id == id) {
            break;  /* size id of helper was considered as id of the buffer */
        } 
        next = buffer->next;
        buffer = next;
    }
    if (buffer) { 
        if (!buffer->helper_for_read) {   
            printf("ERROR: null beffer->helper_for_read in io_set_this_buffer\n");
            return -1;
        }
        if (!buffer->names)  {  
             printf("ERROR: null buffer->names in io_set_this_buffer\n");
             return -1;
        }   
        if (buffer->tot_size_in_data == 0) {  
            printf("ERROR: buffer->tot_size = 0 in io_set_this_buffer\n");
            return -1;
        }
        helper   = buffer->helper_for_read; 
        namelist = buffer->names;
        narrays_pe_buf = buffer->narrays_pe_buf;
        sizes_pe  = buffer->sizes_pe; 
        data_pe_buf = buffer->data_pe_buf;
        datatypes_pe_buf = buffer->datatypes_pe_buf;
    }  
    else { 
        if (collective) { 
            sprintf(helpername, "%s%s", io_bufhelper_pre, purename);
            sprintf(name,       "%s%s", io_bufnames_pre,  purename);
        }
        else {  
            sprintf(helpername, "%s%s", io_ibufhelper_pre, purename);
            sprintf(name,       "%s%s", io_ibufnames_pre,  purename);
        }
        if (file->id != parentid) { 
            tree = io_find_tree(file, parentid);
        }
        else {
            tree = file->root;
        }
        obj_idx = tree->idx;
        a = (file->object)->attr[obj_idx];

        nattr = 0;
        attr_list = 0;
        io_get_attr_list(a, helpername, &nattr, &attr_list);
//      bio_attr_list(file->id, parentid, helpername, &nattr, &attr_list); 
   
        if (nattr != 1) { 
            printf("ERROR: nattr = %d in io_set_this_buffer\n", nattr);
        } 
        helper = (long long *) attr_list->value;
        free(attr_list->name);
        free(attr_list);
   
        nattr = 0;
        attr_list = NULL;
        io_get_attr_list(a, name, &nattr, &attr_list);
//      bio_attr_list(file->id, parentid, name, &nattr, &attr_list); 
        if (nattr != 1) { 
            printf("ERROR: nattr = %d in io_set_this_buffer\n", nattr);
        } 
        namelist = (char *)attr_list->value;
        free(attr_list->name);
        free(attr_list);
    }   
    io_list_null(1, list);

    if (collective) { 

    /*  structure of helper:
           header                                    \ 
           narrays numbers for datatypes              \ 
           npes numbers for size_in_name_each_pe       \
           pe0: narrays_buf0, narrays_buf1, ....        \  m as follows
           pe1: narrays_buf0, narrays_buf1, ....        /
           ....                                        /
           lastpe: narrays_buf0, narrays_buf1, ...    /
    
           pe0:  (offset, size) of array0, (offset, size) of array1, ....
           pe1:  (offset, size) of array0, (offset, size) of array1, ....
           ... 
           pelast: (offset, size) of array0, (offset, size) of array1, ....
    */
        ncpu = (int) helper[0];
        coll = (int) helper[1];
        narrays  = (int) helper[2];
        nbuf     = (int) helper[3];
        tot_size = helper[4];
                   
        lldatatypes = helper + io_nshared_in_buf;
        size_in_name_pe = helper + (io_nshared_in_buf + narrays);
   
        names  = (char **)malloc((size_t)narrays * sizeof(char*));
        myname = namelist;
        na = 0;
        while (na < narrays) { 
            slen = strlen(myname);
            names[na] = (char *) malloc((slen + 1));
            strcpy(names[na], myname);
            myname += (slen + 1 + io_szint); /* there is an int after each nname, which is (bio_Object_Type) obj_type */
            na++;
        }
        if (!buffer) { 
            io_get_buffer_info(helper, &narrays_pe_buf, &sizes_pe, &data_pe_buf, &datatypes_pe_buf);
        } 
        sizes = (long long *) malloc((size_t)(narrays * io_szllong));
        for (i = 0; i < narrays; i++) {
             sizes[i] = 0;
        }
        for (k = 0; k < ncpu; k++) {
            offsets_array = sizes_pe[k];
            llp = sizes;
            for (b = 0; b < nbuf; b++) {
                na = narrays_pe_buf[k][b];
                for (i = 0; i < na; i++) {
                    llp[i] += offsets_array[i+i+1];
                }
                llp += na;
             } 
        }
        slen = strlen(purename);
        list->name = (char *) malloc((size_t)(slen + 1));
        strcpy(list->name, purename);
        list->id = id;  

        list->type = bio_buffer;
        list->value = NULL;

        list->size    = tot_size;    /* (poffset[1] * io_szchar)/io_sizeof(*pd); */
        list->sizes   = sizes;
        list->narrays = narrays;
        list->names   = names;
        datatypes = (bio_Data_Type *)malloc(narrays * sizeof(bio_Data_Type));
        memcpy(datatypes, datatypes_pe_buf[0][0], (size_t)(narrays * sizeof(bio_Data_Type)));
        list->datatypes = datatypes;

        if (!buffer) {  
            buffer = (io_Buffer *) malloc(sizeof(io_Buffer));
            io_buffer_null(buffer);
   
            buffer->name = (char *)malloc((size_t)(slen + 1));
            strcpy(buffer->name, purename);
                       
            buffer->id = id;
            buffer->grp_id = parentid;

            buffer->file = file;
            file->buffer = buffer;

            buffer->collective = collective;
            buffer->narrays = narrays; 
            buffer->narrays_buf = NULL;
                       
            buffer->npes = ncpu;
            buffer->nbufs_written = nbuf;
            buffer->names  = namelist;
            buffer->tot_size_in_data = tot_size;
            buffer->data  = NULL;
                       
            buffer->helper = NULL; 
            buffer->helper_for_read = helper;
   
            buffer->sizes_pe = sizes_pe;          /* offset, size */      
            buffer->narrays_pe_buf = narrays_pe_buf;
            buffer->datatypes = datatypes_pe_buf[0][0]; 
            buffer->datatypes_pe_buf = datatypes_pe_buf;
            buffer->data_pe_buf = data_pe_buf; 
              
            buffer->next = io_buffers;
            io_buffers = buffer;
         } 
    }                                /*  collective      */   
    else  {                          /*  collective == 0 */
        
         /*  structure of helper:
             header (ie, npes, collective, tot_nbufs, tot_size_in_data,
                     tot_size_in_name, tot_narrays)  

         part1:
             nbufs_pe0,        nbufs_pe1,        ...
             narrays_pe0,      narrays_pe1,      ...
             size_in_data_pe0, size_in_data_pe1, ,,,
             size_in_name_pe0, size_in_name_pe1, ...

             pe0:    narrays_buf0, narrays_buf1, ....
             pe1:    narrays_buf0, narrays_buf1, ....
             ....
             npes-1: narrays_buf0, narrays_buf1,

             pe0:     datatypes of narrays,   size_array0, size_array1, .....
             pe1:     datatypes of narrays,   size_array1, size_array1, .....
             ....
             last_pe: datatypes of narrays,   size_array0, aize_array1, .....
         */
         ncpu = helper[0];
         tot_nbufs = helper[2]; 
         tot_size_in_data = helper[3];
         tot_size_in_name = helper[4];
         tot_narrays      = helper[5];

         nbufs_pe = helper + io_nshared_in_ibuf;
         narrays_pe = nbufs_pe + ncpu;

         if (!buffer) { 
             /* datatypes_pe_buf[0][0] is the datatypes  for all arrays */

             io_get_buffer_info_nonc(helper, &narrays_pe_buf, &sizes_pe, &data_pe_buf, &datatypes_pe_buf);
         } 
         names = (char **)malloc(tot_narrays * sizeof(char*));
         myname = namelist;
         for (i = 0; i < tot_narrays; i++) { 
             slen = strlen(myname);
             names[i] = (char *) malloc((slen + 1));
             c = names[i];
             strcpy(c, myname);
             for (k = slen - 1; k >= 0; k--) { 
                 if (c[k] == '$') { 
                     c[k] = '\0';
                     break;
                 }
             } 
             myname += (slen + 1 + io_szint); /* there is an int after each nname, 
                                                 which is (bio_Object_Type) obj_type. */ 
         } 
         sizes = (long long *) malloc((size_t)(tot_narrays * io_szllong));
         datatypes = (bio_Data_Type *) malloc(tot_narrays * sizeof(bio_Data_Type)); 
         list->datatypes = datatypes;  
         memcpy(datatypes, datatypes_pe_buf[0][0], (size_t)(tot_narrays * sizeof(bio_Data_Type)));

         offset = 0;
         for (k = 0; k < ncpu; k++) {  /* sizes_pe[0], sizes_pe[1], .. not continguous im memory */ 
             narrays = narrays_pe[k];
             memcpy(sizes + offset,sizes_pe[k], (size_t)(narrays * io_szllong)); 
             offset += narrays;
         }       
         slen = strlen(purename);
         list->name = (char *) malloc((size_t)(slen + 1));
         strcpy(list->name, purename);

         list->id       = id;          /* the id stored for the first block of the buffer */
         list->type     = bio_buffer;  /* (bio_Object_Type)(*po); */
         
         list->sizes   = sizes;
         list->value   = NULL;
         list->narrays = tot_narrays;
         list->names   = names;

         if (!buffer) { 
             buffer = (io_Buffer *) malloc(sizeof(io_Buffer));
             io_buffer_null(buffer);

             buffer->name = (char *)malloc((size_t)(slen + 1));
             strcpy(buffer->name, purename);

             buffer->id = id;
             buffer->grp_id = parentid;

             buffer->file = file;
             file->buffer = buffer;

             buffer->collective = 0;
             buffer->narrays = tot_narrays;
             buffer->datatypes = datatypes_pe_buf[0][0];
             buffer->datatypes_pe_buf = datatypes_pe_buf;

             buffer->narrays_buf = NULL;

             buffer->npes = ncpu;
             buffer->nbufs_written = tot_nbufs;
             buffer->names  = namelist;
             buffer->tot_size_in_data = tot_size_in_data;
             buffer->data  = NULL;

             buffer->helper = NULL;
             buffer->helper_for_read = helper;

             buffer->sizes_pe = sizes_pe;;
             buffer->narrays_pe_buf = narrays_pe_buf;
             buffer->data_pe_buf = data_pe_buf;

             buffer->next = io_buffers;
             io_buffers = buffer;
         }                         /* !buffer    */  
    }

    return 0;
}
/****************************************************************************************/
int io_get_attr_list(io_Attr *a, char *filter, int *num, bio_List_Struct **list)
{          
    int attr_num, m, i, n, d, slen, slen_stord;
    int szdatatype_stord;
    char *pc, *c, *myname;
    long long size, a2[2];
    bio_Data_Type datatype_stord; 
    bio_List_Struct *mylist;

    *num = 0;
    if (!a) return 0;

    attr_num = *(a->num);
    
    d = 3 * io_szllong;
    pc = (char *) a->value + (3 * io_szllong);

    if (filter)
       { slen = strlen(filter);
         m = 0;
         for (i = 0; i < attr_num; i++)
           { 
             /****
             pll = (long long *) pc;
             slen = *pll;
             c = (char *)(pll + 1);
             pd = (long long *) (c + slen);
             mydatatype = (bio_Data_Type) ((int)*pd);
             pll = pd + 1;
             n = (int) *pll;
             if (!strcmp(c, filter))
                { m++;
                 }
             *****/

             memcpy(a2, pc, (size_t)io_szllong);
             slen_stord = (int) a2[0];
             myname = pc + io_szllong;

             c = myname + slen_stord;
             memcpy(a2, c, (size_t) (io_szllong + io_szllong));
             datatype_stord = (bio_Data_Type) ((int)a2[0]);
             n = (int) a2[1];
             if (!strncmp(myname, filter, (size_t)slen) && (slen == slen_stord) )
                { m++;
                 }
             pc += (slen_stord + n + d);
            }
        }
    else 
       { m = attr_num;
        }  
    if (!m) 
       { *num  = 0;
         return 0;
        }
    *num = m;
    if (!(*list))
       { mylist = (bio_List_Struct *) malloc((size_t)m * sizeof(bio_List_Struct));
         *list = mylist;
        }
    else 
       { mylist = *list; 
        }

    pc = (char *) a->value + (3 * io_szllong);
    if (filter)
       { m = 0; 
         for (i = 0; i < attr_num; i++)
           { /****
             slen = *pll;
             pc = (char *)(pll + 1);
             pd = (long long *) (pc + slen);
             pll = pd + 1;
             mydatatype = (bio_Data_Type) ((int)*pd);
             n = (int) *pll;
             if (!strcmp(pc, filter))
                { mylist[m].name = (char *) malloc((size_t)(slen + 1));
                  memcpy(mylist[m].name, pc, (size_t)slen);
                  mylist[m].name[slen] = '\0';
                  mylist[m].id = -1;
                  mylist[m].type = bio_object_type_undefined;
                  mylist[m].datatype = mydatatype;
                  mylist[m].size     = n/io_sizeof(mydatatype);
                  mylist[m].value = malloc((size_t)n);
                  memcpy(mylist[m].value, pll+1, (size_t)n);
                  
                  m++;
                 }
             pc += (slen + n + io_szllong + io_szllong);
             pll = (long long *)pc;
             *****/

             memcpy(a2, pc, (size_t)io_szllong);
             slen_stord = (int) a2[0];
             myname = pc + io_szllong;
             c = myname + slen_stord;
             memcpy(a2, c, (size_t) (io_szllong + io_szllong));
             c += (io_szllong + io_szllong);
             datatype_stord = (bio_Data_Type) ((int)a2[0]);
             n = (int) a2[1];
             if (!strncmp(myname, filter, (size_t)slen) && slen == slen_stord )
                { mylist[m].names = NULL;
                  mylist[m].sizes = NULL;
                  mylist[m].datatypes = NULL;
                  mylist[m].name = (char *) malloc((size_t)(slen_stord + 1));
                  memcpy(mylist[m].name, myname, (size_t)slen_stord);
                  mylist[m].name[slen_stord] = '\0';
                  mylist[m].id = -1;
                  mylist[m].type = bio_object_type_undefined;
                  mylist[m].datatype = datatype_stord;
                  szdatatype_stord = io_sizeof(datatype_stord);
                  size = n /io_sizeof(datatype_stord);
                  mylist[m].size  = size;
                  mylist[m].value = malloc((size_t)n);
                  memcpy(mylist[m].value, c, (size_t)n);
                  m++;
                 }
             pc += (slen_stord + n + d);
            }
        }
    else 
       { for (i = 0; i < attr_num; i++)
           { 
             memcpy(a2, pc, (size_t)io_szllong);
             slen_stord = (int) a2[0];
             myname = pc + io_szllong;
             c = myname + slen_stord;
             memcpy(a2, c, (size_t) (io_szllong + io_szllong));
             c += (io_szllong + io_szllong);
             datatype_stord = (bio_Data_Type) ((int)a2[0]);
             n = (int) a2[1];

             mylist[i].names = NULL;
             mylist[i].sizes = NULL;
             mylist[i].datatypes = NULL;
             mylist[i].name = (char *) malloc((size_t)(slen_stord + 1));
             memcpy(mylist[i].name, myname, (size_t)slen_stord);
             mylist[i].name[slen_stord] = '\0';
             mylist[i].id = -1;
             mylist[i].type = bio_object_type_undefined;
             mylist[i].datatype = datatype_stord;
             szdatatype_stord = io_sizeof(datatype_stord);
             size = n /io_sizeof(datatype_stord);
             mylist[i].size  = size;
             mylist[i].value = malloc((size_t)n);
             memcpy(mylist[i].value, c, (size_t)n);
             pc += (slen_stord + n + d);
            }
        }
    return 0;
   }
/****************************************************************************************/
int bio_count(int fileid, int grp_id, char *filter, int *n)
{
    char *c, *myname;
    int slen, myslen, obj_idx, id_stord, parentid;
    bio_Object_Type type;
    bio_Data_Type dtype;
    long long a6[6], offset, gsize;
    io_Obj *obj;
    io_File *file;
    io_Tree *tree;
    io_Child *child, *next;
 
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != grp_id) { 
        tree = io_find_tree(file, grp_id);
    }
    if (!tree || !file) {
        printf("ERROR: grp_id not found in bio_count\n");
        return -1;
    }
    c = (char *)tree->p;
    memcpy(a6, c, (size_t)io_szllong);
    myslen = a6[0];
    myname = c + io_szllong;
    c = myname + myslen;
    memcpy(a6, c, (size_t)(6 * io_szllong));
    id_stord = (int) a6[0];
    type = (bio_Object_Type) ((int)a6[2]);

    if (type != bio_group) {
        printf("ERROR: grp_id is not for a group in bio_count\n");
        return -1;
    }
    io_get_list(1, file, tree, 0, filter, n, NULL);

    return 0;
   }
/****************************************************************************************/
int bio_list(int fileid,
             int grp_id, char *filter, int read_value_too, int *n, bio_List_Struct **list)
{ 
    char *c, *myname;
    int slen, myslen, obj_idx, id_stord, parentid;
    bio_Object_Type type;
    bio_Data_Type dtype;
    long long a6[6], offset, gsize; 
    io_Obj *obj;
    io_File *file;
    io_Tree *tree;
    io_Buffer *buffer, *next;
    
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != grp_id) { 
        tree = io_find_tree(file, grp_id);
    }  
    if (!tree || !file) {
        printf("ERROR: grp_id not found in bio_list\n");
        return -1;
    }
    c = (char *)tree->p;
    memcpy(a6, c, (size_t)io_szllong);
    myslen = a6[0];
    myname = c + io_szllong;
    c = myname + myslen;
    memcpy(a6, c, (size_t)(6 * io_szllong));
    id_stord = (int) a6[0];
    type = (bio_Object_Type) ((int)a6[2]);

    if (type != bio_group) {

//      buffer = io_buffers;      This is more safe 
        buffer = file->buffer;
        while (buffer) {
            if (buffer->id == grp_id) {
                break;
            }
            next   = buffer->next;
            buffer = next;
        }
        if (!buffer) { 
            // printf("Warning: grp_id is not for a group in bio_list\n");
            return -1;
        }
        else { 
            io_get_list_from_buffer(file, 0, buffer, read_value_too, filter, n, list);
        }
    }
    else { 
        io_get_list(0, file, tree, read_value_too, filter, n, list);
    }
    return 0;
   } 
/****************************************************************************************/
int bio_attr_list(int fileid, int id, char *filter, int *n, bio_List_Struct **list)
{ 
    int obj_idx, parentid;
    bio_Object_Type type;
    bio_Data_Type dtype;
    long long offset, gsize; 
    io_Obj *obj;
    io_File *file;
    io_Attr *a;
    io_Tree *tree;

    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != id) {
        tree = io_find_tree(file, id);
    }
    if (!tree || !file) { 
        if (!mype) printf("ERROR: id not found in bio_attr_list\n");
        return -1; 
    }
    obj_idx = tree->idx;

    a = (file->object)->attr[obj_idx];
    io_get_attr_list(a, filter, n, list);
   
    return 0;
   }
/****************************************************************************************/
int io_swap(void *buffer, int unitsize, long long size) 
{ 
    long long i;
    int j;
    char *s, t[64];

/*****
    int is_little_endian;
    int n = 1;
    char carray[sizeof(int)];

    memcpy(carray, &n, sizeof(int));
    if (carray[0] == 0x01) {
        is_little_endian = 1;
    }
    else {
        is_little_endian = 0;
    }

*****/

    assert(unitsize <= 64);
    s = (char *) buffer;


    
/***
    if (unitsize == 8) { 
        for (i = 0; i < size; i++) { 
            t[0] = s[7];
            t[1] = s[6];
            t[2] = s[5];
            t[3] = s[4];
            t[4] = s[3];
            t[5] = s[2];
            t[6] = s[1];
            t[7] = s[0];
            
            s[0] = t[0];
            s[1] = t[1];
            s[2] = t[2];
            s[3] = t[3];
            s[4] = t[4];
            s[5] = t[5];
            s[6] = t[6];
            s[7] = t[7];
            s += 8;
        }
    }
    else if (unitsize == 4) { 
        for (i = 0; i < size; i++) { 
            t[0] = s[3];
            t[1] = s[2];
            t[2] = s[1];
            t[3] = s[0];
                   
            s[0] = t[0];
            s[1] = t[1];
            s[2] = t[2];
            s[3] = t[3];
            s += 4;
        }
    }
    else if (unitsize == 2) {
        for (i = 0; i < size; i++) {
            t[0] = s[1];
            t[1] = s[0];

            s[0] = t[0];
            s[1] = t[1];
            s += 2;
        }
    }
    else if (unitsize > 0) { 
*****/

        for (i = 0; i < size; i++) { 
            for (j = 0; j < unitsize; j++) { 
                t[j] = s[unitsize - 1 - j];
            }
            for (j = 0; j < unitsize; j++) { 
                s[j] = t[j];
            }
            s += unitsize;
        }
    return 0;
  } 
/****************************************************************************************/
int io_swap_metadata(void *metadata)
{ 
   int ii, i, j, k, num_obj, slen, num_attr;
   char *s, *ptr, t[32];
   long long size_meta_stord, size_data, size_value, size;  
   long long a6ll[6];
   long a6l[6];
   int  a6i[6];
   bio_Data_Type datatype_stord;
/*
    size_meta,          szllong
    fid,                szllong
    num,                szllong
    size_data,          szllong
    size of name,       szllong
    name,               slen*szchar
    id,                 szllong
    parent_id,          szllong
    obj_type,           szllong
    datatype,           szllong
    offset,             szllong
    size,               szllong
*/
   s = (char *) metadata;
   for (i = 0; i < 4; i++) {
       for (j = 0; j < io_szllong_stord; j++) { 
           t[j] = s[io_szllong_stord - 1 - j];
       }
       for (j = 0; j < io_szllong_stord; j++) { 
           s[j] = t[j];
       }
       s += io_szllong_stord;
   }
   if (io_szllong_stord == io_szllong) { 
       memcpy(a6ll, metadata, (size_t)(4 * io_szllong_stord));
       size_meta_stord = a6ll[0];
    /* fileid    = a6ll[1];  */
       num_obj   = a6ll[2];
       size_data = a6ll[3];
   }
   else if (io_szllong_stord == io_szint) { 
       memcpy(a6i, metadata, (size_t)(4 * io_szllong_stord));
       size_meta_stord = a6i[0];
   /*  fileid    = a6i[1];   */  
       num_obj   = a6i[2];
       size_data = a6i[3];
   }
   else if (io_szllong_stord == io_szlong) {  
       memcpy(a6l, metadata, (size_t)(4 * io_szllong_stord));
       size_meta_stord = a6l[0];
   /*  fileid    = a6l[1];   */
       num_obj   = a6l[2];
       size_data = a6l[3];
   }
   else { 
       printf("ERROR: io_szllong_stord is out the scope in io_swap_metadata\n");
       return -1;
   } 
   s = (char *) metadata + (4 * io_szllong_stord);
   for (k = 0; k < num_obj; k++) { 
       for (j = 0; j < io_szllong_stord; j++) { 
           t[j] = s[io_szllong_stord - 1 - j];
       }   
       for (j = 0; j < io_szllong_stord; j++) { 
           s[j] = t[j];
       }
       if (io_szllong_stord == io_szllong) {
           memcpy(a6ll, s, (size_t)io_szllong_stord);
           slen = a6ll[0];
       }
       else if (io_szllong_stord == io_szint) {
           memcpy(a6i, s, (size_t)io_szllong_stord);
           slen = a6i[0];
       }
       else if (io_szllong_stord == io_szlong) {
           memcpy(a6l, s, (size_t)io_szllong_stord);
           slen = a6l[0];
       }
       s += (io_szllong_stord + slen);
       
       for (i = 0; i < 6; i++) { 
           for (j = 0; j < io_szllong_stord; j++) {
               t[j] = s[io_szllong_stord - 1 - j];
           }  
           for (j = 0; j < io_szllong_stord; j++) {
               s[j] = t[j];
           }
          s += io_szllong_stord;
       }
   }
/****
    for attrs

    size,            szllong
    id,              szllong
    num,             szllong
    size of name,    szllong
    name,            slen*szchar
    datatype,        szllong
    size of value,   szllong
    values,          n*io_sizeof(datatype)
***/
 
   s = (char *) metadata + size_meta_stord;
   for (ii = 0; ii < num_obj; ii++) {

       ptr = s;                    /* size, id, num_attr */ 
       for (i = 0; i < 3; i++) { 
           for (j = 0; j < io_szllong_stord; j++) {
               t[j] = s[io_szllong_stord - 1 - j];
           }
           for (j = 0; j < io_szllong_stord; j++) {
               s[j] = t[j];
           }
           s += io_szllong_stord;
       }
       if (io_szllong_stord == io_szllong) { 
           memcpy(a6ll, ptr, (size_t)(3 * io_szllong_stord));
           num_attr = a6ll[2];
       }
       else if (io_szllong_stord == io_szint) { 
           memcpy(a6i, ptr, (size_t)(3 * io_szllong_stord));
           num_attr = a6i[2];
       }
       else if (io_szllong_stord == io_szlong) {  
           memcpy(a6l, ptr, (size_t)(3 * io_szllong_stord));
           num_attr = a6l[2];
       }
       for (k = 0; k < num_attr; k++) { 
           for (j = 0; j < io_szllong_stord; j++) {
               t[j] = s[io_szllong_stord - 1 - j];
           }
           for (j = 0; j < io_szllong_stord; j++) {
               s[j] = t[j];
           }
           if (io_szllong_stord == io_szllong) { 
               memcpy(a6ll, s, (size_t)io_szllong_stord);
               slen = a6ll[0];
           }
           else if (io_szllong_stord == io_szint) { 
               memcpy(a6i, s, (size_t)io_szllong_stord);
               slen = a6i[0];
           }
           else if (io_szllong_stord == io_szlong) { 
               memcpy(a6l, s, (size_t)io_szllong_stord);
               slen = a6l[0];
           }
           s += (io_szllong_stord + slen);
           ptr = s;                  /* datatype, size of value */  
           for (i = 0; i < 2; i++) { 
               for (j = 0; j < io_szllong_stord; j++) {
                   t[j] = s[io_szllong_stord - 1 - j];
               }
               for (j = 0; j < io_szllong_stord; j++) {
                   s[j] = t[j];
               }
               s += io_szllong_stord;
           }
           if (io_szllong_stord == io_szllong) {
               memcpy(a6ll, ptr, (size_t)(io_szllong_stord + io_szllong_stord));
               datatype_stord = (bio_Data_Type) ((int)a6ll[0]);
               size_value = a6ll[1];
           }
           else if (io_szllong_stord == io_szint) {
               memcpy(a6i, ptr, (size_t)(io_szllong_stord + io_szllong_stord));
               datatype_stord = (bio_Data_Type) ((int)a6i[0]);
               size_value = a6i[1];
           }
           else if (io_szllong_stord == io_szlong) {
               memcpy(a6l, ptr, (size_t)(io_szllong_stord + io_szllong_stord));
               datatype_stord = (bio_Data_Type) ((int)a6l[0]);
               size_value = a6l[1];
           }
           if (datatype_stord == bio_int) {
               size = size_value /io_szint_stord;
               io_swap(s, (int)io_szint_stord, size);
           }
           else if (datatype_stord == bio_long) { 
               size = size_value /io_szlong_stord;
               io_swap(s, (int)io_szlong_stord, size);
           }
           else if (datatype_stord == bio_long_long) { 
               size = size_value /io_szllong_stord;
               io_swap(s, (int)io_szllong_stord, size);
           }
           else if (datatype_stord == bio_float) { 
               size = size_value /io_szfloat_stord;
               io_swap(s, (int)io_szfloat_stord, size);
           }
           else if (datatype_stord == bio_double) {    
               size = size_value /io_szdouble_stord;
               io_swap(s, (int)io_szdouble_stord, size);
           }
           else if (datatype_stord == bio_long_double) {
               size = size_value /io_szldouble_stord;
               io_swap(s, (int)io_szldouble_stord, size);
           }
           s += size_value;
       } 
   }

   return 0;
 }

void io_list_null(int n, bio_List_Struct *list)
{
     int i;
     bio_List_Struct *obj;
  
     for (i = 0; i < n; i++) { 
         obj = list + i;
         obj->narrays = 0;
         obj->names = NULL;
         obj->datatypes = NULL;
         obj->sizes = NULL;

         obj->name = NULL;
         obj->id = -1;
         obj->type = bio_object_type_undefined; 
         obj->datatype = bio_datatype_invalid; 
         obj->value   = NULL;
     }
     return;
 } 

int bio_free_list(int n, bio_List_Struct *list) 
{
    int i;
    bio_List_Struct *obj;

    for (i = 0; i < n; i++) {
        obj = list + i;
        if (obj->names) free(obj->names);
        obj->names = NULL;
        if (obj->datatypes) free(obj->datatypes);
        obj->datatypes = NULL;
        if (obj->sizes) free(obj->sizes);
        obj->sizes = NULL;

        if (obj->name) free(obj->name);
        obj->name = NULL;
        if (obj->value) free(obj->value);
        obj->value = NULL;
    }
    return 0;
 } 

/****************************************************************************************/
char *mystrchar(char *s, char separator)
{
    char *c;
    
    c = s;
    while (*c != '\0') {
       if (*c == separator) { 
           return c;
       }
       c++;
    }
    return NULL;
 }  
/****************************************************************************************/

int bio_debug(int fileid, int id)
{
    io_Obj *obj;
    io_File *file;
    io_Tree *tree;

#ifdef MPI 
    MPI_Barrier(bio_comm);
#endif
    
    file = NULL;
    tree = io_find_file(fileid, &file);
    if (fileid != id) {
        tree = io_find_tree(file, id);
    }
    if (!tree || !file) { 
        if (!mype) printf("ERROR: null file/tree in bio_debug\n");
    }
    else { 
       obj = file->object;
    }
#ifdef MPI
    MPI_Barrier(bio_comm);
#endif

    return 0;
} 
     
