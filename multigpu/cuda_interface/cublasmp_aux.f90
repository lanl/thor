!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file cublasmp_aux.f90
!! @author Ismael Djibrilla Boureima, Oleg Korobkin
!! @date  December 2023
!! @brief Fortran interface to C wrappers for cuBLASMp library functions
!!
module cublasmp_aux
use iso_c_binding
implicit none

   type(c_ptr)      :: cublasmp_handle
   logical, private :: cublasmp_initialized = .false.

   private :: cublasMp_stop_if_error, &
              show_cublasMpHandle_c, &
              show_cublasMpGrid_c, &
              show_cublasMpMatrixDescriptor_c, &
              set_cublasMpHandle_c, &
              destroy_cublasMpHandle_c, &
              set_cublasMpGrid_c, &
              destroy_cublasMpGrid_c, &
              set_cublasMpMatrixDesc_c, &
              destroy_cublasMpMatrixDesc_c, &
              cublasMpNumroc_int64_c, &
              redistribute_rectMatrix_c, &
              cublasMpGeadd_c, &
              cublasMpGemm_c

interface
!///////////////////////////////////////////////////////////////
!// SHOW DSTRUC
!//////////////////////////////////////////////////////////////
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-Prototype: int show_cublasMpHandle(cublasMpHandle_t *handle);
   integer(c_int) function show_cublasMpHandle_c(handle) bind(C,name="show_cublasMpHandle")
   use iso_c_binding
     type(c_ptr)             :: handle
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int show_cublasMpGrid(cublasMpGrid_t *grid);
   integer(c_int) function show_cublasMpGrid_c(grid) bind(C,name="show_cublasMpGrid")
   use iso_c_binding
     type(c_ptr)             :: grid
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int show_cublasMpMatrixDescriptor(cublasMpMatrixDescriptor_t *descriptor);;
   integer(c_int) function show_cublasMpMatrixDescriptor_c(descriptor)&
                  bind(C,name="show_cublasMpMatrixDescriptor")
   use iso_c_binding
     type(c_ptr)             :: descriptor
   end function
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   !//                                                   SET/DESTROY  DSTRUCT
   !///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ! C-PROTOTYPE: int set_cublasMpHandle(cublasMpHandle_t *handle, cudaStream_t *stream, int VRBZ);
   integer(c_int) function set_cublasMpHandle_c(handle, stream, VRBZ) bind(C,name="set_cublasMpHandle")
   use iso_c_binding
     type(c_ptr)             :: handle, stream
     integer(c_int), value   :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int destroy_cublasMpHandle(cublasMpHandle_t *handle, int VRBZ);
   integer(c_int) function destroy_cublasMpHandle_c(handle, VRBZ) bind(C,name="destroy_cublasMpHandle")
   use iso_c_binding
     type(c_ptr)             :: handle
     integer(c_int), value   :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int set_cublasMpGrid(cublasMpHandle_t *handle, cublasMpGrid_t *grid, cal_comm_t *cal_comm,
   !                      int64_t numRowDevices, int64_t numColDevices, int VRBZ);
   integer(c_int) function set_cublasMpGrid_c(grid, cal_comm, numRowDevices, numColDevices,&
                                 VRBZ) bind(C,name="set_cublasMpGrid")
   use iso_c_binding
     type(c_ptr)               :: grid, cal_comm
     integer(c_int64_t), value :: numRowDevices, numColDevices
     integer(c_int), value     :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int destroy_cublasMpGrid(cublasMpHandle_t *handle, cublasMpGrid_t *grid,  int VRBZ);
   integer(c_int) function destroy_cublasMpGrid_c(grid, VRBZ) bind(C,name="destroy_cublasMpGrid")
   use iso_c_binding
     type(c_ptr)             :: grid
     integer(c_int), value   :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int set_cublasMpMatrixDesc(cublasMpHandle_t *handle, cublasMpMatrixDescriptor_t *descrA,
   !                        cublasMpGrid_t *gridA, int64_t M, int64_t N, int64_t MB_A, int64_t NB_A,
   !                        int64_t numRowDevices, int64_t numColDevices, int64_t dev_row_idx,
   !                        int64_t dev_col_idx, int VRBZ);
   integer(c_int) function  set_cublasMpMatrixDesc_c(handle, descrA, gridA, M, N, MB_A, NB_A, numRowDevices,&
                                 dev_row_idx, VRBZ ) bind(C,name="set_cublasMpMatrixDesc")
   use iso_c_binding
     type(c_ptr)               :: handle, descrA, gridA
     integer(c_int64_t), value :: M,N,MB_A,NB_A, numRowDevices, dev_row_idx
     integer(c_int), value     :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int destroy_cublasMpMatrixDesc(cublasMpHandle_t *handle, cublasMpMatrixDescriptor_t *descrA, int VRBZ);
   integer(c_int) function destroy_cublasMpMatrixDesc_c(handle, descrA, VRBZ) bind(C,name="destroy_cublasMpMatrixDesc")
   use iso_c_binding
     type(c_ptr)             :: handle, descrA
     integer(c_int), value   :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int64_t cublasMpNumroc(int64_t M, int64_t MB_A,
   !                                     int64_t dev_idx, uint64_t RC_SRCA, int64_t numDevices);
   pure integer(c_int64_t) function cublasMpNumroc_int64_c(M, MB_A, dev_idx, RC_SRCA, numDevices) bind(C,name="cublasMpNumroc")
   use iso_c_binding
     integer(c_int64_t), value :: M, MB_A, dev_idx, RC_SRCA, numDevices
   end function

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int redistribute_rectMatrix(cublasMpHandle_t *handle,
   !         cal_comm_t *cal_comm, cudaStream_t *stream,
   !         int64_t M, int64_t N,
   !         cublasMpMatrixDescriptor_t *descrA, double* d_A,
   !         cublasMpMatrixDescriptor_t *descrB, double* d_B, int VRBZ) {
   integer(c_int) function redistribute_rectMatrix_c (handle, cal_comm, stream, &
        M, N, descrA, d_A, descrB, d_B, VRBZ) bind(C,name="redistribute_rectMatrix")
   use iso_c_binding
     type(c_ptr)               :: handle, descrA, descrB, cal_comm, stream
     type(c_ptr), value        :: d_A, d_B
     integer(c_int64_t), value :: M, N
     integer(c_int), value     :: VRBZ
   end function
   !//////////////////////////////////////////////////////////////////////////
   !//                            MATRIX ADDITION
   !//////////////////////////////////////////////////////////////////////////
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int cublasMpGeadd_wrapper(cublasMpHandle_t *handle,
   !                        cal_comm_t *cal_comm,
   !                        cudaStream_t *stream,
   !                        int64_t M, int64_t N, int trA, double alpha, double beta,
   !                        cublasMpMatrixDescriptor_t *descrA, double* d_A,
   !                        cublasMpMatrixDescriptor_t *descrC, double* d_C, int VRBZ);
   integer(c_int) function cublasMpGeadd_c(handle, cal_comm, stream, &
        M, N, trA, alpha, beta, descrA, d_A, descrC, d_C, VRBZ) &
        bind(C,name="cublasMpGeadd_wrapper")
   use iso_c_binding
   type(c_ptr)               :: handle, descrA, descrC, cal_comm, stream
   type(c_ptr), value        :: d_A, d_C
   integer(c_int64_t), value :: M, N
   double precision, value   :: alpha, beta
   integer(c_int), value     :: trA, VRBZ
   end function
   !/////////////////////////////////////////////////////////////////////////
   !//                      MATRIX MULTIPLICATION
   !/////////////////////////////////////////////////////////////////////////
   !> C-PROTOTYPE: int cublasMpGemm_wrapper(cublasMpHandle_t *handle,
   !!     cal_comm_t *cal_comm, cudaStream_t *stream,
   !!     int64_t M, int64_t N, int64_t K, int trA, int trB,
   !!     cublasMpMatrixDescriptor_t *descrA, double* d_A,
   !!     cublasMpMatrixDescriptor_t *descrB, double* d_B,
   !!     cublasMpMatrixDescriptor_t *descrC, double* d_C, int VRBZ)
   integer(c_int) function cublasMpGemm_c(handle, cal_comm, stream, &
        M, N, K, trA, trB, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ) &
        bind(C,name="cublasMpGemm_wrapper")
   use iso_c_binding
     type(c_ptr)               :: handle, cal_comm, stream
     type(c_ptr)               :: descrA, descrB, descrC
     type(c_ptr), value        :: d_A, d_B, d_C
     integer(c_int64_t), value :: M, N, K
     integer(c_int), value     :: trA, trB, VRBZ
   end function

end interface
contains

   subroutine cublasMp_stop_if_error(err, msg)
   use, intrinsic :: iso_fortran_env, only: error_unit
   use iso_c_binding
   integer(c_int), intent(in) :: err
   character(*), intent(in)   :: msg
      if (err.ne.0) then
         write(error_unit,'(A,I8)') msg, err
         error stop
      endif
   end subroutine


   subroutine cublasmp_init(VRBZ_)
   use iso_c_binding
   use cuda_comm_aux, only: local_stream
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= set_cublasMpHandle_c(cublasmp_handle, local_stream, VRBZ)
      call cublasMp_stop_if_error(err, "[!][cublasmp_init], err=")
      cublasmp_initialized= .true.
   end subroutine


   logical function cublasmp_is_initialized()
   cublasmp_is_initialized= cublasmp_initialized
   end function


   subroutine show_cublasmp_handle()
   use iso_c_binding
   !
   integer(c_int) :: err
      err= show_cublasMpHandle_c(cublasmp_handle)
      call cublasMp_stop_if_error(err, "[!][show_cublasmp_handle], err=")
   end subroutine


   subroutine cublasmp_destroy(VRBZ_)
   use iso_c_binding
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= destroy_cublasMpHandle_c(cublasmp_handle, VRBZ)
      call cublasMp_stop_if_error(err, "[!][cublasmp_destroy], err=")
      cublasmp_initialized= .false.
   end subroutine


   type(c_ptr) function set_cublasMpGrid(nRowDevs, nColDevs, VRBZ_) result(grid)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle
   integer, intent(IN)           :: nRowDevs, nColDevs
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err

      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= set_cublasMpGrid_c(grid, cal_handle, &
                            int(nRowDevs,c_int64_t), &
                            int(nColDevs,c_int64_t), VRBZ)
      call cublasMp_stop_if_error(err, "[!][set_cublasMpGrid], err=")
   end function set_cublasMpGrid


   subroutine show_cublasMpGrid(grid)
   use iso_c_binding
   type(c_ptr), intent(IN) :: grid
   !
   integer(c_int) :: c_err
      c_err = show_cublasMpGrid_c(grid)
      call cublasMp_stop_if_error(c_err, '[!][show_cublasMpGrid], err=')
   end subroutine


   subroutine destroy_cublasMpGrid(grid, VRBZ_)
   use iso_c_binding
   type(c_ptr), intent(IN) :: grid
   integer, optional       :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= destroy_cublasMpGrid_c(grid, VRBZ)
      call cublasMp_stop_if_error(err, "[!][destroy_cublasMpGrid], err=")
   end subroutine destroy_cublasMpGrid


   type(c_ptr) function  set_cublasMpMatrixDesc(grid, &
                         M, N, MB, NB, numRowDevices, &
                         dev_row_idx, VRBZ_ ) result(descr)
   use iso_c_binding
   type(c_ptr), intent(IN) :: grid
   integer, intent(IN)     :: M, N, MB, NB
   integer, intent(IN)     :: numRowDevices
   integer, intent(IN)     :: dev_row_idx
   integer, optional       :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= set_cublasMpMatrixDesc_c(cublasmp_handle, descr, grid, &
           int(M, c_int64_t), int(N, c_int64_t), &
           int(MB,c_int64_t), int(NB,c_int64_t), &
           int(numRowDevices, c_int64_t), &
           int(dev_row_idx, c_int64_t), VRBZ)
      call cublasMp_stop_if_error(err, "[!][set_cublasMpMatrixDesc], err=")
   end function


   subroutine show_cublasMpMatrixDesc(descriptor)
   use iso_c_binding
   type(c_ptr), intent(IN) :: descriptor
   !
   integer(c_int) :: err
      err= show_cublasMpMatrixDescriptor_c(descriptor)
      call cublasMp_stop_if_error(err, "[!][show_cublasMpMatrixDesc], err=")
   end subroutine



   subroutine destroy_cublasMpMatrixDesc(descr, VRBZ_)
   use iso_c_binding
   type(c_ptr), intent(IN) :: descr
   integer, optional       :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= destroy_cublasMpMatrixDesc_c(cublasmp_handle, descr, VRBZ)
      call cublasMp_stop_if_error(err, "[!][destroy_cublasMpMatrixDesc], err=")
   end subroutine


   pure elemental integer function &
   cublasMpNumroc(M, MB_A, dev_idx, RC_SRCA, numDevices) result(retval)
   use iso_c_binding
   integer, intent(IN) :: M, MB_A, dev_idx, RC_SRCA, numDevices

     retval = int(cublasMpNumroc_int64_c( &
                  int(M, c_int64_t), int(MB_A, c_int64_t), &
                  int(dev_idx, c_int64_t), int(RC_SRCA, c_int64_t), &
                  int(numDevices, c_int64_t)))
   end function


   !> Number of rows or columns for a matrix size M distributed with blocks MB
   pure elemental integer function numroc(M, MB, rank, nranks)
   integer, intent(IN) :: M, MB, rank, nranks
      numroc = M/(MB*nranks)*MB &
             + min(MB, max(0, mod(M, MB*nranks) - mod(rank,nranks)*MB))
   end function


   !> Use local index to compute global index in a distributed matrix
   pure elemental integer function iloc2glob(j, MB, rank, nranks)
   integer, intent(IN) :: j, MB, rank, nranks
      iloc2glob = 1 + MB*(mod(rank,nranks) + ((j - 1)/MB)*nranks) &
                + mod(j - 1, MB)
   end function


   !> Use global index to compute on which rank an element resides
   pure elemental integer function iglob2rank(J, MB, nranks)
   integer, intent(IN) :: J, MB, nranks
      iglob2rank = mod((J - 1)/MB, nranks)
   end function


   !> Use global index and rank to compute local intex of an element
   pure elemental integer function iglob2loc(J, MB, rank, nranks)
   integer, intent(IN) :: J, MB, rank, nranks
      iglob2loc = 1 + ((J - 1 - mod(rank,nranks)*MB)/(MB*nranks))*MB &
                + mod(J - 1, MB)
   end function


   integer function &
   AxB(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_) result(retval)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN):: descrA, descrB, descrC
   double precision, device, target, intent(IN):: d_A(:), d_B(:)
   double precision, device, target, intent(INOUT):: d_C(:)
   integer, intent(IN) :: M, N, K
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ

      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      retval = int(cublasMpGemm_c(cublasmp_handle, cal_handle, local_stream, &
                  int(M,c_int64_t),int(N,c_int64_t),int(K,c_int64_t), 0, 0,   &
                  descrA, C_LOC(d_A), descrB, C_LOC(d_B), descrC, C_LOC(d_C), &
                  int(VRBZ,c_int)))
   end function


   integer function &
   ATxB(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_) result(retval)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN):: descrA, descrB, descrC
   double precision, device, target, intent(IN):: d_A(:), d_B(:)
   double precision, device, target, intent(INOUT):: d_C(:)
   integer, intent(IN) :: M, N, K
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ

      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      retval = int(cublasMpGemm_c(cublasmp_handle, cal_handle, local_stream, &
                  int(M,c_int64_t),int(N,c_int64_t),int(K,c_int64_t), 1, 0,   &
                  descrA, C_LOC(d_A), descrB, C_LOC(d_B), descrC, C_LOC(d_C), &
                  int(VRBZ,c_int)))
   end function


   integer function &
   AxBT(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_) result(retval)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN):: descrA, descrB, descrC
   double precision, device, target, intent(IN):: d_A(:), d_B(:)
   double precision, device, target, intent(INOUT):: d_C(:)
   integer, intent(IN) :: M, N, K
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ

      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      retval = int(cublasMpGemm_c(cublasmp_handle, cal_handle, local_stream, &
                  int(M,c_int64_t),int(N,c_int64_t),int(K,c_int64_t), 0, 1,   &
                  descrA, C_LOC(d_A), descrB, C_LOC(d_B), descrC, C_LOC(d_C), &
                  int(VRBZ,c_int)))
   end function


   integer function &
   ATxBT(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_) result(retval)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN):: descrA, descrB, descrC
   double precision, device, target, intent(IN):: d_A(:), d_B(:)
   double precision, device, target, intent(INOUT):: d_C(:)
   integer, intent(IN) :: M, N, K
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ

      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      retval = int(cublasMpGemm_c(cublasmp_handle, cal_handle, local_stream, &
                  int(M,c_int64_t),int(N,c_int64_t),int(K,c_int64_t), 1, 1,   &
                  descrA, C_LOC(d_A), descrB, C_LOC(d_B), descrC, C_LOC(d_C), &
                  int(VRBZ,c_int)))
   end function


   integer function &
   redistribute_rectMatrix(M, N, descrA, d_A, descrB, d_B, VRBZ_) result(retval)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   integer, intent(IN)       :: M, N
   type(c_ptr), intent(IN)   :: descrA, descrB
   double precision, device, target, intent(IN)   :: d_A(:)
   double precision, device, target, intent(INOUT):: d_B(:)
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ

      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      retval = int(redistribute_rectMatrix_c(cublasmp_handle, cal_handle, &
                   local_stream, int(M,c_int64_t), int(N,c_int64_t),  &
                   descrA, C_LOC(d_A), descrB, C_LOC(d_B), VRBZ))
   end function

   !> General matrix addition: C = alpha*OP(A) + beta*C
   integer function &
   cublasMpGeadd(M, N, trA, alpha, beta, descrA, d_A, descrC, d_C, VRBZ_) result(err)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN) :: descrA, descrC
   double precision, pointer, device, intent(IN)    :: d_A(:)
   double precision, pointer, device, intent(INOUT) :: d_C(:)
   integer, intent(IN)           :: M, N, trA
   double precision, intent(IN)  :: alpha, beta
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err = int(cublasMpGeadd_c(&
                cublasmp_handle, cal_handle, local_stream, &
                int(M,c_int64_t), int(N,c_int64_t), &
                int(trA,c_int), alpha, beta, &
                descrA, C_LOC(d_A), &
                descrC, C_LOC(d_C), VRBZ))
   end function


   !> Matrix addition: C = C + A
   integer function &
   CaddA(M, N, descrA, d_A, descrC, d_C, VRBZ_) result(err)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN) :: descrA, descrC
   double precision, pointer, device, intent(IN)    :: d_A(:)
   double precision, pointer, device, intent(INOUT) :: d_C(:)
   integer, intent(IN)                              :: M, N
   integer, intent(IN), optional                    :: VRBZ_
   integer(c_int) :: VRBZ
   integer(c_int), parameter :: TRANS_OP_N = 0
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err = int(cublasMpGeadd_c(&
                cublasmp_handle, cal_handle, local_stream, &
                int(M,c_int64_t), int(N,c_int64_t), &
                TRANS_OP_N, 1d0, 1d0, &
                descrA, C_LOC(d_A), &
                descrC, C_LOC(d_C), VRBZ))
   end function


   !> Matrix addition with a transposed A: C = C + A^T
   integer function &
   CaddAT(M, N, descrA, d_A, descrC, d_C, VRBZ_) result(err)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN) :: descrA, descrC
   double precision, pointer, device, intent(IN)    :: d_A(:)
   double precision, pointer, device, intent(INOUT) :: d_C(:)
   integer, intent(IN)                              :: M, N
   integer, intent(IN), optional                    :: VRBZ_
   integer(c_int) :: VRBZ
   integer(c_int), parameter :: TRANS_OP_T = 1
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err = int(cublasMpGeadd_c(&
                cublasmp_handle, cal_handle, local_stream, &
                int(M,c_int64_t), int(N,c_int64_t), &
                TRANS_OP_T, 1d0, 1d0, &
                descrA, C_LOC(d_A), &
                descrC, C_LOC(d_C), VRBZ))
   end function


   !> Matrix subtraction: C = C - A
   integer function &
   CsubA(M, N, descrA, d_A, descrC, d_C, VRBZ_) result(err)
   use iso_c_binding
   use cuda_comm_aux, only: cal_handle, local_stream
   type(c_ptr), intent(IN) :: descrA, descrC
   double precision, pointer, device, intent(IN)    :: d_A(:)
   double precision, pointer, device, intent(INOUT) :: d_C(:)
   integer, intent(IN)                              :: M, N
   integer, intent(IN), optional                    :: VRBZ_
   integer(c_int) :: VRBZ
   integer(c_int), parameter :: TRANS_OP_N = 0
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err = int(cublasMpGeadd_c(&
                cublasmp_handle, cal_handle, local_stream, &
                int(M,c_int64_t), int(N,c_int64_t), &
                TRANS_OP_N, -1d0, 1d0, &
                descrA, C_LOC(d_A), &
                descrC, C_LOC(d_C), VRBZ))
   end function

end module cublasmp_aux
