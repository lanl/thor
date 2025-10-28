!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file cusolvermp_aux.f90
!! @author Ismael Djibrilla Boureima, Oleg Korobkin
!! @date  November 2023
!! @brief Fortran interface to C wrappers for cuSolverMp library functions
!!
module cusolvermp_aux
use iso_c_binding
implicit none

   type(c_ptr)      :: cusolvermp_handle
   logical, private :: cusolvermp_initialized = .false.

   private :: cusolverMp_stop_if_error, &
              show_device_int_array_c, &
              show_device_real_array_c, &
              show_device_double_array_c, &
              print_host_matrix_c, &
              show_cusolverMpHandle_c, &
              show_cusolverMpGrid_c, &
              show_cusolverMpMatrixDescriptor_c, &
              set_dist_mat_row_col_params_c, &
              set_cusolverMpHandle_c, &
              destroy_cusolverMpHandle_c, &
              set_cusolverMpGrid_c, &
              destroy_cusolverMpGrid_c, &
              set_cusolverMpMatrixDesc_c, &
              destroy_cusolverMpMatrixDesc_c, &
              check_qr_buffsize_c, &
              dist_geqrf_c, &
              dist_ormqr_c, &
              dist_sytrd_c, &
              dist_stedc_c


interface
!![!] Notes:
!! Wrapping C   ->   Fortran:
!! TYPE C       -> TYPE FORTRAN
!! void*        -> type(c_ptr), value  // See handling of 'void* dst' and 'void* src) in the wrapping of cudaMemcpy()
!! void**       -> type(c_ptr)         // See handling of 'void** devPtr' in the wrapping of cudaMalloc
!! double*      -> type(c_ptr), value  // See handling of 'float* A'' in the wrapping of cusolverDnSgetrf()

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-Prototype: int show_device_int_array(int *d_array, int size_array);
   integer(c_int) function show_device_int_array_c(d_array, size_array)&
                  bind(C,name="show_device_int_array")
   use iso_c_binding
   type(c_ptr), value      :: d_array
   integer(c_int), value   :: size_array
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int show_device_real_array(float *d_array, int size_array);
   integer(c_int) function show_device_real_array_c(d_array, size_array)&
                  bind(C,name="show_device_real_array")
   use iso_c_binding
   type(c_ptr), value      :: d_array
   integer(c_int), value   :: size_array
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int show_device_double_array(double *d_array, int size_array);
   integer(c_int) function show_device_double_array_c(d_array, size_array)&
                  bind(C,name="show_device_double_array")
   use iso_c_binding
   type(c_ptr), value      :: d_array
   integer(c_int), value   :: size_array
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: static void print_host_matrix(int64_t M, int64_t N, double* A, int64_t lda, const char* msg);
   subroutine print_host_matrix_c(M,N,A,lda, msg) bind(C,name="print_host_matrix")
   use iso_c_binding
   integer(c_int64_t) :: M,N,lda
   real(c_double)     :: A(*)
   character(c_char)  :: msg(*)
   end subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int show_cusolverMpHandle(cusolverMpHandle_t *cusolverMpHandle);;
   integer(c_int) function show_cusolverMpHandle_c(cusolverMpHandle) bind(C,name="show_cusolverMpHandle")
   use iso_c_binding
   type(c_ptr)              :: cusolverMpHandle
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: show_cusolverMpGrid(cusolverMpGrid_t *cusolverMpGrid);;
   integer(c_int) function show_cusolverMpGrid_c(cusolverMpGrid) bind(C,name="show_cusolverMpGrid")
   use iso_c_binding
   type(c_ptr)              :: cusolverMpGrid
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: show_cusolverMpMatrixDescriptor(cusolverMpMatrixDescriptor_t *cusolverMpMatrixDescriptor);;
   integer(c_int) function show_cusolverMpMatrixDescriptor_c(cusolverMpMatrixDescriptor)&
                  bind(C,name="show_cusolverMpMatrixDescriptor")
   use iso_c_binding
   type(c_ptr)              :: cusolverMpMatrixDescriptor
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int set_dist_mat_row_col_params(int64_t M, int64_t N, int64_t MA, int64_t NA,
   !              int rankId, int64_t numRowDevices, int64_t numColDevices,
   !              int *ML, int *NL, int VRBZ){
   integer(c_int) function set_dist_mat_row_col_params_c(M,N,MA,NA, rankId, &
                                  numRowDevices, numColDevices, ML, NL, VRBZ) &
                                  bind(C,name="set_dist_mat_row_col_params")
   use iso_c_binding
   type(c_ptr), value           :: ML, NL
   integer(c_int), value        :: rankId, VRBZ
   integer(c_int64_t), value    :: M, N, MA, NA, numRowDevices, numColDevices
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int set_cusolverMpHandle(cusolverMpHandle_t *cusolverMpHandle, int localDeviceId, cudaStream_t *stream, int VRBZ);
   integer(c_int) function set_cusolverMpHandle_c(cusolverMpHandle, localDeviceId, stream, VRBZ) bind(C,name="set_cusolverMpHandle")
   use iso_c_binding
   type(c_ptr)           :: stream, cusolverMpHandle
   integer(c_int), value :: localDeviceId, VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int destroy_cusolverMpHandle(cusolverMpHandle_t *cusolverMpHandle, int VRBZ);
   integer(c_int) function destroy_cusolverMpHandle_c(cusolverMpHandle, VRBZ) bind(C,name="destroy_cusolverMpHandle")
   use iso_c_binding
   type(c_ptr)           :: cusolverMpHandle
   integer(c_int), value :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int set_cusolverMpGrid(cusolverMpHandle_t *cusolverMpHandle, cusolverMpGrid_t *gridA,
   !                     cal_comm_t *cal_comm, int64_t numRowDevices, int64_t numColDevices, int VRBZ);
   integer(c_int) function set_cusolverMpGrid_c(cusolverMpHandle, gridA, cal_comm, numRowDevices, numColDevices,&
                                              VRBZ) bind(C,name="set_cusolverMpGrid")
   use iso_c_binding
   type(c_ptr)                :: cusolverMpHandle, cal_comm, gridA
   integer(c_int64_t), value  :: numRowDevices, numColDevices
   integer(c_int), value      ::  VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int destroy_cusolverMpGrid(cusolverMpGrid_t *gridA, int VRBZ);
   integer(c_int) function destroy_cusolverMpGrid_c(gridA, VRBZ) bind(C,name="destroy_cusolverMpGrid")
   use iso_c_binding
   type(c_ptr)           :: gridA
   integer(c_int), value :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int set_cusolverMpMatrixDesc(cusolverMpMatrixDescriptor_t *descrA,
   !                        cusolverMpGrid_t *gridA, int64_t M, int64_t N,
   !                        int64_t MA, int64_t NA,  int64_t numRowDevices, int dev_row_idx, int VRBZ)
   integer(c_int) function set_cusolverMpMatrixDesc_c(descrA,gridA,M,N,MA,NA,numRowDevices,&
                              dev_row_idx, VRBZ) bind(C,name="set_cusolverMpMatrixDesc")
   use iso_c_binding
   type(c_ptr)               :: descrA, gridA
   integer(c_int64_t), value :: M,N,MA,NA,numRowDevices
   integer(c_int), value     :: dev_row_idx, VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int destroy_cusolverMpMatrixDesc(cusolverMpMatrixDescriptor_t *descrA, int VRBZ);;
   integer(c_int) function destroy_cusolverMpMatrixDesc_c(descrA, VRBZ) bind(C,name="destroy_cusolverMpMatrixDesc")
   use iso_c_binding
   type(c_ptr)    :: descrA
   integer(c_int), value :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int check_qr_buffsize(cusolverMpHandle_t *cusolverMpHandle,
   !                   const int64_t M, const int64_t N,
   !                   const int64_t IA, const int64_t JA,
   !                   cusolverMpMatrixDescriptor_t *descrA, void* d_A);
   integer(c_int) function check_qr_buffsize_c(cusolverMpHandle, M, N, IA, JA, descrA, d_A) bind(C,name="check_qr_buffsize")
   use iso_c_binding
   type(c_ptr)               :: cusolverMpHandle, descrA
   type(c_ptr)               :: d_A
   integer(c_int64_t), value :: M,N,IA,JA
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int dist_geqrf(cal_comm_t *cal_comm, cudaStream_t *stream,
   !                          cusolverMpHandle_t *cusolverMpHandle,
   !                          cusolverMpMatrixDescriptor_t *descrA, double* d_A,
   !                          const int64_t M, const int64_t N, double* d_tau)
   integer(c_int) function dist_geqrf_c(stream, cusolverMpHandle, &
                                      descrA, d_A, M, N, d_tau, VRBZ) bind(C,name="dist_geqrf")
   use iso_c_binding
   type(c_ptr)               :: stream, cusolverMpHandle, descrA
   integer(c_int64_t), value :: M, N
   type(c_ptr), value        :: d_A, d_tau
   integer(c_int), value     :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE:
   ! int dist_ormqr(cal_comm_t *cal_comm, cudaStream_t *stream, cusolverMpHandle_t *cusolverMpHandle,
   !        cusolverMpMatrixDescriptor_t *descA, double* d_A, const int iside, const int itrans,
   !        const int64_t M, const int64_t N, cusolverMpMatrixDescriptor_t *descC, double* d_C,
   !        const int64_t K, const double *d_tau) {
   integer(c_int) function dist_ormqr_c(stream, cusolverMpHandle, &
                                   descrA, d_A, iside, itrans, M, N, &
                                   descrC, d_C, K, d_tau, verb) &
                                   bind(C,name="dist_ormqr")
   use iso_c_binding
   type(c_ptr)               :: cusolverMpHandle, stream
   type(c_ptr)               :: descrA, descrC
   integer(c_int64_t), value :: M, N, K
   type(c_ptr), value        :: d_A, d_C, d_tau
   integer(c_int), value     :: iside, itrans
   integer(c_int), value     :: verb
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE:
   ! int dist_mpsytrd(cal_comm_t *cal_comm, cudaStream_t *stream, cusolverMpHandle_t *cusolverMpHandle,
   !        const int64_t N, cusolverMpMatrixDescriptor_t *descA, double* d_A, double* d_d, double* d_e,
   !        double* d_tau)
   integer(c_int) function dist_sytrd_c(stream, cusolverMpHandle, &
                                   N, descrA, d_A, d_d, d_e, d_tau) &
                                   bind(C,name="dist_mpsytrd")
   use iso_c_binding
   type(c_ptr)               :: cusolverMpHandle, stream
   type(c_ptr)               :: descrA
   integer(c_int64_t), value :: N
   type(c_ptr), value        :: d_A, d_d, d_e, d_tau
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE:
   ! int dist_mpstedc(cusolverMpHandle_t *cusolverMpHandle,
   !        const int64_t N, double* d_d, double* d_e,
   !        cusolverMpMatrixDescriptor_t *descQ, double* d_Q)
   integer(c_int) function dist_stedc_c(cusolverMpHandle, N, d_d, d_e, descrQ, d_Q, verb) &
                                   bind(C,name="dist_mpstedc")
   use iso_c_binding
   type(c_ptr)               :: cusolverMpHandle
   type(c_ptr)               :: descrQ
   integer(c_int64_t), value :: N
   type(c_ptr), value        :: d_d, d_e, d_Q
   integer(c_int), value     :: verb
   end function

end interface

contains

   subroutine cusolverMp_stop_if_error(err, msg)
   use iso_c_binding
   use, intrinsic :: iso_fortran_env, only: error_unit
   integer(c_int), intent(in) :: err
   character(*), intent(in)   :: msg
      if (err.ne.0) then
         write(error_unit,'(A,I8)') msg, err
         error stop
      endif
   end subroutine


   subroutine  cusolvermp_init(localDeviceId, VRBZ_) result(handle)
   use iso_c_binding
   use cuda_comm_aux, only: local_stream
   integer, intent(IN)     :: localDeviceId
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if(present(VRBZ_)) VRBZ= int(VRBZ_, c_int)
      err= set_cusolverMpHandle_c(cusolvermp_handle, int(localDeviceId,c_int), &
                                local_stream, VRBZ)
      call cusolverMp_stop_if_error(err, "[!][set_cusolverMpHandle], err=")
   end subroutine cusolvermp_init


   logical function cusolvermp_is_initialized()
   cusolvermp_is_initialized= cusolvermp_initialized
   end function


   subroutine show_cusolverMpHandle()
   use iso_c_binding
   !
   integer(c_int) :: err
      err= show_cusolverMpHandle_c(cusolvermp_handle)
      call cusolverMp_stop_if_error(err, "[!][show_cusolvermp_handle], err=")
   end subroutine


   subroutine cusolvermp_destroy(VRBZ_)
   use iso_c_binding
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= destroy_cusolverMpHandle_c(cusolvermp_handle, VRBZ)
      call cusolverMp_stop_if_error(err, "[!][cusolvermp_destroy], err=")
      cusolvermp_initialized= .false.
   end subroutine


   !> Compute local matrix size
   subroutine set_dist_mat_row_col(M, N, MA, NA, rank, nRowDevs, nColDevs, &
                                   ML, NL, VRBZ_)
   use iso_c_binding
   integer, intent(IN) :: M, N               !< global matrix dimensions: M x N
   integer, intent(IN) :: MA, NA             !< blocking dimensions
   integer, intent(IN) :: rank               !< current rank
   integer, intent(IN) :: nRowDevs, nColDevs !< processor grid dimensions
   integer, intent(OUT)  :: ML, NL   !< [OUT] local matrix dimensions: ML x NL
   integer, intent(IN), optional :: VRBZ_
   !
   integer :: VRBZ
   integer, target :: LDA, colsA
   integer(c_int) :: err
   integer(c_int64_t), target :: ML64, NL64

      VRBZ=0; if (present(VRBZ_)) VRBZ= int(VRBZ_,c_int)
      err= set_dist_mat_row_col_params_c( &
               int(M, c_int64_t), int(N, c_int64_t), &
               int(MA,c_int64_t), int(NA,c_int64_t), int(rank, c_int), &
               int(nRowDevs, c_int64_t), int(nColDevs, c_int64_t), &
               C_LOC(ML64),  C_LOC(NL64), VRBZ)
      ML = int(ML64)
      NL = int(NL64)
      call cusolverMp_stop_if_error(err, "[!][set_dist_mat_row_col_params], err=")

   end subroutine set_dist_mat_row_col


   type(c_ptr) function set_cusolverMpGrid(numRowDevs, numColDevs, VRBZ_) result(grid)
   use iso_c_binding
   use cuda_comm_aux, only: nccl_handle
   integer, intent(IN)           :: numRowDevs, numColDevs
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= set_cusolverMpGrid_c(cusolvermp_handle, grid, nccl_handle, &
                            int(numRowDevs,c_int64_t), &
                            int(numColDevs,c_int64_t), VRBZ)
      call cusolverMp_stop_if_error(err, "[!][set_cusolverMpGrid], err=")
   end function set_cusolverMpGrid


   subroutine show_cusolverMpGrid(grid)
   use iso_c_binding
   type(c_ptr), intent(IN) :: grid
   !
   integer(c_int) :: c_err
      c_err = show_cusolverMpGrid_c(grid)
      call cusolverMp_stop_if_error(c_err, '[!][show_cusolverMpGrid], err=')
   end subroutine


   subroutine destroy_cusolverMpGrid(grid, VRBZ_)
   use iso_c_binding
   type(c_ptr), intent(IN) :: grid
   integer, optional       :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ= int(VRBZ_,c_int)
      err= destroy_cusolverMpGrid_c(grid, VRBZ)
      call cusolverMp_stop_if_error(err, "[!][destroy_cusolverMpGrid], err=")
   end subroutine destroy_cusolverMpGrid


   type(c_ptr) function  set_cusolverMpMatrixDesc(grid, M, N, MB, NB, &
                         numRowDevs, dev_row_idx, VRBZ_) result(descr)
   use iso_c_binding
   type(c_ptr), intent(IN) :: grid
   integer, intent(IN)     :: M, N, MB, NB
   integer, intent(IN)     :: numRowDevs, dev_row_idx
   integer, optional       :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= set_cusolverMpMatrixDesc_c(descr, grid, &
           int(M, c_int64_t), int(N, c_int64_t), &
           int(MB,c_int64_t), int(NB,c_int64_t), &
           int(numRowDevs, c_int64_t), int(dev_row_idx, c_int), VRBZ)
      call cusolverMp_stop_if_error(err, "[!][set_cusolverMpMatrixDesc], err=")
   end function set_cusolverMpMatrixDesc


   subroutine show_cusolverMpMatrixDesc(descriptor)
   use iso_c_binding
   type(c_ptr), intent(IN) :: descriptor
   !
   integer(c_int) :: err
      err= show_cusolverMpMatrixDescriptor_c(descriptor)
      call cusolverMp_stop_if_error(err, "[!][show_cusolverMpMatrixDesc], err=")
   end subroutine


   subroutine destroy_cusolverMpMatrixDesc(descr, VRBZ_)
   use iso_c_binding
   type(c_ptr), intent(IN) :: descr
   integer, optional       :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= destroy_cusolverMpMatrixDesc_c(descr, VRBZ)
      call cusolverMp_stop_if_error(err, "[!][destroy_cusolverMpMatrixDesc], err=")
   end subroutine destroy_cusolverMpMatrixDesc


   integer function dist_geqrf(descrA, d_A, M, N, d_tau, VRBZ_) result(err)
   use iso_c_binding
   use cuda_comm_aux, only: local_stream
   type(c_ptr), intent(IN) :: descrA
   integer, intent(IN)     :: M, N
   double precision, device, target, intent(INOUT) :: d_A(:), d_tau(:)
   integer, optional       :: VRBZ_
   !
   integer(c_int) :: VRBZ
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= int(dist_geqrf_c(local_stream, cusolvermp_handle, descrA, &
                          C_LOC(d_A), int(M, c_int64_t), int(N, c_int64_t), &
                          C_LOC(d_tau), VRBZ))
   end function dist_geqrf


   integer function check_qr_buffsize(M, N, descrA, d_A) result(err)
   use iso_c_binding
   type(c_ptr), intent(IN) :: descrA
   integer, intent(IN)     :: M,N
   double precision, device, target, intent(INOUT) :: d_A(:)
      err= int(check_qr_buffsize_c(cusolvermp_handle, &
               int(M, c_int64_t), int(N, c_int64_t), &
               int(1, c_int64_t), int(1, c_int64_t), &
               descrA, C_LOC(d_A)))
      call cusolverMp_stop_if_error(err, "[!][check_qr_buffsize], err=")
   end function check_qr_buffsize


   integer function dist_ormqr(descrA, d_A, iside, itrans, M, N, &
                               descrC, d_C, K, d_tau, VRBZ_) result(err)
   use iso_c_binding
   use cuda_comm_aux, only: local_stream
   type(c_ptr), intent(IN) :: descrA, descrC
   integer, intent(IN) :: M, N, K, iside, itrans
   double precision, device, target, intent(IN) :: d_A(:), d_tau(:)
   double precision, device, target, intent(INOUT) :: d_C(:)
   integer, intent(IN), optional :: VRBZ_
   !
   integer :: VRBZ
      VRBZ= 0; if (present(VRBZ_)) VRBZ= VRBZ_
      err= int(dist_ormqr_c(local_stream, cusolvermp_handle, &
                   descrA, C_LOC(d_A), int(iside, c_int), int(itrans, c_int), &
                   int(M, c_int64_t),  int(N, c_int64_t), &
                   descrC, C_LOC(d_C), int(K, c_int64_t), C_LOC(d_tau), int(VRBZ,c_int)))
      call cusolverMp_stop_if_error(err, "[!][dist_ormqr], err=")
   end function


   integer function dist_sytrd(N, descrA, d_A, d_d, d_e, d_tau) result(err)
   use iso_c_binding
   use cuda_comm_aux, only: local_stream
   type(c_ptr), intent(IN) :: descrA
   integer, intent(IN) ::  N
   double precision, device, target, intent(INOUT) :: d_A(:), d_d(:), d_e(:), d_tau(:)
      err= int(dist_sytrd_c(local_stream, cusolvermp_handle, &
                   int(N, c_int64_t), descrA, C_LOC(d_A), &
                   C_LOC(d_d), C_LOC(d_e), C_LOC(d_tau)))
      call cusolverMp_stop_if_error(err, "[!][dist_sytrd], err=")
   end function


   integer function dist_stedc(N, d_d, d_e, descrQ, d_Q, VRBZ_) result(err)
   use iso_c_binding
   type(c_ptr), intent(IN) :: descrQ
   integer, intent(IN) ::  N
   integer, intent(IN), optional :: VRBZ_
   !
   double precision, device, target, intent(INOUT) :: d_d(:), d_e(:), d_Q(:)
   integer :: VRBZ
      VRBZ= 0; if (present(VRBZ_)) VRBZ= VRBZ_
      err= int(dist_stedc_c(cusolvermp_handle, int(N, c_int64_t), &
               C_LOC(d_d), C_LOC(d_e), descrQ, C_LOC(d_Q), int(VRBZ,c_int)))
      call cusolverMp_stop_if_error(err, "[!][dist_stedc], err=")
   end function

end module cusolvermp_aux

