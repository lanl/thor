!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file cublas_aux.f90
!! @author Ismael Djibrilla, Oleg Korobkin
!! @date  May 2025
!! @brief Fortran interface to CUDA cuBLAS library wrappers (see cuda_cublas_aux.c)
!!
module cublas_aux
use iso_c_binding
implicit none

type(c_ptr) :: cublas_handle

interface
   ! C PROTOTYPE: int set_cublasHandle(cublasHandle_t *handle, int VRBZ)
   integer(c_int) function set_cublasHandle(handle, VRBZ) &
              bind(C,name="set_cublasHandle")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: VRBZ
   end function


   ! C PROTOTYPE: int destroy_cublasMpHandle(cublasMpHandle_t *handle, int VRBZ);
   integer(c_int) function destroy_cublasHandle(handle, VRBZ) &
              bind(C,name="destroy_cublasHandle")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: VRBZ
   end function


   ! C PROTOTYPE: int compute_cublasIdamax(cublasHandle_t *handle, int n,
   !                                       const double *x, int incx, int VRBZ)
   integer(c_int) function compute_cublasIdamax(handle, n, d_X, incx, VRBZ) &
              bind(C,name="compute_cublasIdamax")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incx, VRBZ
     type(c_ptr), value    :: d_X
   end function


   ! C PROTOTYPE: int compute_cublasIdamin(cublasHandle_t *handle, int n,
   !                                       const double *x, int incx, int VRBZ)
   integer(c_int) function compute_cublasIdamin(handle, n, d_X, incx, VRBZ) &
              bind(C,name="compute_cublasIdamin")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incx, VRBZ
     type(c_ptr), value    :: d_X
   end function


   ! C PROTOTYPE: double compute_cublasDasum(cublasHandle_t *handle, int n,
   !                                         const double *x, int incx, int VRBZ)
   real(c_double) function compute_cublasDasum(handle, n, d_X, incx, VRBZ) &
              bind(C,name="compute_cublasDasum")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incx, VRBZ
     type(c_ptr), value    :: d_X
   end function


   ! C PROTOTYPE: int compute_cublasDaxpy(cublasHandle_t *handle, int n,
   !                                      const double *alpha,
   !                                      const double *x, int incx,
   !                                      double       *y, int incy,
   !                                      int VRBZ)
   integer(c_int) function compute_cublasDaxpy(handle, n, &
              alpha, d_X, incX, d_Y, incY, VRBZ) &
              bind(C,name="compute_cublasDaxpy")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incX, incY, VRBZ
     type(c_ptr), value    :: d_X, d_Y
     real(c_double), value :: alpha
   end function


   ! C PROTOTYPE: int compute_cublasDcopy(cublasHandle_t *handle, int n,
   !                                      const double *x, int incx,
   !                                      double       *y, int incy,
   !                                      int VRBZ)
   integer(c_int) function compute_cublasDcopy(handle, n, d_X, incX, d_Y, incY, VRBZ) &
              bind(C,name="compute_cublasDcopy")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incX, incY, VRBZ
     type(c_ptr), value    :: d_X, d_Y
   end function


   ! C PROTOTYPE:  double compute_cublasDdot(cublasHandle_t *handle, int n,
   !                                         const double *x, int incx,
   !                                         const double *y, int incy, int VRBZ)
   real(c_double) function compute_cublasDdot(handle, n, d_X, incX, d_Y, incY, VRBZ) &
              bind(C,name="compute_cublasDdot")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incX, incY, VRBZ
     type(c_ptr), value    :: d_X, d_Y
   end function


   ! C PROTOTYPE:  double compute_cublasDnrm2(cublasHandle_t *handle, int n,
   !                                          const double *x, int incx, int VRBZ)
   real(c_double) function compute_cublasDnrm2(handle, n, d_X, incx, VRBZ) &
              bind(C,name="compute_cublasDnrm2")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incx, VRBZ
     type(c_ptr), value    :: d_X
   end function


   ! C PROTOTYPE:  double compute_cublasDrot(cublasHandle_t *handle, int n,
   !                                         const double *x, int incx,
   !                                         const double *y, int incy,
   !                                         const double *c, const double *s,
   !                                         int VRBZ)
   integer(c_int) function compute_cublasDrot(handle, n, d_X, incX, d_Y, incY, VRBZ) &
              bind(C,name="compute_cublasDrot")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incX, incY, VRBZ
     type(c_ptr), value    :: d_X, d_Y
     real(c_double), value :: c, s
   end function


   ! C PROTOTYPE: int compute_cublasDscal(cublasHandle_t *handle, int n,
   !                                      const double *alpha,
   !                                      double *x, int incx,
   !                                      int VRBZ)
   integer(c_int) function compute_cublasDscal(handle, n, alpha, d_X, incX, VRBZ) &
              bind(C,name="compute_cublasDscal")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incX, VRBZ
     type(c_ptr), value    :: d_X
     real(c_double), value :: alpha
   end function


   ! C PROTOTYPE: int compute_cublasDswap(cublasHandle_t *handle, int n,
   !                                      double *x, int incx,
   !                                      double *y, int incy,
   !                                      int VRBZ)
   integer(c_int) function compute_cublasDswap(handle, n, d_X, incX, d_Y, incY, VRBZ) &
              bind(C,name="compute_cublasDswap")
   use iso_c_binding
     type(c_ptr)           :: handle
     integer(c_int), value :: n, incX, incY, VRBZ
     type(c_ptr), value    :: d_X, d_Y
   end function

end interface

contains

   !> cuBLAS library initialization routine: setup handle, etc.
   subroutine setup_cublas(VRBZ_)
   use iso_c_binding
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= set_cublasHandle(cublas_handle, VRBZ)
      call cublas_stop_if_error(err, "[!][setup_cublas], err=")
   end subroutine


   !> Library cleanup routine: destroy handle, etc.
   subroutine cleanup_cublas(VRBZ_)
   use iso_c_binding
   integer, intent(IN), optional :: VRBZ_
   !
   integer(c_int) :: VRBZ, err
      VRBZ=0; if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      err= destroy_cublasHandle(cublas_handle, VRBZ)
      call cublas_stop_if_error(err, "[!][cleanup_cublas], err=")
   end subroutine


   !> Returns index of the maximum element in an array (base 1)
   integer function cublas_iamax(d_X, incx_, VRBZ_) result(idx)
   use iso_c_binding
   double precision, device, intent(IN) :: d_X(:)
   integer, intent(IN), optional        :: incx_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, N
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      N= int(size(d_X),c_int)
      idx= int(compute_cublasIdamax(cublas_handle, N, c_loc(d_X), incx, VRBZ))
   end function


   !> Returns index of the minimum element in an array (base 1)
   integer function cublas_iamin(d_X, incx_, VRBZ_) result(idx)
   use iso_c_binding
   double precision, device, intent(IN) :: d_X(:)
   integer, intent(IN), optional        :: incx_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, N
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      N= int(size(d_X),c_int)
      idx= int(compute_cublasIdamin(cublas_handle, N, c_loc(d_X), incx, VRBZ))
   end function


   !> Computes the sum of the absolute values of the elements of vector
   double precision function cublas_asum(d_X, incx_, VRBZ_) result(asum)
   use iso_c_binding
   double precision, device, intent(IN) :: d_X(:)
   integer, intent(IN), optional        :: incx_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, N
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      N= int(size(d_X),c_int)
      asum= compute_cublasDasum(cublas_handle, N, c_loc(d_X), incx, VRBZ)
   end function


   !> Computes y <- a*x + y
   subroutine cublas_axpy(alpha, d_X, d_Y, incx_, incy_, VRBZ_)
   use iso_c_binding
   double precision, intent(IN)            :: alpha
   double precision, device, intent(IN)    :: d_X(:)
   double precision, device, intent(INOUT) :: d_Y(:)
   integer, intent(IN), optional           :: incx_, incy_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, incy, N, err
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      incy=one; if (present(incy_)) incy = int(incy_,c_int)
      N= int(size(d_X),c_int)
      err= compute_cublasDaxpy(cublas_handle, N, alpha, c_loc(d_X), incx, c_loc(d_Y), incy, VRBZ)
      call cublas_stop_if_error(err, "[!][cublas_axpy], err=")
   end subroutine


   !> Copies x-> y
   subroutine cublas_copy(d_X, d_Y, incx_, incy_, VRBZ_)
   use iso_c_binding
   double precision, device, intent(IN)    :: d_X(:)
   double precision, device, intent(INOUT) :: d_Y(:)
   integer, intent(IN), optional           :: incx_, incy_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, incy, N, err
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      incy=one; if (present(incy_)) incy = int(incy_,c_int)
      N= int(size(d_X),c_int)
      err= compute_cublasDcopy(cublas_handle, N, c_loc(d_X), incx, c_loc(d_Y), incy, VRBZ)
      call cublas_stop_if_error(err, "[!][cublas_copy], err=")
   end subroutine


   !> Copies 2D (sub)matrix x-> y
   !!
   !! Inputs:
   !! ------
   !!    d_X2, d_Y2 : source & destination arrays, resp.
   !!    M_, N_     : (optional) size of the 2D matrix to copy
   !!    xi0_, xj0_ : (optional) top-left corner point of the copied block in d_X2
   !!    yi0_, yj0_ : (optional) top-left corner point of where to copy in d_Y2
   !!
   subroutine cublas_copy_2d(d_X2, d_Y2, M_, N_, xi0_, xj0_, yi0_, yj0_, VRBZ_)
   use iso_c_binding
   double precision, device, intent(IN)    :: d_X2(:,:)
   double precision, device, intent(INOUT) :: d_Y2(:,:)
   integer, intent(IN), optional :: M_, N_, xi0_, xj0_, yi0_, yj0_, VRBZ_
   !
   integer(c_int) :: VRBZ, M, N, xi0, xj0, yi0, yj0, incr, err
   integer(c_int), parameter :: one = int(1, c_int)
   integer :: j
      incr= one
      VRBZ=0;  if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      xi0=one; if (present(xi0_)) xi0 = int(xi0_,c_int)
      xj0=one; if (present(xj0_)) xj0 = int(xj0_,c_int)
      yi0=one; if (present(yi0_)) yi0 = int(yi0_,c_int)
      yj0=one; if (present(yj0_)) yj0 = int(yj0_,c_int)
      M= int(size(d_X2, 1),c_int); if (present(M_)) M = int(M_,c_int)
      N= int(size(d_X2, 2),c_int); if (present(N_)) N = int(N_,c_int)
      do j=1,int(N)
         ! d_Y2(yi0:yi0+M-1,yj0+j-1)= d_X2(xi0:xi0+M-1,xj0+j-1)
         err= compute_cublasDcopy(cublas_handle, M, c_loc(d_X2(xi0,xj0+j-1)), incr, &
                                                    c_loc(d_Y2(yi0,yj0+j-1)), incr, VRBZ)
         call cublas_stop_if_error(err, "[!][cublas_copy], err=")
      enddo
   end subroutine


   !> Scalar product (x, y)
   double precision function cublas_dot(d_X, d_Y, incx_, incy_, VRBZ_) result(xy)
   use iso_c_binding
   double precision, device, intent(IN) :: d_X(:)
   double precision, device, intent(IN) :: d_Y(:)
   integer, intent(IN), optional        :: incx_, incy_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, incy, N
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      incy=one; if (present(incy_)) incy = int(incy_,c_int)
      N= int(size(d_X),c_int)
      xy= compute_cublasDdot(cublas_handle, N, c_loc(d_X), incx, c_loc(d_Y), incy, VRBZ)
   end function


   !> Computes maximum value in an array
   double precision function cublas_maxval(d_X, incx_, VRBZ_) result(mval)
   double precision, device, intent(IN) :: d_X(:)
   integer, intent(IN), optional        :: incx_, VRBZ_
   !
   integer :: idx
      idx= cublas_iamax(d_X, incx_, VRBZ_)
      mval= d_X(idx)
   end function


   !> Computes minimum value in an array
   double precision function cublas_minval(d_X, incx_, VRBZ_) result(mval)
   double precision, device, intent(IN) :: d_X(:)
   integer, intent(IN), optional        :: incx_, VRBZ_
   !
   integer :: idx
      idx= cublas_iamin(d_X, incx_, VRBZ_)
      mval= d_X(idx)
   end function


   !> Computes the Euclidean norm of the vector x (= sqrt(sum[|x_i|^2]))
   double precision function cublas_nrm2(d_X, incx_, VRBZ_) result(nrm2)
   use iso_c_binding
   double precision, device, intent(IN) :: d_X(:)
   integer, intent(IN), optional        :: incx_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, N
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      N= int(size(d_X),c_int)
      nrm2= compute_cublasDnrm2(cublas_handle, N, c_loc(d_X), incx, VRBZ)
   end function


   !> Scales a vector by the specified factor
   subroutine cublas_scal(alpha, d_X, incx_, VRBZ_)
   use iso_c_binding
   double precision, intent(IN)            :: alpha
   double precision, device, intent(INOUT) :: d_X(:)
   integer, intent(IN), optional           :: incx_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, N, err
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      N= int(size(d_X),c_int)
      err= compute_cublasDscal(cublas_handle, N, alpha, c_loc(d_X), incx, VRBZ)
      call cublas_stop_if_error(err, "[!][cublas_axpy], err=")
   end subroutine


   !> Swaps two vectors x <-> y
   subroutine cublas_swap(d_X, d_Y, incx_, incy_, VRBZ_)
   use iso_c_binding
   double precision, device, intent(INOUT) :: d_X(:)
   double precision, device, intent(INOUT) :: d_Y(:)
   integer, intent(IN), optional           :: incx_, incy_, VRBZ_
   !
   integer(c_int) :: VRBZ, incx, incy, N, err
   integer(c_int), parameter :: one = int(1, c_int)
      VRBZ=0;   if (present(VRBZ_)) VRBZ = int(VRBZ_,c_int)
      incx=one; if (present(incx_)) incx = int(incx_,c_int)
      incy=one; if (present(incy_)) incy = int(incy_,c_int)
      N= int(size(d_X),c_int)
      err= compute_cublasDswap(cublas_handle, N, c_loc(d_X), incx, c_loc(d_Y), incy, VRBZ)
      call cublas_stop_if_error(err, "[!][cublas_swap], err=")
   end subroutine


   !> local error handler
   subroutine cublas_stop_if_error(err, msg)
   use, intrinsic :: iso_fortran_env, only: error_unit
   use iso_c_binding
   integer(c_int), intent(in) :: err
   character(*), intent(in)   :: msg
      if (err.ne.0) then
         write(error_unit,'(A,I8)') msg, err
         error stop
      endif
   end subroutine

end module cublas_aux
