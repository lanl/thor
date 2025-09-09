!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file qtt.f90
!! @author Oleg Korobkin
!! @date   December 2024
!! @brief  Declaration of quantic tensor train types (qtt_tensor and qtt_matrix)
!!
!#define VERBOSE4

module qtt_lib
  use thor_lib
  implicit none
  type, extends(dtt_tensor), public :: qtt_tensor
    integer :: d
  contains
    procedure, pass(self) :: qtt_tensor_assign
    generic :: assignment(=) => qtt_tensor_assign
    final :: dealloc_qtt_tensor
  end type

  type, extends(dtt_matrix), public :: qtt_matrix
    integer :: d
  contains
    procedure, pass(self) :: qtt_matrix_assign
    generic :: assignment(=) => qtt_matrix_assign
    procedure :: full_matrix => qtt_full_matrix
    final :: dealloc_qtt_matrix
  end type qtt_matrix

  !! CONSTRUCTOR\DESTRUCTOR
  interface qtt_tensor ! constructor interface
     procedure :: empty_qtt_tensor
     procedure :: qtt_tensor_from_1darray
     procedure :: qtt_tensor_from_2darray !! TODO: unfinished
  end interface

  interface qtt_matrix ! constructor interface
     procedure :: empty_qtt_matrix
     procedure :: qttm_from_dtt_tensor
     procedure :: qttm_from_square_matrix
  end interface 

  interface kron
     module procedure qtt_kron, qttm_kron
  end interface

  interface core2cell
     module procedure core2cell_qtt, core2cell_qttm
  end interface

  !!TODO !! T_IJK: ACCESS ELEMENT IJK (VALUE)
  !! interface tijk         ! particular element from tt-tensor as function of index
  !!     module procedure dtt_ijk,ztt_ijk
  !! end interface
  !! !! VALUE
  !! interface value        ! value of tt-tensor as function of its coordinates on [0:1]
  !!     module procedure dtt_value,dtt_value0,ztt_value,ztt_value0
  !! end interface

  !! OPERATORS / OVERLOADING
  interface operator (+)
     module procedure qtt_add
  end interface

  interface operator (-)
     module procedure qtt_sub
  end interface


  contains

   !!
   !! TENSOR CONSTRUCTOR \ DESTRUCTOR
   !!
   !> Creates empty qtt-tensor with the given ranks
   !!
   !! Input:
   !!  * d    : number of dimensions
   !!  * r_(:): array of ranks; assumed zero-based [1]
   !!  * l_   : starting index: n(l:) and r(l-1:) [1]
   !!  * m_   : ending index  : n(:m_) and r(:m_) [size(n)]
   !!
   function empty_qtt_tensor(d, r_, l_, m_) result(this)
   type(qtt_tensor):: this
   integer, intent(in) :: d
   integer, intent(in), optional :: r_(0:d)
   integer, intent(in), optional :: l_, m_
   !
   integer, allocatable :: n(:)
#  ifdef VERBOSE4
   print '("[+][empty_qtt_tensor] entry")'
#  endif

      allocate(n(d))
      n(:) = 2
      this= dtt_tensor(n, r_, l_, m_)
      this% d = d
      deallocate(n)

#  ifdef VERBOSE4
   print '("[+][empty_qtt_tensor] exit")'
#  endif
   end function empty_qtt_tensor

   !!
   !! ONES
   !!

   !> Creates a qtt_tensor of of ones of size n = 2^d
   function qtt_ones(d, l_, m_) result(this)
   type(qtt_tensor)    :: this
   integer, intent(in) :: d
   integer, intent(in), optional :: l_, m_
   !
   integer :: k, l, m
#  ifdef VERBOSE4
   print '("[+][qtt_ones] entry")'
#  endif
      l= 1; if (present(l_)) l = l_
      m= d; if (present(m_)) m = m_
      this = empty_qtt_tensor(d, l_=l, m_=m)
      do k=l,m
         this% u(k)% p= 1d0
      end do
#  ifdef VERBOSE4
   print '("[+][qtt_ones] exit")'
#  endif
   end function qtt_ones


   !> from 1D array
   function qtt_tensor_from_1darray(arr, eps_, rmax_) result(this)
   type(qtt_tensor):: this
   double precision, intent(in):: arr(:)
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding
   !
   integer :: d, N, N2
   integer, allocatable :: sh(:)
   double precision, allocatable:: arr2(:)

      N = size(arr)
      d = int(ceiling(log(dble(N))/log(2d0)))
      N2 = 2**d
      allocate(arr2(N2), sh(d))
      arr2(1:N2)= 0d0
      arr2(1:N)= arr(1:N)
      sh(1:d)= 2
      this= dtt_from_darray(arr2, sh, eps_, rmax_)
      this% d= d
      deallocate(sh)
      deallocate(arr2)

   end function qtt_tensor_from_1darray


   !> creates quantic tt-tensor from 2D array
   function qtt_tensor_from_2darray(arr, eps_, rmax_) result(this)
   type(qtt_tensor):: this
   double precision, intent(in):: arr(:,:)
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding
   !
   integer :: d, k
   integer, allocatable :: sh(:), sh2(:)
   double precision, allocatable:: arr2(:,:)
      
   !!! TODO
   !!   sh = shape(arr)
   !!   d = size(sh) ! it's two
   !!   allocate(sh2(d))
   !!   do k=1,d
   !!      sh2(k)= int(log(dble(sh(k)))/log(2d0)) + 1
   !!   enddo
   !!   allocate(arr2(sh2(1),sh2(2)))
   !!   arr2(1:sh(1),1:sh(2))= arr(1:sh(1),1:sh(2))

   !!   deallocate (sh, sh2, arr2)

   !!   !! TODO


   end function qtt_tensor_from_2darray

   !> the qtt_tensor assignment operator, A = B
   subroutine qtt_tensor_assign(self, other)
   use thor_lib, only: dtt_tensor_baseclass_assignment
   class(qtt_tensor), intent(inout) :: self
   type(qtt_tensor), intent(in), target :: other
#  ifdef VERBOSE4
   print '("[+][qtt_tensor_assign] entry")'
#  endif

      call dtt_tensor_baseclass_assignment(self, other)
      self% d = other% d

#  ifdef VERBOSE4
   print '("[+][qtt_tensor_assign] exit")'
#  endif
   end subroutine qtt_tensor_assign


   !> destructor (= deallocator)
   subroutine dealloc_qtt_tensor(this)
   type(qtt_tensor), intent(inout):: this
   integer:: k
#  ifdef VERBOSE4
   print '("[+][dealloc_qtt_tensor] entry")'
#  endif

      do k=this% l, this% m
         if (associated(this% u(k)% p)) deallocate(this% u(k)% p)
         this% u(k)% p=> null()
      enddo
      this% l= 1
      this% m= 0
      this% is_allocated= .false.

#  ifdef VERBOSE4
   print '("[+][dealloc_qtt_tensor] exit")'
#  endif
   end subroutine dealloc_qtt_tensor

!!
!! QTT ARITHMETICS
!!

   function qtt_add(x,y) result(z)
   type(qtt_tensor),intent(in),target :: x,y
   !
   type(qtt_tensor), target :: z
   class(dtt_tensor), pointer :: pt_z
#  ifdef VERBOSE4
   print '("[+][qtt_add] entry")'
#  endif

      pt_z=> z
      call dtt_add_baseclass(pt_z, x, y)
      z% d = x% d

#  ifdef VERBOSE4
   print '("[+][qtt_add] exit")'
#  endif
   end function qtt_add


   function qtt_sub(x,y) result(z)
   type(qtt_tensor),intent(in),target :: x,y
   !
   type(qtt_tensor), target :: z
   class(dtt_tensor), pointer :: pt_z
#  ifdef VERBOSE4
   print '("[+][qtt_sub] entry")'
#  endif

      pt_z=> z
      call dtt_sub_baseclass(pt_z, x, y)
      z% d = x% d

#  ifdef VERBOSE4
   print '("[+][qtt_sub] exit")'
#  endif
   end function qtt_sub

!!
!! MATRIX CONSTRUCTOR \ DESTRUCTOR
!!

   !> creates empty qtt matrix 2^d x 2^d with some ranks (optional)
   function empty_qtt_matrix(d, r_, l_, m_) result(this)
   type(qtt_matrix):: this
   integer, intent(in) :: d
   integer, intent(in), optional :: r_(0:d)
   integer, intent(in), optional :: l_, m_
   !
   integer, allocatable :: q(:), s(:)
#  ifdef VERBOSE4
   print '("[+][empty_qtt_matrix] entry")'
#  endif

      allocate(q(d),s(d))
      q(:) = 2
      s(:) = 2
      this= empty_dtt_matrix(q, s, r_, l_, m_)
      this% d = d
      deallocate(q)
      deallocate(s)

#  ifdef VERBOSE4
   print '("[+][empty_qtt_matrix] exit")'
#  endif
   end function empty_qtt_matrix


   !> deep copy constructor
   function qttm_from_dtt_tensor(dtt) result(this)
   class(dtt_tensor), intent(in):: dtt
   type(qtt_matrix):: this
   !
   integer:: d, l, m, k, sh(3)
#  ifdef VERBOSE4
   print '("[+][qttm_from_dtt_tensor] entry")'
#  endif

      l= dtt% l
      m= dtt% m
      d= m - l + 1
      if (any(dtt% n(l:m).ne.4)) error stop "[!][qttm_from_dtt_tensor]: n/=4"
      this= empty_qtt_matrix(d, dtt% r(l-1:m), l_=l, m_=m) 
      do k=l,m
         !sh= shape(dtt% u(k)% p)
         !call realloc_to_shape(this% u(k), sh)
         this% u(k)% p(:,:,:)= dtt% u(k)% p(:,:,:)
         !this% u4(k)% p(1:this% r(k-1), 1:2, 1:2,  &
         !               1:this% r(k))=> this% u(k)% p
      enddo
      this% is_allocated= .true.
      this% nu_core= dtt% nu_core

#  ifdef VERBOSE4
   print '("[+][qttm_from_dtt_tensor] exit")'
#  endif
   end function qttm_from_dtt_tensor


   !> construction from 2^d x 2^d matrix
   function qttm_from_square_matrix(A, eps_) result(this)
   use matrix_util, only: interlace_2d
   double precision, intent(IN):: A(:,:)
   type(qtt_matrix):: this
   double precision, optional:: eps_
   !
   double precision, allocatable:: q(:)
   integer:: d, l, m, k, sh(2)
   integer, allocatable:: nn(:)
   type(dtt_tensor):: tt
#  ifdef VERBOSE4
   print '("[+][qttm_from_square_matrix] entry")'
#  endif

      sh= shape(A)
      if (sh(1) /= sh(2)) &
         error stop "[!][qttm_from_square_matrix]: matrix not square!"

      q= interlace_2d(A)
      d = int(log(dble(sh(1)))/log(2d0) + 1d-11)
      allocate(nn(d)); nn= 4
      this= qttm_from_dtt_tensor(dtt_tensor(q, nn, eps_))
      !this= qttm_from_dtt_tensor(tt)
      deallocate(q, nn)

#  ifdef VERBOSE4
   print '("[+][qttm_from_square_matrix] exit")'
#  endif
   end function qttm_from_square_matrix


   !> creates the identity matrix 2^d x 2^d in QTT form
   function qtt_eye(d) result(this)
   type(qtt_matrix):: this
   integer, intent(in):: d
   !
   integer :: k
   type(qtt_matrix):: I2
#  ifdef VERBOSE4
   print '("[+][qtt_eye] entry")'
#  endif

      this= empty_qtt_matrix(1)
      this% u(1)% p(1,:,1)= [1d0, 0d0, 0d0, 1d0] 
      if (d.gt.1) then
         I2= this
         do k=2,d
            this= kron(this, I2)
         enddo
      endif

#  ifdef VERBOSE4
   print '("[+][qtt_eye] exit")'
#  endif
   end function qtt_eye


   !> the qtt_matrix assignment operator, A = B
   subroutine qtt_matrix_assign(self, other)
   use thor_lib, only: dtt_tensor_baseclass_assignment
   class(qtt_matrix), intent(inout) :: self
   type(qtt_matrix), intent(in), target :: other
   !
   integer :: k, k1, l, m, sh(4)
#  ifdef VERBOSE4
   print '("[+][qtt_matrix_assign] entry")'
#  endif

      l= other% l
      m= other% m
      call dtt_tensor_baseclass_assignment(self, other)
      self% d = other% d
      self% q(self% l:self% m)= other% q(l:m)
      self% s(self% l:self% m)= other% s(l:m)

      do k=l,m
         k1= k - l + self% l
         sh= shape(other% u4(k)% p)
         self% u4(k1)% p => null()
         call realloc_to_shape(self% u4(k1), sh)
         self% u4(k1)% p(1:self% r(k1-1), &
                         1:self% q(k1),   &
                         1:self% s(k1),   &
                         1:self% r(k1))=> self% u(k1)% p
      enddo

#  ifdef VERBOSE4
   print '("[+][qtt_matrix_assign] exit")'
#  endif
   end subroutine qtt_matrix_assign


   !!
   !! FULL: THOR-VECT -> ARRAY
   !!
   !> returns a full matrix
   function qtt_full_matrix(this) result(a)
   use matrix_util, only: binary_repr, integer_from_binary
   implicit none
   class(qtt_matrix), intent(in) :: this
   double precision, allocatable :: a(:,:)
   !
   integer :: l, i, j, d, ij
   integer, allocatable :: ind1(:), ind2(:), ind12(:)
   double precision, allocatable :: arr(:)
#  ifdef VERBOSE4
   print '("[+][qtt_full_matrix] entry")'
#  endif

      d= this% d
      l = 2**d
      allocate(a(l,l),ind1(d),ind2(d),ind12(2*d))

      arr = dtt_full(this)
      ind1(:)= 0
      ind2(:)= 0
      do i = 1,l
         ind1= binary_repr(l+i-1)
         ind12(1::2) = ind1(1:d)
         do j = 1,l
            ind2= binary_repr(l+j-1)
            ind12(2::2) = ind2(1:d)
            ij = integer_from_binary(ind12) + 1
            a(i,j)= arr(ij)
         enddo
      enddo

      deallocate(arr)
      deallocate(ind1)
      deallocate(ind2)
      deallocate(ind12)

#  ifdef VERBOSE4
   print '("[+][qtt_full_matrix] exit")'
#  endif
   end function qtt_full_matrix


   !> Kronecker product of two qtt-tensors
   !!
   !! a, b : qtt_tensor
   !!     operands
   !!
   function qtt_kron(a, b) result(c)
   type(qtt_tensor), intent(in) :: a, b
   type(qtt_tensor) :: c
   !
   integer :: d1, d2, d, k
   integer, allocatable :: r(:)

      d1 = a%m - a%l + 1
      d2 = b%m - b%l + 1
      d = d1 + d2
      allocate(r(0:d))

      r(0:d1)   = a% r(0:d1)
      r(d1:d)   = b% r(0:d2)

      c = empty_qtt_tensor(d, r_=r)
      do k=1,d1
         c% u(k)% p= a% u(k)% p
      enddo
      do k=d1+1,d1+d2
         c% u(k)% p= b% u(k-d1)% p
      enddo
    
   end function qtt_kron


   !> Kronecker product of two dtt matrices
   !!
   !!   A, B : dtt_tensor
   !!       operands
   !!
   function qttm_kron(A, B) result(C)
   type(qtt_matrix), intent(in) :: A, B
   type(qtt_matrix) :: C
   !
   integer :: d1, d2, d, k
   integer, allocatable :: r(:)

      d1 = A%m - A%l + 1
      d2 = B%m - B%l + 1
      d = d1 + d2
      allocate(r(0:d))

      r(0:d1)   = A% r(0:d1)
      r(d1:d)   = B% r(0:d2)

      c = empty_qtt_matrix(d, r_=r)
      do k=1,d1
         c% u4(k)% p= a% u4(k)% p
      enddo
      do k=d1+1,d1+d2
         c% u4(k)% p= b% u4(k-d1)% p
      enddo
    
   end function qttm_kron


   !> matrix destructor (= deallocator)
   subroutine dealloc_qtt_matrix(this)
   type(qtt_matrix), intent(inout):: this
   integer:: k
#  ifdef VERBOSE4
   print '("[+][dealloc_qtt_matrix] entry")'
#  endif

      do k=this% l, this% m
         if (associated(this% u4(k)% p)) deallocate(this% u4(k)% p)
         this% u4(k)% p=> null()
         this% u(k)% p=> null()
      enddo
      this% l= 1
      this% m= 0
      this% is_allocated= .false.

#  ifdef VERBOSE4
   print '("[+][dealloc_qtt_matrix] exit")'
#  endif
   end subroutine dealloc_qtt_matrix


   !!
   !! CORE2CELL FUNCTIONS
   !!
   function core2cell_qtt(tt) result(cc)
   use matlab_struct_module, only: array3d
   implicit none
   type(qtt_tensor), intent(in):: tt
   type(array3d) :: cc(tt% m - tt% l + 1)
   !
   integer :: i, d

       d = tt% m - tt% l + 1 ! number of dimensions of tt
       do i=1,d
          cc(i)% arr = tt% u(i)% p
       enddo

   end function core2cell_qtt


   function core2cell_qttm(tt) result(cc)
   use matlab_struct_module, only: array4d
   implicit none
   type(qtt_matrix), intent(in):: tt
   type(array4d) :: cc(tt% m - tt% l + 1)
   !
   integer :: i, d

       d = tt% m - tt% l + 1 ! number of dimensions of tt
       do i=1,d
          cc(i)% arr = tt% u4(i)% p
       enddo

   end function core2cell_qttm


   !!
   !! HELPER FUNCTIONS
   !!
end module qtt_lib

