!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file thor.f90
!! @author Ismael Djibrilla Boureima, Oleg Korobkin
!! @date   April 2023
!! @brief  Declaration of dtt_tensor and dtt_matrix types
!!
!#define VERBOSE4

module thor_lib
  use string_lib, only: str
  use mat_lib, only: thor_matmul
  use thor_pointers, only: pointd,pointz, pointd3, pointz3, pointd4, pointz4
  implicit none
  integer,parameter :: tt_size=2048
  type, public  :: dtt_tensor        ! Double precision tensor train tt-dtt_tensor
    integer        :: l = 1          ! index of the leftmost core
    integer        :: m = 0          ! index of the rightmost core (dimension)
    integer        :: n(tt_size)=0   ! mode sizes
    integer        :: r(0:tt_size)=0 ! tensor ranks
    type(pointd3)  :: u(tt_size)     ! Host/CPU tensor cores
    logical        :: is_allocated=.false. ! Allocation status
    logical        :: is_rhs=.false. ! is a temporary, "righ-hand side" object?
    integer        :: nu_core= 0     ! the non-unitary core; zero before SVD
  contains
    procedure, pass(self) :: dtt_tensor_assign
    generic :: assignment(=) => dtt_tensor_assign
    procedure      :: full => dtt_full
    procedure      :: svd => dtt_svd
    procedure      :: mround => matlab_round
    procedure      :: rounding => dtt_svd
    procedure      :: round => dtt_round
    procedure      :: size => dtt_size
    procedure      :: log10_fullsize => dtt_log10_fullsize
    procedure      :: fullsize_string => dtt_fullsize_string
    procedure      :: norm => dtt_norm
    procedure      :: normb => dtt_bounded_norm
    procedure      :: grand_sum => dtt_grand_sum
    procedure      :: flat_cores => dtt_flat_cores
    procedure      :: say => dtt_say
    procedure      :: pprint => dtt_pprint_aux
    final          :: dealloc_dtt_tensor
  end type

  type, extends(dtt_tensor), public :: dtt_matrix  ! Double precision tt-matrix
    integer        :: q(tt_size)=0    ! first mode sizes (for matrices)
    integer        :: s(tt_size)=0    ! second mode sizes (for matrices)
    type(pointd4)  :: u4(tt_size)     ! Host/CPU cores
  contains
    procedure, pass(self) :: dtt_matrix_assign
    generic :: assignment(=) => dtt_matrix_assign
    procedure      :: full2d => dttm_full
    procedure      :: say => dttm_say
    procedure      :: pprint => dttm_pprint_aux
    final :: dealloc_dtt_matrix
  end type dtt_matrix

  !! CONSTRUCTOR\DESTRUCTOR
  interface dtt_tensor ! constructor interface
     procedure :: empty_dtt_tensor_int4
     procedure :: empty_dtt_tensor_int8
     procedure :: dtt_from_dtt_tensor
     procedure :: dtt_from_cell_array
     procedure :: dtt_from_dttm
     procedure :: dtt_from_darray
     procedure :: dtt_from_flat_corearr
     procedure :: dtt_from_2darray
     procedure :: dtt_from_3darray
     procedure :: dtt_from_4darray
     procedure :: dtt_from_5darray
     procedure :: dtt_from_6darray
     procedure :: dtt_from_7darray
  end interface dtt_tensor
  interface dtt_matrix ! constructor interface
     procedure :: empty_dtt_matrix
     procedure :: dttm_from_dtt_matrix
     procedure :: dttm_from_dtt_tensor
     procedure :: dttm_from_2dmatrix
     procedure :: dttm_from_cell_array
     procedure :: dttm_from_flat_corearr
  end interface dtt_matrix

  !! ALLOC\DEALLOC\READY
  interface alloc
      module procedure alloc_d
  end interface
  interface dealloc
      module procedure dealloc_d
  end interface
  interface ready
      module procedure dtt_ready
  end interface
  !! REALLOC
  interface realloc_to_shape ! reallocate to size
     module procedure realloc_pointd3, realloc_pointd4, &
                      realloc_pointz3, realloc_pointz4
  end interface realloc_to_shape
  !! MEMORY
  interface memory          ! memory to keep all cores
      module procedure dtt_mem
  end interface
  !! SAY
  interface sayfull
      module procedure dtt_sayfull
  end interface
  !! PRETTY PTINT
  interface pprint
      module procedure dtt_pprint, dttm_pprint
  end interface
  !! RANK
  interface erank        ! effective rank of thor tensor
      module procedure dtt_rank
  end interface
  !! ORT
  interface ort          ! tt-orthogonalization
      module procedure dtt_ort
  end interface
  !! DOT
  interface dot_prod
      module procedure dtt_dot
  end interface
  !! T_IJK: ACCESS ELEMENT IJK (VALUE)
  interface tijk         ! particular element from tt-tensor as function of index
      module procedure dtt_ijk
  end interface
  !! RANDOM ORTHOGONAL TENSORS
  interface dtt_random_ortho
      module procedure dtt_random_ortho_full, dtt_random_ortho_short
  end interface
  !! OPERATORS / OVERLOADING
  interface operator (+)
     module procedure dtt_add, dtt_plus_d, d_plus_dtt, dttm_add
  end interface
  interface operator (-)
     module procedure dtt_sub, dtt_minus_d, d_minus_dtt, minus_dtt, &
                      dttm_sub
  end interface
  interface operator (*)
      module procedure dtt_mul_d, d_mul_dtt, dtt_mult
  end interface
  interface operator (/)
      module procedure dtt_div_d
  end interface
  interface operator (**)
      module procedure dtt_pwr_i
  end interface
  interface core2cell
     module procedure core2cell_dtt, core2cell_dttm
  end interface
  !! KRONECKER PRODUCT
  interface kron
      module procedure dtt_kron, dttm_kron
  end interface


  contains


   !!
   !! TENSOR CONSTRUCTOR \ DESTRUCTOR
   !!
   !> Creates empty tensor from arrays of modes and ranks
   !!
   !! Input:
   !!  * n(:) : array of modes
   !!  * r_(:): array of ranks; assumed zero-based [1]
   !!  * l_   : starting index: n(l:) and r(l-1:) [1]
   !!  * m_   : ending index  : n(:m_) and r(:m_) [size(n)]
   !!
   function empty_dtt_tensor_int4(n, r_, l_, m_) result(this)
   type(dtt_tensor):: this
   integer(kind=4), intent(in) :: n(:)
   integer(kind=4), intent(in), optional :: r_(0:size(n))
   integer(kind=4), intent(in), optional :: l_, m_
   !
   integer(kind=4) d, m, k, sh(3), l
   integer(kind=4), allocatable :: r(:)
#  ifdef VERBOSE4
   print '("[+][empty_dtt_tensor_int4] entry")'
#  endif
      l= 1; if (present(l_)) l= l_
      m= size(n); if (present(m_)) m= m_
      this% l = l; this% m = m
      d= m - l + 1
      if (l.gt.size(n)) error stop "in empty_dtt_tensor_int4: l > |n|"
      if (m.gt.size(n)) error stop "in empty_dtt_tensor_int4: m > |n|"
      if (d.lt.0) error stop "in empty_dtt_tensor_int4: m + 1 < l"
      allocate(r(0:size(n))); r = 1
      if (present(r_)) r(l-1:m) = r_(l-1:m)
      if(any(n(l:m).le.0)) then
         write(*,*) 'empty_dtt_tensor_int4: bad n: ',n(l:m)
         error stop
      endif
      if(any(r(l:m).le.0)) then
         write(*,*) 'empty_dtt_tensor_int4: bad r: ',r(l:m)
         error stop
      endif
      this% n(l:m)= n(l:m)
      this% r(l-1:m)= 1
      if (present(r_)) this% r(l-1:m)= r_(l-1:m)
      do k=l,m
         sh= [this% r(k-1), this% n(k), this% r(k)]
         call realloc_to_shape(this% u(k), int(sh))
         this% u(k)% p= 0d0
      enddo
      this% nu_core= 0
      this% is_allocated= .true.
      deallocate(r)
#  ifdef VERBOSE4
   print '("[+][empty_dtt_tensor_int4] exit")'
#  endif
   end function empty_dtt_tensor_int4


   function empty_dtt_tensor_int8(n, r_, l_, m_) result(this)
   type(dtt_tensor):: this
   integer(kind=8), intent(in) :: n(:)
   integer(kind=8), intent(in), optional :: r_(0:size(n))
   integer(kind=8), intent(in), optional :: l_, m_
   !
   integer(kind=8) :: d, m, k, sh(3), l
   integer(kind=8), allocatable :: r(:)
#  ifdef VERBOSE4
   print '("[+][empty_dtt_tensor_int8] entry")'
#  endif
      l= 1; if (present(l_)) l= l_
      m= size(n); if (present(m_)) m= m_
      this% l = l; this% m = m
      d= m - l + 1
      if (l.gt.size(n)) error stop "in empty_dtt_tensor_int8: l > |n|"
      if (m.gt.size(n)) error stop "in empty_dtt_tensor_int8: m > |n|"
      if (d.lt.0) error stop "in empty_dtt_tensor_int8: m + 1 < l"
      allocate(r(0:size(n))); r = 1
      if (present(r_)) r(l-1:m) = r_(l-1:m)
      if(any(n(l:m).le.0)) then
         write(*,*) 'empty_dtt_tensor_int8: bad n: ',n(l:m)
         error stop
      endif
      if(any(r(l:m).le.0)) then
         write(*,*) 'empty_dtt_tensor_int8: bad r: ',r(l:m)
         error stop
      endif
      this% n(l:m)= n(l:m)
      this% r(l-1:m)= 1
      if (present(r_)) this% r(l-1:m)= r_(l-1:m)
      do k=l,m
         sh= [this% r(k-1), this% n(k), this% r(k)]
         call realloc_to_shape(this% u(k), int(sh))
         this% u(k)% p= 0d0
      enddo
      this% nu_core= 0
      this% is_allocated= .true.
      deallocate(r)
#  ifdef VERBOSE4
   print '("[+][empty_dtt_tensor_int8] exit")'
#  endif
   end function empty_dtt_tensor_int8


   !> deep copy constructor
   function dtt_from_dtt_tensor(other) result(this)
   type(dtt_tensor):: this
   type(dtt_tensor), intent(in):: other
   integer:: l, m, k, sh(3)

      call dealloc(this)

      l= other% l
      m= other% m
      this% l= l
      this% m= m
      this% n(l:m)= other% n(l:m)
      this% r(l-1:m)= other% r(l-1:m)
      do k=l,m
         sh= shape(other% u(k)% p)
         call realloc_to_shape(this% u(k), sh)
         this% u(k)% p(:,:,:)= other% u(k)% p(:,:,:)
      enddo
      this% nu_core= other% nu_core
      this% is_allocated= .true.

   end function dtt_from_dtt_tensor


   !> constructor from a cell array
   function dtt_from_cell_array(cells) result(this)
   use matlab_struct_module
   type(dtt_tensor):: this
   type(array3d), intent(in):: cells(:)
   integer:: l, m, d, k, sh(3)

      call dealloc(this)

      d= size(cells)
      m= d
      this% l= 1
      this% m= m
      this% r(0)= 1
      do k=1,m
         sh= shape(cells(k)% arr)
         this% n(k)= sh(2)
         this% r(k)= sh(3)
         call realloc_to_shape(this% u(k), sh)
         this% u(k)% p(:,:,:)= cells(k)% arr(:,:,:)
      enddo
      this% nu_core= 0 ! do not assume cores are orthogonzalized
      this% is_allocated= .true.

   end function dtt_from_cell_array


   !> converts dtt_matrix into dtt_tensor
   function dtt_from_dttm(this) result(dtt)
   type(dtt_matrix), intent(in):: this
   type(dtt_tensor):: dtt
   integer k, l, m, s
#  ifdef VERBOSE4
   print '("[+][dtt_from_dttm] entry")'
#  endif

      l= this% l
      m= this% m
      dtt= dtt_tensor(this% n(l:m), l_=l, r_=this% r(l-1:m))
      do k=l,m
         s= this% r(k-1)*this% n(k)*this% r(k)
         call dcopy(s, this% u(k)% p, 1, dtt% u(k)% p, 1)
      enddo
      dtt% nu_core= this% nu_core

#  ifdef VERBOSE4
   print '("[+][dtt_from_dttm] exit")'
#  endif
   end function dtt_from_dttm


   !> create dtt_tensor from an assumed-shape array of doubles
   function dtt_from_darray(a,n,eps_,rmax_) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: a(*)  !< array to convert to tt-format
   integer,intent(in)          :: n(:)  !< modes of a (typically from shape(a))
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding
#  ifdef VERBOSE4
   print '("[+][dtt_from_darray] entry")'
#  endif

      call dtt_svd0(n, a, tt, eps_, rmax_)

#  ifdef VERBOSE4
   print '("[+][dtt_from_darray] exit")'
#  endif
   end function dtt_from_darray


   !> constructor from a flattened array of the sequence of cores
   function dtt_from_flat_corearr(cr, n, r) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: cr(:) !< flattened sequence of cores to convert from
   integer,intent(in) :: n(:)           !< modes of tt
   integer,intent(in) :: r(0:size(n))   !< ranks of tt (mandatory argument!)
   !
   integer :: d, k, ofs, cs, sh(3)
   double precision, pointer :: p1(:)
#  ifdef VERBOSE4
   print '("[+][dtt_from_flat_corearr] entry")'
#  endif

      d = size(n)
      tt = dtt_tensor(n, r_=r)
      ofs = 0
      do k=1,d
         sh = [r(k-1), n(k), r(k)]
         cs = product(sh)
         p1(1:cs)=> tt% u(k)% p
         p1(1:cs) = cr(ofs+1:ofs+cs)
         ofs = ofs + cs
      enddo

#  ifdef VERBOSE4
   print '("[+][dtt_from_flat_corearr] exit")'
#  endif
   end function dtt_from_flat_corearr


   function dtt_from_2darray(a,eps_,rmax_) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: a(:,:)  !< array to convert to tt-format
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding

      tt = dtt_from_darray(a, shape(a), eps_, rmax_)

   end function dtt_from_2darray


   function dtt_from_3darray(a,eps_,rmax_) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: a(:,:,:)  !< array to convert to tt-format
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding

      tt = dtt_from_darray(a, shape(a), eps_, rmax_)

   end function dtt_from_3darray


   function dtt_from_4darray(a,eps_,rmax_) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: a(:,:,:,:)  !< array to convert to tt-format
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding

      tt = dtt_from_darray(a, shape(a), eps_, rmax_)

   end function dtt_from_4darray


   function dtt_from_5darray(a,eps_,rmax_) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: a(:,:,:,:,:)  !< array to convert to tt-format
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding

      tt = dtt_from_darray(a, shape(a), eps_, rmax_)

   end function dtt_from_5darray


   function dtt_from_6darray(a,eps_,rmax_) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: a(:,:,:,:,:,:)  !< array to convert to tt-format
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding

      tt = dtt_from_darray(a, shape(a), eps_, rmax_)

   end function dtt_from_6darray


   function dtt_from_7darray(a,eps_,rmax_) result(tt)
   type(dtt_tensor) :: tt
   double precision,intent(in) :: a(:,:,:,:,:,:,:)  !< array to convert to tt-format
   double precision,intent(in),optional :: eps_ !< approximation tolerance
   integer,intent(in),optional :: rmax_ !< maximum rank for the tt-rounding

      tt = dtt_from_darray(a, shape(a), eps_, rmax_)

   end function dtt_from_7darray

   !!
   !! MATRIX CONSTRUCTOR \ DESTRUCTOR
   !!
   !> creates empty matrix from arrays of dimensions and ranks (optional)
   function empty_dtt_matrix(q, s, r_, l_, m_) result(this)
   type(dtt_matrix):: this
   integer, intent(in) :: q(:)
   integer, intent(in) :: s(size(q))
   integer, intent(in), optional :: r_(0:size(q))
   integer, intent(in), optional :: l_, m_
   !
   integer:: d, k, l, m, sh(4)
   integer, allocatable :: r(:)
#  ifdef VERBOSE4
   print '("[+][empty_dtt_matrix] entry")'
#  endif

      l = 1; if (present(l_)) l = l_
      m = size(q); if (present(m_)) m = m_
      this% l = l; this% m = m
      d= m - l + 1
      if (l.gt.size(q)) error stop "in empty_dtt_matrix: l > |q|"
      if (m.gt.size(q)) error stop "in empty_dtt_matrix: m > |q|"
      if (d.lt.0) error stop "in empty_dtt_matrix: m + 1 < l"
      allocate(r(0:size(q))); r = 1
      if (present(r_)) r(l-1:m) = r_(l-1:m)
      if(any(q(l:m).le.0)) then
         write(*,*) 'empty_dtt_matrix: bad q: ',q(l:m)
         error stop
      endif
      if(any(s(l:m).le.0)) then
         write(*,*) 'empty_dtt_matrix: bad s: ',s(l:m)
         error stop
      endif
      if(any(r(l:m).le.0)) then
         write(*,*) 'empty_dtt_matrix: bad r: ',r(l:m)
         error stop
      endif
      this% q(l:m)= q(l:m)
      this% s(l:m)= s(l:m)
      this% n(l:m)= q(l:m)*s(l:m)
      this% r(l-1:m)= 1
      if (present(r_)) this% r(l-1:m)= r_(l-1:m)
      do k=l,m
         sh= [this% r(k-1), this% q(k), this% s(k), this% r(k)]
         call realloc_to_shape(this% u4(k), sh)
         this% u(k)% p(1:this% r(k-1), &
                       1:this% n(k),   &
                       1:this% r(k))=> this% u4(k)% p
         this% u(k)% p= 0d0
      enddo
      this% is_allocated= .true.
      deallocate(r)

#  ifdef VERBOSE4
   print '("[+][empty_dtt_matrix] exit")'
#  endif
   end function empty_dtt_matrix


   !> deep copy constructor
   function dttm_from_dtt_matrix(other) result(this)
   type(dtt_matrix):: this
   type(dtt_matrix), intent(in):: other
   integer:: d, l, m, k, sh(4)
#  ifdef VERBOSE4
   print '("[+][dttm_from_dtt_matrix] entry")'
#  endif

      l= other% l
      m= other% m
      this% l= l
      this% m= m
      this% q(l:m)= other% q(l:m)
      this% s(l:m)= other% s(l:m)
      this% n(l:m)= other% n(l:m)
      this% r(l-1:m)= other% r(l-1:m)
      do k=l,m
         sh= shape(other% u4(k)% p)
         call realloc_to_shape(this% u4(k), sh)
         this% u4(k)% p(:,:,:,:)= other% u4(k)% p(:,:,:,:)
         this% u(k)% p(1:this% r(k-1), &
                       1:this% n(k),   &
                       1:this% r(k))=> this% u4(k)% p
      enddo
      this% is_allocated= .true.
      this% nu_core= other% nu_core

#  ifdef VERBOSE4
   print '("[+][dttm_from_dtt_matrix] exit")'
#  endif
   end function dttm_from_dtt_matrix


   !> constructor from the base class object: tensor
   function dttm_from_dtt_tensor(dtt, q, s) result(this)
   type(dtt_tensor), intent(in):: dtt
   integer, intent(in) :: q(:)
   integer, intent(in) :: s(size(q))
   type(dtt_matrix):: this
   integer:: d, l, m, k, sh(3)
#  ifdef VERBOSE4
   print '("[+][dttm_from_dtt_tensor] entry")'
#  endif

      l= dtt% l
      m= dtt% m
      d= m - l + 1
      this% l= l
      this% m= m
      this% n(l:m)= dtt% n(l:m)
      this% q(l:m)= q(1:d)
      this% s(l:m)= s(1:d)
      this% r(l-1:m)= dtt% r(l-1:m)
      if (any(this% n(l:m).ne.q(1:d)*s(1:d))) then
         stop "[!][dttm_from_dtt_tensor] incompatible dimensions"
      endif
      do k=l,m
         sh= shape(dtt% u(k)% p)
         call realloc_to_shape(this% u(k), sh)
         this% u(k)% p(:,:,:)= dtt% u(k)% p(:,:,:)
         this% u4(k)% p(1:this% r(k-1), &
                        1:this% q(k), 1:this% s(k),  &
                        1:this% r(k))=> this% u(k)% p
      enddo
      this% is_allocated= .true.
      this% nu_core= dtt% nu_core

#  ifdef VERBOSE4
   print '("[+][dttm_from_dtt_tensor] exit")'
#  endif
   end function dttm_from_dtt_tensor


   !> dconstructor from a matrix
   function dttm_from_2dmatrix(A2, q, s, eps_, rmax_) result(this)
   use matrix_util, only: unravel_index, ravel_multi_index
   type(dtt_matrix):: this
   double precision,           intent(IN):: A2(:,:)   !< matrix to convert to tt-format
   integer,                    intent(IN):: q(:)      !< row dimensions
   integer,                    intent(IN):: s(size(q))!< col dimensions
   double precision, optional, intent(IN):: eps_      !< approximation tolerance
   integer,          optional, intent(IN):: rmax_     !< maximum rank for the tt-rounding
   !
   character(*), parameter :: subnam = "[dttm_from_2dmatrix]"
   double precision, allocatable :: a(:)
   integer :: i,j,k,ij,l,m,d,MA,NA,sh(2)
   integer,allocatable :: i_mi(:),j_mi(:),ij_mi(:)
#  ifdef VERBOSE4
   print '("[+]'//subnam//' entry")'
#  endif

      ! transpose array a(:) by interlacing indices
      d= size(q)
      MA = product(q)
      NA = product(s)
      sh = shape(A2)
      if (MA/=sh(1) .or. NA/=sh(2)) error stop '[!]'//subnam//': dimensions mismatch'
      allocate(a(MA*NA),ij_mi(d))
      do j = 1, NA
         j_mi = unravel_index(j - 1, s)
         do i = 1, MA
            i_mi = unravel_index(i - 1, q)
            ij_mi = 0
            do k = 1,size(i_mi)
               ij_mi(k) = i_mi(k)
            enddo
            do k = 1,size(j_mi)
               ij_mi(k) = ij_mi(k) + j_mi(k)*q(k) !! TODO: check
            enddo
            ij = ravel_multi_index(ij_mi, q*s) + 1
            a(ij) = A2(i,j)
         enddo
      enddo

      this = dttm_from_dtt_tensor(dtt_from_darray(a, q*s, eps_, rmax_), q, s)

#  ifdef VERBOSE4
   print '("[+]'//subnam//' exit")'
#  endif
   end function dttm_from_2dmatrix


   !> constructor from a cell array
   function dttm_from_cell_array(cell_array) result(this)
   use matlab_struct_module
   type(dtt_matrix):: this
   type(cell4d_array), intent(in):: cell_array
   integer:: l, m, d, k, sh(4)
#  ifdef VERBOSE4
   print '("[+][dttm_from_cell array] entry")'
#  endif

      d= size(cell_array% cells)
      m= d
      this% l= 1
      this% m= m

      do k=1,m
         sh= shape(cell_array% cells(k)% arr)
         this% r(k-1)= sh(1) ! TODO assert r(k) == sh(1)
         this% q(k)= sh(2)
         this% s(k)= sh(3)
         this% r(k)= sh(4)
      enddo
      this% n(1:m)= this% q(1:m)*this% s(1:m)

      do k=1,m
         sh= shape(cell_array% cells(k)% arr)
         call realloc_to_shape(this% u4(k), sh)
         this% u4(k)% p(:,:,:,:)= cell_array% cells(k)% arr(:,:,:,:)
         this% u(k)% p(1:this% r(k-1), &
                       1:this% n(k),   &
                       1:this% r(k))=> this% u4(k)% p
      enddo

      this% is_allocated= .true.
      this% nu_core= 0

#  ifdef VERBOSE4
   print '("[+][dttm_from_cell array] exit")'
#  endif
   end function dttm_from_cell_array


   !> constructor from a flattened array of the sequence of cores
   function dttm_from_flat_corearr(cr, q, s, r) result(ttm)
   type(dtt_matrix) :: ttm
   double precision,intent(in) :: cr(:) !< flattened sequence of cores to convert from
   integer,intent(in) :: q(:), s(:)     !< modes of ttm
   integer,intent(in) :: r(0:size(q))   !< ranks of ttm (mandatory argument!)
   !
   character(len=*),parameter    :: subnam='[dttm_from_flat_corearr]:'
   integer :: d, k, ofs, cs, sh(4)
   double precision, pointer :: p1(:)
#  ifdef VERBOSE4
   print '(A)' "[+]"//subnam//" entry"
#  endif

      d = size(s)
      if (size(q)/=d)   error stop "[!]"//subnam//": size(q) /= size(s)"
      ttm = dtt_matrix(q, s, r_=r)
      ofs = 0
      do k=1,d
         sh = [r(k-1), q(k), s(k), r(k)]
         cs = product(sh)
         p1(1:cs)=> ttm% u(k)% p
         p1(1:cs) = cr(ofs+1:ofs+cs)
         ofs = ofs + cs
      enddo

#  ifdef VERBOSE4
   print '(A)' "[+]"//subnam//" exit"
#  endif
   end function dttm_from_flat_corearr

   !!
   !! ALLOC \ DEALLOC
   !!
   subroutine alloc_d(this)
   class(dtt_tensor), intent(inout)  :: this
   character(len=*),parameter    :: subnam='[!][TENSOR][ALLOC]:'
   integer :: i,m, l, info
      l = this%l;  m = this%m
      if(m.lt.l) return
      if(l.le.0) error stop subnam//' %l should be > 0'
      if(m.gt.tt_size) error stop subnam//' %m exceeds tt_size'
      do i=l,m
         if (associated(this%u(i)%p)) deallocate(this%u(i)%p)
         allocate(this%u(i)%p(this%r(i-1), this%n(i), this%r(i)), stat=info)
         if (info.ne.0) error stop 'dtt_tensor allocate fail: no memory'
      end do
   end subroutine


   subroutine dealloc_d(this)
     class(dtt_tensor), intent(inout) :: this
     character(len=*),parameter   :: subnam='[!][TENSOR][DEALLOC]:'
     integer :: i
     do i=1,tt_size
       if(associated(this%u(i)%p))deallocate(this%u(i)%p)
     end do
   end subroutine


   !> dtt_tensor destructor (= deallocator)
   subroutine dealloc_dtt_tensor(this)
   type(dtt_tensor), intent(inout):: this
   integer:: k
#  ifdef VERBOSE4
   print '("[+][dealloc_dtt_tensor] entry")'
#  endif

      do k=this% l, this% m
         if (associated(this% u(k)% p)) deallocate(this% u(k)% p)
      enddo
      this% l= 1
      this% m= 0
      this% is_allocated= .false.

#  ifdef VERBOSE4
   print '("[+][dealloc_dtt_tensor] exit")'
#  endif
   end subroutine dealloc_dtt_tensor


   !> dtt_matrix destructor (= deallocator)
   subroutine dealloc_dtt_matrix(this)
   type(dtt_matrix), intent(inout):: this
   integer:: k
#  ifdef VERBOSE4
   print '("[+][dealloc_dtt_matrix] entry")'
   print '("[+][dealloc_dtt_matrix] l=",I2,", m=",I2)', this%l, this%m
#  endif

      do k=this% l, this% m
         if (associated(this% u4(k)% p)) deallocate(this% u4(k)% p)
         this% u(k)% p=> null()
         this% u4(k)% p=> null()
      enddo
      this% l= 1
      this% m= 0
      this% is_allocated= .false.

#  ifdef VERBOSE4
   print '("[+][dealloc_dtt_matrix] exit")'
#  endif
   end subroutine dealloc_dtt_matrix

   !!
   !! REALLOC
   !!
   subroutine realloc_pointd3(ptr, sh)
   type(pointd3), intent(inout):: ptr
   integer, intent(in):: sh(3)
   !
   integer:: s2(3)
      if (associated(ptr% p)) then
         s2= shape(ptr% p)
         if (sum((s2 - sh)**2).ne.0) then
            deallocate(ptr% p)
            allocate(ptr% p(sh(1), sh(2), sh(3)))
         endif
      else
         allocate(ptr% p(sh(1), sh(2), sh(3)))
      endif
   end subroutine realloc_pointd3


   subroutine realloc_pointd4(ptr, sh)
   type(pointd4), intent(inout):: ptr
   integer, intent(in):: sh(4)
   !
   integer:: s2(4)
      if (associated(ptr% p)) then
         s2= shape(ptr% p)
         if (sum((s2 - sh)**2).ne.0) then
            if(associated(ptr% p)) deallocate(ptr% p)
            allocate(ptr% p(sh(1), sh(2), sh(3), sh(4)))
         endif
      else
         allocate(ptr% p(sh(1), sh(2), sh(3), sh(4)))
      endif
   end subroutine realloc_pointd4


   subroutine realloc_pointz3(ptr, sh)
   type(pointz3), intent(inout):: ptr
   integer, intent(in):: sh(3)
   !
   integer:: s2(3)
      if (associated(ptr% p)) then
         s2= shape(ptr% p)
         if (sum((s2 - sh)**2).ne.0) then
            if(associated(ptr% p)) deallocate(ptr% p)
            allocate(ptr% p(sh(1), sh(2), sh(3)))
         endif
      else
         allocate(ptr% p(sh(1), sh(2), sh(3)))
      endif
   end subroutine realloc_pointz3


   subroutine realloc_pointz4(ptr, sh)
   type(pointz4), intent(inout):: ptr
   integer, intent(in):: sh(4)
   !
   integer:: s2(4)
      if (associated(ptr% p)) then
         s2= shape(ptr% p)
         if (sum((s2 - sh)**2).ne.0) then
            deallocate(ptr% p)
            allocate(ptr% p(sh(1), sh(2), sh(3), sh(4)))
         endif
      else
         allocate(ptr% p(sh(1), sh(2), sh(3), sh(4)))
      endif
   end subroutine realloc_pointz4

   !!
   !! Pointers mem address
   !!

   subroutine compute_ps(m,r,n,ps)
     implicit none
     integer, intent(in) :: m, r(*), n(*)
     integer, intent(out) :: ps(*)
     integer i
     ps(1)=1;
     do i=1,m
       ps(i+1) = ps(i) + r(i)*n(i)*r(i+1)
     end do
   end subroutine compute_ps
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_vec_ps(v,ps)
     implicit none
     type(dtt_tensor),intent(in)    :: v
     integer, intent(out) :: ps(*)
     integer i
     ps(1)=1;
     do i=1,v%m
       ps(i+1) = ps(i) + v%r(i)*v%n(i)*v%r(i+1)
     end do
   end subroutine get_vec_ps

   !!
   !! ASSIGNEMENT (OPERATOR =)
   !!

   subroutine dtt_assign(b,a)
   implicit none
   type(dtt_tensor),intent(inout) :: b
   type(dtt_tensor),intent(in)    :: a
   !
   integer                 :: k,l,m

       l=a%l; m=a%m
       b%l=l; b%m=m; b%n(l:m)=a%n(l:m); b%r(l-1:m)=a%r(l-1:m); call alloc(b)
       do k=l,m
          call dcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(k)%p,1)
       end do
       b% nu_core= a% nu_core

   end subroutine


   !> the dtt_tensor assignment operator, A = B
   subroutine dtt_tensor_assign(self, other)
   class(dtt_tensor), intent(inout) :: self
   type(dtt_tensor), intent(in), target :: other
#  ifdef VERBOSE4
   print '("[+][dtt_tensor_assign] entry")'
#  endif

      call dtt_tensor_baseclass_assignment(self, other)

#  ifdef VERBOSE4
   print '("[+][dtt_tensor_assign] exit")'
#  endif
   end subroutine dtt_tensor_assign

   !> the dtt_tensor assignment operator, A = B
   subroutine dtt_tensor_baseclass_assignment(self, other)
   class(dtt_tensor), intent(inout) :: self
   class(dtt_tensor), intent(in), target :: other
   integer:: k, k1, d, l, m, sh(3)
#  ifdef VERBOSE4
   print '("[+][dtt_tensor_baseclass_assignment] entry")'
#  endif

      l= other% l
      m= other% m
      d= m - l + 1
      if ((self% m - self% l + 1).ne.d) then
         if (self% is_allocated) call dealloc_dtt_tensor(self)
         self% l= l
         self% m= m
      endif
      self% n(self% l:self% m)= other% n(l:m)
      self% r(self% l-1:self% m)= other% r(l-1:m)

      if (other% is_rhs) then
         do k=l,m
            k1= k - l + self% l
            sh= shape(other% u(k)% p)
            self% u(k1)% p(1:self% r(k1-1), &
                           1:self% n(k1),   &
                           1:self% r(k1))=> other% u(k)% p
         enddo
      else
         do k=l,m
            k1= k - l + self% l
            sh= shape(other% u(k)% p)
            call realloc_to_shape(self% u(k1), sh)
            self% u(k1)% p(:,:,:)= other% u(k)% p(:,:,:)
         enddo
      endif
      self% is_allocated= .true.
      self%nu_core= other% nu_core

#  ifdef VERBOSE4
   print '("[+][dtt_tensor_baseclass_assignment] exit")'
#  endif
   end subroutine dtt_tensor_baseclass_assignment

   !> the dtt_matrix assignment operator, A = B
   subroutine dtt_matrix_assign(self, other)
   class(dtt_matrix), intent(inout) :: self
   type(dtt_matrix), intent(in)    :: other
   integer:: k, k1, d, l, m, sh(4)
#  ifdef VERBOSE4
   print '("[+][dtt_matrix_assign] entry")'
#  endif

      l= other% l
      m= other% m
      d= m - l + 1
      if ((self% m - self% l + 1).ne.d) then
         if (self% is_allocated) call dealloc_dtt_matrix(self)
         self% l= l
         self% m= m
      endif
      self% q(self% l:self% m)= other% q(l:m)
      self% s(self% l:self% m)= other% s(l:m)
      self% n(self% l:self% m)= other% n(l:m)
      self% r(self% l-1:self% m)= other% r(l-1:m)

      do k=l,m
         k1= k - l + self% l
         sh= shape(other% u4(k)% p)
         call realloc_to_shape(self% u4(k1), sh)
         self% u4(k1)% p(:,:,:,:)= other% u4(k)% p(:,:,:,:)
         self% u(k1)% p(1:self% r(k1-1), &
                        1:self% n(k1),   &
                        1:self% r(k1))=> self% u4(k1)% p
      enddo
      self% is_allocated= .true.
      self% nu_core = other% nu_core

#  ifdef VERBOSE4
   print '("[+][dtt_matrix_assign] exit")'
#  endif
   end subroutine dtt_matrix_assign

   !!
   !! COPY
   !!

   subroutine dtt_copy(a,b,low)
       implicit none
       type(dtt_tensor),intent(in) :: a
       type(dtt_tensor),intent(inout) :: b
       integer,intent(in),optional :: low
       integer :: k,l,m,ll,mm
       l=a%l; m=a%m
       if( present(low) ) then
           ll=low
       else
           ll=b%l
       end if
       mm=ll-l+m
       b%l=ll; b%m=mm; b%n(ll:mm)=a%n(l:m); b%r(ll-1:mm)=a%r(l-1:m);
       if(.not.all(a%n(l:m)>0))return;if(.not.all(a%r(l-1:m)>0))return;call alloc(b)
       do k=l,m; call dcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(ll-l+k)%p,1); end do
   end subroutine


   !!
   !! READY
   !!

   logical function dtt_ready(arg) result(l)
     implicit none
     type(dtt_tensor),intent(in) :: arg
     integer :: i
     l=(arg%l .le. arg%m)
     if(.not.l)return
     l=all(arg%n(arg%l:arg%m)>0)
     if(.not.l)return
     l=all(arg%r(arg%l-1:arg%m)>0)
     if(.not.l)return
     do i=arg%l,arg%m; l=l.and.associated(arg%u(i)%p); enddo
     if(.not.l)return
     do i=arg%l,arg%m
       l=l.and.(size(arg%u(i)%p,1).eq.arg%r(i-1))
       l=l.and.(size(arg%u(i)%p,2).eq.arg%n(i))
       l=l.and.(size(arg%u(i)%p,3).eq.arg%r(i))
     enddo
     if(.not.l)return
     return
   end function

   !!
   !! MEMORY
   !!

   integer function dtt_mem(arg) result (sz)
     implicit none
     type(dtt_tensor),intent(in) :: arg
     integer :: i
     sz=0
     do i=arg%l,arg%m
       sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
     end do
   end function


   double precision function dtt_mb(arg) result (sz)
     implicit none
     type(dtt_tensor),intent(in) :: arg
     integer :: i
     sz=0.d0
     do i=arg%l,arg%m
       sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
     end do
     sz=sz*8.d0/(2**20)
   end function


   !!
   !! MULT
   !!

   type(dtt_tensor) function dtt_mul_d(b,a) result(c)
   class(dtt_tensor),intent(in) :: b
   double precision,intent(in) :: a
   integer :: k,l,m
      l=b%l; m=b%m; c%l=l; c%m=m; c%n=b%n; c%r=b%r; call alloc(c)
      do k=l,m
         call dcopy(b%r(k-1)*b%n(k)*b%r(k),b%u(k)%p,1,c%u(k)%p,1)
      end do
      call dscal(b%r(l-1)*b%n(l)*b%r(l),a,c%u(l)%p,1)
      c%nu_core= b% nu_core
   end function dtt_mul_d


   type(dtt_tensor) function d_mul_dtt(a,b) result(c)
   double precision,intent(in) :: a
   class(dtt_tensor),intent(in) :: b
   integer :: k,l,m
      l=b%l; m=b%m; c%l=l; c%m=m; c%n=b%n; c%r=b%r; call alloc(c)
      do k=l,m
         call dcopy(b%r(k-1)*b%n(k)*b%r(k),b%u(k)%p,1,c%u(k)%p,1)
      end do
      call dscal(b%r(l-1)*b%n(l)*b%r(l),a,c%u(l)%p,1)
   end function d_mul_dtt


   function dtt_mult(x,y) result(z)
   class(dtt_tensor), intent(in), target :: x,y
   type(dtt_tensor)                      :: z
   !
   character(len=*),parameter      :: subnam='[!][TENSOR][dtt_mult]'
   integer                         :: l,m,k,i,j,ii,jj,kk,p
   integer,pointer                 :: r(:),q(:),n(:)
   logical                         :: zx,zy
   integer, allocatable            :: rr(:)
#  ifdef VERBOSE4
   print '("[+][dtt_mult] entry")'
#  endif

      !! TODO: handle the case where dimensions are not the same
      if (x%l.gt.x%m.or.y%l.gt.y%m) stop subnam//": l > m"
      if (x%l.ne.y%l.or.x%m.ne.y%m) stop subnam//": length mismatch"
      l=x%l; m=x%m
      if (any(x%n(l:m).le.0).or.any(y%n(l:m).le.0)) stop subnam//": degenerate ranks"
      if (any(x%n(l:m).ne.y%n(l:m))) stop subnam//": size mismatch"
      n=>x%n;r=>x%r;q=>y%r
      if(r(l-1).ne.q(l-1) .or. r(m).ne.q(m)) stop subnam//": border ranks mismatch"
      if (m.eq.l) stop subnam//": special code required for m=l!"

      ! initialize empty dtt_matrix
      allocate(rr(0:size(n)))
      rr(l-1:m)= r(l-1:m)*q(l-1:m)
      z= dtt_tensor(n, l_=l, m_=m, r_=rr)
      z%r(l-1)= r(l-1); z%r(m)= r(m)

      forall(i=1:r(l-1),j=1:n(l),k=1:r(l),kk=1:q(l)) &
         z%u(l)%p(i,j,k+(kk-1)*r(l))= x%u(l)%p(i,j,k)*y%u(l)%p(i,j,kk)
      do p=l+1,m-1
         forall(i=1:r(p-1),ii=1:q(p-1),j=1:n(p),k=1:r(p),kk=1:q(p)) &
            z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))= x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
      end do
      forall(i=1:r(m-1),ii=1:q(m-1),j=1:n(m),k=1:r(m)) &
         z%u(m)%p(i+(ii-1)*r(m-1),j,k)= x%u(m)%p(i,j,k)*y%u(m)%p(ii,j,k)
      deallocate(rr)
      z% nu_core = 0

#  ifdef VERBOSE4
   print '("[+][dtt_mult] exit")'
#  endif
   end function dtt_mult

   !!
   !! ADD z=x+y
   !!

   type(dtt_tensor) function dtt_plus_d(a,b) result(c)
     class(dtt_tensor),intent(in),target :: a
     double precision,intent(in) :: b
     integer :: i,j,k,l,m,p
     integer,pointer :: n(:),r(:)
     if(.not.(ready(a)))return
     l=a%l; m=a%m; n=>a%n; r=>a%r
     c%l=a%l; c%m=a%m; c%n=a%n; c%r=0; c%r(l-1)=a%r(l-1); c%r(m)=a%r(m); c%r(l:m-1)=a%r(l:m-1)+1
     call alloc(c)
     forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) c%u(l)%p(i,j,     k)=a%u(l)%p(i,j,k)
     forall(i=1:r(l-1),j=1:n(l))          c%u(l)%p(i,j,r(l)+1)=b
     do p=l+1,m-1
      forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) c%u(p)%p(       i,j,     k)=a%u(p)%p(i,j,k)
      forall(i=1:r(p-1),j=1:n(p))          c%u(p)%p(       i,j,r(p)+1)=0.d0
      forall(           j=1:n(p),k=1:r(p)) c%u(p)%p(r(p-1)+1,j,     k)=0.d0
      forall(           j=1:n(p))          c%u(p)%p(r(p-1)+1,j,r(p)+1)=1.d0
     end do
     forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) c%u(m)%p(       i,j,k)=a%u(m)%p(i,j,k)
     forall(           j=1:n(m),k=1:r(m)) c%u(m)%p(r(m-1)+1,j,k)=1.d0
     c% nu_core= 0
   end function dtt_plus_d


   type(dtt_tensor) function d_plus_dtt(b,a) result(c)
   double precision,intent(in) :: b
   class(dtt_tensor),intent(in),target :: a

      c = dtt_plus_d(a,b)

   end function


   !> Addition of two dtt_tensor base-class objects
   subroutine dtt_add_baseclass(z,x,y)
     implicit none
     type(dtt_tensor), intent(inout) :: z
     class(dtt_tensor),intent(in),target :: x,y
     double precision,parameter     :: a=1.d0, b=1.d0 !< actually, y <- b*y + a*x
     character(len=*),parameter     :: subnam='[!][TENSOR][dtt_add_baseclass]'
     integer                        :: l,m,k,i,j,p
     integer,pointer                :: r(:),q(:),n(:)
     logical                        :: zx,zy
     zx=.true.; if(x%m.ge.x%l)zx=.not.(all(x%n(x%l:x%m)>0).and.all(x%r(x%l-1:x%m)>0))
     zy=.true.; if(y%m.ge.y%l)zy=.not.(all(y%n(y%l:y%m)>0).and.all(x%r(x%l-1:x%m)>0))

     if(zx.and.zy)return
     if(zx)then;call dscal(y%r(l-1)*y%n(l)*y%r(l),b,y%u(l)%p,1);return;endif
     if(zy)then;z=a*x;return;endif

     if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': length mismatch';stop;endif
     l=x%l; m=x%m
     if(.not.all(x%n(l:m)==y%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
     n=>x%n;r=>x%r;q=>y%r
     if(r(l-1).ne.q(l-1) .or. r(m).ne.q(m))then;write(*,*)subnam,': border ranks mismatch';stop;endif

     z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1)=r(l-1); z%r(m)=r(m); z%r(l:m-1)=r(l:m-1)+q(l:m-1)
     if(l.eq.m)then;z%u(l)%p=a*x%u(l)%p + b*y%u(l)%p;return;endif
     call alloc(z)

     forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) z%u(l)%p(i,j,     k)=a*(x%u(l)%p(i,j,k))
     forall(i=1:r(l-1),j=1:n(l),k=1:q(l)) z%u(l)%p(i,j,r(l)+k)=b*(y%u(l)%p(i,j,k))
     do p=l+1,m-1
       forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(       i,j,     k)=x%u(p)%p(i,j,k)
       forall(i=1:r(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(       i,j,r(p)+k)=0.d0
       forall(i=1:q(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(r(p-1)+i,j,     k)=0.d0
       forall(i=1:q(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(r(p-1)+i,j,r(p)+k)=y%u(p)%p(i,j,k)
     end do
     forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(       i,j,k)=x%u(m)%p(i,j,k)
     forall(i=1:q(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(r(m-1)+i,j,k)=y%u(m)%p(i,j,k)
     z% nu_core= 0
   end subroutine dtt_add_baseclass


   type(dtt_tensor) function dtt_add(x,y) result(z)
   type(dtt_tensor),intent(in),target :: x,y
#  ifdef VERBOSE4
   print '("[+][dtt_add] entry")'
#  endif

      call dtt_add_baseclass(z, x, y)

#  ifdef VERBOSE4
   print '("[+][dtt_add] exit")'
#  endif
   end function dtt_add


   !> addition of two dtt_matrix objects
   function     dttm_add(left, right) result(total)
   class(dtt_matrix), intent(in) :: left, right
   type(dtt_matrix) :: total
   !
   integer:: k, l, m, d, q, s, r1, r2
   integer, allocatable:: r(:)
#  ifdef VERBOSE4
   print '("[+][dttm_add] entry")'
#  endif

      ! check that the dimensions match
      l= left% l
      m= left% m
      d= m - l + 1
      if (m - l.ne.right% m - right% l) then
         stop "ERROR: in dttm_add: incompatible operators"
      else
         if (any(left% q(l:m).ne.right% q(right% l:right% m)) &
         .or.any(left% s(l:m).ne.right% s(right% l:right% m))) then
            stop "ERROR: in dttm_add: incompatible dimensions"
         endif
      endif

      ! add ranks
      allocate(r(0:m))
      r(1:m-1)= left% r(l:m - 1) + right% r(right% l:right% m - 1)
      r(0)= 1
      r(m)= 1

      ! initialize empty dtt_matrix
      total= dtt_matrix(left% q(l:m), left% s(l:m), l_=l, r_=r, m_=m)

      ! write cores
      q=  total% q(1)
      s=  total% s(1)
      r2= r(1)

      total%    u4(l)% p(1,1:q,1:s,1:left% r(1)) &
        = left% u4(l)% p(1,1:q,1:s,1:left% r(1))
      total%    u4(l)% p(1,1:q,1:s,r2-right% r(1)+1:r2) &
        = right%u4(l)% p(1,1:q,1:s, 1:right% r(1))

      do k=l+1,m-1
         r1= total% r(k-1)
         q=  total% q(k)
         s=  total% s(k)
         r2= total% r(k)

         total% u4(k)% p= 0d0
         total%   u4(k)% p(1:left% r(k-1),      1:q,1:s, 1:left% r(k))&
          = left% u4(k)% p(1:left% r(k-1),      1:q,1:s, 1:left% r(k))
         total%   u4(k)% p(r1-right%r(k-1)+1:r1,1:q,1:s,r2-right%r(k)+1:r2)&
          = right%u4(k)% p( 1:right%r(k-1),     1:q,1:s, 1:right%r(k))
      enddo

      r1= r(m-1)
      q=  total% q(m)
      s=  total% s(m)

      total%   u4(m)% p(1:left% r(m-1),       1:q,1:s,1) &
       = left% u4(m)% p(1:left% r(m-1),       1:q,1:s,1)
      total%   u4(m)% p(r1-right% r(m-1)+1:r1,1:q,1:s,1) &
       = right%u4(m)% p( 1:right% r(m-1),     1:q,1:s,1)

      total% nu_core= 0

      ! cleanup
      deallocate(r)

#  ifdef VERBOSE4
   print '("[+][dttm_add] exit")'
#  endif
   end function dttm_add

   !!
   !! SUBTRACT z=x-y
   !!

   !> Subtraction of two dtt_tensor base-class objects
   subroutine dtt_sub_baseclass(z,x,y)
     ! y=b*y+a*x
     implicit none
     type(dtt_tensor), intent(inout) :: z
     class(dtt_tensor),intent(in),target :: x,y
     double precision,parameter     :: a=1.d0, b=-1.d0
     character(len=*),parameter     :: subnam='[!][TENSOR][dtt_sub_baseclass]'
     integer                        :: l,m,k,i,j,p
     integer,pointer                :: r(:),q(:),n(:)
     logical                        :: zx,zy
     zx=.true.; if(x%m.ge.x%l)zx=.not.(all(x%n(x%l:x%m)>0).and.all(x%r(x%l-1:x%m)>0))
     zy=.true.; if(y%m.ge.y%l)zy=.not.(all(y%n(y%l:y%m)>0).and.all(x%r(x%l-1:x%m)>0))

     if (zx.and.zy) return
     if (zx) then
        z= -y
        return
     endif
     if (zy) then
        z= x
        return
     endif

     if (x%l.ne.y%l .or. x%m.ne.y%m) error stop subnam//': length mismatch'
     l=x%l
     m=x%m

     if (.not.all(x%n(l:m)==y%n(l:m))) error stop subnam//': size mismatch'
     n=>x%n
     r=>x%r
     q=>y%r

     if (r(l-1).ne.q(l-1) .or. r(m).ne.q(m)) error stop subnam//': border ranks mismatch'
     z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1)=r(l-1); z%r(m)=r(m); z%r(l:m-1)=r(l:m-1)+q(l:m-1)
     call alloc(z)

     if (l.eq.m) then
        z%u(l)%p= x%u(l)%p - y%u(l)%p
        return
     endif

     forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) z%u(l)%p(i,j,     k)=a*(x%u(l)%p(i,j,k))
     forall(i=1:r(l-1),j=1:n(l),k=1:q(l)) z%u(l)%p(i,j,r(l)+k)=b*(y%u(l)%p(i,j,k))
     do p=l+1,m-1
       forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(       i,j,     k)=x%u(p)%p(i,j,k)
       forall(i=1:r(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(       i,j,r(p)+k)=0.d0
       forall(i=1:q(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(r(p-1)+i,j,     k)=0.d0
       forall(i=1:q(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(r(p-1)+i,j,r(p)+k)=y%u(p)%p(i,j,k)
     end do
     forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(       i,j,k)=x%u(m)%p(i,j,k)
     forall(i=1:q(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(r(m-1)+i,j,k)=y%u(m)%p(i,j,k)
   end subroutine dtt_sub_baseclass


   type(dtt_tensor) function dtt_sub(x,y) result(z)
   type(dtt_tensor),intent(in),target :: x,y
#  ifdef VERBOSE4
   print '("[+][dtt_sub] entry")'
#  endif

      call dtt_sub_baseclass(z, x, y)

#  ifdef VERBOSE4
   print '("[+][dtt_sub] exit")'
#  endif
   end function dtt_sub


   !> subtract double from dtt_tensor: a - b
   type(dtt_tensor) function dtt_minus_d(a,b) result(c)
   class(dtt_tensor),intent(in),target :: a
   double precision,intent(in) :: b

      c = dtt_plus_d(a, -b)

   end function dtt_minus_d


   !> subtract dtt_tensor from a double: b - a
   type(dtt_tensor) function d_minus_dtt(b,a) result(c)
   double precision,intent(in) :: b
   class(dtt_tensor),intent(in),target :: a

      c = dtt_plus_d(-a, b)

   end function d_minus_dtt


   !> subtraction of two dtt_matrix objects
   function dttm_sub(left, right) result(difference)
   class(dtt_matrix), intent(in) :: left, right
   type(dtt_matrix) :: difference
   !
   integer:: q, s, rr, rl
#  ifdef VERBOSE4
   print '("[+][dttm_sub] entry")'
#  endif

      ! add, then negate the first sub-core
      difference= left + right

      ! write cores
      q=  difference% q(1)
      s=  difference% s(1)
      rl= left% r(1)
      rr= right% r(1)

      difference%    u4(1)% p(1,1:q,1:s,rl+1:rr+rl) &
        =     -right%u4(1)% p(1,1:q,1:s,   1:rr)

#  ifdef VERBOSE4
   print '("[+][dttm_sub] exit")'
#  endif
   end function dttm_sub

   !> division operator a/b
   type(dtt_tensor) function dtt_div_d(a,b) result(c)
   class(dtt_tensor),intent(in)  :: a
   double precision,intent(in) :: b
   integer :: k,l,m
   double precision :: d

      c = dtt_mul_d(a, 1d0/b)

   end function


   !> unary minus operator: a -> (-a)
   type(dtt_tensor) function minus_dtt(a) result(b)
   class(dtt_tensor),intent(in) :: a

      b = (-1d0)*a

   end function

   !> integer power: c = a**n
   type(dtt_tensor) function dtt_pwr_i(a,n) result (c)
   class(dtt_tensor),intent(in)  :: a
   integer,        intent(in)  :: n
   !
   integer :: i

   if (n.eq.0) then
      c= dtt_tensor_ones_like(a)
   else if (n.eq.1) then
      c= a
   else if (n.lt.0) then
      error stop "dtt_pwr_i: not implemented!"
   else
      c= dtt_mult(a, a)
      do i=3,n
        c= dtt_mult(c, a)
      enddo
   endif
   end function



   !!
   !! DOT PRODUCT
   !!

   double precision function dtt_dot(x,y) result(dot)
     implicit none
     type(dtt_tensor),intent(in) :: x,y
     character(len=*),parameter :: subnam='dtt_dot'
     integer :: i,l,m,info,rx(0:tt_size),ry(0:tt_size),n(tt_size)
     double precision, allocatable :: phi(:), res(:)
     if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': dimensions not match';stop;endif
     if(.not.all(x%n(x%l:x%m)==y%n(y%l:y%m)))then;write(*,*)subnam,': sizes not match';stop;endif
     l=x%l;m=x%m; rx=x%r;ry=y%r; n=x%n
     !!
     allocate(phi(maxval(rx(l-1:m)*ry(l-1:m))), res(maxval(rx(l-1:m-1)*n(l:m)*ry(l:m))), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
     !!
     phi(1)=1.d0
     do i = l,m
       call dgemm('n','n', rx(i-1),n(i)*ry(i),ry(i-1), 1d0,phi,rx(i-1), y%u(i)%p,ry(i-1), 0d0,res,rx(i-1))
       call dgemm('t','n', rx(i),ry(i),rx(i-1)*n(i), 1d0,x%u(i)%p,rx(i-1)*n(i), res,rx(i-1)*n(i), 0d0,phi,rx(i))
     end do
     dot=phi(1)
     deallocate(phi,res)
   end function

   !!
   !! ONES
   !!

   !> Function for creating a dtt_tensor of array of ones
   function dtt_tensor_ones(n, l_, m_) result(this)
   type(dtt_tensor)    :: this
   integer, intent(in) :: n(:)
   integer, intent(in), optional :: l_, m_
   !
   integer :: k, l, m
#  ifdef VERBOSE4
   print '("[+][dtt_tensor_ones] entry")'
#  endif
      l= 1; if (present(l_)) l = l_
      m= size(n); if (present(m_)) m = m_
      this = dtt_tensor(n, l_=l, m_=m)
      do k=l,m
         this% u(k)% p= 1d0
      end do
      this% nu_core= 0
#  ifdef VERBOSE4
   print '("[+][dtt_tensor_ones] exit")'
#  endif
   end function dtt_tensor_ones


   !> Function for creating a dtt_tensor of array of ones like a different tensor
   function dtt_tensor_ones_like(v) result(this)
   type(dtt_tensor) :: this
   type(dtt_tensor),intent(in) :: v

      this = dtt_tensor_ones(v%n, l_=v%l, m_=v%m)

   end function dtt_tensor_ones_like


   !> Function for creating a dtt_matrix of array of ones
   function dtt_matrix_ones(q, s, l_, m_) result(this)
   type(dtt_matrix)    :: this
   integer, intent(in) :: q(:), s(:)
   integer, intent(in), optional :: l_, m_
   !
   integer :: k, l, m
#  ifdef VERBOSE4
   print '("[+][dtt_matrix_ones] entry")'
#  endif
      l= 1; if (present(l_)) l = l_
      m= size(q); if (present(m_)) m = m_
      this = empty_dtt_matrix(q, s, l_=l, m_=m)
      do k=l,m
         this% u(k)% p= 1d0
      end do
      this% nu_core= 0
#  ifdef VERBOSE4
   print '("[+][dtt_matrix_ones] exit")'
#  endif
   end function dtt_matrix_ones


   !> Function for creating a dtt_matrix of array of ones like a different tensor
   function dtt_matrix_ones_like(A) result(this)
   type(dtt_matrix) :: this
   type(dtt_matrix),intent(in) :: A

      this = dtt_matrix_ones(A%q, A%s, l_=A%l, m_=A%m)

   end function dtt_matrix_ones_like

   !!
   !!RANDOM
   !!

   !> Create a dtt_tensor with specified dims and ranks with random cores
   !!
   !! Parameters:
   !! * n(:)  : tensor modes
   !! * r_(:) : tensor ranks [1]
   !! * l_    : starting index in n(:), l_-1 in r(:)
   !! * m_    : ending index (inclusive) in n(:) and in r(:)
   !! * randl_, randr_ : generate random numbers in the interval [randl_, randr_]
   !!
   function dtt_tensor_rand(n, r_, l_, m_, randl_, randr_) result(tens)
   use rnd_lib, only : random
   integer, intent(in) :: n(:)
   integer, intent(in), optional :: r_(0:)
   integer, intent(in), optional :: l_, m_
   double precision, optional :: randr_, randl_
   type(dtt_tensor)    :: tens
   !
   integer :: k, l, m, d
   integer, allocatable :: r(:)
   double precision :: randr, randl
#  ifdef VERBOSE4
   print '("[+][dtt_tensor_rand] entry")'
#  endif

     l = 1; if (present(l_)) l = l_
     m = size(n); if (present(m_)) m = m_
     d = m - l + 1
     if (d.lt.0) return

     randr =-1d0; if (present(randr_)) randr = randr_
     randl = 1d0; if (present(randl_)) randl = randl_

     allocate(r(0:size(n))); r = 1
     if (present(r_)) r(l-1:m) = r_(l-1:m)
     tens = dtt_tensor(n, r_=r, l_=l, m_=m)
     do k=l,m
       call random(tens%u(k)%p)
       tens%u(k)%p = randr + (randl - randr)*tens%u(k)%p
     end do
     tens% nu_core= 0
     deallocate(r)

#  ifdef VERBOSE4
   print '("[+][dtt_tensor_rand] exit")'
#  endif
   end function dtt_tensor_rand


   function dtt_tensor_rand_like(v) result(tens)
   !
   ! Create a tensor of the same shape as v with random cores
   !
   type(dtt_tensor), intent(in)    :: v
   type(dtt_tensor)                :: tens

     tens = dtt_tensor_rand(v% n, r_=v% r, l_=v% l, m_=v% m)

   end function dtt_tensor_rand_like


   !> Create a dtt_tensor with specified dims and ranks with random cores
   function dtt_matrix_rand(q, s, r_, l_, m_) result(this)
   !
   ! Parameters:
   ! * q(:), s(:) : tt-matrix modes
   ! * r_(:) : tensor ranks [1]
   ! * l_    : starting index in n(:), l_-1 in r(:)
   ! * m_    : ending index (inclusive) in n(:) and in r(:)
   !
   use rnd_lib, only: random
   type(dtt_matrix):: this
   integer, intent(in) :: q(:)
   integer, intent(in) :: s(size(q))
   integer, intent(in), optional :: r_(0:size(q))
   integer, intent(in), optional :: l_, m_
   !
   integer:: k, l, m, d
   integer, allocatable :: r(:)
#  ifdef VERBOSE4
   print '("[+][dtt_matrix_rand] entry")'
#  endif

      l = 1; if (present(l_)) l = l_
      m = size(q); if (present(m_)) m = m_
      d = m - l + 1
      if (d.lt.0) return

      allocate(r(0:size(q)))
      if (present(r_)) r(l-1:m) = r_(l-1:m)
      this = empty_dtt_matrix(q, s, r_=r, l_=l, m_=m)
      do k=l,m
         call random(this% u4(k)% p)
      enddo
      this% nu_core= 0
      deallocate(r)

#  ifdef VERBOSE4
   print '("[+][dtt_matrix_rand] exit")'
#  endif
   end function dtt_matrix_rand


   !> Create random matrix with shape like a give matrix A
   function dtt_matrix_rand_like(A) result(this)
   type(dtt_matrix), intent(in)    :: A
   type(dtt_matrix)                :: this

     this = dtt_matrix_rand(A% q, A% s, r_=A% r, l_=A% l, m_=A% m)

   end function dtt_matrix_rand_like


   !!
   !! ZEROS
   !!

   !> function for creating a dtt_tensor of array of zeros
   !!
   !! Parameters:
   !! * n(:)  : tensor modes
   !! * r_(:) : tensor ranks [1]
   !! * l_    : starting index in n(:), l_-1 in r(:)
   !! * m_    : ending index (inclusive) in n(:) and in r(:)
   !!
   function dtt_tensor_zeros(n,r_,l_,m_) result(this)
   integer, intent(in) :: n(:)
   integer, intent(in), optional :: r_(0:size(n)), l_, m_
   type(dtt_tensor)            :: this
   !
   integer :: k, l, m, d
   integer, allocatable :: r(:)
#  ifdef VERBOSE4
   print '("[+][dtt_tensor_zeros] entry")'
#  endif

     l = 1; if (present(l_)) l = l_
     m = size(n); if (present(m_)) m = m_
     d = m - l + 1
     if (d.lt.0) return
     allocate(r(0:size(n)))
     r = 1; if(present(r_)) r = r_
     this = dtt_tensor(n, r_=r, l_=l, m_=m)
     do k=l,m
        this%u(k)%p = 0d0
     enddo
     this% nu_core= 0
     deallocate(r)

#  ifdef VERBOSE4
   print '("[+][dtt_tensor_zeros] exit")'
#  endif
   end function dtt_tensor_zeros


   !> Function for creating a dtt_tensor of array of zeros like v
   function dtt_tensor_zeros_like(v) result(tens)
   type(dtt_tensor), intent(in)    :: v
   type(dtt_tensor)                :: tens
#  ifdef VERBOSE4
   print '("[+][dtt_tensor_zeros_like] entry")'
#  endif

      tens= dtt_tensor_zeros(v%n, r_=v%r, l_=v%l, m_=v%m)

#  ifdef VERBOSE4
   print '("[+][dtt_tensor_zeros_like] exit")'
#  endif
   end function dtt_tensor_zeros_like


   !> function for creating a dtt_matrix of zeros
   !!
   !! Parameters:
   !! * iq(:), s(:) : matrix modes
   !! * r_(:) : tensor ranks [1]
   !! * l_    : starting index in n(:), l_-1 in r(:)
   !! * m_    : ending index (inclusive) in n(:) and in r(:)
   !!
   function dtt_matrix_zeros(q, s, r_, l_, m_) result(this)
   integer, intent(in) :: q(:), s(size(q))
   integer, intent(in), optional :: r_(0:size(q)), l_, m_
   type(dtt_matrix)            :: this
   !
   integer :: k, l, m, d
   integer, allocatable :: r(:)
#  ifdef VERBOSE4
   print '("[+][dtt_matrix_zeros] entry")'
#  endif

     l = 1; if (present(l_)) l = l_
     m = size(q); if (present(m_)) m = m_
     d = m - l + 1
     if (d.lt.0) return
     allocate(r(0:size(q)))
     r = 1; if(present(r_)) r = r_
     this = empty_dtt_matrix(q, s, r_=r, l_=l, m_=m)
     do k=l,m
        this%u(k)%p = 0d0
     enddo
     this% nu_core= 0
     deallocate(r)

#  ifdef VERBOSE4
   print '("[+][dtt_matrix_zeros] exit")'
#  endif
   end function dtt_matrix_zeros


   !> Function for creating a dtt_matrix of zeros like v
   function dtt_matrix_zeros_like(A) result(this)
   type(dtt_matrix), intent(in)    :: A
   type(dtt_matrix)                :: this
#  ifdef VERBOSE4
   print '("[+][dtt_matrix_zeros_like] entry")'
#  endif

      this= dtt_matrix_zeros(A%q, A%s, r_=A%r, l_=A%l, m_=A%m)

#  ifdef VERBOSE4
   print '("[+][dtt_matrix_zeros_like] exit")'
#  endif
   end function dtt_matrix_zeros_like


   !!
   !! ORT:
   !!
   subroutine dtt_ort(arg)
   ![VEC] ortogonalize from left
     implicit none
     type(dtt_tensor),intent(inout),target :: arg
     character(len=*),parameter :: subnam='dtt_ort'
     integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
     integer,pointer :: r(:),n(:)
     double precision,allocatable :: work(:),tau(:),mat(:),u(:)
     double precision :: err,nrm,lognrm
     double precision,external :: dnrm2
     !!
     l=arg%l; m=arg%m
     if(m.lt.l)return
     r=>arg%r; n=>arg%n
     !!
     nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
     lwork=128*nn*rr
     allocate(work(lwork),tau(nn*rr), mat(rr*nn*rr),u(rr*nn*rr), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif
     !!
     lognrm=0.d0
     do k=l,m-1
       mm=r(k-1)*n(k); nn=r(k); mn=min(mm,nn); kk=n(k+1)*r(k+1)
       call dcopy(mm*nn, arg%u(k)%p,1,u,1)
       call dgeqrf(mm,nn, u,mm,tau,work,lwork,info)
       if(info.ne.0)then; write(*,*) subnam,': dgeqrf info: ',info; stop; end if
       do j=1,nn
         forall(i=1:min(j,mm))    mat(i+(j-1)*mn)=u(i+(j-1)*mm)
         forall(i=min(j,mm)+1:mn) mat(i+(j-1)*mn)=0.d0
       end do
       nrm=dnrm2(mn*nn,mat,1)
       if(nrm.ne.0.d0)then
         call dscal(mn*nn,1.d0/nrm,mat,1)
         lognrm=lognrm+dlog(nrm)
       endif
       !!
       call dorgqr(mm,mn,mn,u,mm,tau,work,lwork,info)
       if(info.ne.0)then; write(*,*) subnam,': dorgqr info: ',info; stop; end if
       !!
       call dcopy(mm*mn, u,1,arg%u(k)%p,1)
       call dgemm('n','n',mn,kk,nn,1.d0,mat,mn,arg%u(k+1)%p,nn,0.d0,u,mn)
       if(r(k).ne.mn)then
         call dcopy(mm*mn, arg%u(k)%p,1,mat,1)
         deallocate(arg%u(k)%p,arg%u(k+1)%p)
         r(k)=mn
         allocate(arg%u(k)%p(r(k-1),n(k),r(k)),arg%u(k+1)%p(r(k),n(k+1),r(k+1)))
         call dcopy(mm*mn, mat,1,arg%u(k)%p,1)
       end if
       call dcopy(mn*kk, u,1,arg%u(k+1)%p,1)
     end do
     deallocate(work,tau,mat,u) ! this may malfunction
     !!
     nrm=dnrm2(r(m-1)*n(m)*r(m),arg%u(m)%p,1)
     if(nrm.ne.0.d0)then
       call dscal(r(m-1)*n(m)*r(m),1.d0/nrm,arg%u(m)%p,1)
       lognrm=lognrm+dlog(nrm)
     endif
     !!
     lognrm=lognrm/(m-l+1)
     nrm=dexp(lognrm)
     do k=l,m
       call dscal(r(k-1)*n(k)*r(k),nrm,arg%u(k)%p,1)
     end do
     arg% nu_core= m
     end subroutine


   !!
   !! SVD
   !!
   function dtt_round(tt,tol_,rmax_) result(newtt)
   implicit none
   class(dtt_tensor),intent(in),target :: tt
   double precision,intent(in),optional :: tol_
   integer,intent(in),optional :: rmax_
   type(dtt_tensor) :: newtt

      newtt = tt
      call dtt_svd(newtt, tol_, rmax_)

   end function dtt_round


   subroutine dtt_svd(tt,tol_,rmax_)
   use mat_lib, only : svd
   implicit none
   class(dtt_tensor),intent(inout),target :: tt
   double precision,intent(in),optional :: tol_
   integer,intent(in),optional :: rmax_
   !
   character(len=*),parameter :: subnam='[!][TENSOR][dtt_svd]'
   integer :: l,m,k,i,j,info,nn,rr,mn,mm,kk
   integer,pointer :: r(:),n(:)
   double precision,pointer,contiguous :: mat(:,:),u(:,:),s(:),v(:,:),tmp(:)
   double precision :: err,nrm,lognrm,tol
   double precision,external :: dnrm2
#  ifdef VERBOSE4
   print '("[+][dtt_svd] entry")'
#  endif
   !!
      l=tt%l; m=tt%m
      if(m.le.l)return
      r=>tt%r; n=>tt%n
      nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
      allocate(tmp(rr*nn*rr), stat=info)
      if(info.ne.0) error stop subnam//': no memory'
      tol=1d-14;if(present(tol_)) tol=tol_
      !!
      call dtt_ort(tt)
      lognrm=0.d0
      !!
      do k=m,l+1,-1
         mm=r(k-1); nn=n(k)*r(k); mn=min(mm,nn); kk=r(k-2)*n(k-1)
         mat(1:mm,1:nn) => tt%u(k)%p
         call svd(mat,u,s,v,tol,rmax_,err,info)
         if(info.ne.0) then
            write(*,*) subnam,': svd info: ',info
            error stop
         endif
         rr=size(s); nrm=dnrm2(rr,s,1)
         if(nrm.ne.0.d0)then
            call dscal(rr,1.d0/nrm,s,1)
            lognrm=lognrm+dlog(nrm)
         end if
         do j=1,rr; call dscal(mm,s(j),u(1:mm,j),1); enddo
         !!
         call dgemm('n','n',kk,rr,mm,1.d0,tt%u(k-1)%p,kk,u,mm,0.d0,tmp,kk)
         if(r(k-1).ne.rr)then
            deallocate(tt%u(k-1)%p,tt%u(k)%p)
            r(k-1)=rr
            allocate(tt%u(k-1)%p(r(k-2),n(k-1),r(k-1)),&
                       tt%u(k)%p(r(k-1),n(k),r(k)))
         end if
         call dcopy(kk*rr,tmp,1,tt%u(k-1)%p,1)
         call dcopy(rr*nn,v,1,tt%u(k)%p,1)
         deallocate(u,s,v)
         nullify(mat)
      end do
      !!
      nrm=dnrm2(r(l-1)*n(l)*r(l),tt%u(l)%p,1)
      if(nrm.ne.0.d0)then
         call dscal(r(l-1)*n(l)*r(l),1.d0/nrm,tt%u(l)%p,1)
         lognrm=lognrm+dlog(nrm)
      endif
      !!
      lognrm=lognrm/(m-l+1)
      nrm=dexp(lognrm)
      do k=l,m
         call dscal(r(k-1)*n(k)*r(k),nrm,tt%u(k)%p,1)
      end do
      deallocate(tmp)
      tt% nu_core= 0
#  ifdef VERBOSE4
   print '("[+][dtt_svd] exit")'
#  endif
   end subroutine


   subroutine dtt_svd0(n,a,tt,tol_,rmax_)
   use mat_lib, only : svd
   implicit none
   type(dtt_tensor),intent(inout),target :: tt
   integer,intent(in) :: n(:)
   double precision,intent(in) :: a(*)
   double precision,intent(in),optional :: tol_
   integer,intent(in),optional :: rmax_
   character(len=*),parameter :: subnam='[!][TENSOR][dtt_svd0]'
   integer :: i,l,m,k,nn,mm,mn,info
   double precision,pointer,contiguous :: b(:,:),u(:,:),v(:,:),tmp(:)
   double precision,pointer :: s(:)
   integer,pointer :: r(:)
   double precision,external :: dnrm2
#  ifdef VERBOSE4
   print '("[+][dtt_svd0] entry")'
#  endif
   !!
      l=tt%l; m=l+size(n)-1;tt%m=m
      tt%n(l:m)=n
      r=>tt%r; r(l-1)=1;r(m)=1
      nn=product(n)
      if (8*nn.le.0) error stop subnam//': nn too big'
      allocate(tmp(nn),stat=info)
      if (info.ne.0) error stop subnam//': no memory'
      if (dnrm2(nn,a,1).eq.0.d0) then
         tt= dtt_tensor_zeros(n)
         return
      endif
      call dcopy(nn,a,1,tmp,1)
      !!
      do k=m,l+1,-1 ! right-to-left SVD
         mm=r(l-1)*product(tt%n(l:k-1)); nn=tt%n(k)*r(k); mn=min(mm,nn)
         b(1:mm,1:nn)=>tmp
         call svd(b,u,s,v,tol_,rmax_,info=info)
         if (info.ne.0) error stop subnam//': svd error '//str(info)
         r(k-1)=size(s)
         do i=1,r(k-1)
            call dscal(mm,s(i),u(1:mm,i),1)
         enddo
         if(associated(tt%u(k)%p)) deallocate(tt%u(k)%p)
         allocate(tt%u(k)%p(r(k-1),tt%n(k),r(k)))
         call dcopy(r(k-1)*tt%n(k)*r(k),v,1,tt%u(k)%p,1)
         call dcopy(mm*r(k-1),u,1,tmp,1)
         deallocate(u,s,v)
      end do
      if(associated(tt%u(l)%p))deallocate(tt%u(l)%p)
      allocate(tt%u(l)%p(r(l-1),tt%n(l),r(l)))
      call dcopy(r(l-1)*tt%n(l)*r(l),tmp,1,tt%u(l)%p,1)
      deallocate(tmp)
      tt% nu_core= 1
      tt% is_allocated= .true.
#  ifdef VERBOSE4
   print '("[+][dtt_svd0] exit")'
#  endif
   end subroutine


   !> Returns total size of all cores
   integer function dtt_size(tt) result(s)
   class(dtt_tensor), intent(in):: tt
   integer :: i

      s= 0
      do i=tt%l ,tt%m
         s= s + tt%r(i-1)*tt%n(i)*tt%r(i)
      enddo

   end function dtt_size


   !> Returns total size of full tensor (log10 of it)
   double precision function dtt_log10_fullsize(tt)
   class(dtt_tensor), intent(in):: tt

      dtt_log10_fullsize= sum(log10(dble(tt%n(tt%l:tt%m))))

   end function dtt_log10_fullsize


   !> Returns printable fullsize string
   function dtt_fullsize_string(tt) result(str)
   character(len=24) :: str
   class(dtt_tensor), intent(in):: tt
   double precision :: lgfs, mantis, fs
   integer :: expon
   character(len=1) :: ed

      lgfs= dtt_log10_fullsize(tt)
      expon= int(floor(lgfs))
      mantis= lgfs - dble(expon)
      if (expon.gt.7) then
         write(ed, '(i1)') int(log10(dble(expon))) + 1
         write(str, '(f14.7,"e+",i'//ed//')') 10**mantis, expon
      else
         write(ed, '(i1)') int(lgfs)+1
         fs= 10**lgfs
         if (fs - floor(fs) .le. ceiling(fs) - fs) then
            write(str, '(i'//ed//')') int(10**lgfs)
         else
            write(str, '(i'//ed//')') int(10**lgfs) + 1
         endif
      endif

   end function dtt_fullsize_string


   !> Compute the "grand sum" - sum of all elements in the matrix
   function dtt_grand_sum(this) result(gsum)
   class(dtt_tensor), intent(in) :: this
   double precision :: gsum
   double precision, allocatable :: s1(:,:), s2(:,:)
   integer :: k, l, m

      l = this%l; m = this% m
      s1 = sum(this% u(l)% p, 2)/this%n(1)
      do k=l+1,m
         s2 = thor_matmul(s1, sum(this% u(k)% p, 2))
         s1 = s2/this% n(k)
      enddo
      gsum = s1(1,1)

   end function dtt_grand_sum


   !> Compute Frobenius norm
   !!
   !! Translated from tt-toolbox @tt_tensor/norm.m
   function dtt_norm(this) result(norm)
   use mat_lib, only: qr => matlab_qr, normfro
   class(dtt_tensor), intent(in) :: this
   double precision :: norm
   !
   integer :: i, d, l, m, sh(2)
   integer, allocatable :: nn(:), rr(:)
   double precision, allocatable :: core0(:,:), core1(:,:)
   double precision, allocatable :: Q(:,:), R(:,:), nrm(:)

      l= this% l
      m= this% m
      d= m - l + 1
      allocate(nn(d), rr(0:d), nrm(d))
      nn(1:d)= this% n(l:m)
      rr(0:d)= this% r(l-1:m)

      core1= reshape(this% u(l)% p, [1, rr(0)*nn(1)*rr(1)])
      ortho_l2r: do i = 1,d-1
         core0= reshape(core1, [rr(i-1)*nn(i), rr(i)])
         call qr(Q, R, core0)
         nrm(i)= normfro(R)
         R = R/(nrm(i) + tiny(1d0))
         core1= reshape(this% u(i+1)% p, [rr(i),nn(i+1)*rr(i+1)])
         core1= thor_matmul(R, core1)
      enddo ortho_l2r
      nrm(d)= normfro(core1)
      norm= product(nrm(1:d))
      deallocate (nn, rr, core0, core1, Q, R, nrm)

   end function dtt_norm


   !> Compute bounded Frobenius norm (divided by total number of elements)
   function dtt_bounded_norm(this) result(norm)
   use mat_lib, only: qr => matlab_qr, normfro
   class(dtt_tensor), intent(in) :: this
   double precision :: norm
   !
   integer :: i, d, l, m, sh(2)
   integer, allocatable :: nn(:), rr(:)
   double precision, allocatable :: core0(:,:), core1(:,:)
   double precision, allocatable :: Q(:,:), R(:,:), nrm(:)

      if(this% nu_core.gt.0) then
         l= this% nu_core
         norm= normfro(this% u(l)% p)
      else
         l= this% l
         m= this% m
         d= m - l + 1
         allocate(nn(d), rr(0:d), nrm(d))
         nn(1:d)= this% n(l:m)
         rr(0:d)= this% r(l-1:m)

         core1= reshape(this% u(l)% p, [rr(0)*nn(1)*rr(1), 1])
         ortho_l2r: do i = 1,d-1
            core0= reshape(core1, [rr(i-1)*nn(i), rr(i)])
            call qr(Q, R, core0)
            nrm(i)= normfro(R)
            R = R/(nrm(i) + tiny(1d0))
            core1= reshape(this% u(i+1)% p, [rr(i),nn(i+1)*rr(i+1)])
            core1= thor_matmul(R, core1)
            sh = shape(Q)
            rr(i)= sh(2)
         enddo ortho_l2r
         nrm(d)= normfro(core1)
         norm= product(nrm(1:d)/sqrt(dble(nn(1:d))))
         deallocate (nn, rr, core0, core1, Q, R, nrm)
      endif

   end function dtt_bounded_norm


   !!
   !! SAY: PRINT/WRITE THOR TENSOR
   !!

   subroutine dtt_say(arg)
   implicit none
   class(dtt_tensor),intent(in) :: arg
   !
   character(len=1) d
   integer l, m
      l = arg%l; m = arg%m
      write(*,'(a,i4,a,i4,a,f6.2)') 'dtt[',l,':', m,']: rank ',erank(arg)
      write (d, '(I1)') int(log10(dble(maxval(arg%n(l:m))))) + 2
      if(all(arg%n(l+1:m).eq.arg%n(l)))then
        write(*,'(a,1x,i'//d//',a)') 'n: ',arg%n(l),' for all modes'
      else
        write(*,'(a,1x,2048i'//d//')') 'n: ',arg%n(l:m)
      end if
      if (maxval(arg%r).gt.maxval(arg%n)) &
        write (d, '(I1)') int(log10(dble(maxval(arg%r(l:m))))) + 2
      write(*,'("r: ",i'//d//')',advance='no') arg%r(l-1)
      if(all(arg%r(l+1:m-1).eq.arg%r(l)).and.m-l.gt.5) then
        write(*,'(1x,i'//d//'" ... "i'//d//')',advance='no') arg%r(l),arg%r(m-1)
      else
        write(*,'(2048i'//d//')',advance='no') arg%r(l:m-1)
      end if
      write(*,'(i'//d//')') arg%r(m)
      write(*,'("size in tt-format: ", I15)') arg% size()
      write(*,'("size in full format: ", A)') trim(arg% fullsize_string())
   end subroutine


   subroutine dttm_say(arg)
   implicit none
   class(dtt_matrix),intent(in) :: arg
   character(len=1) d
   integer :: l, m
      l = arg%l; m = arg%m
      write(*,'(a,i4,a,i4,a,f6.2)') 'dttm[',l,':', m,']: rank ',erank(arg)
      write (d, '(I1)') int(log10(dble(maxval(arg%n(l:m))))) + 2
      if(all(arg%q(l+1:m).eq.arg%q(l)))then
         write(*,'(a,1x,i'//d//',a)') 'q: ',arg%q(l),' for all modes'
      else
         write(*,'(a,1x,2048i'//d//')') 'q: ',arg%q(l:m)
      end if
      if(all(arg%s(l+1:m).eq.arg%s(l)))then
         write(*,'(a,1x,i'//d//',a)') 's: ',arg%s(l),' for all modes'
      else
         write(*,'(a,1x,2048i'//d//')') 's: ',arg%s(l:m)
      end if
      if (maxval(arg%r).gt.maxval(arg%n)) &
         write (d, '(I1)') int(log10(dble(maxval(arg%r(l:m))))) + 2
      write(*,'("r: ",i'//d//')',advance='no') arg%r(l-1)
      if(all(arg%r(l+1:m-1).eq.arg%r(l)).and.m-l.gt.5) then
         write(*,'(1x,i'//d//'" ... "i'//d//')',advance='no') arg%r(l),arg%r(m-1)
      else
         write(*,'(2048i'//d//')',advance='no') arg%r(l:m-1)
      end if
      write(*,'(i'//d//')') arg%r(m)
      write(*,'("size in tt-format: ", I15)') arg% size()
      write(*,'("size in full format: ", A)') trim(arg% fullsize_string())
   end subroutine


   subroutine dtt_sayfull(arg)
     implicit none
     class(dtt_tensor),intent(in) :: arg
     character(len=*),parameter :: subnam='dtt_sayfull'
     integer :: l,m,n,i,info
     double precision,allocatable :: a(:)
     l=arg%l; m=arg%m; if(l.gt.m)return
     n=arg%r(l-1)*product(arg%n(l:m))*arg%r(m); if(n.le.0)return
     if(n.gt.2**20)then;write(*,*)subnam,': too much to say';return;endif
     allocate(a(n),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
     a= arg% full()
     do i=1,n
       write(*,'(es22.14)') a(i)
     end do
     deallocate(a)
   end subroutine


   subroutine dtt_pprint(this, label_, line_len_)
   use matlab_struct_module, only: pprint_matrix3d
   implicit none
   type(dtt_tensor),intent(in) :: this
   character(len=*),intent(in),optional :: label_
   integer,intent(in),optional :: line_len_
   !
   integer :: k, line_len
   character(len=20) :: label

      label= "dtt_tensor: "; if (present(label_)) label= label_
      line_len= 120; if (present(line_len_)) line_len= line_len_
      print '(A)', label
      do k=this% l,this% m
         print '("-- core #",I4,": [",I6,", ",I6,", ",I6,"] ---")', &
               k, this% r(k-1), this% n(k), this% r(k)
         call pprint_matrix3d (this% u(k)% p, &
              flip23_=.true., frmt_='(ES10.3)', line_len_=line_len)
      enddo
      print '("=====")'

   end subroutine


   !> helper function to use as tt% pprint
   subroutine dtt_pprint_aux(this, label_, line_len_)
   implicit none
   class(dtt_tensor),intent(in) :: this
   character(len=*),intent(in),optional :: label_
   integer,intent(in),optional :: line_len_

      call dtt_pprint(this, label_, line_len_)

   end subroutine dtt_pprint_aux


   subroutine dttm_pprint(this, label_, line_len_)
   use matlab_struct_module, only: pprint_matrix3d
   implicit none
   class(dtt_matrix),intent(in) :: this
   character(len=*),intent(in),optional :: label_
   integer,intent(in),optional :: line_len_
   !
   integer :: k, line_len
   character(len=20) :: label

      label= "dtt_matrix: "; if (present(label_)) label= label_
      line_len= 120; if (present(line_len_)) line_len= line_len_
      print '(A)', label
      do k=this% l,this% m
         print '("-- core #"I4": ["I6", "I6", "I6", "I6"] ---")', &
               k, this% r(k-1), this% q(k), this% s(k), this% r(k)
         call pprint_matrix3d (this% u(k)% p, &
              flip23_=.true., frmt_='(ES10.3)', line_len_=line_len)
      enddo
      print '("=====")'

   end subroutine


   !> helper function to use as mA% pprint
   subroutine dttm_pprint_aux(this, label_, line_len_)
   implicit none
   class(dtt_matrix),intent(in) :: this
   character(len=*),intent(in),optional :: label_
   integer,intent(in),optional :: line_len_

      call dttm_pprint(this, label_, line_len_)

   end subroutine dttm_pprint_aux


   !!
   !! RANK
   !!

   double precision function dtt_rank(arg) result (r)
     implicit none
     class(dtt_tensor),intent(in) :: arg
     integer :: l,m,i,a,b,d
     l=arg%l;m=arg%m;d=m-l+1
     if(d.le.0)then;r=-1.d0;return;endif
     if(d.eq.1)then;r= 0.d0;return;endif
     r=0.d0
     do i=arg%l,arg%m
       r=r+arg%r(i-1)*arg%n(i)*arg%r(i)
     end do
     if(r.eq.0.d0)return
     b=arg%r(l-1)*arg%n(l) + arg%n(m)*arg%r(m)
     if(d.eq.2)then;r=r/b;return;endif
     a=sum(arg%n(l+1:m-1))
     r=(dsqrt(b*b+4.d0*a*r)-b)/(2.d0*a)
     return
   end function


   !!
   !! FULL: THOR-VECT -> ARRAY
   !!

   function dtt_full(arg,part,ind) result(a)
   ! A = [beta]*A + [alpha]*FULL(dtt_tensor), dtt_tensor = arg([l:m]), l,m=[part]
   ! A size r(l-1) *n(l)*...*n(m)* r(m)
   implicit none
   double precision,allocatable :: a(:)
   class(dtt_tensor),intent(in),target :: arg
   integer,intent(in),optional :: part(2),ind(:)
   !
   character(len=*),parameter :: subnam='[!][TENSOR][dtt_full]'
   double precision :: alp,bet
   type(pointd) :: p(0:1)
   integer,pointer :: r(:),n(:)
   integer :: l,m,na,nb,mem,mem1,rr,info,i,j,pp
   integer,allocatable :: ii(:),nn(:)
   double precision,allocatable :: q(:)
#  ifdef VERBOSE4
   print '("[+][dtt_full] entry")'
#  endif

      l=arg%l; m=arg%m; if (present(part)) then; l=part(1); m=part(2); endif
      if(l.gt.m) error stop subnam//": empty input"
      r=>arg%r; n=>arg%n

      allocate(ii(l:m),nn(l:m)); ii=0; nn(l:m)=arg%n(l:m)
      if (present(ind)) then
         do i=l,m
            ii(i)=ind(i-l+1)
            if(ii(i).lt.0 .or. ii(i).gt.n(i)) stop subnam//': invalid ind'
            if(ii(i).ne.0) nn(i)=1
         end do
      end if
      na= r(l-1) * product(nn(l:m)) * r(m)
      mem=na; rr=1
      do i=l,m
         mem1= r(l-1)*product(nn(l:i))*r(i)
         if (mem1.le.0) stop subnam//': oversized mem'
         mem= max(mem, mem1)
         rr=  max(rr,r(i-1)*r(i))
      end do

      allocate(p(0)%p(mem),p(1)%p(mem),q(rr),stat=info)
      if(info.ne.0) stop subnam//': could not allocate'
      if(ii(l).eq.0)then
         call dcopy(r(l-1)*n(l)*r(l), arg%u(l)%p, 1, p(0)%p, 1)
      else
         do j=1,r(l)
            call dcopy(r(l-1), arg%u(l)%p(1:r(l-1),ii(l),j), 1, &
                                   p(0)%p(1+(j-1)*r(l-1):j*r(l-1)), 1)
         enddo
      end if

      pp=0
      do i=l+1,m
         nb=r(l-1) * product(nn(l:i-1))
         if(ii(i).eq.0)then
            if(nb*n(i)*r(i).gt.mem) stop subnam//': nb-by-n-by-r > mem'
            call dgemm('n','n',nb,n(i)*r(i),r(i-1),1.d0, &
                       p(pp)%p,nb,arg%u(i)%p,r(i-1),0.d0,p(1-pp)%p,nb)
         else
            if(nb*r(i).gt.mem) stop subnam//': nb-by-r > mem'
            do j=1,r(i)
               call dcopy(r(i-1), arg%u(i)%p(1:r(i-1),ii(i),j), 1, &
                                           q(1+(j-1)*r(i-1):j*r(i-1)), 1)
            enddo
            call dgemm('n','n',nb,r(i),r(i-1),1.d0, p(pp)%p,nb, &
                        q,r(i-1),0.d0,p(1-pp)%p,nb)
         end if
         pp= 1 - pp
      end do
      a= p(pp)%p(1:na)
      deallocate(p(0)%p, p(1)%p, q, ii, nn)
#  ifdef VERBOSE4
   print '("[+][dtt_full] exit")'
#  endif
   end function dtt_full


   !> Convert tt-matrix to a regular matrix
   function dttm_full(ttm) result(A2)
   use matrix_util, only: unravel_index, ravel_multi_index
   implicit none
   double precision,allocatable :: A2(:,:)
   class(dtt_matrix),intent(IN) :: ttm
   !
   character(len=*),parameter :: subnam='[dttm_full]'
   double precision, allocatable :: a(:)
   integer :: i,j,k,ij,l,m,d,MA,NA
   integer,allocatable :: i_mi(:),j_mi(:),ij_mi(:)
   double precision,allocatable :: q(:)
#  ifdef VERBOSE4
   print '("[+]'//subnam//' entry")'
#  endif

      a = dtt_full(ttm)
      l = ttm%l; m = ttm%m
      d = m - l + 1
      MA = product(ttm%q(l:m))
      NA = product(ttm%s(l:m))
      allocate(A2(MA,NA),ij_mi(d))
      do j = 1, NA
         j_mi = unravel_index(j - 1, ttm%s(l:m))
         do i = 1, MA
            i_mi = unravel_index(i - 1, ttm%q(l:m))
            ij_mi = 0
            do k = 1,size(i_mi)
               ij_mi(k) = i_mi(k)
            enddo
            do k = 1,size(j_mi)
               ij_mi(k) = ij_mi(k) + j_mi(k)*ttm%q(k) !! TODO: check
            enddo
            ij = ravel_multi_index(ij_mi, ttm%n(l:m)) + 1
            A2(i,j) = a(ij)
         enddo
      enddo

      deallocate(a,i_mi,j_mi,ij_mi)

#  ifdef VERBOSE4
   print '("[+][dttm_full] exit")'
#  endif
   end function dttm_full


   !!
   !! T_IJK: ACCESS ELEMENT IJK (VALUE)
   !!

   double precision function dtt_ijk(arg,ind) result (a)
     implicit none
     type(dtt_tensor),intent(in) :: arg
     integer,intent(in) :: ind(:)
     character(len=*),parameter :: subnam='dtt_ijk'
     integer :: info,i,l,m,n(tt_size),r(0:tt_size)
     double precision,pointer :: x(:,:),y(:,:),z(:,:)
     l=arg%l;m=arg%m;n=arg%n;r=arg%r
     if(size(ind).lt.m-l+1) then;a=-2.d0;return;endif
     if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>n(l:m)))then;a=-3.d0;return;endif
     if(r(l-1).ne.1 .or. r(m).ne.1)then;a=-4.d0;return;endif
     allocate(x(r(m-1),r(m)),stat=info)
     if(info.ne.0)then;a=-1.d0;return;endif
     x=arg%u(m)%p(:,ind(m-l+1),:)
     do i=m-1,l,-1
       allocate(y(r(i-1),r(i)),z(r(i-1),r(m)),stat=info)
       if(info.ne.0)then;a=-2.d0;return;endif
       y=arg%u(i)%p(:,ind(i-l+1),:)
       z=thor_matmul(y,x)
       deallocate(x,y); x=>z; nullify(z)
     end do
     a=x(1,1)
     deallocate(x)
   end function


  !> Returns all the cores as a single flat 1D array
  function dtt_flat_cores(tt) result(flatcores)
  implicit none
  class(dtt_tensor), intent(in) :: tt
  double precision, allocatable :: flatcores(:)
  double precision, pointer :: p1(:)
  !
  integer :: k, sz, ofs

     allocate(flatcores(tt% size()))

     ! copy cores
     ofs = 1
     do k=tt% l, tt% m
        sz = tt% r(k-1) * tt% n(k) * tt% r(k)
        p1(1:sz)=> tt% u(k)% p
        flatcores(ofs:ofs+sz-1)= p1(1:sz)
        ofs = ofs + sz
     enddo

  end function dtt_flat_cores


  !> MATLAB-like rounding
  subroutine matlab_round(tt,tol_,rmax_)
  use mat_lib, only : d2submat, my_chop2, plain_svd
  use nan_lib, only : nan
  implicit none
  class(dtt_tensor),intent(inout),  target :: tt
  double precision, intent(in),   optional :: tol_
  integer,          intent(in),   optional :: rmax_
  !
  character(len=*),parameter               :: subnam='[!][TENSOR][matlab_round]'
  integer                                  :: l,m,k,i,j,info,nn,rr,mn,mm,kk,rmax
  integer, pointer                         :: r(:),n(:)
  double precision,pointer,contiguous      :: mat(:,:),u(:,:),s(:),v(:,:),tmp(:)
  double precision, allocatable            :: ss(:),uu(:,:),vv(:,:),nrm(:)
  double precision                         :: err,nrm0,tol
  double precision,external                :: dnrm2
  integer                                  :: lwork
  double precision, allocatable            :: work(:),tau(:),ru(:),core0(:)

     l = tt%l
     m = tt%m
     if(m.le.l)return

     tol= 1d-12; if (present(tol_)) tol = tol_/sqrt(dble(m-l+1))
     rmax= product(tt%n(l:m)); if (present(rmax_)) rmax= rmax_

     r => tt%r
     n => tt%n
     !
     nn = maxval(n(l:m))
     rr = maxval(r(l-1:m))
     !
     allocate(tmp(rr*nn*rr), stat=info)
     if(info.ne.0) error stop subnam//': no memory for allocating tmp'
     !
     allocate(nrm(l:m))
     !
     ! Begin: Orthogonalization from left-to-right
     !
     mm = 0
     do k=l,m-1
       mm = max(mm,r(k-1)*n(k))
     end do
     !
     lwork = -1
     allocate(work(1),tau(nn*rr), ru(rr*nn*rr),core0(rr*nn*rr), stat=info)
     if(info.ne.0) error stop subnam//': no memory for allocating work,tau,ru,core0'
     !
     call dgeqrf(mm,nn,core0,mm,tau,work,lwork,info)
     lwork = int(work(1)) + 1
     !
     deallocate(work)
     allocate(work(lwork),stat=info)
     !
     if(info.ne.0) error stop subnam//': no memory for allocating work'
     !
     do k=l,m-1
       mm = r(k-1)*n(k)
       nn = r(k)
       mn = min(mm,nn)
       kk = n(k+1)*r(k+1)
       !
       call dcopy(mm*nn,tt%u(k)%p,1,core0,1)
       call dgeqrf(mm,nn,core0,mm,tau,work,lwork,info) ! to calculate: ru in round.m
       if(info.ne.0) error stop subnam//': dgeqrf info: '//str(info)
       !
       do j=1,nn
         forall(i=1:min(j,mm))    ru(i+(j-1)*mn)=core0(i+(j-1)*mm)
         forall(i=min(j,mm)+1:mn) ru(i+(j-1)*mn)=0.d0
       end do
       !
       nrm(k+1)=dnrm2(mn*nn,ru,1)
       !
       if(nrm(k+1).ne.0.d0)then ! line 42-44 in round.m
         call dscal(mn*nn,1.d0/nrm(k+1),ru,1)
       endif
       !
       call dorgqr(mm,mn,mn,core0,mm,tau,work,lwork,info) ! to calculate: core0 in round.m
       if(info.ne.0) error stop subnam//': dorgqr info: '//str(info)

       ! line 49 in round.m
       call dcopy(mm*mn, core0,1,tt%u(k)%p,1)
       ! line 47 in round.m: core0 will be used again later to set core0=ru*core1 & core1=core0
       call dgemm('n','n',mn,kk,nn,1.d0,ru,mn,tt%u(k+1)%p,nn,0.d0,core0,mn)
       ! line 48 and required adjustments in case rank r(k) is truncated
       if(r(k).ne.mn)then
         call dcopy(mm*mn, tt%u(k)%p,1,ru,1)
         deallocate(tt%u(k)%p,tt%u(k+1)%p)
         r(k)=mn
         allocate(tt%u(k)%p(r(k-1),n(k),r(k)),tt%u(k+1)%p(r(k),n(k+1),r(k+1)))
         call dcopy(mm*mn, ru,1,tt%u(k)%p,1)
       end if
       ! line 47 in round.m: now set core1=ru*core1
       call dcopy(mn*kk, core0,1,tt%u(k+1)%p,1)
     end do
     !
     deallocate(work,tau,ru,core0)
     !
     ! End: Orthogonalization from left-to-right
     !
     do k=m,l+1,-1
       !
       mm = r(k-1)
       nn = n(k)*r(k)
       mn = min(mm,nn)
       kk = r(k-2)*n(k-1)
       !
       mat(1:mm,1:nn) => tt%u(k)%p
       !
       allocate(uu(mm,mn),vv(mn,nn),ss(mn))
       !
       call plain_svd(mat,uu,ss,vv,info) ! line 65 in round.m
       error stop subnam//': plain_svd dgesvd error: '//str(info)
       !
       rr = my_chop2(ss,tol*dnrm2(mn,ss,1))
       rr = min(rr,rmax)
       !
       allocate(u(mm,rr),s(rr),v(rr,nn))
       !
       call dcopy(mm*rr,uu,1,u,1)   ! line 80 in round.m
       call dcopy(rr,ss,1,s,1)      ! line 80 in round.m
       call d2submat(rr,nn,vv,mn,v) ! line 80 in round.m
       !
       deallocate(uu,vv,ss)
       !
       do j=1,rr
         call dscal(mm,s(j),u(1:mm,j),1) ! line 81 in round.m
       enddo
       ! line 83 in round.m
       call dgemm('n','n',kk,rr,mm,1.d0,tt%u(k-1)%p,kk,u,mm,0.d0,tmp,kk)
       ! line 82 in round.m and required adjustments in case rank r(i) is modified
       if(r(k-1).ne.rr)then
         deallocate(tt%u(k-1)%p,tt%u(k)%p)
         r(k-1)=rr
         allocate(tt%u(k-1)%p(r(k-2),n(k-1),r(k-1)),&
                    tt%u(k)%p(r(k-1),n(k),r(k)))
       end if
       !
       call dcopy(kk*rr,tmp,1,tt%u(k-1)%p,1) ! line 86 in round.m
       call dcopy(rr*nn,v,1,tt%u(k)%p,1)     ! line 85 in round.m
       !
       deallocate(u,s,v)
       nullify(mat)
       !
     end do
     ! line 95 in round.m
     nrm(1)=dnrm2(r(l-1)*n(l)*r(l),tt%u(l)%p,1)

     ! line 96-98 in round.m
     if(nrm(1).ne.0.d0)then
        call dscal(r(l-1)*n(l)*r(l),1.d0/nrm(1),tt%u(l)%p,1)
     endif
     ! line 102-103 in round.m
     nrm0 = sum(log(abs(nrm)))
     nrm0 = nrm0/dble(m-l+1)
     nrm0 = dexp(nrm0)
     ! line 105-107 in round.m
     do k=l,m
        call dscal(r(k-1)*n(k)*r(k),nrm0,tt%u(k)%p,1)
     end do
     !
     deallocate(tmp,nrm)
     !
  end subroutine matlab_round


  !!
  !! CORE2CELL FUNCTIONS
  !!
  function core2cell_dtt(tt) result(cc)
  use matlab_struct_module, only: array3d
  implicit none
  type(dtt_tensor), intent(in):: tt
  type(array3d) :: cc(tt% m - tt% l + 1)
  !
  integer :: i, d

      d = tt% m - tt% l + 1 ! number of dimensions of tt
      do i=1,d
         cc(i)% arr = tt% u(i)% p
      enddo

  end function core2cell_dtt


  function core2cell_dttm(tt) result(cc)
  use matlab_struct_module, only: array4d
  implicit none
  type(dtt_matrix), intent(in):: tt
  type(array4d) :: cc(tt% m - tt% l + 1)
  !
  integer :: i, d

      d = tt% m - tt% l + 1 ! number of dimensions of tt
      do i=1,d
         cc(i)% arr = tt% u4(i)% p
      enddo

  end function core2cell_dttm


  !> Generates a random tt-matrix with orthogonal cores
  function dttm_random_ortho(q, s, r_, rall_, lr_) result(ttm)
  type(dtt_matrix)    :: ttm
  integer, intent(in) :: q(:), s(:)
  integer, intent(in), optional :: r_(size(q)+1)
  integer, intent(in), optional :: rall_
  logical, intent(in), optional :: lr_
  !
  integer :: d, k

     d= size(q)
     if (size(s) /= d) error stop "[!][dttm_random_ortho]: q =/= s"
     ttm= dtt_random_ortho(q*s, r_, rall_, lr_)
     ttm% q(1:d)= q(1:d)
     ttm% s(1:d)= s(1:d)
     do k=1,d
        ttm% u4(k)% p(1:ttm% r(k-1), &
                      1:ttm% q(k),   &
                      1:ttm% s(k),   &
                      1:ttm% r(k))=> ttm% u(k)% p
     enddo

  end function dttm_random_ortho


  !> Generates a random tensor with orthogonal cores
  !!
  !! A wrapper for dtt_random_ortho_full with scalar arguments,
  !! that produces a random tensor with orth. cores of the size
  !! n x n x ... x n (d times) and core rands (1, r_, r_, ..., r_, 1)
  !! n : integer
  !!     size of the tensor in each dimension
  !!
  !! d : integer
  !!     number of dimensions
  !!
  !! r_: integer (optional, default = 1)
  !!     core ranks
  !!
  !! lr_: logical (optional, default = true)
  !!     is direction of orthogonalization left-to-right?
  !!
  function dtt_random_ortho_short(n,d,r_,lr_) result(tt)
  type(dtt_tensor)    :: tt
  integer, intent(in) :: n, d
  integer, intent(in), optional :: r_
  logical, intent(in), optional :: lr_
  !
  integer, allocatable :: nn(:), rr(:)
  logical :: lr

     allocate(nn(d), rr(d+1))
     nn = n
     rr = 1; if (present(r_)) rr(2:d) = r_
     lr = .true.; if (present(lr_)) lr = lr_
     tt = dtt_random_ortho_full(nn, r_=rr, lr_=lr)

  end function dtt_random_ortho_short


  !> Kronecker product of two dtt-tensors
  !!
  !! a, b : dtt_tensor
  !!     operands
  !!
  function dtt_kron(a, b) result(c)
  type(dtt_tensor), intent(in) :: a, b
  type(dtt_tensor) :: c
  !
  integer :: d1, d2, d, k, sh(3)
  integer, allocatable :: r(:), n(:)

     d1 = a%m - a%l + 1
     d2 = b%m - b%l + 1
     d = d1 + d2
     allocate(n(d), r(0:d))

     n(1:d1)   = a% n(1:d1)
     n(d1+1:d) = b% n(1:d2)
     r(0:d1)   = a% r(0:d1)
     r(d1:d)   = b% r(0:d2)

     c = dtt_tensor(n, r_=r)
     do k=1,d1
        c% u(k)% p= a% u(k)% p
     enddo
     do k=d1+1,d1+d2
        c% u(k)% p= b% u(k-d1)% p
     enddo

  end function dtt_kron


  !> Kronecker product of two dtt matrices
  !!
  !!   A, B : dtt_tensor
  !!       operands
  !!
  function dttm_kron(A, B) result(C)
  type(dtt_matrix), intent(in) :: A, B
  type(dtt_matrix) :: C
  !
  integer :: d1, d2, d, k
  integer, allocatable :: r(:), q(:), s(:)

     d1 = A%m - A%l + 1
     d2 = B%m - B%l + 1
     d = d1 + d2
     allocate(q(d), s(d), r(0:d))

     q(1:d1)   = A% q(1:d1)
     q(d1+1:d) = B% q(1:d2)

     s(1:d1)   = A% s(1:d1)
     s(d1+1:d) = B% s(1:d2)

     r(0:d1)   = A% r(0:d1)
     r(d1:d)   = B% r(0:d2)

     C = empty_dtt_matrix(q, s, r_=r)
     do k=1,d1
        C% u4(k)% p= A% u4(k)% p
        C% u(k)% p(1:r(k-1),1:s(k)*r(k),1:r(k))=> C% u4(k)% p
     enddo
     do k=d1+1,d1+d2
        C% u4(k)% p= B% u4(k-d1)% p
        C% u(k)% p(1:r(k-1),1:s(k)*r(k),1:r(k))=> C% u4(k)% p
     enddo

  end function dttm_kron


  !> Diagnoal dtt-matrix from a dtt-tensor
  function dttm_diag(x) result(A)
  type(dtt_tensor), intent(in) :: x
  type(dtt_matrix) :: A
  !
  integer :: i, k, l, m
  integer, allocatable :: n(:), r(:)

     l = x%l
     m = x%m
     n = x%n(l:m)
     r = x%r(l-1:m)

     A = dtt_matrix_zeros(n, n, r_=r, l_=l, m_=m)
     do k = l, m
        do i = 1, n(k)
           A% u4(k)% p(:,i,i,:)= x% u(k)% p(:,i,:)
        enddo
     enddo

  end function dttm_diag


  !!
  !! RANDOM TENSORS WITH ORTHOGONAL CORES
  !!

  !> Generates a random tensor with orthogonal cores
  !!
  !! Produces a random tensor with cores that are orthogonal matrices, i.e.:
  !!  - if dir>0 (left-to-right): for the j-th core of size r(j)*n(j)*r(j+1),
  !!    its matrix reshape into [r(j)*n(j)] X [r(j+1)] matrix is orthogonal:
  !!    C2' * C2 = I
  !!  - if dir<0 (right-to-left): for the j-th core of size r(j)*n(j)*r(j+1),
  !!    the transpose of its matrix reshape into [r(j)] X [n(j)*r(j+1)] matrix
  !!    is orthogonal:
  !!    C2 * C2' = I
  !!
  !! n(:) : array of integers
  !!    array of dimensions for this tensor
  !!
  !! r_(:) : array of integers (optional, default - array of 1s)
  !!    array of ranks for the cores
  !!
  !! rall_ : integer (optional, default - 1)
  !!    use constant rank for all ranks
  !!
  !! lr_: logical (optional, default = true)
  !!     is direction of orthogonalization left-to-right?
  !!
  !! ---------------------------
  !!
  !!  Adopted from: core/tt_rand.m in TT-Toolbox 2.2, 2009-2012
  !!
  !!---------------------------
  function dtt_random_ortho_full(n,r_,rall_,lr_) result(tt)
  use mat_lib, only: qr => matlab_qr, cumsum
  use rnd_lib, only: rndmat
  type(dtt_tensor)    :: tt
  integer, intent(in) :: n(:)
  integer, intent(in), optional :: r_(size(n)+1)
  integer, intent(in), optional :: rall_
  logical, intent(in), optional :: lr_
  !
  integer :: dir, d, i, sh(3)
  integer, allocatable :: r(:)
  double precision, allocatable :: cr(:,:),cr1(:,:),rv(:,:)
  logical :: lr

     ! we expect size(n) = size(r) - 1 > 1
     d = size(n)
     allocate(r(d+1))
     if (present(r_)) then
        r(:) = r_(:)
     else if (present(rall_)) then
        r = rall_
     else
        r = 1
     endif

     lr = .true.; if (present(lr_)) lr = lr_
     tt= dtt_tensor(n, r_=r)

     if (lr) then ! left-to-right orthogonality
        do i=1,d
           cr= rndmat(r(i)*n(i),r(i+1))
           call qr(cr1, rv, cr)
           r(i+1)= size(cr1, 2)
           sh= [r(i), n(i), r(i+1)]
           call realloc_to_shape(tt% u(i), sh)
           tt% u(i)% p= reshape(cr1, sh)
        enddo
     else         ! right-to-left orthogonality
        do i=d,1,-1
           cr= rndmat(n(i)*r(i+1), r(i))
           call qr(cr1, rv, cr)
           cr1= transpose(cr1)
           r(i)= size(cr1, 1)
           sh= [r(i), n(i), r(i+1)]
           call realloc_to_shape(tt% u(i), sh)
           tt% u(i)% p= reshape(cr1, sh)
        enddo
     endif

     ! cleanup
     deallocate(r,cr,cr1,rv)

  end function dtt_random_ortho_full

end module thor_lib

