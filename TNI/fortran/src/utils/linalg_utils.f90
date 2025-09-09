!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file linalg_utils.f90
!! @author Ismael Boureima
!! @date  October 2023
!! @brief Data structures and functions inspired Matlab linear algebra toolbox

module linalg_module
use string_lib, only: str
implicit none
! --- module Interfaces  ------------------------------------
interface zeros
      module procedure zeros1, zeros2, zeros3, zeros4
end interface

interface zeros_like
      module procedure zeros1_like, zeros2_like, zeros3_like, zeros4_like
end interface

interface ones
      module procedure ones1, ones2, ones3, ones4
end interface

interface ones_like
      module procedure ones1_like, ones2_like, ones3_like, ones4_like
end interface

interface diag
      module procedure diag_ivect, diag_imat, diag_vect, diag_mat
end interface

interface kron
      module procedure kron_11, kron_22, kron_12, kron_21
end interface

interface colstack
      module procedure colstack_12, colstack_21, colstack_22
end interface

interface rowstack
      module procedure rowstack_12, rowstack_21
end interface

interface flip
      module procedure flip_vect, flip_col_row, flip_ivect, flip_icol_row
end interface

interface blkdiag
      module procedure blkdiag_22
end interface

interface lin_solve
      module procedure sgesv_21, sgesv_22, dgesv_21, dgesv_22
end interface

interface pprint
      module procedure pprint_f, pprint_d
end interface

contains ! --- module  ------------------------------------


!
! ASSERT
!

  subroutine assert(condition, message_)
  logical, intent(in)     :: condition
  character, intent(in), optional   :: message_
  if (present(message_)) then
     if (.not.condition) error stop 'Assertion failure: '//message_
  else
     if (.not.condition) error stop 'Assertion failure'
  endif
  end subroutine assert

!
! ONES
!
  function ones1(n) result(arr)
    integer, intent(in) :: n
    real(kind=kind(0.0d0)), allocatable :: arr(:)
    integer             :: info
    allocate(arr(n), stat=info)
    if (info.ne.0) error stop '[!][ones1]: Unable to allocate arr'//str(info)
    arr = 1d0
  end function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ones2(n,m) result(arr)
    integer, intent(in)                 :: n,m
    !real(kind=kind(0.0d0))             :: arr(n,m)
    !real(kind=kind(0.0d0)), allocatable :: arr(:,:)
    real(8), allocatable                :: arr(:,:)
    integer                :: info
    allocate(arr(n,m), stat=info)
    if (info.ne.0) error stop '[!][ones2]: Unable to allocate arr'//str(info)
    arr = 1.d0
    !print*, "shape(arr) = ", shape(arr)
    !print*, "arr = ", arr
  end function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ones3(n,m,o) result(arr)
    integer, intent(in) :: n,m,o
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:)
    integer             :: info
    allocate(arr(n,m,o), stat=info)
    if (info.ne.0) error stop '[!][ones3]: Unable to allocate arr'//str(info)
    arr = 1d0
  end function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ones4(n,m,o,p) result(arr)
    integer, intent(in) :: n,m,o,p
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:,:)
    integer             :: info
    allocate(arr(n,m,o,p), stat=info)
    if (info.ne.0) error stop '[!][ones3]: Unable to allocate arr'//str(info)
    arr = 1d0
  end function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ones1_like(v) result(arr)
    real(8), intent(in)                 :: v(:)
    real(kind=kind(0.0d0)), allocatable :: arr(:)
    arr = ones1(size(v,1))
  end function ones1_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ones2_like(a) result(arr)
    real(8), intent(in)                 :: a(:,:)
    real(kind=kind(0.0d0)), allocatable :: arr(:,:)
    arr = ones2(size(a,1), size(a,2))
  end function ones2_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ones3_like(a) result(arr)
    real(8), intent(in)                 :: a(:,:,:)
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:)
    arr = ones3(size(a,1), size(a,2), size(a,3))
  end function ones3_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ones4_like(a) result(arr)
    real(8), intent(in)                 :: a(:,:,:,:)
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:,:)
    arr = ones4(size(a,1), size(a,2), size(a,3), size(a,3))
  end function ones4_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! ZEROS
!

  function zeros1(n) result(arr)
    integer, intent(in) :: n
    real(kind=kind(0.0d0)), allocatable :: arr(:)
    integer             :: info
    allocate(arr(n), stat=info)
    if (info.ne.0) error stop '[!][zeros1]: Unable to allocate arr'//str(info)
    arr = 0d0
  end function zeros1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function zeros2(n,m) result(arr)
    integer, intent(in) :: n,m
    real(kind=kind(0.0d0)), allocatable :: arr(:,:)
    integer             :: info
    allocate(arr(n,m), stat=info)
    if (info.ne.0) error stop '[!][zeros2]: Unable to allocate arr'//str(info)
    arr = 0d0
  end function zeros2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function zeros3(n,m,o) result(arr)
    integer, intent(in) :: n,m,o
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:)
    integer             :: info
    allocate(arr(n,m,o), stat=info)
    if (info.ne.0) error stop '[!][zeros3]: Unable to allocate arr'//str(info)
    arr = 0d0
  end function zeros3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function zeros4(n,m,o,p) result(arr)
    integer, intent(in) :: n,m,o,p
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:,:)
    integer             :: info
    allocate(arr(n,m,o,p), stat=info)
    if (info.ne.0) error stop '[!][zeros4]: Unable to allocate arr'//str(info)
    arr = 0d0
  end function zeros4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function zeros1_like(v) result(arr)
    real(8), intent(in)                 :: v(:)
    real(kind=kind(0.0d0)), allocatable :: arr(:)
    arr = zeros1(size(v,1))
  end function zeros1_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function zeros2_like(a) result(arr)
    real(8), intent(in)                 :: a(:,:)
    real(kind=kind(0.0d0)), allocatable :: arr(:,:)
    arr = zeros2(size(a,1), size(a,2))
  end function zeros2_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function zeros3_like(a) result(arr)
    real(8), intent(in)                 :: a(:,:,:)
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:)
    arr = zeros3(size(a,1), size(a,2), size(a,3))
  end function zeros3_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function zeros4_like(a) result(arr)
    real(8), intent(in)                 :: a(:,:,:,:)
    real(kind=kind(0.0d0)), allocatable :: arr(:,:,:,:)
    arr = zeros4(size(a,1), size(a,2), size(a,3), size(a,3))
  end function zeros4_like
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! EYE
!

  function eye(n) result(res)
    integer, intent(in) :: n
    real(kind=kind(0.0d0)), allocatable :: res(:,:)
    integer :: i
    allocate(res(n,n))
    res = 0d0
    forall(i=1:n) res(i,i)=1d0
  end function
!
! DIAG
!
  function diag_imat(A) result(d)
    integer, intent(in) :: A(:,:)
    integer             :: d(size(A,1))
    integer             :: i
    do i = 1, size(A,1)
        d(i) = A(i,i)
    end do
  end function diag_imat
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function diag_ivect(v) result(d)
    integer, intent(in) :: v(:)
    integer             :: d(size(v,1), size(v,1))
    integer             :: i, m, info
    m = size(v,1)
    !allocate(d(m, m), stat=info)
    !if(info.ne.0)then;write(*,*)'[!][diag_ivect] problem allocating d',info;stop;endif
    d = 0
    do i = 1, m
        d(i,i) = v(i)
    end do
  end function diag_ivect
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function diag_mat(A) result(d)
    real(8), intent(in) :: A(:,:)
    real(8), allocatable :: d(:,:)
    integer             :: i,j, sz1, sz2, info
    sz1 = size(A,1) ;  sz2 = size(A,2)
    call assert( .not. ((sz1==1) .and. (sz2==1)), "size(A,1) = size(A,2) = 1")
    if (sz1 == 1) then
       allocate(d(sz2, sz2), stat=info)
       if(info.ne.0) error stop '[!][diag_mat] problem allocating d/1'//str(info)
       do i = 1, size(A,2)
          d(i,i) = A(1,i)
       end do
    else if (sz2 == 1) then
       allocate(d(sz1, sz1), stat=info)
       if(info.ne.0) error stop '[!][diag_mat] problem allocating d/2'//str(info)
       do i = 1, size(A,1)
          d(i,i) = A(i,1)
       end do
    else if (sz1==sz2) then
       allocate(d(sz1, 1), stat=info)
       if(info.ne.0) error stop '[!][diag_mat] problem allocating d/3'//str(info)
       do i = 1, size(A,1)
          d(i, 1) = A(i,i)
       end do
    else
       write(*,*) '[!][FEAUTURE]: diag(reanctangualr mat) is not supported yet'
       call exit(1)
    end if
  end function diag_mat
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function diag_vect(v) result(d)
    real(8), intent(in) :: v(:)
    real(8)             :: d(size(v,1), size(v,1))
    integer             :: i, m, info
    m = size(v,1)
    !allocate(d(m, m), stat=info)
    !if(info.ne.0)then;write(*,*)'[!][diag_vect] problem allocating d',info;stop;endif
    d = 0.d0
    do i = 1, m
        d(i,i) = v(i)
    end do
  end function diag_vect
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! KERON
!

  function kron_11(A,B) result(C)
    IMPLICIT NONE
    real(8), dimension (:)  , intent(in)  :: A, B
    real(8), dimension (:,:), allocatable :: C
    integer                               :: i, j, k, l, m, n, p, q
    integer                               :: info
    allocate(C(size(A,1)*size(B,1), 1), stat=info)
    if(info.ne.0) error stop '[!][kron_11]'//str(info)
    C = 0
    do i = 1,size(A,1)
      !do j = 1,size(A,2)
        n=(i-1)*size(B,1) + 1
        m=n+size(B,1) - 1
        p=1 !(j-1)*size(B,2) i+ 1
        q=p !+size(B,2) - 1
        C(n:m,1) = A(i)*B
      !enddo ! j
    enddo  ! i
  end function kron_11
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function kron_22(A,B) result(C)
    IMPLICIT NONE
    real(8), dimension (:,:), intent(in)  :: A, B
    real(8), dimension (:,:), allocatable :: C
    integer :: i = 0, j = 0, k = 0, l = 0
    integer :: m = 0, n = 0, p = 0, q = 0
    integer             :: info
    allocate(C(size(A,1)*size(B,1),size(A,2)*size(B,2)), stat=info)
    if(info.ne.0) error stop '[!][kron_22]'//str(info)
    C = 0
    do i = 1,size(A,1)
      do j = 1,size(A,2)
        n=(i-1)*size(B,1) + 1
        m=n+size(B,1) - 1
        p=(j-1)*size(B,2) + 1
        q=p+size(B,2) - 1
        C(n:m,p:q) = A(i,j)*B
      enddo ! j
    enddo  ! i
  end function kron_22
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function kron_12(A,B) result(C)
    IMPLICIT NONE
    real(8), dimension (:)  , intent(in)  :: A
    real(8), dimension (:,:), intent(in)  :: B
    real(8), dimension (:,:), allocatable :: C
    integer                               :: i, j, k, l, m, n, p, q
    integer                               :: info
    allocate(C(size(A,1)*size(B,1), size(B,2)), stat=info)
    if(info.ne.0) error stop '[!][kron_12]'//str(info)
    C = 0
    do i = 1,size(A,1)
        n=(i-1)*size(B,1) + 1
        m=n+size(B,1) - 1
        p=1 !(j-1)*size(B,2) + 1
        q=p+size(B,2) - 1
        C(n:m,p:q) = A(i)*B
    enddo  ! i
  end function kron_12
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function kron_21(A,B) result(C)
    IMPLICIT NONE
    real(8), dimension (:,:), intent(in)  :: A
    real(8), dimension (:)  , intent(in)  :: B
    real(8), dimension (:,:), allocatable :: C
    integer                               :: i, j, k, l, m, n, p, q
    integer                               :: info
    allocate(C(size(A,1)*size(B,1),size(A,2)), stat=info)
    if(info.ne.0) error stop '[!][kron_21]'//str(info)
    C = 0
    do i = 1,size(A,1)
      do j = 1,size(A,2)
        n=(i-1)*size(B,1) + 1
        m=n+size(B,1) - 1
        p=(j-1) + 1
        q=p !+size(B,2) - 1
        C(n:m,p) = A(i,j)*B
      enddo ! j
    enddo  ! i
  end function kron_21
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
! STACKING VECTORS AND MATRICES
!
  function colstack_21(A,B) result(C)
  implicit none
  double precision, intent(in)  :: A(:,:), B(:)
  double precision, allocatable :: C(:,:)
  integer :: info
     call assert(size(A,1)==size(B,1), "size(A,1) != size(B,1)")
     allocate(C(size(A,1), size(A,2)+1), stat=info)
     if(info.ne.0) error stop '[!][colstack_21]: Unable to allocate C'//str(info)
     C(:, 1:size(A,2)) = A
     C(:, size(A,2)+1) = B
  end function colstack_21


  function colstack_12(A,B) result(C)
  implicit none
  double precision, intent(in)  :: A(:), B(:,:)
  double precision, allocatable :: C(:,:)
  integer :: info
     call assert(size(A,1)==size(B,1), "size(A,1) != size(B,1)")
     allocate(C(size(A,1),size(B,2)+1), stat=info)
     if(info.ne.0) error stop '[!][colstack_12]: Unable to allocate C'//str(info)
     C(:, 2:size(B,2)+1) = B
     C(:, 1) = A
  end function colstack_12


  function colstack_22(A,B) result(C)
  implicit none
  double precision, intent(in)  :: A(:,:), B(:,:)
  double precision, allocatable :: C(:,:)
  integer :: info
     call assert(size(A,1)==size(B,1), "size(A,1) != size(B,1)")
     allocate(C(size(A,1), size(A,2)+size(B,2)), stat=info)
     if(info.ne.0) error stop '[!][colstack_22]: Unable to allocate C'//str(info)
     C(:, 1:size(A,2))  = A
     C(:, size(A,2)+1:) = B
  end function colstack_22


  function rowstack_21(A,B) result(C)
  implicit none
  double precision, intent(in)  :: A(:,:), B(:)
  double precision, allocatable :: C(:,:)
  integer :: info
     call assert(size(A,2)==size(B,1), "size(A,2) != size(B,1)")
     allocate(C(size(A,1)+1, size(A,2)), stat=info)
     if(info.ne.0) error stop '[!][rowstack_21]: Unable to allocate C'//str(info)
     C(1:size(A,1), :) = A
     C(size(A,1)+1, :) = B
  end function rowstack_21


  function rowstack_12(A,B) result(C)
  implicit none
  double precision, intent(in)  :: A(:), B(:,:)
  double precision, allocatable :: C(:,:)
  integer :: info
     call assert(size(B,2)==size(A,1), "size(B,2) != size(A,1)")
     allocate(C(size(B,1)+1, size(B,2)), stat=info)
     if(info.ne.0) error stop '[!][rowstack_12]: Unable to allocate C'//str(info)
     C(1, :) = A
     C(2:size(B,1)+1, :) = B
  end function rowstack_12


!
! SUM: Takes in Matrices A(i,j),B(k,l), assumed 2D, returns Direct sum C(i+k,j+l)
!
  function DirSum(A,B) result(C)
  double precision, intent(in)  :: A(:,:), B(:,:)
  double precision, allocatable :: C(:,:)
  integer :: p = 0, q = 0
     allocate(C(size(A,1)+size(B,1),size(A,2)+size(B,2)))
     C = 0
     p = size(A,1) + size(B,1)
     q = size(A,2) + size(B,2)
     C(1:size(A,1),1:size(A,2)) = A
     C(size(A,1)+1:p,size(A,2)+1:q) = B
  end function DirSum

  function VecDirSum(A,B) result(C)
  ! Takes 2 vectors, A(i),B(j), returns Direct Sum C(i+j)
  double precision, intent(in)  :: A(:), B(:)
  double precision, allocatable :: C(:)
     allocate(C(size(A)+size(B)))
     C = 0
     C(1:size(A)) = A
     C(size(A)+1:size(A)+size(B)) = B
  end function VecDirSUm


!
! FLIP
!
  function flip_vect(v) result(C)
  double precision, intent(in)  :: v(:)
  double precision, allocatable :: C(:)
  integer :: i,n, info
     n = size(v)
     allocate(C(n), stat=info)
     if(info.ne.0) error stop '[!][flip_vect()]: Unable to allocate C'//str(info)
     do i=1,n
        C(i) = v(n+1-i)
     end do
  end function flip_vect

  function flip_ivect(v) result(C)
  integer, intent(in)  :: v(:)
  integer, allocatable :: C(:)
  integer :: i,n, info
     n = size(v)
     allocate(C(n), stat=info)
     if(info.ne.0) error stop '[!][flip_ivect()]: Unable to allocate C'//str(info)
     do i=1,n
        C(i) = v(n+1-i)
     end do
  end function flip_ivect


  function flip_col_row(A) result(C)
  ! Takes in column or row vector v and returns v flipped
  double precision, intent(in)  :: A(:,:)
  double precision, allocatable :: C(:,:)
  integer :: i, m, n, info
     m = size(A,1); n = size(A,2)
     if(m>n.and.n.ne.1) error stop '[!][flip_col_row()]: m='//str(m)//' > n=' &
                                   //str(n)//' and n .ne.1 '//str(info)
     if(m<n.and.m.ne.1) error stop '[!][flip_col_row()]: m='//str(m)//' < n=' &
                                   //str(n)//' and m .ne.1 '//str(info)
     allocate(C(m,n), stat=info)
     if(info.ne.0) error stop '[!][flip_mat()]: Unable to allocate C '//str(info)
     C = 0
     if (m>n) then
        do i=1,n
           C(i,1) = A(n+1-i,1)
        end do
     else
        do i=1,n
           C(1,i) = A(1,n+1-i)
        end do
     end if
  end function flip_col_row


  function flip_icol_row(A) result(C)
  ! Takes in column or row vector v and returns v flipped
  integer, intent(in)  :: A(:,:)
  integer, allocatable :: C(:,:)
  integer :: i, m, n, info
     m = size(A,1); n = size(A,2)
     if(m>n.and.n.ne.1) error stop '[!][flip_icol_row()]: m='//str(m)//' > n=' &
                                   //str(n)//' and n .ne.1 '//str(info)
     if(m<n.and.m.ne.1) error stop '[!][flip_icol_row()]: m='//str(m)//' < n=' &
                                   //str(n)//' and m .ne.1 '//str(info)
     allocate(C(m,n), stat=info)
     if(info.ne.0) error stop '[!][flip_mat()]: Unable to allocate C'//str(info)
     C = 0
     if (m>n) then
        do i=1,m
           C(i,1) = A(n+1-i,1)
        end do
     else
        do i=1,n
           C(1,i) = A(1,n+1-i)
        end do
     end if
  end function flip_icol_row


!
! BLKDIAG
!

  function blkdiag_22(A,B) result(C)
  implicit none
  double precision, intent(in)  :: A(:,:), B(:,:)
  double precision, allocatable :: C(:,:)
  integer                               :: info
     allocate(C(size(A,1)+size(B,1), size(A,2)+size(B,2)), stat=info)
     if(info.ne.0) error stop '[!][blkdiag_22]: Unable to allocate C'//str(info)
     C(1:size(A,1),  1:size(A,2))  = A
     C(size(A,1)+1:, size(A,2)+1:) = B
  end function blkdiag_22


!
! SOLVE LINEAR SYSTEM: A @ x = B
!
  function sgesv_21(A,B)  result(X)
  external :: sgesv
  real(4), intent(in)  :: A(:,:), B(:)
  real(4), allocatable :: LU(:,:), X(:)
  integer, allocatable :: pivot(:) ! indices that define the permutation matrix
  integer  :: info, n,nrhs,lda,ldb
  character(*), parameter :: subname= '[sgesv_21]: '
     info=0;
     call assert(size(A,2)==size(B,1), subname//"size(A,2) != size(B,1)")
     n=size(A,2); lda=size(A,1); nrhs=1; ldb=size(B,1);
     allocate(X(ldb), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate X'//str(info)
     allocate(LU(lda, n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate LU'//str(info)
     allocate(pivot(n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate pivot'//str(info)
     X = B*1.0; LU = A*1.0
     call sgesv(n, nrhs, LU, lda, pivot, X, ldb, info)
     if(info.ne.0) error stop '[!]'//subname//'Error calling sgesv: '//str(info)
  end function sgesv_21


  function sgesv_22(A,B)  result(X)
  external :: sgesv
  real(4), intent(in)  :: A(:,:), B(:,:)
  real(4), allocatable :: LU(:,:), X(:,:)
  integer, allocatable :: pivot(:) ! indices that define the permutation matrix
  integer  :: info, n,nrhs,lda,ldb
  character(*), parameter :: subname= '[sgesv_22]: '
     info=0;
     call assert(size(A,2)==size(B,1), subname//"size(A,2) != size(B,1)")
     n=size(A,2); lda=size(A,1); nrhs=size(B,2); ldb=size(B,1);
     allocate(X(ldb, nrhs), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate X'//str(info)
     allocate(LU(lda, n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate LU'//str(info)
     allocate(pivot(n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate pivot'//str(info)
     X = B*1.0; LU = A*1.0
     call sgesv(n, nrhs, LU, lda, pivot, X, ldb, info)
     if(info.ne.0) error stop '[!]'//subname//'Error calling sgesv: '//str(info)
  end function sgesv_22


  function dgesv_21(A,B)  result(X)
  external :: dgesv
  double precision, intent(in)  :: A(:,:), B(:)
  double precision, allocatable :: LU(:,:), X(:)
  integer, allocatable :: pivot(:) ! indices that define the permutation matrix
  integer  :: info, n,nrhs,lda,ldb
  character(*), parameter :: subname= '[dgesv_21]: '
     info=0;
     call assert(size(A,2)==size(B,1), subname//"size(A,2) != size(B,1)")  ! Check A@B shape match
     n=size(A,2); lda=size(A,1); nrhs=1; ldb=size(B,1);
     allocate(X(ldb), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate X'//str(info)
     allocate(LU(lda, n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate LU'//str(info)
     allocate(pivot(n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate pivot'//str(info)
     X = B*1.0; LU = A*1.0
     call dgesv(n, nrhs, LU, lda, pivot, X, ldb, info)
     if(info.ne.0) error stop '[!]'//subname//'Error calling sgesv: '//str(info)
  end function dgesv_21


  function dgesv_22(A,B)  result(X)
  external :: dgesv
  double precision, intent(in)  :: A(:,:), B(:,:)
  double precision, allocatable :: LU(:,:), X(:,:)
  integer, allocatable :: pivot(:) ! indices that define the permutation matrix
  integer  :: info, n,nrhs,lda,ldb
  character(*), parameter :: subname= '[dgesv_22]: '
     info=0;
     call assert(size(A,2)==size(B,1), subname//"size(A,2) != size(B,1)")
     n=size(A,2); lda=size(A,1); nrhs=size(B,2); ldb=size(B,1);
     allocate(X(ldb, nrhs), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate X'//str(info)
     allocate(LU(lda, n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate LU'//str(info)
     allocate(pivot(n), stat=info)
     if(info.ne.0) error stop '[!]'//subname//'Unable to allocate pivot'//str(info)
     X = B*1.0; LU = A*1.0
     call dgesv(n, nrhs, LU, lda, pivot, X, ldb, info)
     if(info.ne.0) error stop '[!]'//subname//'Error calling sgesv: '//str(info)
  end function dgesv_22


!
!  PRETTY PRINTING
!

  subroutine pprint_f(arr)
    real, dimension(:,:), intent(in) :: arr
    integer :: i
    do i = 1, size(arr,dim=1)
      write (*,"(*(f8.4,:,','))") arr(i,:) ! unlimited format
    end do
  end subroutine pprint_f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pprint_d(arr)
    real(8), dimension(:,:), intent(in) :: arr
    integer :: i
    do i = 1, size(arr,dim=1)
      write (*,"(*(f8.4,:,','))") arr(i,:) ! unlimited format
    end do
  end subroutine pprint_d

end module linalg_module
