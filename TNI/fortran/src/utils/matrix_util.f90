!!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file matrix_util.f90
!! @author Ismael Djibrilla Boureima, Oleg Korobkin
!! @date  October 2023
!! @brief  General full matrix and tensor routines
!!
module matrix_util
  interface transp
     module procedure dtransp, ztransp, d1_transp, z1_transp
  end interface transp

  interface perm3d
     module procedure dperm_arr3d, dperm_3d, zperm_arr3d
  end interface perm3d

contains

  subroutine d1copy(A, B, N)
  double precision, intent(in):: A(*)
  double precision, intent(INOUT):: B(:)
  integer, intent(in) :: N
     B(1:N)= A(1:N)
  end subroutine d1copy


  subroutine eye(n,a)
    integer, intent(in) :: n
    double precision, intent(inout) :: a(n,n)
    integer i
    a(:,:) = 0d0
    do i = 1,n
       a(i,i) = 1d0
    end do
  end subroutine eye


  subroutine dqr(n,m, A, R)
    integer, intent(in):: n,m
    real(8), intent(inout) :: A(n,m), R(min(n,m),m)
    integer, parameter :: nb = 256
    integer :: lwork
    double precision :: tau(min(n,m))
    double precision :: work(nb*n)
    integer info
    integer rnew, k,j
    lwork = nb*n
    call dgeqrf(n, m, A, n, tau, work,lwork,info)
    if (info.ne.0) then
       print *, 'qr: dgeqrf failed'
    end if
    rnew = min(n,m)
    R(:,:)=0d0
    do j=1,m
       R(1:min(j,n),j)=A(1:min(j,n),j)
    end do
    call dorgqr(n,rnew,rnew,A,n,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'qr: dorgqr failed'
    end if

  end subroutine dqr


  subroutine zqr(n,m, A, R)
    integer, intent(in):: n,m
    complex(8), intent(inout) :: A(n,m), R(min(n,m),m)
    integer, parameter :: nb = 256
    integer :: lwork
    complex(8) :: tau(min(n,m))
    complex(8) :: work(nb*n)
    complex(8) :: ZERO
    parameter( ZERO=(0.0d0,0.0d0) )
    integer info
    integer rnew, k,j
    lwork = nb*n
    call zgeqrf(n, m, A, n, tau, work,lwork,info)
    if (info.ne.0) then
       print *, 'zqr: zgeqrf failed'
    end if
    rnew = min(n,m)
    R(:,:)=ZERO
    do j=1,m
       R(1:min(j,n),j)=A(1:min(j,n),j)
    end do
    call zungqr(n,rnew,rnew,A,n,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'zqr: zungqr failed'
    end if

  end subroutine zqr


  subroutine drow_add(m,n,k,A,B)
    integer, intent(in) :: m, n, k
    real(8), intent(inout) :: A(*)
    real(8), intent(in) :: B(*)
    real(8) swp(m)
    integer i

    do i=n,1,-1
       call dcopy(m, A(1+(i-1)*m), 1, swp, 1)
       call dcopy(m, swp, 1, A(1+(i-1)*(m+k)), 1)
       call dcopy(k, B(1+(i-1)*k), 1, A(m+1+(i-1)*(m+k)), 1)
    end do
  end subroutine drow_add

  subroutine zrow_add(m,n,k,A,B)
    integer, intent(in) :: m, n, k
    complex(8), intent(inout) :: A(*)
    complex(8), intent(in) :: B(*)
    complex(8) swp(m)
    integer i

    do i=n,1,-1
       call zcopy(m, A(1+(i-1)*m), 1, swp, 1)
       call zcopy(m, swp, 1, A(1+(i-1)*(m+k)), 1)
       call zcopy(k, B(1+(i-1)*k), 1, A(m+1+(i-1)*(m+k)), 1)
    end do
  end subroutine zrow_add

  subroutine drow_cut(m,n,k,A)
    integer, intent(in) :: m, n, k
    real(8), intent(inout) :: A(*)
    integer i

    do i=2,n
        call dcopy(k, A(1+(i-1)*m), 1, A(1+(i-1)*k), 1)
    end do
  end subroutine drow_cut

  subroutine zrow_cut(m,n,k,A)
    integer, intent(in) :: m, n, k
    complex(8), intent(inout) :: A(*)
    integer i

    do i=2,n
       call zcopy(k, A(1+(i-1)*m), 1, A(1+(i-1)*k), 1)
    end do
  end subroutine zrow_cut

  subroutine d1_transp(n, m, A, B)
    integer, intent(in):: n,m
    real(8), intent(in):: A(:)
    real(8), intent(inout), optional, target :: B(:)
    if ( present(B) ) then
        call dtransp(n, m, A, B)
    else
        call dtransp(n, m, A)
    end if
  end subroutine d1_transp

  subroutine dtransp(n, m, A, B)
    integer, intent(in):: n,m
    real(8), intent(in):: A(n,m)
    real(8), intent(inout), optional, target :: B(m,n)
    double precision, pointer :: C(:,:)
    integer i,j
    if ( present(B) ) then
       C => B
    else
       allocate(C(m,n))
    end if
    do i=1,n
       call dcopy(m, A(i,1), n, C(1,i),1)
    end do
    if ( .not. present(B) ) then
       call dcopy(n*m, C, 1, A, 1)
       deallocate(C)
    end if
  end subroutine dtransp

  subroutine ztransp(n, m, A, B)
    integer, intent(in):: n,m
    complex(8), intent(in):: A(n,m)
    complex(8), intent(inout), optional, target :: B(m,n)
    complex(8), pointer :: C(:,:)
    integer i,j
    if ( present(B) ) then
       C => B
    else
       allocate(C(m,n))
    end if
    do i=1,n
       call zcopy(m, A(i,1), n, C(1,i),1)
    end do
    if ( .not. present(B) ) then
       call zcopy(n*m, C, 1, A, 1)
       deallocate(C)
    end if
  end subroutine ztransp

  subroutine z1_transp(n, m, A, B)
    integer, intent(in):: n,m
    complex(8), intent(in):: A(:)
    complex(8), intent(inout), optional, target :: B(:)
    if ( present(B) ) then
        call ztransp(n, m, A, B)
    else
        call ztransp(n, m, A)
    end if

  end subroutine z1_transp

  subroutine dperm1324(n1,n2,n3,n4, A, B)
    integer, intent(in) :: n1,n2,n3,n4
    real(8), intent (in) :: A(n1,n2,n3,n4)
    real(8), intent (inout) :: B(n1,n3,n2,n4)

    integer i2,i3,i4

    do i4=1,n4
       do i3=1,n3
          do i2=1,n2
             call dcopy(n1, A(1,i2,i3,i4), 1, B(1,i3,i2,i4), 1)
          end do
       end do
    end do

  end subroutine dperm1324

  subroutine zperm1324(n1,n2,n3,n4, A, B)
    integer, intent(in) :: n1,n2,n3,n4
    complex(8), intent (in) :: A(n1,n2,n3,n4)
    complex(8), intent (inout) :: B(n1,n3,n2,n4)

    integer i2,i3,i4

    do i4=1,n4
       do i3=1,n3
          do i2=1,n2
             call zcopy(n1, A(1,i2,i3,i4), 1, B(1,i3,i2,i4), 1)
          end do
       end do
    end do

  end subroutine zperm1324

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! TT and svd stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_ps(d,r,n,ps)
    integer, intent(in) :: d,r(*),n(*)
    integer, intent(out) :: ps(*)
    integer i

    ps(1)=1;
    do i=1,d
       ps(i+1) = ps(i) + r(i)*n(i)*r(i+1)
    end do
  end subroutine compute_ps

  ! RELATIVE accuracy
    integer function my_chop3(n, s, eps)
    real(8), intent(in) :: s(*), eps
    integer, intent(in) :: n
    real(8) cursum, nrm
    integer i
    real(8) dnrm2

    nrm = dnrm2(n,s,1)
    nrm = (nrm*eps)*(nrm*eps);

    cursum = 0d0
    i = n;
    do while (i>0)
       cursum = cursum+(s(i)*s(i))
       if (cursum>nrm) then
          exit
       end if
       i = i-1
    end do

    my_chop3 = min(i+1, n)
  end function my_chop3


  recursive subroutine dperm_arr3d(A, pmt)
  double precision, allocatable, intent(INOUT) :: A(:,:,:)
  integer, intent(IN) :: pmt ! \in {123, 132, 213, 231, 312, 321}'
  !
  integer :: sh(3), n1,n2,n3,i,j,k
  double precision, allocatable :: aux(:), B(:,:,:)
     sh = shape(A)
     n1 = sh(1); n2 = sh(2); n3 = sh(3)
     select case (pmt)
     case (123)
        ! do nothing
     case (132) ! IX
        if (n2 == n3) then
           do j = 1,n2
              do k = j+1,n3
                 aux = A(:,k,j)
                 A(:,k,j) = A(:,j,k)
                 A(:,j,k) = aux
              enddo
           enddo
        else
           B = A
           deallocate(A)
           allocate(A(n1,n3,n2))
           forall(i=1:n1,j=1:n2,k=1:n3) A(i,k,j)= B(i,j,k)
           deallocate(B)
        endif
     case (213) ! XI
        if (n1 == n2) then
           do i = 1,n1
              do j = i+1,n2
                 aux = A(j,i,:)
                 A(j,i,:) = A(i,j,:)
                 A(i,j,:) = aux
              enddo
           enddo
        else
           B = A
           deallocate(A)
           allocate(A(n2,n1,n3))
           forall(i=1:n1,j=1:n2,k=1:n3) A(j,i,k)= B(i,j,k)
           deallocate(B)
        endif
     case (231) ! '//, == XI * IX
        call dperm_arr3d(A, 213)
        call dperm_arr3d(A, 132)
     case (312) ! ,\\' == IX * XI
        call dperm_arr3d(A, 132)
        call dperm_arr3d(A, 213)
     case (321) ! >|<
        if (n1 == n3) then
           do k = 1,n3
              do i = k+1,n1
                 aux = A(k,:,i)
                 A(k,:,i) = A(i,:,k)
                 A(i,:,k) = aux
              enddo
           enddo
        else
           B = A
           deallocate(A)
           allocate(A(n3,n2,n1))
           forall(i=1:n1,j=1:n2,k=1:n3) A(k,j,i)= B(i,j,k)
           deallocate(B)
        endif
     case default
        error stop "[!][dperm_arr3d]: invalid permutation"
     end select
  end subroutine dperm_arr3d


  recursive subroutine dperm_3d(A, n1, n2, n3, pmt)
  double precision, target, intent(INOUT) :: A(:)
  integer, intent(IN) :: n1,n2,n3
  integer, intent(IN) :: pmt ! \in {123, 132, 213, 231, 312, 321}'
  !
  integer :: i,j,k
  double precision, pointer :: p(:,:,:), q(:,:,:)
  double precision, allocatable :: aux(:), B(:,:,:)
     p(1:n1,1:n2,1:n3)=> A
     select case (pmt)
     case (123)
        ! do nothing
        p=> null()
        return
     case (132) ! IX
        q(1:n1,1:n3,1:n2)=> A
        if (n2 == n3) then
           do j = 1,n2
              do k = j+1,n3
                 aux = p(:,k,j)
                 p(:,k,j) = q(:,j,k)
                 q(:,j,k) = aux
              enddo
           enddo
        else
           B = p
           forall(i=1:n1,j=1:n2,k=1:n3) q(i,k,j)= B(i,j,k)
           deallocate(B)
        endif
     case (213) ! XI
        q(1:n2,1:n1,1:n3)=> A
        if (n1 == n2) then
           do i = 1,n1
              do j = i+1,n2
                 aux = p(j,i,:)
                 p(j,i,:) = q(i,j,:)
                 q(i,j,:) = aux
              enddo
           enddo
        else
           B = p
           forall(i=1:n1,j=1:n2,k=1:n3) q(j,i,k)= B(i,j,k)
           deallocate(B)
        endif
     case (231) ! '//, == XI * IX
        call dperm_3d(A, n1, n2, n3, 213)
        call dperm_3d(A, n2, n1, n3, 132)
     case (312) ! ,\\' == IX * XI
        call dperm_3d(A, n1, n2, n3, 132)
        call dperm_3d(A, n1, n3, n2, 213)
     case (321) ! >|<
        q(1:n3,1:n2,1:n1)=> A
        if (n1 == n3) then
           do k = 1,n3
              do i = k+1,n1
                 aux = p(k,:,i)
                 p(k,:,i) = q(i,:,k)
                 q(i,:,k) = aux
              enddo
           enddo
        else
           B = p
           forall(i=1:n1,j=1:n2,k=1:n3) q(k,j,i)= B(i,j,k)
           deallocate(B)
        endif
     case default
        error stop "[!][dperm_3d]: invalid permutation"
     end select
     p=> null()
     q=> null()
  end subroutine dperm_3d


  recursive subroutine zperm_arr3d(A, pmt)
  double complex, allocatable, intent(INOUT) :: A(:,:,:)
  integer, intent(IN) :: pmt ! \in {123, 132, 213, 231, 312, 321}'
  !
  integer :: sh(3), n1,n2,n3,i,j,k
  double complex, allocatable :: aux(:), B(:,:,:)
     sh = shape(A)
     n1 = sh(1); n2 = sh(2); n3 = sh(3)
     select case (pmt)
     case (123)
        ! do nothing
     case (132) ! IX
        if (n2 == n3) then
           do j = 1,n2
              do k = j+1,n3
                 aux = A(:,k,j)
                 A(:,k,j) = A(:,j,k)
                 A(:,j,k) = aux
              enddo
           enddo
        else
           B = A
           deallocate(A)
           allocate(A(n1,n3,n2))
           forall(i=1:n1,j=1:n2,k=1:n3) A(i,k,j)= B(i,j,k)
           deallocate(B)
        endif
     case (213) ! XI
        if (n1 == n2) then
           do i = 1,n1
              do j = i+1,n2
                 aux = A(j,i,:)
                 A(j,i,:) = A(i,j,:)
                 A(i,j,:) = aux
              enddo
           enddo
        else
           B = A
           deallocate(A)
           allocate(A(n2,n1,n3))
           forall(i=1:n1,j=1:n2,k=1:n3) A(j,i,k)= B(i,j,k)
           deallocate(B)
        endif
     case (231) ! '//, == XI * IX
        call zperm_arr3d(A, 213)
        call zperm_arr3d(A, 132)
     case (312) ! ,\\' == IX * XI
        call zperm_arr3d(A, 132)
        call zperm_arr3d(A, 213)
     case (321) ! >|<
        if (n1 == n3) then
           do k = 1,n3
              do i = k+1,n1
                 aux = A(k,:,i)
                 A(k,:,i) = A(i,:,k)
                 A(i,:,k) = aux
              enddo
           enddo
        else
           B = A
           deallocate(A)
           allocate(A(n3,n2,n1))
           forall(i=1:n1,j=1:n2,k=1:n3) A(k,j,i)= B(i,j,k)
           deallocate(B)
        endif
     case default
        error stop "[!][zperm_arr3d]: invalid permutation"
     end select
  end subroutine zperm_arr3d


  !> helper function: compute a binary represenation of integer n
  !!
  !! Examples:
  !!   2 -> {0, 1, 0, 0, 0, 0, 0, 0}
  !!  15 -> {1, 1, 1, 1, 0, 0, 0, 0}
  function binary_repr(n) result(ind)
#ifdef __NVCOMPILER
  ! NVIDIA compiler does not define shiftl and shiftr intrinsics (yet - O.K. 09/2025)
  use mat_lib, only: shiftr
#endif
  implicit none
  integer, intent(in) :: n
  integer, allocatable :: ind(:)
  !
  integer :: tmp, d, pos

     allocate(ind(128))
     ind = 0
     tmp = n
     d = 0
     do while (tmp > 0)
        d = d + 1
        ind(d) = iand(tmp, 1)
        tmp = shiftr(tmp, 1)
     enddo
     ind = ind(1:d)
     if (d.eq.0) ind = [0]

  end function binary_repr


  !> helper function to compute n from its binary representation
  !!
  !! Examples:
  !!  {0, 1, 0, 0, 0, 0, 0, 0} ->  2
  !!  {1, 1, 1, 1, 0, 0, 0, 0} -> 15
  function integer_from_binary(ind) result(n)
  implicit none
  integer, intent(in) :: ind(:)
  integer :: n
  !
  integer :: k, fof2

     n = 0
     fof2 = 1
     do k = 1, size(ind)
        n = n + ind(k)*fof2
        fof2 = fof2*2
     enddo

  end function integer_from_binary


  !> similar to numpy.unravel_index but with column-major ordering
  pure function unravel_index(k, ishape) result(ind)
  implicit none
  integer, intent(IN) :: k, ishape(:)
  integer, allocatable :: ind(:)
  !
  integer :: i, d, k1

     d = size(ishape)
     allocate(ind(d))
     k1 = k
     do i = 1, d
        ind(i) = mod(k1, ishape(i))
        k1 = k1/ishape(i)
     enddo

  end function unravel_index


  !> similar to numpy.ravel_multi_index but with column-major ordering
  pure function ravel_multi_index(ind, ishape) result(k)
  implicit none
  integer:: k
  integer, intent(IN) :: ind(:), ishape(:)
  !
  integer :: i, d

     d = size(ishape)
     k = ind(d)
     do i = d-1,1,-1
        k = k*ishape(i) + ind(i)
     enddo

  end function ravel_multi_index


  !> Transpose elements from a 2D matrix in a 1D array by interlaced indices
  !!
  !! Example:
  !!
  !!     i\j  |   0-0   0-1   1-0   1-1                  
  !!      ----|------------------------        0  1  4  5
  !!      0-0 |  0000  0001  0100  0101  --->  2  3  6  7
  !!      0-1 |  0010  0011  0110  0111        8  9 12 13
  !!      1-0 |  1000  1001  1100  1101       10 11 14 15
  !!      1-1 |  1010  1011  1110  1111 
  !!
  !!     --> { 0, 3, 8, 10, 1, 3, 9, 11, 4, 6, 12, 14, 5, 7, 13, 15 }
  !!           |
  !!           V
  !!         { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}
  !!
  function interlace_2d(A2) result(A)
  double precision, intent(IN) :: A2(:,:)
  double precision, allocatable :: A(:)
  !
  integer :: i, j, ij, N, sh(2), d
  integer, allocatable :: ii(:), jj(:), ind(:)

     sh = shape(A2)
     if (sh(1)/=sh(2)) &
        error stop "[!][interlace_2d] A2 is not square"
     N = sh(1)
     d = int(log(dble(N))/log(2d0) + 1d-11)
     allocate(A(N*N),ind(2*d))
     do i=1,N
        ind= 0
        ii= binary_repr(i-1)
        ind(2:2*size(ii):2)= ii
        do j=1,N
           jj= binary_repr(j-1)
           ind(1:2*size(jj)-1:2)= jj
           ij= integer_from_binary(ind) + 1
           A(ij)= A2(j, i)
        enddo
     enddo

  end function interlace_2d


  !> Same as interlace_2d but with 3D array
  function interlace_3d(A3) result(A)
  double precision, intent(IN) :: A3(:,:,:)
  double precision, allocatable :: A(:)
  !
  integer :: i, j, k, ijk, N, sh(3), d
  integer, allocatable :: ii(:), jj(:), kk(:), ind(:)

     sh = shape(A3)
     if (sh(1)/=sh(2) .or. sh(1)/=sh(3)) &
        error stop "[!][interlace_3d] A3 is not a cube"
     N = sh(1)
     d = int(log(dble(N))/log(2d0) + 1d-11)
     allocate(A(N*N*N),ind(3*d))
     do i=1,N
        ind= 0
        ii= binary_repr(i-1)
        ind(3:3*size(ii):3)= ii
        do j=1,N
           jj= binary_repr(j-1)
           ind(2:3*size(jj)-1:3)= jj
           do k=1,N
              kk= binary_repr(k-1)
              ind(1:3*size(kk)-2:3)= kk
              ijk= integer_from_binary(ind) + 1
              A(ijk)= A3(k, j, i)
           enddo
        enddo
     enddo

  end function interlace_3d

end module matrix_util


