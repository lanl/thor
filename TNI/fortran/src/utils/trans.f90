module trans_lib
 implicit none
!---------------------------------------------------------------
!     per = 1 -> 1 2 3 (No transpose, simple copy)
!     per = 2 -> 2 3 1
!     per = 3 -> 3 1 2
!     per = 4 -> 3 2 1
!     per = 5 -> 1 3 2
!     per = 6 -> 2 1 3
!---------------------------------------------------------------

 interface trans
  module procedure d3_trans,z3_trans,d2_trans,z2_trans,trans2d,trans2z
 end interface

 private :: d3_trans,d2_trans
contains

 subroutine d2_trans(a,b)
  implicit none
  double precision,intent(in) :: a(:,:)
  double precision,intent(out):: b(:,:)
  character(len=*),parameter :: subnam='d2_trans'
  integer :: i,j, m,n
  m=size(a,1); n=size(a,2)
  if(size(b,1).ne.n .or. size(b,2).ne.m)then
   write(*,*)subnam,': size mismatch: a[',shape(a),'] b[',shape(b),']'
   stop
  end if
  IF(m.gt.n)then
!$OMP PARALLEL DO SHARED(a,b)
   do j=1,n
    b(j,:)=a(:,j)
   end do
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SHARED(a,b)
   do i=1,m
    b(:,i)=a(i,:)
   end do
!$OMP END PARALLEL DO
  END IF
  return
 end subroutine
 subroutine z2_trans(a,b)
  implicit none
  double complex,intent(in) :: a(:,:)
  double complex,intent(out):: b(:,:)
  character(len=*),parameter :: subnam='z2_trans'
  integer :: i,j, m,n
  m=size(a,1); n=size(a,2)
  if(size(b,1).ne.n .or. size(b,2).ne.m)then
   write(*,*)subnam,': size mismatch: a[',shape(a),'] b[',shape(b),']'
   stop
  end if
  IF(m.gt.n)then
!$OMP PARALLEL DO SHARED(a,b)
   do j=1,n
    b(j,:)=a(:,j)
   end do
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SHARED(a,b)
   do i=1,m
    b(:,i)=a(i,:)
   end do
!$OMP END PARALLEL DO
  END IF
  return
 end subroutine

 subroutine d3_trans(p,a,b)
  implicit none
  integer,intent(in) :: p
  double precision,intent(in)  :: a(:,:,:)
  double precision,intent(out) :: b(:,:,:)
  character(len=*),parameter :: subnam='d3_trans'
  integer :: n(3),m(3),mm(3),i,j,k,q(3)
  if(p.lt.1 .or. p.gt.6)then;write(*,*)subnam,': illegal p: ',p;stop;endif
  m=shape(a); n=shape(b)
  mm=prm3(p,m)
  if(.not.all(mm==n))then
   write(*,'(2a,i1,a,3i8,a,3i8,a)')subnam,'(',p,') size mismatch: a[',m,'] b[',n,']'
   return
  end if

  select case(p)
   case(1)
!$OMP PARALLEL WORKSHARE SHARED(a,b)
    b=a
!$OMP END PARALLEL WORKSHARE
   case(2)
!$OMP PARALLEL DO SHARED(a,b)
    do k=1,m(3)
     do i=1,m(1)
      b(:,k,i)=a(i,:,k)
     end do
    end do
!$OMP END PARALLEL DO
   case(3)
!$OMP PARALLEL DO SHARED(a,b)
    do j=1,m(2)
     do i=1,m(1)
      b(:,i,j)=a(i,j,:)
     end do
    end do
!$OMP END PARALLEL DO
   case(4)
!$OMP PARALLEL DO SHARED(a,b)
    do j=1,m(2)
     do i=1,m(1)
      b(:,j,i)=a(i,j,:)
     end do
    end do
!$OMP END PARALLEL DO
   case(5)
!$OMP PARALLEL DO SHARED(a,b)
    do k=1,m(3)
     do j=1,m(2)
      b(:,k,j)=a(:,j,k)
     end do
    end do
!$OMP END PARALLEL DO
   case(6)
!$OMP PARALLEL DO SHARED(a,b)
    do k=1,m(3)
     do j=1,m(2)
      b(j,:,k)=a(:,j,k)
     end do
    end do
!$OMP END PARALLEL DO
   case default
    write(*,*)subnam,': illegal(!) p: ',p; stop
  end select
  return
 end subroutine
 subroutine z3_trans(p,a,b)
  implicit none
  integer,intent(in) :: p
  double complex,intent(in)  :: a(:,:,:)
  double complex,intent(out) :: b(:,:,:)
  character(len=*),parameter :: subnam='z3_trans'
  integer :: n(3),m(3),mm(3),i,j,k,q(3)
  if(p.lt.1 .or. p.gt.6)then;write(*,*)subnam,': illegal p: ',p;stop;endif
  m=shape(a); n=shape(b)
  mm=prm3(p,m)
  if(.not.all(mm==n))then
   write(*,'(2a,i1,a,3i8,a,3i8,a)')subnam,'(',p,') size mismatch: a[',m,'] b[',n,']'
   return
  end if

  select case(p)
   case(1)
!$OMP PARALLEL WORKSHARE SHARED(a,b)
    b=a
!$OMP END PARALLEL WORKSHARE
   case(2)
!$OMP PARALLEL DO SHARED(a,b)
    do k=1,m(3)
     do i=1,m(1)
      b(:,k,i)=a(i,:,k)
     end do
    end do
!$OMP END PARALLEL DO
   case(3)
!$OMP PARALLEL DO SHARED(a,b)
    do j=1,m(2)
     do i=1,m(1)
      b(:,i,j)=a(i,j,:)
     end do
    end do
!$OMP END PARALLEL DO
   case(4)
!$OMP PARALLEL DO SHARED(a,b)
    do j=1,m(2)
     do i=1,m(1)
      b(:,j,i)=a(i,j,:)
     end do
    end do
!$OMP END PARALLEL DO
   case(5)
!$OMP PARALLEL DO SHARED(a,b)
    do k=1,m(3)
     do j=1,m(2)
      b(:,k,j)=a(:,j,k)
     end do
    end do
!$OMP END PARALLEL DO
   case(6)
!$OMP PARALLEL DO SHARED(a,b)
    do k=1,m(3)
     do j=1,m(2)
      b(j,:,k)=a(:,j,k)
     end do
    end do
!$OMP END PARALLEL DO
   case default
    write(*,*)subnam,': illegal(!) p: ',p; stop
  end select
  return
 end subroutine

 subroutine trans2d(m,n,a,b)
  implicit none
  integer,intent(in) :: m,n
  double precision,intent(in) :: a(m,n)
  double precision,intent(out):: b(n,m)
  character(len=*),parameter :: subnam='trans2d'
  integer :: i,j
!$OMP PARALLEL DO SHARED(a,b)
  do i=1,m
    call dcopy(n, a(i,1), m, b(1,i), 1)
  end do
!$OMP END PARALLEL DO
 end subroutine
 subroutine trans2z(m,n,a,b)
  implicit none
  integer,intent(in) :: m,n
  double complex,intent(in) :: a(m,n)
  double complex,intent(out):: b(n,m)
  character(len=*),parameter :: subnam='trans2z'
  integer :: i,j
!$OMP PARALLEL DO SHARED(a,b)
  do i=1,m
    call zcopy(n, a(i,1), m, b(1,i), 1)
  end do
!$OMP END PARALLEL DO
 end subroutine

 pure function prm3(p,a)
  implicit none
  integer,intent(in) :: p
  integer,dimension(3),intent(in)  :: a
  integer,dimension(3) :: prm3
  integer,dimension(3,6),parameter :: prm = reshape( &
      source=(/   1,2,3,  2,3,1,   3,1,2,  3,2,1,  1,3,2 , 2,1,3 /), &
      shape=(/ 3 ,6 /) )
  prm3(:)=a(prm(:,p))
  return
 end function
end module
