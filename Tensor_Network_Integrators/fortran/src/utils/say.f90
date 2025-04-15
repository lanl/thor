module say_lib
 interface saynnz
  module procedure saynnz_d2,saynnz_d3,saynnz_z2,saynnz_z3
 end interface
 interface say
  module procedure say_d1,say_z1, say_d2,say_z2, say_d3, say_i2,say_i3
 end interface
contains
 subroutine say_d1(a)
  implicit none
  double precision,intent(in) :: a(:)
  character(len=*),parameter :: subnam='say_d1'
  integer :: n, i,j
  character(len=32) :: frm
  n=size(a,1)
  do i=1,n
   write(*,'(e22.12)') a(i)
  end do
  return
 end subroutine
 subroutine say_z1(a)
  implicit none
  double complex,intent(in) :: a(:)
  character(len=*),parameter :: subnam='say_z1'
  integer :: n, i,j
  character(len=32) :: frm
  n=size(a,1)
  do i=1,n
   write(*,'(e22.12,1x,e22.12)') a(i)
  end do
  return
 end subroutine
 subroutine say_d2(a)
  implicit none
  double precision,intent(in) :: a(:,:)
  character(len=*),parameter :: subnam='say_d2'
  integer :: m,n, i,j
  character(len=32) :: frm
  m=size(a,1); n=size(a,2)
  
  write(frm,'(a,i4,a)') '(a,',n,'(e9.3,1x),a)'
  !write(*,*)'frm: ',frm
  do i=1,m
   write(*,fmt=frm) '[ ',(a(i,j),j=1,n),']'
  end do
  return
 end subroutine
 subroutine say_z2(a)
  implicit none
  double complex,intent(in) :: a(:,:)
  character(len=*),parameter :: subnam='say_z2'
  integer :: m,n, i,j
  character(len=64) :: frm
  m=size(a,1); n=size(a,2)
  
  write(frm,'(a,i4,a)') '(a,',n,'(e9.3,1x,e9.3,2x),a)'
  !write(*,*)'frm: ',frm
  do i=1,m
   write(*,fmt=frm) '[ ',(a(i,j),j=1,n),']'
  end do
  return
 end subroutine
 subroutine say_d3(a)
  implicit none
  double precision,intent(in) :: a(:,:,:)
  integer :: m,n,p, i,j,k
  character(len=*),parameter :: subnam='say_d3'
  character(len=32) :: frm
  m=size(a,1); n=size(a,2); p=size(a,3)
  !write(*,*)subnam,':',m,n,p

  !write(frm,'(a,i3,a)') '(a,',n,'(e9.3,1x),a)'
  write(frm,'(a,i3,a)') '(a,',n,'(f5.3,1x),a)'
!  !write(*,*)'frm: ',frm
!  do k=1,p
!   do i=1,m
!    write(*,fmt=frm) '[ ',(a(i,j,k),j=1,n),']'
!   end do
!   write(*,*)
!  end do
  
  do i=1,m
   do k=1,p
    write(*,fmt=frm) '[ ',(a(i,j,k),j=1,n),']'
   end do
   write(*,*)
  end do
  return
 end subroutine
 
 subroutine say_i2(a)
  implicit none
  integer,intent(in) :: a(:,:)
  integer :: m,n, i,j
  character(len=*),parameter :: subnam='say_i2'
  character(len=32) :: frm
  m=size(a,1); n=size(a,2)
  write(frm,'(a,i3,a)') '(a,',n,'(i2,1x),a)'
  do i=1,m
   write(*,fmt=frm) '[ ',(a(i,j),j=1,n),']'
  end do
  return
 end subroutine
 subroutine say_i3(a)
  implicit none
  integer,intent(in) :: a(:,:,:)
  integer :: m,n,p, i,j,k
  character(len=*),parameter :: subnam='say_i3'
  character(len=32) :: frm
  m=size(a,1); n=size(a,2); p=size(a,3)
  write(frm,'(a,i3,a)') '(a,',n,'(i2,1x),a)'
  do k=1,p
   do i=1,m
    write(*,fmt=frm) '[ ',(a(i,j,k),j=1,n),']'
   end do
   write(*,*)
  end do
  return
 end subroutine
 
 subroutine saynnz_d2(x,tol)
  implicit none
  double precision,intent(in) :: x(:,:)
  double precision,intent(in),optional :: tol
  double precision :: eps,nrm
  integer :: i,j
  if(present(tol))then;eps=tol;else;eps=1.d-13;endif
  nrm=maxval(abs(x));eps=eps*nrm
  do j=1,size(x,2)
   do i=1,size(x,1)
    if(dabs(x(i,j)).gt.eps)write(*,'(a,2i6,a,e12.5)')'(',i,j,') ',x(i,j)
   enddo
  enddo
 end subroutine
 subroutine saynnz_d3(x,tol)
  implicit none
  double precision,intent(in) :: x(:,:,:)
  double precision,intent(in),optional :: tol
  double precision :: eps,nrm
  integer :: i,j,k
  if(present(tol))then;eps=tol;else;eps=1.d-13;endif
  nrm=maxval(abs(x));eps=eps*nrm
  do k=1,size(x,3)
   do j=1,size(x,2)
    do i=1,size(x,1)
    if(dabs(x(i,j,k)).gt.eps)write(*,'(a,3i6,a,e12.5)')'(',i,j,k,') ',x(i,j,k)
    enddo
   enddo
  enddo 
 end subroutine
 
 subroutine saynnz_z2(x,tol)
  implicit none
  double complex,intent(in) :: x(:,:)
  double precision,intent(in),optional :: tol
  double precision :: eps,nrm
  integer :: i,j
  if(present(tol))then;eps=tol;else;eps=1.d-13;endif
  nrm=maxval(cdabs(x));eps=eps*nrm
  do j=1,size(x,2)
   do i=1,size(x,1)
    if(cdabs(x(i,j)).gt.eps)write(*,'(a,2i6,a,2e12.5)')'(',i,j,') ',x(i,j)
   enddo
  enddo
 end subroutine
 subroutine saynnz_z3(x,tol)
  implicit none
  double complex,intent(in) :: x(:,:,:)
  double precision,intent(in),optional :: tol
  double precision :: eps,nrm
  integer :: i,j,k
  if(present(tol))then;eps=tol;else;eps=1.d-13;endif
  nrm=maxval(cdabs(x));eps=eps*nrm
  do k=1,size(x,3)
   do j=1,size(x,2)
    do i=1,size(x,1)
    if(cdabs(x(i,j,k)).gt.eps)write(*,'(a,3i6,a,2e12.5)')'(',i,j,k,') ',x(i,j,k)
    enddo
   enddo
  enddo 
 end subroutine
end module 

