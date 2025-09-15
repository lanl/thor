module rnd_lib
 implicit none
 interface random
  module procedure d1rnd,z1rnd, d2rnd,z2rnd, d3rnd,z3rnd, d4rnd,z4rnd, &
                   d5rnd,z5rnd, d6rnd,z6rnd, d7rnd,z7rnd
 end interface
 interface rndmat
  module procedure random_d1,random_d2,random_d3,random_d4
 end interface
 interface randn
  module procedure randn_d1, randn_d2, randn_d3, randn_d4, &
                   randn_d5, randn_d6
 end interface

 ! these functions are only available via interface
 private randn_d1, randn_d2, randn_d3, randn_d4, randn_d5, randn_d6, &
         d1rnd,z1rnd, d2rnd,z2rnd, d3rnd,z3rnd, d4rnd,z4rnd, &
         d5rnd,z5rnd, d6rnd,z6rnd, d7rnd,z7rnd, &
         random_d1,random_d2,random_d3,random_d4
contains

 subroutine arnd()
  implicit none
  integer s,clock,crate,cmax
  integer,allocatable :: seed(:)
  call random_seed(size=s)
  allocate(seed(s))
  call random_seed(get=seed)
  write(*,*) 'oldseed: ',seed
  call system_clock(clock,crate,cmax)
  !write(*,*)clock,crate,cmax
  seed(1)=clock
  call random_seed(put=seed)
  write(*,*) 'newseed: ',seed
 end subroutine

 double precision function drnd( )
  implicit none
  call random_number(drnd)
 return
 end function

 subroutine d1rnd(d)
  double precision,intent(out)  :: d(:)
  call random_number(d)
 end subroutine
 subroutine z1rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:)
  character(len=*),parameter :: subnam='z1rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine

 subroutine d2rnd(d)
  double precision,intent(out)  :: d(:,:)
  call random_number(d)
 end subroutine
 subroutine z2rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:)
  character(len=*),parameter :: subnam='z2rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine

 subroutine d3rnd(d)
  double precision,intent(out)  :: d(:,:,:)
  call random_number(d)
 end subroutine
 subroutine z3rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:,:)
  character(len=*),parameter :: subnam='z3rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine


 subroutine d4rnd(d)
  double precision,intent(out)  :: d(:,:,:,:)
  call random_number(d)
 end subroutine
 subroutine z4rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:,:,:)
  character(len=*),parameter :: subnam='z4rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine
 
 function random_d1(m) result(a)
  integer, intent(in) :: m
  double precision :: a(m)
  call random_number(a)
 end function random_d1
 function random_d2(m, n) result(a)
  integer, intent(in) :: m, n
  double precision :: a(m,n)
  call random_number(a)
 end function random_d2
 function random_d3(m, n, p) result(a)
  integer, intent(in) :: m, n, p
  double precision :: a(m,n,p)
  call random_number(a)
 end function random_d3
 function random_d4(m, n, p, q) result(a)
  integer, intent(in) :: m, n, p, q
  double precision :: a(m,n,p,q)
  call random_number(a)
 end function random_d4


 !> generate an array of normally distributed (mean=0, stddev=1) rn
 !!
 !! Uses Box-Muller transform
 !! Parameters:
 !!  * a(*): array to fill with random normally distributed numbers
 !!  * m   : array length
 !!
 subroutine d1randn(a, m)
 double precision, intent(INOUT) :: a(*)
 integer, intent(IN) :: m
 !
 integer :: m1, m2
 double precision, allocatable :: u1(:), u2(:)
 double precision, parameter :: twopi = 8d0*atan(1d0)
    
    m2 = m/2
    m1 = m - m2
    allocate(u1(m1),u2(m1))
    call random_number(u2)
    u1 = dsqrt(-2d0*dlog(u2))
    call random_number(u2)
    a(1:m1)   = u1*dcos(twopi*u2)
    a(m1+1:m) = u1(1:m2)*dsin(twopi*u2(1:m2))
    deallocate(u1,u2)

 end subroutine


 function randn_d1(m) result(a)
 integer, intent(IN) :: m
 double precision :: a(m)
    call d1randn(a, m)
 end function
 function randn_d2(m, n) result(a)
 integer, intent(IN) :: m, n
 double precision :: a(m, n)
    call d1randn(a, m*n)
 end function
 function randn_d3(m, n, p) result(a)
 integer, intent(IN) :: m, n, p
 double precision :: a(m, n, p)
    call d1randn(a, m*n*p)
 end function
 function randn_d4(m, n, p, q) result(a)
 integer, intent(IN) :: m, n, p, q
 double precision :: a(m, n, p, q)
    call d1randn(a, m*n*p*q)
 end function
 function randn_d5(m, n, p, q, s) result(a)
 integer, intent(IN) :: m, n, p, q, s
 double precision :: a(m, n, p, q, s)
    call d1randn(a, m*n*p*q*s)
 end function
 function randn_d6(m, n, p, q, s, t) result(a)
 integer, intent(IN) :: m, n, p, q, s, t
 double precision :: a(m, n, p, q, s, t)
    call d1randn(a, m*n*p*q*s*t)
 end function


 subroutine d5rnd(d)
  double precision,intent(out)  :: d(:,:,:,:,:)
  call random_number(d)
 end subroutine
 subroutine z5rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:,:,:,:)
  character(len=*),parameter :: subnam='z5rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine


 subroutine d6rnd(d)
  double precision,intent(out)  :: d(:,:,:,:,:,:)
  call random_number(d)
 end subroutine
 subroutine z6rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:,:,:,:,:)
  character(len=*),parameter :: subnam='z6rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine


 subroutine d7rnd(d)
  double precision,intent(out)  :: d(:,:,:,:,:,:,:)
  call random_number(d)
 end subroutine
 subroutine z7rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:,:,:,:,:,:)
  character(len=*),parameter :: subnam='z7rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine

 integer function irnd( maxi )
  integer,intent(in) :: maxi
  double precision :: d
  call random_number(d)
  irnd=int(d*maxi)+1
 return
 end function
 subroutine irand(maxi,ix)
  integer,intent(in)  :: maxi
  integer,intent(out) :: ix(:)
  integer :: i,n
  double precision,allocatable :: d(:)
  n=size(ix)
  allocate(d(n))
  call random_number(d)
  ix=int(d*maxi)+1
  deallocate(d)
 return
 end subroutine


 subroutine lottery2(npnt,m,n,wcol,wrow,points)
  implicit none
  integer,intent(in) :: npnt,m,n
  double precision,intent(in) :: wcol(m),wrow(n)
  integer,intent(out) :: points(npnt,2)
  character(len=*),parameter :: subnam='lottery2'
  double precision,allocatable :: pcol(:),prow(:),d(:,:)
  double precision :: scol,srow
  integer :: ipnt,i,j,info
  double precision,external :: dasum
  allocate(pcol(0:m),prow(0:n),d(npnt,2),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  scol=dasum(m,wcol,1); srow=dasum(n,wrow,1)
  if(dabs(srow) < 1d-16) srow=1d-16
  if(dabs(scol) < 1d-16) scol=1d-16
  pcol(0)=0.d0; do i=1,m;pcol(i)=pcol(i-1)+dabs(wcol(i))/scol;enddo
  prow(0)=0.d0; do j=1,n;prow(j)=prow(j-1)+dabs(wrow(j))/srow;enddo
  call random_number(d)
  do ipnt=1,npnt
   points(ipnt,1)=find_d(m+1,pcol,d(ipnt,1)); if(points(ipnt,1).gt.m)points(ipnt,1)=m
   points(ipnt,2)=find_d(n+1,prow,d(ipnt,2)); if(points(ipnt,2).gt.n)points(ipnt,2)=n
  end do
  deallocate(pcol,prow,d)
 end subroutine

 pure integer function find_d(n,x,y) result (pos)
  ! for sorted vector x(1) <= x(2) <= ... <= x(n) and value y find pos, s.t. x(pos) <= y < x(pos+1)
  implicit none
  integer,intent(in) :: n
  double precision,intent(in) :: x(n),y
  integer :: s,t,i
  logical :: key
  if(n.eq.0)then;pos=0;return;endif
  if(y.lt.x(1))then;pos=0;return;endif
  if(x(n).le.y)then;pos=n;return;endif
  s=1;t=n;pos=(t+s)/2
  do while(t-s.gt.1)
   if(y.lt.x(pos))then;t=pos;else;s=pos;end if
   pos=(s+t)/2
  enddo
  return
 end function


 !*****************************************
 !                                        *
 !  Initializes random number generator   *
 !                                        *
 !*****************************************
 subroutine init_random_seed(rseed)
 integer, intent(in):: rseed
 integer :: i, n, clock
 integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))
    clock= 0
    IF (rseed.EQ.0) call system_clock(count=clock)
    seed= clock + 37*(/ (i-1, i=rseed+1, rseed+n) /)
    call random_seed(put = seed)
    deallocate(seed)

 end subroutine init_random_seed


 subroutine init_random_seed_mpi(rseed)
 use mpi_f08
 integer, intent(in):: rseed
 !
 integer(MPI_INTEGER_KIND) :: rank, ierr
 integer, dimension(:), allocatable :: seed

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call init_random_seed(rseed + rank)

 end subroutine init_random_seed_mpi


end module
