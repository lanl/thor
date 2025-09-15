module lr_lib
 use ort_lib
 use default_lib
 implicit none

 interface lr
  module procedure lr_d2
 end interface

 contains
 subroutine lr_d2(a,u,b,tol,maxr,err)
  implicit none
  double precision,intent(in) :: a(:,:)
  double precision,pointer :: u(:,:),b(:,:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: maxr
  double precision,intent(out),optional :: err
  character(len=*),parameter :: subnam='lr_d2'
  double precision,allocatable :: x(:,:),y(:,:),z(:,:),g(:,:),v(:)
  integer,allocatable :: q(:)
  integer :: m,n,mn,i,j,jj,ij(2),rmax,r,info
  double precision :: er,nrm,nr1,er1,zz,xx
  logical :: done
  double precision,external :: dnrm2
  integer,external :: idamax

  m=size(a,1); n=size(a,2); mn=min(m,n)
  if(present(maxr))then;if(maxr.ge.0)then;rmax=min(maxr,mn);else;rmax=mn;endif; else;rmax=mn;end if

  allocate(x(m,rmax),y(n,rmax),z(m,n),v(n),q(n), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory!';stop;endif
  
  call dcopy(m*n,a,1,z,1)
  nrm=dnrm2(m*n,a,1); er=nrm; nr1=-1.d0
  r=0; done=(r.ge.rmax)
  
  do while(.not.done)
!   ij=maxloc(abs(z)); i=ij(1); j=ij(2)
!...OMP parallel do shared(z,q,v,m,n)
   do jj=1,n
    q(jj)=idamax(m,z(1,jj),1); v(jj)=abs(z(q(jj),jj))
   end do
!...OMP END parallel do
   j=maxloc(v,1); i=q(j); zz=z(i,j); er1=dabs(zz)
   if(r.eq.0)nr1=er1
   r=r+1
   call dcopy(m,z(1,j),1,x(1,r),1); xx=dnrm2(m,x(1,r),1); call dscal(m,1.d0/xx,x(1,r),1)
   call dcopy(n,z(i,1),m,y(1,r),1); call dscal(n,xx/zz,y(1,r),1)
   call dger(m,n,-1.d0,x(1,r),1,y(1,r),1,z,m)
   er=dnrm2(m*n,z,1)
   done=r.eq.rmax
   if(present(tol))done=done.or.er.le.tol*nrm

   !write(*,'(i3,a,2i6,a,2e10.3)')r,' @ ',i,j,' er ',er/nrm,er1/nr1
  end do
  
  allocate(u(m,r),b(r,n),g(r,r), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory!';stop;endif
  
  call ort0(x(:,1:r),out=u,mat=g)
  call dgemm('n','t',r,n,r,1.d0,g,r,y,n,0.d0,b,r)
!-  
!  call dcopy(m*n,a,1,z,1)
!  call dgemm('n','n',m,n,r,-1.d0,u,m,b,r,1.d0,z,m)
!  er=dnrm2(m*n,z,1)
!  write(*,*) 'truerr : ',er/nrm
!- 
  deallocate(x,y,z,v,q,g)
  if(present(err))err=er/nrm
 end subroutine


 subroutine d2_lrg(m,n,r,a,u,v,ind,jnd)
 implicit none
  integer,intent(in) :: m,n,r
  double precision,intent(in) :: a(m,n)
  double precision,intent(out) :: u(m,r),v(r,n)
  integer,intent(out) :: ind(r),jnd(r)
  character(len=*),parameter :: subnam='d2_lrg'
  double precision,allocatable :: e(:,:)
  integer :: info,p,ij,i,j
  integer,external :: idamax
  allocate(e(m,n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': alloc fail';stop;endif
  call dcopy(m*n,a,1,e,1)
  do p=1,r
   ij=idamax(m*n,e,1)
   j=(ij-1)/m+1; i=ij-(j-1)*m
   ind(p)=i; jnd(p)=j
   call dcopy(m,e(1,j),1,u(1,p),1)
   call dcopy(n,e(i,1),m,v(p,1),r)
   call dscal(n,1.d0/e(i,j),v(p,1),r)
   call dger(m,n,-1.d0,u(1,p),1,v(p,1),r,e,m)
  end do
  deallocate(e)
 end subroutine
 
 subroutine d2_lug(r,a,g)
 implicit none
  integer,intent(in) :: r
  double precision,intent(in) :: a(r,r)
  double precision,intent(out) :: g(r*r)
  character(len=*),parameter :: subnam='d2_lug'
  double precision,allocatable :: e(:,:)
  integer :: info,p,pos,i
  double precision :: aa
  allocate(e(r,r),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': alloc fail';stop;endif
  call dcopy(r*r,a,1,e,1)
  pos=0
  do p=1,r
   aa=1.d0/e(p,p)
   g(pos+1)=aa; pos=pos+1
   call dscal(r-p,aa,e(p,p+1),r)
   if(p.gt.1)then
    forall(i=1:p-1)g(pos+i)=e(i,p); pos=pos+p-1
    forall(i=1:p-1)g(pos+i)=e(p,i); pos=pos+p-1
   end if
   call dger(r-p,r-p,-1.d0,e(p+1,p),1,e(p,p+1),r,e(p+1,p+1),r)
  end do
  deallocate(e)
 end subroutine
 
 subroutine d2_lual(m,r,g,col,from)
 implicit none
  integer,intent(in) :: m,r
  double precision,intent(in) :: g(r*r)
  double precision,intent(inout) :: col(m,r)
  integer,intent(in),optional :: from
  character(len=*),parameter :: subnam='d2_lual'
  integer :: info,p,i,p0
  p0=default(1,from)
  do p=p0,r
   if(p.gt.1)then
    call dgemv('n',m,p-1,-1.d0,col,m,g(p**2-p+1:),1,1.d0,col(1,p),1)
   end if
   if (g(p**2).ne.0d0) call dscal(m,1.d0/g(p**2),col(1,p),1)
  end do
 end subroutine
 subroutine d2_luar(n,r,g,row,from)
 implicit none
  integer,intent(in) :: r,n
  double precision,intent(in) :: g(r*r)
  double precision,intent(inout) :: row(r,n)
  integer,intent(in),optional :: from
  character(len=*),parameter :: subnam='d2_luar'
  integer :: info,p,i,p0
  p0=default(1,from)
  do p=p0,r
   if(p.gt.1)then
    call dgemv('t',p-1,n,-1.d0,row,r,g(p**2-2*p+2:),1,1.d0,row(p,1),r)
   end if
  end do
 end subroutine
 
 subroutine d2_luil(m,r,g,col,opt)
 implicit none
  integer,intent(in) :: m,r
  double precision,intent(in)    :: g(r*r)
  double precision,intent(inout) :: col(m,r)
  integer,intent(in)             :: opt ! 1 or 2
  character(len=*),parameter :: subnam='d2_luil'
  integer :: p,k
  double precision :: x(r)
  if(.not. (opt.eq.1 .or. opt.eq.2))then; write(*,*)subnam,': unknown opt:',opt;stop;endif
   do p=r,1,-1  !ll
    if(p.gt.1)then
     call dgemv('n',m,p-1,-1.d0,col,m,g(p**2-(p-1)),1,1.d0,col(1,p),1)
    end if
    !call dscal(m,dsqrt(dabs(g(p**2)))/g(p**2),col(1,p),1)
    call dscal(m,1.d0/g(p**2),col(1,p),1)
   end do
  if(opt.eq.2)then 
   do p=1,r    !rl
    !call dscal(m,1.d0/dsqrt(dabs(g(p**2))),col(1,p),1)
    if(p.lt.r)then
     do k=p+1,r
      !x(k)=g(k**2-2*(k-1)+(p-1)) / dsqrt(dabs(g(k**2)))
      x(k)=g(k**2-2*(k-1)+(p-1)) 
     end do 
     call dgemv('n',m,r-p,-1.d0,col(1,p+1),m,x(p+1),1,1.d0,col(1,p),1)
    end if
   end do
  end if 
 end subroutine
 subroutine d2_luir(n,r,g,row,opt)
 implicit none
  integer,intent(in) :: n,r
  double precision,intent(in)    :: g(r*r)
  double precision,intent(inout) :: row(r,n)
  integer,intent(in)             :: opt
  character(len=*),parameter :: subnam='d2_luir'
  integer :: p,k
  double precision :: x(r)
  if(.not. (opt.eq.1 .or. opt.eq.2))then; write(*,*)subnam,': unknown opt:',opt;stop;endif
   do p=r,1,-1 ! rr
    if(p.gt.1)then
     call dgemv('t',p-1,n,-1.d0,row,r,g(p**2-2*(p-1)),1,1.d0,row(p,1),r)
    end if
    !call dscal(n,1.d0/dsqrt(dabs(g(p**2))),row(p,1),r)
   end do
  if(opt.eq.2)then 
   do p=1,r !lr
    !call dscal(n,dsqrt(dabs(g(p**2)))/g(p**2),row(p,1),r)
    call dscal(n,1.d0/g(p**2),row(p,1),r)
    if(p.lt.r)then
     do k=p+1,r
      !x(k)=g(k**2-(k-1)+(p-1)) * dsqrt(dabs(g(k**2)))/g(k**2)
      x(k)=g(k**2-(k-1)+(p-1)) / g(k**2)
     end do 
     call dgemv('t',r-p,n,-1.d0,row(p+1,1),r,x(p+1),1,1.d0,row(p,1),r)
    end if
   end do
  end if
 end subroutine
end module
 
