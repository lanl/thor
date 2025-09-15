module ort_lib
 use nan_lib
 interface ort1
  module procedure ort1_d
 end interface
 interface ort0
  module procedure ort0_d,ort0_z
 end interface
 interface orto
  module procedure orto_d,orto_z
 end interface

 logical,private,parameter :: loud=.false.
 
 !private :: ort0_d, orto_d
contains
 subroutine ort0_d(inp,out,mat)
  implicit none
  double precision,intent(inout),target :: inp(:,:)
  double precision,intent(out),optional,target :: out(size(inp,1),size(inp,2))
  double precision,intent(out),optional :: mat(size(inp,2),size(inp,2))
  character(len=*),parameter  :: subnam='ort0_d'
  integer :: lwork,info,m,n,k,i,j
  double precision, allocatable :: work(:),tau(:)
  double precision,dimension(:,:),pointer :: y
!-for check
!  double precision :: nrm,err
!  double precision,allocatable :: tmp(:,:),mnn(:,:)
!  double precision,external :: dnrm2
!-      
  m=size(inp,1); n=size(inp,2); k=max(m,n)
  if(m.lt.n)then 
   !write(*,*)subnam,': less rows than columns: ',shape(inp)
   if(present(mat))then
    forall(i=1:m,j=1:n)mat(i,j)=inp(i,j)
    forall(i=m+1:n,j=1:n)mat(i,j)=0.d0
   end if
   if(present(out))then
    forall(i=1:m,j=1:n)out(i,j)=0.d0
    forall(i=1:m)out(i,i)=1.d0
   else 
    forall(i=1:m,j=1:n)inp(i,j)=0.d0
    forall(i=1:m)inp(i,i)=1.d0
   end if 
   return
  end if
!-for check- 
!  allocate(tmp(m,n), mnn(n,n))
!  call dcopy(m*n,inp,1,tmp,1); nrm=dnrm2(m*n,inp,1)
!-
  lwork=64*k
  allocate(tau(k),work(lwork),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  if(present(out))then;y=>out;call dcopy(m*n,inp,1,out,1);else;y=>inp;endif 

  call dgeqrf(m,n,y,m,tau,work,lwork,info)
  if(info.ne.0)then; write(*,*) subnam,': dgeqrf info: ',info; stop; endif
  if(nan(y))then;write(*,*)subnam,': dgeqrf returns NaNs, are you using Intel MKL 2011?';stop;endif
  
  if(present(mat))then
   forall(j=1:n) 
    forall(i=1:j)   mat(i,j)=y(i,j)
    forall(i=j+1:n) mat(i,j)=0.d0
   end forall 
  end if 
!-for check
!   forall(j=1:n) 
!    forall(i=1:j)   mnn(i,j)=y(i,j)
!    forall(i=j+1:n) mnn(i,j)=0.d0
!   end forall 
!-  
  call dorgqr(m,n,n,y,m,tau,work,lwork,info)
  if(info.ne.0)then; write(*,*) subnam,': dorgqr info: ',info; stop; end if
!-check-  
!  call dgemm('n','n',m,n,n,-1.d0,y,m,mnn,n,1.d0,tmp,m)
!  err=dnrm2(m*n,tmp,1)
!  write(*,*)subnam,': err: ',err/nrm
!-
  nullify(y)
  deallocate(tau,work)
 end subroutine

 subroutine ort0_z(inp,out,mat)
  implicit none
  double complex,intent(inout),target :: inp(:,:)
  double complex,intent(out),optional,target :: out(size(inp,1),size(inp,2))
  double complex,intent(out),optional :: mat(size(inp,2),size(inp,2))
  character(len=*),parameter  :: subnam='ort0_z'
  double complex,parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
  integer :: lwork,info,m,n,k,i,j
  double complex, allocatable :: work(:),tau(:)
  double complex,dimension(:,:),pointer :: y
!-for check
!  double precision :: nrm,err
!  double precision,allocatable :: tmp(:,:),mnn(:,:)
!  double precision,external :: dnrm2
!-      
  m=size(inp,1); n=size(inp,2); k=max(m,n)
  if(m.lt.n)then 
   !write(*,*)subnam,': less rows than columns: ',shape(inp)
   if(present(mat))then
    forall(i=1:m,j=1:n)mat(i,j)=inp(i,j)
    forall(i=m+1:n,j=1:n)mat(i,j)=zero
   end if
   if(present(out))then
    forall(i=1:m,j=1:n)out(i,j)=zero
    forall(i=1:m)out(i,i)=one
   else 
    forall(i=1:m,j=1:n)inp(i,j)=zero
    forall(i=1:m)inp(i,i)=one
   end if 
   return
  end if
!-for check- 
!  allocate(tmp(m,n), mnn(n,n))
!  call dcopy(m*n,inp,1,tmp,1); nrm=dnrm2(m*n,inp,1)
!-
  lwork=64*k
  allocate(tau(k),work(lwork),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  if(present(out))then;y=>out;call zcopy(m*n,inp,1,out,1);else;y=>inp;end if 
  
  call zgeqrf(m,n,y,m,tau,work,lwork,info)
  if(info.ne.0)then; write(*,*) subnam,': zgeqrf info: ',info; stop; end if
  if(nan(y))then;write(*,*)subnam,': zgeqrf returns NaNs, are you using Intel MKL 2011?';stop;endif

  if(present(mat))then
   forall(j=1:n) 
    forall(i=1:j)   mat(i,j)=y(i,j)
    forall(i=j+1:n) mat(i,j)=(0.d0,0.d0)
   end forall 
  end if 

!-for check
!   forall(j=1:n) 
!    forall(i=1:j)   mnn(i,j)=y(i,j)
!    forall(i=j+1:n) mnn(i,j)=0.d0
!   end forall 
!-  
  call zungqr(m,n,n,y,m,tau,work,lwork,info)
  if(info.ne.0)then; write(*,*) subnam,': zungqr info: ',info; stop; end if
!-check-  
!  call dgemm('n','n',m,n,n,-1.d0,y,m,mnn,n,1.d0,tmp,m)
!  err=dnrm2(m*n,tmp,1)
!  write(*,*)subnam,': err: ',err/nrm
!-
  nullify(y)
  deallocate(tau,work)
 end subroutine
   

 subroutine ort1_d(u,v,y,b)
  implicit none
  double precision,intent(in) :: u(:,:)
  double precision,intent(inout),target :: v(:)
  double precision,intent(out),optional,target :: y(size(v))
  double precision,intent(out),optional,target :: b(size(u,2)+1)
  character(len=*),parameter :: subnam='ort1_d'
  integer :: nu,nv,n,r,pass
  double precision,allocatable :: g(:)
  double precision,allocatable,target :: btmp(:)
  double precision,pointer :: bb(:),yy(:)
  double precision :: nrm1,nrm2
  logical :: reort
  double precision,external :: dnrm2
  
  r=size(u,2); nu=size(u,1); nv=size(v)
  if(nu.ne.nv)then
   write(*,*)subnam,': size mismatch: u[',nu,r,'] v[',nv,']'
   stop
  end if
  n=nu
  if(present(y))then
   yy=>y; y=v
  else
   yy=>v
  end if
  if(r.ge.n)then
   !write(*,'(2a,2i10)')subnam,': warning: n,r: ',n,r
   if(present(b))then
    call dgemv('t',n,n,1.d0,u,n,v,1, 0.d0,b(1:n),1)
    b(n+1:r+1)=0.d0
   end if 
   yy=0.d0
   return
  end if
  if(r.eq.0)then
   nrm2=dnrm2(n,v,1)
   if(nrm2.le.0.d0)then
    if(loud)write(*,'(2a,e10.3)')subnam,': small norm: ',nrm2
    yy=0.d0
   else
    yy=yy/nrm2
   end if
   if(present(b))b(1)=nrm2
   !write(*,*)subnam,': return since r=0'
   return
  end if

  if(present(b))then
   bb=>b
  else
   allocate(btmp(r+1)); bb=>btmp
  end if
  allocate(g(r))
      
  nrm1=dnrm2(n,yy,1)
  call dgemv('t',n,r,1.d0,u,n,yy,1, 0.d0,bb,1)
  call dgemv('n',n,r,-1.d0,u,n,bb,1, 1.d0,yy,1)
  nrm2=dnrm2(n,yy,1); reort=nrm2.lt.nrm1/2; pass=0
  
  do while(reort.and.pass.lt.3)
   call dgemv('t',n,r,1.d0,u,n,yy,1, 0.d0,g,1)
   call dgemv('n',n,r,-1.d0,u,n,g,1, 1.d0,yy,1)
   bb(1:r)=bb(1:r)+g
   nrm1=nrm2; nrm2=dnrm2(n,yy,1); reort=nrm2.lt.nrm1/2; pass=pass+1
  end do
  if(nrm2.lt.1.d-15)then
   if(loud)write(*,*)subnam,': small norm: ',nrm2
   yy=0.d0; bb(r+1)=0.d0
  else
   yy=yy/nrm2; bb(r+1)=nrm2
  end if
  
  nullify(bb,yy)
  deallocate(g)
  if(allocated(btmp))deallocate(btmp)
 end subroutine


 subroutine orto_d(u,v,out,mat)
  implicit none
  double precision,intent(in) :: u(:,:)
  double precision,intent(inout),target :: v(:,:)
  double precision,intent(out),optional,target :: out(size(v,1),size(v,2))
  double precision,intent(out),optional,target :: mat(size(u,2)+size(v,2),size(v,2))
  character(len=*),parameter :: subnam='orto_d'
  integer :: ru,rv,nu,nv,n,j, r,pass
  double precision,allocatable :: x(:),gu(:),gv(:)
  double precision,allocatable,target :: btmp(:,:)
  double precision,pointer :: bb(:,:),yy(:,:)
  double precision :: nrm1,nrm2
  logical :: reort
  double precision,external :: dnrm2
  
  !write(*,'(a,2(a,2i5),a)')subnam,': u[',size(u),'] v[',size(v),']'
  ru=size(u,2); rv=size(v,2); nu=size(u,1); nv=size(v,1)
  if(nu.ne.nv)then
   write(*,*)subnam,': size mismatch: u[',nu,ru,'] v[',nv,rv,']'
   stop
  end if
  n=nu
  if(present(out))then
   yy=>out; out=v
  else
   yy=>v
  end if
  if(present(mat))then
   bb=>mat
  else
   allocate(btmp(ru+rv,rv)); bb=>btmp
  end if
  allocate(x(n),gu(ru),gv(rv))
  
  do j=1,rv
   r=j-1 
   x=v(:,j); nrm1=dnrm2(n,x,1)
   call dgemv('t',n,ru,1.d0,u,n,x,1, 0.d0,gu,1)
   call dgemv('t',n,r,1.d0,yy,n,x,1, 0.d0,gv,1)
   call dgemv('n',n,ru,-1.d0,u,n,gu,1, 1.d0,x,1)
   call dgemv('n',n,r,-1.d0,yy,n,gv,1, 1.d0,x,1)
   bb(1:ru,j)=gu; bb(ru+1:ru+r,j)=gv(1:r)
   nrm2=dnrm2(n,x,1); reort=nrm2.lt.nrm1/2; pass=0
   
   do while(reort.and.pass.lt.3)
    call dgemv('t',n,ru,1.d0,u,n,x,1, 0.d0,gu,1)
    call dgemv('t',n,r,1.d0,yy,n,x,1, 0.d0,gv,1)
    call dgemv('n',n,ru,-1.d0,u,n,gu,1, 1.d0,x,1)
    call dgemv('n',n,r,-1.d0,yy,n,gv,1, 1.d0,x,1)
    bb(1:ru,j)=bb(1:ru,j)+gu; bb(ru+1:ru+r,j)=bb(ru+1:ru+r,j)+gv(1:r)
    nrm1=nrm2; nrm2=dnrm2(n,x,1); reort=nrm2.lt.nrm1/2; pass=pass+1
   end do
   if(nrm2.lt.1.d-15)then
    if(loud)write(*,*)subnam,': small norm: ',nrm2
    yy(:,j)=0.d0; bb(ru+j,j)=0.d0
   else
    yy(:,j)=x/nrm2; bb(ru+j,j)=nrm2
   end if
   bb(ru+j+1:ru+rv,j)=0.d0
  end do
  
  nullify(bb,yy)
  deallocate(x,gu,gv)
  if(allocated(btmp))deallocate(btmp)
 end subroutine
 subroutine orto_z(u,v,out,mat)
  implicit none
  double complex,intent(in) :: u(:,:)
  double complex,intent(inout),target :: v(:,:)
  double complex,intent(out),optional,target :: out(size(v,1),size(v,2))
  double complex,intent(out),optional,target :: mat(size(u,2)+size(v,2),size(v,2))
  character(len=*),parameter :: subnam='orto_d'
  double complex,parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
  integer :: ru,rv,nu,nv,n,j, r,pass
  double complex,allocatable :: x(:),gu(:),gv(:)
  double complex,allocatable,target :: btmp(:,:)
  double complex,pointer :: bb(:,:),yy(:,:)
  double precision :: nrm1,nrm2
  logical :: reort
  double precision,external :: dznrm2
  
  !write(*,'(a,2(a,2i5),a)')subnam,': u[',size(u),'] v[',size(v),']'
  ru=size(u,2); rv=size(v,2); nu=size(u,1); nv=size(v,1)
  if(nu.ne.nv)then
   write(*,*)subnam,': size mismatch: u[',nu,ru,'] v[',nv,rv,']'
   stop
  end if
  n=nu
  if(present(out))then
   yy=>out; out=v
  else
   yy=>v
  end if
  if(present(mat))then
   bb=>mat
  else
   allocate(btmp(ru+rv,rv)); bb=>btmp
  end if
  allocate(x(n),gu(ru),gv(rv))
  
  do j=1,rv
   r=j-1 
   x=v(:,j); nrm1=dznrm2(n,x,1)
   call zgemv('c',n,ru,one,u,n,x,1, zero,gu,1)
   call zgemv('c',n,r,one,yy,n,x,1, zero,gv,1)
   call zgemv('n',n,ru,-one,u,n,gu,1, one,x,1)
   call zgemv('n',n,r,-one,yy,n,gv,1, one,x,1)
   bb(1:ru,j)=gu; bb(ru+1:ru+r,j)=gv(1:r)
   nrm2=dznrm2(n,x,1); reort=nrm2.lt.nrm1/2; pass=0
   
   do while(reort.and.pass.lt.3)
    call zgemv('c',n,ru,one,u,n,x,1, zero,gu,1)
    call zgemv('c',n,r,one,yy,n,x,1, zero,gv,1)
    call zgemv('n',n,ru,-one,u,n,gu,1, one,x,1)
    call zgemv('n',n,r,-one,yy,n,gv,1, one,x,1)
    bb(1:ru,j)=bb(1:ru,j)+gu; bb(ru+1:ru+r,j)=bb(ru+1:ru+r,j)+gv(1:r)
    nrm1=nrm2; nrm2=dznrm2(n,x,1); reort=nrm2.lt.nrm1/2; pass=pass+1
   end do
   if(nrm2.lt.1.d-15)then
    if(loud)write(*,*)subnam,': small norm: ',nrm2
    yy(:,j)=zero; bb(ru+j,j)=zero
   else
    yy(:,j)=x/nrm2; bb(ru+j,j)=nrm2
   end if
   bb(ru+j+1:ru+rv,j)=zero
  end do
  
  nullify(bb,yy)
  deallocate(x,gu,gv)
  if(allocated(btmp))deallocate(btmp)
 end subroutine
end module
