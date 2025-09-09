module mat_lib
 use nan_lib
 implicit none

 interface matinv
  module procedure matinv_d,matinv_z
 end interface
 private :: matinv_d,matinv_z,matinv_svd_d,matinv_svd_z

 interface svd
  ! A = USV
  module procedure d_svd, z_svd
 end interface
 interface svd_old
  ! A = USV
  module procedure d_svd_old, z_svd_old
 end interface
 interface eye
  module procedure eye_d2,eye_z2,d2eye,z2eye
 end interface
 interface matlab_eye
  module procedure f_d2eye
 end interface
 interface laplace
  module procedure laplace_d2,laplace_z2
 end interface
 interface cumsum
  module procedure i_cumsum,d_cumsum
 end interface
 interface normfro
  module procedure norm1fro_d, norm2fro_d, norm3fro_d, norm4fro_d
 end interface
 interface normfro_bounded
  module procedure norm1fro_bounded_d, norm2fro_bounded_d, norm3fro_bounded_d
 end interface
 interface matlab_zeros
  module procedure dzeros_1d, dzeros_2d, dzeros_3d, dzeros_4d
 end interface
 interface matlab_diag
  module procedure diag_m2v, diag_v2m
 end interface
 interface matlab_qr
  module procedure matlab_qr_d
 end interface

contains
!MATINV
 subroutine matinv_d(a,ainv,alg,tol)
  implicit none
  double precision,intent(inout)  :: a(:,:)
  double precision,intent(out)    :: ainv(size(a,2),size(a,1))
  character(len=1),intent(in),optional :: alg
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_d'
  if(.not.present(alg) .or. size(a,1).ne.size(a,2))then
   call matinv_svd_d(a,ainv,tol)
  else
   select case(alg)
    case('s','S')
     call matinv_svd_d(a,ainv,tol)
    case('t','T')
     call matinv_lu_d(a,ainv)
    case default
     write(*,'(3a)')subnam,': unknown alg: ',alg
     stop
   end select
  end if
  return
 end subroutine
 subroutine matinv_z(a,ainv,alg,tol)
  implicit none
  double complex,intent(inout)  :: a(:,:)
  double complex,intent(out)    :: ainv(size(a,2),size(a,1))
  character(len=1),intent(in),optional :: alg
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_z'
  if(.not.present(alg) .or. size(a,1).ne.size(a,2))then
   call matinv_svd_z(a,ainv,tol)
  else
   select case(alg)
    case('s','S')
     call matinv_svd_z(a,ainv,tol)
    case('t','T')
     call matinv_lu_z(a,ainv)
    case default
     write(*,'(3a)')subnam,': unknown alg: ',alg
     stop
   end select
  end if
  return
 end subroutine

 subroutine matinv_svd_d(a,ainv,tol)
  implicit none
  double precision,intent(inout)  :: a(:,:)
  double precision,intent(out),optional :: ainv(size(a,2),size(a,1))
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_svd_d'
  double precision,allocatable :: b(:,:),u(:,:),v(:,:),s(:),work(:)
  integer :: r,lwork,info,i,ntrunc,m,n
  double precision :: si,s1,small
  character(len=180) :: str

  m=size(a,1); n=size(a,2)
  if(present(tol))then;small=tol;else;small=1.d-14;endif
  write(str,'(a,a,i8,1x,i8)')subnam,': m,n:',m,n
  !call plog(2,str)
  r=min(m,n); lwork=64*max(m,n)
  allocate(b(m,n),u(m,r),v(r,n),s(r),work(lwork))
  b=a
  call dgesvd('s','s',m,n,b,m,s,u,m,v,r, work,lwork,info)
  if(info.ne.0)then
   write(*,'(a,a,i10)')subnam,': dgesvd info: ',info
   stop
  end if

  s1=s(1)
  !write(str,'(a,80(e9.3,1x))')'sv:',(fmem(s+i)/s1,i=1,r)
  !call plog(1,str)

  ntrunc=0
  do i=1,r
   si=s(i)
    if(si.lt.small*s1)then
     s(i)=0.d0; ntrunc=ntrunc+1
    else
     s(i)=1.d0/s(i)
    end if
  end do

  forall(i=1:r)u(:,i)=u(:,i)*s(i)
  if(present(ainv))then
   call dgemm('t','t',n,m,r,1.d0,v,r,u,m,0.d0,ainv,n)
  else
   call dgemm('t','t',n,m,r,1.d0,v,r,u,m,0.d0,a,n)
  end if
  if(ntrunc.gt.0)then
   write(*,'(a,a,i5,a,i5)')subnam,': truncate: ',ntrunc,' of: ',r
   !call plog(1,str)
  end if

  deallocate(b,u,v,s,work)
  return
 end subroutine
 subroutine matinv_svd_z(a,ainv,tol)
  implicit none
  double complex,intent(in)  :: a(:,:)
  double complex,intent(out) :: ainv(size(a,2),size(a,1))
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_svd_z'
  double complex,allocatable   :: b(:,:),u(:,:),v(:,:),work(:)
  double precision,allocatable :: s(:),rwork(:)
  integer :: r,lwork,lrwork,info,i,ntrunc,m,n
  double precision :: small=1.d-14
  double precision :: si,s1
  character(len=180) :: str

  m=size(a,1); n=size(a,2)
  if(present(tol))small=tol
  write(str,'(a,a,i8,1x,i8)')subnam,': m,n:',m,n
  !call plog(2,str)
  r=min(m,n); lwork=64*max(m,n); lrwork=5*min(m,n)
  allocate(b(m,n),u(m,r),v(n,r),s(r),work(lwork),rwork(lrwork))
  b=a
  call zgesvd('s','s',m,n,b,m,s,u,m,v,n,work,lwork,rwork,info)
  if(info.ne.0)then
   write(*,'(a,a,i10)')subnam,': zgesvd info: ',info
   stop
  end if

  s1=s(1)
  !write(str,'(a,80(e9.3,1x))')'sv:',(fmem(s+i)/s1,i=1,r)
  !call plog(1,str)

  ntrunc=0
  do i=1,r
   si=s(i)
    if(si.lt.small*s1)then
     s(i)=0.d0; ntrunc=ntrunc+1
    else
     s(i)=1.d0/s(i)
    end if
  end do

  forall(i=1:r)u(:,i)=u(:,i)*s(i)
  call zgemm('c','c',n,m,r,(1.d0,0.d0),v,r,u,m,(0.d0,0.d0),ainv,n)
  if(ntrunc.gt.0)then
   write(str,'(a,a,i5,a,i5)')subnam,': truncate: ',ntrunc,' of: ',r
   !call plog(1,str)
  end if

  deallocate(b,u,v,s,work,rwork)
  return
 end subroutine

 subroutine matinv_lu_d(a,ainv)
  implicit none
  double precision,intent(inout),target  :: a(:,:)
  double precision,intent(out),optional,target :: ainv(size(a,2),size(a,1))
  character(len=*),parameter   :: subnam='matinv_lu_d'
  double precision,allocatable :: work(:)
  double precision,dimension(:,:),pointer :: aa
  integer :: lwork,m,n,info
  integer,allocatable :: piv(:)
  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   write(*,*)subnam,': matrix not square: ',m,n
   stop
  end if
  if(present(ainv))then
   call dcopy(m*n,a,1,ainv,1)
   aa=>ainv
  else
   aa=>a
  end if

  lwork=64*n
  allocate(work(lwork),piv(n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call dgetrf(n,n,aa,n,piv,info)
  if(info.ne.0)then; write(*,*)subnam,': dgertf: info: ',info; stop; end if
  call dgetri(n,aa,n,piv,work,lwork,info)
  if(info.ne.0)then; write(*,*)subnam,': dgerti: info: ',info; stop; end if

  nullify(aa)
  deallocate(work,piv)
  return
 end subroutine
 subroutine matinv_lu_z(a,ainv)
  implicit none
  double complex,intent(inout),target  :: a(:,:)
  double complex,intent(out),optional,target :: ainv(size(a,2),size(a,1))
  character(len=*),parameter   :: subnam='matinv_lu_z'
  double complex,allocatable :: work(:)
  double complex,dimension(:,:),pointer :: aa
  integer :: lwork,m,n,info
  integer,allocatable :: piv(:)
  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   write(*,*)subnam,': matrix not square: ',m,n
   stop
  end if
  if(present(ainv))then
   call zcopy(m*n,a,1,ainv,1)
   aa=>ainv
  else
   aa=>a
  end if

  lwork=64*n
  allocate(work(lwork),piv(n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif

  call zgetrf(n,n,aa,n,piv,info)
  if(info.ne.0)then; write(*,*)subnam,': dgertf: info: ',info; stop; end if
  call zgetri(n,aa,n,piv,work,lwork,info)
  if(info.ne.0)then; write(*,*)subnam,': dgerti: info: ',info; stop; end if

  deallocate(work,piv)
  return
 end subroutine

!EYE
 subroutine eye_d2(a)
  implicit none
  double precision,intent(inout) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  a=0.d0
  forall(i=1:min(m,n))a(i,i)=1.d0
  return
 end subroutine
 subroutine eye_z2(a)
  implicit none
  double complex,intent(inout) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  a=(0.d0,0.d0)
  forall(i=1:min(m,n))a(i,i)=(1.d0,0.d0)
  return
 end subroutine
 subroutine d2eye(a,n)
  implicit none
  double precision,intent(out) :: a(n,n)
  integer,intent(in) :: n
  integer :: i
  call dscal(n*n,0.d0,a,1)
  forall(i=1:n)a(i,i)=1.d0
 end subroutine
 subroutine z2eye(n,a)
  implicit none
  integer,intent(in) :: n
  double complex,intent(out) :: a(n,n)
  integer :: i
  call zdscal(n*n,0.d0,a,1)
  forall(i=1:n)a(i,i)=(1.d0,0.d0)
 end subroutine

 function f_d2eye(n) result(A)
  implicit none
  integer,intent(in) :: n
  double precision,allocatable :: A(:,:)
  allocate(A(n,n))
  call d2eye(A,n)
 end function f_d2eye

!LAPLACE
 subroutine laplace_d2(a)
  implicit none
  double precision,intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  call dscal(m*n,0.d0,a,1)
  forall(i=1:min(m,n))a(i,i)=2.d0
  forall(i=1:min(m-1,n))a(i+1,i)=-1.d0
  forall(i=1:min(m,n-1))a(i,i+1)=-1.d0
  return
 end subroutine
 subroutine laplace_z2(a)
  implicit none
  double complex,intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  call zdscal(m*n,0.d0,a,1)
  forall(i=1:min(m,n))a(i,i)=(2.d0,0.d0)
  forall(i=1:min(m-1,n))a(i+1,i)=-(1.d0,0.d0)
  forall(i=1:min(m,n-1))a(i,i+1)=-(1.d0,0.d0)
  return
 end subroutine


!SUBMAT
 subroutine d2submat(m,n,a,lda,b)
  implicit none
  integer,intent(in) :: m,n,lda
  double precision,intent(in) :: a(lda,n)
  double precision,intent(out) :: b(m,n)
  integer :: i,j
  if(m.gt.lda)then;write(*,*)'d2submat: lda,m: ',lda,m;stop;endif
  if(m.eq.lda)then;call dcopy(m*n,a,1,b,1);return;endif
  forall(i=1:m,j=1:n)b(i,j)=a(i,j)
 end subroutine
 subroutine z2submat(m,n,a,lda,b)
  implicit none
  integer,intent(in) :: m,n,lda
  double complex,intent(in) :: a(lda,n)
  double complex,intent(out) :: b(m,n)
  integer :: i,j
  if(m.gt.lda)then;write(*,*)'d2submat: lda,m: ',lda,m;stop;endif
  if(m.eq.lda)then;call zcopy(m*n,a,1,b,1);return;endif
  forall(i=1:m,j=1:n)b(i,j)=a(i,j)
 end subroutine
 subroutine d2subset(m,n,r,ind,jnd,a,b)
  implicit none
  integer,intent(in) :: m,n,r
  integer,intent(in) :: ind(r),jnd(r)
  double precision,intent(in) :: a(m,n)
  double precision,intent(out) :: b(r,r)
  integer :: i,j
  forall(i=1:r,j=1:r)b(i,j)=a(ind(i),jnd(j))
 end subroutine
 subroutine z2subset(m,n,r,ind,jnd,a,b)
  implicit none
  integer,intent(in) :: m,n,r
  integer,intent(in) :: ind(r),jnd(r)
  double complex,intent(in) :: a(m,n)
  double complex,intent(out) :: b(r,r)
  integer :: i,j
  forall(i=1:r,j=1:r)b(i,j)=a(ind(i),jnd(j))
 end subroutine

! SVD
 subroutine d_svd(a,u,s,v,tol,rmax,err,info)
  implicit none
  double precision,intent(in) :: a(:,:)
  double precision,pointer,intent(out) :: u(:,:),s(:),v(:,:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='d_svd'
  double precision,allocatable :: b(:,:),work(:),ss(:),uu(:,:),vv(:,:)
  integer,allocatable :: iwork(:)
  integer :: m,n,mn,mx,lwork,ierr,r,i,j
  double precision,external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); mx=max(m,n); lwork=-1
  allocate(work(1), b(m,n),uu(m,mn),vv(mn,n),ss(mn))
  ! call dgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,iwork,ierr)
  call dgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,ierr)
  lwork=int(work(1))+1
  deallocate(work)
  allocate(work(lwork),iwork(8*mn), stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call dcopy(m*n,a,1,b,1)
  call dgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,ierr)
  !call dgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,iwork,ierr)
  if(present(info))info=ierr
  deallocate(b,work,iwork)
  if(ierr.ne.0)then
   write(*,*)subnam,': dgesvd info: ',ierr
   if(ierr.lt.0)stop
   if(nan(a))then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   else
    write(*,*) subnam,': min/max element of input: ',minval(a),maxval(a)
   end if
   u=>null(); v=>null(); s=>null()
  else
   r=chop(ss,tol,rmax,err)
   if(present(err)) err= err/max(dnrm2(mn,ss,1),1d-16)
   allocate(u(m,r),s(r),v(r,n))
   call dcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   call d2submat(r,n,vv,mn,v)
  end if
  deallocate(uu,vv,ss)
 end subroutine
 subroutine z_svd(a,u,s,v,tol,rmax,err,info)
  implicit none
  double complex,intent(in) :: a(:,:)
  double complex,pointer :: u(:,:),v(:,:)
  double precision,pointer :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='z_svd'
  double complex,allocatable :: b(:,:),work(:),uu(:,:),vv(:,:)
  double precision,allocatable :: ss(:),rwork(:)
  integer,allocatable :: iwork(:)
  integer :: m,n,mn,mx,lwork,lrwork,ierr,r,i,j
  double precision,external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); mx=max(m,n); lwork=-1; lrwork=mn*max(5*(mn+1),2*(mn+mx+1))
  allocate(work(1),rwork(lrwork),b(m,n), uu(m,mn),vv(mn,n),ss(mn))
  ! call zgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,iwork,ierr)
  call zgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,ierr)
  lwork=int(real(work(1)))+1
  deallocate(work)
  allocate(work(lwork),iwork(8*mn), stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call zcopy(m*n,a,1,b,1)
  call zgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,ierr)
  !call zgesdd('s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,iwork,ierr)
  if(present(info))info=ierr
  deallocate(b,work,rwork,iwork)
  if(ierr.ne.0)then
   write(*,*)subnam,': zgesvd info: ',ierr
   if(ierr.lt.0)stop
   if (nan(a)) then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   end if
   u=>null(); v=>null(); s=>null()
  else
   r=chop(ss,tol,rmax,err)
   if(present(err))err=err/dnrm2(mn,ss,1)
   allocate(u(m,r),s(r),v(r,n))
   call zcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   call z2submat(r,n,vv,mn,v)
  end if
  deallocate(uu,vv,ss)
 end subroutine

 integer function chop(s,tol,rmax,err) result (r)
  implicit none
  double precision,intent(in) :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  double precision :: nrm,er,er2,bound
  double precision,external :: dnrm2
  r=size(s); er2=0.d0
  if(present(rmax))then
   if(rmax.lt.r)then
    er2=dot_product(s(rmax+1:r),s(rmax+1:r))
    r=rmax
   end if
  end if
  if(present(tol))then
   nrm=dnrm2(size(s),s,1)
   bound=tol*tol*nrm*nrm
   er=er2+s(r)*s(r)
   do while(er.lt.bound)
    er2=er; r=r-1; er=er+s(r)*s(r)
   end do
  end if
  if(present(err))err=dsqrt(er2)
  return
 end function


! NORM2
 double precision function norm2_d(a) result(nrm)
  implicit none
  double precision,intent(in)  :: a(:,:)
  character(len=*),parameter :: subnam='norm2_d'
  double precision,pointer :: u(:,:),v(:,:),s(:)
  integer :: info
  call d_svd(a,u,s,v,tol=1.d-3,info=info)
  if(info.ne.0)then;write(*,'(a,a,i10)')subnam,': svd info: ',info;stop;end if
  nrm=s(1)
  deallocate(u,v,s)
 end function

 double precision function norm2p_d(a,x) result(nrm)
  implicit none
  double precision,intent(in)  :: a(:,:)
  double precision,intent(inout) :: x(size(a,2))
  character(len=*),parameter :: subnam='norm2p_d'
  double precision,allocatable :: y(:)
  double precision :: xnrm,ynrm
  integer :: m,n,i,info
  double precision,external :: dnrm2

  m=size(a,1); n=size(a,2)
  allocate(y(m),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif

  xnrm=dnrm2(n,x,1)
  if(xnrm.le.1.d-14)then
   write(*,*)subnam,': please provide non-zero x!'
   nrm=-1
   return
  end if
  if(xnrm.ne.1.d0)call dscal(n,1.d0/xnrm,x,1)

  do i=1,32
   call dgemv('n',m,n,1.d0,a,m,x,1,0.d0,y,1)
   ynrm=dnrm2(m,y,1)
   if(ynrm.ne.1.d0)call dscal(m,1.d0/ynrm,y,1)
   call dgemv('t',m,n,1.d0,a,m,y,1,0.d0,x,1)
   xnrm=dnrm2(n,x,1)
   if(xnrm.ne.1.d0)call dscal(n,1.d0/xnrm,x,1)
   nrm=xnrm
  end do

  deallocate(y)
 end function

 !> SVD translated from from MatLab
 subroutine plain_svd(a,u,s,v,info)
   implicit none
   ! argument parameters
   double precision,intent(in)    :: a(:,:)
   double precision,intent(inout) :: u(:,:),s(:),v(:,:)
   integer,         intent(out)   :: info
   ! local parameters
   character(len=*), parameter    :: subnam='plain_svd'
   double precision, allocatable  :: b(:,:),work(:)
   integer,          allocatable  :: iwork(:)
   integer                        :: m,n,mn,mx,lwork,r,i,j
   double precision, external     :: dnrm2
   !
   m  = size(a,1)
   n  = size(a,2)
   mn = min(m,n)
   mx = max(m,n)
   !
   allocate(work(1), b(m,n))
   !
   lwork=-1
   call dgesvd('s','s',m,n,b,m,s,u,m,v,mn,work,lwork,info)
   !
   lwork=int(work(1))+1
   deallocate(work)
   allocate(work(lwork),iwork(8*mn), stat=info)
   !
   if(info.ne.0) then
     write(*,*) subnam,': cannot allocate'
     stop
   endif
   !
   call dcopy(m*n,a,1,b,1)
   call dgesvd('s','s',m,n,b,m,s,u,m,v,mn,work,lwork,info)
   !
   deallocate(b,work,iwork)
   !
 end subroutine plain_svd
 

 !> Frobenius norm of a matrix
 double precision function norm1fro_d(a) result(nrm)
 double precision,intent(in)  :: a(:)
   nrm= sqrt(sum(a*a))
 end function


 double precision function norm2fro_d(a) result(nrm)
 double precision,intent(in)  :: a(:,:)
   nrm= sqrt(sum(a*a))
 end function


 double precision function norm3fro_d(a) result(nrm)
 double precision,intent(in)  :: a(:,:,:)
   nrm= sqrt(sum(a*a))
 end function


 double precision function norm4fro_d(a) result(nrm)
 double precision,intent(in)  :: a(:,:,:,:)
   nrm= sqrt(sum(a*a))
 end function


 !> Bounded Frobenius norm (divided by the size)
 double precision function norm1fro_bounded_d(a) result(nrm)
 double precision,intent(in)  :: a(:)
   nrm= sqrt(sum(a*a)/dble(size(a)))
 end function


 double precision function norm2fro_bounded_d(a) result(nrm)
 double precision,intent(in)  :: a(:,:)
   nrm= sqrt(sum(a*a)/dble(size(a)))
 end function


 double precision function norm3fro_bounded_d(a) result(nrm)
 double precision,intent(in)  :: a(:,:,:)
   nrm= sqrt(sum(a*a)/dble(size(a)))
 end function




 !> QR decomposition, Matlab-style
 !!
 !! qr(A, Q, R) is equivalent to matlab command `[Q, R] = qr(A, 0)`, which in
 !! turn is the same as `[Q, R, P] = qr(A, "vector", "econ")
 !! - P is the permutation vector;
 !! - the "vector" option means that A(:,P) = Q*R (as opposed to A*P = Q*R
 !!   if P were a matrix);
 !! - the "econ" option produces an economy-size decomposition, which means
 !!   that if A is an (m,n)-size matrix, and m >= n, qr computes only the
 !!   first n columns of Q and the first n rows of R:
 !!         / R \           / R \
 !!   A = Q |   | = (Q1 Q2) |   | = Q1 * R (we also denote Q1 as C)
 !!         \ 0 /           \ 0 /
 !!
 !! Matlab example:
 !! >> A = rand(4, 4)
 !!
 !! A =
 !!
 !!     0.8147    0.6324    0.9575    0.9572
 !!     0.9058    0.0975    0.9649    0.4854
 !!     0.1270    0.2785    0.1576    0.8003
 !!     0.9134    0.5469    0.9706    0.1419
 !!
 !! >> [Q, R, P] = qr(A, "vector", "econ")
 !!
 !! Q =
 !!
 !!    -0.5707    0.4307    0.2564   -0.6504
 !!    -0.5751   -0.0867   -0.8029    0.1308
 !!    -0.0939    0.7694    0.0861    0.6259
 !!    -0.5785   -0.4636    0.5313    0.4100
 !!
 !!
 !! R =
 !!
 !!    -1.6777   -0.9827   -0.7595   -1.5263
 !!          0    0.9201    0.2246   -0.0534
 !!          0         0    0.3983   -0.0222
 !!          0         0         0    0.0425
 !!
 !!
 !! P =
 !!
 !!      3     4     2     1
 !!
 !!>> norm( A(:,P))
 !!
 !!ans =
 !!
 !!   4.8600e-16
 !!
 subroutine matlab_qr_d(Q, R, A)
 implicit none
 double precision, intent(inout), allocatable :: Q(:,:)
 double precision, intent(inout), allocatable :: R(:,:)
 double precision, intent(in) :: A(:,:)
 !
 integer :: m, n, i, mn, lda, lwork, info
 double precision, allocatable :: tau(:), work(:)

    ! in LAPAQK, the routine DGEQRF computes the QR factorization
    ! (https://netlib.org/lapack/lug/node40.html)
    m = size(A,1)
    n = size(A,2)
    mn = min(m,n)
    R = A

    ! The routine  xGEQRF computes  the QR factorization.  The matrix Q is not
    ! formed explicitly,  but is represented  as a product  of elementary  ref-
    ! lectors.  Associated  routines are provided  to work  with Q: xORGQR (or 
    ! xUNGQR in  the complex  case) can generate  all  or  part  of  Q,  while 
    ! xORMQR  (or xUNMQR) can pre- or post-multiply a given matrix by Q or QT 
    ! (QH if complex).

    ! determine the work size
    lwork = -1
    allocate(work(1))
    call dgeqrf(m,n,R,max(1,m),tau,work,lwork,info)
    if(info.ne.0) error stop info

    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork), tau(mn))
    call dgeqrf(m,n,R,max(1,m),tau,work,lwork,info)
    if(info.ne.0) error stop info

    if (allocated(Q)) deallocate(Q)
    allocate(Q(m,mn)); Q = 0d0; Q(1:m,1:mn) = R(1:m,1:mn)
    call dorgqr(m,mn,mn,Q,m,tau,work,lwork,info)
    if(info.ne.0) error stop info
    deallocate(work,tau)
    R = R(1:mn,1:n)

    do i=2,mn
       R(i,1:i-1)= 0d0
    enddo

 end subroutine matlab_qr_d


 !> Cholesky factorization (matlab version)
 !!
 !! U = chol(A) factorizes symmetric positive definite matrix A into
 !! an upper triangular U that satisfies A = U'*U. If A is nonsymmetric,
 !! then chol treats the matrix as symmetric and uses only the diagonal
 !! and upper triangle of A.
 !!
 !! flag_: an optional integer indicating whether A is symmetric positive
 !!        definite. If the flag_ is supplied, the function does not raise
 !!        an error when A is not symmetric positive-definite.
 !!
 !!  * flag-> 0: A is symmetric positive-definite, factorization successful;
 !!  * flag=/=0: A is not SPD, the value of the flag indicates the index of
 !!              the pivot position where the factorization failed.
 !!
 function d_chol(A, flag_) result(U)
 implicit none
 double precision, intent(in) :: A(:,:)
 integer, intent(out), optional :: flag_
 double precision, allocatable :: U(:,:)
 !
 integer :: info, i, N

    U = A
    N = size(A,1)
    call dpotrf('U', N, U, max(1,N), info)
    if (present(flag_)) flag_= info
    do i=2,N
       U(i,1:i-1)= 0d0
    enddo

 end function d_chol

 !> stack two 2D matrices together
 !!
 !! Equivalent to using matlab expression: c = [a,b]
 !! If d_ = 1, stack them vertically (default d_= 2)
 !!
 function stack(A, B, d_) result(C)
 implicit none
 double precision, allocatable:: C(:,:)
 integer, intent(in), optional :: d_
 double precision, intent(in) :: A(:,:), B(:,:)
 !
 integer :: i, d, m1, m2, n1, n2

    d = 2; if (present(d_)) d = d_
    m1 = size(A, 1); n1 = size(A, 2)
    m2 = size(B, 1); n2 = size(B, 2)
    if (d.eq.1) then
       allocate(C(m1+m2, n1))
       C(1:m1,1:n1) = A(1:m1,1:n1)
       C(m1+1:m1+m2,1:n1) = B(1:m2,1:n2) ! n1 == n2
    else if (d.eq.2) then
       allocate(C(m1, n1+n2))
       C(1:m1,1:n1) = A(1:m1,1:n1)
       C(1:m1,n1+1:n1+n2) = B(1:m2,1:n2) ! m1 == m2
    endif

 end function stack


 !> Creates a diagonal matrix from a vector
 !!
 !! Takes a vector v and creates a matrix with d being its diagonal
 !! TODO: generalize to use other matlab interfaces
 function diag_v2m(v) result (A)
 implicit none
 double precision, intent(in)  :: v(:)
 double precision, allocatable :: A(:,:)
 !
 integer :: i, d

    d = size(v)
    allocate(A(d,d))
    A = 0d0
    do i=1,d
       A(i,i) = v(i)
    enddo

 end function diag_v2m

 !> Takes a matrix and extracts the diagonal
 !!
 !! like in matlab
 function diag_m2v(A) result (v)
 implicit none
 double precision, intent(in)  :: A(:,:)
 double precision, allocatable :: v(:)
 !
 integer :: i, d

    d = min(size(A,1),size(A,2))
    allocate(v(d))
    do i=1,d
       v(i) = A(i,i)
    enddo

 end function diag_m2v


 !> Highly experimental acceleration of SVD/QR using Gram matrix.
 !  Use with caution for m>>n only!
 !    call svdgram(u,s,r,A,[tol])
 !    u is the left singular factor of A,
 !    s is the singular values (vector!),
 !    r has the meaning of diag(s)*v'.
 !    if tol is given, performs the truncation with Fro-threshold.
 subroutine svdgram(up,sp,rp,A,tol)
 implicit none
 double precision, pointer, intent(out):: up(:,:), sp(:), rp(:,:)
 double precision, intent(in):: A(:,:)
 double precision, intent(in), optional :: tol
 !
 double precision :: tol_
 double precision, pointer :: vp(:,:)
 double precision, allocatable :: R2(:,:), iR(:,:), RR(:,:)
 double precision, allocatable :: u(:,:), s(:), v(:,:), r(:,:)
 integer :: p, sh(2)

 ! WIP
   ! TODO
   R2 = matmul(transpose(A), A)   ! R2 = A'*A
   call d_svd(R2,up,sp,vp)        ! [u,s,v]=svd(R2, 'econ');
   deallocate(R2)
   v = transpose(vp)              ! v needs to be transposed for compat. w/matlab
   u = matmul(A, v)               ! u = A*v
   s = sum(u**2, 1)               ! s = sum(u.^2, 1);
   s = sqrt(s)                    ! s = sqrt(s.') "s.'" is a non-conjugate transpose
                                  !               (simple transpose)

   if (present(tol)) then         ! if (nargin>1)&&(~isempty(tol))
      p = my_chop2(s, sqrt(sum(s**2))*tol)!     p = my_chop2(s, norm(s)*tol);
      u = u(:,1:p) !!! but u is tr!     u = u(:,1:p);
      s = s(1:p)                  !     s = s(1:p);
      v = v(:,1:p)                !     v = v(:,1:p);
   else                           ! end;
      p = size(s)
   endif

   u = matmul(u, matlab_diag(1d0/s))  !  u = u*spdiags(1./s, 0, numel(s), numel(s));
   !                                     if (issparse(u))
   !                                         u = full(u);
   !                                     end;
                                      !  r = diag(s)*v';
   r = matmul(matlab_diag(s), vp)

   !% Run chol for reortogonalization.
   !% It will stop if the matrix will be singular.
   !% Fortunately, it means rank truncation with eps in our business.
   if (s(1)/s(p).gt.1d7) then         !  if (s(1)/s(end)>1e7)
      p = 1                           !      p = 1;
      do while(p.gt.0d0)              !      while (p>0)
         R2 = matmul(transpose(u), u) !          R2 = u'*u;
         RR = d_chol(R2, flag_=p)     !          [R,p] = chol(R2);
         deallocate(R2)
         if (p.gt.0) then             !          if (p>0)
            u = u(:, 1:p-1)           !              u = u(:,1:p-1);
            s = s(1:p-1)              !              s = s(1:p-1);
            r = r(1:p-1,:)            !              r = r(1:p-1,:);
         endif                        !          end;
         call matinv(RR, iR)          !          iR = inv(R);
         u = matmul(u, iR)            !          u = u*iR;
         r = matmul(RR, r)            !          r = R*r;
         deallocate(RR,iR)
      enddo                           !      end;
   endif                              !  end;

   sh = shape(u)
   deallocate(up)
   allocate(up(sh(1),sh(2)))
   up = u

   deallocate(sp)
   allocate(sp(size(s)))
   sp = s

   sh = shape(r)
   allocate(rp(sh(1),sh(2)))
   rp = r

 end subroutine svdgram

 !> Truncation by absolution precision in Frobenius norm
 !!   [R]=MY_CHOP2(SV,EPS) truncation by absolute precision in the Frobenius
 !!   norm. Finds minimal possible R such that \sqrt(sum(SV(R:n)^2)) < EPS
 !!
 !!  Adopted from TT-Toolbox 2.2, 2009-2012
 !!
 !! ---------------------------
 integer function my_chop2(sv,eps) result(r)
 double precision, intent(in):: sv(:)
 double precision, intent(in):: eps
 !
 double precision, allocatable:: sv0(:)
 double precision :: eps2
 integer :: i, sz

    sz = size(sv)
    eps2 = eps*eps
    if (sum(sv**2).eq.0d0) then
       r = 1
    else if (eps.le.0d0) then
       ! Check for zero tolerance
       r = sz
    else
       sv0 = cumsum(sv(sz:1:-1)**2)
       r = sz                         !ff = find(sv0 < eps.^2)
       do i=1,sz                      !if (isempty(ff) )
          if (sv0(i).ge.eps2) then    !   r=numel(sv);
             r = sz - i + 1; exit     !else
          endif                       !   r=numel(sv)-ff(end);
       enddo                          !end
       deallocate(sv0)
    endif
    ! GPU-firendly code
    !do i=1,sz
    !   if (sv0(i).lt.eps2) exit
    !enddo
    !r = sz - i*int(i.le.sz) ! GPU-friendly
 end function

 !> Function to compute cumulative sum
 !! TODO: optimize (can be O(log(N))
 function i_cumsum(arr) result(carr)
 integer, dimension(:), intent(in) :: arr
 integer, dimension(size(arr)) :: carr
 integer :: i

   carr(1) = arr(1)
   do i = 2, size(arr)
      carr(i) = carr(i-1) + arr(i)
   end do
 end function i_cumsum

 !> Function to compute cumulative sum
 !! TODO: optimize (can be O(log(N))
 function d_cumsum(arr) result(carr)
 double precision, dimension(:), intent(in) :: arr
 double precision, dimension(size(arr)) :: carr
 integer :: i

   carr(1) = arr(1)
   do i = 2, size(arr)
      carr(i) = carr(i-1) + arr(i)
   end do
 end function d_cumsum

 !> adopted from matlab
 function bfun3(Phi1, A, Phi2, x) result (y)
 double precision, intent(in):: Phi1(:,:,:)
 double precision, intent(in):: A(:,:)
 double precision, intent(in):: Phi2(:,:,:)
 double precision, intent(in):: x(:,:)
 double precision, allocatable :: y(:,:)
 !
 integer :: b, rx1,ry1,ra1,rx2,ry2,ra2,n,m
 double precision, allocatable :: P1(:,:), P2(:,:)

    b= size(x, 2)
    ! Phi1: ry1, rx1, ra1
    ry1 = size(Phi1,1)
    rx1 = size(Phi1,2)
    ra1 = size(Phi1,3)

    ! % Phi2: rx2, ra2, ry2
    ry2 = size(Phi2,3)
    rx2 = size(Phi2,1)
    ra2 = size(Phi2,2)
    n = size(A,1)/ra1;
    m = size(A,2)/ra2;

    y = reshape(x, (/ b*rx1*m, rx2/))
    P2 = reshape(Phi2, (/rx2, ra2*ry2/))
    y = thor_matmul(y, P2)
    y = transpose(reshape(y, (/ b*rx1, m*ra2*ry2 /)))
    y = reshape(y, (/ m*ra2, ry2*b*rx1 /))
    y = thor_matmul(A, y)
    y = transpose(reshape(y, (/ ra1*n*ry2*b, rx1 /)))
    y = reshape(y, (/rx1*ra1, n*ry2*b/))
    P1 = reshape(Phi1, (/ ry1, rx1*ra1 /))
    y = thor_matmul(P1, y)
    y = reshape(y, (/ ry1*n*ry2, b/))

    deallocate(P1, P2)

 end function bfun3


 function dzeros_1d(m) result(a)
 integer, intent(in) :: m
 double precision, allocatable :: a(:)
 allocate(a(m)); a = 0d0
 end function dzeros_1d
 function dzeros_2d(m,n) result(a)
 integer, intent(in) :: m, n
 double precision, allocatable :: a(:,:)
 allocate(a(m,n)); a = 0d0
 end function dzeros_2d
 function dzeros_3d(m,n,p) result(a)
 integer, intent(in) :: m,n,p
 double precision, allocatable :: a(:,:,:)
 allocate(a(m,n,p)); a = 0d0
 end function dzeros_3d
 function dzeros_4d(m,n,p,q) result(a)
 integer, intent(in) :: m,n,p,q
 double precision, allocatable :: a(:,:,:,:)
 allocate(a(m,n,p,q)); a = 0d0
 end function dzeros_4d

 subroutine d_svd_old(a,u,v,s,tol,rmax,err,info)
  implicit none
  real(8),intent(in) :: a(:,:)
  real(8),pointer :: u(:,:),v(:,:),s(:)
  real(8),intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  real(8),intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='d_svd'
  real(8),allocatable :: b(:,:),work(:),ss(:),uu(:,:),vv(:,:)
  integer :: m,n,mn,lwork,ierr,r,i,j
  real(8),external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); lwork=256*max(m,n)
  allocate(uu(m,mn),vv(mn,n),ss(mn),b(m,n),work(lwork),stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call dcopy(m*n,a,1,b,1)
  call dgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,ierr)
  if(present(info))info=ierr
  if(ierr.ne.0)then
   write(*,*)subnam,': dgesvd info: ',ierr
   if(ierr.lt.0)stop
   if(nan(a))then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   else
    write(*,*) subnam,': min/max element of input: ',minval(a),maxval(a)
   end if
   u=>null(); v=>null(); s=>null()
  else
   r=chop(ss,tol,rmax,err)
   if(present(err))err=err/dnrm2(mn,ss,1)
   allocate(u(m,r),v(n,r),s(r))
   call dcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   forall(i=1:n,j=1:r)v(i,j)=vv(j,i)
  end if
  deallocate(uu,vv,ss,b,work)
 end subroutine
 subroutine z_svd_old(a,u,v,s,tol,rmax,err,info)
  implicit none
  double complex,intent(in) :: a(:,:)
  double complex,pointer :: u(:,:),v(:,:)
  real(8),pointer :: s(:)
  real(8),intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  real(8),intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='z_svd'
  double complex,allocatable :: b(:,:),work(:),uu(:,:),vv(:,:)
  real(8),allocatable :: ss(:),rwork(:)
  integer :: m,n,mn,lwork,ierr,r,i,j
  real(8),external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); lwork=256*max(m,n)
  allocate(uu(m,mn),vv(mn,n),ss(mn),b(m,n),work(lwork),rwork(lwork),stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call zcopy(m*n,a,1,b,1)
  call zgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,ierr)
  if(present(info))info=ierr
  if(ierr.ne.0)then
   write(*,*)subnam,': zgesvd info: ',ierr
   if(ierr.lt.0)stop
   if (nan(a)) then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   end if
   u=>null(); v=>null(); s=>null()
  else
   r=chop(ss,tol,rmax,err)
   if(present(err))err=err/dnrm2(mn,ss,1)
   allocate(u(m,r),v(n,r),s(r))
   call zcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   forall(i=1:n,j=1:r)v(i,j)=vv(j,i)
  end if
  deallocate(uu,vv,ss,b,work)
 end subroutine

 !> Converts a sequence of piv indices to sequence for indexing arrays
 !!
 !! Converts a sequence of piv indices to a sequence of indices useful for indexing
 !! arrays. Piv indicies come out of lapack with lu_factor, and is a sequence of
 !! indices where the i-th row should be swapped with piv[i].
 !! Swapping should be done in-place.
 !! input:
 !!   piv: 1d array, a sequence of piv indices
 !! output:
 !!   perm: 1d array, a sequence of indices useful for indexing arrays
 !!
 !! Adopted from: https://gitlab.lanl.gov/ewskau/python_amen_cross, 'piv_to_py()'
 function piv_to_perm(piv) result(perm) ! a.k.a piv_to_py
 integer, intent(in) :: piv(:)
 integer :: perm(size(piv))
 !
 integer :: i, j, k
 do i=1,size(piv)
    perm(i)= i
 enddo
 do i=1,size(piv)
    ! in perm, swap i and piv(j)
    j = piv(i)
    k = perm(i); perm(i) = perm(j); perm(j) = k
 enddo
 end function piv_to_perm


 !> Matrix-matrix multiplication with BLAS dgemm function
 function thor_matmul(A, B, tsp_) result(C)
 implicit none
 double precision, allocatable :: C(:,:)
 double precision, intent(IN):: A(:,:), B(:,:)
 character(len=2), optional :: tsp_  !< are matrices transposed (NN/NT/TN/TT)
 !
 character(len=2) :: tsp
 character(len=1) :: tspA, tspB
 logical :: A_transposed, B_transposed
 integer :: m, n, k, k1, ldA, ldB, ldC, shA(2), shB(2)

    tsp= 'nn'; if (present(tsp_)) tsp= tsp_
    tspA= tsp(1:1); tspB= tsp(2:2)
    A_transposed = (tspA.eq.'t')
    B_transposed = (tspB.eq.'t')

    shA = shape(A)
    shB = shape(B)
    
    if (A_transposed) then
       k = shA(1); m = shA(2)
    else
       m = shA(1); k = shA(2)
    endif

    if (B_transposed) then
       n = shB(1); k1 = shB(2)
    else
       n = shB(2); k1 = shB(1)
    endif
    allocate(C(m,n))

    call dgemm(tsp(1:1),  & ! transa
               tsp(2:2),  & ! transb
               m, n, k,   & ! A(m,k)*B(k,n) -> C(m,n)
               1d0,       & ! alpha
               A, shA(1), & ! A(lda,*)
               B, shB(1), & ! B(ldb,*)
               0d0,       & ! beta
               C, m)        ! C(ldc,*)
 end function thor_matmul


 !> Matrix-matrix multiplication with BLAS dgemm function
 function thor_matmul_mnk(m, n, k, A, B, tsp_) result(C)
 implicit none
 double precision, allocatable :: C(:,:)
 double precision, intent(IN):: A(*), B(*)
 integer, intent(IN) :: m, n, k
 character(len=2), optional :: tsp_  !< are matrices transposed (NN/NT/TN/TT)
 !
 character(len=2) :: tsp
 character(len=1) :: tspA, tspB
 logical :: A_transposed, B_transposed
 integer :: k1, ldA, ldB, ldC, shA(2), shB(2)

    tsp= 'nn'; if (present(tsp_)) tsp= tsp_
    tspA= tsp(1:1); tspB= tsp(2:2)
    A_transposed = (tspA.eq.'t')
    B_transposed = (tspB.eq.'t')

    if (A_transposed) then
       shA = [k, m]
    else
       shA = [m, k]
    endif

    if (B_transposed) then
       shB = [n, k]
    else
       shB = [k, n]
    endif
    allocate(C(m,n))

    call dgemm(tsp(1:1),  & ! transa
               tsp(2:2),  & ! transb
               m, n, k,   & ! A(m,k)*B(k,n) -> C(m,n)
               1d0,       & ! alpha
               A, shA(1), & ! A(lda,*)
               B, shB(1), & ! B(ldb,*)
               0d0,       & ! beta
               C, m)        ! C(ldc,*)
 end function thor_matmul_mnk


 !> Returns block of n row indices with maxvol square submatrix
 !!
 !! Takes a $m \times n$ matrix with $n \leq m$, and returns $n$ row indices that
 !! provide the maximum volume square submatrix.
 !!
 !!     Parameters:
 !!         X (array): An $m \times n$ array.
 !!         inds (list of ints): The starting collection of indices, default = None.
 !!         tol (double): a tolerance paramter, default = 5e-2.
 !!
 !!     Returns:
 !!         inds (array): A list of $n$ row indices that provide the maximum volume submatrix.
 !!
 !! Adopted from Erik Skau's code
 !! https://gitlab.lanl.gov/ewskau/python_amen_cross/-/blob/main/amen_cross/maxvol.py?ref_type=heads
 !!
 function maxvol(X, inds_, tol_) result(inds)
 use linalg_module, only: assert
 implicit none
 double precision, intent(inout) :: X(:,:)
 integer, intent(in), optional :: inds_(:)
 double precision, intent(in), optional :: tol_
 integer, allocatable :: inds(:)
 !
 integer :: i, sh(2)
 double precision :: tol, rcond
 double precision, allocatable :: Y(:,:), A(:,:), B(:,:), wrk(:)
 double precision, allocatable :: A1(:,:), A2(:,:), A12(:,:)
 integer, allocatable :: ipiv(:), jpvt(:), mxl(:)
 integer :: info, N, mn, Nrhs, lwrk, mi, mi1, mi2, rankA

    tol=5d-2; if (present(tol_)) tol= tol_
    allocate(inds(0)); if (present(inds_)) inds= inds_
    Y = X ! create copy of X so it doesn't get destroyed in the function

    sh = shape(X)
    if (sh(1).eq.sh(2)) then
       allocate(inds(sh(1)))
       do i=1,sh(1)
          inds(i) = i
       enddo
       return
    endif

    call assert (sh(1).ge.sh(2), &
                 "Maxvol works for m by n matrices where m >= n.")

    if (size(inds).ne.sh(2))then
       allocate(ipiv(sh(2)))

       ! lu,piv = scipy.linalg.lu_factor(X)
       call dgetrf(sh(1), sh(2), Y, sh(1), ipiv, info)

       ! inds = piv_to_py(piv)
       inds = piv_to_perm(ipiv)
    endif

    !#A = np.linalg.solve(X[inds, :].T, X.T).T
    !A = np.linalg.lstsq(X[inds,:].T, X.T, rcond=None)[0].T
    N = size(inds); mn = min(N,sh(2)); Nrhs = sh(1)
    lwrk = max(mn + 3+1, 2*mn + Nrhs)
    allocate(A(N,sh(2)),jpvt(N),wrk(lwrk))
    B = transpose(X)
    rcond = 1d0
    call dgelsd(sh(1),N,sh(1),A,sh(1),B,sh(2),jpvt, rcond,rankA,wrk,lwrk,info)
    A = B

    maxvol_loop: do i=1,1000
       !mi = np.unravel_index(np.argmax(np.abs(A)), X.shape)
       mxl = maxloc(reshape(dabs(A), (/size(A)/) ))
       mi = mxl(1)
       mi1 = mi / sh(1) + 1
       mi2 = mod(mi, sh(1)) + 1
       if (A(mi1,mi2).le.1d0 + tol) exit maxvol_loop

       ! construct new A:
       ! A = A + A[:,[mi[1]]]@(A[[inds[mi[1]]], :] - A[[mi[0]],:])/A[mi[0],mi[1]]
       ! construct A[:,[mi[1]]]
       if (allocated(A1)) deallocate(A1)
       allocate(A1(size(A,1),1))
       A1(:,1) = A(:,mi2)

       ! construct A[[inds[mi[1]]],:] - A[[mi[0]],:]
       if (allocated(A2)) deallocate(A2)
       allocate(A2(1,size(A,2)))
       A2(1,:) = A(inds(mi2),:) - A(mi1, :)

       ! tensor product A1'@A2
       A12 = matmul(A1,A2)/A(mi1,mi2)

       ! finally
       A = A - A12

       !# A = np.linalg.lstsq(X[inds,:].T, X.T, rcond=None)[0].T
       inds(mi2) = mi1
    enddo maxvol_loop

    !return np.sort(inds)
    !! TODO

    if (allocated(A1))   deallocate(A1)
    if (allocated(A2))   deallocate(A2)
    if (allocated(A12))  deallocate(A12)
    if (allocated(A))    deallocate(A)
    if (allocated(B))    deallocate(B)
    if (allocated(wrk))  deallocate(wrk)
    if (allocated(ipiv)) deallocate(ipiv)
    if (allocated(jpvt)) deallocate(jpvt)
    if (allocated(mxl))  deallocate(mxl)
    deallocate(Y)

end function maxvol



end module
