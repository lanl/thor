module dmrgg_lib
 use thor_lib
 use lr_lib
 use rnd_lib
 use default_lib
 use time_lib
 use thor_pointers, only : pointi2, pointd3, pointz3
 implicit none

contains

 subroutine dtt_dmrgg(arg,fun,par,accuracy,maxrank,pivoting,neval,quad,tru)
 ! approximate [fun] in TT format using dmrg greedy cross interpolation
 implicit none
 type(dtt_tensor),intent(inout),target :: arg
 double precision,external :: fun
 double precision,intent(in),optional :: par(*)
 double precision,intent(in),optional :: accuracy
 integer,intent(in),optional :: maxrank
 integer,intent(in),optional :: pivoting
 integer(kind=8),intent(out),optional :: neval
 type(dtt_tensor),intent(in),optional :: quad
 double precision,intent(in),optional :: tru

 character(len=*),parameter :: subnam='[dtt_dmrgg]'
 integer,parameter :: lft=1,rgt=2,smin=8
 double precision :: small_element, small_pivot
 integer,pointer :: r(:),n(:)
 integer :: l,m,i,j,k,p,q,ii,jj,kk,pp,qq,nn,s,t,it,crs,dir,ijkq,ij,jk,kq,nlot,ilot,piv,snum
 integer :: ihave,ineed,strike, first,last
 integer :: nproc, info, me, they
 integer(kind=8) :: nevalloc,nevalall
 double precision :: t1,t2,t3,amax,pivot,pivotmax,pivotmin,val,err, pivotmax_prev,val_prev
 double precision,allocatable :: a(:,:,:,:),b(:),c(:),d(:),bcol(:,:,:),brow(:,:,:),bcol1(:,:),brow1(:,:),acol1(:,:),arow1(:,:)
 double precision :: bb(4),cc(4)
 integer,allocatable :: lot(:,:),shifts(:),own(:)
 type(dtt_tensor) :: col,row,ttqq
 type(pointd) :: inv(0:tt_size)
 type(pointi2) :: vip(0:tt_size)
 integer :: ind(tt_size),rr(0:tt_size),tape(4,0:tt_size),tmpp(4,0:tt_size)
 logical :: ready,done,haverow,havecol,skipcol,needsend,needrecv, upd(0:tt_size)
 character(len=2) :: sdir
 character(len=1024) :: str,stmp,sfmt
 integer,external :: idamax
 double precision,external :: dnrm2,ddot

  t1=timef()
  piv=default(3,pivoting)
  ! check how our double precision type compares to hard-coded MPI types
  select case(storage_size(1.d0))
   case(32);  small_element= 5*epsilon(1.d0);  small_pivot=1.d-3; sfmt='e14.8'
   case(64);  small_element=10*epsilon(1.d0);  small_pivot=1.d-5; sfmt='e20.14'
   case(128); small_element=50*epsilon(1.d0);  small_pivot=1.d-7; sfmt='e39.33'
   case default
    write(*,'(a,a,i5,a)')subnam,': unknown storage_size(1.d0): ',storage_size(1.d0),' guessing 64'
    small_element=5*epsilon(1.d0);  small_pivot=1.d-5; sfmt='e20.14'
  end select

  ! initialise dimensions, mode sizes and bond ranks
  if(arg%l.gt.arg%m)then;write(*,*)subnam,': l,m: ',arg%l,arg%m;stop;endif
  l=arg%l; m=arg%m; r=>arg%r; n=>arg%n
  !if(any(r(l-1:m).gt.1))then;if(me.eq.0)write(*,*)'arg came in non-empty, wiping it clear';endif
  r(l-1:m)=1; call alloc(arg)
  nproc = 1; me = 0

  ! distribute bonds (ranks) between procs
  allocate(own(0:nproc), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate own bonds';stop;endif
  own(0) = l; own(1) = m

  ! allocating vip: i=vip(p)%p(1,r); j=vip(p)%p(2,r); k=vip(p)%p(3,r); q=vip(p)%p(4,r)
  allocate(shifts(0:nproc), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate init shifts';stop;endif
  do p=l-1,m
   allocate(inv(p)%p(1),vip(p)%p(4,1), stat=info)
   if(info.ne.0)then;write(*,*)subnam,': cannot allocate init inv and vip';stop;endif
   inv(p)%p(1)=1.d0
  end do

  !if(me.eq.0)write(*,*) 'locating initial cross'
  snum=max(smin,nproc) ! number of shifts
  do p=0,nproc-1
   shifts(p) = int(dble(snum)*dble(p)/nproc)
  end do
  shifts(nproc) = snum
  !write(*,'(a,i2,a,32i4)')'[',me,']: shifts: ',shifts(0:nproc)
  nn=minval(n(l:m))
  nlot=nn*(shifts(me+1)-shifts(me))
  allocate(b(nlot), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate b';stop;endif
  ihave=shifts(me+1)-shifts(me)
  do s=0,ihave-1
!$OMP PARALLEL DO PRIVATE(p,ind)
   do k=1,nn
    forall(p=l:m)ind(p)=mod(k-1+(s+shifts(me))*(p-1),n(p))+1
    b(k+s*nn)=fun(m,ind,n,par)
   end do
!$OMP END PARALLEL DO
  end do
  ilot=idamax(nlot,b,1); amax=dabs(b(ilot))
  nevalloc=nlot
  ilot=ilot+nn*shifts(me)
  deallocate(b)

  s=(ilot-1)/nn; k=mod(ilot-1,nn)+1; forall(p=l:m)ind(p)=mod(k-1+s*(p-1),n(p))+1
  s=(ilot-1)/nn; k=mod(ilot-1,nn)+1; forall(p=l:m)ind(p)=mod(k-1+s*(p-1),n(p))+1
  !write(*,'(a,i2,a,e10.3,a,512i4)')'[',me,'] global max  ',amax,' at ',ind(l:m)
  !write(*,'(a,2(1x,i4),a,e20.13)')'PIVOT ',ind(l:m),' : ',amax
  vip(l-1)%p(:,1)=(/1,1,1,1/)
  forall(p=l:m-1)vip(p)%p(:,1)=(/ 1, ind(p),ind(p+1), 1  /)
  vip(m)%p(:,1)=(/1,1,1,1/)

  !if(me.eq.0)write(*,*)' computing initial cross'
  do p=own(me),own(me+1)
!$OMP PARALLEL DO
   do j=1,n(p)
     arg%u(p)%p(1,j,1)=dmrgg_fun(1,j,ind(p+1),1, p, l,m,vip,r,n, fun, par)
   enddo
!$OMP END PARALLEL DO
   nevalloc=nevalloc+n(p)
   do j=1,n(p); amax=dmax1(amax,dabs(arg%u(p)%p(1,j,1))); enddo
   !write(*,'(a,i2,a,i2,a,50e10.3)')'[',me,']{',p,'}: fiber: ', arg%u(p)%p(1,:,1)
  end do
  pivotmax_prev = amax
  do p=own(me),own(me+1)-1
   pivot=arg%u(p)%p(1,ind(p),1)
   inv(p)%p(1)=pivot
   !write(*,'(a,i2,a,i2,a,e20.15)')'[',me,']{',p,'}: cross entry: ', arg%u(p)%p(1,ind(p),1)
  end do

  !if(me.eq.0)write(*,*)'initialise cols and rows of approximation'
  col=arg; row=arg
  do p=own(me),own(me+1)-1
   call d2_lual(n(p),r(p),inv(p)%p,col%u(p)%p)
   call d2_luar(n(p),r(p),inv(p)%p,row%u(p+1)%p)
  end do

  if(present(quad))then
   !if(me.eq.0)write(*,*)'computing quadrature'
   val=1.d0
   do p=own(me),own(me+1)-1
    val=val*ddot(n(p),arg%u(p)%p,1,quad%u(p)%p,1) / inv(p)%p(1)
   end do
   val=val*ddot(n(m),arg%u(m)%p,1,quad%u(m)%p,1)
   val_prev=val
   !if(me.eq.0)write(*,'(a,e20.15)')'initial value:',val
  end if

  nevalall = nevalloc

  !if(me.eq.0)write(*,*)'clearing the tape'
  forall(p=l-1:m)tape(:,p)=(/ -1, -1, -1, -1 /)
  forall(p=l-1:m)tmpp(:,p)=(/ -2, -2, -2, -2 /)
  forall(p=l-1:m)upd(p)=.false.

  it=0; strike=0; ready=.false.; if(present(maxrank))ready=(it+1.ge.maxrank)
  do while(.not.ready) ! main loop increasing the ranks and hopefully accuracy
   it=it+1
   dir=2-mod(it,2); if(dir.eq.1)then;sdir='>>';else if(dir.eq.2)then;sdir='<<';else;sdir='??';endif
   rr(l-1:m)=r(l-1:m);  pivotmax=-1.d0;pivotmin=-1.d0

   do pp=1,own(me+1)-own(me)  ! sweep over processor own's bonds
    p=own(me)+pp-1; if(dir.eq.2)p=own(me+1)-pp
    !write(*,'(a,i2,a,i3,a2,a,i2)') '[,me,']: ' it,sdir,' bond ',p

    allocate(acol1(r(p-1),n(p)),arow1(n(p+1),r(p+1)),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot allocate acol/arow';stop;endif

    if(piv.eq.-1)then ! full pivoting
     allocate(a(r(p-1),n(p),n(p+1),r(p+1)),b(r(p-1)*n(p)*n(p+1)*r(p+1)), bcol(r(p-1),n(p),r(p)),brow(r(p),n(p+1),r(p+1)), &
             stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate superblock';stop;endif
      ! compute long indices and form rhs
!$OMP PARALLEL DO private(i,j,k,q)
     do ijkq=1,r(p-1)*n(p)*n(p+1)*r(p+1)
      i=ijkq-1
      q=i/(r(p-1)*n(p)*n(p+1)); i=mod(i,r(p-1)*n(p)*n(p+1))
      k=i/(r(p-1)*n(p));        i=mod(i,r(p-1)*n(p))
      j=i/r(p-1);               i=mod(i,r(p-1))
      i=i+1;j=j+1;k=k+1;q=q+1
      a(i,j,k,q)=dmrgg_fun(i,j,k,q, p, l,m,vip,r,n, fun, par)
     end do
!$OMP END PARALLEL DO
     nevalloc=nevalloc+r(p-1)*n(p)*n(p+1)*r(p+1)
     ijkq=idamax(r(p-1)*n(p)*n(p+1)*r(p+1),a,1)-1
     q=ijkq/(r(p-1)*n(p)*n(p+1))+1; ijkq=mod(ijkq,r(p-1)*n(p)*n(p+1))
     k=ijkq/(r(p-1)*n(p))+1;        ijkq=mod(ijkq,r(p-1)*n(p))
     j=ijkq/r(p-1)+1;               i=mod(ijkq,r(p-1))+1
     amax=dmax1(amax,dabs(a(i,j,k,q)))

     ! compute residual and do full pivoting
     call dcopy(r(p-1)*n(p)*n(p+1)*r(p+1),a,1,b,1)
     call dgemm('n','n',r(p-1)*n(p),n(p+1)*r(p+1),r(p),-1.d0,col%u(p)%p,r(p-1)*n(p),row%u(p+1)%p,r(p),1.d0,b,r(p-1)*n(p))

     ijkq=idamax(r(p-1)*n(p)*n(p+1)*r(p+1),b,1)-1
     qq=ijkq/(r(p-1)*n(p)*n(p+1))+1; ijkq=mod(ijkq,r(p-1)*n(p)*n(p+1))
     kk=ijkq/(r(p-1)*n(p))+1;        ijkq=mod(ijkq,r(p-1)*n(p))
     jj=ijkq/r(p-1)+1;               ii=mod(ijkq,r(p-1))+1
     pivot=b(ii + r(p-1)*( jj-1 + n(p)*( kk-1 + n(p+1)*( qq-1 ) ) ) )
     !write(*,'(a,i2,a,i2,a,4i5,a,e10.3)')'[',me,']{',p,'} pivot at ',ii,jj,kk,qq,' : ',pivot

     ! copy column and row to containers
     forall(i=1:r(p-1),j=1:n(p))  acol1(i,j)=a(i,j,kk,qq)
     forall(k=1:n(p+1),q=1:r(p+1))arow1(k,q)=a(ii,jj,k,q)
     !do j=1,n(p);do i=1,r(p-1);   acol1(i,j)=a(i,j,kk,qq); enddo;enddo
     !do q=1,r(p+1);do k=1,n(p+1); arow1(k,q)=a(ii,jj,k,q); enddo;enddo

     ! deallocate superblocks
     deallocate(a,b,bcol,brow)

    else if(piv.ge.0)then ! partial pivoting on (r*n+n*r) random elements
     !nlot=r(p-1)*n(p)+n(p+1)*r(p+1)
     nlot=r(p-1)+n(p)+n(p+1)+r(p+1)
     !nlot=1

     !write(*,'(a,i2,a,i2,a,4i4)')'[',me,']{',p,'}: superblock ',r(p-1),n(p),n(p+1),r(p+1)
     !write(*,'(a,i2,a,i2,a,i4,a,i7)')'[',me,']{',p,'}: pivoting ',piv,' nlot ',nlot
     allocate(bcol1(r(p-1),n(p)),brow1(n(p+1),r(p+1)),lot(nlot,4),b(nlot), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate bcol1/brow1 for pivoting';stop;endif

     ! distribute lottery weights for indices columns and rows
     forall(i=1:r(p-1),j=1:n(p))bcol1(i,j)=1.d0
     forall(k=1:n(p+1),q=1:r(p+1))brow1(k,q)=1.d0
     do s=1,r(p)
      i=vip(p)%p(1,s);j=vip(p)%p(2,s);k=vip(p)%p(3,s);q=vip(p)%p(4,s)
      bcol1(i,j)=0.d0; brow1(k,q)=0.d0
     end do
     ! get indices from lottery
     call lottery2(nlot,r(p-1)*n(p),n(p+1)*r(p+1),bcol1,brow1,lot)
     forall(i=1:nlot)lot(i,3)=lot(i,2)
     do ilot=1,nlot
      lot(ilot,2)=(lot(ilot,1)-1)/r(p-1)+1; lot(ilot,1)=mod(lot(ilot,1)-1,r(p-1))+1  ! ij
      lot(ilot,4)=(lot(ilot,3)-1)/n(p+1)+1; lot(ilot,3)=mod(lot(ilot,3)-1,n(p+1))+1  ! kq
     end do

     ! get matrix elements and residuals for these indices
!$OMP PARALLEL DO private(i,j,k,q)
     do ilot=1,nlot
      i=lot(ilot,1); j=lot(ilot,2); k=lot(ilot,3); q=lot(ilot,4)
      b(ilot)=dmrgg_fun(i,j,k,q, p, l,m,vip,r,n, fun, par)
      !write(*,'(a,i2,a,i2,a,4i4,a,4i4,a,e10.3)')'[',me,']{',p,'} elt: ',i,j,k,q,' : ',i,j,k,q,' : ',b(ilot)
     end do
!$OMP END PARALLEL DO
     nevalloc=nevalloc+nlot
     ilot=idamax(nlot,b,1); amax=dmax1(amax,dabs(b(ilot)))
     do ilot=1,nlot
      i=lot(ilot,1); j=lot(ilot,2); k=lot(ilot,3); q=lot(ilot,4)
      !b(ilot)=b(ilot)-ddot(r(p),col%u(p)%p(i,j,1),r(p-1)*n(p),row%u(p+1)%p(1,k,q),1)
      b(ilot)=b(ilot)-sum(col%u(p)%p(i,j,1:r(p))*row%u(p+1)%p(1:r(p),k,q))
      !write(*,'(a,i2,a,i2,a,4i4,a,4i4,a,e10.3)')'[',me,']{',p,'} err: ',i,j,k,q,' : ',i,j,k,q,' : ',b(ilot)
     end do

     ! find maximum pivot
     ilot=idamax(nlot,b,1)
     ii=lot(ilot,1); jj=lot(ilot,2); kk=lot(ilot,3); qq=lot(ilot,4)
     pivot=b(ilot)
     !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)') '[',me,']{',p,'} rnd pivot ',ii,jj,kk,qq,' : ',pivot

     done=.false.; havecol=.false.; haverow=.false.;
     if(piv.eq.0)then ! compute row and column
!$OMP PARALLEL private(i,j,k,q)
!$OMP DO
       do ij=0,r(p-1)*n(p)-1
        j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1
        acol1(i,j)=dmrgg_fun(i,j,kk,qq, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END DO
!$OMP DO
       do kq=0,n(p+1)*r(p+1)-1
        q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1
        arow1(k,q)=dmrgg_fun(ii,jj,k,q, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END DO
!$OMP END PARALLEL
      nevalloc=nevalloc+r(p-1)*n(p)+n(p+1)*r(p+1)
      done=.true.; havecol=.true.; haverow=.true.;
     end if

     ! rook pivoting to increase abs(pivot)
     crs=0; skipcol=(dir.eq.2)
     do while(.not.done)
      if(.not.skipcol)then
!$OMP PARALLEL DO private(i,j)
       do ij=0,r(p-1)*n(p)-1
        j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1
        acol1(i,j)=dmrgg_fun(i,j,kk,qq, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END PARALLEL DO
       nevalloc=nevalloc+r(p-1)*n(p)
       ij=idamax(r(p-1)*n(p),acol1,1)-1; j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1; amax=dmax1(amax,dabs(acol1(i,j)))
       havecol=.true.; crs=crs+1; done=havecol.and.haverow.and.(crs.ge.2*piv)
       if(.not.done)then
        call dcopy(r(p-1)*n(p),acol1,1,bcol1,1)
        call dgemv('n',r(p-1)*n(p),r(p),-1.d0,col%u(p)%p,r(p-1)*n(p),row%u(p+1)%p(1:r(p),kk,qq),1,1.d0,bcol1,1)
        ij=idamax(r(p-1)*n(p),bcol1,1)-1; j=ij/r(p-1)+1; i=mod(ij,r(p-1))+1
        done=havecol.and.haverow.and. (i.eq.ii .and. j.eq.jj)
        ii=i;jj=j;pivot=bcol1(ii,jj)
       end if
       !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)')'[',me,']{',p,'} col pivot ',ii,jj,kk,qq,' : ',pivot
      end if
      skipcol=.false.
      if(.not.done)then
!$OMP PARALLEL DO private(k,q)
       do kq=0,n(p+1)*r(p+1)-1
        q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1
        arow1(k,q)=dmrgg_fun(ii,jj,k,q, p, l,m,vip,r,n, fun, par)
       end do
!$OMP END PARALLEL DO
       nevalloc=nevalloc+n(p+1)*r(p+1)
       kq=idamax(n(p+1)*r(p+1),arow1,1)-1; q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1; amax=dmax1(amax,dabs(arow1(k,q)))
       haverow=.true.; crs=crs+1; done=havecol.and.haverow.and.(crs.ge.2*piv)
       if(.not.done)then
        call dcopy(n(p+1)*r(p+1),arow1,1,brow1,1)
        call dgemv('t',r(p),n(p+1)*r(p+1),-1.d0,row%u(p+1)%p,r(p),col%u(p)%p(ii,jj,1:r(p)),1,1.d0,brow1,1)
        kq=idamax(n(p+1)*r(p+1),brow1,1)-1; q=kq/n(p+1)+1; k=mod(kq,n(p+1))+1
        done=havecol.and.haverow.and. (k.eq.kk .and. q.eq.qq)
        qq=q;kk=k;pivot=brow1(kk,qq)
       end if
       !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)')'[',me,']{',p,'} row pivot ',ii,jj,kk,qq,' : ',pivot
      end if
     end do
     deallocate(bcol1,brow1,lot,b, stat=info)
     if(info.ne.0)then;write(*,*)'fail to deallocate after pivoting';stop;endif
    else
     write(*,*) subnam,': unknown pivoting: ',piv; stop
    end if !(pivoting)
    !write(*,'(a,i2,a,i4,a,4i4,a,e20.13)')'[',me,']{',p,'} fin pivot ',ii,jj,kk,qq,' : ',pivot
    !write(*,'(a,2(1x,i4),a,e20.13)')'PIVOT ',jj,kk,' : ',pivot

    tape(:,p)=(/-1,-1,-1,-1/)
    upd(p)=(dabs(pivot).gt.small_element*amax) .and. (dabs(pivot).gt.small_pivot*pivotmax_prev)
    if( upd(p) )then
     ! put new pivot in sets
     tape(:,p)=(/ ii,jj,kk,qq /)
     allocate(lot(4,r(p)+1),stat=info); if(info.ne.0)then;write(*,*)'allocate lot fail:',info;stop;endif
     forall(i=1:4,s=1:r(p))lot(i,s)=vip(p)%p(i,s)
     deallocate(vip(p)%p); allocate(vip(p)%p(4,r(p)+1),stat=info); if(info.ne.0)then;write(*,*)'reallocate vip fail:',info;stop;endif
     forall(i=1:4,s=1:r(p))vip(p)%p(i,s)=lot(i,s)
     vip(p)%p(:,r(p)+1)=tape(:,p)
     deallocate(lot)
     if(pivotmax.lt.0.d0)then;pivotmax=dabs(pivot);else;pivotmax=dmax1(pivotmax,dabs(pivot));endif
     if(pivotmin.lt.0.d0)then;pivotmin=dabs(pivot);else;pivotmin=dmin1(pivotmin,dabs(pivot));endif

      !write(*,'(a,i2,a,i2,a,4i5)')'[',me,']{',p,'}: reallocating superblock with sizes: ',r(p-1),n(p),n(p+1),r(p+1)
     allocate(bcol(r(p-1),n(p),r(p)+1),brow(r(p)+1,n(p+1),r(p+1)), bcol1(r(p-1),n(p)),brow1(n(p+1),r(p+1)),b(r(p)**2), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot reallocate blocks for update';stop;endif

     !if(me.eq.0)write(*,*)'updating inverse..'
     call dcopy(r(p)**2,inv(p)%p,1,b,1)
     deallocate(inv(p)%p); allocate(inv(p)%p((r(p)+1)**2),stat=info)
     call dcopy(r(p)**2,b,1,inv(p)%p,1)
     if(info.ne.0)then;write(*,'(a,i2,a,i5)')'[',me,']: cannot reallocate inverse: ',p;stop;endif
     call dcopy(r(p),col%u(p)%p(ii,jj,1:r(p)),1,inv(p)%p(r(p)*r(p)+1:r(p)*(r(p)+1)),1)
     call dcopy(r(p),row%u(p+1)%p(1:r(p),kk,qq),1,inv(p)%p(r(p)*(r(p)+1)+1:r(p)*(r(p)+2)),1)
     inv(p)%p((r(p)+1)**2)=pivot

     !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(arg): copying blocks'
     forall(i=1:r(p-1),j=1:n(p),s=1:r(p))bcol(i,j,s)=arg%u(p)%p(i,j,s)
     forall(i=1:r(p-1),j=1:n(p)         )bcol(i,j,r(p)+1)=acol1(i,j)
     forall(s=1:r(p),k=1:n(p+1),q=1:r(p+1))brow(s,k,q)=arg%u(p+1)%p(s,k,q)
     forall(         k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=arow1(k,q)
     !do s=1,r(p); do j=1,n(p); do i=1,r(p-1); bcol(i,j,s)=arg%u(p)%p(i,j,s); enddo; enddo; enddo
     !             do j=1,n(p); do i=1,r(p-1); bcol(i,j,r(p)+1)=acol1(i,j); enddo; enddo
     !do q=1,r(p+1);do k=1,n(p+1);do s=1,r(p); brow(s,k,q)=arg%u(p+1)%p(s,k,q);enddo;enddo;enddo
     !do q=1,r(p+1);do k=1,n(p+1);             brow(r(p)+1,k,q)=arow1(k,q);enddo;enddo
     deallocate(arg%u(p)%p, arg%u(p+1)%p); allocate(arg%u(p)%p(r(p-1),n(p),r(p)+1),arg%u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': not enough memory for new blocks';stop;endif
     call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,arg%u(p)%p,1)
     call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,arg%u(p+1)%p,1)

     !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(uv): updating factors'
     forall(i=1:r(p-1),j=1:n(p),s=1:r(p))bcol(i,j,s)=col%u(p)%p(i,j,s)
     forall(i=1:r(p-1),j=1:n(p)         )bcol(i,j,r(p)+1)=acol1(i,j)
     forall(s=1:r(p),k=1:n(p+1),q=1:r(p+1))brow(s,k,q)=row%u(p+1)%p(s,k,q)
     forall(         k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=arow1(k,q)
     call d2_lual(r(p-1)*n(p),  r(p)+1,inv(p)%p,bcol,from=r(p)+1)
     call d2_luar(n(p+1)*r(p+1),r(p)+1,inv(p)%p,brow,from=r(p)+1)
     deallocate(col%u(p)%p,row%u(p+1)%p); allocate(col%u(p)%p(r(p-1),n(p),r(p)+1),row%u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': not enough memory for new factors';stop;endif
     call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,col%u(p)%p,1)
     call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,row%u(p+1)%p,1)

     if(p.gt.own(me))then
      !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(v): updating left rows'
      call dcopy(r(p-1)*n(p)*r(p),row%u(p)%p,1,bcol,1)
      call dcopy(r(p-1)*n(p),arg%u(p)%p(:,:,r(p)+1),1,bcol1,1)
      call d2_luar(n(p),r(p-1),inv(p-1)%p,bcol1)
      call dcopy(r(p-1)*n(p),bcol1,1,bcol(1,1,r(p)+1),1)
      deallocate(row%u(p)%p); allocate(row%u(p)%p(r(p-1),n(p),r(p)+1), stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for left rows';stop;endif
      call dcopy(r(p-1)*n(p)*(r(p)+1),bcol,1,row%u(p)%p,1)
     end if
     if(p.lt.own(me+1)-1)then
      !write(*,'(a,i2,a,i2,a)') '[',me,']{',p,'}(u): updating right cols'
      forall(s=1:r(p),k=1:n(p+1),q=1:r(p+1))brow(s,k,q)=col%u(p+1)%p(s,k,q)
      forall(         k=1:n(p+1),q=1:r(p+1))brow1(k,q)=arg%u(p+1)%p(r(p)+1,k,q)
      call d2_lual(n(p+1),r(p+1),inv(p+1)%p,brow1)
      forall(         k=1:n(p+1),q=1:r(p+1))brow(r(p)+1,k,q)=brow1(k,q)
      deallocate(col%u(p+1)%p); allocate(col%u(p+1)%p(r(p)+1,n(p+1),r(p+1)),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': not enough memory for right cols';stop;endif
      call dcopy((r(p)+1)*n(p+1)*r(p+1),brow,1,col%u(p+1)%p,1)
     end if

     ! increase the rank and move on
     r(p)=r(p)+1
     deallocate(bcol,brow,bcol1,brow1,b,stat=info)
     if(info.ne.0)then;write(*,*)'fail to deallocate after update';stop;endif
    end if ! (update)
    deallocate(acol1,arow1)
   end do !(loop over my cores)

   pivotmax_prev = pivotmax
   nevalall = nevalloc

   ! REPORT current progress
   t2 = timef()
!   if(me.eq.0)write(str,'(i3,a2,a,f5.1,a,2f6.3)') &
!     it,sdir,' rank', erank(arg),&
!     ' log10t,n ',dlog(t2-t1)/dlog(10.d0),dlog(dble(nevalall))/dlog(10.d0)

   if(present(quad))then
    ttqq%l=l;ttqq%m=m;ttqq%r(l-1:m)=1;ttqq%n(l:m)=1
    ttqq%r(own(me)-1:own(me+1))=r(own(me)-1:own(me+1))
    call alloc(ttqq)
    first=own(me); last=own(me+1)-1; if(me.eq.nproc-1)last=m
    do p=first,last
     do k=1,r(p)
      call dgemv('n',r(p-1),n(p),1.d0,arg%u(p)%p(:,:,k),r(p-1),quad%u(p)%p,1,0.d0,ttqq%u(p)%p(1:r(p-1),1,k),1)
     end do
    end do
    call dtt_lua(ttqq,inv,own)
    val=dtt_quad(ttqq,mybonds=own)
    call dealloc(ttqq)

    ! print error or internal conv
    if(me.eq.0)then
     if(present(tru))then
      write(stmp,'(a,e8.3,a,'//sfmt//')')' err ',dabs(1.d0-val/tru),' val ',val
     else
      write(stmp,'(a,e8.3,a,'//sfmt//')')' cnv ',dabs(1.d0-val/val_prev),' val ',val
     end if
     str=trim(str)//trim(stmp)
    end if
    val_prev=val
   end if !(quad)

   ! check exit conditions
   if(present(maxrank))ready=ready.or.(it+1.ge.maxrank)
   if(present(accuracy))then
    if(pivotmax.le.accuracy*amax)then;strike=strike+1;else;strike=0;endif
    ready=ready.or.(strike.ge.3)
   end if
  end do  !(loop over ranks)

  call dtt_lua(arg,inv,own)

  do p=l-1,m
   deallocate(inv(p)%p,vip(p)%p)
  end do
  call dealloc(col); call dealloc(row)

  if(present(neval))then
    neval=nevalloc
  end if
 end subroutine dtt_dmrgg


 double precision function dmrgg_fun(i,j,k,q, p, l,m,vip,r,n, fun, par) result(f)
  implicit none
  integer,intent(in) :: i,j,k,q, p, l,m
  integer,intent(in) :: r(l-1:m),n(l:m)
  type(pointi2),intent(in) :: vip(0:tt_size)
  double precision,external :: fun
  double precision,intent(in),optional :: par(*)
  integer :: t,s, ind(tt_size)
  t=i; do s=p-1,l,-1; ind(s-l+1)=vip(s)%p(2,t); t=vip(s)%p(1,t); enddo
  ind(p  )=j; ind(p+1)=k
  t=q; do s=p+1,m-1,+1; ind(s+1-l+1)=vip(s)%p(3,t); t=vip(s)%p(4,t); enddo
  f=fun(m,ind,n,par)
 end function


 subroutine dtt_lua(arg,inv,own)
  implicit none
  type(dtt_tensor),intent(in) :: arg
  type(pointd) :: inv(0:tt_size)
  integer,intent(in) :: own(0:)

  character(len=*),parameter :: subnam='[dtt_lua]'
  integer,parameter :: tag=4
  integer :: me,nproc
  integer :: p,q,m, typed,info,r(0:tt_size),n(1:tt_size)

  ! MPI leftovers
  nproc = 1; me = 0

  r=arg%r; n=arg%n
 !if(me.eq.0)write(*,*)'applying inverses'
  do p=own(me),own(me+1)-1
   !write(*,'(a,i2,a,i2,a)')'[',me,']{',p,'} applying inverses from both sides'
   call d2_luar(n(p)*r(p),r(p-1),inv(p-1)%p,arg%u(p)%p)
   call d2_lual(r(p-1)*n(p),r(p),inv(p)%p,arg%u(p)%p)
  end do
  if(me.eq.nproc-1)then
   m=own(me+1)
   !write(*,'(a,i2,a,i2,a)')'[',me,']{',m,'} applying inverses from the left side'
   call d2_luar(n(m)*r(m),r(m-1),inv(m-1)%p,arg%u(m)%p)
  end if
 end subroutine dtt_lua

 double precision function dtt_quad(arg,quad,mybonds) result(val) ! val is returned only on proc=0
  implicit none
  type(dtt_tensor),intent(in),target :: arg
  type(dtt_tensor),intent(in),optional :: quad
  integer,intent(in),optional,target :: mybonds(0:)

  character(len=*),parameter :: subnam='[dtt_quad]'
  integer,parameter :: tagsize=11,tagdata=12
  integer,pointer :: r(:),n(:),own(:)
  integer :: me,her,him,nproc,info,typed, mym,myn,herm,hern
  integer :: first,last, l,m,p,q,i,j,k, mn(2)
  double precision,pointer :: prev(:,:),curr(:,:),next(:,:)

  val=0.d0
  nproc = 1; me = 0


  l=arg%l; m=arg%m; r=>arg%r; n=>arg%n
  ! distribute bonds (ranks) between procs
  if(present(mybonds))then
   own=>mybonds
  else
   allocate(own(0:nproc), stat=info)
   if(info.ne.0)then;write(*,*)subnam,': cannot allocate own bonds';stop;endif
   own(0) = l
   own(1) = m
  end if
  !write(*,'(a,i3,a,32i4)')'[',me,']: own: ',own(0:nproc)

  !if(me.eq.0)write(*,*)'compute local products'
  first=own(me)
  last=own(me+1)-1; if(me.eq.nproc-1)last=m
  do p=first,last
   allocate( curr(r(p-1), r(p)) )
   if(present(quad))then
    do k=1,r(p)
     call dgemv('n',r(p-1),n(p),1.d0,arg%u(p)%p(:,:,k),r(p-1),quad%u(p)%p,1,0.d0,curr(1:r(p-1),k),1)
    end do
   else
    forall(i=1:r(p-1),k=1:r(p))curr(i,k)=sum(arg%u(p)%p(i,:,k))
   end if
   if(p.eq.first)then
    prev=>curr; nullify(curr)
   else
    allocate(next(r(first-1),r(p)))
    call dgemm('n','n',r(first-1),r(p),r(p-1),1.d0,prev,r(first-1),curr,r(p-1),0.d0,next,r(first-1))
    deallocate(prev,curr)
    prev=>next; nullify(next)
   end if
   !write(*,'(a,i2,a,i2,a,2i5,a,2i5)')'[',me,']{',p,'} block shape ',shape(prev),' r(first-1) r(p): ',r(first-1),r(p)
  end do

  !write(*,'(a,i2,a)')'[',me,'] assembling along the binary tree'
  mym=r(first-1); myn=r(last)
  if(size(prev,1).ne.mym .or. size(prev,2).ne.myn)then
   write(*,'(a,i2,a,2i5,a,2i5)')'[',me,'] size mismatch: mym myn: ',mym,myn,' shape: ',shape(prev)
   stop
  endif
  q=1
  val=prev(1,1)
  deallocate(prev)
end function dtt_quad


end module
