!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file ttnodeop.f90
!! @note   Adopted from tt_fort library by Sergey Dolgov, sergey.v.dolgov@gmail.com 
!!         and Dmitry Savostyanov, dmitry.savostyanov@gmail.com
!! @contributor Ismael Djibrilla Boureima, Oleg Korobkin
!! @date   December 2024
!! @brief  ALS-related routines. Projections, truncations, enrichment, etc.
!!
module ttnodeop_lib
  use thor_lib, only : dtt_tensor
  use thor_pointers, only : pointd, pointd3, pointz, pointz3
  use ttlocsolve_lib
  use default_lib
  use trans_lib
  use mat_lib, only : chop, svd => svd_old
  implicit none

!! Author: Sergey Dolgov, sergey.v.dolgov@gmail.com
!! Contributor: Dmitry Savostyanov, dmitry.savostyanov@gmail.com
contains

! QR
  subroutine dtt_qr(arg,k,dir)
    use ort_lib
    implicit none
  ! dir=+1: QR left-to-right, u(k)->u(k+1) [default]
  ! dir=-1: QR right-to-left, u(k)->u(k-1)
    type(dtt_tensor),intent(inout) :: arg
    integer,intent(in) :: k
    integer,intent(in),optional :: dir
    character(len=*),parameter :: subnam='dtt_qr'
    integer :: info, m,p, rr, rnew, i,j, step
    double precision,allocatable :: R(:,:), U(:,:), nnod(:,:,:)

    step=default(+1,dir)
    if(abs(step).ne.1)then;write(*,*)subnam,': illegal dir: ',step;stop;endif

    ! allocate the auxiliary memory
    if(step.eq.+1)then
     m=arg%r(k-1)*arg%n(k); rr=arg%r(k); p=arg%n(k+1)*arg%r(k+1)
    else
     m=arg%n(k)*arg%r(k); rr=arg%r(k-1); p=arg%r(k-2)*arg%n(k-1)
    end if
    allocate(U(m,rr), R(rr,rr), nnod(arg%r(k+step-1),arg%n(k+step),arg%r(k+step)), stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
    call dcopy(p*rr,arg%u(k+step)%p,1,nnod,1)

    ! orthogonalize the left node u(k)
    if(step.eq.+1)then
     call dcopy(m*rr, arg%u(k)%p, 1, U, 1)
    else
     call trans2d(rr,m, arg%u(k)%p, U)
    end if
    call ort0(U,mat=R)

    ! compress nodes if necessary
    rnew=min(rr,m)
    if(rnew.lt.rr)then
     deallocate(arg%u(k)%p,arg%u(k+step)%p)
     if(step.eq.+1)then
      allocate(arg%u(k)%p(arg%r(k-1), arg%n(k), rnew), arg%u(k+1)%p(rnew, arg%n(k+1), arg%r(k+1)))
      arg%r(k)=rnew
     else
      allocate(arg%u(k-1)%p(arg%r(k-2), arg%n(k-1), rnew), arg%u(k)%p(rnew, arg%n(k), arg%r(k)))
      arg%r(k-1)=rnew
     end if
    end if

    ! copy orthogonal factor to node
    if(step.eq.+1)then
     call dcopy(m*rnew, U, 1, arg%u(k)%p, 1)
    else
     call trans2d(m,rnew, U, arg%u(k)%p)
    end if

    ! Apply R to the next node
    if(step.eq.+1)then
     call dgemm('n','n', rnew,p,rr, 1d0, R,rr, nnod,rr, 0d0, arg%u(k+step)%p,rnew)
    else
     call dgemm('n','t', p,rnew,rr, 1d0, nnod,p, R,rr, 0d0, arg%u(k+step)%p,p)
    end if
    deallocate(R,U,nnod)
  end subroutine

!RICH
  subroutine dtt_enrich(arg, k, rho, left,right)
   implicit none
   ! Increase r(k) -> r(k)+rho, add left and right arrays to cores k, k+1
   type(dtt_tensor),intent(inout) :: arg
   integer :: k,rho
   double precision,optional :: left(*),right(*)
   character(len=*),parameter :: subnam='dtt_enrich'
   integer :: a,b,c,info,rold,rnew
   double precision, allocatable :: knod(:,:,:),nnod(:,:,:)
   ! hold two nodes in a safe place
   allocate(knod(arg%r(k-1),arg%n(k),arg%r(k)), nnod(arg%r(k),arg%n(k+1),arg%r(k+1)),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
   call dcopy(arg%r(k-1)*arg%n(k)*arg%r(k),arg%u(k)%p,1,knod,1)
   call dcopy(arg%r(k)*arg%n(k+1)*arg%r(k+1),arg%u(k+1)%p,1,nnod,1)

   ! increase the rank and reallocate structures
   rold=arg%r(k); rnew=arg%r(k)+rho
   deallocate(arg%u(k)%p, arg%u(k+1)%p)
   allocate(arg%u(k)%p(arg%r(k-1),arg%n(k),rnew), arg%u(k+1)%p(rnew,arg%n(k+1),arg%r(k+1)),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': cannot allocate wider cores';stop;endif

   ! put old nodes in a new structure
   forall(a=1:arg%r(k-1),b=1:arg%n(k),c=1:arg%r(k)) arg%u(k)%p(a,b,c)=knod(a,b,c)
   forall(a=1:arg%r(k),b=1:arg%n(k+1),c=1:arg%r(k+1)) arg%u(k+1)%p(a,b,c)=nnod(a,b,c)

   ! put enrichments in a new structure
   if(present(left))then
    forall(a=1:arg%r(k-1),b=1:arg%n(k),c=1:rho) arg%u(k)%p(a,b,rold+c)=left(a+(b-1)*arg%r(k-1)+(c-1)*arg%r(k-1)*arg%n(k))
   else
    forall(a=1:arg%r(k-1),b=1:arg%n(k),c=1:rho) arg%u(k)%p(a,b,rold+c)=0.d0
   end if
   if(present(right))then
    forall(a=1:rho,b=1:arg%n(k+1),c=1:arg%r(k+1)) arg%u(k+1)%p(rold+a,b,c)=right(a+(b-1)*rho+(c-1)*rho*arg%n(k+1))
   else
    forall(a=1:rho,b=1:arg%n(k+1),c=1:arg%r(k+1)) arg%u(k+1)%p(rold+a,b,c)=0.d0
   end if

   arg%r(k)=rnew
   deallocate(nnod,knod)
  end subroutine


! TRUNC
  subroutine dtt_trunc(arg, k, tol, dir, rmax, A, phi, rhs, xnew)
    ! svd of u(k) to k+1 (dir==1), or k-1 (dir==-1)
    type(dtt_tensor),intent(inout) :: arg
    integer,intent(in) :: k, dir
    double precision,intent(in) :: tol
    integer,intent(in),optional :: rmax
    type(dtt_tensor),intent(in),optional :: A
    type(pointd),optional :: phi(0:) ! l-1:d
    double precision, optional :: rhs(:)
    double precision, pointer, optional :: xnew(:)
    character(len=*),parameter :: subnam='dtt_trunc'
    integer :: info, m,n,i, rr,rmin,rnew
    double precision,allocatable :: knod(:,:),nnod(:,:)
    double precision,pointer :: u(:,:), v(:,:), s(:)
    double precision :: err, norm_rhs
    double precision,external :: dnrm2

    if (dir.eq.+1) then
      m=arg%r(k-1)*arg%n(k)
      n=arg%n(k+1)*arg%r(k+1)
      rr=arg%r(k)
      rmin=min(rr,m)

      allocate(knod(m,rr),nnod(rr,n),stat=info)
      if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
      
      call dcopy(m*rr,arg%u(k)%p,1,knod,1)
      call svd(knod,u,v,s,info=info)
      if(info.ne.0)then 
       write(*,*)subnam,': svd info: ',info
       write(*,*)subnam,': return the system untouched, no truncation made'
       if(present(xnew))call dcopy(m*rr,arg%u(k)%p,1,xnew,1)
      else
       ! Multiply V'=SV'
       do i=1,rmin
         call dscal(rr,s(i),v(1,i),1)
       end do
      
       if ((present(A)).and.(present(phi))) then
       ! Residual truncation
         norm_rhs = dnrm2(m*rr,rhs,1)
      
         ! in the AMRs, the rank cannot change too much => its faster to start from end
         do rnew=rmin-1,1,-1
           call dgemm('n','t', m,rr,rnew, 1d0, u,m, v,rr, 0d0, arg%u(k)%p,m)
           ! A tilde x
           call d2d_mv(arg%r(k-1),arg%n(k),arg%r(k), arg%r(k-1),arg%n(k),arg%r(k), A%r(k-1), A%r(k), &
                       phi(k-1)%p, A%u(k)%p, phi(k)%p, arg%u(k)%p, knod)
           call daxpy(m*rr, -1d0, rhs,1, knod,1)
           err = dnrm2(m*rr, knod,1)/norm_rhs
           if (err>tol)exit
         end do
         rnew=rnew+1
       else
         ! Frobenius-norm truncation
         rnew = chop(S,tol)
       end if
      
       if(present(rmax))rnew=min(rnew,rmax)
      
       ! Reallocate two cores and fill with data
       call dcopy(rr*n,arg%u(k+1)%p,1,nnod,1)
       deallocate(arg%u(k)%p,arg%u(k+1)%p)
       allocate(arg%u(k)%p(arg%r(k-1),arg%n(k),rnew), arg%u(k+1)%p(rnew,arg%n(k+1),arg%r(k+1)))
       call dcopy(m*rnew, u, 1, arg%u(k)%p, 1)
       call dgemm('t','n', rnew,n,rr, 1d0, v,rr, nnod,rr, 0d0, arg%u(k+1)%p,rnew)
       if(present(xnew))call dgemm('n','t', m,rr,rnew, 1d0, u,m, v,rr, 0d0, xnew,m)
       arg%r(k)=rnew
       deallocate(u,v,s)
      end if
      deallocate(knod,nnod)
    else
     write(*,*) subnam,': not yet implemented dir: ',dir
     stop
    end if
  end subroutine



  !!!!!
  !! Partial dot products (projections) for ALS
  !!!!!
  subroutine dtt_YAX(phi_old, y, x, k, dir, phi_new, A)
  ! Computes the projection Y_{<=k}' [A_{<=k}] X_{<=k} (dir>0), or
  ! Y_{>=k}' [A_{>=k}] X_{>=k} (dir<0)
  ! A is optional
   implicit none
   double precision,intent(in) :: phi_old(*)
   type(dtt_tensor),intent(in) :: y, x
   type(dtt_tensor),intent(in),optional :: A
   integer,intent(in) :: k, dir
   double precision,pointer :: phi_new(:)
   character(len=*),parameter :: subnam='dtt_YAX'
   double precision,allocatable :: res1(:), res2(:)
   integer rx1, m, rx2, ry1, n, ry2, ra1, ra2

   rx1 = x%r(k-1); m = x%n(k); rx2 = x%r(k)
   ry1 = y%r(k-1); n = y%n(k); ry2 = y%r(k)
   if (present(A)) then
    ra1 = A%r(k-1); ra2 = A%r(k)
   else
    ra1 = 1; ra2 = 1;
   end if

   if (dir>0) then
     !  phi is of sizes ry1, rx1, ra1
     allocate(res1(max(rx1,rx2)*max(ra1*n, m*ra2)*ry2), res2(max(rx1,rx2)*max(ra1*n, m*ra2)*ry2))
     if (associated(phi_new))deallocate(phi_new)
     allocate (phi_new(ry2*rx2*ra2))

     call dgemm('T', 'N', rx1*ra1, n*ry2, ry1, 1d0, phi_old, ry1, y%u(k)%p, ry1, 0d0, res1, rx1*ra1)
     !    res1: rx1,ra1,n,ry2
     if (present(A)) then
      call trans2d(rx1, ra1*n*ry2, res1, res2)
      call dgemm('T', 'N', ry2*rx1, m*ra2, ra1*n, 1d0, res2, ra1*n, A%u(k)%p, ra1*n, 0d0, res1, ry2*rx1)
      !     res2: ry2,rx1,m,ra2
      call trans2d(ry2, rx1*m*ra2, res1, res2)
      call dgemm('T', 'N', ra2*ry2, rx2, rx1*m, 1d0, res2, rx1*m, x%u(k)%p, rx1*m, 0d0, res1, ra2*ry2)
      call trans2d(ra2, ry2*rx2, res1, phi_new)
      ! phi1: ra2, ry2, rx2
     else
      call dgemm('T', 'N', ry2, rx2, rx1*n, 1d0, res1, rx1*n, x%u(k)%p, rx1*n, 0d0, phi_new, ry2)
     end if

   else
     !  phi is of sizes rx2, ra2, ry2
     allocate(res1(rx1*max(ra1*n, m*ra2)*max(ry2,ry1)), res2(rx1*max(ra1*n, m*ra2)*max(ry2,ry1)))
     if (associated(phi_new))deallocate(phi_new)
     allocate (phi_new(rx1*ra1*ry1))

     call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x%u(k)%p, rx1*m, phi_old, rx2, 0d0, res1, rx1*m)
     !    res1: rx1,m,ra2,ry2
     if (present(A)) then
      call trans2d(rx1, m*ra2*ry2, res1, res2)
      call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A%u(k)%p, ra1*n, res2, m*ra2, 0d0, res1, ra1*n)
      !     res2: ra1,n,ry2,rx1
      call trans2d(ra1, n*ry2*rx1, res1, res2)
      call dgemm('T', 'T', rx1*ra1, ry1, n*ry2, 1d0, res2, n*ry2, y%u(k)%p, ry1, 0d0, phi_new, rx1*ra1)
     else
      call dgemm('N', 'T', rx1, ry1, n*ry2, 1d0, res1, rx1, y%u(k)%p, ry1, 0d0, phi_new, rx1)
     end if

   end if

   deallocate(res1, res2)
  end subroutine

end module
