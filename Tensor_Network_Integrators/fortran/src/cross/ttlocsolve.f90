module ttlocsolve_lib
  use trans_lib
  implicit none

  !!!!!
  !! Standard matrix routines (transposed, slices), and
  !! ALS-related local operations: matvecs, gmres, precs, etc.
  !!!!!

contains

  !!!!!
  !! Local ALS matvec using given phi1, A, phi2
  !!!!!
  subroutine d2d_mv(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y, work1, work2)
  ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
  integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
  double precision, intent(in) :: phi1(*), A(*), phi2(*), x(*)
  double precision, intent(out) :: y(*)
  double precision, intent(inout),pointer,optional :: work1(:), work2(:)
  !
  character(*),parameter :: subnam='[d2d_mv]'
  double precision, pointer :: res1(:), res2(:)
  integer :: info
  logical :: loc1,loc2

    loc1=.true.; if(present(work1)) loc1=.not.associated(work1)
    if (loc1) then
       allocate(res1(rx1*max(m*ra2, ra1*n)*ry2),stat=info)
       if(info.ne.0) error stop '[!]'//subnam//': cannot allocate res1'
    else
       res1=>work1
    endif

    loc2=.true.; if(present(work2)) loc2=.not.associated(work2)
    if(loc2)then
       allocate(res2(rx1*max(m*ra2, ra1*n)*ry2),stat=info)
       if(info.ne.0) error stop '[!]'//subnam//': cannot allocate res2'
    else
     res2=>work2
    endif

    !   phi2: rx2, ra2, ry2
    !  x: rx1,m,rx2
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2
    call trans2d(rx1, m*ra2*ry2, res1, res2)
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res2, m*ra2, 0d0, res1, ra1*n)
    !     res2: ra1,n,ry2,rx1
    call trans2d(ra1*n*ry2, rx1, res1, res2)
    !    phi1: ry1, rx1, ra1
    call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res2, rx1*ra1, 0d0, y, ry1)

    if(loc1)deallocate(res1)
    if(loc2)deallocate(res2)
  end subroutine


!!! Generate full local matrix
  subroutine d2d_fullmat(rx1, m, rx2, ra1, ra2, phi1, A, phi2, B)
  integer, intent(in) :: rx1, m, rx2, ra1, ra2
  double precision, intent(in) :: phi1(*), A(*), phi2(*)
  double precision, intent(inout) :: B(*)
  !
  double precision, allocatable :: res1(:), res2(:), phi2t(:)

    allocate(res1(rx1*m*rx1*m*max(ra2, rx2*rx2)), &
             res2(rx1*m*rx1*m*max(ra2, rx2*rx2)), phi2t(ra2*rx2*rx2))

    ! phi1: ry1,rx1,ra1
    call dgemm('N', 'N', rx1*rx1, m*m*ra2, ra1, 1d0, phi1, rx1*rx1, A, ra1, 0d0, res1, rx1*rx1)
    ! res1: ry1,rx1,n,m,ra2
    call dperm1324(rx1, rx1, m, m*ra2, res1, res2)
    ! phi2: rx2,ra2,ry2
    call trans2d(rx2, ra2*rx2, phi2, phi2t)
    call dgemm('N', 'N', rx1*m*rx1*m, rx2*rx2, ra2, 1d0, res2, rx1*m*rx1*m, phi2t, ra2, 0d0, res1, rx1*m*rx1*m);
    call dperm1324(rx1*m, rx1*m, rx2, rx2, res1, B)
    deallocate(res1, res2, phi2t)
  end subroutine


  subroutine d2d_jac_gen(ptype, rx1, n, rx2, ra1, ra2, Phi1, A, Phi2, jacs)
  ! sizes of work1, work2:
  !   center:
  !       rx*ra, rx1*n*n*rx2, rx2*ra2*rx2
  !   left:
  !       rx*ra, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2, rx2*ra2*rx2
  !   right:
  !       rx*ra, rx1*n*n*rx2*rx2, rx2*ra2*rx2
  ! sizes of jacs:
  !   center:
  !       rx1*n*n*ra2,  rx1*n*n*rx2
  !   left:
  !       rx1*rx1*n*n*ra2,  rx1*rx1*n*n*rx2
  !   right:
  !       n*n*rx2*rx2*rx1, rx1*n*n*ra2
  character, intent(in):: ptype
  integer, intent(in) :: rx1, n, rx2, ra1, ra2
  double precision, intent(in) :: Phi1(*), A(*), Phi2(*)
  double precision, intent(inout),pointer :: jacs(:) !, work1(*), work2(*)
  !
  double precision, allocatable :: work1(:), work2(:)
  integer, allocatable :: ipiv(:)
  integer :: i, info

    if ((ptype=='c').or.(ptype=='C')) then

       i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2)
       allocate(work1(i))
       allocate(work2(i))
       allocate(ipiv(n))

       call trans2d(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do
       do i=0,ra1-1
          call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do

       call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
       call dgemm('N', 'T', rx1*n*n, rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2, 0d0, work1, rx1*n*n)
       call trans2d(rx1, n*n*rx2, work1, work2)
       ! inversion
       do i=1,n*n
          work1(i)=0d0
       end do
       do i=1,n
          work1(i+(i-1)*n)=1d0
       end do
       do i=0,rx2*rx1-1
          call dcopy(n*n, work1, 1, jacs(i*n*n+1:(i+1)*n*n), 1)
          call dgesv(n, n, work2(i*n*n+1), n, ipiv, jacs(i*n*n+1), n, info)
       end do
    end if

    if ((ptype=='l').or.(ptype=='L')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx1*n
       allocate(ipiv(i))

       call trans2d(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do

       call dgemm('N', 'N', rx1*rx1, n*n*ra2, ra1, 1d0, Phi1, rx1*rx1, A, ra1, 0d0, work1, rx1*rx1)
       call dperm1324(rx1, rx1, n, n*ra2, work1, jacs)
       call dgemm('N', 'T', rx1*n*rx1*n, rx2, ra2, 1d0, jacs, rx1*n*rx1*n, work2, rx2, 0d0, work1, rx1*n*rx1*n)
       ! inversion
       do i=1,rx1*n*rx1*n
          work2(i)=0d0
       end do
       do i=1,rx1*n
          work2(i+(i-1)*rx1*n)=1d0
       end do
       do i=0,rx2-1
          call dcopy(rx1*n*rx1*n, work2, 1, jacs(i*rx1*n*rx1*n+1:(i+1)*rx1*n*rx1*n), 1)
          call dgesv(rx1*n, rx1*n, work1(i*rx1*n*rx1*n+1), rx1*n, ipiv, jacs(i*rx1*n*rx1*n+1), rx1*n, info)
          if(info.ne.0)then;write(*,*)': dgesv problem ', info;endif
       end do
    end if

    if ((ptype=='r').or.(ptype=='R')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx2*n
       allocate(ipiv(i))

       call trans2d(rx2*ra2, rx2, Phi2, work2)
       do i=0,ra1-1
          call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do

       call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
       call dgemm('N', 'T', rx1*n*n, rx2*rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2*rx2, 0d0, work1, rx1*n*n)
       call dperm1324(rx1*n, n, rx2, rx2, work1, work2)
       call trans2d(rx1, n*rx2*n*rx2, work2, work1)
       ! inversion
       do i=1,rx2*n*rx2*n
          work2(i)=0d0
       end do
       do i=1,rx2*n
          work2(i+(i-1)*rx2*n)=1d0
       end do
       do i=0,rx1-1
          call dcopy(rx2*n*rx2*n, work2, 1, jacs(i*rx2*n*rx2*n+1:(i+1)*rx2*n*rx2*n), 1)
          call dgesv(rx2*n, rx2*n, work1(i*rx2*n*rx2*n+1), rx2*n, ipiv, jacs(i*rx2*n*rx2*n+1), rx2*n, info)
       end do
    end if

    deallocate(work1,work2, ipiv)
  end subroutine


  subroutine d2d_jac_apply(ptype, rx1, n, rx2, jacs, x, y, res1)
    ! sizes of work1: rx1*n*rx2
    character, intent(in) :: ptype
    integer, intent(in) :: rx1, n, rx2
    double precision, intent(in) :: jacs(*), x(*)
    double precision, intent(inout) :: y(*)
    double precision, intent(inout),pointer,optional :: res1(:)
    double precision, pointer :: work1(:)
    integer i, info
    logical :: loc1

    loc1=.true.; if(present(res1))then;if(associated(res1))loc1=.false.;endif
    if(loc1)then
     allocate(work1(rx1*n*rx2),stat=info)
     if(info.ne.0)then;write(*,*)': cannot allocate work1';stop;endif
    else
     work1=>res1
    endif

    if ((ptype=='c').or.(ptype=='C')) then
       ! jacs is n,n,rx2,rx1
       call trans2d(rx1, n*rx2, x, work1)
       do i=0,(rx2*rx1-1)
          call dgemv('N', n, n, 1d0, jacs(i*n*n+1), n, work1(i*n+1), 1, 0d0, y(i*n+1), 1)
       end do
       call trans2d(n*rx2, rx1, y, work1)
       call dcopy(rx1*n*rx2, work1, 1, y, 1)
    end if

    if ((ptype=='l').or.(ptype=='L')) then
       ! jacs is rx1*n,rx1*n, rx2
       call dcopy(rx1*n*rx2, x, 1, work1, 1)
       do i=0,rx2-1
          call dgemv('N', rx1*n, rx1*n, 1d0, jacs(i*rx1*n*rx1*n+1), rx1*n, work1(i*rx1*n+1), 1, 0d0, y(i*rx1*n+1), 1)
       end do
    end if

    if ((ptype=='r').or.(ptype=='R')) then
       ! jacs is n*rx2, n*rx2, rx1
       call trans2d(rx1, n*rx2, x, work1)
       do i=0,rx1-1
          call dgemv('N', n*rx2, n*rx2, 1d0, jacs(i*n*rx2*n*rx2+1), n*rx2, work1(i*n*rx2+1), 1, 0d0, y(i*n*rx2+1), 1)
       end do
       call trans2d(n*rx2, rx1, y, work1)
       call dcopy(rx1*n*rx2, work1, 1, y, 1)
    end if

    if(loc1)deallocate(work1)
  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!
!!!!! GMRES !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
  subroutine d2d_gmresr(Phi1, A, Phi2, rhs, rx1, n, rx2, ra1, ra2, nrestart, tol, niters, ptype, jacs, sol, verb)
    ! right preconditioned - for residual tolerance
    ! This one with Householder tranforms
    double precision,intent(in):: Phi1(*), A(*), Phi2(*), rhs(*)
    double precision,intent(in), pointer :: jacs(:)
    integer, intent(in) :: rx1,n,rx2, ra1,ra2, nrestart, niters, verb
    double precision, intent(in) :: tol
    character, intent(in) :: ptype
    double precision, intent(inout) :: sol(*)
    character(len=*),parameter :: subnam='d2d_gmresr'

    integer :: i,j,it, sz, nrestart0,info
    double precision, pointer :: U(:,:), w(:), R(:,:), JJ(:,:), tau(:), res1(:), res2(:)
    double precision  nrmr, curres, nrmrhs, dbeta, dalpha

    double precision dnrm2
    double precision ddot


    sz = rx1*n*max(ra2,ra1)*rx2;
    allocate(res1(sz), res2(sz),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot alloc res12';stop;endif

    sz = rx1*n*rx2;

    nrestart0 = min(sz, nrestart)

    allocate(U(sz, nrestart0+1), w(sz), tau(nrestart0+1), JJ(2,nrestart0+1), R(nrestart0+1, nrestart0+1))

    do j=1,sz
       sol(j)=0d0
    end do

    do it=0,niters-1
       ! r0
       if (.not.((ptype=='n').or.(ptype=='N'))) then
          call d2d_jac_apply(ptype, rx1, n, rx2, jacs, sol, w, res1=res1)
          call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2)
       end if
       if ((ptype=='n').or.(ptype=='N')) then
          call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, sol, w, work1=res1,work2=res2)
       end if
       call daxpy(sz,-1d0,rhs,1,w,1)

       nrmr = dnrm2(sz,w,1);
       if (it==0) then
          nrmrhs = nrmr
       end if
       if (nrmr==0d0) then
          exit
       end if
       ! initial HHT
       dbeta = nrmr;
       if (w(1)<0d0) then
          dbeta = -dbeta
       end if
       w(1) = w(1)+dbeta
       tau(1) = -dbeta
       nrmr = dnrm2(sz,w,1)
       dbeta = 1d0/nrmr
       call dscal(sz,dbeta,w,1)
       call dcopy(sz,w,1,U,1)

       do j=1,nrestart0
          ! HHT on last U
          call dcopy(sz,U(1:sz,j),1,w,1)
          dbeta = -2d0*U(j,j)
          call dscal(sz,dbeta,w,1)
          w(j) = w(j) + 1d0
          do i=j-1,1,-1
             dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
             call daxpy(sz,dbeta,U(1,i),1,w,1);
          end do
          dbeta = dnrm2(sz,w,1);
          dbeta = 1d0/dbeta;
          call dscal(sz,dbeta,w,1);

          ! precvec, matvec
          if (.not.((ptype=='n').or.(ptype=='N'))) then
             call d2d_jac_apply(ptype, rx1, n, rx2, jacs, w, w, res1=res1);
             call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2);
          end if
          if ((ptype=='n').or.(ptype=='N')) then
             call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2);
          end if

          ! Orthog w to j projectors
          do i=1,j
             dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
             call daxpy(sz,dbeta,U(1,i),1,w,1);
          end do

          ! new P_{j+1}
          if (j<sz) then
             do i=1,j
                U(i,j+1)=0d0;
             end do
             i = sz-j;
             call dcopy(i, w(j+1:j+i), 1, U(j+1:j+i, j+1), 1);
             dalpha = dnrm2(i, U(j+1:j+i, j+1), 1);
             if (.not.(dalpha==0d0)) then
                if (w(j+1)<0d0) then
                   dalpha = -dalpha;
                end if
                U(j+1, j+1) = U(j+1,j+1) + dalpha;
                dbeta = dnrm2(i, U(j+1:j+i, j+1), 1);
                dbeta = 1d0/dbeta;
                call dscal(i,dbeta,U(j+1, j+1),1);

                w(j+1) = -dalpha;
                do i=j+2,sz
                   w(i)=0d0;
                end do
             end if
          end if

          ! Givens rotators to the top of w
          do i=1,j-1
             dbeta = w(i);
             w(i) = JJ(1,i)*w(i) + JJ(2,i)*w(i+1);
             w(i+1) = -JJ(2,i)*dbeta + JJ(1,i)*w(i+1);
          end do

          ! New rotator
          if (j<sz) then
             dalpha = dnrm2(2, w(j:j+1), 1)
             JJ(1,j) = w(j)/dalpha;
             JJ(2,j) = w(j+1)/dalpha;
             tau(j+1) = -JJ(2,j)*tau(j);
             tau(j) = JJ(1,j)*tau(j);
             w(j) = dalpha;
             w(j+1) = 0d0;
          end if

          call dcopy(j, w, 1, R(1:j,j), 1);

          ! residual
          curres = dabs(tau(j+1))/nrmrhs;

          if (curres<tol) then
             exit
          end if
       end do ! inner

       if (j>nrestart0) then
          j = nrestart0
       end if

       call dtrsv('u','n','n',j,R,nrestart0+1,tau,1);

       ! Correction
       call dcopy(sz, U(1:sz,j), 1, w, 1);
       dbeta = -2d0*U(j,j)*tau(j);
       call dscal(sz, dbeta, w, 1);
       w(j) = w(j) + tau(j);
       do i=j-1,1,-1
          w(i) = w(i)+tau(i);
          dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
          call daxpy(sz,dbeta,U(1,i),1,w,1);
       end do
       dalpha=-1d0;
       call daxpy(sz,dalpha,w,1,sol,1);
       if (curres<tol) then
        exit;
       end if
    end do ! iters

    if (verb>0) then
       write(*,"(A,I0,A,I0,A,ES10.3)") 'gmres conducted[', it, ',', j, '] iters to relres ', curres
    end if

    deallocate(U,w,tau,JJ,R,res1,res2)
  end subroutine


end module
