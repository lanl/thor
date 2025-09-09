module ttamen_lib
!
!   This file contains AMEn solver for SPD problems in higher dimensions described in the following reference
!
!   Alternating minimal energy methods for linear systems in higher dimensions.
!   Part I: SPD systems, http://arxiv.org/abs/1301.6068,
!   Part II: Faster algorithm and application to nonsymmetric systems, http://arxiv.org/abs/1304.1222
!   Other experimental/development routines for approximation, eigenproblem, etc. are also included.
!   Ispiered from tt-fort/ttamen.f90
!
  use thor_lib, only : tt_size, dtt_tensor, dtt_random_ortho, dealloc_dtt_tensor
  use thor_pointers, only : pointd, pointz !, pointd3, gpointd3, pointz, pointz3, gpointz3
  use ttlocsolve_lib
  use ttnodeop_lib
  use default_lib

  implicit none
  double precision,parameter :: resid_damp=2d0
  interface amen_solve
      module procedure dtt_amen_solve
  end interface
contains

! AMEn Solver
 function dtt_amen_solve(A, y, tol, kickrank, nswp, local_prec, local_iters, local_restart, trunc_norm, max_full_size, verb_) result(x)
  implicit none
  class(dtt_tensor),intent(in),target :: A
  type(dtt_tensor),intent(in) :: y
  double precision,intent(in) :: tol
  type(dtt_tensor) :: x
  integer,intent(in), optional :: kickrank, local_iters, local_restart, nswp, trunc_norm, verb_, max_full_size ! trunc_norm: 0 - fro, 1 - resid
  character,intent(in), optional :: local_prec ! 'n', 'l', 'c', 'r'
  character(len=*),parameter :: subnam='amen_solve'
  integer :: kickrank0, local_iters0, local_restart0, nswp0, trunc_norm0, verb0, max_full_size0
  character :: local_prec0
  type(dtt_tensor) :: z
  type(pointd) :: phixax(0:tt_size), phixy(0:tt_size), phizax(0:tt_size), phizy(0:tt_size)
  integer :: i, j, l, m, d, swp, rx1, rx2, ra1, ra2, info
  integer, pointer :: n(:)
  integer, allocatable :: ipiv(:), xrr(:)
  double precision,allocatable :: rhs(:), resid(:), B(:)
  double precision,pointer :: jacs(:), sol(:)
  double precision :: err, err_max, res_max, real_tol
  logical :: prc, jacsused
  double precision,external :: dnrm2
  logical,external :: lsame
  target :: y,x,z

  double precision, allocatable :: zarr(:)

  kickrank0 = default(4, kickrank)
  nswp0 = default(20, nswp)
  local_prec0 = default('n',local_prec)
  local_iters0 = default(2, local_iters)
  local_restart0 = default(40, local_restart)
  trunc_norm0 = default(1, trunc_norm)
  max_full_size0 = default(250, max_full_size)
  verb0 = default(0, verb_)
  prc=.not.lsame(local_prec0,'n')

  l= y%l; m= y% m
  d= m - l + 1
  n => y%n
  ! initial guess for x
  allocate(xrr(0:d)); xrr(0)=1; xrr(1:d-1)=2; xrr(d)=1
  x= dtt_random_ortho(A%n(l:m)/y%n(l:m), r_=xrr)
  deallocate(xrr)

  real_tol = tol/dsqrt(1.d0*(d-l))
  real_tol = real_tol/resid_damp

  ! prepare projections
  allocate(phixax(l-1)%p(1), phixy(l-1)%p(1), phixax(d)%p(1), phixy(d)%p(1))
  phixax(l-1)%p(1) = 1d0; phixax(d)%p(1) = 1d0
  phixy(l-1)%p(1) = 1d0; phixy(d)%p(1) = 1d0

  if (kickrank0.gt.0) then
  ! prepare z
    !z%l = l; z%m = d; z%n(l:d)=n(l:d); z%r(l-1)=1; z%r(l:d-1)=kickrank0; z%r(d)=1
    !call random(z)
    z= dtt_random_ortho(n(l:d),rall_=kickrank0)
    allocate(phizax(l-1)%p(1), phizy(l-1)%p(1), phizax(d)%p(1), phizy(d)%p(1))
    phizax(l-1)%p(1) = 1d0; phizax(d)%p(1) = 1d0
    phizy(l-1)%p(1) = 1d0; phizy(d)%p(1) = 1d0
  end if

  ! Main cycle
  do swp=1,nswp0
    ! orthogonalization
    do i=d,l+1,-1
      if ((kickrank0.gt.0).and.(swp.gt.1)) then
        ! update Z block
        call d2d_mv(x%r(i-1),n(i),x%r(i), z%r(i-1),n(i),z%r(i), A%r(i-1),A%r(i), phizax(i-1)%p, &
                    A%u(i)%p, phizax(i)%p, x%u(i)%p, z%u(i)%p)

        allocate(rhs(z%r(i-1)*n(i)*z%r(i)), resid(z%r(i-1)*n(i)*max(z%r(i),y%r(i))))
        call dgemm('n','n', z%r(i-1), n(i)*y%r(i), y%r(i-1), 1d0, phizy(i-1)%p, z%r(i-1), y%u(i)%p, y%r(i-1), 0d0, resid, z%r(i-1))
        call dgemm('n','n', z%r(i-1)*n(i), z%r(i), y%r(i), 1d0, resid, z%r(i-1)*n(i), phizy(i)%p, y%r(i), 0d0, rhs, z%r(i-1)*n(i))
        call daxpy(z%r(i-1)*n(i)*z%r(i), -1d0, rhs, 1, z%u(i)%p, 1)
        ! now z(i) = Z'(Ax-y)
        deallocate(rhs, resid)
      end if

      call dtt_qr(x, i, dir=-1)

      ! update X projections
      call dtt_YAX(phixax(i)%p, x, x, i, -1, phixax(i-1)%p, A=A)
      call dtt_YAX(phixy(i)%p, x, y, i, -1, phixy(i-1)%p)

      if (kickrank0>0) then
        call dtt_qr(z, i, dir=-1)
        ! update Z projections
        call dtt_YAX(phizax(i)%p, z, x, i, -1, phizax(i-1)%p, A=A)
        call dtt_YAX(phizy(i)%p, z, y, i, -1, phizy(i-1)%p)
      end if
    end do ! ort

    res_max = 0d0; err_max = 0d0
    ! optimization
    do i=l,d
      ! extract sizes
      rx1 = x%r(i-1); rx2 = x%r(i)
      ra1 = A%r(i-1); ra2 = A%r(i)

      ! prepare rhs
      ! first, memory allocations
      if (kickrank0.gt.0) then
        allocate(rhs(max(rx1,z%r(i-1))*n(i)*max(rx2,z%r(i))), resid(max(rx1,z%r(i-1))*n(i)*max(rx2,y%r(i),z%r(i))))
        allocate(sol(max(rx1,z%r(i-1))*n(i)*max(rx2,z%r(i))))
      else
        allocate(rhs(rx1*n(i)*rx2), resid(rx1*n(i)*max(rx2,y%r(i))))
        allocate(sol(rx1*n(i)*rx2))
      end if

      call dgemm('n','n', rx1, n(i)*y%r(i), y%r(i-1), 1d0, phixy(i-1)%p, rx1, y%u(i)%p, y%r(i-1), 0d0, resid, rx1)
      call dgemm('n','n', rx1*n(i), rx2, y%r(i), 1d0, resid, rx1*n(i), phixy(i)%p, y%r(i), 0d0, rhs, rx1*n(i))

      ! res_prev
      call d2d_mv(rx1, n(i), rx2, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, x%u(i)%p, resid)
      call daxpy(rx1*n(i)*rx2, -1d0, rhs, 1, resid, 1)
      err = dnrm2(rx1*n(i)*rx2, resid, 1)/dnrm2(rx1*n(i)*rx2, rhs, 1)
      res_max = dmax1(res_max,err)
      jacsused = .false.

      if (err.gt.real_tol) then
        if (rx1*n(i)*rx2<max_full_size0) then ! full solution
          allocate(B(rx1*n(i)*rx2*rx1*n(i)*rx2), ipiv(rx1*n(i)*rx2),stat=info)
          if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
          call d2d_fullmat(rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, B)
          call dcopy(rx1*n(i)*rx2, resid, 1, sol, 1)
          call dgesv(rx1*n(i)*rx2, 1, B, rx1*n(i)*rx2, ipiv, sol, rx1*n(i)*rx2, info)
          if(info.ne.0)then; write(*,*)subnam,': dgesv info: ',info; stop; end if
          if (verb0>1) write(*,"(A,I0,A)") 'amen_solve: block=', i, ', direct local solver'
          deallocate(B, ipiv)
        else ! iter. solution
          ! allocate and compute jac prec, if necessary
          if (lsame(local_prec0,'l')) then
            allocate(jacs(rx1*rx1*n(i)*n(i)*max(ra2, rx2)))
            jacsused = .true.
          else if (lsame(local_prec0,'c'))then
            allocate(jacs(rx1*n(i)*n(i)*max(ra2, rx2)))
            jacsused = .true.
          else if (lsame(local_prec0,'r'))then
            allocate(jacs(rx1*n(i)*n(i)*rx2*max(ra2, rx2)))
            jacsused = .true.
          end if
          if (prc) call d2d_jac_gen(local_prec0, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, jacs)

          ! Run the solver
          call d2d_gmresr(phixax(i-1)%p, A%u(i)%p, phixax(i)%p, resid, rx1, n(i), rx2, ra1, ra2, &
                          local_restart0, real_tol/err, local_iters0, local_prec0, jacs, sol, max(verb0-1,0))

          if (prc) call d2d_jac_apply(local_prec0, rx1, n(i), rx2, jacs, sol, sol)
        end if
        ! Now, sol is the -correction, since the rhs was resid = Ax-y
        call daxpy(rx1*n(i)*rx2, -1d0, sol, 1, x%u(i)%p, 1)
        ! the updated core is ready
        err = dnrm2(rx1*n(i)*rx2, sol, 1)/dnrm2(rx1*n(i)*rx2, x%u(i)%p, 1)
        err_max = dmax1(err_max,err)
        ! compute res_new
        call d2d_mv(rx1, n(i), rx2, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, x%u(i)%p, resid)
        call daxpy(rx1*n(i)*rx2, -1d0, rhs, 1, resid, 1)
        err = dnrm2(rx1*n(i)*rx2, resid, 1)/dnrm2(rx1*n(i)*rx2, rhs, 1)
      end if


      if (i.lt.d) then
        if (kickrank0.gt.0) then
          if (trunc_norm0.eq.1) then ! truncation in the A-norm of error (i.e. residual)
            err = dmax1(err,real_tol*resid_damp)
            call dtt_trunc(x, i, err, 1, A=A, phi=phixax, rhs=rhs, xnew=sol)
          else  ! truncation in the frobenius norm of solution
            call dtt_trunc(x, i, real_tol*resid_damp, 1, xnew=sol)
          end if
          call dcopy(rx1*n(i)*rx2, sol, 1, rhs, 1)

          ! prepare the enrichments: both for X and for Z
          call d2d_mv(rx1, n(i), rx2, rx1, n(i), z%r(i), ra1, ra2, phixax(i-1)%p, A%u(i)%p, phizax(i)%p, rhs, sol)
          call d2d_mv(rx1,n(i),rx2, z%r(i-1),n(i),z%r(i), ra1,ra2, phizax(i-1)%p, A%u(i)%p, phizax(i)%p, rhs, z%u(i)%p)

          call dgemm('n','n', rx1, n(i)*y%r(i), y%r(i-1), 1d0, phixy(i-1)%p, rx1, y%u(i)%p, y%r(i-1), 0d0, resid, rx1)
          call dgemm('n','n', rx1*n(i), z%r(i), y%r(i), 1d0, resid, rx1*n(i), phizy(i)%p, y%r(i), 0d0, rhs, rx1*n(i))
          call daxpy(rx1*n(i)*z%r(i), -1d0, rhs, 1, sol, 1)
          ! now sol is the enrichment to x
          ! update the block for z
          call dgemm('n','n', z%r(i-1), n(i)*y%r(i), y%r(i-1), 1d0, phizy(i-1)%p, z%r(i-1), &
                              y%u(i)%p, y%r(i-1), 0d0, resid, z%r(i-1))
          call dgemm('n','n', z%r(i-1)*n(i), z%r(i), y%r(i), 1d0, resid, z%r(i-1)*n(i), phizy(i)%p, y%r(i), 0d0, rhs, z%r(i-1)*n(i))
          call daxpy(z%r(i-1)*n(i)*z%r(i), -1d0, rhs, 1, z%u(i)%p, 1)
          ! now z(i) = Z'(Ax-y)

          ! enrichment
          call dtt_enrich(x, i, z%r(i), left=sol)
          ! orthogonalization
          call dtt_qr(x, i, dir=+1)
          call dtt_qr(z, i, dir=+1)
          ! update projections
          call dtt_YAX(phizax(i-1)%p, z, x, i, 1, phizax(i)%p, A=A)
          call dtt_YAX(phizy(i-1)%p, z, y, i, 1, phizy(i)%p)
        else
          call dtt_qr(x, i, dir=+1)
        end if

        ! update projections
        call dtt_YAX(phixax(i-1)%p, x, x, i, 1, phixax(i)%p, A=A)
        call dtt_YAX(phixy(i-1)%p, x, y, i, 1, phixy(i)%p)

      end if ! i<d

      ! deallocate all work arrays
      if (prc)then
        if (jacsused)deallocate(jacs)
      end if
      deallocate(rhs, resid, sol)
    end do ! i

    if (verb0>0) then
!#ifdef MEXWRITE
!      write(mxline,"(A,I0,A,ES10.3,A,ES10.3,A,I0,A)"), 'amen_solve: swp=', swp, ', max_dx=', err_max, ', max_res=', res_max, ', max_rank=', maxval(x%r(l-1:d)), char(0)
!      call mexputstr()
!#else
      write(*,"(A,I0,A,ES10.3,A,ES10.3,A,I0)") 'amen_solve: swp=', swp, ', max_dx=', err_max, &
                                               ', max_res=', res_max, ', max_rank=', maxval(x%r(l-1:d))
!#endif
    end if

    if (trunc_norm0.eq.1) then
      if (res_max.lt.tol) exit
    else
      if (err_max.lt.tol) exit
    end if

  end do ! swp

  do i=l-1,d
    deallocate(phixax(i)%p, phixy(i)%p)
    if (kickrank0>0)deallocate(phizax(i)%p, phizy(i)%p)
  end do
  if (kickrank0>0)call dealloc_dtt_tensor(z)
 end function dtt_amen_solve



 !> Approximate the matrix-by-vector via the AMEn iteration
 !!
 !!   [y,z]=amen_mv(A, x, tol, varargin)
 !!   Attempts to approximate the y = A*x
 !!   with accuracy TOL using the AMEn+ALS iteration.
 !!   Matrix A has to be given in the TT-format, right-hand side x should be
 !!   given in the TT-format also.
 !!
 !!   Options are provided in form
 !!   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
 !!   on. The parameters are set to default (in brackets in the following)
 !!   The list of option names and default values are:
 !!       o y0 - initial approximation to Ax [rand rank-2]
 !!       o nswp - maximal number of sweeps [20]
 !!       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
 !!       o kickrank - compression rank of the error,
 !!         i.e. enrichment size [3]
 !!       o init_qr - perform QR of the input (save some time in ts, etc) [true]
 !!       o renorm - Orthog. and truncation methods: direct (svd,qr) or gram
 !!         (apply svd to the gram matrix, faster for m>>n) [direct]
 !!       o fkick - Perform solution enrichment during forward sweeps [false]
 !!         (rather questionable yet; false makes error higher, but "better
 !!         structured": it does not explode in e.g. subsequent matvecs)
 !!       o z0 - initial approximation to the error Ax-y [rand rank-kickrank]
 !!
 !!
 !!********
 !!   For description of adaptive ALS please see
 !!   Sergey V. Dolgov, Dmitry V. Savostyanov,
 !!   Alternating minimal energy methods for linear systems in higher dimensions.
 !!   Part I: SPD systems, http://arxiv.org/abs/1301.6068,
 !!   Part II: Faster algorithm and application to nonsymmetric systems, http://arxiv.org/abs/1304.1222
 !!
 !!   Use {sergey.v.dolgov, dmitry.savostyanov}@gmail.com for feedback
 !!********
 !!
 !!  Adopted from: core/amen_mv.m in TT-Toolbox 2.2, 2009-2012
 !!
 !!---------------------------
 !! TODO:
 !!  - estimate size of the matrices and allocate things in advance;
 !!  - add arguments to thor_matmul to allow multiplication without reshaping (saves on copy)
 subroutine amen_mv(y, A, x, tol, &
    y0_, nswp_, verb_, kickrank_, init_qr_, renorm_, fkick_)
 use matlab_struct_module, only: array2d, array3d, array4d, rand_cell3d
 use mat_lib, only: svd, svdgram, bfun3, my_chop2, stack, &
                    diag => matlab_diag, zeros => matlab_zeros, &
                    qr => matlab_qr, norm => normfro
 use rnd_lib, only: rndmat
 use thor_lib, only: dtt_tensor, dtt_matrix, core2cell, dtt_random_ortho
 use mat_lib, only: thor_matmul
 use time_lib
 implicit none
 class(dtt_tensor), intent(INOUT) :: y
 class(dtt_matrix), intent(IN) :: A
 class(dtt_tensor), intent(IN) :: x
 double precision, intent(IN):: tol
 class(dtt_tensor), intent(IN), optional:: y0_!< y0 - initial approximation to Ax [rand rank-2]
 integer, intent(IN),optional :: nswp_        !< maximal number of sweeps [20]
 integer, intent(IN),optional :: verb_        !< verbosity level, 0-silent, 1-sweep info, 2-block info [1]
 integer, intent(IN),optional :: kickrank_    !< compression rank of the error,
                                              !  i.e. enrichment size [3]
 logical, intent(IN),optional :: init_qr_     !< perform QR of the input (save some time in ts, etc) [true]
 integer, intent(IN),optional :: renorm_      !< renorm - Orthog. and truncation methods: 1:direct (svd,qr) or 2:gram
                                              !  (apply svd to the gram matrix, faster for m>>n) [direct]
 logical, intent(IN),optional :: fkick_       !< Perform solution enrichment during forward sweeps [false]
                                              !  (rather questionable yet; false makes error higher, but "better
                                              !  structured": it does not explode in e.g. subsequent matvecs)
 !
 type(array3d), allocatable:: z(:) !< z0 - initial approximation to the error Ax-y [rand rank-kickrank]
 integer, allocatable :: m(:), rx(:), n(:), ra(:), ry(:), rz(:), r2(:)
 integer :: d,dir,i,r,swp
 integer :: nswp,verb,kickrank,renorm
 logical :: init_qr,fkick
 integer :: kickrank2, Ntim
 type(array3d), allocatable :: phiax(:), phiyax(:), phizax(:), phizy(:)
 type(array2d), allocatable :: y2(:), A2(:)
 double precision :: dx,max_dx,nrm,nrm_s,nrmr,nrmz
 double precision, allocatable :: s2(:,:)
 double precision, allocatable, target :: cr2(:,:),crx(:,:),cry(:,:),dum__(:,:)
 double precision, allocatable :: nrms(:),crs(:,:),ys(:,:),yz(:,:),cr1(:,:)
 double precision, allocatable, target :: cr(:,:),s(:),RR(:,:),u(:,:),v(:,:),crz(:,:)
 double precision, pointer :: crp(:,:),up(:,:),sp(:),vp(:,:),RRp(:,:)

    nswp = 20;         if(present(nswp_))     nswp = nswp_
    kickrank  = 4;     if(present(kickrank_)) kickrank= kickrank_
    kickrank2 = 0
    verb = 0;          if(present(verb_))     verb = verb_
    init_qr = .TRUE.;  if(present(init_qr_))  init_qr = init_qr_
    renorm = 1;        if(present(renorm_))   renorm = renorm_! 1:'direct', 2:'gram'
    fkick = .FALSE.;   if(present(fkick_))    fkick = fkick_

    allocate(z(0),dum__(0,0))

    !! else : x is cell3d_array (always)
    d = x%m - x%l + 1
    allocate(m(d), rx(d+1))
    m(1:d) = x%n(1:d)
    rx(1:d+1) = x%r(0:d)
    allocate(n(d), ra(d+1))
    n(1:d)= A% q(1:d)
    ra(1:d+1)= A% r(0:d)

    !! A is always a cell4d
    allocate(A2(d))
    do i=1,d
       A2(i)%arr = reshape(A% u4(i)% p, [ra(i)*n(i), m(i)*ra(i+1)])
    enddo

    if (.not.present(y0_)) then
       allocate(r2(d+1)); r2=2; r2(1)=1; r2(d+1)=1
       y = dtt_random_ortho(n,r_=r2)
    else
       y = y0_
    endif



    allocate(y2(d)) ! needed for reshaping y
    allocate(ry(d+1))
    ry(1:d+1) = y% r(0:d)



    if (kickrank+kickrank2.gt.0) then
       if (size(z).eq.0) then
          allocate(rz(d+1))
          rz = kickrank + kickrank2
          rz(1) = 1; rz(d+1) = 1
          z = core2cell(dtt_random_ortho(n,r_=rz))
       else
          allocate(rz(d+1))
       endif
       rz(1)= 1
       do i=1,d
          rz(i+1) = size(z(i)% arr, 3)
       enddo

       allocate(phizax(d+1))
       !! assume atype to be always 1
       allocate(phizax(1)% arr(1,ra(1),1))
       allocate(phizax(d+1)% arr(1,ra(d+1),1))
       phizax(1)% arr= 1d0; phizax(d+1)% arr= 1d0
       allocate(phizy(d+1))
       allocate(phizy(1)%arr(1,1,1))
       allocate(phizy(d+1)% arr(1,1,1))
       phizy(1)%arr = 1d0; phizy(d+1)%arr=1d0
    endif ! kickrank_+ kickrank2 > 0

    allocate(phiyax(d+1))
    !! atype is always 1
    allocate(phiyax(1)% arr(1,ra(1),1))
    allocate(phiyax(d+1)% arr(1,ra(1),1))
    phiyax(1)% arr= 1d0; phiyax(d+1)% arr= 1d0


    allocate(nrms(d)); nrms = 1

    !% Initial ort
    do i=1,d-1
       if (init_qr) then
          cr= reshape(y%u(i)%p,[ry(i)*n(i),ry(i+1)])
          if (renorm.eq.2.and.ry(i)*n(i).gt.5*ry(i+1)) then
             call svdgram(crp,sp,RRp, cr)
             cr = crp; s = sp; RR = RRp
          else
             call qr(cr1, RR, cr)
             cr = cr1
          endif
          nrmr = norm(RR)
          if (nrmr.gt.0) then
             RR = RR/nrmr
          endif
          cr2 = reshape(y%u(i+1)%p,[ry(i+1),n(i+1)*ry(i+2)])
          cr2 = thor_matmul(RR, cr2)
          ry(i+1) = size(cr, 2)
          deallocate(y%u(i)%p, y%u(i+1)%p)
          allocate(y%u(i)%p(ry(i),n(i),ry(i+1)))
          allocate(y%u(i+1)%p(ry(i+1),n(i+1),ry(i+2)))
          y%r(i-1:i+1)= ry(i:i+2)
          y%u(i)% p= reshape(cr, [ry(i),n(i),ry(i+1)])
          y%u(i+1)% p= reshape(cr2, [ry(i+1),n(i+1),ry(i+2)])
       endif

       call compute_next_Phi(phiyax(i+1)%arr, nrms(i), phiyax(i)%arr, y%u(i)%p, A2(i)%arr, x%u(i)%p, 'lr')

       if (kickrank + kickrank2 .gt. 0) then
          cr = reshape(z(i)% arr, [rz(i)*n(i), rz(i+1)])
          if (renorm.eq.2.and.rz(i)*n(i).gt.5*rz(i+1)) then
             call svdgram(crp,sp,RRp, cr)
             cr = crp; s = sp; RR = RRp
          else
             call qr(cr1, RR, cr)
             cr = cr1
          endif
          nrmr = norm(RR)
          if (nrmr.gt.0) then
              RR = RR/nrmr
          endif
          cr2 = reshape(z(i+1)%arr, [rz(i+1),n(i+1)*rz(i+2)])
          cr2 = thor_matmul(RR, cr2)
          rz(i+1) = size(cr, 2)
          z(i)%arr = reshape(cr, [rz(i),n(i),rz(i+1)])
          z(i+1)%arr = reshape(cr2,[rz(i+1),n(i+1),rz(i+2)])

          call compute_next_Phi(phizax(i+1)%arr, nrm, phizax(i)%arr, z(i)%arr, A2(i)%arr, x%u(i)%p, 'lr', extnrm=nrms(i))
          call compute_next_Phi(phizy(i+1)%arr, nrm, phizy(i)%arr, z(i)%arr, dum__, y%u(i)%p, 'lr')
       endif
    enddo

    i = d
    dir = -1
    swp = 1
    max_dx = 0

    swaploop: do while (swp.le.nswp)
       !% Project the MatVec generating vector
       crx = reshape(x%u(i)%p, [rx(i)*m(i)*rx(i+1), 1])
       cry = bfun3(phiyax(i)%arr,A2(i)%arr,phiyax(i+1)%arr,crx)
       nrms(i) = norm(cry)
       !% The main goal is to keep y{i} of norm 1
       if (nrms(i).gt.0d0) then
          cry = cry/nrms(i)
       else
          nrms(i) = 1d0
       endif
       y2(i)%arr = reshape(y%u(i)%p, [ry(i)*n(i)*ry(i+1),1])
       dx = sqrt(sum((cry - y2(i)%arr)**2))
       max_dx = max(max_dx, dx)

       !% Truncation and enrichment
       if (dir.gt.0 .and. i.lt.d) then
          cry = reshape(cry, [ry(i)*n(i), ry(i+1)])
          if (renorm .eq. 2) then
              call svdgram(up,sp,vp,cry,tol/sqrt(dble(d)))
              v = transpose(vp)
              r = size(up,2)
              u = up; s = sp
          else
              call svd(cry,up,sp,vp)
              s = sp
              r = my_chop2(s, tol*norm(s)/sqrt(dble(d)))
              u = up(:,1:r)
              !!!TODO v = thor_matmul(conjg(v(:,1:r)),diag(s(1:r)))
              v = thor_matmul(vp(1:r,:),diag(s(1:r)),'tn')
          endif
          up=> null(); sp=> null(); vp=> null()

          !% Prepare enrichment, if needed
          if (kickrank + kickrank2 .gt. 0) then
             cry = thor_matmul(u, v, tsp_='nt')
             cry = reshape(cry, [ry(i)*n(i), ry(i+1)])
             !% For updating z
             crz = bfun3(phizax(i)%arr, A2(i)%arr, phizax(i+1)%arr, crx)
             crz = reshape(crz, [rz(i)*n(i), rz(i+1)])
             ys = thor_matmul(cry, phizy(i+1)%arr(:,:,1))
             yz = reshape(ys, [ry(i), n(i)*rz(i+1)])
             yz = thor_matmul(phizy(i)%arr(:,:,1),yz)
             yz = reshape(yz, [rz(i)*n(i), rz(i+1)])
             crz = crz/nrms(i) - yz
             nrmz = norm(crz)
             if (kickrank2 .gt. 0d0) then
                call svd(crz,up,sp,vp)
                crz = up(:, 1:min(size(crz,2), kickrank))
                crz = stack(crz, rndmat(rz(i)*n(i), kickrank2))
                up=> null(); sp=> null(); vp=> null()
             end if
             !% For adding into solution
             if (fkick) then
                crs = bfun3(phiyax(i)%arr, A2(i)%arr, phizax(i+1)%arr, crx)
                crs = reshape(crs, [ry(i)*n(i), rz(i+1)])
                crs = crs/nrms(i) - ys
                u = stack(u, crs)

                if (renorm.eq.2.and.ry(i)*n(i).gt.5*(ry(i+1)+rz(i+1))) then
                   call svdgram(up, sp, RRp, u)
                   u = up; s = sp; RR = RRp
                else
                   call qr(cr1, RR, u)
                   u = cr1
                endif
                v = stack(v, zeros(ry(i+1), rz(i+1)))
                v = thor_matmul(v, RR, tsp_='nt')
                r = size(u, 2)
             endif
          endif
          deallocate(y%u(i)%p)
          allocate(y%u(i)%p(ry(i), n(i), r))
          y%u(i)%p = reshape(u, [ry(i), n(i), r])
          y%r(i-1) = ry(i); y%r(i) = r

          cr2 = reshape(y%u(i+1)%p,[ry(i+1),n(i+1)*ry(i+2)])
          v = reshape(v, [ry(i+1), r])
          cr2 = thor_matmul(v, cr2, tsp_='tn')
          deallocate(y%u(i+1)%p)
          allocate(y%u(i+1)%p(r, n(i+1), ry(i+2)))
          y%u(i+1)%p = reshape(cr2, [r, n(i+1), ry(i+2)])
          y%r(i) = r; y%r(i+1) = ry(i+2)

          ry(i+1) = r


          call compute_next_Phi(phiyax(i+1)%arr, nrms(i), phiyax(i)%arr, y%u(i)%p, A2(i)%arr, x%u(i)%p, 'lr')

          if (kickrank + kickrank2 .gt. 0) then
             if (renorm.eq.2.and.rz(i)*n(i).gt.5*rz(i+1)) then
                call svdgram(crp, sp, RRp, crz)
                crz = crp; s = sp; RR = RRp
             else
                call qr(cr1, RR, crz)
                crz = cr1
             endif
             rz(i+1) = size(crz, 2)
             z(i)%arr = reshape(crz, [rz(i),n(i),rz(i+1)])
             !% z{i+1} will be recomputed from scratch in the next step
             call compute_next_Phi(phizax(i+1)%arr, nrm, phizax(i)%arr, z(i)%arr, A2(i)%arr, x%u(i)%p, 'lr', extnrm=nrms(i))
             call compute_next_Phi(phizy(i+1)%arr, nrm, phizy(i)%arr, z(i)%arr, dum__, y%u(i)%p, 'lr')
          endif
       else if (dir.lt.0.and.i.gt.1) then
          cry = reshape(cry,[ry(i), n(i)*ry(i+1)])
          if (renorm.eq.2) then
             call svdgram(vp,sp,up,transpose(cry),tol/sqrt(dble(d)))
             u = transpose(up)
             r = size(vp, 2)
             s = sp; v = vp
          else
             call svd(cry,up,sp,vp)
             s = sp
             r = my_chop2(s, tol*norm(s)/sqrt(dble(d)))
             !!!TODO v = conjg(v(:,1:r))
             v = transpose(vp(1:r,:))
             u = thor_matmul(up(:,1:r),diag(s(1:r)))
          endif
          up=> null(); sp=> null(); vp=> null()

          !% Prepare enrichment, if needed
          if (kickrank + kickrank2 .gt. 0) then
             cry = thor_matmul(u, v, tsp_='nt')
             cry = reshape(cry, [ry(i), n(i)*ry(i+1)])
             !% For updating z

             crz = bfun3(phizax(i)%arr, A2(i)%arr, phizax(i+1)%arr, crx)
             crz = reshape(crz, [rz(i), n(i)*rz(i+1)])
             ys = thor_matmul(phizy(i)%arr(:,:,1), cry)
             yz = reshape(ys, [rz(i)*n(i), ry(i+1)])
             yz = thor_matmul(yz, phizy(i+1)%arr(:,:,1))
             yz = reshape(yz, [rz(i), n(i)*rz(i+1)])
             crz = crz/nrms(i) - yz
             nrmz = norm(crz)
             if (kickrank2.gt.0) then
                call svd(crz,up,sp,vp)
                crz = transpose(vp(:,1:min(size(crz,2), kickrank)))
                crz = stack(crz, rndmat(kickrank2, n(i)*rz(i+1)), 1)
             endif
             !% For adding into solution
             crs = bfun3(phizax(i)%arr, A2(i)%arr, phiyax(i+1)%arr, crx)
             crs = reshape(crs, [rz(i), n(i)*ry(i+1)])
             crs = crs/nrms(i) - ys
             v = stack(v, transpose(crs))

             if (renorm.eq.2.and.n(i)*ry(i+1)>5*(ry(i)+rz(i))) then
                call svdgram(vp,sp,RRp,v)
                v = vp; s = sp; RR = RRp
             else
                call qr(cr1,RR,v)
                v = cr1
             endif
             u = stack(u, zeros(ry(i), rz(i)))
             u = thor_matmul(u, RR, tsp_='nt')
             r = size(v, 2)
          endif
          cr2 = reshape(y%u(i-1)%p, [ry(i-1)*n(i-1), ry(i)])
          cr2 = thor_matmul(cr2, u)
          deallocate(y%u(i-1)%p, y%u(i)%p)
          allocate(y%u(i-1)%p(ry(i-1), n(i-1), r))
          allocate(y%u(i)%p(r,n(i),ry(i+1)))
          y%u(i-1)%p = reshape(cr2, [ry(i-1), n(i-1), r])
          y%u(i)%p = reshape(transpose(v), [r,n(i),ry(i+1)])
          y%r(i-2) = ry(i-1); y%r(i-1) = r; y%r(i) = ry(i+1)

          ry(i) = r


          call compute_next_Phi(phiyax(i)%arr, nrms(i), phiyax(i+1)%arr, y%u(i)%p, A2(i)%arr, x%u(i)%p, 'rl')

          if (kickrank + kickrank2 .gt. 0) then

              if (renorm.eq.2 .and. n(i)*rz(i+1).gt.5*rz(i)) then
                  call svdgram(crp,sp,RRp,transpose(crz))
                  crz = crp; s = sp; RR = RRp
              else
                  call qr(cr1,RR,transpose(crz))
                  crz = cr1
              endif
              rz(i) = size(crz, 2)

              z(i)%arr = reshape(transpose(crz), [rz(i), n(i), rz(i+1)])
              !% don't update z{i-1}, it will be recomputed from scratch
              call compute_next_Phi(phizax(i)%arr, nrm, phizax(i+1)%arr, z(i)%arr, A2(i)%arr, x%u(i)%p, 'rl', extnrm=nrms(i))
              call compute_next_Phi(phizy(i)%arr, nrm, phizy(i+1)%arr, z(i)%arr, dum__, y%u(i)%p, 'rl')
          endif
       endif

       if (verb.gt.1) then

          print '("amen-mv: swp=["I4","I4"], dx="ES10.3", r="I4", |y|="ES10.3", |z|="ES10.3)', &
                swp, i, dx, r, norm(cry), nrmz
       endif

       !% Stopping or reversing

       if ((dir.gt.0.and.i.eq.d).or.(dir.lt.0.and.i.eq.1)) then
          if (verb.gt.0) then
             print '("amen-mv: swp="I4"{"I4"}, max_dx="ES10.3", max_r="I4)', &
                   swp, (1-dir)/2, max_dx, maxval(ry)
          endif
          if ((max_dx.lt.tol .or. swp.eq.nswp).and.(dir.gt.0)) then
             exit swaploop
          else
             ! We are at the terminal block
             deallocate(y%u(i)%p)
             allocate(y%u(i)%p(ry(i),n(i),ry(i+1)))
             y%u(i)%p= reshape(cry, [ry(i),n(i),ry(i+1)])
             y%r(i-1)= ry(i); y%r(i)= ry(i+1)
             if (dir.gt.0) swp = swp + 1
          endif
          max_dx = 0d0
          dir = -dir
       else
          i = i + dir
       endif
    enddo swaploop

    !% Distribute norms equally...
    nrms = dexp(sum(dlog(nrms))/dble(d))
    !% ... and plug them into y
    do i=1,d
       y%u(i)%p = y%u(i)%p * nrms(i)
    enddo

    ! cleanup
    if (allocated(s2)) deallocate(s2)
    if (allocated(cr1)) deallocate(cr1)
    if (allocated(cr2)) deallocate(cr2)
    if (allocated(crx)) deallocate(crx)
    if (allocated(cry)) deallocate(cry)
    if (allocated(dum__)) deallocate(dum__)
    if (allocated(nrms)) deallocate(nrms)
    if (allocated(crs)) deallocate(crs)
    if (allocated(ys)) deallocate(ys)
    if (allocated(yz)) deallocate(yz)
    if (allocated(cr)) deallocate(cr)
    if (allocated(s)) deallocate(s)
    if (allocated(RR)) deallocate(RR)
    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)
    if (allocated(crz)) deallocate(crz)

    if (allocated(z)) then
       do i=1,size(z); deallocate(z(i)% arr); enddo
       deallocate(z)
    endif

 end subroutine amen_mv


 !> Performs the recurrent Phi (or Psi) matrix computation
 !! Phi = Phi_prev * (x'Ay).
 !! If direction is 'lr', computes Psi
 !! if direction is 'rl', computes Phi
 !! A can be empty, then only x'y is computed.
 !!  Phi1: rx1, ry1, ra1, or {rx1, ry1}_ra, or rx1, ry1
 !!  Phi2: ry2, ra2, rx2, or {ry2, rx2}_ra, or ry2, rx2
 subroutine compute_next_Phi(Phi, nrm, Phi_prev, x, A, y, direction, extnrm)
 use mat_lib, only: normfro, thor_matmul
 use time_lib
 implicit none
 double precision, intent(inout), allocatable :: Phi(:,:,:)
 double precision, intent(inout) :: nrm
 double precision, intent(in)  :: Phi_prev(:,:,:)
 double precision, intent(in)  :: x(:,:,:)
 double precision, intent(in)  :: A(:,:)
 double precision, intent(in)  :: y(:,:,:)
 character(len=2), intent(in)  :: direction
 double precision, intent(in), optional :: extnrm

 double precision, allocatable :: x2(:,:)   ! 2D-reshaped version of x
 double precision, allocatable :: y2(:,:)   ! 2D-reshaped version of y
 double precision, allocatable :: Phi2(:,:) ! 2D-reshaped version of Phi_prev
 integer :: rx1,ry1,ra1,rx2,ry2,ra2,m,n
 logical :: isempty_A

    !
    !!! TODO: WIP
    !
    !! if (nargin<6)
    !!     extnrm = [];
    !! end
    !!
    !!
    rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
    ry1 = size(y,1); m = size(y,2); ry2 = size(y,3);
    ra1 = size(A,1)/n; ra2 = size(A,2)/m;
    isempty_A = ra1.eq.0 .and. ra2.le.0
    if (isempty_A) then
       ra1 = 1; ra2 = 1
    endif

    !!                                                if (isa(Phi_prev, 'cell'))
    !!                                                 . . . < skipped >
    !!                                                else (Phi_prev is not a cell)
    if (direction.eq.'lr') then
       !%lr: Phi1
       x2 = reshape(x, [rx1, n*rx2])              !! x = reshape(x, rx1, n*rx2);
       Phi2 = reshape(Phi_prev, [rx1,ry1*ra1])    !! Phi = reshape(Phi_prev, rx1, ry1*ra1);
       Phi2 = thor_matmul(x2, Phi2, tsp_='tn')    !! Phi = x'*Phi;
       if (.not.isempty_A) then                   !! if (~isempty(A))
          Phi2 = reshape(Phi2, [n*rx2*ry1, ra1])  !!     Phi = reshape(Phi, n*rx2*ry1, ra1);
          Phi2 = transpose(Phi2)                  !!     Phi = Phi.';
          Phi2 = reshape(Phi2, [ra1*n, rx2*ry1])  !!     Phi = reshape(Phi, ra1*n, rx2*ry1);
          Phi2 = thor_matmul(A, Phi2, tsp_='tn')  !!     Phi = A.'*Phi;
          Phi2 = reshape(Phi2, [m, ra2*rx2*ry1])  !!     Phi = reshape(Phi, m, ra2*rx2*ry1);
       else                                       !! else
          Phi2 = reshape(Phi2, [n, rx2*ry1])      !!     Phi = reshape(Phi, n, rx2*ry1);
       endif                                      !! end;
       Phi2 = transpose(Phi2)                     !! Phi = Phi.';
       Phi2 = reshape(Phi2, [ra2*rx2, ry1*m])     !! Phi = reshape(Phi, ra2*rx2, ry1*m);
       y2 = reshape(y, [ry1*m, ry2])              !! y = reshape(y, ry1*m, ry2);
       Phi2 = thor_matmul(Phi2, y2)               !! Phi = Phi*y;
       if (.not.isempty_A) then                   !! if (~isempty(A))
          Phi2 = reshape(Phi2, [ra2, rx2*ry2])    !!     Phi = reshape(Phi, ra2, rx2*ry2
          Phi2 = transpose(Phi2)                  !!     Phi = Phi.';
       endif                                      !! end;
       Phi = reshape(Phi2, [rx2, ry2, ra2])       !! Phi = reshape(Phi, rx2, ry2, ra2);
    else
       !! %rl: Phi2
       y2 = reshape(y, [ry1*m, ry2])              !! y = reshape(y, ry1*m, ry2);
       Phi2 = reshape(Phi_prev, [ry2,ra2*rx2])    !! Phi = reshape(Phi_prev, ry2, ra2*rx2);
       Phi2 = thor_matmul(y2, Phi2)               !! Phi = y*Phi;
       if(.not.isempty_A) then                    !! if (~isempty(A))
          Phi2 = reshape(Phi2, [ry1, m*ra2*rx2])  !!     Phi = reshape(Phi, ry1, m*ra2*rx2);
          Phi2 = transpose(Phi2)                  !!     Phi = Phi.';
          Phi2 = reshape(Phi2, [m*ra2, rx2*ry1])  !!     Phi = reshape(Phi, m*ra2, rx2*ry1);
          Phi2 = thor_matmul(A, Phi2)             !!     Phi = A*Phi;
          Phi2 = reshape(Phi2, [ra1*n*rx2, ry1])  !!     Phi = reshape(Phi, ra1*n*rx2, ry1);
          Phi2 = transpose(Phi2)                  !!     Phi = Phi.';
       endif                                      !! end
       Phi2 = reshape(Phi2, [ry1*ra1, n*rx2])     !! Phi = reshape(Phi, ry1*ra1, n*rx2);
       x2 = reshape(x, [rx1, n*rx2])              !! x = reshape(x, rx1, n*rx2);
       Phi2 = thor_matmul(Phi2, x2, tsp_='nt')    !! Phi = Phi*x';
       if (.not.isempty_A) then                   !! if (~isempty(A))
          Phi = reshape(Phi2, [ry1, ra1, rx1])    !!     Phi = reshape(Phi, ry1, ra1, rx1);
       else                                       !! else
          Phi = reshape(Phi2, [ry1, rx1, 1])      !!     Phi = reshape(Phi, ry1, rx1);
       endif                                      !! end
    endif

    if (present(extnrm)) then                  !! if (~isempty(extnrm))
       !% Override the normalization by the external one
       Phi = Phi/extnrm                        !!     Phi = Phi/extnrm;
    else                                       !! elseif (nargout>1)
       ! Extract the scale to prevent overload
       nrm = normfro(Phi)                      !!     nrm = norm(Phi(:), 'fro');
       if (nrm.gt.0) then                      !!     if (nrm>0)
           Phi = Phi/nrm                       !!         Phi = Phi/nrm;
       else                                    !!     else
           nrm = 1d0                           !!         nrm=1;
       endif                                   !!     end
    endif                                      !! end
    !!                                            end %  if (isa(Phi_prev, 'cell'))
    ! cleanup
    if (allocated(x2)) deallocate(x2)
    if (allocated(y2)) deallocate(y2)
    if (allocated(Phi2)) deallocate(Phi2)
 end subroutine compute_next_Phi


 !> function [X]=amen_mm(A, Y, tol, varargin)
 !! Approximate the matrix-by-matrix via the AMEn iteration
 !!    [X]=amen_mv(A, Y, tol, varargin)
 !!    Attempts to approximate the X = A*Y with accuracy TOL using the
 !!    AMEn+ALS iteration. A is a n x m matrix, Y is a m x k matrix.
 !!    Matrices A,Y can be given in tt_matrix format, the output is tt_matrix.
 !!    Y can be tt_tensor, it's considered as column, the output is tt_tensor.
 !!    A and Y can be a {d,R} cell array. However, X is always a "single" TT
 !!    (no tensor chain), since that's how ALS works. Generally, X has the
 !!    same form as Y, except that it's always {d,1} in essense. X and Y can't
 !!    be sparse (SVD will destroy it anyways), but A can be.
 !!
 !!    Options are provided in form
 !!    'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
 !!    on. The parameters are set to default (in brackets in the following)
 !!    The list of option names and default values are:
 !!        o x0 - initial approximation to AY [rand with ranks of Y(:,1)]
 !!        o nswp - maximal number of sweeps [20]
 !!        o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
 !!        o kickrank - compression rank of the error,
 !!                     i.e. enrichment size [4]
 !!
 !!
 !! ********
 !!    For description of adaptive ALS please see
 !!    Sergey V. Dolgov, Dmitry V. Savostyanov,
 !!    Alternating minimal energy methods for linear systems in higher dimensions.
 !!    Part I: SPD systems, http://arxiv.org/abs/1301.6068,
 !!    Part II: Faster algorithm and application to nonsymmetric systems, http://arxiv.org/abs/1304.1222
 !!
 !! Call tree:
 !! +- amen_mm
 !!    +-- leftreduce_matrix
 !!    +-- leftreduce_vector
 !!    +-- local_matvec
 !!    +-- rightreduce_matrix
 !!    +-- rightreduce_vector

 subroutine amen_mm(X, A, Y, tol, x0_, nswp_, verb_, kickrank_, init_qr_)
 use matlab_struct_module, only: array2d, array3d, array4d, rand_cell3d, pprint_matrix
 use mat_lib, only: svd, svdgram, bfun3, my_chop2, stack, &
                    diag => matlab_diag, zeros => matlab_zeros, &
                    qr => matlab_qr, norm => normfro
 use rnd_lib, only: randn
 use thor_lib, only: dtt_tensor, dtt_matrix, core2cell, dttm_random_ortho
 use mat_lib, only: thor_matmul, thor_matmul_mnk
 use thorio_lib, only: dttm_read_ascii_file, dttm_write_ascii_file
 use time_lib
 implicit none
 class(dtt_matrix),intent(INOUT) :: X
 class(dtt_matrix),intent(IN) :: A, Y
 double precision, intent(IN) :: tol
 class(dtt_matrix),intent(IN), optional:: x0_ !< initial approximation to A@Y
 integer, intent(IN),optional :: nswp_        !< maximal number of sweeps [20]
 integer, intent(IN),optional :: verb_        !< verbosity level, 0-silent, 1-sweep info, 2-block info [1]
 integer, intent(IN),optional :: kickrank_    !< compression rank of the error,
                                              !  i.e. enrichment size [3]
 logical, intent(IN),optional :: init_qr_     !< perform QR of the input (save some time in ts, etc) [true]
 !
 type(array3d), allocatable:: z(:) !< z0 - initial approximation to the error Ax-y [rand rank-kickrank]
 integer, allocatable :: m(:), k(:), rx(:),n(:),ra(:),ry(:),rz(:),r2(:)
 integer :: d,dir,i,j,swp,r,Rry,Rra, sh(2)
 integer :: nswp,verb,kickrank,renorm
 logical :: init_qr,fkick
 integer :: kickrank2, Ntim
 type(array3d), allocatable :: phiax(:), phiyax(:), phizax(:), phizy(:)
 type(array2d), allocatable :: y2(:), A2(:)
 double precision :: dx,max_dx,nrm,nrm_s,nrmr,nrmz, aux
 double precision, allocatable :: s2(:,:), cr4(:,:,:,:)
 double precision, allocatable, target :: cr2(:,:),crx(:,:),cry(:,:),dum__(:,:)
 double precision, allocatable :: nrms(:),crs(:,:),ys(:,:),yz(:,:),cr1(:,:)
 double precision, allocatable, target :: cr(:,:),s(:),RR(:,:),u(:,:),v(:,:),crz(:,:)
 double precision, pointer   :: crp(:,:),up(:,:),sp(:),vp(:,:),RRp(:,:)
 double precision, parameter :: one2(1,1) = reshape([1d0], [1,1])
 double precision, allocatable :: XAY2(:,:)
 type(array2d), allocatable  :: XAY(:), ZAY(:), ZX(:)

 nswp= 20;       if(present(nswp_)) nswp= nswp_             ! nswp = 20;
 kickrank= 4;    if(present(kickrank_)) kickrank= kickrank_ ! kickrank = 4;
 verb= 1;        if(present(verb_)) verb= verb_             ! verb = 1;
 init_qr=.true.; if(present(init_qr_)) init_qr= init_qr_    ! init_qr = true;

 d= Y% m - Y% l + 1
 allocate(m(d),n(d),k(d),ry(d+1))
 m= Y% q(1:d)
 k= Y% s(1:d)
 ry=Y% r(0:d)
 Rry= 1

 !% Grumble matrix
 if (A%m-A%l+1.ne.d) error stop "[!][amen_mm]: mismatched orders in A and Y"
 allocate(ra(d+1)); ra(1:d+1)= A% r(0:d)
 if (any(A%s(1:d).ne.m(1:d))) &
    error stop "[!][amen_mm]: Mismatching TT dimensionalities in A and Y"
 n= A% q(1:d)

 allocate(A2(d))
 do i=1,d
    A2(i)%arr = reshape(A% u4(i)% p, [ra(i)*n(i), m(i)*ra(i+1)])
 enddo

 if (.not.present(x0_)) then
    X = dttm_random_ortho(n(1:d), k(1:d), r_=ry(1:d+1))
    allocate(rx(d+1)); rx(1:d+1)= X% r(0:d)
    init_qr= .false.
 else
    X = x0_
    if (X%m-X%l+1 /= d .or. any(X%q(1:d)/=n) .or. any(X%s(1:d)/=k)) &
       error stop "[!][amen_mm]: mismatched dimensions in X & Y"
    allocate(rx(d+1)); rx(1:d+1)= X% r(0:d)
 endif

 !% Reductions
 allocate(XAY(d+1))
 allocate(ZAY(d+1))
 allocate(ZX(d+1))
 ZX(1)%arr = one2; ZX(d+1)%arr = one2
 XAY(1)%arr= one2; XAY(d+1)%arr= one2
 ZAY(1)%arr= one2; ZAY(d+1)%arr= one2

 !% Residual rank
 allocate(rz(d+1)); rz= kickrank; rz(1)=1; rz(d+1)=1
 allocate(nrms(d)); nrms= 1d0

 !% Initial ort
 do i=1,d-1
    if (init_qr) then
       cr= reshape(X%u4(i)%p, [rx(i)*n(i)*k(i), rx(i+1)])
       call qr(cr1, RR, cr)
       cr = cr1
       nrmr= norm(RR)
       if (nrmr > 0d0) then
          RR= RR/nrmr
       endif
       cr2= thor_matmul_mnk(rx(i+1), rx(i+1), n(i+1)*k(i+1)*rx(i+2), &
                            RR, X%u4(i+1)%p)
       rx(i+1)= size(cr, 2)
       deallocate(X% u4(i)% p)
       allocate(X% u4(i)% p(rx(i), n(i), k(i), rx(i+1)))
       X% u4(i)% p= reshape(cr, [rx(i), n(i), k(i), rx(i+1)])
       deallocate(X% u4(i+1)% p)
       allocate(X% u4(i+1)% p(rx(i+1), n(i+1), k(i+1), rx(i+2)))
       X% u4(i+1)% p= reshape(cr2, [rx(i+1), n(i+1), k(i+1), rx(i+2)])
    endif

    !% Reduce
    XAY(i+1)%arr= leftreduce_matrix(nrms(i), XAY(i)%arr, X% u4(i)% p, &
                           A2(i)%arr, Y% u4(i)% p, rx(i),n(i),k(i),rx(i+1),      &
                           ra(i),ra(i+1),ry(i),m(i),ry(i+1))

    !% Residual reductions
    if (kickrank > 0) then
       cr= randn(rz(i)*n(i)*k(i), rz(i+1))
       call qr(cr1, RR, cr)
       rz(i+1)= size(cr1, 2)
       ZAY(i+1)%arr= leftreduce_matrix(aux, ZAY(i)% arr, &
                           reshape(cr1,[size(cr1),1,1,1]), &
                           A2(i)%arr, Y% u4(i)% p, rz(i),n(i),k(i),rz(i+1), &
                           ra(i),ra(i+1), ry(i),m(i),ry(i+1), extnrm_=nrms(i))
       ZX(i+1)%arr= leftreduce_vector(ZX(i)%arr, cr1, X%u4(i)%p, &
                       rz(i),n(i),k(i),rz(i+1), rx(i),rx(i+1))
    endif
 enddo

 !% Normalize the last block for small error on good initial guess
 nrmr= norm(X% u4(d)% p)
 if (nrmr > 0d0) then
    X% u4(d)% p= X% u4(d)% p/nrmr
 endif

 i = d
 dir = -1
 swp = 1
 max_dx = 0d0

 !% Iteration
 sweeps: do while(swp <= nswp)
    !% Project the MatVec generating vector
    cr= local_matvec(Y% u4(i)% p, ry(i), m(i), k(i), ry(i+1), &
            rx(i),n(i),rx(i+1), XAY(i)%arr, A2(i)%arr, XAY(i+1)%arr, &
            ra(i),ra(i+1))

    nrms(i)= norm(cr)
    !% The main goal is to keep y{i} of norm 1
    if (nrms(i) > 0d0) then
       cr= cr/nrms(i)
    else
       nrms(i)= 1d0
    endif

    dx= norm(cr - reshape(X% u4(i)% p, [rx(i)*n(i)*k(i)*rx(i+1), 1]))
    max_dx= max(dx, max_dx)

    !% Truncation and enrichment
    if (dir > 0 .and. i < d) then
       cr= reshape(cr, [rx(i)*n(i)*k(i), rx(i+1)])
       call svd(cr,up,sp,vp)
       s= sp
       r= my_chop2(s, tol*norm(s)/sqrt(dble(d)))
       u= up(:,1:r)
       v= thor_matmul(vp(1:r,:),diag(s(1:r)),'tn')

       !% Prepare enrichment, if needed
       if (kickrank > 0) then
          cr= thor_matmul(u, v, tsp_='nt')
          cr= reshape(cr, [rx(i)*n(i)*k(i), rx(i+1)])
          !% For updating z
          crz= local_matvec(Y% u4(i)% p, ry(i),m(i),k(i),ry(i+1), &
                   rz(i),n(i),rz(i+1), ZAY(i)%arr, A2(i)%arr, &
                   ZAY(i+1)%arr, ra(i),ra(i+1))
          crz= reshape(crz, [rz(i)*n(i)*k(i), rz(i+1)])
          ys= thor_matmul(cr, ZX(i+1)%arr)
          yz= reshape(ys, [rx(i), n(i)*k(i)*rz(i+1)])
          yz= thor_matmul(ZX(i)%arr, yz)
          yz= reshape(yz, [rz(i)*n(i)*k(i), rz(i+1)])
          crz= crz/nrms(i) - yz
          nrmz= norm(crz)
          !% For adding into solution
          crs= local_matvec(Y% u4(i)% p, ry(i),m(i),k(i),ry(i+1), &
                   rx(i),n(i),rz(i+1), XAY(i)%arr, A2(i)%arr, &
                   ZAY(i+1)%arr, ra(i),ra(i+1))
          crs= reshape(crs, [rx(i)*n(i)*k(i), rz(i+1)])
          crs= crs/nrms(i) - ys
          u= stack(u, crs)
          call qr(cr1, RR, u)
          u= cr1
          v= stack(v, zeros(rx(i+1), rz(i+1)))
          v= thor_matmul(v, RR, tsp_='nt')
          r= size(u, 2)
       endif
       deallocate(X% u4(i)% p)
       allocate(X% u4(i)% p(rx(i),n(i),k(i),r))
       call dcopy(rx(i)*n(i)*k(i)*r, u, 1, X% u4(i)% p, 1)
       cr2= reshape(X% u4(i+1)% p, [rx(i+1), n(i+1)*k(i+1)*rx(i+2)])
       v= reshape(v, [rx(i+1),r])
       cr2= thor_matmul(v, cr2, tsp_='tn')
       deallocate(X% u4(i+1)% p)
       allocate(X% u4(i+1)% p(r,n(i+1),k(i+1),rx(i+2)))
       call dcopy(r*n(i+1)*k(i+1)*rx(i+2), cr2, 1, X% u4(i+1)% p, 1)

       rx(i+1)= r

       nrms(i+dir)= nrms(i)
       !% Reduce
       XAY(i+1)%arr= leftreduce_matrix(nrms(i), XAY(i)%arr, X%u4(i)%p, A2(i)%arr, &
                     Y%u4(i)%p, rx(i),n(i),k(i),rx(i+1),ra(i),ra(i+1),ry(i),m(i),ry(i+1))
       !% Enrichment
       if (kickrank > 0) then
          call qr(cr1, RR, crz); crz= cr1
          rz(i+1)= size(crz, 2)
          ZAY(i+1)%arr= leftreduce_matrix(aux, ZAY(i)%arr, &
                        reshape(crz,[size(crz),1,1,1]), A2(i)%arr, &
                        Y%u4(i)%p, rz(i),n(i),k(i),rz(i+1), ra(i),ra(i+1), &
                        ry(i),m(i),ry(i+1), nrms(i))
          ZX(i+1)%arr= leftreduce_vector(ZX(i)%arr, crz, X%u4(i)%p, &
                       rz(i),n(i),k(i),rz(i+1), rx(i),rx(i+1))
       endif

    else if (dir < 0 .and. i > 1) then
       cr= reshape(cr, [rx(i), n(i)*k(i)*rx(i+1)])
       call svd(cr,up,sp,vp)
       s= sp
       r= my_chop2(s, tol*norm(s)/sqrt(dble(d)))
       v= transpose(vp(1:r,:))
       u= thor_matmul(up(:,1:r),diag(s(1:r)))

       !% Prepare enrichment, if needed
       if (kickrank > 0) then
          cr= thor_matmul(u, v, tsp_='nt')
          cr= reshape(cr, [rx(i), n(i)*k(i)*rx(i+1)])
          !% For updating z
          crz= local_matvec(Y%u4(i)%p, ry(i),m(i),k(i),ry(i+1), &
               rz(i),n(i),rz(i+1), ZAY(i)%arr, A2(i)%arr, ZAY(i+1)%arr, &
               ra(i),ra(i+1))
          crz= reshape(crz, [rz(i), n(i)*k(i)*rz(i+1)])
          ys= thor_matmul(ZX(i)%arr, cr)
          yz= reshape(ys, [rz(i)*n(i)*k(i), rx(i+1)])
          yz= thor_matmul(yz, ZX(i+1)%arr)
          yz= reshape(yz, [rz(i), n(i)*k(i)*rz(i+1)])
          crz= crz/nrms(i) - yz
          nrmz= norm(crz)
          !% For adding into solution
          crs= local_matvec(Y%u4(i)%p, ry(i),m(i),k(i),ry(i+1), &
               rz(i),n(i),rx(i+1), ZAY(i)%arr, A2(i)%arr, XAY(i+1)%arr, &
               ra(i),ra(i+1))
          crs= reshape(crs, [rz(i), n(i)*k(i)*rx(i+1)])
          crs= crs/nrms(i) - ys
          v= stack(v, transpose(crs))
          call qr(cr1, RR, v); v= cr1
          u= stack(u, zeros(rx(i), rz(i)))
          u= thor_matmul(u, RR, tsp_='nt')
          r= size(v, 2)
       endif
       cr2= reshape(X% u4(i-1)% p, [rx(i-1)*n(i-1)*k(i-1), rx(i)])
       cr2= thor_matmul(cr2, u)
       deallocate(X% u4(i-1)% p)
       allocate(X% u4(i-1)% p(rx(i-1),n(i-1),k(i-1),r))
       call dcopy(rx(i-1)*n(i-1)*k(i-1)*r, cr2, 1, X% u4(i-1)% p, 1)
       v= transpose(v)
       deallocate(X% u4(i)% p)
       allocate(X% u4(i)% p(r,n(i),k(i),rx(i+1)))
       call dcopy(r*n(i)*k(i)*rx(i+1), v, 1, X% u4(i)% p, 1)

       rx(i)= r

       nrms(i+dir) = nrms(i) !% this is the norm of my block, save it
       !% Reduce
       XAY(i)%arr= rightreduce_matrix(nrms(i), XAY(i+1)%arr, X%u4(i)%p, &
                   A2(i)%arr, Y%u4(i)%p, rx(i),n(i),k(i),rx(i+1), &
                   ra(i),ra(i+1), ry(i),m(i),ry(i+1))
       !% Enrich
       if (kickrank > 0) then
          call qr(cr1, RR, transpose(crz)); crz= cr1
          rz(i)= size(crz, 2)
          ZAY(i)%arr= rightreduce_matrix(aux, ZAY(i+1)%arr, &
                      reshape(crz, [size(crz),1,1,1]), A2(i)%arr, &
                      Y%u4(i)%p, rz(i),n(i),k(i),rz(i+1), ra(i),ra(i+1), &
                      ry(i),m(i),ry(i+1), nrms(i))
          ZX(i)%arr= rightreduce_vector(ZX(i+1)%arr, crz, X%u4(i)%p, &
                     rz(i),n(i),k(i),rz(i+1), rx(i),rx(i+1))
       endif
    else
          deallocate(X% u4(i)% p)
          allocate(X% u4(i)% p(rx(i),n(i),k(i),rx(i+1)))
          call dcopy(rx(i)*n(i)*k(i)*rx(i+1), cr, 1, X%u4(i)%p, 1)

    endif

    if (verb > 1) then
       print '("amen-mm: swp=["I3","I3"], dx="ES14.7", r="I5", '// &
             '|X|="ES14.7", |z|="ES14.7)', swp, i, dx, r, norm(cr), nrmz
    endif

    !% Stopping or reversing
    if ((dir > 0.and.i == d) .or. (dir < 0.and.i == 1)) then
       if (verb > 0) then
          print '("amen-mm: swp=",I3,"{",I2,"}, max_dx=",ES14.7,", max_r=",I5)', &
                swp, (1-dir)/2, max_dx, maxval(rx)
       endif
       if ((max_dx < tol .or. swp==nswp).and. dir > 0) then
          exit sweeps
       endif
       swp= swp + 1
       max_dx= 0d0
       dir= -dir
    else
       i= i + dir
    endif
 enddo sweeps

 !% Distribute norms equally...
 nrmz= dexp(sum(dlog(nrms))/dble(d))
 !% ... and plug them into y
 X% r(0:d)= rx(1:d+1)
 X% n(1:d)= X% q(1:d) * X% s(1:d)
 do i=1,d
    X% u4(i)% p= X% u4(i)% p * nrmz
    !!! TODO: make sure other fields of X are correctly filled out!!!
    X% u(i)% p(1:rx(i),1:X%n(i),1:rx(i+1))=> X% u4(i)% p
 enddo
 end subroutine amen_mm


 !% Accumulates the left reduction W{1:k}'*A{1:k}*X{1:k}
 function leftreduce_matrix(nrm, WAX1, w, A2, x, rw1, n, k, rw2, &
                            ra1, ra2, rx1, m, rx2, extnrm_) result(WAX2)
 use matrix_util, only: perm3d
 use mat_lib, only: norm=> normfro, thor_matmul_mnk
 use time_lib
 implicit none
 double precision, allocatable, target   :: WAX2(:,:), aux(:,:)
 double precision, intent(INOUT) :: nrm
 double precision, intent(IN)    :: WAX1(:,:), w(:,:,:,:), A2(:,:), x(:,:,:,:)
 integer, intent(IN)             :: rw1, n, k, rw2, m, ra1, ra2, rx1, rx2
 double precision, intent(IN), optional :: extnrm_
 !
 double precision, allocatable :: w2(:,:), w3(:,:,:), xc(:,:,:), xc2(:,:)
 double precision, pointer :: wp(:)

    !% Left WAX has the form of the first matrix TT block, i.e. [rw, rx, ra]
    nrm= 0d0
    w3 = reshape(w, [rw1*n, k, rw2])
    call perm3d(w3, 132)
    xc = reshape(x, [rx1*m, k, rx2])
    call perm3d(xc, 321)
    WAX2= thor_matmul_mnk(n*rw2*k,rx1*ra1,rw1,w3,WAX1,tsp_='tn')
    WAX2= transpose(reshape(WAX2, [n, rw2*k*rx1*ra1]))
    WAX2= thor_matmul_mnk(rw2*k*rx1, m*ra2, n*ra1, WAX2, A2)
    WAX2= transpose(reshape(WAX2, [rw2, k*rx1*m*ra2]))
    WAX2= thor_matmul_mnk(rx2, ra2*rw2, k*rx1*m, xc, WAX2)
    WAX2= transpose(reshape(WAX2, [rx2*ra2, rw2]))
    nrm= max(nrm, norm(WAX2))

    !% Extract the scale to prevent overload
    if (.not.present(extnrm_)) then
       if (nrm > 0d0) then
          WAX2= WAX2/nrm
       else
          nrm = 1d0
       endif
    else
       !% Override the normalization
       WAX2= WAX2/extnrm_
    endif
 end function leftreduce_matrix


 !% Accumulates the right reduction W{k:d}'*A{k:d}*X{k:d}
 function rightreduce_matrix(nrm, WAX2, w, A2, x, rw1, n, k, rw2, &
                             ra1, ra2, rx1, m, rx2, extnrm_) result(WAX1)
 use matrix_util, only: perm3d
 use mat_lib, only: norm=> normfro, thor_matmul_mnk
 implicit none
 double precision, allocatable   :: WAX1(:,:)
 double precision, intent(INOUT) :: nrm
 double precision, intent(IN)    :: WAX2(:,:), w(:,:,:,:), A2(:,:), x(:,:,:,:)
 integer, intent(IN)             :: rw1, n, k, rw2, m, ra1, ra2, rx1, rx2
 double precision, intent(IN), optional :: extnrm_
 !
 double precision, allocatable :: w3(:,:,:), xc(:,:,:), xc2(:,:)

    !% Right WAX has the form of the last matrix TT block, i.e. [ra, rw, rx]
    nrm= 0d0
    w3 = reshape(w, [rw1*n, k, rw2])
    call perm3d(w3, 132)
    xc = reshape(x, [rx1*m, k, rx2])
    call perm3d(xc, 213)
    WAX1= thor_matmul_mnk(k*rx1*m, ra2*rw2, rx2, xc, WAX2, tsp_='nt')
    WAX1= transpose(reshape(WAX1, [k*rx1, m*ra2*rw2]))
    WAX1= thor_matmul_mnk(ra1*n, rw2*k*rx1, m*ra2, A2, WAX1)
    WAX1= transpose(reshape(WAX1, [ra1, n*rw2*k*rx1]))
    WAX1= thor_matmul_mnk(rw1, rx1*ra1, n*rw2*k, w3, WAX1)
    WAX1= transpose(reshape(WAX1, [rw1*rx1, ra1]))
    nrm= max(nrm, norm(WAX1))

    !% Extract the scale to prevent overload
    if (.not.present(extnrm_)) then
       if (nrm > 0d0) then
           WAX1= WAX1/nrm
       else
          nrm = 1d0
       endif
    else
       !% Override the normalization
       WAX1= WAX1/extnrm_
    endif
 end function rightreduce_matrix


 !% Accumulates the left reduction W{1:k}'*X{1:k}
 function leftreduce_vector(WX1,w,x,rw1,n,k,rw2,rx1,rx2,extnrm_) result(WX2)
 use mat_lib, only: thor_matmul, thor_matmul_mnk
 implicit none
 double precision, allocatable :: WX2(:,:)
 double precision, intent(IN)  :: WX1(:,:),w(:,:),x(:,:,:,:)
 integer, intent(IN) :: rw1, n, k, rw2, rx1, rx2
 double precision, intent(IN), optional :: extnrm_
 !
 double precision, allocatable :: w2(:,:), tmp(:,:)
    WX2= thor_matmul_mnk(n*k*rw2,rx1,rw1,w,WX1,tsp_='tn')
    WX2= transpose(reshape(WX2, [n*k, rw2*rx1]))
    WX2= thor_matmul_mnk(rw2, rx2, rx1*n*k, WX2, x)
    if (present(extnrm_)) then
       WX2= WX2/extnrm_
    endif
 end function leftreduce_vector


 !% Accumulates the right reduction W{k:d}'*X{k:d}
 function rightreduce_vector(WX2,w,x,rw1,n,k,rw2,rx1,rx2,extnrm_) result(WX1)
 use mat_lib, only: thor_matmul_mnk
 implicit none
 double precision, allocatable :: WX1(:,:)
 double precision, intent(IN)  :: WX2(:,:),w(:,:),x(:,:,:,:)
 integer, intent(IN) :: rw1, n, k, rw2, rx1, rx2
 double precision, intent(IN), optional :: extnrm_
    WX1= thor_matmul_mnk(rx1*n*k, rw2, rx2, x, WX2)
    WX1= thor_matmul_mnk(rx1, rw1, n*k*rw2, WX1, w, tsp_='nt')                               !    WX1{1,i} = WX1{1,i}*wc'; % size rx1, rw1
    if (present(extnrm_)) WX1= WX1/extnrm_
 end function rightreduce_vector


 !% A matrix-matrix product for the matrix in the 3D TT (WAX1-A-WAX2), and
 !% full matrix of size (rx1*m*k*rx2). Returns (rw1*n*k*rw2)
 function local_matvec(x, rx1,m,k,rx2, rw1,n,rw2, WAX1, A, WAX2, ra1,ra2) result(w)
 use mat_lib, only: thor_matmul, thor_matmul_mnk
 use matrix_util, only: perm3d
 implicit none
 double precision, allocatable :: w(:,:)
 double precision, intent(IN) :: x(:,:,:,:), WAX1(:,:), A(:,:), WAX2(:,:)
 integer, intent(IN) :: rx1, m, k, rx2, rw1, n, rw2, ra1, ra2
 !
 double precision, allocatable :: w2(:,:), w3(:,:,:), xc(:,:,:), xc2(:,:), tmp(:,:), wk(:,:)
    allocate(w2(rw1*n*rw2, k)); w2= 0d0
    xc= reshape(x, [rx1*m, k, rx2])
    call perm3d(xc, 213)
    wk= thor_matmul_mnk(k*rx1*m, ra2*rw2, rx2, xc, WAX2, tsp_='nt')
    wk= transpose(reshape(wk, [k*rx1, m*ra2*rw2]))
    wk= thor_matmul_mnk(ra1*n, rw2*k*rx1, m*ra2, A, wk)
    wk= transpose(reshape(wk, [ra1*n*rw2*k, rx1]))
    wk= thor_matmul_mnk(rw1, n*rw2*k, rx1*ra1, WAX1, wk)
    wk= reshape(wk, [rw1*n*rw2, k])
    w2= w2 + wk
    w3= reshape(w2, [rw1*n, rw2, k])
    call perm3d(w3, 132)
    w= reshape(w3, [rw1*n*rw2*k, 1])
 end function local_matvec


end module ttamen_lib
