!> @file dtt_amen.f90
!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin (korobkin@lanl.gov)
!! @date   February 2025
!! @brief  Benchmarking amen_mv and amen_solve
!!

program main
use rnd_lib, only: random, init_random_seed
implicit none
character(len=128) arg
!
integer,parameter :: Ntim = 10
double precision  :: timers(Ntim)
integer :: i, err
integer :: rseed = 0 ! random seed (0 = from clock)
double precision  :: epsi = 1d-8
integer           :: batch_size = 10
integer           :: low_rank = 5
integer           :: num_sweeps = 1
integer           :: d = 4
integer           :: num_threads = 1

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-e') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) epsi
        if (err.ne.0) error stop "ERROR: when parsing -e <???>"
        i= i + 1
     elseif (trim(arg).eq.'-l') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) low_rank
        if (err.ne.0) error stop "ERROR: when parsing -l <???>"
        i= i + 1
     elseif (trim(arg).eq.'-d') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) d
        if (err.ne.0) error stop "ERROR: when parsing -d <???>"
        i= i + 1
     elseif (trim(arg).eq.'-n') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) num_sweeps
        if (err.ne.0) error stop "ERROR: when parsing -n <???>"
        i= i + 1
     elseif (trim(arg).eq.'-b') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) batch_size
        if (err.ne.0) error stop "ERROR: when parsing -b <???>"
        i= i + 1
     elseif (trim(arg).eq.'-t') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) num_threads
        if (err.ne.0) error stop "ERROR: when parsing -t <???>"
        i= i + 1
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     elseif (trim(arg).eq.'-h') then
        call help_msg
        stop
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call init_random_seed(rseed)
  !call bench_amen(d, low_rank, num_sweeps, batch_size, epsi, num_threads)
  call bench_amen_mm(d, low_rank, num_sweeps, batch_size, epsi, num_threads)

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Benchmark amen_mv and amen_solve", &
   "Usage: ./dtt_amen.exe [-h]", &
   "Options:", &
   " -d <val>   : tensor order (dimensionality) [4]", &
   " -l <val>   : maximal rank of the generated tensor [5]", &
   " -n <val>   : number of sweeps [1]", &
   " -e <val>   : rounding tolerance epsilon [1e-8]", &
   " -b <val>   : batch size [10]", &
   " -t <val>   : number of threads [1]", &
   " --rseed <num> : supply random seed [none]", &
   " -h         : print this help message"
  end subroutine


  subroutine bench_amen(d, low_rank, nswp, batch_size, epsi, num_threads)
  use time_lib
  use thor_lib, only: dtt_tensor, dtt_matrix, dtt_tensor_rand, dtt_matrix_rand
  use ttamen_lib, only: amen_solve, amen_mv
  implicit none
  integer, intent(IN) :: d          !< tensor order
  integer, intent(IN) :: low_rank   !< maximal rank of the generated tensor
  integer, intent(IN) :: nswp       !< number of sweeps to use
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  integer, intent(IN) :: num_threads
  !
  integer, allocatable :: nn(:), qq(:), ss(:), rr(:)
  integer :: mode_size, i, test
  type(dtt_tensor) :: x, b
  type(dtt_matrix) :: A
  double precision :: t_amen_mv, t_amen_solve

     call openblas_set_num_threads(num_threads)
     allocate(nn(d), qq(d), ss(d), rr(0:d))
     rr(0)= 1; rr(1:d-1)= low_rank; rr(d)= 1

     call system_timer_init()
     print '("# THOR: amen_mv optimizing")'
     print '("# Tensor order = "I5", low rank = "I7)', d, low_rank
     print '("# Modes = (M, M, ..., M)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:amen-mv[s] 3:amen-solve[s]")'
     mode_size = 2
     scaling_loop: do i=1,12
        qq(1:d) = mode_size
        ss(1:d) = mode_size
        nn = qq*ss
        x = dtt_tensor_rand(qq, r_=rr)
        A = dtt_matrix_rand(qq, ss, r_=rr)
        !call x% pprint(label_="###  x  ###")
        timers(1:Ntim)= 0d0
        call system_timer_start
        do test = 1, batch_size
           call amen_mv(b, A, x, epsi, nswp_=nswp, verb_=0)
        enddo
        call system_timer_stop
        t_amen_mv = system_dt/dble(batch_size)
        timers(1:Ntim)= timers(1:Ntim)/dble(batch_size)

        !call system_timer_start
        !do test = 1, batch_size
        !   x = amen_solve(A, b, tol=epsi, nswp=nswp, verb_=0)
        !enddo
        !call system_timer_stop
        !t_amen_solve = system_dt/dble(batch_size)
        print '(I10, 20(1X,ES10.3))', mode_size, t_amen_mv, t_amen_solve
        !                            , timers(1:Ntim)
        mode_size = 2*mode_size
     enddo scaling_loop

     deallocate(nn, qq, ss, rr)

  end subroutine bench_amen


  subroutine bench_amen_mm(d, low_rank, nswp, batch_size, epsi, num_threads)
  use time_lib
  use thor_lib, only: dtt_matrix, dttm_random_ortho
  use thorio_lib, only: dttm_read_ascii_file, dttm_write_ascii_file
  use ttamen_lib, only: amen_mm
  implicit none
  integer, intent(IN) :: d          !< tensor order
  integer, intent(IN) :: low_rank   !< maximal rank of the generated tensor
  integer, intent(IN) :: nswp       !< number of sweeps to use
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  integer, intent(IN) :: num_threads
  !
  integer, allocatable :: nn(:), qq(:), ss(:), rr(:)
  integer :: mode_size, i, test
  type(dtt_matrix) :: A, B, C
  double precision :: t_amen_mm

     call openblas_set_num_threads(num_threads)
     allocate(qq(d), ss(d), rr(0:d))
     rr(0)= 1; rr(1:d-1)= low_rank; rr(d)= 1

     call system_timer_init()
     print '("# THOR: amen_mm, number of threads: "I4)', num_threads
     print '("# Tensor order = "I5", low rank = "I7)', d, low_rank
     print '("# Modes = (MxM, MxM, ..., MxM)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:amen-mm[s]")'
     mode_size = 2
     scaling_loop: do i=1,12
        qq(1:d) = mode_size
        ss(1:d) = mode_size
        A = dttm_random_ortho(qq, ss, r_=rr)
        B = dttm_random_ortho(qq, ss, r_=rr)
!call dttm_write_ascii_file(A, 'bench/A.dat', verb_=.true.)
!call dttm_write_ascii_file(B, 'bench/B.dat', verb_=.true.)
!A = dttm_read_ascii_file('/tmp/A.dat', verb_=.true.)
!B = dttm_read_ascii_file('/tmp/B.dat', verb_=.true.)
        call system_timer_start
        do test = 1, batch_size
           call amen_mm(C, A, B, epsi, nswp_=nswp, verb_=0)
        enddo
        call system_timer_stop
        t_amen_mm = system_dt/dble(batch_size)

        print '(I10, 20(1X,ES10.3))', mode_size, t_amen_mm
        mode_size = 2*mode_size
     enddo scaling_loop

     deallocate(qq, ss, rr)

  end subroutine bench_amen_mm



!!  !> function [X]=amen_mm(A, Y, tol, varargin)
!!  !! Approximate the matrix-by-matrix via the AMEn iteration
!!  !!    [X]=amen_mv(A, Y, tol, varargin)
!!  !!    Attempts to approximate the X = A*Y with accuracy TOL using the
!!  !!    AMEn+ALS iteration. A is a n x m matrix, Y is a m x k matrix.
!!  !!    Matrices A,Y can be given in tt_matrix format, the output is tt_matrix.
!!  !!    Y can be tt_tensor, it's considered as column, the output is tt_tensor.
!!  !!    A and Y can be a {d,R} cell array. However, X is always a "single" TT
!!  !!    (no tensor chain), since that's how ALS works. Generally, X has the
!!  !!    same form as Y, except that it's always {d,1} in essense. X and Y can't
!!  !!    be sparse (SVD will destroy it anyways), but A can be.
!!  !!
!!  !!    Options are provided in form
!!  !!    'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
!!  !!    on. The parameters are set to default (in brackets in the following)
!!  !!    The list of option names and default values are:
!!  !!        o x0 - initial approximation to AY [rand with ranks of Y(:,1)]
!!  !!        o nswp - maximal number of sweeps [20]
!!  !!        o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
!!  !!        o kickrank - compression rank of the error,
!!  !!                     i.e. enrichment size [4]
!!  !!
!!  !!
!!  !! ********
!!  !!    For description of adaptive ALS please see
!!  !!    Sergey V. Dolgov, Dmitry V. Savostyanov,
!!  !!    Alternating minimal energy methods for linear systems in higher dimensions.
!!  !!    Part I: SPD systems, http://arxiv.org/abs/1301.6068,
!!  !!    Part II: Faster algorithm and application to nonsymmetric systems, http://arxiv.org/abs/1304.1222
!!  !!
!!  !! Call tree:
!!  !! +- amen_mm
!!  !!    +-- leftreduce_matrix
!!  !!    +-- leftreduce_vector
!!  !!    +-- local_matvec
!!  !!    +-- rightreduce_matrix
!!  !!    +-- rightreduce_vector
!!
!!  subroutine amen_mm(X, A, Y, tol, x0_, nswp_, verb_, kickrank_, init_qr_)
!!  use matlab_struct_module, only: array2d, array3d, array4d, rand_cell3d, pprint_matrix
!!  use mat_lib, only: svd, svdgram, bfun3, my_chop2, stack, &
!!                     diag => matlab_diag, zeros => matlab_zeros, &
!!                     qr => matlab_qr, norm => normfro
!!  use rnd_lib, only: randn
!!  use thor_lib, only: dtt_tensor, dtt_matrix, core2cell, dttm_random_ortho
!!  use mat_lib, only: thor_matmul, thor_matmul_mnk
!!  use thorio_lib, only: dttm_read_ascii_file, dttm_write_ascii_file
!!  use time_lib
!!  implicit none
!!  class(dtt_matrix),intent(INOUT) :: X
!!  class(dtt_matrix),intent(IN) :: A, Y
!!  double precision, intent(IN) :: tol
!!  class(dtt_matrix),intent(IN), optional:: x0_ !< initial approximation to A@Y
!!  integer, intent(IN),optional :: nswp_        !< maximal number of sweeps [20]
!!  integer, intent(IN),optional :: verb_        !< verbosity level, 0-silent, 1-sweep info, 2-block info [1]
!!  integer, intent(IN),optional :: kickrank_    !< compression rank of the error,
!!                                               !  i.e. enrichment size [3]
!!  logical, intent(IN),optional :: init_qr_     !< perform QR of the input (save some time in ts, etc) [true]
!!  !
!!  type(array3d), allocatable:: z(:) !< z0 - initial approximation to the error Ax-y [rand rank-kickrank]
!!  integer, allocatable :: m(:), k(:), rx(:),n(:),ra(:),ry(:),rz(:),r2(:)
!!  integer :: d,dir,i,j,swp,r,Rry,Rra, sh(2)
!!  integer :: nswp,verb,kickrank,renorm
!!  logical :: init_qr,fkick
!!  integer :: kickrank2, Ntim
!!  type(array3d), allocatable :: phiax(:), phiyax(:), phizax(:), phizy(:)
!!  type(array2d), allocatable :: y2(:), A2(:)
!!  double precision :: dx,max_dx,nrm,nrm_s,nrmr,nrmz, aux
!!  double precision, allocatable :: s2(:,:), cr4(:,:,:,:)
!!  double precision, allocatable, target :: cr2(:,:),crx(:,:),cry(:,:),dum__(:,:)
!!  double precision, allocatable :: nrms(:),crs(:,:),ys(:,:),yz(:,:),cr1(:,:)
!!  double precision, allocatable, target :: cr(:,:),s(:),RR(:,:),u(:,:),v(:,:),crz(:,:)
!!  double precision, pointer   :: crp(:,:),up(:,:),sp(:),vp(:,:),RRp(:,:)
!!  double precision, parameter :: one2(1,1) = reshape([1d0], [1,1])
!!  double precision, allocatable :: XAY2(:,:)
!!  type(array2d), allocatable  :: XAY(:), ZAY(:), ZX(:)
!! !!double precision :: dt1, dt2, dt3, dt3a, dt3b, dt3c, dt3d, dt4, dt5, dt6, dt7, dt8, dt9, dt10
!!
!! !!call system_timer_start
!!  nswp= 20;       if(present(nswp_)) nswp= nswp_             ! nswp = 20;
!!  kickrank= 4;    if(present(kickrank_)) kickrank= kickrank_ ! kickrank = 4;
!!  verb= 1;        if(present(verb_)) verb= verb_             ! verb = 1;
!!  init_qr=.true.; if(present(init_qr_)) init_qr= init_qr_    ! init_qr = true;
!!                              ! X = [];
!!                              !
!!                              ! for i=1:2:length(varargin)-1
!!                              !     switch lower(varargin{i})
!!                              !         case 'nswp'
!!                              !             nswp=varargin{i+1};
!!                              !         case 'x0'
!!                              !             X=varargin{i+1};
!!                              !         case 'verb'
!!                              !             verb=varargin{i+1};
!!                              !         case 'kickrank'
!!                              !             kickrank=varargin{i+1};
!!                              !         case 'init_qr'
!!                              !             init_qr = varargin{i+1};
!!                              !
!!                              !         otherwise
!!                              !             warning('Unrecognized option: %s\n',varargin{i});
!!                              !     end
!!                              ! end
!!                              !
!!                              ! % Grumble inputs
!!                              ! % First, the vector
!!                              ! if (isa(Y,'tt_tensor'))
!!                              !     d = Y.d;
!!                              !     m = Y.n;
!!                              !     k = ones(d,1);
!!                              !     ry = Y.r;
!!                              !     Ry = 1;
!!                              !     Y = core2cell(Y);
!!                              !     outtype = 0;
!!                              ! elseif (isa(Y,'tt_matrix'))
!!  d= Y% m - Y% l + 1          !     d = Y.d;
!!  allocate(m(d),n(d),k(d),ry(d+1))
!!  m= Y% q(1:d)                !     m = Y.n;
!!  k= Y% s(1:d)                !     k = Y.m;
!!  ry=Y% r(0:d)                !     ry = Y.r;
!!  Rry= 1                      !     Ry = 1;
!!
!!                              !     Y = core2cell(Y);
!!                              !     outtype = 1;
!!                              ! else % {d,R}
!!                              !     [d,Ry] = size(Y);
!!                              !     m = zeros(d,1);
!!                              !     k = zeros(d,1);
!!                              !     ry = ones(d+1,Ry);
!!                              !     for j=1:Ry
!!                              !         for i=1:d
!!                              !             [r1,n1,n2,r2]=size(Y{i,j});
!!                              !             if (r1~=ry(i,j))
!!                              !                 error('Inconsistent ry(%d)', i);
!!                              !             end;
!!                              !             m(i) = n1;
!!                              !             k(i) = n2;
!!                              !             ry(i+1,j) = r2;
!!                              !         end;
!!                              !     end;
!!                              !     outtype = 2;
!!                              ! end;
!!                              !
!! ! % Grumble matrix
!!                                                                        ! if (isa(A, 'tt_matrix'))
!!                                                                        !     % Copy TT ranks and blocks from a tt_matrix
!!  if (A%m-A%l+1.ne.d) error stop "[!][amen_mm]: mismatched orders in A and Y"
!!  allocate(ra(d+1)); ra(1:d+1)= A% r(0:d)                               !     ra = A.r;
!!  if (any(A%s(1:d).ne.m(1:d))) &                                        !     if (A.d~=d)||(any(A.m~=m))
!!     error stop "[!][amen_mm]: Mismatching TT dimensionalities in A and Y" !      error('Mismatching TT dimensionalities in A and Y');
!!                                                                        !     end;
!!  n= A% q(1:d)                                                          !     n = A.n;
!!  ! Ra = 1 always for tt-matrix input (not sparse)                      !     Ra = 1;
!!
!!  allocate(A2(d))                                                       !     A = core2cell(A);
!!  do i=1,d                                                              !     for i=1:d
!!     A2(i)%arr = reshape(A% u4(i)% p, [ra(i)*n(i), m(i)*ra(i+1)])       !         A{i} = reshape(A{i}, ra(i)*n(i), m(i)*ra(i+1));
!!  enddo                                                                 !     end;
!!                                                                        ! else
!!                                                                        !     % The matrix is given as a {d,R} cell array, meaning a sum
!!                                                                        !     % of R TT formats (needed for sparse blocks, in which case
!!                                                                        !     % it's the sum of R Kronecker products)
!!                                                                        !     Ra = size(A,2);
!!                                                                        !     if (size(A,1)~=d)
!!                                                                        !         error('Mismatching TT dimensionalities in A and Y');
!!                                                                        !     end;
!!                                                                        !     ra = ones(d+1,Ra);
!!                                                                        !     n = ones(d,1);
!!                                                                        !     for j=1:Ra
!!                                                                        !         for i=1:d
!!                                                                        !             if (issparse(A{i,j}))
!!                                                                        !                 % Sparse matrix can only have TT ranks 1
!!                                                                        !                 [n1,n2]=size(A{i,j});
!!                                                                        !                 r2 = n2/m(i);
!!                                                                        !                 n1 = n1/ra(i,j);
!!                                                                        !                 if(abs(r2-round(r2))>sqrt(eps))
!!                                                                        !                     error('A is sparse, but the column size are not divisible by m');
!!                                                                        !                 end;
!!                                                                        !                 n2 = m(i);
!!                                                                        !             else
!!                                                                        !                 [r1,n1,n2,r2]=size(A{i,j});
!!                                                                        !                 if (r1~=ra(i,j))
!!                                                                        !                     error('Inconsistent ra(%d)', i);
!!                                                                        !                 end;
!!                                                                        !                 A{i,j} = reshape(A{i,j}, r1*n1, n2*r2);
!!                                                                        !             end;
!!                                                                        !             if (n2~=m(i))
!!                                                                        !                 error('Mismatching m in A');
!!                                                                        !             end;
!!                                                                        !             n(i) = n1;
!!                                                                        !             ra(i+1,j) = r2;
!!                                                                        !         end;
!!                                                                        !     end;
!!                                                                        ! end;
!!                                                                        !
!! !!call system_timer_stop
!! !!dt1 = system_dt
!! !!call system_timer_start
!!  !% Initial guess
!!  if (.not.present(x0_)) then                                           ! if (isempty(X))
!!     X = dttm_random_ortho(n(1:d), k(1:d), r_=ry(1:d+1))                !     X = tt_rand(n.*k, d, ry(:,1), 1);
!! !!call dttm_write_ascii_file(X, 'bench/X.dat', verb_=.true.)
!! !!!DEBUG!X = dttm_read_ascii_file('bench/X.dat', verb_=.true.)
!! !!!DEBUG!
!!     allocate(rx(d+1)); rx(1:d+1)= X% r(0:d)                            !     rx = X.r;
!!                                                                        !     X = core2cell(X);
!!     init_qr= .false.                                                   !     init_qr = false;
!!  else                                                                  ! else
!!     X = x0_
!!                                                                        !     % X was given, parse it
!!                                                                        !     if (isa(X,'tt_tensor'))
!!                                                                        !         if (X.d~=d)||(any(X.n~=n.*k))
!!                                                                        !             error('Mismatching dimensions in x0');
!!                                                                        !         end;
!!                                                                        !         rx = X.r;
!!                                                                        !         X = core2cell(X);
!!                                                                        !     elseif (isa(X,'tt_matrix'))
!!     if (X%m-X%l+1 /= d .or. any(X%q(1:d)/=n) .or. any(X%s(1:d)/=k)) &  !         if (X.d~=d)||(any(X.n~=n))||(any(X.m~=k))
!!        error stop "[!][amen_mm]: mismatched dimensions in X & Y"       !             error('Mismatching dimensions in x0');
!!                                                                        !         end;
!!     allocate(rx(d+1)); rx(1:d+1)= X% r(0:d)                            !         rx = X.r;
!!                                                                        !         X = core2cell(X);
!!                                                                        !     else % {d,R}
!!                                                                        !         X = X(:,end); % There could be several TC columns
!!                                                                        !         rx = ones(d+1,1);
!!                                                                        !         for i=1:d
!!                                                                        !             [r1,n1,n2,r2]=size(X{i});
!!                                                                        !             if (r1~=rx(i))
!!                                                                        !                 error('Inconsistent rx(%d)', i);
!!                                                                        !             end;
!!                                                                        !             if (n1~=n(i))||(n2~=k(i))
!!                                                                        !                 error('Inconsistent n/k in x0');
!!                                                                        !             end;
!!                                                                        !             rx(i+1,j) = r2;
!!                                                                        !         end;
!!                                                                        !     end;
!!  endif                                                                 ! end;
!!                                                                        !
!!  !% Reductions
!!  allocate(XAY(d+1))                                                    ! XAY = cell(d+1,Ra,Ry);
!!  allocate(ZAY(d+1))                                                    ! ZAY = cell(d+1,Ra,Ry);
!!  allocate(ZX(d+1))                                                     ! ZX = cell(d+1,1);
!!  ZX(1)%arr = one2; ZX(d+1)%arr = one2                                  ! ZX{1} = 1; ZX{d+1} = 1;
!!                                                                        ! for j=1:Ry
!!                                                                        !     for i=1:Ra
!!  XAY(1)%arr= one2; XAY(d+1)%arr= one2                                  !         XAY{1,i,j} = 1; XAY{d+1,i,j} = 1;
!!  ZAY(1)%arr= one2; ZAY(d+1)%arr= one2                                  !         ZAY{1,i,j} = 1; ZAY{d+1,i,j} = 1;
!!                                                                        !     end;
!!                                                                        ! end;
!!  !% Residual rank
!!  allocate(rz(d+1)); rz= kickrank; rz(1)=1; rz(d+1)=1                   ! rz = [1; kickrank*ones(d-1,1); 1];
!!  allocate(nrms(d)); nrms= 1d0                                          ! nrms = ones(d,1); % To prevent doulbe overflow
!!
!! !!call system_timer_stop
!! !!dt2 = system_dt; dt3a=6.66d0; dt3b=6.77d0
!! !!call system_timer_start
!!  !% Initial ort
!!  do i=1,d-1                                                            ! for i=1:d-1
!!     if (init_qr) then                                                  !     if (init_qr)
!!        cr= reshape(X%u4(i)%p, [rx(i)*n(i)*k(i), rx(i+1)])              !         cr = reshape(X{i}, rx(i)*n(i)*k(i), rx(i+1));
!!        call qr(cr1, RR, cr)                                            !         [cr,R]=qr(cr, 0);
!!        cr = cr1
!!        nrmr= norm(RR)                                                  !         nrmr = norm(R, 'fro');
!!        if (nrmr > 0d0) then                                            !         if (nrmr>0)
!!           RR= RR/nrmr                                                  !             R = R/nrmr;
!!        endif                                                           !         end;
!!        !cr2= reshape(X%u4(i+1)%p, [rx(i+1), n(i+1)*k(i+1)*rx(i+2)])     !         cr2 = reshape(X{i+1}, rx(i+1), n(i+1)*k(i+1)*rx(i+2));
!!        cr2= thor_matmul_mnk(rx(i+1), rx(i+1), n(i+1)*k(i+1)*rx(i+2), &
!!                             RR, X%u4(i+1)%p)                           !         cr2 = R*cr2;
!!        rx(i+1)= size(cr, 2)                                            !         rx(i+1) = size(cr, 2);
!!        deallocate(X% u4(i)% p)
!!        allocate(X% u4(i)% p(rx(i), n(i), k(i), rx(i+1)))
!!        X% u4(i)% p= reshape(cr, [rx(i), n(i), k(i), rx(i+1)])          !         X{i} = reshape(cr, rx(i), n(i), k(i), rx(i+1));
!!        deallocate(X% u4(i+1)% p)
!!        allocate(X% u4(i+1)% p(rx(i+1), n(i+1), k(i+1), rx(i+2)))
!!        X% u4(i+1)% p= reshape(cr2, [rx(i+1), n(i+1), k(i+1), rx(i+2)]) !         X{i+1} = reshape(cr2, rx(i+1), n(i+1), k(i+1), rx(i+2));
!!     endif                                                              !     end;
!!     !% Reduce
!! !!call system_timer_stop
!! !!dt3a = system_dt
!! !!call system_timer_start
!!     ![XAY(i+1,:,:),nrms(i)] = leftreduce_matrix(XAY(i,:,:), X{i},
!!     !                         A(i,:), Y(i,:), rx(i),n(i),k(i),rx(i+1),
!!     !                         Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:));
!!     XAY(i+1)%arr= leftreduce_matrix(nrms(i), XAY(i)%arr, X% u4(i)% p, &
!!                            A2(i)%arr, Y% u4(i)% p, rx(i),n(i),k(i),rx(i+1),      &
!!                            ra(i),ra(i+1),ry(i),m(i),ry(i+1))
!! !DEBUG!call pprint_matrix(XAY(i+1)%arr)
!!     !% Residual reductions
!! !!call system_timer_stop
!! !!dt3b = system_dt
!! !!call system_timer_start
!!     if (kickrank > 0) then                                             !     if (kickrank>0)
!!        cr= randn(rz(i)*n(i)*k(i), rz(i+1))                            !         cr = randn(rz(i)*n(i)*k(i), rz(i+1));
!!        call qr(cr1, RR, cr)                                            !         [cr,~]=qr(cr, 0);
!!        rz(i+1)= size(cr1, 2)                                           !         rz(i+1) = size(cr, 2);
!!        ! ZAY(i+1,:,:) = leftreduce_matrix(ZAY(i,:,:), cr,
!!        !                A(i,:), Y(i,:), rz(i),n(i),k(i),rz(i+1),
!!        !                Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:), nrms(i));
!!        ZAY(i+1)%arr= leftreduce_matrix(aux, ZAY(i)% arr, &
!!                            reshape(cr1,[size(cr1),1,1,1]), &
!!                            A2(i)%arr, Y% u4(i)% p, rz(i),n(i),k(i),rz(i+1), &
!!                            ra(i),ra(i+1), ry(i),m(i),ry(i+1), extnrm_=nrms(i))
!! !!call system_timer_stop
!! !!dt3c = system_dt
!! !!call system_timer_start
!! !DEBUG!call pprint_matrix(ZAY(i+1)%arr)
!!        ! ZX(i+1) = leftreduce_vector(ZX(i), cr, X(i),
!!        !           rz(i),n(i),k(i),rz(i+1), 1,rx(i),rx(i+1));
!!        ZX(i+1)%arr= leftreduce_vector(ZX(i)%arr, cr1, X%u4(i)%p, &
!!                        rz(i),n(i),k(i),rz(i+1), rx(i),rx(i+1))
!! !DEBUG!call pprint_matrix(ZX(i+1)%arr)
!!     endif                                                              !     end;
!!  enddo                                                                 ! end;
!!                                                                        !
!!  !% Normalize the last block for small error on good initial guess
!!  nrmr= norm(X% u4(d)% p)                                               ! nrmr = norm(X{d}(:));
!!  if (nrmr > 0d0) then                                                  ! if (nrmr>0)
!!     X% u4(d)% p= X% u4(d)% p/nrmr                                      !     X{d} = X{d}/nrmr;
!!  endif                                                                 ! end
!!                                                                        !
!!  i = d                                                                 ! i = d;
!!  dir = -1                                                              ! dir = -1;
!!  swp = 1                                                               ! swp = 1;
!!  max_dx = 0d0                                                          ! max_dx = 0;
!!  !% Iteration
!! !!call system_timer_stop
!! !!dt3 = system_dt
!! !!call system_timer_start
!!  sweeps: do while(swp <= nswp)                                         ! while (swp<=nswp)
!!     !% Project the MatVec generating vector
!!     !  cr = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:),
!!     !       rx(i),n(i),rx(i+1), XAY(i,:,:), A(i,:), XAY(i+1,:,:),
!!     !       Ra,ra(i,:),ra(i+1,:));
!!     cr= local_matvec(Y% u4(i)% p, ry(i), m(i), k(i), ry(i+1), &
!!             rx(i),n(i),rx(i+1), XAY(i)%arr, A2(i)%arr, XAY(i+1)%arr, &
!!             ra(i),ra(i+1))
!!
!!     nrms(i)= norm(cr)                                                  !     nrms(i) = norm(cr, 'fro');
!!     !% The main goal is to keep y{i} of norm 1
!!     if (nrms(i) > 0d0) then                                            !     if (nrms(i)>0)
!!        cr= cr/nrms(i)                                                  !         cr = cr/nrms(i);
!!     else                                                               !     else
!!        nrms(i)= 1d0                                                    !         nrms(i)=1;
!!     endif                                                              !     end;
!!                                                                        !     X{i} = reshape(X{i}, rx(i)*n(i)*k(i)*rx(i+1), 1);
!!     dx= norm(cr - reshape(X% u4(i)% p, [rx(i)*n(i)*k(i)*rx(i+1), 1]))  !     dx = norm(cr-X{i});
!!     max_dx= max(dx, max_dx)                                            !     max_dx = max(max_dx, dx);
!! !!call system_timer_stop
!! !!dt4 = system_dt
!! !!call system_timer_start
!!
!!     !% Truncation and enrichment
!!     if (dir > 0 .and. i < d) then                                      !     if ((dir>0)&&(i<d))
!!        cr= reshape(cr, [rx(i)*n(i)*k(i), rx(i+1)])                     !         cr = reshape(cr, rx(i)*n(i)*k(i), rx(i+1));
!!        call svd(cr,up,sp,vp)                                           !         [u,s,v]=svd(cr, 'econ');
!!        s= sp                                                           !         s = diag(s);
!!        r= my_chop2(s, tol*norm(s)/sqrt(dble(d)))                       !         r = my_chop2(s, tol*norm(s)/sqrt(d));
!!        u= up(:,1:r)                                                    !         u = u(:,1:r);
!!        v= thor_matmul(vp(1:r,:),diag(s(1:r)),'tn')                     !         v = conj(v(:,1:r))*diag(s(1:r));
!!
!!        !% Prepare enrichment, if needed
!!        if (kickrank > 0) then                                          !         if (kickrank>0)
!!           cr= thor_matmul(u, v, tsp_='nt')                             !             cr = u*v.';
!!           cr= reshape(cr, [rx(i)*n(i)*k(i), rx(i+1)])                  !             cr = reshape(cr, rx(i)*n(i)*k(i), rx(i+1));
!!           !% For updating z
!!           !  crz = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:),
!!           !        rz(i),n(i),rz(i+1), ZAY(i,:,:), A(i,:),
!!           !        ZAY(i+1,:,:), Ra,ra(i,:),ra(i+1,:));
!!           crz= local_matvec(Y% u4(i)% p, ry(i),m(i),k(i),ry(i+1), &
!!                    rz(i),n(i),rz(i+1), ZAY(i)%arr, A2(i)%arr, &
!!                    ZAY(i+1)%arr, ra(i),ra(i+1))
!!           crz= reshape(crz, [rz(i)*n(i)*k(i), rz(i+1)])                !             crz = reshape(crz, rz(i)*n(i)*k(i), rz(i+1));
!!           ys= thor_matmul(cr, ZX(i+1)%arr)                             !             ys = cr*ZX{i+1};
!!           yz= reshape(ys, [rx(i), n(i)*k(i)*rz(i+1)])                  !             yz = reshape(ys, rx(i), n(i)*k(i)*rz(i+1));
!!           yz= thor_matmul(ZX(i)%arr, yz)                               !             yz = ZX{i}*yz;
!!           yz= reshape(yz, [rz(i)*n(i)*k(i), rz(i+1)])                  !             yz = reshape(yz, rz(i)*n(i)*k(i), rz(i+1));
!!           crz= crz/nrms(i) - yz                                        !             crz = crz/nrms(i) - yz;
!!           nrmz= norm(crz)                                              !             nrmz = norm(crz,'fro');
!!           !% For adding into solution
!!           !  crs = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:),
!!           !        rx(i),n(i),rz(i+1), XAY(i,:,:), A(i,:),
!!           !        ZAY(i+1,:,:), Ra,ra(i,:),ra(i+1,:));
!!           crs= local_matvec(Y% u4(i)% p, ry(i),m(i),k(i),ry(i+1), &
!!                    rx(i),n(i),rz(i+1), XAY(i)%arr, A2(i)%arr, &
!!                    ZAY(i+1)%arr, ra(i),ra(i+1))
!!           crs= reshape(crs, [rx(i)*n(i)*k(i), rz(i+1)])                !             crs = reshape(crs, rx(i)*n(i)*k(i), rz(i+1));
!!           crs= crs/nrms(i) - ys                                        !             crs = crs/nrms(i) - ys;
!!           u= stack(u, crs)                                             !             u = [u,crs];
!!           call qr(cr1, RR, u)                                          !             [u,R]=qr(u, 0);
!!           u= cr1
!!           v= stack(v, zeros(rx(i+1), rz(i+1)))                         !             v = [v, zeros(rx(i+1), rz(i+1))];
!!           v= thor_matmul(v, RR, tsp_='nt')                             !             v = v*R.';
!!           r= size(u, 2)                                                !             r = size(u, 2);
!!        endif                                                           !         end;
!!        deallocate(X% u4(i)% p)
!!        allocate(X% u4(i)% p(rx(i),n(i),k(i),r))
!!        call dcopy(rx(i)*n(i)*k(i)*r, u, 1, X% u4(i)% p, 1)             !         X{i} = reshape(u, rx(i), n(i), k(i), r);
!!        cr2= reshape(X% u4(i+1)% p, [rx(i+1), n(i+1)*k(i+1)*rx(i+2)])   !         cr2 = reshape(X{i+1}, rx(i+1), n(i+1)*k(i+1)*rx(i+2));
!!        v= reshape(v, [rx(i+1),r])                                      !         v = reshape(v, rx(i+1), r);
!!        cr2= thor_matmul(v, cr2, tsp_='tn')                             !         cr2 = v.'*cr2;
!!        deallocate(X% u4(i+1)% p)
!!        allocate(X% u4(i+1)% p(r,n(i+1),k(i+1),rx(i+2)))
!!        call dcopy(r*n(i+1)*k(i+1)*rx(i+2), cr2, 1, X% u4(i+1)% p, 1)   !         X{i+1} = reshape(cr2, r, n(i+1), k(i+1), rx(i+2));
!!
!!        rx(i+1)= r                                                      !         rx(i+1) = r;
!!
!!        nrms(i+dir)= nrms(i)                                            !         nrms(i+dir) = nrms(i); % this is the norm of my block, save it
!!        !% Reduce
!!        !  [XAY(i+1,:,:),nrms(i)] = leftreduce_matrix(XAY(i,:,:), X{i}, A(i,:),
!!        !                Y(i,:), rx(i),n(i),k(i),rx(i+1), Ra,ra(i,:),ra(i+1,:),
!!        !                                            Ry,ry(i,:),m(i),ry(i+1,:));
!!        XAY(i+1)%arr= leftreduce_matrix(nrms(i), XAY(i)%arr, X%u4(i)%p, A2(i)%arr, &
!!                      Y%u4(i)%p, rx(i),n(i),k(i),rx(i+1),ra(i),ra(i+1),ry(i),m(i),ry(i+1))
!!        !% Enrichment
!!        if (kickrank > 0) then                                          !         if (kickrank>0)
!!           call qr(cr1, RR, crz); crz= cr1                              !             [crz,~]=qr(crz, 0);
!!           rz(i+1)= size(crz, 2)                                        !             rz(i+1) = size(crz, 2);
!!           ! ZAY(i+1,:,:) = leftreduce_matrix(ZAY(i,:,:), crz, A(i,:),
!!           !                Y(i,:), rz(i),n(i),k(i),rz(i+1), Ra,ra(i,:),ra(i+1,:),
!!           !                Ry,ry(i,:),m(i),ry(i+1,:), nrms(i));
!!           ZAY(i+1)%arr= leftreduce_matrix(aux, ZAY(i)%arr, &
!!                         reshape(crz,[size(crz),1,1,1]), A2(i)%arr, &
!!                         Y%u4(i)%p, rz(i),n(i),k(i),rz(i+1), ra(i),ra(i+1), &
!!                         ry(i),m(i),ry(i+1), nrms(i))
!!           ! ZX(i+1) = leftreduce_vector(ZX(i), crz, X(i),
!!           !           rz(i),n(i),k(i),rz(i+1), 1,rx(i),rx(i+1));
!!           ZX(i+1)%arr= leftreduce_vector(ZX(i)%arr, crz, X%u4(i)%p, &
!!                        rz(i),n(i),k(i),rz(i+1), rx(i),rx(i+1))
!!        endif                                                           !         end;
!!
!!     else if (dir < 0 .and. i > 1) then                                 !     elseif ((dir<0)&&(i>1))
!!        cr= reshape(cr, [rx(i), n(i)*k(i)*rx(i+1)])                     !         cr = reshape(cr, rx(i), n(i)*k(i)*rx(i+1));
!!        call svd(cr,up,sp,vp)                                           !         [u,s,v]=svd(cr, 'econ');
!!        s= sp                                                           !         s = diag(s);
!!        r= my_chop2(s, tol*norm(s)/sqrt(dble(d)))                       !         r = my_chop2(s, tol*norm(s)/sqrt(d));
!!        v= transpose(vp(1:r,:))                                         !         v = conj(v(:,1:r));
!!        u= thor_matmul(up(:,1:r),diag(s(1:r)))                          !         u = u(:,1:r)*diag(s(1:r));
!!
!!        !% Prepare enrichment, if needed
!!        if (kickrank > 0) then                                          !         if (kickrank>0)
!!           cr= thor_matmul(u, v, tsp_='nt')                             !             cr = u*v.';
!!           cr= reshape(cr, [rx(i), n(i)*k(i)*rx(i+1)])                  !             cr = reshape(cr, rx(i), n(i)*k(i)*rx(i+1));
!!           !% For updating z
!!           !  crz = local_matvec(Y(i,:), ry(i,:),m(i),k(i),ry(i+1,:),
!!           !        rz(i),n(i),rz(i+1), ZAY(i,:,:), A(i,:), ZAY(i+1,:,:),
!!           !        ra(i,:),ra(i+1,:));
!!           crz= local_matvec(Y%u4(i)%p, ry(i),m(i),k(i),ry(i+1), &
!!                rz(i),n(i),rz(i+1), ZAY(i)%arr, A2(i)%arr, ZAY(i+1)%arr, &
!!                ra(i),ra(i+1))
!!           crz= reshape(crz, [rz(i), n(i)*k(i)*rz(i+1)])                !             crz = reshape(crz, rz(i), n(i)*k(i)*rz(i+1));
!!           ys= thor_matmul(ZX(i)%arr, cr)                               !             ys = ZX{i}*cr;
!!           yz= reshape(ys, [rz(i)*n(i)*k(i), rx(i+1)])                  !             yz = reshape(ys, rz(i)*n(i)*k(i), rx(i+1));
!!           yz= thor_matmul(yz, ZX(i+1)%arr)                             !             yz = yz*ZX{i+1};
!!           yz= reshape(yz, [rz(i), n(i)*k(i)*rz(i+1)])                  !             yz = reshape(yz, rz(i), n(i)*k(i)*rz(i+1));
!!           crz= crz/nrms(i) - yz                                        !             crz = crz/nrms(i) - yz;
!!           nrmz= norm(crz)                                              !             nrmz = norm(crz,'fro');
!!           !% For adding into solution
!! !!call system_timer_stop
!! !!dt5 = system_dt
!! !!call system_timer_start
!!           !  crs = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:),
!!           !        rz(i),n(i),rx(i+1), ZAY(i,:,:), A(i,:), XAY(i+1,:,:),
!!           !        Ra,ra(i,:),ra(i+1,:));
!!           crs= local_matvec(Y%u4(i)%p, ry(i),m(i),k(i),ry(i+1), &
!!                rz(i),n(i),rx(i+1), ZAY(i)%arr, A2(i)%arr, XAY(i+1)%arr, &
!!                ra(i),ra(i+1))
!!           crs= reshape(crs, [rz(i), n(i)*k(i)*rx(i+1)])                !             crs = reshape(crs, rz(i), n(i)*k(i)*rx(i+1));
!!           crs= crs/nrms(i) - ys                                        !             crs = crs/nrms(i) - ys;
!!           v= stack(v, transpose(crs))                                  !             v = [v,crs.'];
!!           call qr(cr1, RR, v); v= cr1                                  !             [v,R]=qr(v, 0);
!!           u= stack(u, zeros(rx(i), rz(i)))                             !             u = [u, zeros(rx(i), rz(i))];
!!           u= thor_matmul(u, RR, tsp_='nt')                             !             u = u*R.';
!!           r= size(v, 2)                                                !             r = size(v, 2);
!! !!call system_timer_stop
!! !!dt6 = system_dt
!! !!call system_timer_start
!!        endif                                                           !         end;
!!        cr2= reshape(X% u4(i-1)% p, [rx(i-1)*n(i-1)*k(i-1), rx(i)])     !         cr2 = reshape(X{i-1}, rx(i-1)*n(i-1)*k(i-1), rx(i));
!!        cr2= thor_matmul(cr2, u)                                        !         cr2 = cr2*u;
!!        deallocate(X% u4(i-1)% p)
!!        allocate(X% u4(i-1)% p(rx(i-1),n(i-1),k(i-1),r))
!!        call dcopy(rx(i-1)*n(i-1)*k(i-1)*r, cr2, 1, X% u4(i-1)% p, 1)   !         X{i-1} = reshape(cr2, rx(i-1), n(i-1), k(i-1), r);
!!        v= transpose(v)
!!        deallocate(X% u4(i)% p)
!!        allocate(X% u4(i)% p(r,n(i),k(i),rx(i+1)))
!!        call dcopy(r*n(i)*k(i)*rx(i+1), v, 1, X% u4(i)% p, 1)           !         X{i} = reshape(v.', r, n(i), k(i), rx(i+1));
!!
!!        rx(i)= r                                                        !         rx(i) = r;
!!
!!        nrms(i+dir) = nrms(i) !% this is the norm of my block, save it  !         nrms(i+dir) = nrms(i); % this is the norm of my block, save it
!!        !% Reduce
!!        !  [XAY(i,:,:),nrms(i)] = rightreduce_matrix(XAY(i+1,:,:), X{i}, A(i,:),
!!        !           Y(i,:), rx(i),n(i),k(i),rx(i+1), Ra,ra(i,:),ra(i+1,:),
!!        !           Ry,ry(i,:),m(i),ry(i+1,:));
!!        XAY(i)%arr= rightreduce_matrix(nrms(i), XAY(i+1)%arr, X%u4(i)%p, &
!!                    A2(i)%arr, Y%u4(i)%p, rx(i),n(i),k(i),rx(i+1), &
!!                    ra(i),ra(i+1), ry(i),m(i),ry(i+1))
!! !!call system_timer_stop
!! !!dt7 = system_dt
!! !!call system_timer_start
!!        !% Enrich
!!        if (kickrank > 0) then                                          !         if (kickrank>0)
!!           call qr(cr1, RR, transpose(crz)); crz= cr1                   !             [crz,~]=qr(crz.', 0);
!!           rz(i)= size(crz, 2)                                          !             rz(i) = size(crz, 2);
!!           !  ZAY(i,:,:) = rightreduce_matrix(ZAY(i+1,:,:), crz, A(i,:),
!!           !     Y(i,:), rz(i),n(i),k(i),rz(i+1), Ra,ra(i,:),ra(i+1,:),
!!           !     Ry,ry(i,:),m(i),ry(i+1,:), nrms(i));
!!           ZAY(i)%arr= rightreduce_matrix(aux, ZAY(i+1)%arr, &
!!                       reshape(crz, [size(crz),1,1,1]), A2(i)%arr, &
!!                       Y%u4(i)%p, rz(i),n(i),k(i),rz(i+1), ra(i),ra(i+1), &
!!                       ry(i),m(i),ry(i+1), nrms(i))
!!           !  ZX(i) = rightreduce_vector(ZX(i+1), crz, X(i),
!!           !          rz(i),n(i),k(i),rz(i+1), 1,rx(i),rx(i+1));
!!           ZX(i)%arr= rightreduce_vector(ZX(i+1)%arr, crz, X%u4(i)%p, &
!!                      rz(i),n(i),k(i),rz(i+1), rx(i),rx(i+1))
!!        endif                                                           !         end;
!! !!call system_timer_stop
!! !!dt8 = system_dt
!! !!call system_timer_start
!!     else                                                               !     else
!!           deallocate(X% u4(i)% p)
!!           allocate(X% u4(i)% p(rx(i),n(i),k(i),rx(i+1)))
!!           call dcopy(rx(i)*n(i)*k(i)*rx(i+1), cr, 1, X%u4(i)%p, 1)     !         X{i} = reshape(cr, rx(i), n(i), k(i), rx(i+1));
!!
!!     endif                                                              !     end;
!!
!!     if (verb > 1) then                                                 !     if (verb>1)
!!        !fprintf('amen-mm: swp=[%d,%d], dx=%3.3e, r=%d, |X|=%3.3e, |z|=%3.3e\n', swp, i, dx, r, norm(cr,'fro'), nrmz);
!!        print '("amen-mm: swp=["I3","I3"], dx="ES14.7", r="I5", '// &
!!              '|X|="ES14.7", |z|="ES14.7)', swp, i, dx, r, norm(cr), nrmz
!!     endif                                                              !     end;
!!
!!     !% Stopping or reversing
!!     if ((dir > 0.and.i == d) .or. (dir < 0.and.i == 1)) then           !     if ((dir>0)&&(i==d))||((dir<0)&&(i==1))
!!        if (verb > 0) then                                              !         if (verb>0)
!!           ! fprintf('amen-mm: swp=%d{%d}, max_dx=%3.3e, max_r=%d\n', swp, (1-dir)/2, max_dx, max(rx));
!!           print '("amen-mm: swp="I3"{"I2"}, max_dx="ES14.7", max_r="I5)', &
!!                 swp, (1-dir)/2, max_dx, maxval(rx)
!!        endif                                                           !         end;
!!        if ((max_dx < tol .or. swp==nswp).and. dir > 0) then            !         if ((max_dx<tol)||(swp==nswp))&&(dir>0)
!!           exit sweeps                                                  !             break;
!!        endif                                                           !         end;
!!        swp= swp + 1                                                    !         swp = swp+1;
!!        max_dx= 0d0                                                     !         max_dx = 0;
!!        dir= -dir                                                       !         dir = -dir;
!!     else                                                               !     else
!!        i= i + dir                                                      !         i = i+dir;
!!     endif                                                              !     end;
!!  enddo sweeps                                                          ! end;
!! !!call system_timer_stop
!! !!dt9 = system_dt
!! !!call system_timer_start
!!
!!  !% Distribute norms equally...
!!  nrmz= dexp(sum(dlog(nrms))/dble(d))                                   ! nrms = exp(sum(log(nrms))/d);
!!  !% ... and plug them into y
!!  X% r(0:d)= rx(1:d+1)
!!  X% n(1:d)= X% q(1:d) * X% s(1:d)
!!  do i=1,d                                                              ! for i=1:d
!!     X% u4(i)% p= X% u4(i)% p * nrmz                                    !     X{i} = X{i}*nrms;
!!     !!! TODO: make sure other fields of X are correctly filled out!!!
!!     X% u(i)% p(1:rx(i),1:X%n(i),1:rx(i+1))=> X% u4(i)% p
!!  enddo                                                                 ! end;
!!
!!                                                                        ! % Return the correct form
!!                                                                        ! if (outtype==0)
!!                                                                        !     for i=1:d
!!                                                                        !         X{i} = reshape(X{i}, rx(i), n(i)*k(i), rx(i+1));
!!                                                                        !     end;
!!                                                                        !     X = cell2core(tt_tensor, X);
!!                                                                        ! elseif (outtype==1)
!!                                                                        !     X = cell2core(tt_matrix, X);
!!                                                                        ! end;
!! !!call system_timer_stop
!! !!dt10 = system_dt
!! !!!print '(10(F9.5,1X))', dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10
!! !!print '(10(F9.5,1X))', dt1,dt2,dt3,dt3a,dt3b,dt3c,dt3d
!!  end subroutine amen_mm
!!
!!
!!  !% Accumulates the left reduction W{1:k}'*A{1:k}*X{1:k}
!!  !function [WAX2,nrm] =
!!  !             leftreduce_matrix(WAX2, w, A,  x, rw1, n, k, rw2, Ra, ra1, ra2, Rx, rx1, m, rx2, extnrm)
!!  function leftreduce_matrix(nrm, WAX1, w, A2, x, rw1, n, k, rw2,     ra1, ra2,     rx1, m, rx2, extnrm_) result(WAX2)
!!  use matrix_util, only: perm3d
!!  use mat_lib, only: norm=> normfro, thor_matmul_mnk
!!  use time_lib
!!  implicit none
!!  double precision, allocatable, target   :: WAX2(:,:), aux(:,:)
!!  double precision, intent(INOUT) :: nrm
!!  double precision, intent(IN)    :: WAX1(:,:), w(:,:,:,:), A2(:,:), x(:,:,:,:)
!!  integer, intent(IN)             :: rw1, n, k, rw2, m, ra1, ra2, rx1, rx2
!!  double precision, intent(IN), optional :: extnrm_
!!  !
!!  double precision, allocatable :: w2(:,:), w3(:,:,:), xc(:,:,:), xc2(:,:)
!!  double precision, pointer :: wp(:)
!!
!! !!double precision:: dt1,dt2,dt3,dt4,dt_mtm,dt_tsp,dt_prm
!! !!integer :: ii,iimax
!!
!!     !% Left WAX has the form of the first matrix TT block, i.e. [rw, rx, ra]
!!                                                                        !if (nargin<16)
!!                                                                        !    extnrm = [];
!!                                                                        !end;
!!     nrm= 0d0                                                           !if (nargout>1)
!!                                                                        !    nrm = 0;
!!                                                                        !end;
!!                                                                        !
!! !!dt_mtm=0d0;dt_tsp=0d0;dt_prm=0d0
!! !!iimax = 1
!! !!do ii=1,iimax
!! !!call system_timer_start
!! !!!DEBUG!print '(A,3(I5,1X))', "rw1*n, k, rw2 = ", rw1*n, k, rw2
!!     w3 = reshape(w, [rw1*n, k, rw2])                                   !w = reshape(w, rw1*n, k, rw2);
!!     call perm3d(w3, 132)                                               !w = permute(w, [1,3,2]);
!!     w2 = reshape(w3, [rw1, n*rw2*k])                                   !w = reshape(w, rw1, n*rw2*k);
!!
!!     ! Rx == 1                                                          !for j=1:Rx
!!     xc = reshape(x, [rx1*m, k, rx2])                                   !    xc = reshape(x{j}, rx1*m, k, rx2);
!!     call perm3d(xc, 321)                                               !    xc = permute(xc, [3,2,1]);
!!     xc2 = reshape(xc, [rx2, k*rx1*m])                                  !    xc = reshape(xc, rx2, k*rx1*m);
!! !!call system_timer_stop
!! !!dt_prm= dt_prm + system_dt
!! !!!print '("MTM: w^T*WAX2: ["I8","I8"]x["I8","I8"]")', n*rw2*k, rw1, rw1, rx1*ra1
!! !!call system_timer_start
!!     ! Ra == 1                                                          !    for i=1:Ra
!! !!    WAX2= reshape(WAX1, [rw1, rx1*ra1])                                !        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rw1, rx1*ra1);
!!     WAX2= thor_matmul_mnk(n*rw2*k,rx1*ra1,rw1,w2,WAX1,tsp_='tn')       !        WAX2{1,i,j} = w'*WAX2{1,i,j}; % size n rw2 x rx1 ra1
!! !!    WAX2= aux
!! !!call system_timer_stop
!! !!dt_mtm= dt_mtm + system_dt
!! !print '("TSP: WAX2: ["I8","I8"]")', n, rw2*k*rx1*ra1
!! !!call system_timer_start
!!     WAX2= transpose(reshape(WAX2, [n, rw2*k*rx1*ra1]))                 !        WAX2{1,i,j} = reshape(WAX2{1,i,j}, n, rw2*k*rx1*ra1);
!! !!call system_timer_stop
!! !!dt_tsp= dt_tsp + system_dt
!! !!call system_timer_start
!! !!!wp(1:n * rw2*k*rx1 * ra1)=> WAX2
!! !!!call perm3d(wp,n,rw2*k*rx1,ra1, 132)
!! !!!aux= thor_matmul_mnk(rw2*k*rx1, m*ra2, n*ra1, wp, A2, tsp_='tn')
!!     !WAX2= transpose(WAX2)                                              !        WAX2{1,i,j} = WAX2{1,i,j}.';
!! !!!!    WAX2= reshape(WAX2, [rw2*k*rx1, ra1*n])                            !        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rw2*k*rx1, ra1*n);
!!     WAX2= thor_matmul_mnk(rw2*k*rx1, m*ra2, n*ra1, WAX2, A2)                                         !        WAX2{1,i,j} = WAX2{1,i,j}*A{i}; % size rw2 rx1 m ra2
!!     !WAX2= aux
!! !!call system_timer_stop
!! !!dt_mtm= dt_mtm + system_dt
!! !print '("TSP: WAX2: ["I8","I8"]")', rw2, k*rx1*m*ra2
!! !!call system_timer_start
!!     WAX2= transpose(reshape(WAX2, [rw2, k*rx1*m*ra2]))                 !        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rw2, k*rx1*m*ra2);
!! !!call system_timer_stop
!! !!dt_tsp= dt_tsp + system_dt
!! !!call system_timer_start
!!     !WAX2= transpose(WAX2)                                              !        WAX2{1,i,j} = WAX2{1,i,j}.';
!! !!!!    WAX2= reshape(WAX2, [k*rx1*m, ra2*rw2])                            !        WAX2{1,i,j} = reshape(WAX2{1,i,j}, k*rx1*m, ra2*rw2);
!!     WAX2= thor_matmul_mnk(rx2, ra2*rw2, k*rx1*m, xc2, WAX2)                                        !        WAX2{1,i,j} = xc*WAX2{1,i,j}; % size rx2, ra2 rw2
!!     !WAX2= aux
!! !!call system_timer_stop
!! !!dt_mtm= dt_mtm + system_dt
!! !print '("TSP: WAX2: ["I8","I8"]")', rx2*ra2, rw2
!! !!call system_timer_start
!!     WAX2= transpose(reshape(WAX2, [rx2*ra2, rw2]))                     !        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rx2*ra2, rw2);
!! !!call system_timer_stop
!! !!dt_tsp= dt_tsp + system_dt
!!     !WAX2= transpose(WAX2)                                              !        WAX2{1,i,j} = WAX2{1,i,j}.';
!! !!enddo
!! !DEBUG!call system_timer_stop
!! !DEBUG!dt4= system_dt
!! !!print '("leftreduce_matrix timings: ",12(F9.5))', &
!! !!      dt_prm/iimax,dt_mtm/iimax,dt_tsp/iimax
!!                                                                        !%         WAX2{i,j} = reshape(WAX2{i,j}, rw2, rx2(j), ra2(i));
!!     nrm= max(nrm, norm(WAX2))                                          !              if (nargout>1)
!!                                                                        !            nrm = max(nrm, norm(WAX2{1,i,j}, 'fro'));
!!                                                                        !        end;
!!                                                                        !    end;
!!                                                                        !end;
!!                                                                        !if (nargout>1)
!!     !% Extract the scale to prevent overload
!!     if (.not.present(extnrm_)) then                                    !    if (nrm>0)
!!        if (nrm > 0d0) then                                             !        for i=1:Rx*Ra
!!                                                                        !            WAX2{i} = WAX2{i}/nrm;
!!           WAX2= WAX2/nrm                                               !        end;
!!
!!        else                                                            !    else
!!           nrm = 1d0                                                    !        nrm=1;
!!        endif                                                           !    end;
!!     else                                                               !elseif (~isempty(extnrm))
!!        !% Override the normalization
!!                                                                        !    for i=1:Rx*Ra
!!        WAX2= WAX2/extnrm_                                              !        WAX2{i} = WAX2{i}/extnrm;
!!                                                                        !    end;
!!     endif                                                              !end;
!!  end function leftreduce_matrix                                        !end
!!
!!
!!  !% Accumulates the right reduction W{k:d}'*A{k:d}*X{k:d}
!!  !function [WAX1,nrm] =
!!  !             rightreduce_matrix(WAX1, w, A,  x, rw1, n, k, rw2, Ra,ra1, ra2, Rx,rx1, m, rx2, extnrm)
!!  function rightreduce_matrix(nrm, WAX2, w, A2, x, rw1, n, k, rw2,    ra1, ra2,    rx1, m, rx2, extnrm_) result(WAX1)
!!  use matrix_util, only: perm3d
!!  use mat_lib, only: norm=> normfro, thor_matmul_mnk
!!  implicit none
!!  double precision, allocatable   :: WAX1(:,:)
!!  double precision, intent(INOUT) :: nrm
!!  double precision, intent(IN)    :: WAX2(:,:), w(:,:,:,:), A2(:,:), x(:,:,:,:)
!!  integer, intent(IN)             :: rw1, n, k, rw2, m, ra1, ra2, rx1, rx2
!!  double precision, intent(IN), optional :: extnrm_
!!  !
!!  double precision, allocatable :: w2(:,:), w3(:,:,:), xc(:,:,:), xc2(:,:)
!!
!!     !% Right WAX has the form of the last matrix TT block, i.e. [ra, rw, rx]
!!                                                                        !if (nargin<16)
!!                                                                        !    extnrm = [];
!!                                                                        !end;
!!     nrm= 0d0                                                           !if (nargout>1)
!!                                                                        !    nrm = 0;
!!                                                                        !end;
!!                                                                        !
!!     w3 = reshape(w, [rw1*n, k, rw2])                                   !w = reshape(w, rw1*n, k, rw2);
!!     call perm3d(w3, 132)                                               !w = permute(w, [1,3,2]);
!!     w2 = reshape(w3, [rw1, n*rw2*k])                                   !w = reshape(w, rw1, n*rw2*k);
!!                                                                        !w = conj(w);
!!     ! j == 1                                                           !for j=1:Rx
!!     xc = reshape(x, [rx1*m, k, rx2])                                   !    xc = reshape(x{j}, rx1(j)*m, k, rx2(j));
!!     call perm3d(xc, 213)                                               !    xc = permute(xc, [2,1,3]);
!!     !xc2 = reshape(xc, [k*rx1*m, rx2])                                  !    xc = reshape(xc, k*rx1(j)*m, rx2(j));
!!     ! i == 1                                                           !    for i=1:Ra
!!     !WAX1= reshape(WAX2, [ra2*rw2, rx2])                                !      WAX1{1,i,j} = reshape(WAX1{1,i,j}, ra2(i)*rw2, rx2(j));
!!     !WAX1= thor_matmul(xc2, WAX1, tsp_='nt')                            !      WAX1{1,i,j} = xc*WAX1{1,i,j}.'; % size rx1 m x ra2 rw2
!!     WAX1= thor_matmul_mnk(k*rx1*m, ra2*rw2, rx2, xc, WAX2, tsp_='nt')                            !      WAX1{1,i,j} = xc*WAX1{1,i,j}.'; % size rx1 m x ra2 rw2
!!     WAX1= transpose(reshape(WAX1, [k*rx1, m*ra2*rw2]))                 !      WAX1{1,i,j} = reshape(WAX1{1,i,j}, k*rx1(j), m*ra2(i)*rw2);
!!     !WAX1= transpose(WAX1)                                              !      WAX1{1,i,j} = WAX1{1,i,j}.';
!!     !WAX1= reshape(WAX1, [m*ra2, rw2*k*rx1])                            !      WAX1{1,i,j} = reshape(WAX1{1,i,j}, m*ra2(i), rw2*k*rx1(j));
!!     !WAX1= thor_matmul(A2, WAX1)                                        !      WAX1{1,i,j} = A{i}*WAX1{1,i,j}; % size ra1(k)*n, rw2*rx1
!!     WAX1= thor_matmul_mnk(ra1*n, rw2*k*rx1, m*ra2, A2, WAX1)                                        !      WAX1{1,i,j} = A{i}*WAX1{1,i,j}; % size ra1(k)*n, rw2*rx1
!!     WAX1= transpose(reshape(WAX1, [ra1, n*rw2*k*rx1]))                 !      WAX1{1,i,j} = reshape(WAX1{1,i,j}, ra1(i), n*rw2*k*rx1(j));
!!     !WAX1= transpose(WAX1)                                              !      WAX1{1,i,j} = WAX1{1,i,j}.';
!!     !WAX1= reshape(WAX1, [n*rw2*k, rx1*ra1])                            !      WAX1{1,i,j} = reshape(WAX1{1,i,j}, n*rw2*k, rx1(j)*ra1(i));
!!     !WAX1= thor_matmul(w2, WAX1)                                        !      WAX1{1,i,j} = w*WAX1{1,i,j}; % size rw1, rx1 ra1
!!     WAX1= thor_matmul_mnk(rw1, rx1*ra1, n*rw2*k, w3, WAX1)             !      WAX1{1,i,j} = w*WAX1{1,i,j}; % size rw1, rx1 ra1
!!     WAX1= transpose(reshape(WAX1, [rw1*rx1, ra1]))                     !      WAX1{1,i,j} = reshape(WAX1{1,i,j}, rw1*rx1(j), ra1(i));
!!     !WAX1= transpose(WAX1)                                              !      WAX1{1,i,j} = WAX1{1,i,j}.';
!!                                                                        !%         WAX1{i,j} = reshape(WAX1{i,j}, ra1(i), rw1, rx1(j));
!!     nrm= max(nrm, norm(WAX1))                                          !        if (nargout>1)
!!                                                                        !            nrm = max(nrm, norm(WAX1{1,i,j}, 'fro'));
!!                                                                        !        end;
!!                                                                        !    end;
!!                                                                        !end;
!!                                                                        !if (nargout>1)
!!     !% Extract the scale to prevent overload
!!     if (.not.present(extnrm_)) then                                    !    if (nrm>0)
!!        if (nrm > 0d0) then                                             !        for i=1:Rx*Ra
!!                                                                        !            WAX1{i} = WAX1{i}/nrm;
!!            WAX1= WAX1/nrm                                              !        end;
!!
!!        else                                                            !    else
!!           nrm = 1d0                                                    !        nrm=1;
!!        endif                                                           !    end;
!!     else                                                               !elseif (~isempty(extnrm))
!!        !% Override the normalization
!!                                                                        !    for i=1:Rx*Ra
!!        WAX1= WAX1/extnrm_                                              !        WAX1{i} = WAX1{i}/extnrm;
!!                                                                        !    end;
!!     endif                                                              !end;
!!  end function rightreduce_matrix                                       !end
!!
!!
!!  !% Accumulates the left reduction W{1:k}'*X{1:k}
!!  !% function [WX2] = leftreduce_vector(WX1, w, x, rw1,n,k,rw2, Rx,rx1,rx2, extnrm)
!!  function leftreduce_vector(WX1,w,x,rw1,n,k,rw2,rx1,rx2,extnrm_) result(WX2)
!!  use mat_lib, only: thor_matmul, thor_matmul_mnk
!!  implicit none
!!  double precision, allocatable :: WX2(:,:)
!!  double precision, intent(IN)  :: WX1(:,:),w(:,:),x(:,:,:,:)
!!  integer, intent(IN) :: rw1, n, k, rw2, rx1, rx2
!!  double precision, intent(IN), optional :: extnrm_
!!  !
!!  double precision, allocatable :: w2(:,:), tmp(:,:)
!!
!!                                                                        !if (nargin<11)
!!                                                                        !    extnrm = [];
!!                                                                        !end;
!!
!!     !% Left WX has the form of the first vector TT block, i.e. [rw, rx]
!!     !WX2= WX1                                                           !WX2 = WX1;
!!     !w2= reshape(w, [rw1, n*k*rw2])                                     !wc = reshape(w, rw1, n*k*rw2);
!!     ! i == 1                                                           !for i=1:Rx
!!     !WX2= thor_matmul(w2,WX2,tsp_='tn')                                 !    WX2{1,i} = wc'*WX2{1,i}; % size n rw2 x rx1
!!     WX2= thor_matmul_mnk(n*k*rw2,rx1,rw1,w,WX1,tsp_='tn')              !    WX2{1,i} = wc'*WX2{1,i}; % size n rw2 x rx1
!!     WX2= transpose(reshape(WX2, [n*k, rw2*rx1]))                       !    WX2{1,i} = reshape(WX2{1,i}, n*k, rw2*rx1(i));
!!     !WX2= transpose(WX2)                                                !    WX2{1,i} = WX2{1,i}.';
!!     !WX2= reshape(WX2, [rw2, rx1*n*k])                                  !    WX2{1,i} = reshape(WX2{1,i}, rw2, rx1(i)*n*k);
!!     !tmp= reshape(x, [rx1*n*k, rx2])                                    !    tmp = reshape(x{i}, rx1(i)*n*k, rx2(i));
!!     !WX2= thor_matmul(WX2, tmp)                                         !    WX2{1,i} = WX2{1,i}*tmp; % size rw2, rx2
!!     WX2= thor_matmul_mnk(rw2, rx2, rx1*n*k, WX2, x)                    !    WX2{1,i} = WX2{1,i}*tmp; % size rw2, rx2
!!                                                                        !end;
!!     if (present(extnrm_)) then                                         !if (~isempty(extnrm))
!!        !% Override the normalization
!!        WX2= WX2/extnrm_                                                !    for i=1:Rx
!!                                                                        !        WX2{i} = WX2{i}/extnrm;
!!     endif                                                              !    end;
!!  end function leftreduce_vector                                        !end;
!!
!!
!!  !% Accumulates the right reduction W{k:d}'*X{k:d}
!!  !% function [WX1] = rightreduce_vector(WX2, w, x, rw1,n,k,rw2, Rx,rx1,rx2, extnrm)
!!  function rightreduce_vector(WX2,w,x,rw1,n,k,rw2,rx1,rx2,extnrm_) result(WX1)
!!  use mat_lib, only: thor_matmul_mnk
!!  implicit none
!!  double precision, allocatable :: WX1(:,:)
!!  double precision, intent(IN)  :: WX2(:,:),w(:,:),x(:,:,:,:)
!!  integer, intent(IN) :: rw1, n, k, rw2, rx1, rx2
!!  double precision, intent(IN), optional :: extnrm_
!!  !
!!  double precision, allocatable :: w2(:,:), tmp(:,:)
!!
!!                                                                        !if (nargin<11)
!!                                                                        !    extnrm = [];
!!                                                                        !end;
!!                                                                        !
!!     !% Right WX has the form of the last vector TT block, i.e. [rx, rw]
!!     !WX1= WX2                                                           !WX1 = WX2;
!!     !w2 = reshape(w, [rw1, n*k*rw2])                                    !wc = reshape(w, rw1, n*k*rw2);
!!     ! i == 1                                                           !for i=1:Rx
!!     !tmp= reshape(x, [rx1*n*k, rx2])                                    !    tmp = reshape(x{i}, rx1(i)*n*k, rx2(i));
!!     !WX1= thor_matmul(tmp, WX2)                                         !    WX1{1,i} = tmp*WX1{1,i}; % size rx1 n x rw2
!!     WX1= thor_matmul_mnk(rx1*n*k, rw2, rx2, x, WX2)                    !    WX1{1,i} = tmp*WX1{1,i}; % size rx1 n x rw2
!!     !WX1= reshape(WX1, [rx1, n*k*rw2])                                  !    WX1{1,i} = reshape(WX1{1,i}, rx1(i), n*k*rw2);
!!     WX1= thor_matmul_mnk(rx1, rw1, n*k*rw2, WX1, w, tsp_='nt')                               !    WX1{1,i} = WX1{1,i}*wc'; % size rx1, rw1
!!                                                                        !end;
!!     if (present(extnrm_)) then                                         !if (~isempty(extnrm))
!!        !    % Override the normalization
!!                                                                        !    for i=1:Rx
!!        WX1= WX1/extnrm_                                                !        WX1{i} = WX1{i}/extnrm;
!!                                                                        !    end;
!!     endif                                                              !end;
!!  end function rightreduce_vector
!!
!!
!!  !% A matrix-matrix product for the matrix in the 3D TT (WAX1-A-WAX2), and
!!  !% full matrix of size (rx1*m*k*rx2). Returns (rw1*n*k*rw2)
!!  !function [w]=local_matvec(x, Rx,rx1,m,k,rx2, rw1,n,rw2, WAX1, A, WAX2, Ra,ra1,ra2)
!!  function local_matvec(x,         rx1,m,k,rx2, rw1,n,rw2, WAX1, A, WAX2,    ra1,ra2) result(w)
!!  use mat_lib, only: thor_matmul, thor_matmul_mnk
!!  use matrix_util, only: perm3d
!!  implicit none
!!  double precision, allocatable :: w(:,:)
!!  double precision, intent(IN) :: x(:,:,:,:), WAX1(:,:), A(:,:), WAX2(:,:)
!!  integer, intent(IN) :: rx1, m, k, rx2, rw1, n, rw2, ra1, ra2
!!  !
!!  double precision, allocatable :: w2(:,:), w3(:,:,:), xc(:,:,:), xc2(:,:), tmp(:,:), wk(:,:)
!!
!!     allocate(w2(rw1*n*rw2, k)); w2= 0d0                                !w = zeros(rw1*n*rw2,k);
!!     ! j == 1 always, because Rx == 1                                   !for j=1:Rx
!!     xc= reshape(x, [rx1*m, k, rx2])                                    !    xc = reshape(x{j}, rx1(j)*m, k, rx2(j));
!!     call perm3d(xc, 213)                                               !    xc = permute(xc, [2,1,3]);
!!     !xc2= reshape(xc, [k*rx1*m, rx2])                                   !    xc = reshape(xc, k*rx1(j)*m, rx2(j));
!!     !! i == 1 always                                                    !    for i=1:Ra
!!     !tmp= reshape(WAX2, [ra2*rw2, rx2])                                 !        tmp = reshape(WAX2{1,i,j}, ra2(i)*rw2, rx2(j));
!!     !wk= thor_matmul(xc2, tmp, tsp_='nt')                               !        wk = xc*tmp.';
!!     wk= thor_matmul_mnk(k*rx1*m, ra2*rw2, rx2, xc, WAX2, tsp_='nt')    !        wk = xc*tmp.';
!!     wk= transpose(reshape(wk, [k*rx1, m*ra2*rw2]))                     !        wk = reshape(wk, k*rx1(j), m*ra2(i)*rw2);
!!     !wk= reshape(wk, [k*rx1, m*ra2*rw2])                                !        wk = reshape(wk, k*rx1(j), m*ra2(i)*rw2);
!!     !wk= transpose(wk)                                                  !        wk = wk.';
!!     !wk= reshape(wk, [m*ra2, rw2*k*rx1])                                !        wk = reshape(wk, m*ra2(i), rw2*k*rx1(j));
!!     wk= thor_matmul_mnk(ra1*n, rw2*k*rx1, m*ra2, A, wk)                !        wk = A{i}*wk;
!!     wk= transpose(reshape(wk, [ra1*n*rw2*k, rx1]))                     !        wk = reshape(wk, ra1(i)*n*rw2*k, rx1(j));
!!     !wk= transpose(wk)                                                  !        wk = wk.';
!!     !wk= reshape(wk, [rx1*ra1, n*rw2*k])                                !        wk = reshape(wk, rx1(j)*ra1(i), n*rw2*k);
!!     !tmp= reshape(WAX1, [rw1, rx1*ra1])                                 !        tmp = reshape(WAX1{1,i,j}, rw1, rx1(j)*ra1(i));
!!     !wk= thor_matmul(tmp, wk)                                           !        wk = tmp*wk;
!!     wk= thor_matmul_mnk(rw1, n*rw2*k, rx1*ra1, WAX1, wk)               !        wk = tmp*wk;
!!     wk= reshape(wk, [rw1*n*rw2, k])                                    !        wk = reshape(wk, rw1*n*rw2, k);
!!     w2= w2 + wk                                                        !        w = w+wk;
!!                                                                        !    end;
!!                                                                        !end;
!!     w3= reshape(w2, [rw1*n, rw2, k])                                   !w = reshape(w, rw1*n, rw2, k);
!!     call perm3d(w3, 132)                                               !w = permute(w, [1,3,2]);
!!     w= reshape(w3, [rw1*n*rw2*k, 1])                                   !w = reshape(w,[],1);
!!
!!  end function local_matvec                                             !end


end program
