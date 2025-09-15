!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file low_rank.f90
!! @author [IDB], [OGK]
!! @date  March 2025
!! @brief Testing how much a 2D matrix can be compressed
!!

program main
use rnd_lib, only: random, init_random_seed
implicit none
integer i, ierr
character(len=128) arg
!
integer :: d = 3 ! binary log of the tensor order
logical :: lsave = .false.
logical :: lbrief = .false.
double precision :: tol = 1d-12
character(100) :: tf = "circle"

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     elseif (trim(arg).eq.'-d') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=ierr) d
        if (ierr.ne.0) error stop "ERROR: when parsing -d <???>"
     elseif (trim(arg).eq.'-t') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=ierr) tol
        if (ierr.ne.0) error stop "ERROR: when parsing -t <???>"
     elseif (trim(arg).eq.'--save'.or.trim(arg).eq.'-s') then
        lsave= .true.
     elseif (trim(arg).eq.'--brief'.or.trim(arg).eq.'-b') then
        lbrief= .true.
     elseif (trim(arg).eq.'--diag') then
        tf= "diag"
     elseif (trim(arg).eq.'--line') then
        tf= "line"
     elseif (trim(arg).eq.'--radiant') then
        tf= "radiant"
     elseif (trim(arg).eq.'--clocks') then
        tf= "clocks"
     elseif (trim(arg).eq.'--patrick') then
        tf= "patrick"
     elseif (trim(arg).eq.'--permute') then
        tf= "permute"
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call test_low_rank(d, tol, lsave, .not.lbrief, tf)

contains

  subroutine help_msg()
  print '(11(/,A))', &
   "Testing how much a 2D matrix can be compressed", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>   : tensor order <ord> [3]", &
   " -t <tol>   : use this tolerance in SVD step [1e-12]", &
   " -b|--brief : only output the ranks [F]", &
   " -s|--save  : save generated tensors to /tmp [F]", &
   " --diag|--line|--circle|--radiant"// &
   "|--clocks|--patrick|--permute : test function [circle]"
  end subroutine

  subroutine test_low_rank(d, tol, lsave, verbose, test_function)
  use time_lib
  use qtt_lib
  use matlab_struct_module, only : pprint_matrix
  use thor_lib
  use matrix_util, only : d1copy
  use mat_lib, only: normfro
  implicit none
  integer, intent(in) :: d
  double precision, intent(in) :: tol
  logical, intent(in) :: lsave
  logical, intent(in) :: verbose
  character(100), intent(in) :: test_function
  !
  type(dtt_tensor) :: t1
  type(qtt_tensor) :: q1
  type(qtt_matrix) :: qm
  integer :: failed, i, j, MA, fd
  integer :: mr_tt, mr_qtt_dumb, mr_qtt2, mr_qtt4
  integer, allocatable :: nn(:)
  double precision :: dt, size_tt, size_dumb_qtt, size_qtt2, size_qtt4, err
  double precision, allocatable:: A(:), A2(:,:), AQ(:,:)

     ! initialize system timer
     call system_timer_init()

     if (verbose) print '("Ranks of 2D matrix in TT and QTT representation")'

     if (verbose) print '("==="/"-- initializing a matrix")'
     MA= 2**d
     select case(trim(test_function))
     case ("diag")
        A2= Aij_diag(d)
     case ("line")
        A2= Aij_line(d)
     case ("circle")
        A2= Aij_circle(d)
     case ("radiant")
        A2= Aij_radiant(d)
     case ("clocks")
        A2= Aij_clocks(d)
     case ("patrick")
        A2= Aij_patrick(d)
     case ("permute")
        A2= Aij_permute(d)
     end select

     if (verbose) call pprint_matrix(A2)

     if (lsave .and. d <= 12) then
        if (verbose) print '(/,"-- saving the matrix")'
        open (fd, file='/tmp/A2.dat', status='replace')
        write (fd, '("# matrix A2 from low_rank.exe")')
        do i=1,MA
           write (fd, '(4096(ES14.7,1X))') A2(i,:)
        enddo
        print '("file /tmp/A2.dat written")'
     endif
     if (lsave .and. d > 12) print '("WARNING: matrix too big, not saving")'

     !!if (verbose) print '(/,"-- computing 2D TT tensor")'
     !!t1 = dtt_tensor(A2, eps_=tol)
     !!size_tt= t1% size()
     !!mr_tt= maxval(t1% r(0:2))
     !!if (verbose) call t1% say()

     !if (verbose) print '(/,"-- computing with dumb QTT")'
     !allocate(A(MA*MA))
     !call d1copy(A2, A, MA*MA)
     !q1= qtt_tensor(A, eps_=tol)
     !size_dumb_qtt= q1% size()
     !if (verbose) call q1% say()

     if (verbose) print '(/,"-- computing with qtt")'
     qm= qtt_matrix(A2, eps_=tol)
     if (verbose) call qm% say()
     mr_qtt4= maxval(qm% r(0:d))
     size_qtt4= qm% size()

     if (verbose) print '(/,"-- recovering the full matrix")'
     AQ= qm% full_matrix()
     err= normfro(AQ - A2)

     if (verbose) print '(/,"-- d, maxranks:tt,qtt; sizes:full,qtt; error")'
     print '(I2,2(1X,I5),5(1X,ES14.7))', &
             d, mr_tt,mr_qtt4, dble(MA)**2,size_tt,size_qtt4, err

     if (lsave .and. d <= 12) then
        if (verbose) print '(/,"-- saving the qtt-approximated matrix")'
        open (fd, file='/tmp/A2qtt.dat', status='replace')
        write (fd, '("# matrix A2 after QTT, tol = ",ES14.7)') tol
        do i=1,MA
           write (fd, '(4096(ES14.7,1X))') AQ(i,:)
        enddo
        print '("file /tmp/A2qtt.dat written")'
        print '("norm of the difference with the orig.matrix = ",ES14.7)',err
     endif

  end subroutine test_low_rank


  !> simple unit diagonal
  function Aij_diag(d) result(A2)
  double precision, allocatable :: A2(:,:)
  integer, intent(IN) :: d
  integer :: i, N
     N= 2**d
     allocate(A2(N, N))
     A2= 0d0
     do i=1,N
        A2(i,i)= 1d0
     enddo
  end function Aij_diag


  !> arbitrary straight line
  function Aij_line(d) result(A2)
  double precision, allocatable :: A2(:,:)
  double precision, parameter :: a= -0.18, b= 0.63
  integer, intent(IN) :: d
  integer :: i, j, N
     N= 2**d
     allocate(A2(N, N))
     A2= 0d0
     do i=1,N
        j= int(N*(a*(i - 1)/dble(N) + b))
        if (j <= N .and. j > 0) A2(i,j)= 1d0
     enddo
  end function Aij_line


  !> pixelated filled "circle"
  function Aij_circle(d) result(A2)
  double precision, allocatable :: A2(:,:)
  integer, intent(IN) :: d
  integer :: i, j, N, M
     N= 2**(d-1)
     M = N*8/10 !< adjust the "radius" here
     allocate(A2(2*N, 2*N))
     A2= 0d0
     do i=1,2*N; do j=1,2*N
        if ((i - N)**2 + (j - N)**2 < M*M) A2(i,j)= 1d0
     enddo; enddo
  end function Aij_circle


  !> Energy distribution in a Green functin for a 2D wave eqn
  function Aij_radiant(d) result(A2)
  double precision, allocatable :: A2(:,:)
  integer, intent(IN) :: d
  !
  double precision :: x, y, zeta
  double precision, parameter :: a = 1.1d0
  integer :: i, j, N

     N= 2**d

     allocate(A2(N, N))
     A2= 0d0

     ! radiative flux energy
     do i=1,N; do j=1,N
        x= a*(2*dble(i)/dble(N) - 1d0)
        y= a*(2*dble(j)/dble(N) - 1d0)
        if (x*x + y*y < a*a) then
           zeta= dsqrt(x*x + y*y)
           A2(i,j)= datan2(zeta, 1d0 - zeta*zeta)
        endif
     enddo; enddo

  end function Aij_radiant


  !> this is an experiment to try to find low-rank representaiton
  !!      of a matrix A_ij = exp(2*pi*i*j/N)
  function Aij_clocks(d) result(A2)
  double precision, allocatable :: A2(:,:)
  integer, intent(IN) :: d
  !
  double precision :: x, y, zeta
  double precision, parameter :: a = 1.1d0
  double precision, parameter :: pi = 4d0*datan(1d0)
  integer :: i, j, N

     N= 2**d

     allocate(A2(N, N))
     A2= 0d0

     ! radiative flux energy
     do i=1,N; do j=1,N
        !A2(i,j)= dsin(2*pi*dble((i-1)*(j-1))/dble(N))
        !A2(i,j)= dexp(-2*pi*dble((i-1)*(j-1))/dble(N))
        A2(i,j)= dexp(-.999d0*dble((i-1)*(j-1))/dble(N))
     enddo; enddo

  end function Aij_clocks


  function Aij_patrick(d) result(A2)
  double precision, allocatable :: A2(:,:)
  integer, intent(IN) :: d
  integer :: i, j, N, M
  integer, parameter :: p = 5
  double precision, parameter :: r0=.8d0, amp=.2d0
  double precision, parameter :: pi= 4d0*datan(1d0)
  double precision :: r, phi
     N= 2**(d-1)
     M = N*8/10 !< adjust the "radius" here
     allocate(A2(2*N, 2*N))
     A2= 0d0
     do i=1,2*N; do j=1,2*N
        r = sqrt(dble((i - N)**2 + (j - N)**2))
        phi = atan2(dble(j-N),dble(i-n))
        if (N*(r0 + amp*cos(p*phi + pi*0.13)) > r) A2(i,j)= 1d0
     enddo; enddo
  end function Aij_patrick


  function Aij_permute(d) result(A2)
  use qtt_lib, only: binary_repr, integer_from_binary
  double precision, allocatable :: A2(:,:)
  integer, intent(IN) :: d
  integer :: i, j, ij, N, d2
  integer, allocatable :: ii(:), jj(:), ind(:)
     d2= d / 2
     if (d2*2 /= d) error stop "In Aij: d must be even!"

     N= 2**d2
     allocate(A2(N*N, N*N), ind(d))
     A2= 0d0
     do i=1,N
        ind= 0
        ii= binary_repr(i-1)
        ind(2:2*size(ii):2)= ii
        do j=1,N
           jj= binary_repr(j-1)
           ind(1:2*size(jj)-1:2)= jj
           ij= integer_from_binary(ind) + 1
           A2(j + (i - 1)*N, ij)= 1d0
           !print '(14(1X,I4))', i,j,ij,-1,ii(:),-1,jj(:),-1,ind(:)
        enddo
        !print '(/)'
     enddo
     !stop
  end function Aij_permute

end program
