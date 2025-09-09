!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_amen_mm.f90
!! @author [IDB], [OGK]
!! @date  October 2023
!! @brief Testing amen_mm
!!
program main
use rnd_lib, only: random, init_random_seed
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 4     ! tensor order
integer :: n = 3     ! tensor dimension size
integer :: r = 2     ! tensor rank
integer :: rseed = 0 ! random seed (0 = from clock)
integer, allocatable :: q(:), s(:), rr(:)
logical :: lrandomize_nr = .false.
logical :: lsave_tt = .true.
logical :: lread_tt = .false.
double precision, allocatable :: rnd(:)

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     elseif (trim(arg).eq.'-d') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) d
        if (err.ne.0) error stop "ERROR: when parsing -d <???>"
        i= i + 1
     elseif (trim(arg).eq.'-n') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) n
        if (err.ne.0) error stop "ERROR: when parsing -n <???>"
        i= i + 1
     elseif (trim(arg).eq.'-r') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) r
        if (err.ne.0) error stop "ERROR: when parsing -r <???>"
        i= i + 1
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     elseif (trim(arg).eq.'--save-tt') then
        lsave_tt = .true.
     elseif (trim(arg).eq.'--randomize-nr') then
        lrandomize_nr= .true.
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call init_random_seed(rseed)
  allocate(q(d), s(d), rr(d+1), rnd(d))
  rr(1)= 1; rr(d+1)= 1
  if (lrandomize_nr) then
     call random(rnd)
     q= 1 + int(rnd*n)
     call random(rnd)
     s= 1 + int(rnd*n)
     call random(rnd)
     rr(2:d)= 1 + int(rnd(2:d)*r)
     do i=1,d-1
        rr(i)= min(rr(i), max(q(i)*s(i), q(i+1)*s(i+1)))
     enddo
  else
     !qq = [3, 2, 3, 1]
     !ss = [3, 4, 2, 7]
     q = n; s = n + 1; rr = r
     rr(1) = 1; rr(d+1) = 1
  endif

  !! call permute ! test permutation function
  call test_amen_mm(q, s, rr, lsave_tt, lread_tt)
  deallocate(q, s, rr, rnd)

contains

  subroutine help_msg()
  print '(11(/,A))', &
   "Test amen-mv and amen-solve functions", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : matrix order <ord> [4]", &
   " -n <dim>  : max dimension for each mode [3]", &
   " -r <rank> : rank of the cores [2]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]", &
   " --save-tt : save files to /tmp/* [false]", &
   " --read-tt : read files mA.dat, mB.dat, mC.dat (from /tmp) [false]"
   !" -l <num>  : starting index [1]", &
   !" -m <num>  : ending index [4]", &
  end subroutine


  subroutine permute
  use matlab_struct_module, only: pprint_matrix
  use matrix_util, only: perm3d
  use rnd_lib, only: randn
  use mat_lib, only: thor_matmul
  implicit none
  double precision, allocatable, target :: A(:,:), B(:,:)
  double precision, pointer :: wp(:)
  integer :: n1, n2, n3, i, j, l
     n1 = 3
     n2 = 5
     n3 = 2
     allocate(A(n1,n2*n3))
     do i = 1,n1
        do l = 1,n3
           do j = 1,n2
              A(i,j+(l-1)*n2) = dble(i) + 0.1d0*dble(j) + 0.01d0*dble(l)
           enddo
        enddo
     enddo
     print '("Initial matrix:")'
     call pprint_matrix(A, line_len_=160)
     wp(1:n1*n2*n3)=> A
     call perm3d(wp,n1,n2,n3,132)
     print '("After 132 permutation:")'
     call pprint_matrix(A, line_len_=160)
  end subroutine permute


  subroutine test_amen_mm(qq, ss, r, lsave_tt, lread_tt)
  use time_lib
  use matlab_struct_module, only : pprint_matrix, array3d
  use thorio_lib, only: dttm_write_ascii_file
  use thor_lib
  use ttamen_lib, only: amen_mm
  implicit none
  integer, intent(in) :: qq(:), ss(:), r(:)
  logical, intent(in) :: lsave_tt, lread_tt
  !
  double precision, parameter :: err_tol = 1d-15
  double precision :: err
  integer :: i
  integer, allocatable :: nn(:), rr(:)
  character(:), allocatable :: strs(:)
  type(dtt_tensor) :: tx, ty
  type(dtt_matrix) :: mA, mB, mC
  type(array3d), allocatable :: cy(:)
  double precision, allocatable :: y(:,:,:,:)
  character(len=32) :: fname

     ! initialize system timer
     call system_timer_init()

     print '("------------------------------")'
     print '("generating arbitrary diagonally dominant matrix A:")'
     mA= dtt_matrix_rand(qq, ss, r_=r)
     if (lsave_tt) then
        call dttm_write_ascii_file(mA,"/tmp/mA.dat",verb_=.true.)
     else
        print '("matrix mA:")'
        call mA% say
     endif

     mB= dtt_matrix_rand(ss, qq+2, r_=r)
     if (lsave_tt) then
        call dttm_write_ascii_file(mB,"/tmp/mB.dat",verb_=.true.)
     else
        print '("matrix mB:")'
        call mB% say
     endif

     print '("------------------------------")'
     print '("testing amen_mm(mA, mB):")'
     call system_timer_start
     call amen_mm(mC, mA, mB, 1d-10)
     call system_timer_stop
     if (lsave_tt) then
        call dttm_write_ascii_file(mC,"/tmp/mC.dat",verb_=.true.)
        print '("You can now check the result with Matlab:")'
        print '("$ matrun.sh tests/testdata/test_amen_mm.m")'
     else
        print '("resulting matrix mC = mA @ mB:")'
        call mC% say
     endif

  end subroutine test_amen_mm

end program
