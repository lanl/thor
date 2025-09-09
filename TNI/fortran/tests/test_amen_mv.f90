!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_amen_mv.f90
!! @author [IDB], [OGK]
!! @date  October 2023
!! @brief Testing amen_mv
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
  else
     !qq = [3, 2, 3, 1]
     !ss = [3, 4, 2, 7]
     q = n; s = n + 1; rr = r
     rr(1) = 1; rr(d+1) = 1
  endif

  !! call permute ! test permutation function
  call test_amen_mv(q, s, rr, lsave_tt)
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
   " --save-tt : save files to /tmp/* [false]"
  end subroutine


  subroutine test_amen_mv(qq, ss, r, lsave_tt)
  use time_lib
  use matlab_struct_module, only : pprint_matrix, array3d
  use thorio_lib, only: dtt_write_ascii_file, dttm_write_ascii_file
  use thor_lib
  use ttamen_lib, only: amen_mv
  implicit none
  integer, intent(in) :: qq(:), ss(:), r(:)
  logical, intent(in) :: lsave_tt
  !
  double precision, parameter :: err_tol = 1d-15
  double precision :: err
  integer :: i
  integer, allocatable :: nn(:), rr(:)
  character(:), allocatable :: strs(:)
  type(dtt_tensor) :: tx, ty
  type(dtt_matrix) :: mA
  character(len=32) :: fname

     ! initialize system timer
     call system_timer_init()

     print '("------------------------------")'
     print '("generating arbitrary tt-matrix A:")'
     mA= dtt_matrix_rand(qq, ss, r_=r)
     if (lsave_tt) then
        call dttm_write_ascii_file(mA,"/tmp/mA.dat",verb_=.true.)
     else
        print '("matrix mA:")'
        call mA% say
     endif

     print '("------------------------------")'
     print '("generating arbitrary tt-tensor x:")'
     tx = dtt_tensor_rand(ss, r_=r)
     if (lsave_tt) then
        call dtt_write_ascii_file(tx, '/tmp/tx.dat', verb_=.true.)
     else
        print '("tensor tx:")'
        call tx% say
     endif

     print '("------------------------------")'
     print '("testing ty = amen_mv(mA, tx):")'
     call system_timer_start
     call amen_mv(ty, mA, tx, 1d-12)
     call system_timer_stop
     if (lsave_tt) then
        call dtt_write_ascii_file(ty, "/tmp/ty.dat", verb_=.true.)
        print '("You can now check the result with Matlab:")'
        print '("$ matrun.sh tests/testdata/test_amen_mv.m")'
     else
        print '("resulting tt-tensor ty = mA * tx:")'
        call ty% say
     endif

  end subroutine test_amen_mv

end program
