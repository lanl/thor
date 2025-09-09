!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_amen_ops.f90
!! @author [IDB], [OGK]
!! @date  October 2023
!! @brief Testing amen_mv and amen_solve
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
integer :: l = 1     ! starting index
integer :: m         ! ending index
integer :: rseed = 0 ! random seed (0 = from clock)
integer, allocatable :: q(:), s(:), rr(:)
logical :: lrandomize_nr = .false.
logical :: m_given = .false.
logical :: lsave_tt = .true.
logical :: lread_tt = .false.
double precision, allocatable :: rnd(:)

  m = d
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
     elseif (trim(arg).eq.'-l') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) l
        if (err.ne.0) error stop "ERROR: when parsing -l <???>"
        i= i + 1
     elseif (trim(arg).eq.'-m') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) m
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
        m_given = .true.
        i= i + 1
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     elseif (trim(arg).eq.'--read-tt') then
        lread_tt = .true.
     elseif (trim(arg).eq.'--dontsave-tt') then
        lsave_tt = .false.
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
     q = n; s = n; rr = r
     rr(1) = 1; rr(d+1) = 1
  endif

  if (.not.m_given) m = d
  call test_amen_mv(q, s, rr, l, m, lsave_tt, lread_tt)
  deallocate(q, s, rr, rnd)

contains

  subroutine help_msg()
  print '(14(/,A))', &
   "Test amen-mv and amen-solve functions", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : matrix order <ord> [4]", &
   " -n <dim>  : max dimension for each mode [3]", &
   " -r <rank> : rank of the cores [2]", &
   " -l <num>  : starting index [1]", &
   " -m <num>  : ending index [4]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]", &
   " --dontsave-tt : don't save files [saves to /tmp/* by default]", &
   " --read-tt : read files tx.sdv, mA.sdv and ty.sdv (from tests/testdata) [false]"
  end subroutine


  subroutine test_amen_mv(qq, ss, r, l, m, lsave_tt, lread_tt)
  use time_lib
  use matlab_struct_module, only : pprint_matrix, array3d
  use thorio_lib, only: dtt_write_sdv_file, dtt_read_sdv_file, &
                     dtt_write_ascii_file
  use thor_lib
  use ttamen_lib, only: amen_solve, amen_mv
  implicit none
  integer, intent(in) :: qq(:), ss(:), r(:), l, m
  logical, intent(in) :: lsave_tt, lread_tt
  !
  double precision, parameter :: err_tol = 1d-15
  double precision :: err
  integer :: i
  integer, allocatable :: nn(:), rr(:)
  character(:), allocatable :: strs(:)
  type(dtt_tensor) :: ta, tx, ty, tx1, dtx
  type(dtt_matrix) :: mA
  type(array3d), allocatable :: cy(:)
  double precision, allocatable :: y(:,:,:,:)
  character(len=32) :: fname

     ! initialize system timer
     call system_timer_init()

     print '("-----")'
     print '("testing amen_mv")'
     !! print '("matlab code :")'
     !! print '(">>  ta = tt_ones(9,4)")'
     !! print '(">>  mA = tt_matrix(a, 3, 3)")'
     !! print '(">>  tx = 5*tt_ones(3,4)")'
     !! print '(">>  [y, z] = amen_mv( A, x, 1e-8)")'

     !! ta = dtt_tensor_ones([9, 9, 9, 9])
     !! mA = dtt_matrix(ta, [3,3,3,3], [3,3,3,3])
     !! !print '("ta = ")'
     !! !call pprint(ta)
     !! !print '("mA = ")'
     !! !call pprint(mA)
     !! tx = 5d0*dtt_tensor_ones([3, 3, 3, 3])
     !! ty = amen_mv(mA, tx, 1d-8)
     !! print '("ty = ")'
     !! call pprint(ty)
     !! print '("expected matlab result:")'
     !! print '(">> y.core")'
     !! print '("")'
     !! print '("ans =")'
     !! print '("")'
     !! print '("    4.4860")'
     !! print '("    4.4860")'
     !! print '("    4.4860")'
     !! print '("   -4.4860")'
     !! print '("   -4.4860")'
     !! print '("   -4.4860")'
     !! print '("   -4.4860")'
     !! print '("   -4.4860")'
     !! print '("   -4.4860")'
     !! print '("    4.4860")'
     !! print '("    4.4860")'
     !! print '("    4.4860")'
     !! y = reshape(ty% full(), [3,3,3,3])
     !! print '("y(1,3,2,1) = ", F4.0)', y(1,3,2,1)
     !! call pprint_matrix(reshape(y, [9, 9]), frmt="(F4.0)")

     nn = qq*ss
     print '("------------------------------")'
     if (lread_tt) then
        print '("reading tensor tx from file:")'
        tx= dtt_read_sdv_file('tests/testdata/tx.sdv')
     else
        print '("generating arbitrary tensor x:")'
        tx = dtt_tensor_rand(ss, r_=r)
     endif
     print '("tensor tx:")'
     call tx% say
     if (lsave_tt) then
        call dtt_write_sdv_file(tx, '/tmp/tx.sdv')
        call dtt_write_ascii_file(tx, '/tmp/tx.dat', verb_=.true.)
     endif

     print '("------------------------------")'
     if (lread_tt) then
        print '("reading matrix mA from file:")'
        ta= dtt_read_sdv_file('tests/testdata/mA.sdv')
        mA= dtt_matrix(ta, ta%n(l:m)/tx%n(l:m), tx%n(l:m))
     else
        print '("generating arbitrary diagonally dominant matrix A:")'
        mA= dtt_matrix_rand(qq, ss, r_=r)
!        mA= mA + dttm_diag(100d0*dtt_tensor_ones(qq) + dtt_tensor_rand(qq, r_=r))
     endif
     print '("matrix mA:")'
     call mA% say
     if (lsave_tt) then ! TODO: have dttm_write as well
        call dtt_write_sdv_file(mA, '/tmp/mA.sdv')
!!        call dttm_write_ascii_file(mA, '/tmp/mA.dat', verb_=.true.)
     endif

     print '("------------------------------")'
     print '("testing amen_mv(mA, tx):")'
     call system_timer_start
     call amen_mv(ty, mA, tx, 1d-12)
     call system_timer_stop
     print '("tensor ty = mA*tx:")'
     call ty% say
     print '("ty = amen_mv(mA, tx) timing: "F9.4" s")', system_dt
     if (lsave_tt) then ! TODO: have dttm_write as well
        call dtt_write_sdv_file(ty, '/tmp/ty.sdv')
        call dtt_write_ascii_file(ty, '/tmp/ty.dat', verb_=.true.)
     endif
     if (lread_tt) then
        print '("reading tensor ty1 from file:")'
        tx1= dtt_read_sdv_file('tests/testdata/ty.sdv')
        dtx= tx1 - ty
        print '("|ty - ty1|/|ty| = "ES14.7)', dtx% normb()/ty% normb()
     endif

     print '("------------------------------")'
     print '("testing tx = amen_solve(mA, ty, tol):")'
     if (lread_tt) then
        print '("reading initial guess tx0 from file:")'
        tx1= dtt_read_sdv_file('tests/testdata/tx0.sdv')
     else
        if (allocated(rr)) deallocate(rr)
        allocate(rr(0:d))
        rr = 2; rr(0) = 1; rr(d) = 1
        tx1= dtt_random_ortho(qq, r_=rr, lr_=.false.) ! TODO: make more intutivie "left-right"
     endif
     print '("initial guess for tx (tensor tx1):")'
     call tx1% say
     if (lsave_tt) call dtt_write_ascii_file(tx1, '/tmp/tx0.dat', verb_=.true.)
     call system_timer_start
     tx1= amen_solve(mA, ty, 1d-10, verb_=1)
     call system_timer_stop
     if (lsave_tt) call dtt_write_ascii_file(tx1, '/tmp/tx1.dat', verb_=.true.)
     call tx1% say
     print '("tx1 = amen_solve(mA, ty) timing: "F9.4" s")', system_dt
     dtx= tx1 - tx
     print '("|tx - tx1|/|tx| = "ES14.7)', dtx% normb()/tx% normb()

  end subroutine test_amen_mv


end program
