!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_thor_tensor_ones.f90
!! @author [IDB], [OGK]
!! @date  March 2024
!! @brief Testing "ones" tensor generator
!!

program main
use rnd_lib, only: random, init_random_seed
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 3  ! tensor order
integer :: n = 4  ! tensor dimension size
integer :: r = 2  ! tensor dimension size
integer :: l = 1  ! starting index
integer :: m = 3  ! ending index
integer :: rseed = 0 ! random seed (0 = from clock)
integer, allocatable :: nn(:)
logical :: lrandomize_n = .false.
logical :: m_given = .false.
double precision :: tol = 1d-12
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
     elseif (trim(arg).eq.'-l') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) l
        if (err.ne.0) error stop "ERROR: when parsing -l <???>"
        i= i + 1
     elseif (trim(arg).eq.'-m') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) m
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
        i= i + 1
     elseif (trim(arg).eq.'-t') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) tol
        if (err.ne.0) error stop "ERROR: when parsing -t <???>"
        m_given = .true.
        i= i + 1
     elseif (trim(arg).eq.'--randomize-n') then
        lrandomize_n= .true.
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call init_random_seed(rseed)
  allocate(nn(d), rnd(d))
  if (lrandomize_n) then
     call random(rnd)
     nn= 1 + int(rnd*n)
  else
     nn= n
  endif

  if (.not.m_given) m = d
  call test_tensor_ones(nn, tol, l, m)
  deallocate(nn, rnd)

contains

  subroutine help_msg()
  print '(11(/,A))', &
   "Test dtt_tensor_ones", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : tensor order <ord> [3]", &
   " -n <dim>  : tensor dimensions <dim> [4]", &
   " -l <num>  : starting index [1]", &
   " -m <num>  : ending index [3]", &
   " -t <tol>  : use this tolerance in SVD step [1e-12]", &
   " --randomize-n : generate arrays of random dimensions <= n [F]", &
   " --rseed <num> : supply random seed [none]"
  end subroutine

  subroutine test_tensor_ones(n, tol, l, m)
  use time_lib
  use thor_lib, only: dtt_tensor,dtt_tensor_ones,dtt_tensor_ones_like, &
                      operator(-)
  implicit none
  integer, intent(in) :: n(:), l, m
  double precision, intent(in) :: tol
  !
  type(dtt_tensor) :: v1, v2, vd
  integer :: info, sa, failed
  double precision, external :: dnrm2
  double precision, allocatable :: a1(:)
  double precision :: err, dt1, dt2, lg_ftsize
  double precision, parameter :: lg_ftsize_max = 9.0d0

     ! initialize system timer
     call system_timer_init()

     print '(/,"Test 01: Illustration of using tensor_ones",/)'
     failed= 0

     print '("-- [1-1] testing: v1 = dtt_tensor_ones(n, l_=l, m_=m)")'
     call system_timer_start
     v1 = dtt_tensor_ones(n, l_=l, m_=m)
     call system_timer_stop
     dt1 = system_dt
     lg_ftsize = v1% log10_fullsize()
     call v1% say
     print '("|v1|_2 = ",ES12.5)', v1% normb()

     if (lg_ftsize.lt.lg_ftsize_max) then
        a1 = v1% full()
        sa = size(a1)
        call system_timer_start
        err= dnrm2(sa, a1 - 1d0, 1)
        call system_timer_stop
        print '("|v1.full() - 1|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     else
        print '("v1.full() is too big: using normb for checks")'
        vd = v1 - 1d0
        call system_timer_start
        err= vd% normb()
        call system_timer_stop
        print '("|v1 - 1|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
        err= 0d0
     endif
     write (*,'("[+][TEST01-1]["F7.2" s]:&
        & tensor of ones v1 created: ")',&
        advance='no') dt1
     if (dabs(err).lt.tol) then
        write(*,'("PASS",/)')
     else
        write(*,'("FAIL",/)')
        failed= failed + 1
     endif

     print '("-- [1-2] testing: v2 = dtt_tensor_ones_like(v1)")'
     call system_timer_start
     v2 = dtt_tensor_ones_like(v1)
     call system_timer_stop
     dt1= system_dt
     if (lg_ftsize.lt.lg_ftsize_max) then
        a1 = v2% full()
        sa = size(a1)
        call system_timer_start
        err= dnrm2(sa, a1 - 1d0, 1)
        call system_timer_stop
        print '("|v2.full() - 1|_2 = ",ES12.5," ["F7.2" s]")', err,system_dt
     else
        vd = v2 - 1d0
        call system_timer_start
        err= vd% normb()
        call system_timer_stop
        print '("|v2 - 1|_2 = ",ES12.5," ["F7.2" s]")', err,system_dt
     endif
     call v2% say
     write (*,'("[+][TEST01-2]["F7.2" s]:&
        & tensor of ones v2 like v1 created: ")',&
        advance='no') dt1
     if (dabs(err).lt.tol) then
        write(*,'("PASS",/)')
     else
        write(*,'("FAIL",/)')
        failed= failed + 1
     endif

     print '("-- [1-3] testing: v2% svd(rmax=1, tol="ES12.5")")', tol
     call system_timer_start
     call v2% svd(rmax_=1, tol_=tol)
     call system_timer_stop
     dt1= system_dt

     if (lg_ftsize.lt.lg_ftsize_max) then
        a1 = v2% full()
        sa = size(a1)
        call system_timer_start
        err= dnrm2(sa, a1 - 1d0, 1)
        call system_timer_stop
        print '("|v2.full() - 1|_2 = ",ES12.5," ["F7.2" s]")', err,system_dt
     else
        vd = v2 - 1d0
        call system_timer_start
        err= vd% normb()
        call system_timer_stop
        print '("|v2 - 1|_2 = ",ES12.5," ["F7.2" s]")', err,system_dt
     endif
     write(*,'("[+][TEST01-3][",F7.2," s]: &
        &tensor of ones SVD compression (tol=",ES14.7,"):")',&
        advance='no') dt1, tol
     if (dabs(err).lt.tol) then
        write(*,'("PASS",/)')
     else
        write(*,'("FAIL",/)')
        failed= failed + 1
     endif

     if(failed.eq.0) then
        print '("[+][TEST01]: PASSED")'
     else
        print '("[-][TEST01]: FAILED with ",I2," error(s)")', failed
        error stop failed
     end if

  end subroutine test_tensor_ones

end program
