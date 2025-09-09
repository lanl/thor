!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_thor_tensor_rand.f90
!! @author [IDB], [OGK]
!! @date  March 2024
!! @brief Testing random tensor generator
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
integer :: r = 2  ! tensor rank
integer :: l = 1  ! starting index
integer :: m      ! ending index
integer :: rseed = 0 ! random seed (0 = from clock)
integer, allocatable :: nn(:), rr(:)
logical :: lrandomize_nr = .false.
logical :: m_given = .false.
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
     elseif (trim(arg).eq.'--randomize-nr') then
        lrandomize_nr= .true.
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
  allocate(nn(d), rr(d+1), rnd(d))
  rr(1)= 1; rr(d+1)= 1
  if (lrandomize_nr) then
     call random(rnd)
     nn= 1 + int(rnd*n)
     call random(rnd)
     rr(2:d)= 1 + int(rnd(2:d)*r)
     do i=1,d-1
        rr(i)= min(rr(i), max(nn(i),nn(i+1)))
     enddo
  else
     nn= n; rr(2:d)= r
  endif
  if (.not.m_given) m = d
  call test_tensor_rand(nn, rr, l, m)
  deallocate(nn, rr, rnd)

contains

  subroutine help_msg()
  print '(11(/,A))', &
   "Test a random tensor generator", &
   "Usage: ./test.exe [-h] [-d <ord>] [-n <dim>] [-r <rank>] [--randomize-nr]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : tensor order <ord> [3]", &
   " -n <dim>  : tensor dimensions <dim> [4]", &
   " -r <rank> : rank of the cores [2]", &
   " -l <num>  : starting index in the array n [1]", &
   " -m <num>  : ending index in the array n [d]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]"
  end subroutine

  subroutine test_tensor_rand(n, r, l, m)
  use time_lib
  use thor_lib, only: dtt_tensor, dtt_tensor_rand, &
      dtt_tensor_rand_like, dealloc
  implicit none
  integer, intent(in):: n(:), r(:), l, m
  !
  type(dtt_tensor) :: v1, v2

     ! initialize system timer
     call system_timer_init()

     print '(/"Test 04: Testing a random tensor generator"/)'

     print '("-- [4-1] create random tensor v1")'
     call system_timer_start()
     v1 = dtt_tensor_rand(n, r_=r, l_=l, m_=m)
     call system_timer_stop()
     print '("v1 info:")'
     call v1% say
     print '("[+][TEST04-1][",F7.2,&
            &"s]: random vector generated OK: PASS"/)', system_dt

     print '("-- [4-2] Create random tensor v2")'
     call system_timer_start()
     v2 = dtt_tensor_rand_like(v1)
     call system_timer_stop()
     print '("v2 info:")'
     call v2% say
     print '("[+][TEST04-2][",F7.2,&
            &"s]: random vector v2 like v1 generated OK: PASS"/)', system_dt

     print '("-- [4-3] Deallocate random tensors v1 and v2")'
     call system_timer_start()
     call dealloc(v1)
     call dealloc(v2)
     call system_timer_stop()
     print '("[+][TEST04-3][",F7.2,&
            &"s]: random vectors deallocated OK: PASS"/)', system_dt

     print '("[+][TEST04]: PASS")'

 end subroutine test_tensor_rand


end program
