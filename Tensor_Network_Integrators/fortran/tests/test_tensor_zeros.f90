!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_thor_tensor_zeros.f90
!! @author [IDB], [OGK]
!! @date  March 2024
!! @brief Testing "zeros" tensor generators
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
integer :: m = 3  ! ending index
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
  call test_zero_tensor_generator(nn, rr, l, m)
  deallocate(nn, rr, rnd)

contains

  subroutine help_msg()
  print '(12(/,A))', &
   "Test zero tensor constructors", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : tensor order <ord> [3]", &
   " -n <dim>  : tensor dimensions <dim> [4]", &
   " -r <rank> : rank of the cores [2]", &
   " -l <num>  : starting index [1]", &
   " -m <num>  : ending index [3]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]"
  end subroutine

  subroutine test_zero_tensor_generator(n, r, l, m)
  use time_lib
  use thor_lib, only: dtt_tensor, dtt_tensor_zeros, &
                      dtt_tensor_zeros_like, dealloc
  implicit none
  integer, intent(in):: n(:), r(:), l, m
  !
  type(dtt_tensor) :: v0, v1
  double precision :: nrm, tol
  integer :: failed

     ! initialize system timer
     call system_timer_init()
     failed = 0
     tol = 1d-12

     print '(/,"Test 00: Routines to create zero tt-tensor",/)'

     print '("-- [0-1] testing dtt_tensor_zeros(n,r_=r,l_=l,m_=m)")'
     call system_timer_start
     v0 = dtt_tensor_zeros(n,r_=r,l_=l,m_=m)
     call system_timer_stop
     call v0% say
     nrm = v0% normb()
     print '("|A|_2 = "F3.1)', nrm
     write(*, '("[+][TEST00-1]["F7.2" s]: &
        &zero tensor v0 using dtt_tensor_zeros: ")', &
        advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS",/)')
     else
        write(*,'("FAIL",/)')
        failed= failed + 1
     endif

     print '("-- [0-2] testing dtt_tensor_zeros_like(v0)")'
     call system_timer_start
     v1 = dtt_tensor_zeros_like(v0)
     call system_timer_stop
     call v1% say
     nrm = v1% normb()
     print '("|A|_2 = "F3.1)', nrm
     write(*, '("[+][TEST00-2]["F7.2" s]: &
        &zero tensor v1 using dtt_tensor_zeros_like(v0): ")', &
        advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS",/)')
     else
        write(*,'("FAIL",/)')
        failed= failed + 1
     endif

     print '("-- [0-3] testing dtt_tensor(n,l_=l,m_=m); v0 = 0")'
     call system_timer_start
     v0 = dtt_tensor(n,l_=l,m_=m)
     call system_timer_stop
     call v0% say
     nrm = v0% normb()
     print '("|A|_2 = "F3.1)', nrm
     write(*, '("[+][TEST00-3]["F7.2" s]: &
        &zero tensor v0 using tt constructor dtt_tensor(): ")', &
        advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS",/)')
     else
        write(*,'("FAIL",/)')
        failed= failed + 1
     endif

     print '("-- [0-4] testing deallocation:")'
     ![0-4] Deallocated all tensors
     call system_timer_start
     call dealloc(v0)
     call dealloc(v1)
     call system_timer_stop
     print '("[+][TEST00-4][",F7.2," s]: &
           &tensor v0 deallocated: PASS")', system_dt

     if(failed.eq.0) then
        print '("[+][TEST00]: PASSED")'
     else
        print '("[-][TEST00]: FAILED with ",I2," error(s)")', failed
        error stop failed
     end if

  end subroutine test_zero_tensor_generator

end program
