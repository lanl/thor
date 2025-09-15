!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_rnd_lib.f90
!! @author Oleg Korobkin
!! @date   July 2025
!! @brief Testing rnd_lib
!!
program main
use rnd_lib, only: random, init_random_seed
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 2     ! array rank
integer :: n = 100   ! array size
integer :: rseed = 0 ! random seed (0 = from clock)

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
  call test_rnd(d, n)

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Testing rnd_lib", &
   "Usage: ./test_rnd_lib.exe [OPTIONS]", &
   "Options:", &
   " -d <val>   : array dimensionality [2]", &
   " -n <val>   : array size [100]", &
   " --rseed <num> : supply random seed [none]", &
   " -h         : print this help message"
  end subroutine


  subroutine test_rnd(d, n)
  use rnd_lib, only: randn
  use time_lib
  implicit none
  integer, intent(IN) :: d, n
  !
  integer :: nd
  double precision, allocatable :: A1(:), A2(:,:), A3(:,:,:), &
         A4(:,:,:,:), A5(:,:,:,:,:), A6(:,:,:,:,:,:)
  double precision :: A_mean, A_stddev, A_max, A_min

     ! initialize system timer
     call system_timer_init()

     print '("--------------------------------------")'
     print '("Testing normal distribution generator:")'
     print '(" - "I1"-dimensional array, "I10" along each side")', d, n
     nd = n**d
     print '(" - total elements: "I12)', nd
     print '("Generating . . .")'
     select case(d)
     case(1)
        call system_timer_start
        A1 = randn(n)
        call system_timer_stop
        A_mean = sum(A1)/dble(n)
        A_stddev = dsqrt(sum(A1**2)/dble(n) - A_mean**2)
        A_min = minval(A1)
        A_max = maxval(A1)

     case(2)
        call system_timer_start
        A2 = randn(n,n)
        call system_timer_stop
        A_mean = sum(A2)/dble(nd)
        A_stddev = dsqrt(sum(A2**2)/dble(nd) - A_mean**2)
        A_min = minval(A2)
        A_max = maxval(A2)

     case(3)
        call system_timer_start
        A3 = randn(n,n,n)
        call system_timer_stop
        A_mean = sum(A3)/dble(nd)
        A_stddev = dsqrt(sum(A3**2)/dble(nd) - A_mean**2)
        A_min = minval(A3)
        A_max = maxval(A3)

     case(4)
        call system_timer_start
        A4 = randn(n,n,n,n)
        call system_timer_stop
        A_mean = sum(A4)/dble(nd)
        A_stddev = dsqrt(sum(A4**2)/dble(nd) - A_mean**2)
        A_min = minval(A4)
        A_max = maxval(A4)

     case(5)
        call system_timer_start
        A5 = randn(n,n,n,n,n)
        call system_timer_stop
        A_mean = sum(A5)/dble(nd)
        A_stddev = dsqrt(sum(A5**2)/dble(nd) - A_mean**2)
        A_min = minval(A5)
        A_max = maxval(A5)

     case(6)
        call system_timer_start
        A6 = randn(n,n,n,n,n,n)
        call system_timer_stop
        A_mean = sum(A6)/dble(nd)
        A_stddev = dsqrt(sum(A6**2)/dble(nd) - A_mean**2)
        A_min = minval(A6)
        A_max = maxval(A6)

     end select

     print '("Elapsed time: ",ES10.3)', system_dt
     print '("Mean = ", ES14.7)', A_mean
     print '("Stddev = ", ES14.7)', A_stddev
     print '("Min = ", ES14.7)', A_min
     print '("Max = ", ES14.7)', A_max

  end subroutine test_rnd


end program
