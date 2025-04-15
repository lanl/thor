!> @file b_randn.f90
!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin
!! @date   July 2025
!! @brief  Benchmarking functions from rnd_lib
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
integer :: b = 10    ! batch size
integer :: rseed = 0 ! random seed (0 = from clock)

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     elseif (trim(arg).eq.'-b') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) b
        if (err.ne.0) error stop "ERROR: when parsing -b <???>"
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
  call bench__randn(b)

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Benchmark functions from rnd_lib", &
   "Usage: ./b_randn.exe [OPTIONS]", &
   "Options:", &
   " -b <num>      : batch size for timing [10]", &
   " --rseed <num> : supply random seed [none]", &
   " -h            : print this help message"
  end subroutine


  subroutine bench__randn(batch_size)
  use rnd_lib, only: rndmat, randn
  use time_lib
  implicit none
  integer, intent(IN) :: batch_size
  !
  double precision, allocatable :: v(:)
  integer :: i, n, test
  double precision :: dt_rndmat, dt_randn

     ! initialize system timer
     call system_timer_init()

     print '("# THOR test for rndmat and randn from rnd_lib")'
     print '("# Batch size: b = "I6)', batch_size
     print '("# 1:size 2:rndmat[s] 3:randn[s]")'
     n = 1000
     scaling_loop: do i = 3,8
        call system_timer_start
        test_loop: do test = 1, batch_size
           v = randn(n)
        enddo test_loop
        call system_timer_stop
        dt_randn = system_dt/dble(batch_size)

        call system_timer_start
        rndmat_loop: do test = 1, batch_size
           v = rndmat(n)
        enddo rndmat_loop
        call system_timer_stop
        dt_rndmat = system_dt/dble(batch_size)

        print '(I10, 2(1X,ES12.5))', n, dt_rndmat, dt_randn
        n = n*10
     enddo scaling_loop

  end subroutine bench__randn

end program

