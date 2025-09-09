!> @file dtt_ones.f90
!!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin (korobkin@lanl.gov)
!! @date   February 2025
!! @brief  Benchmarking "dtt_ones" tensor generator
!!

program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg
!
integer :: err

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call bench__dtt_tensor_ones()

contains

  subroutine help_msg()
  print '(10(/,A))', &
   "Benchmark dtt_tensor_ones", &
   "Usage: ./dtt_tensor_ones.exe [-h]", &
   "Options:", &
   " -h : print this help message"
  end subroutine

  subroutine bench__dtt_tensor_ones()
  use time_lib
  use thor_lib, only: dtt_tensor, dtt_tensor_ones
  implicit none
  integer, allocatable :: n(:)
  integer, parameter :: batch_size = 1000
  integer, parameter :: mode_size  = 100  !< mode size
  integer :: d                            !< tensor dimension
  integer :: i, test
  type(dtt_tensor) :: v

     call system_timer_init()

     print '("# THOR test: dtt_tensor_ones, mode size = 100")'
     print '("# 1:dimension 2:time[s]")'
     d = 2
     scaling_loop: do i = 1, 11
        allocate(n(d))
        n = mode_size
        call system_timer_start
        test_loop: do test = 1, batch_size
           v = dtt_tensor_ones(n)
        enddo test_loop
        call system_timer_stop
        print '(I10, 1X, ES12.5)', d, system_dt/dble(batch_size)
        d = d*2
        deallocate(n)
     enddo scaling_loop
     

  end subroutine bench__dtt_tensor_ones

end program
