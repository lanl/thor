!> @file matmuls.f90
!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin (korobkin@lanl.gov)
!! @date   February 2025
!! @brief  Benchmarking matrix multiplication and transpose
!!

program main
use rnd_lib, only: random, init_random_seed
implicit none
character(len=128) arg
!
integer,parameter :: Ntim = 10
double precision  :: timers(Ntim)
integer :: i, err
integer :: rseed = 0 ! random seed (0 = from clock)
double precision  :: epsi = 1d-8
integer           :: batch_size = 10
integer           :: n = 20
integer           :: num_threads = 1

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-n') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) n
        if (err.ne.0) error stop "ERROR: when parsing -n <???>"
        i= i + 1
     elseif (trim(arg).eq.'-b') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) batch_size
        if (err.ne.0) error stop "ERROR: when parsing -b <???>"
        i= i + 1
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        call init_random_seed(rseed)
        i= i + 1
     elseif (trim(arg).eq.'-h') then
        call help_msg
        stop
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call bench__matmuls(n, batch_size, num_threads)

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Benchmark matrix multiplication and transpose", &
   "Usage: ./matmuls.exe [-h]", &
   "Options:", &
   " -n <val>   : matrix dimension [100]", &
   " -b <val>   : batch size [10]", &
   " --rseed <num> : supply random seed [none]", &
   " -h         : print this help message"
  end subroutine


  subroutine bench__matmuls(n, batch_size, num_threads)
  use time_lib
  use mat_lib, only: thor_matmul
  use rnd_lib, only: randn
  implicit none
  integer, intent(IN) :: n          !< matrix size
  integer, intent(IN) :: batch_size
  integer, intent(IN) :: num_threads
  !
  integer :: mode_size, i, test
  double precision, allocatable :: A(:,:), B(:,:), C(:,:)
  double precision :: t_matmul, t_thmatmul, t_transpose

     call openblas_set_num_threads(num_threads)

     call system_timer_init()
     print '("# THOR: matrix multiplication timing")'
     print '("# 1:matrix size 2:matmul[s] 3:thor_matmul[s] 4:transpose[s]")'
     mode_size = n
     scaling_loop: do i=1,8
        A = randn(mode_size,2*mode_size)
        B = randn(2*mode_size,3*mode_size)
        call system_timer_start
        do test = 1, batch_size
           C = matmul(A, B)
        enddo
        call system_timer_stop
        t_matmul = system_dt/dble(batch_size)

        call system_timer_start
        do test = 1, batch_size
           C = thor_matmul(A, B)
        enddo
        call system_timer_stop
        t_thmatmul = system_dt/dble(batch_size)

        call system_timer_start
        do test = 1, batch_size
           C = transpose(C)
        enddo
        call system_timer_stop
        t_transpose = system_dt/dble(batch_size)

        print '(I10, 20(1X,ES10.3))', mode_size, t_matmul, t_thmatmul, t_transpose
        mode_size = 2*mode_size
        deallocate(A, B, C)

     enddo scaling_loop


  end subroutine bench__matmuls

end program
