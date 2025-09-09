!> @file dtt_matmul.f90
!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin (korobkin@lanl.gov)
!! @date   February 2025
!! @brief  Benchmarking FORTRAN matmul
!!

program main
use rnd_lib, only: random
implicit none
character(len=128) arg
!
integer :: i, err
integer :: batch_size = 10

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-b') then
        i=i+1;call getarg(i, arg)
        read(arg,*,iostat=err) batch_size
        if (err.ne.0) error stop "ERROR: when parsing -b <???>"
     elseif (trim(arg).eq.'-h') then
        call help_msg
        stop
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call bench_matmul(batch_size)

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Benchmark FORTRAN matmul", &
   "Usage: ./dtt_matmul.exe [-h]", &
   "Options:", &
   " -b <val>   : batch size [10]", &
   " -h         : print this help message"
  end subroutine


  subroutine bench_matmul(batch_size)
  use time_lib
  use thor_lib, only: dtt_tensor
  use mat_lib, only: thor_matmul
  implicit none
  integer, intent(IN) :: batch_size
  !
  double precision, allocatable :: A(:,:), B(:,:), C(:,:)
  double precision :: dt_matmul, dt_matmul_tp, dt_dgemm, dt_dgemm_tp
  integer :: mode_size, i, test

     ! call openblas_set_num_threads(32)

     call system_timer_init()
     print '("# THOR: matmul operation")'
     print '("# 1:mode_size 2:matmul[s] 3:matmul(transpose(A),B)[s]")'
     print '("# 4:dgemm(N,N)[s] 5:dgemm(T,N)[s]")'
     !mode_size = 2
     scaling_loop: do mode_size=1000,10000,1000
        allocate(A(mode_size,mode_size), B(mode_size,mode_size/100))
        call random(A)
        call random(B)
        call system_timer_start
        do test = 1, batch_size
           C = matmul(A, B)
        enddo
        call system_timer_stop
        dt_matmul= system_dt/batch_size
        call system_timer_start
        do test = 1, batch_size
           C = matmul(transpose(A), B)
        enddo
        call system_timer_stop
        dt_matmul_tp= system_dt/batch_size

        call system_timer_start
        do test = 1, batch_size
           C = thor_matmul(A, B)
        enddo
        call system_timer_stop
        dt_dgemm= system_dt/batch_size

        call system_timer_start
        do test = 1, batch_size
           C = thor_matmul(A, B, tsp_='tn')
        enddo
        call system_timer_stop
        dt_dgemm_tp= system_dt/batch_size

        print '(I10, 20(1X,ES10.3))', mode_size, &
                dt_matmul, dt_matmul_tp, dt_dgemm, dt_dgemm_tp
        deallocate(A, B, C)
     enddo scaling_loop

  end subroutine bench_matmul

end program
