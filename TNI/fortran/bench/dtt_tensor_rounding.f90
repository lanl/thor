!> @file dtt_tensor_rounding.f90
!!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin (korobkin@lanl.gov)
!! @date   February 2025
!! @brief  Benchmarking dtt_tensor rounding operation
!!

program main
use rnd_lib, only: random
implicit none
character(len=128) arg
!
integer :: i, err
double precision  :: epsi = 1d-8
logical           :: epsi_given = .false.
integer           :: batch_size = 100
logical           :: batch_size_given = .false.
integer           :: low_rank = 5
logical           :: low_rank_given = .false.
character(len=20) :: which_test

  which_test = "default"

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-e') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) epsi
        if (err.ne.0) error stop "ERROR: when parsing -e <???>"
        epsi_given = .true.
     elseif (trim(arg).eq.'-b') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) epsi
        read(arg,*,iostat=err) batch_size
        if (err.ne.0) error stop "ERROR: when parsing -b <???>"
        batch_size_given = .true.
     elseif (trim(arg).eq.'-r') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) epsi
        read(arg,*,iostat=err) low_rank
        if (err.ne.0) error stop "ERROR: when parsing -r <???>"
        low_rank_given = .true.
     elseif (trim(arg).eq.'-h') then
        call help_msg
        stop
     endif
  enddo

  select case (trim(which_test))
  case ("default")
     if (.not.epsi_given) epsi = 1e-8
     if (.not.batch_size_given) batch_size = 10
     if (.not.low_rank_given) low_rank = 5
     call bench__dtt_tensor_4d_round(batch_size, epsi, low_rank)
  endselect

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Benchmark dtt_tensor rounding operation", &
   "Usage: ./dtt_tensor_rounding.exe [-h] [options]", &
   "Options:", &
   " -e <val>   : rounding tolerance epsilon [1e-8]", &
   " -b <val>   : batch size [10]", &
   " -r <val>   : low-rank value [5]", &
   " -h         : print this help message"
  end subroutine

  subroutine bench__dtt_tensor_4d_round(batch_size, epsi, low_rank)
  use time_lib
  use thor_lib, only: dtt_tensor, dtt_random_ortho
  implicit none
  integer, intent(IN) :: batch_size, low_rank
  double precision, intent(IN) :: epsi
  !
  double precision, allocatable :: A4(:,:,:,:)
  integer :: mode_size, n(4), r(0:4) !< mode size
  integer :: i, test
  type(dtt_tensor) :: vr

     call system_timer_init()
     print '("# THOR test: round a random low-rank 4D tensor")'
     print '("# Modes = (M, 2M, 3M, 4M), ranks = (1, "I4","I4",..., 1)")', &
           low_rank, low_rank
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:time[s] 3-7:ranks")'

     r =  [1, low_rank, low_rank, low_rank, 1]
     mode_size = 125
     scaling_loop: do i = 1,15
        n = [mode_size, mode_size*2, mode_size*3, mode_size*4]
        vr = dtt_random_ortho(n, r_=r)
        call system_timer_start
        test_loop: do test = 1, batch_size
           call vr% mround(epsi)
        enddo test_loop
        call system_timer_stop
        print '(I10, 1X, ES12.5, 14(I6,1X))', &
              mode_size, system_dt/dble(batch_size), vr%r(0:4)
        mode_size = mode_size*2
     enddo scaling_loop

  end subroutine bench__dtt_tensor_4d_round

end program
