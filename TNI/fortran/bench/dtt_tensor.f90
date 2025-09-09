!> @file dtt_tensor.f90
!!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin (korobkin@lanl.gov)
!! @date   February 2025
!! @brief  Benchmarking "dtt_tensor" constructor by approximating full data
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
character(len=20) :: which_test

  which_test = "a2"

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
        read(arg,*,iostat=err) batch_size
        if (err.ne.0) error stop "ERROR: when parsing -b <???>"
        batch_size_given = .true.
     elseif (trim(arg).eq.'-a2') then
        which_test = "a2"
     elseif (trim(arg).eq.'-a3') then
        which_test = "a3"
     elseif (trim(arg).eq.'-a4') then
        which_test = "a4"
     elseif (trim(arg).eq.'-a4-lrank') then
        which_test = "a4-lrank"
     elseif (trim(arg).eq.'-a4-mround') then
        which_test = "a4-mround"
     elseif (trim(arg).eq.'-h') then
        call help_msg
        stop
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  select case (trim(which_test))
  case ("a2")
     if (.not.epsi_given) epsi = 1e-8
     if (.not.batch_size_given) batch_size = 10
     call bench__dtt_tensor_2d(batch_size, epsi)
  case ("a3")
     if (.not.epsi_given) epsi = 1e-8
     if (.not.batch_size_given) batch_size = 10
     call bench__dtt_tensor_3d(batch_size, epsi)
  case ("a4")
     if (.not.epsi_given) epsi = 1e-8
     if (.not.batch_size_given) batch_size = 4
     call bench__dtt_tensor_4d(batch_size, epsi)
  case ("a4-lrank")
     if (.not.epsi_given) epsi = 1e-8
     if (.not.batch_size_given) batch_size = 10
     call bench__dtt_tensor_lrank_4d(batch_size, epsi)
  case ("a4-mround")
     if (.not.epsi_given) epsi = 1e-8
     if (.not.batch_size_given) batch_size = 10
     call bench__dtt_tensor_mround_4d(batch_size, epsi)
  endselect

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Benchmark dtt_tensor constructor: from various n-dim arrays", &
   "Usage: ./dtt_tensor.exe [-h]", &
   "Options:", &
   " -e <val>   : rounding tolerance epsilon [1e-8]", &
   " -b <val>   : batch size [depens on the test]", &
   " -a2        : create tt from random 2D full-rank tensor A2(:,:)", &
   " -a3        : create tt from random 3D full-rank tensor A3(:,:,:)", &
   " -a4        : create tt from random 4D full-rank tensor A4(:,:,:,:)", &
   " -a4-lrank  : create tt from random 4D low-rank tensor", &
   " -a4-mround : round random low-rank 4D tensor using mround", &
   " -h         : print this help message"
  end subroutine

  subroutine bench__dtt_tensor_2d(batch_size, epsi)
  use time_lib
  use thor_lib, only: dtt_tensor
  implicit none
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  !
  double precision, allocatable :: A2(:,:)
  integer :: mode_size, n(2) !< mode size & dimensions
  integer :: i, test
  type(dtt_tensor) :: v

     call system_timer_init()
     print '("# THOR test: dtt_tensor(A2,epsi)")'
     print '("# Modes = (M, 2M)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:time[s] 3-5:ranks")'
     scaling_loop: do mode_size = 100, 1000, 100
        n = [mode_size, mode_size*2]
        allocate(A2(n(1), n(2)))
        call random(A2)
        call system_timer_start
        test_loop: do test = 1, batch_size
           v = dtt_tensor(A2, eps_=epsi)
        enddo test_loop
        call system_timer_stop
        print '(I10, 1X, ES12.5, 14(I6,1X))', &
              mode_size, system_dt/dble(batch_size), v%r(0:2)
        deallocate(A2)
     enddo scaling_loop

  end subroutine bench__dtt_tensor_2d


  subroutine bench__dtt_tensor_3d(batch_size, epsi)
  use time_lib
  use thor_lib, only: dtt_tensor
  implicit none
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  !
  double precision, allocatable :: A3(:,:,:)
  integer :: mode_size, n(3) !< mode size
  integer :: i, test
  type(dtt_tensor) :: v

     call system_timer_init()
     print '("# THOR test: dtt_tensor(A3,epsi)")'
     print '("# Modes = (M, 2M, 3M)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:time[s] 3-6:ranks")'
     scaling_loop: do mode_size = 10, 100, 10
        n = [mode_size, mode_size*2, mode_size*3]
        allocate(A3(n(1), n(2), n(3)))
        call random(A3)
        call system_timer_start
        test_loop: do test = 1, batch_size
           v = dtt_tensor(A3, eps_=epsi)
        enddo test_loop
        call system_timer_stop
        print '(I10, 1X, ES12.5, 14(I6,1X))', &
              mode_size, system_dt/dble(batch_size), v%r(0:3)
        deallocate(A3)
     enddo scaling_loop

  end subroutine bench__dtt_tensor_3d


  subroutine bench__dtt_tensor_4d(batch_size, epsi)
  use time_lib
  use thor_lib, only: dtt_tensor
  implicit none
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  !
  double precision, allocatable :: A4(:,:,:,:)
  integer :: mode_size, n(4) !< mode size
  integer :: i, test
  type(dtt_tensor) :: v

     call system_timer_init()
     print '("# THOR test: dtt_tensor(A4,epsi)")'
     print '("# Modes = (M, 2M, 3M, 4M)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:time[s] 3-7:ranks")'
     scaling_loop: do mode_size = 10, 50, 5
        n = [mode_size, mode_size*2, mode_size*3, mode_size*4]
        allocate(A4(n(1), n(2), n(3), n(4)))
        call random(A4)
        call system_timer_start
        test_loop: do test = 1, batch_size
           v = dtt_tensor(A4, eps_=epsi)
        enddo test_loop
        call system_timer_stop
        print '(I10, 1X, ES12.5, 14(I6,1X))', &
              mode_size, system_dt/dble(batch_size), v%r(0:4)
        deallocate(A4)
     enddo scaling_loop

  end subroutine bench__dtt_tensor_4d


  subroutine bench__dtt_tensor_lrank_4d(batch_size, epsi)
  use time_lib
  use thor_lib, only: dtt_tensor, dtt_random_ortho
  implicit none
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  !
  double precision, allocatable :: A4(:,:,:,:)
  integer :: mode_size, n(4), r(0:4) !< mode size
  integer :: i, test
  type(dtt_tensor) :: v, vr

     call system_timer_init()
     print '("# THOR test: dtt_tensor(A4,epsi) from low-rank tensor")'
     print '("# Modes = (M, 2M, 3M, 4M), ranks = (1, 5, 5, 5, 1)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:time[s] 3-7:ranks")'

     r =  [1, 5, 5, 5, 1]
     scaling_loop: do mode_size = 10, 50, 5
        n = [mode_size, mode_size*2, mode_size*3, mode_size*4]
        vr = dtt_random_ortho(n, r_=r)
        A4 = reshape(vr% full(), n)
        call system_timer_start
        test_loop: do test = 1, batch_size
           v = dtt_tensor(A4, eps_=epsi)
        enddo test_loop
        call system_timer_stop
        print '(I10, 1X, ES12.5, 14(I6,1X))', &
              mode_size, system_dt/dble(batch_size), v%r(0:4)
        deallocate(A4)
     enddo scaling_loop

  end subroutine bench__dtt_tensor_lrank_4d


  subroutine bench__dtt_tensor_mround_4d(batch_size, epsi)
  use time_lib
  use thor_lib
  implicit none
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  !
  double precision, allocatable :: A4(:,:,:,:)
  integer :: mode_size, n(4), r(0:4) !< mode size
  integer :: i, test
  type(dtt_tensor) :: v, vr, w

     call system_timer_init()
     print '("# THOR test: dtt_tensor(A4,epsi) from low-rank tensor")'
     print '("# Modes = (M, 2M, 3M, 4M), ranks = (1, 3, 9, 7, 1)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:time[s] 3-7:ranks")'

     r =  [1, 3, 9, 7, 1]
     scaling_loop: do mode_size = 10, 50, 5
        n = [mode_size, mode_size*2, mode_size*3, mode_size*4]
        vr = dtt_random_ortho(n, r_=r)
        v  = dtt_random_ortho(n, r_=r)
        w = v + vr
        call system_timer_start
        test_loop: do test = 1, batch_size
           call w% mround(epsi)
        enddo test_loop
        call system_timer_stop
        print '(I10, 1X, ES12.5, 14(I6,1X))', &
              mode_size, system_dt/dble(batch_size), w% r(0:4)
     enddo scaling_loop

  end subroutine bench__dtt_tensor_mround_4d

end program
