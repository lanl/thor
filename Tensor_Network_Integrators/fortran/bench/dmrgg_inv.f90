!> @file dmrgg_inv.f90
!!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @author Oleg Korobkin (korobkin@lanl.gov)
!! @date   February 2025
!! @brief  Benchmarking tensor inversion using DMRG
!!

program main
use rnd_lib, only: random, init_random_seed
use thor_lib, only: dtt_tensor, dtt_tensor_ones, tijk, &
                    dtt_tensor_rand
use dmrgg_lib, only: dtt_dmrgg
implicit none
character(len=128) arg
!
integer :: i, err
double precision  :: epsi = 1d-8
integer           :: batch_size = 10
character(len=20) :: which_test
integer           :: low_rank = 5
integer           :: d = 4
type(dtt_tensor)  :: xg

  which_test = "sin"

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-e') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) epsi
        if (err.ne.0) error stop "ERROR: when parsing -e <???>"
     elseif (trim(arg).eq.'-b') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) batch_size
        if (err.ne.0) error stop "ERROR: when parsing -b <???>"
     elseif (trim(arg).eq.'-l') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) low_rank
        if (err.ne.0) error stop "ERROR: when parsing -l <???>"
     elseif (trim(arg).eq.'-d') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) d
        if (err.ne.0) error stop "ERROR: when parsing -d <???>"
     elseif (trim(arg).eq.'-h') then
        call help_msg
        stop
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call bench_dmrgg_inv(d, low_rank, epsi, batch_size)

contains

  subroutine help_msg()
  print '(100(/,A))', &
   "Benchmark dtt_tensor constructor: from various n-dim arrays", &
   "Usage: ./dtt_tensor.exe [-h]", &
   "Options:", &
   " -d <val>   : tensor order (dimensionality) [4]", &
   " -l <val>   : maximal rank of the generated tensor [5]", &
   " -e <val>   : rounding tolerance epsilon [1e-8]", &
   " -b <val>   : batch size [10]", &
   " -h         : print this help message"
  end subroutine


  subroutine bench_dmrgg_inv(d, low_rank, epsi, batch_size)
  use time_lib
  use dmrgg_interpolate_lib, only    : inv_tens_d
  use matlab_struct_module, only: pprint_matrix
  implicit none
  include 'mpif.h'
  integer, intent(IN) :: d          !< tensor order
  integer, intent(IN) :: low_rank   !< maximal rank of the generated tensor
  integer, intent(IN) :: batch_size
  double precision, intent(IN) :: epsi
  !
  type(dtt_tensor)    :: y
  integer :: mode_size, test, info, nproc, me, failed
  integer, allocatable :: nn(:), rr(:)
  double precision, allocatable      :: arr1(:)
  double precision                   :: eps, err
  double precision :: dt_inv_tens, dt_dtt_dmrgg

     call mpi_init(info)
     call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
     call mpi_comm_rank(MPI_COMM_WORLD,me,info)
     call system_timer_init

     !if (me.eq.0) print '("[MPI]: initialized OK: size = ",I5)', nproc
     call mpi_barrier(  MPI_COMM_WORLD, info)
     allocate(nn(d),rr(0:d))
     rr(0)= 1; rr(1:d-1)= low_rank; rr(d)= 1

     call system_timer_init()
     print '("# THOR DMRGG inversion test: inv(A2,epsi)")'
     print '("# Tensor order = "I5", low rank = "I7)', d, low_rank
     print '("# Modes = (M, M, ... M)")'
     print '("# Rounding tolerance: epsi = ",ES12.5)', epsi
     print '("# 1:mode_size 2:dmrg-inv[s]")'
     mode_size = 4
     scaling_loop: do i = 1,14 !mode_size = 10, 100, 10
        nn(1:d)= mode_size
        xg = dtt_tensor_rand(nn, r_=rr)

        !call system_timer_start
        !test_loop: do test = 1, batch_size
        !   y = inv_tens_d(x, accuracy_=eps)
        !enddo test_loop
        !call system_timer_stop
        !dt_inv_tens = system_dt/dble(batch_size)

        y = dtt_tensor_ones(nn)
        call system_timer_start
        test_loop: do test = 1, batch_size
           call dtt_dmrgg(y, xfun, accuracy=eps)
        enddo test_loop
        call system_timer_stop
        dt_dtt_dmrgg = system_dt/dble(batch_size)
        print '(I10, 2(1X, ES12.5),10(1X,I4))', &
              mode_size, dt_inv_tens, dt_dtt_dmrgg, y%r(0:d)
        mode_size = 2*mode_size
     enddo scaling_loop
     call mpi_finalize(info)

  end subroutine bench_dmrgg_inv


  pure function xfun(d, ind, nn) result(y)
  use thor_lib, only: tt_size
  implicit none
  double precision :: y
  integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
  !double precision, intent(IN), optional :: par(*)
  !
  integer :: l, m
  double precision :: x(tt_size)

    l= xg%l; m= xg% m
    y = -1.37707d0; if (m - l + 1.ne.d) return
    y = -2.37707d0; if (size(ind).ne.tt_size) return
    y = -3.37707d0; if (size(nn).ne.tt_size) return
    y = -4.37707d0; if (any(ind(l:m)<1).or.any(ind(l:m)>nn(l:m))) return
    
    x(l:m)= dble(ind(l:m))/dble(nn(l:m))
    !y = dsin(103*x(l)) + 0.3d0*dcos(20*x(m))
    y = 1d0/(1d-1 + sum(x(l:m)**2))

  end function


  pure function xfun_tt(d, ind, nn) result(y)
  use thor_lib, only: tt_size
  implicit none
  double precision :: y
  integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
  !double precision, intent(IN), optional :: par(*)
  !
  integer :: l, m

    l= xg%l; m= xg% m
    y = -1.37707d0; if (m - l + 1.ne.d) return
    y = -2.37707d0; if (size(ind).ne.tt_size) return
    y = -3.37707d0; if (size(nn).ne.tt_size) return
    y = -4.37707d0; if (any(ind(l:m)<1).or.any(ind(l:m)>nn(l:m))) return
    y = tijk(xg, ind(l:m))
    if (y.gt.0d0) then
       y = 1d0/(y + 1d-15)
    else
       y = 1d0/(y - 1d-15)
    endif

  end function

end program 
