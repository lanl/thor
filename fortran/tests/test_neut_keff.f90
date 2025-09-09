!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_neut_keff.f90
!! @author [IDB], [OGK]
!! @date  October 2023
!! @brief Test neutron criticality calculation
!!
program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg

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

  call test_neut_keff

contains

  subroutine help_msg
  print '(3(/,A))', &
   "Test (DT) using the matrices from neutron transport problem", &
   "To run, launch Thor_1DSlab_Pu239.m in tests/testdata first", &
   "Usage: ./test.exe [-h]"
  end subroutine


  subroutine test_neut_keff
  use time_lib
  use matlab_struct_module, only : array3d
  use thorio_lib, only: dtt_read_sdv_file
  use thor_lib
  use ttamen_lib
  implicit none
  !
  type(dtt_tensor) :: qttPsi0, RHS, qttPsi1, tx, dtx
  type(dtt_matrix) :: qttH, qttHs, qttHf
  integer :: niter, iter
  integer, allocatable :: nn(:)
  type(array3d), allocatable :: tz_c(:)
  double precision :: tol, fixed_point_tol, k, convrate
  character(len=64):: fname

     niter = 100
     tol = 1d-6
     fixed_point_tol = 1d-6
     k = 1d0

     allocate(nn(15))
     nn = 2

     fname = 'tests/testdata/qttH.sdv'
     call check_if_exists(fname)
     tx= dtt_read_sdv_file(fname)
     qttH = dtt_matrix(tx, nn, nn)
     print '("qttH: ")'
     call say(qttH)

     fname= 'tests/testdata/qttHs.sdv'
     call check_if_exists(fname)
     tx= dtt_read_sdv_file(fname)
     qttHs = dtt_matrix(tx, nn, nn)
     print '("qttHs: ")'
     call say(qttHs)
     
     fname= 'tests/testdata/qttHf.sdv'
     call check_if_exists(fname)
     tx= dtt_read_sdv_file(fname)
     qttHf = dtt_matrix(tx, nn, nn)
     print '("qttHf: ")'
     call say(qttHf)

     fname= 'tests/testdata/qttPsi0.sdv'
     call check_if_exists(fname)
     qttPsi0 = dtt_read_sdv_file(fname)
     print '("qttPsi0: ")'
     call say(qttPsi0)

     ! fixed-point iteration to find k-effective
     fixed_point_loop: do iter = 1, niter
       print '("Iteration = ",I3,"/",I3)',iter,niter

       if (allocated(tz_c)) deallocate(tz_c)
       allocate(tz_c(0))

       RHS =       amen_mv(core2cell(qttHs), core2cell(qttPsi0), tz_c, tol) &
           + 1d0/k*amen_mv(core2cell(qttHf), core2cell(qttPsi0), tz_c, tol)

       ! solve a linear system for Psi1 = Htt\RHS;
       qttPsi1 = qttPsi0
       call dtt_amen_solve(qttH, RHS, tol, qttPsi1)

       ! update eigenvalue
       k = k*sumall(amen_mv(core2cell(qttHf), core2cell(qttPsi1), tz_c, tol)) &
            /sumall(amen_mv(core2cell(qttHf), core2cell(qttPsi0), tz_c, tol))
       print '("iter = "I3", k_eff = ", F10.6)', iter, k
      
       dtx = qttPsi1 - qttPsi0
       convrate = sqrt(dabs(dtx% normb()/qttPsi0% normb()))
       if (convrate .lt. fixed_point_tol) then
         print '("The fixed point scheme converged at iter = "I3)', iter
         print '("k-effective = "F10.6)', k
         exit fixed_point_loop
       endif
       qttPsi0 = qttPsi1

     enddo fixed_point_loop

     print '("qttPsi_solution (Thor): ")'
     call say(qttPsi1)

     fname= 'tests/testdata/qttPsi_solution.sdv'
     call check_if_exists(fname)
     qttPsi0= dtt_read_sdv_file(fname)
     print '("qttPsi_solution (Matlab): ")'
     call say(qttPsi0)
     dtx = qttPsi1 - qttPsi0
     print '("difference with Matlab solution: "ES14.7)', &
           dsqrt(dabs(dtx% normb()/qttPsi0% normb()))
  end subroutine test_neut_keff

  subroutine check_if_exists(fname)
  character(*), intent(in) :: fname
  logical :: file_exists

     inquire(file=trim(fname), exist=file_exists)
     if (file_exists) return
     print '("ERROR: File ",A," does not exist")',trim(fname)
     print '("Run script Thor_1DSlab_Pu239.m in tests/testdata first")'
     error stop
  end subroutine check_if_exists

end program
