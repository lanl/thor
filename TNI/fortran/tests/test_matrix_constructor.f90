!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_thor_matrix_constructor.f90
!! @author [IDB], [OGK]
!! @date  October 2023
!! @brief Testing tt_matrix constructors
!!

program main
use rnd_lib, only: random, init_random_seed
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 4  ! tensor order
integer :: n = 3  ! tensor dimension size
integer :: r = 2  ! tensor rank
integer :: rseed = 0 ! random seed (0 = from clock)
integer, allocatable :: q(:), s(:), rr(:)
logical :: lrandomize_nr = .false.
double precision, allocatable :: rnd(:)

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
  allocate(q(d), s(d), rr(d+1), rnd(d))
  rr(1)= 1; rr(d+1)= 1
  if (lrandomize_nr) then
     call random(rnd)
     q= 1 + int(rnd*n)
     call random(rnd)
     s= 1 + int(rnd*n)
     call random(rnd)
     rr(2:d)= 1 + int(rnd(2:d)*r)
     do i=1,d-1
        rr(i)= min(rr(i), max(q(i)*s(i), q(i+1)*s(i+1)))
     enddo
  else
     q = n; s = n; rr = r
     rr(1) = 1; rr(d+1) = 1
  endif

  call test_matrix_constructors(q, s, rr)
  deallocate(q, s, rr, rnd)

contains

  subroutine help_msg()
  print '(12(/,A))', &
   "Test matrix constructors", &
   "Usage: ./test.exe [-h] [-d <ord>] [-n <dim>] [-r <rank>]", &
   "                  [-l <num>] [-m <num>] [--randomize-nr]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : matrix order <ord> [4]", &
   " -n <dim>  : max dimension for each mode [3]", &
   " -r <rank> : rank of the cores [2]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]"
  end subroutine


  subroutine test_matrix_constructors(q, s, r)
  use time_lib
  use rnd_lib, only: rndmat
  use mat_lib, only: normfro
  use thor_lib, only: dtt_matrix, dtt_matrix_zeros, dtt_matrix_zeros_like, &
      dtt_matrix_ones, dtt_matrix_ones_like
  use matrix_util, only: unravel_index, ravel_multi_index
  use matlab_struct_module, only : pprint_matrix
  implicit none
  integer, intent(in) :: q(:), s(:), r(:)
  !
  type(dtt_matrix) :: A, B
  double precision :: nrm, tol
  double precision, allocatable :: A2(:,:), B2(:,:), A1(:)
  integer :: failed, MA, NA
  integer, allocatable :: qq(:), ss(:)

     ! initialize system timer
     call system_timer_init()
     failed = 0
     tol = 1d-12

     print '(/"Test 07: Testing dtt_matrix constructor methods"/)'

     print '("-- [7-1] testing: A = dtt_matrix(q, s, r_=r, l_=1)")'
     call system_timer_start
     A = dtt_matrix(q, s, r_=r, l_=1)
     call system_timer_stop
     call A% say
     nrm = A% normb()
     print '("|A|_2 = "F3.1)', nrm
     write(*, '("[+][TEST07-1]["F7.2" s]:&
        &matrix A(q,s,r) allocated OK and has zero norm: ")', &
        advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif
     
     print '("-- [7-2] testing: A = dtt_matrix(s, q)")'
     call system_timer_start
     A = dtt_matrix(s, q)
     call system_timer_stop
     call A% say
     nrm = A% normb()
     print '("|A|_2 = "F3.1)', nrm
     write (*,'("[+][TEST07-2]["F7.2" s]: matrix A(s,q) allocated: ")' &
             , advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif
     
     print '("-- [7-3] testing: A = dtt_matrix(s, q, l_=1, m_=size(s))")'
     call system_timer_start
     A = dtt_matrix(s, q, l_=1, m_=size(s))
     call system_timer_stop
     call A% say
     nrm = A% normb()
     print '("|A|_2 = "F3.1)', nrm
     write (*,'("[+][TEST07-3]["F7.2" s]: matrix A(s,q,l_,m_) allocated:")' &
             , advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [7-4] testing deep copy constructor: B = dtt_matrix(A)")'
     call system_timer_start
     B = dtt_matrix(A)
     call system_timer_stop
     call B% say
     nrm = B% normb()
     print '("|B|_2 = "F3.1)', nrm
     write (*,'("[+][TEST07-4]["F7.2" s]: deep copy constructor B(<-A): ")' &
             , advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [7-5] testing zeros: A = dtt_matrix(q, s, r_=r)")'
     call system_timer_start
     A = dtt_matrix_zeros(q, s, r_=r)
     call system_timer_stop
     call A% say
     nrm = A% normb()
     print '("|A|_2 = "F3.1)', nrm
     write (*,'("[+][TEST07-5]["F7.2" s]: matrix A(q,s,r) with zeros: ")' &
             , advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [7-6] testing zeros: B = dtt_matrix_zeros_like(A)")'
     call system_timer_start
     B = dtt_matrix_zeros_like(A)
     call system_timer_stop
     call B% say
     nrm = B% normb()
     print '("|B|_2 = "F3.1)', nrm
     write (*,'("[+][TEST07-6]["F7.2" s]: matrix B of zeros like A created: ")' &
             , advance='no') system_dt
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [7-7] testing zeros: A = dtt_matrix_ones(q, s)")'
     call system_timer_start
     A = dtt_matrix_ones(q, s)
     call system_timer_stop
     call A% say
     nrm = A% normb()
     print '("|A|_2 = "F3.1)', nrm
     write (*,'("[+][TEST07-7]["F7.2" s]: matrix ones A(q,s,r) created: ")' &
             , advance='no') system_dt
     if (dabs(nrm - 1d0).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [7-8] testing zeros: B = dtt_matrix_ones_like(A)")'
     call system_timer_start
     B = dtt_matrix_ones_like(A)
     call system_timer_stop
     call B% say
     nrm = B% normb()
     print '("|B|_2 = "F3.1)', nrm
     write (*, '("[+][TEST07-8]["F7.2" s]: matrix B of ones like A created: ")' &
             , advance='no') system_dt
     if (dabs(nrm - 1d0).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [7-9] testing matrix constructor from a 2D array")'
     qq = [4, 6, 3]
     ss = [17, 3, 16]
     MA = product(qq); NA = product(ss)
     print '("   Matrix A2 before tensorization:")'
     A2 = rndmat(MA,NA)
     call pprint_matrix(A2, line_len_=160)
     A = dtt_matrix(A2, qq, ss)
     B2 = A% full2d() 
     nrm = normfro(B2 - A2)
     write (*, '("[+][TEST07-9]: difference after tensorization = "ES9.2": ")' &
             , advance='no') nrm
     if (dabs(nrm).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [7-10] experimenting with rank-1 subsystems")'
     qq = [2, 2, 3]
     ss = [3, 2, 3]
     A  = dtt_matrix_zeros(qq, ss) ! rank 1
     A% u4(1)% p(1,2,3,1) = 1d0
     A% u4(2)% p(1,:,:,1) = 1d0
     A% u4(3)% p(1,:,:,1) = 11d0
     A2 = A% full2d() 
     print '("   Matrix in full form:")'
     call pprint_matrix(A2, frmt_='(F3.0)', max_rows_=14, line_len_=160)

     if(failed.eq.0) then
        print '("[+][TEST07]: PASSED")'
     else
        print '("[-][TEST07]: FAILED with ",I2," error(s)")', failed
        error stop failed
     end if

  end subroutine test_matrix_constructors

end program 
