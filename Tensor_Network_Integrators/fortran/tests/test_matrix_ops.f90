!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_thor_matrix_constructor.f90
!! @author [IDB], [OGK]
!! @date  October 2023
!! @brief Testing dtt_matrix constructors
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
integer :: l = 1  ! starting index
integer :: m      ! ending index
integer :: rseed = 0 ! random seed (0 = from clock)
integer, allocatable :: q(:), s(:), rr(:)
logical :: lrandomize_nr = .false.
logical :: m_given = .false.
double precision, allocatable :: rnd(:)

  m = d
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
     elseif (trim(arg).eq.'-l') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) l
        if (err.ne.0) error stop "ERROR: when parsing -l <???>"
        i= i + 1
     elseif (trim(arg).eq.'-m') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) m
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
        m_given = .true.
        i= i + 1
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     elseif (trim(arg).eq.'--randomize-nr') then
        lrandomize_nr= .true.
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

  if (.not.m_given) m = d
  !call test_tt_matrix_ops(q, s, rr, l, m)
  call test_dttm_add(q, s, rr)
  deallocate(q, s, rr, rnd)

contains

  subroutine help_msg()
  print '(11(/,A))', &
   "Test matrix operations", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : matrix order <ord> [4]", &
   " -n <dim>  : max dimension for each mode [3]", &
   " -r <rank> : rank of the cores [2]", &
   " -l <num>  : starting index [1]", &
   " -m <num>  : ending index [4]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]"
  end subroutine


  subroutine test_tt_matrix_ops(q, s, r, l, m)
  use time_lib
  use thor_lib
  use rnd_lib, only : init_random_seed
  implicit none
  integer, intent(in) :: q(:), s(:), r(:), l, m
  !
  integer:: sa, failed
  double precision, allocatable:: a1(:), a2(:)
  type(dtt_matrix) :: A, B, C, dM
  type(dtt_tensor) :: v
  double precision :: err, lg_ftsize
  double precision, parameter :: tol = 1d-14
  double precision, external  :: dnrm2
  double precision, parameter :: lg_ftsize_max = 9d0

     ! initialize system timer
     call system_timer_init()

     print '(/"Test 08: Testing thor matrix operations"/)'
     failed = 0
     call init_random_seed(1525)

     print '("-- [8-1] testing: A = dtt_matrix_rand(..)")'
     call system_timer_start
     A = dtt_matrix_rand(q, s, l_=1, r_=r)
     call system_timer_stop
     print '("matrix A:")'
     call A% say()
     lg_ftsize = A% log10_fullsize()
     if (lg_ftsize.lt.lg_ftsize_max) then
        a1 = A% full()
     endif
     write (*,'("[+][TEST08-1]["F7.2" s]:&
        & random matrix A created: PASS"/)') system_dt

     print '("-- [8-2] testing: A -> v conversion")'
     call system_timer_start
     v = dtt_tensor(A)
     call system_timer_stop
     print '("tensor v:")'
     call v% say()
     if (lg_ftsize.lt.lg_ftsize_max) then
        a2 = v% full()
        sa = size(a2)
        err= dnrm2(sa, a1 - a2, 1)
        print '("|A - v|_2 = ",ES12.5)', err
     else
        dM= A - dtt_matrix(v, q, s)
        err= dM% normb()
        print '("|A - M(v)|_2 = ",ES12.5)', err
     endif
     write (*,'("[+][TEST08-2]["F7.2" s]:&
        & A -> v conversion correct: ")', &
        advance='no') system_dt
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [8-3] testing: v -> A conversion")'
     call system_timer_start
     B = dtt_matrix(v, q, s)
     call system_timer_stop
     print '("matrix B:")'
     call B% say()
     if (lg_ftsize.lt.lg_ftsize_max) then
        a2 = B% full()
        err= dnrm2(sa, a1 - a2, 1)
        print '("|A - B|_2 = ",ES12.5)', err
     else
        v = v - dtt_tensor(B)
        err= v% normb()
        print '("|v - v(M(v))|_2 = ",ES12.5)', err
     endif
     write (*,'("[+][TEST08-3]["F7.2" s]:&
        & v -> A conversion correct: ")', &
        advance='no') system_dt
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [8-4] testing: C = A - B computation")'
     !B = A + (B - A)
     call system_timer_start
     C = A - B
     call system_timer_stop
     print '("matrix C:")'
     call C% say()
     if (lg_ftsize.lt.lg_ftsize_max) then
        a2 = C% full()
        err= dnrm2(sa, a2, 1)
        print '("|B - A|_2 = ",ES12.5)', err
     else
        err= C% normb()
        print '("|B - A|_2 = ",ES12.5)', err
     endif
     write (*, '("[+][TEST8]["F7.2" s]:&
        & C = (B - A) == 0 computed correctly: ")', &
        advance='no') system_dt
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     if(failed.eq.0) then
        print '("[+][TEST8]: PASSED")'
     else
        print '("[-][TEST8]: FAILED with ",I2," error(s)")', failed
        error stop failed
     end if
  end subroutine test_tt_matrix_ops


  subroutine test_dttm_add(q, s, r)
  use time_lib
  use thor_lib
  use rnd_lib, only : init_random_seed
  implicit none
  integer, intent(in) :: q(:), s(:), r(:)
  !
  integer:: sa, failed
  double precision, allocatable:: a1(:), a2(:)
  type(dtt_matrix) :: A, B, C, C1, dC
  type(dtt_tensor) :: u, v, w, w1, dw
  double precision, parameter :: tol = 1d-14

     print '(/"Testing dttm addition bug"/)'
     failed = 0
     call init_random_seed(1525)

     print '("-- creating matrices A, B:")'
     A = dtt_matrix_rand(q, s, r_=r)
     print '("matrix A:")'
     call A% say()

     B = dtt_matrix_rand(q, s, r_=r)
     print '(/"matrix B:")'
     call B% say()

     print '(/"-- checking A + B == B + A:")'
     C= A + B
     print '("C = |A + B|_2 = ", ES14.5)', C% norm()
     call C% say()
     print '("C% full() = ", 16(F5.3,1X))', C% full()

     C1= B + A
     print '(/"C1= |B + A|_2 = ", ES14.5)', C1% norm()
     call C1% say()
     print '("C1%full() = ", 16(F5.3,1X))', C1% full()
     dC= C - C1
     call dC% pprint(label_="dtt_matrix dC:", line_len_=400)

     print '("|C - C1|_2 = ", ES14.5)', dC% norm()

     return

     print '("-- comparing with dtt matrices converted to dtt vectors:")'
     u= dtt_tensor(A)
     v= dtt_tensor(B)
     w= u + v
     call w% pprint(label_="dtt_tensor w:", line_len_=400)
     call C%  pprint(label_="dtt_matrix C:", line_len_=400)
     call C1% pprint(label_="dtt_matrix C1:", line_len_=400)

     print '("-- creating vectors u, v:")'
     u = dtt_tensor_rand(q, r_=r)
     print '("tensor u:")'
     call u% say()

     v = dtt_tensor_rand(q, r_=r)
     print '(/"matrix v:")'
     call v% say()

     print '(/"-- checking u + v == v + u:")'
     w= u + v
     print '("w = |u + v|_2 = ", ES14.5)', w% norm()
     call w% say()
     print '("w% full() = ", 4(F5.3,1X))', w% full()

     w1= v + u
     print '(/"w1= |v + u|_2 = ", ES14.5)', w1% norm()
     call w1% say()
     print '("C1%full() = ", 16(F5.3,1X))', w1% full()
     dw= w - w1
     print '("|w - w1|_2 = ", ES14.5)', dw% norm()

  end subroutine test_dttm_add

end program
