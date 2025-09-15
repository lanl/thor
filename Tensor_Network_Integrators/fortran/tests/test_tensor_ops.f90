!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_thor_tensor_ops.f90
!! @author [IDB], [OGK]
!! @date  March 2024
!! @brief Testing tensor operations
!!

program main
use rnd_lib, only: random, init_random_seed
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 3  ! tensor order
integer :: n = 4  ! tensor dimension size
integer :: r = 2  ! tensor rank
integer :: l = 1  ! starting index
integer :: m      ! ending index
integer :: rseed = 0 ! random seed (0 = from clock)
integer, allocatable :: nn(:), rr(:)
logical :: lrandomize_nr = .false.
logical :: m_given = .false.
double precision, allocatable :: rnd(:)
double precision :: tol = 1d-12

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
     elseif (trim(arg).eq.'-t') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) tol
        if (err.ne.0) error stop "ERROR: when parsing -t <???>"
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
  allocate(nn(d), rr(d+1), rnd(d))
  rr(1)= 1; rr(d+1)= 1
  if (lrandomize_nr) then
     call random(rnd)
     nn= 1 + int(rnd*n)
     call random(rnd)
     rr(2:d)= 1 + int(rnd(2:d)*r)
     do i=1,d-1
        rr(i)= min(rr(i), max(nn(i),nn(i+1)))
     enddo
  else
     nn= n; rr(2:d)= r
  endif
  if (.not.m_given) m = d
  call test_tensor_ops(nn, rr, tol, l, m)
  deallocate(nn, rr, rnd)


contains

  subroutine help_msg()
  print '(10(/,A))', &
   "Test of tensor operations", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : tensor order <ord> [3]", &
   " -n <dim>  : tensor dimensions <dim> [4]", &
   " -r <rank> : rank of the cores [2]", &
   " -l <num>  : starting index [1]", &
   " -m <num>  : ending index [d]", &
   " -t <tol>  : use this tolerance in SVD step [1e-12]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]"
  end subroutine

  subroutine test_tensor_ops(n, r, tol, l, m)
  use time_lib
  use thor_lib
!  use thor_lib, only : dtt_tensor, dtt_tensor_ones, alloc, dealloc, &
!      dtt_tensor_rand, dtt_full, sayfull, pprint, &
!      operator(*), operator(+), operator(-), operator(/)
  use rnd_lib, only : random
  use thorio_lib
  implicit none
  integer, intent(in) :: n(:), r(:), l, m
  double precision, intent(in) :: tol
  include 'mpif.h'
  type(dtt_tensor) :: v1, v2, v3, v4, v5, dv, Z
  double precision :: err, lg_ftsize, alpha, dt1, dt2
  double precision, allocatable      :: a1(:), a2(:)
  integer                            :: info, nproc, me, failed, sa
  double precision, external :: dnrm2
  double precision, parameter :: lg_ftsize_max = 9d0
  double precision, allocatable :: a3(:,:,:), a4(:,:,:,:), a5(:,:,:,:,:), &
                                   a6(:,:,:,:,:,:), a7(:,:,:,:,:,:,:)

     ! initialize system timer
     call system_timer_init()

     print '(/"Test 05: Thor tensor operations"/)'
     failed= 0
     alpha= 2d0

     call mpi_init(info)
     if(info.ne.0) error stop '[MPI]: init fail'
     call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
     if(info.ne.0) error stop '[MPI]: comm_size fail'
     call mpi_comm_rank(MPI_COMM_WORLD,me,info)
     if(info.ne.0) error stop '[MPI]: comm_rank fail'

     if(me .eq. 0) print '("[MPI]: initialized OK: size = ",I5)', nproc
     call MPI_Barrier(  MPI_COMM_WORLD, info)
     if(info.ne.0) error stop '[MPI]: Barrier fail'

     print '("-- [5-1] testing:  v1 = dtt_tensor_ones(..):")'
     call system_timer_start
     v1 = dtt_tensor_ones(n)
     call system_timer_stop
     dt1= system_dt
     call v1% say()
     lg_ftsize = v1% log10_fullsize()

     if (lg_ftsize.lt.lg_ftsize_max) then
        !call sayfull(v1)
        a1 = v1% full()
        sa = size(a1)
        call system_timer_start
        err= dnrm2(sa, a1 - 1d0, 1)
        call system_timer_stop
        print '("|v1.full() - 1|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     else
        call system_timer_start
        dv = v1 - 1d0
        err= dv% normb()
        call system_timer_stop
        print '("|v1 - 1|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     endif
     write (*,'("[+][TEST05-1]["F7.2" s]:&
        & tensor of ones v1 created (tol=",ES8.1,"): ")',&
        advance='no') dt1, tol
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [5-2] testing:  v2 = alpha * v1")'
     call system_timer_start
     v2 = alpha * v1
     call system_timer_stop
     dt1= system_dt
     call v2% say()
     err= .0d0
     if (lg_ftsize.lt.lg_ftsize_max) then
        !call sayfull(v2)
        sa = size(v2%full())
        call system_timer_start
        err= dnrm2(sa, v2%full() - alpha, 1)
        call system_timer_stop
        print '("|v2.full() - a|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     else
        call system_timer_start
        dv = v2 - alpha
        err= dv% normb()
        call system_timer_stop
        print '("|v2 - a|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     endif
     write (*,'("[+][TEST05-2]["F7.2" s]: &
        &v2 = alpha * v1 computed correctly (tol="ES8.1"): ")', &
        advance='no') dt1, tol
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [5-2-3] testing:  v3 = random 3D tensor and its SVD rep")'
     allocate(a3(n(1),n(2),n(3)))
     print *, "shape of a3 = ", shape(a3)
     call random(a3)
     v3 = dtt_from_3darray(a3)
     call v3% say()
     err= dnrm2(size(a3), v3%full() - reshape(a3, [size(a3)]), 1)
     print '("|v3.full() - a3|_2 = ",ES12.5)', err
     if (err.lt.tol) then
        print '("[+][TEST04]: v3 = random tensor from 3D array created: PASS")'
     else
        print '("[-][TEST04]: v3 = random tensor from 3D array incorrect: FAIL")'
        failed= failed + 1
     end if
     deallocate(a3)

     print '("-- [5-2-4] testing:  v3 = random 4D tensor and its SVD rep")'
     allocate(a4(4,4,4,4))
     call random(a4)
     v3 = dtt_tensor(a4)
     call v3% say()
     err= dnrm2(size(a4), v3%full() - reshape(a4, [size(a4)]), 1)
     print '("|v3.full() - a4|_2 = ",ES12.5)', err
     if (err.lt.tol) then
        print '("[+][TEST04]: v3 = random tensor from 4D array created: PASS")'
     else
        print '("[-][TEST04]: v3 = random tensor from 4D array incorrect: FAIL")'
        failed= failed + 1
     end if
     deallocate(a4)

     print '("-- [5-2-5] testing:  v3 = random 5D tensor and its SVD rep")'
     allocate(a5(4,4,4,4,4))
     call random(a5)
     v3 = dtt_tensor(a5)
     call v3% say()
     err= dnrm2(size(a5), v3%full() - reshape(a5, [size(a5)]), 1)
     print '("|v3.full() - a5|_2 = ",ES12.5)', err
     if (err.lt.tol) then
        print '("[+][TEST04]: v3 = random tensor from 5D array created: PASS")'
     else
        print '("[-][TEST04]: v3 = random tensor from 5D array incorrect: FAIL")'
        failed= failed + 1
     end if
     deallocate(a5)

     print '("-- [5-2-6] testing:  v3 = random 6D tensor and its SVD rep")'
     allocate(a6(3,3,3,3,3,3))
     call random(a6)
     v3 = dtt_tensor(a6)
     call v3% say()
     err= dnrm2(size(a6), v3%full() - reshape(a6, [size(a6)]), 1)
     print '("|v3.full() - a6|_2 = ",ES12.5)', err
     if (err.lt.tol) then
        print '("[+][TEST04]: v3 = random tensor from 6D array created: PASS")'
     else
        print '("[-][TEST04]: v3 = random tensor from 6D array incorrect: FAIL")'
        failed= failed + 1
     end if
     deallocate(a6)

     print '("-- [5-2-7] testing:  v3 = random 7D tensor and its SVD rep")'
     allocate(a7(2,2,2,2,2,2,2))
     call random(a7)
     v3 = dtt_tensor(a7)
     call v3% say()
     err= dnrm2(size(a7), v3%full() - reshape(a7, [size(a7)]), 1)
     print '("|v3.full() - a7|_2 = ",ES12.5)', err
     if (err.lt.tol) then
        print '("[+][TEST04]: v3 = random tensor from 7D array created: PASS")'
     else
        print '("[-][TEST04]: v3 = random tensor from 7D array incorrect: FAIL")'
        failed= failed + 1
     end if
     deallocate(a7)

     print '("-- [5-3] testing:  v3 = random tensor from flattened n-array, and its SVD rep")'
     call random(a1)
     v3 = dtt_tensor(a1, n)
     call v3% say()

     if (lg_ftsize.lt.lg_ftsize_max) then
        ! convert it to a full form and reconstruct using svd
        a1 = v3% full()
        call system_timer_start
        v3 = dtt_tensor(a1, n, eps_=tol)
        call system_timer_stop
        print '("reconstructed from full array using SVD, n(:) and &
            &max(r) ["F7.2" s]:")', system_dt
        call v3% say()
        !call sayfull(v3)
        call system_timer_start
        sa = size(a1)
        err= dnrm2(sa, v3%full() - a1, 1)
        call system_timer_stop
        print '("|v3.full() - a1|_2 = ",ES12.5," ["F9.4" s]")', err, system_dt
     else
        v1 = v3
        call system_timer_start
        call v3% mround(tol_=tol)
        call system_timer_stop
        dv = v3 - v1
        err= dv% normb()
        print '("|v3.mround() - v3|_2 = ",ES10.3," ["F9.4" s]")', err, system_dt
        print '("tensor v3 after rounding with mround(tol="ES8.1"):")', tol
        call v3% say()
     endif
     write (*,'("[+][TEST05-3]["F7.2" s]:&
        & v3 = random tensor created (tol="ES8.1"): ")',&
        advance='no') dt1, tol
     if (dabs(err).lt.1e-10) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif
     !call dtt_write_ascii_file(v1, '/tmp/v1.dat', verb_=.true.)

     print '("-- [5-4] testing:  v4 = random tensor and its SVD rep")'
     call system_timer_start
     v4 = dtt_tensor_rand(n,r)
     call system_timer_stop
     dt1= system_dt
     print '("original random tensor with prescribed n(:) and r(:):")'
     call v4% say()
     if (lg_ftsize.lt.lg_ftsize_max) then
        ! convert it to a full form and reconstruct using svd
        call system_timer_start
        a2 = v4% full()
        v4 = dtt_tensor(a2, n, rmax_=maxval(r), eps_=tol)
        call system_timer_stop
        print '("reconstructed from full array using SVD, n(:) and &
            &max(r) ["F7.2" s]:")', system_dt
        call v4% say()
        !call sayfull(v3)
        call system_timer_start
        sa = size(a2)
        err= dnrm2(sa, v4%full() - a2, 1)
        call system_timer_stop
        print '("|v4.full() - a2|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     else
        call system_timer_start
        v2 = v4
        call v4% svd(rmax_=maxval(r), tol_=tol)
        dv = v4 - v2
        err= dv% normb()
        call system_timer_stop
        print '("|v4.svd() - v4|_2 = ",ES10.3," ["F7.2" s]")',err,system_dt
        print '("tensor v4 after rounding with mround(tol="ES8.1"):")', tol
        call v4% say()
     endif
     write (*,'("[+][TEST05-4]["F7.2" s]:&
        & v4 = random tensor created (tol="ES8.1"): ")', &
        advance='no') dt1, tol
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [5-5] testing summation / difference:  v5 = v3 + v4")'
     call system_timer_start
     v5 = v3 + v4
     call system_timer_stop
     dt1= system_dt
     if (lg_ftsize.lt.lg_ftsize_max) then
        !call sayfull(v5)
        sa = size(a1)
        call system_timer_start
        err= dnrm2(sa, a1 + a2 - v5% full(), 1)
        call system_timer_stop
        print '("|v5 - (v3 + v4)|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     else
        call v5% say()
        dv = (v5 - v4) - v3
        print '("tensor dv = (v5 - v4) - v3: ")'
        call dv% say()
        err= dv% normb()
        print '("|(v5 - v4) - v3|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     endif
     write (*,'("[+][TEST05-5][",F7.2," s]:&
        & v5 = v3 + v4 computed correctly (tol=",ES8.1,"): ")', &
        advance='no') dt1, tol
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif


     print '("-- [5-6] testing multiplication / division:  v5 = v2*v4")'
     print '("         (with v2 = alpha * dtt_ones(n))")'
     call system_timer_start
     v2 = alpha*dtt_tensor_ones(n)
     v5 = v2 * v4
     call system_timer_stop
     dt1= system_dt
     if (lg_ftsize.lt.lg_ftsize_max) then
        !call sayfull(v5)
        sa = size(a1)
        call system_timer_start
        err= dnrm2(sa, alpha*a2 - v5% full(), 1)
        call system_timer_stop
        print '("|v5 - v2*v4|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     else
        call v5% say()
        dv = v5/alpha - v4
        print '("tensor dv = v5/alpha - v4: ")'
        call dv% say()
        err= dv% normb()
        print '("|v5/alpha - v4|_2 = ",ES12.5," ["F7.2" s]")', err, system_dt
     endif
     write (*,'("[+][TEST05-6][",F7.2," s]:&
        & v5 = v2 * v4 computed correctly (tol=",ES8.1,"): ")', &
        advance='no') dt1, tol
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     print '("-- [5-7] testing:  v4**2 - v3**2 = (v4 - v3)*(v4 + v3)")'
     call system_timer_start
     v1 = v4*v4
     v1 = v1% round()
     v2 = v3*v3
     call v1% svd(tol_=1d-12)
     call v2% svd(tol_=1d-12)
     dv = v1 - v2
     v1 = (v4 + v3)*(v4 - v3)
     call v1% svd(tol_=1d-12)
     dv = dv - v1
     call dv% svd(tol_=1d-12)
     !dv = (v4*v4 - v3*v3) - ((v4 - v3)*(v4 + v3))
     call system_timer_stop
     dt1= system_dt
     print '("identity difference tensor:")'
     call dv% say()
     err= dv% normb()
     print '("|v4**2 - v3**2 = (v4 - v3)*(v4 + v3)|_2 = "ES12.5)', err
     write (*,'("[+][TEST05-7]["F7.2" s]:&
        & identity is satisfied: ")', advance='no') dt1
     if (dabs(err).lt.tol) then
        write(*,'("PASS"/)')
     else
        write(*,'("FAIL"/)')
        failed= failed + 1
     endif

     if (allocated(a1)) deallocate(a1)
     if (allocated(a2)) deallocate(a2)

     print '("-- [5-8] test tensor rounding after multiplying by zero")'
     Z = dtt_tensor(n, r_=r)
     Z = Z*v5
     call Z% svd
     call Z% say()
     sa = size(Z% full())
     err= dnrm2(sa, Z% full(), 1)
     print '("|0*v5|_2 = ",ES12.5)', err
     if (err.lt.tol) then
        print '("[+][TEST04]: test rounding after multiplying by zero: PASS")'
     else
        print '("[-][TEST04]: test rounding after multiplying by zero: FAIL")'
        failed= failed + 1
     end if


     call MPI_Finalize(info)
     if(info.ne.0) error stop '[MPI]: Finalize fail'
     if(failed.eq.0) then
        print '("[+][TEST05]: PASSED")'
     else
        print '("[-][TEST05]: FAILED with ",I2," error(s)")', failed
        error stop failed
     end if
  end subroutine test_tensor_ops

end program
