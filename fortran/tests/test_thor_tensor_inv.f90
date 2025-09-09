!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_thor_tensor_inv.f90
!! @author [IDB], [OGK]
!! @date  March 2024
!! @brief Testing tensor element-wise inversion operation
!!

module shared_mod
use thor_lib, only: dtt_tensor
type(dtt_tensor) :: xg
end module shared_mod

program main
use rnd_lib, only: random, init_random_seed
use thor_lib, only: dtt_tensor
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 3  ! tensor order
integer :: n = 4  ! tensor dimension size
integer :: r = 2  ! tensor rank
integer :: l = 1  ! starting index
integer :: m = 3  ! ending index
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
     elseif (trim(arg).eq.'-l') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) l
        if (err.ne.0) error stop "ERROR: when parsing -l <???>"
        i= i + 1
     elseif (trim(arg).eq.'-m') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) m
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
        i= i + 1
     elseif (trim(arg).eq.'-r') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) r
        if (err.ne.0) error stop "ERROR: when parsing -r <???>"
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
  call test_tensor_inv(nn, rr, tol, l, m)
  deallocate(nn, rr, rnd)

contains

  subroutine help_msg()
  print '(12(/,A))', &
   "Test tensor elementwise inversion operation: v -> 1/v", &
   "Usage: ./test.exe [-h] [-options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : tensor order <ord> [3]", &
   " -n <dim>  : tensor dimensions <dim> [4]", &
   " -l <num>  : starting index [1]", &
   " -m <num>  : ending index [d]", &
   " -r <rank> : rank of the cores [2]", &
   " -t <tol>  : use this tolerance in SVD step [1e-12]", &
   " --randomize-nr: generate arrays of random ranks & dimensions <= n, r [F]", &
   " --rseed <num> : supply random seed [none]"
  end subroutine

  subroutine test_tensor_inv(n, r, tol, l, m)
  use time_lib
  use thor_lib, only : dtt_tensor, dtt_tensor_ones, operator(-), &
      dtt_tensor_ones_like, sayfull, dealloc, dtt_tensor_rand
  use matlab_struct_module, only: pprint_matrix
  implicit none
  include 'mpif.h'
  integer, intent(in) :: n(:), r(:)
  double precision, intent(in) :: tol
  type(dtt_tensor)    :: v1, v2, v3
  integer :: i, l, m, info, nElements, nproc,me, failed
  !double precision,external :: dfunc_inv_tt, dfunc_tt_gt !,dfunc_ising_discr
  double precision, allocatable      :: arr1(:)
  double precision                   :: eps, err

  print '("")'
  print '("Test 06: division 1/vec using DMRGG cross-interpolation")'
  print '("")'

  call mpi_init(info)
  call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
  call mpi_comm_rank(MPI_COMM_WORLD,me,info)
  call system_timer_init

  if(me.eq.0) print '("[MPI]: initialized OK: size = ",I5)', nproc
  call MPI_Barrier(  MPI_COMM_WORLD, info)

  eps = 1.0e-14; failed = 0

  print '("-- [6-1] testing: create small preset 4x4 matrix:")'
  arr1 = (/0.8641d0, 0.9178d0, 0.8869d0, 0.8460d0, &
           0.9757d0, 0.7043d0, 1.2966d0, 1.5632d0, &
           1.0847d0, 0.8228d0, 1.3630d0, 1.5189d0, &
           1.0d-2,   0.0001d0, 10.3d0,  -0.1d0/)
  print '("Input matrix arr1 =")'
  call pprint_matrix(reshape(arr1, [4,4]))
  print '("")'

  print '("-- [6-2] converting this array to tensor:")'
  v1= dtt_tensor(arr1,[2,2,2,2],eps)
  print*, "tt-tensor v1="
  call v1% say
  call pprint_matrix(reshape(v1% full(), [4,4]))
  print '("")'

  print*, "-- [6-3] computing elementwise inverse v2=1/v1:"
  v2 = inv_tens_d(v1, accuracy_=eps)
  !call svd(var_n, arr1, var_tt, eps)
  call v2% say
  call pprint_matrix(reshape(v2% full(), [4,4]), frmt_='(ES15.4)')
  print '("")'

  print*, "-- [6-4] create random tensor with specified ranks:"
  call system_timer_start
  v1 = dtt_tensor_rand(n, r)
  call system_timer_stop
  print '("random tensor v1 created ["F7.2" s]: ")',system_dt
  call v1% say
  print '("")'

  print*, "-- [6-5] computing pointwise inverse v2 = 1/v1:"
  call system_timer_start
  v2 = inv_tens_d(v1, accuracy_=eps)
  call system_timer_stop
  print '("tensor v2=1/v1 computed ["F7.2" s]: ")',system_dt
  call v1% say
  print '("")'

  print*, "-- [6-6] computing inverse again v3 = 1/v2:"
  call system_timer_start
  v3 = inv_tens_d(v2, accuracy_=eps)
  call system_timer_stop
  print '("tensor v3=1/v2 computed ["F7.2" s]: ")',system_dt
  call v1% say
  call system_timer_start
  v2 = v1 - v3
  err= v2% normb()
  call system_timer_stop
  print '("|v1 - 1/1/v1|_2 = ",ES10.3," ["F9.2" s]")', err, system_dt
  write (*,'("[+][TEST06-6]:&
     & double inverse v1 = 1/1/v1 identity (tol="ES8.1"): ")', &
     advance='no') tol
  if (dabs(err).lt.tol) then
     write(*,'("PASS",/)')
  else
     write(*,'("FAIL",/)')
     failed= failed + 1
  endif

  call MPI_Finalize(info)

  if(failed.eq.0) then
     print '("[+][TEST01]: PASSED")'
  else
     print '("[-][TEST01]: FAILED with ",I2," error(s)")', failed
     error stop failed
  end if

  end subroutine


  function inv_tens_d(x, accuracy_) result(y)
  use thor_lib, only: dtt_tensor_ones_like
  use dmrgg_lib, only: dtt_dmrgg
  use shared_mod, only: xg

  type(dtt_tensor) :: y
  class(dtt_tensor), intent(IN) :: x
  double precision, optional, intent(IN) :: accuracy_
  !
  double precision :: accuracy

     accuracy=1d-12; if (present(accuracy_)) accuracy= accuracy_
     xg = x  ! copy to module variable so inv_fun_tt can access it
     y = dtt_tensor_ones_like(x)
     call dtt_dmrgg(y, inv_fun_tt, accuracy=accuracy)

  end function inv_tens_d


  !pure 
  function inv_fun_tt(d, ind, nn) result(y)
  use thor_lib, only: tt_size, tijk
  use shared_mod, only: xg
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
