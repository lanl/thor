!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_dthor_constructor.f90
!! @author [IDB], [OGK]
!! @date  August 2025
!! @brief Testing THOR construction in thor.f90
!!
program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 4     ! tensor order
integer :: n = 3     ! tensor dimension size
integer :: r = 2     ! tensor rank
integer :: rseed = 0 ! random seed (0 = from clock)
logical :: save_output = .false.
integer, allocatable :: nn(:), rr(:)

  ! parse command-line arguments
  i = 1
  do while (i <= command_argument_count())
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
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     elseif (trim(arg).eq.'-s' .or. trim(arg).eq.'--save') then
        save_output= .true.
     else
        error stop "ERROR: unknown option"
     endif
     i = i + 1
  enddo

  allocate(nn(d), rr(0:d))
  nn = n
  rr = r; rr(0) = 1; rr(d) = 1

  call test_dthor_rounding(nn, rr, rseed, save_output)
  !call debug_rounding(nn, rr, rseed, save_output)

contains

  subroutine help_msg() !! TODO: add options
  print '(13(/,A))', &
   "Test data structures in dist_thor.f90: distributed THOR", &
   "Usage: ./test.exe [-h] [-s]", &
   "", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : matrix order <ord> [4]", &
   " -n <dim>  : max dimension for each mode [3]", &
   " -r <rank> : rank of the cores [2]", &
   " --rseed <num> : supply random seed [0:use clock]", &
   " -s|--save : save output to files 'test_dray.OUT_??'", &
   "", &
   "Example:", &
   "  mrun -n 2 --host cn136,cn135 ./tests/test_dthor_constructor.exe", &
   ""
  end subroutine


  subroutine test_dthor_rounding(nn, rr, rseed, save_output)
  use thor
  use dthorio_lib, only: dthor_write_ascii_files, dthor_read_ascii_files
  use cublas_aux, only: setup_cublas, cleanup_cublas, cublas_nrm2
  use distributed_comms !, only  : super_comm
  use time_lib, only: timef
  use matrix_util, only: pprint_matrix, pprint_matrix3d
  use string_lib, only: str
  use rnd_lib, only: init_random_seed_mpi
  implicit none
  integer, intent(IN) :: nn(:)
  integer, intent(IN) :: rr(0:size(nn))
  integer, intent(IN) :: rseed
  logical, intent(IN) :: save_output
  !
  type(super_comm)             :: cpu_comm
  integer                      :: err, failed
  double precision             :: t1, t2, dt, nrm
  double precision, parameter  :: tol = 1d-12
  type(dthor) :: S, A, A0, AA
  integer :: i,j,k,rank,nranks

     err  = -1
     cpu_comm = super_comm()
     err = cpu_comm % init_mpi()
     if (err==0) then
        if (cpu_comm% rank == 0) print '(/,"-- super_comm constructor: OK")'
     else
        print '("FAILED TO INITIALIZE: err=", I7)', err
        call MPI_FINALIZE(err)
        error stop
     endif

     rank= cpu_comm% rank
     nranks= cpu_comm% comm_size
     d = size(nn)

     if (rank == 0) then
        print '("=========================")'
        print '("Running dthor_constructor test on ",I5," ranks")', nranks
        print '("Tensor modes and ranks: d = '//str(d)// &
                ', n = '//str(nn(1))//', r = '//str(rr(1))//'")'
        print '("---")'
     endif
     call init_random_seed_mpi(rseed)

     if (rank.eq.0) print '(/,"--- supercomm: init_cuda()")'
     t1=timef()
     err = cpu_comm% init_cuda()
     dt = timef() - t1
     if (err==0) then
         print '(A,F9.2,A)',"    [+] CUDA backend initialized OK: dt= ", dt,'s'
     else
         error stop "ERROR: unable to initialize CUDA backend; exiting"
     endif
     t1=timef()
     err = cpu_comm% init_nccl_comm()
     dt = timef() - t1
     if (err==0) then
         print '(A,F9.2,A)',"    [+] NCCL initialized OK: dt= ", dt,'s'
     else
         error stop "ERROR: unable to initialize NCCL; exiting"
     endif
     call setup_cublas()

     if (rank.eq.0) print '(/,"--- generating random dthor constructor")'
     A = random_dthor(cpu_comm, nn, stacking_='v', r_=rr) !, VRBZ_=3)

     if (rank.eq.0) print '(/,"--- computing tensor A+A+A+A and subtracting 4*A")'
     AA = A + A + A + A
     S = AA - 4*A

     if (rank.eq.0) print '(/,"--- running round_rand_orth")'
     A0 = S% round_rand_orth([(1,i=0,d)]) !, VRBZ_=3)
     if (rank.eq.0) print '(/,"--- A0 = round(AAAA - 4*A):")'
     print '("['//str(rank)//'/'//str(nranks)//'] ==> A0% cores: ")'
     call A0% pprint_cores(frmt_='(ES14.7)')
     print '("['//str(rank)//'/'//str(nranks) // &
             '] ==> norm of the last core of A0 = ",ES9.3)',&
            cublas_nrm2(A0% core(A0%d)% d_A)

     call cleanup_cublas()
     call MPI_FINALIZE(err)

  end subroutine test_dthor_rounding


  subroutine debug_rounding(n_in, r_in, rseed, save_output)
  use matrix_ops
  use cudafor
  use thor
  use distributed_arrays
  use dthorio_lib, only: dthor_write_ascii_files, dthor_read_ascii_files
  use cublas_aux, only: setup_cublas, cleanup_cublas, cublas_nrm2
  use distributed_comms !, only  : super_comm
  use time_lib, only: timef
  use matrix_util, only: pprint_matrix, pprint_matrix3d
  use string_lib, only: str
  use rnd_lib, only: init_random_seed_mpi
  implicit none
  integer, intent(IN) :: n_in(:)
  integer, intent(IN) :: r_in(0:size(nn))
  integer, intent(IN) :: rseed
  logical, intent(IN) :: save_output
  !
  type(super_comm)             :: cpu_comm
  integer                      :: err, failed
  double precision             :: t1, t2, dt, nrm
  double precision, parameter  :: tol = 1d-12

  type(dthor) :: S, A, A0, A2, A2r, Sr, B
  type(dthor) :: AA, BB, AB

  integer :: i,j,k,rank,nranks
  integer, allocatable :: nn(:), rr(:)

     err  = -1
     cpu_comm = super_comm()
     err = cpu_comm % init_mpi()
     if (err==0) then
        if (cpu_comm% rank == 0) print '(/,"-- super_comm constructor: OK")'
     else
        print '("FAILED TO INITIALIZE: err=", I7)', err
        call MPI_FINALIZE(err)
        error stop
     endif

     rank= cpu_comm% rank
     nranks= cpu_comm% comm_size
     call init_random_seed_mpi(rseed)

     if (rank.eq.0) print '(/,"--- supercomm: init_cuda()")'
     t1=timef()
     err = cpu_comm% init_cuda()
     dt = timef() - t1
     if (err==0) then
         print '(A,F9.2,A)',"    [+] CUDA backend initialized OK: dt= ", dt,'s'
     else
         error stop "ERROR: unable to initialize CUDA backend; exiting"
     endif
     t1=timef()
     err = cpu_comm% init_nccl_comm()
     dt = timef() - t1
     if (err==0) then
         print '(A,F9.2,A)',"    [+] NCCL initialized OK: dt= ", dt,'s'
     else
         error stop "ERROR: unable to initialize NCCL; exiting"
     endif
     call setup_cublas()

     if (save_output) then
        nn = n_in
        rr = r_in
        d  = size(n_in)
        if (rank.eq.0) then
           print '(/,"--- generating random thor object")'
           print   '("    nn = ",100(I2,1X))', nn(1:d)
           print   '("    rr = ",100(I2,1X))', rr(0:d)
        endif
        A = random_dthor(cpu_comm, nn, stacking_='h', r_=rr) !, VRBZ_=3)
        call dthor_write_ascii_files(A, "A")
     else
        if (rank.eq.0) print '(/,"--- reading thor object from A_?.dat")'
        A = dthor_read_ascii_files("A", cpu_comm, VRBZ_=3)
        nn = A% n
        rr = A% r
        d  = A% d
        if (rank.eq.0) then
           print '(/,"--- reading thor object from A*.dat files")'
           print   '("    nn = ",100(I2,1X))', nn(1:d)
           print   '("    rr = ",100(I2,1X))', rr(0:d)
        endif
     endif
     call A% to_vstack()
     !call A% pprint_cores()

     if (rank.eq.0) print '(/,"--- testing A -> rounded with rank-1 B")'
     B= empty_dthor(cpu_comm, nn, stacking_='v')
     do k = 1,d
        B% core(k)% d_A = 1.5d0
     enddo

     if (rank.eq.0) print '(/,"--- testing with A2 = A + A:")'
     A2 = A + A
     call dthor_write_ascii_files(A2, "A2", VRBZ_=0)
     call A% to_vstack()
     S = A2 - 2*A
     call S% to_vstack()
     call dthor_write_ascii_files(S, "S", VRBZ_=1)

     if (rank.eq.0) print '(/,"--- computing round(A2 - 2*A):")'
     Sr= S% round_rand_orth([(1,i=0,d)], RTT_=B)
     call dthor_write_ascii_files(Sr, "Sr", VRBZ_=0)
     print '("==> norm of the last core of Sr = ",ES9.3)',&
             cublas_nrm2(Sr% core(S%d)% d_A)
!     call Sr% pprint_cores()

     call cleanup_cublas()
     call MPI_FINALIZE(err)

  end subroutine debug_rounding


end program main

