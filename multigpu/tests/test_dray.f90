!!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_dray.f90
!! @author [IDB], [OGK]
!! @date  December 2024
!! @brief Testing data structures in dray.f90
!!
program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: nx = 4  ! distributed 3D arrays: x-dimension
integer :: ny = 4  ! y-dimension
integer :: nz = 1  ! z-dimension
integer :: nb =-1  ! block size
integer :: nx1 = 3  ! distributed 3D arrays: x-dimension
integer :: ny1 = 6  ! y-dimension
integer :: nz1 = 4  ! z-dimension
integer :: nb1 =-1  ! block size
integer :: nx2 = 2  ! distributed 3D arrays: x-dimension
integer :: ny2 = 5  ! y-dimension
integer :: nz2 = 3  ! z-dimension
integer :: nb2 =-1  ! block size
integer :: M   = 12 ! M in C:MxN matrix
integer :: N   = 5  ! N in C:MxN matrix
integer :: K   = 6  ! contracted dimension
integer :: nzA = 4  ! A stack depth
integer :: nzB = 3  ! B stack depth
logical :: vpart = .true. ! test vertical partitioning [.true.]
logical :: trA = .false. ! is A transposed?
logical :: trB = .false. ! is B transposed?
integer :: what_test = 1 ! test type
integer :: rseed = 0 ! random seed (0 = from clock)
integer :: num_batches = 10

logical :: save_output = .false.

  ! parse command-line arguments
  i = 1
  do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     elseif (trim(arg).eq.'-nx') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) nx
        if (err.ne.0) error stop "ERROR: when parsing -nx <???>"
     elseif (trim(arg).eq.'-ny') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) ny
        if (err.ne.0) error stop "ERROR: when parsing -ny <???>"
     elseif (trim(arg).eq.'-nz') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) nz
        if (err.ne.0) error stop "ERROR: when parsing -nz <???>"
     elseif (trim(arg).eq.'-m') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) M
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
     elseif (trim(arg).eq.'-n') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) N
        if (err.ne.0) error stop "ERROR: when parsing -n <???>"
     elseif (trim(arg).eq.'-k') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) K
        if (err.ne.0) error stop "ERROR: when parsing -k <???>"
     elseif (trim(arg).eq.'-b') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) num_batches
        if (err.ne.0) error stop "ERROR: when parsing -k <???>"
     elseif (trim(arg).eq.'-nza') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) nzA
        if (err.ne.0) error stop "ERROR: when parsing -nza <???>"
     elseif (trim(arg).eq.'-nzb') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) nzB
        if (err.ne.0) error stop "ERROR: when parsing -nzb <???>"
     elseif (trim(arg).eq.'--rseed') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     elseif (trim(arg).eq.'-s' .or. trim(arg).eq.'--save') then
        save_output= .true.
     elseif (trim(arg).eq.'v' .or. trim(arg).eq.'V') then
        vpart = .true.
     elseif (trim(arg).eq.'h' .or. trim(arg).eq.'H') then
        vpart = .false.
     elseif (trim(arg).eq.'-nn') then
        trA = .false.; trB = .false.
     elseif (trim(arg).eq.'-nt') then
        trA = .false.; trB = .true.
     elseif (trim(arg).eq.'-tn') then
        trA = .true.; trB = .false.
     elseif (trim(arg).eq.'-tt') then
        trA = .true.; trB = .true.
     elseif (trim(arg).eq.'--struct') then
        what_test = 0
     elseif (trim(arg).eq.'--qr') then
        what_test = 1
     elseif (trim(arg).eq.'--mult') then
        what_test = 2
     elseif (trim(arg).eq.'--svd') then
        what_test = 3
     elseif (trim(arg).eq.'--sytrd') then
        what_test = 4
     else
        error stop "ERROR: unknown option"
     endif
     i = i + 1
  enddo

  select case(what_test)
  case(1)
     call test_structs(nx, ny, nz, vpart, save_output, rseed)
  case(1)
     call test_qr(nx, ny, nz, vpart, save_output, rseed)
  case(2)
     call test_stack_multiplication(M,N,K,nzA,nzB,-1,-1,vpart,trA,trB,rseed)
  case(3)
     call test_svd(M, N, vpart, num_batches, rseed)
  case(4)
     call test_sytrd(M, N, vpart, rseed)
  end select

contains

  subroutine help_msg()
  print '(28(/,A))', &
   "Test data structures in dray.f90: distributed arrays and stacks", &
   "Usage: ./test [-h] [--test-qr|--test-mult|--test-svd|--test-sytrd]", &
   "              [-nx <nx>] [-ny <ny>] [-nz <nz>] [-s|--save]", &
   "              [v|h] [-m <m>] [-n <n>] [-k <k>] [-nza <nza>] [-nzb <nzb>]", &
   "", &
   "Options:", &
   " -h : print this help message", &
   " --rseed <num> : supply random seed [0:use clock]", &
   " --test-structs : test distributed structures operations [default test]", &
   "   -nx <nx>   : global 3D array size in the x-direction", &
   "   -ny <ny>   : global 3D array size in the y-direction", &
   "   -nz <nz>   : blocking dimension in the z-direction", &
   "   -s|--save  : save output to files 'test_dray.OUT_??'", &
   "", &
   " --test-mult  : test stack multiplication", &
   "   -m  <m>    : 'm' in matrix A size (m x k)^OP", &
   "   -n  <n>    : 'n' in matrix B size (k x n)^OP", &
   "   -k  <k>    : 'k': contracted dimension in the multiplication test", &
   "   -nza <nza> : global z-dimension of matrix A as a stack", &
   "   -nzb <nzb> : global z-dimension of matrix B as a stack", &
   "   [v|h]      : in multiplication test: A and B are (h|v)-Stacks [v]", &
   " --test-svd   : test SVD decomposition", &
   " --test-sytrd : test SYTR decomposition", &
   "", &
   "Example:", &
   "  mrun -n 2 --host cn136,cn135 ./tests/test_dray.exe", &
   ""
  end subroutine


  subroutine test_structs(nxb, nyb, nzb, vpart, save_output, rseed)
  use rnd_lib, only: init_random_seed_mpi
  use matrix_ops
  use cudafor
  use distributed_arrays
  use cublas_aux, only: setup_cublas, cleanup_cublas
  use distributed_comms !, only  : super_comm
  use communicators, only: mpi_global_comm, mpi_sub_comm
  use time_lib, only: timef
  use matrix_util, only: pprint_matrix
  implicit none
  integer, intent(IN) :: nxb,nyb,nzb
  logical, intent(IN) :: vpart, save_output
  integer, intent(IN) :: rseed

  character(*), parameter      :: subnam = "[test_structs]"
  type(super_comm)             :: cpu_comm
  integer                      :: err, failed
  double precision             :: t1, t2, dt, nrm
  double precision, parameter  :: tol = 1d-12

  type(dPart) :: core, hp1, hp2,hp3,hp4, hp5
  type(stack) :: vCore, hCore
  type(stack) :: hst1, hst2, hst3, hst4, hst5, vst1, vst2, vst3, vst4, vst5
  type(stack) :: hQ, hR, vQ, vR, hC, vC

  integer :: i,j,k,nr,rank,nranks,idx,OFFSET,sh(2)
  integer :: i1, nx, ny, nz

  double precision, allocatable     :: h_A(:,:,:), h_A_T(:,:,:), h_A1(:,:,:), h_A2(:,:,:)
  double precision, allocatable     :: X_d_A1(:), X_d_A(:,:), X_d_A2(:,:), Y_d_A(:,:), Y_d_A2(:,:)
  double precision, allocatable, target, device  :: d_A(:,:,:), d_A_T(:,:,:), d_R2(:,:)
  double precision, pointer, device :: d_AV(:), d_AV_T(:)
  double precision, allocatable, target          :: h_QR(:,:,:), v_QR(:,:,:), QR(:,:)
  double precision, pointer                      :: A1d(:)

  ! local output file name
  character(len=8) lrankstr, numstr1, numstr2
  character(len=32) lrank_outfname
  integer :: ofid

     err  = 0
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
     write(lrankstr, '(I8)') rank
     lrankstr= adjustl(lrankstr)

     nx = nxb; ny = nyb; nz = nzb*nranks
     if (rank == 0) then
        print '("=========================")'
        print '("Running dray structures test on ",I5," ranks")', nranks
        print '("Block 3D: (nxb,nyb,nzb) = (",3(I5,1X),")")', nxb, nyb, nzb
        print '("Global 3D:(nx, ny, nz)  = (",3(I5,1X),")")', nx, ny, nz
        print '("hStack: (M, N) = (",I5," rows X ",I5," columns)")', nx, ny*nz
        print '("vStack: (M, N) = (",I5," rows X ",I5," columns)")', nx*nz, ny
        print '("---")'
     endif
     call init_random_seed_mpi(rseed)
     if (save_output) then
        if (rank.eq.0) print '(/,"--- local output files")'
        write(lrank_outfname, '("test_dray.OUT_",A)') trim(lrankstr)
        print '("Rank ",A,": opening """,A,""" for local rank output")', &
              trim(lrankstr), trim(lrank_outfname)
        ofid = 105 + rank
        open (ofid,file=trim(lrank_outfname),action="write")
     endif

     if (rank.eq.0) print '(/,"--- supercomm: init_cuda()")'
     t1=timef()
     err = cpu_comm% init_cuda(VRBZ_=3)
     dt = timef() - t1
     if (err==0) then
         print '(A,F9.2,A)',"    [+] CUDA backend initialized OK: dt= ", dt,'s'
     else
         error stop "ERROR: unable to initialize CUDA backend; exiting"
     endif
     t1=timef()
     call setup_cublas()

     allocate(h_A(nxb,nyb,nzb), h_A_T(nyb,nxb,nzb), h_A1(nxb,nyb,nzb), h_A2(nxb,nyb,nzb))
     allocate(d_A(nxb,nyb,nzb), d_A_T(nyb,nxb,nzb))

     if (rank.eq.0) print '(/,"--- initial data")'
     dt = timef() - t1
     OFFSET = nxb*nyb*nzb*rank
     idx = 1
     do k=1,nzb
        do j=1,nyb
           do i = 1,nxb
              h_A(i,j,k) = OFFSET + idx
              h_A1(i,j,k) = 10.d0*h_A(i,j,k)
              h_A2(i,j,k) = 100.d0*h_A(i,j,k)
              idx = idx + 1
           end do !i
        end do !j
     end do !k
     dt = timef() - t1
     call cpu_comm%barrier()
     print '(A,F9.2)', &
           "    [+]["//trim(lrankstr)//"] initial data finished: dt= ",dt

!!! >>> begin horizontal stacking test
     if (rank.eq.0) print '(/,"--- testing dPart constructor")'
     t1=timef()
     core = dPart(cpu_comm, h_A)
     call cpu_comm%barrier()
     dt = timef() - t1
     print '(A,F9.2)', &
           "    [+]["//trim(lrankstr)//"]core = dPart(..) OK, dt= ",dt

     if (rank == 0) print '(/,"--- testing hstack(dPart)")'
     t1=timef()
     hCore = stack(core, 'h')
     call cpu_comm%barrier()
     dt = timef() - t1
     X_d_A = hCore%d_A2
     sh = shape(X_d_A)
     if (rank == 0) then
        print '(A)', "    [+]["//trim(lrankstr)//"] hCore.d_A2 = "
        call pprint_matrix(X_d_A, frmt_='(F7.0)')
        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
     endif
     call cpu_comm%barrier()
     if (rank /= 0) then
        print '(A)', "    [+]["//trim(lrankstr)//"] hCore.d_A = "
        call pprint_matrix(X_d_A, frmt_='(F7.0)')
     endif
     call cpu_comm%barrier()

     t1=timef()
     hst1 = stack(cpu_comm, h_A1, 'h')
     call print_info(hst1)
     call cpu_comm%barrier()
     dt = timef() - t1
     X_d_A2 = hst1%d_A2
     sh = shape(X_d_A2)
     call cpu_comm%barrier()
     if (rank == 0) then
        print '(A)', "    [+]["//trim(lrankstr)//"] hst1.d_A2 = "
        call pprint_matrix(X_d_A2, frmt_='(F7.0)')
        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
     endif
     call cpu_comm%barrier()

     t1=timef()
     hst2 = hst1 + hCore
     call cpu_comm%barrier()
     dt = timef() - t1
     print '(A,F9.2)', "    [+]["//trim(lrankstr)//"] "// &
           "hst2 = hst1 + hCore computed OK| dt = ",dt
     X_d_A2 = hst2%d_A2
     sh = shape(X_d_A2)
     call cpu_comm%barrier()
     if (rank == 0) then
        print '(A)', "    [+]["//trim(lrankstr)//"] hst2.d_A2 = "
        call pprint_matrix(X_d_A2, frmt_='(F7.0)')
        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
     endif
     call cpu_comm%barrier()
     if (rank /= 0) then
        print '(A)', "    [+]["//trim(lrankstr)//"] hst2.d_A2 = "
        call pprint_matrix(X_d_A2, frmt_='(F7.0)')
     endif
     call cpu_comm%barrier()

     if (rank == 0) print '(/,"--- testing vstack(dPart)")'
     t1=timef()
     call hard_copy_dPart(vCore, core)
     call vstack_dPart(vCore, VRBZ_=0)
     call cpu_comm%barrier()
     dt = timef() - t1
     print '(A,F9.2)', "    [+]["//trim(lrankstr)//"] " // &
           " vCore = vstack(core) constructor OK, dt= ",dt
     call cpu_comm%barrier()
     Y_d_A = vCore%d_A2
     sh = shape(Y_d_A)
     call cpu_comm%barrier()
     if (rank == 0) then
        print '(A)', "    [+]["//trim(lrankstr)//"] vCore.d_A2 = "
        call pprint_matrix(Y_d_A, frmt_='(F7.0)')
        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
     endif
     call cpu_comm%barrier()
     if (rank /= 0) then
        print '(A)', "    [+]["//trim(lrankstr)//"] vCore.d_A2 = "
        call pprint_matrix(Y_d_A, frmt_='(F7.0)')
     endif
     call cpu_comm%barrier()


! !        if (hR% well_stacked()) then
! !           call hR%transpose_slices()
! !           if (save_output) then
! !              write(ofid,'("Matrix hR.T:")')
! !              QR = hR% d_A2
! !              sh = shape(QR)
! !              do i = 1,sh(1)
! !                 do j = 1,sh(2)
! !                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
! !                 enddo
! !                 write(ofid,*)
! !              enddo
! !           endif
! !        else
! !           print '("[w]'//subnam//': hR is not well-stacked, ' // &
! !                 'applying transpose_matrix instead")'
! !           call hR%transpose_matrix()
! !           if (save_output) then
! !              write(ofid,'("Matrix hR.T after transpose_matrix:")')
! !              QR = hR% d_A2
! !              sh = shape(QR)
! !              do i = 1,sh(1)
! !                 do j = 1,sh(2)
! !                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
! !                 enddo
! !                 write(ofid,*)
! !              enddo
! !              close(ofid)
! !           endif
! !        endif

! !        if (vR% well_stacked()) then
! !           call vR%transpose_slices()
! !           if (save_output) then
! !              write(ofid,'("Matrix R.T:")')
! !              QR = vR% d_A2
! !              sh = shape(QR)
! !              do i = 1,sh(1)
! !                 do j = 1,sh(2)
! !                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
! !                 enddo
! !                 write(ofid,*)
! !              enddo
! !              close(ofid)
! !           endif
! !        else
! !           print '("[w]'//subnam//': vR is not well-stacked, ' // &
! !                 'applying transpose_matrix instead")'
! !           call vR%transpose_matrix()
! !           if (save_output) then
! !              write(ofid,'("Matrix vR.T after transpose_matrix:")')
! !              QR = vR% d_A2
! !              sh = shape(QR)
! !              do i = 1,sh(1)
! !                 do j = 1,sh(2)
! !                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
! !                 enddo
! !                 write(ofid,*)
! !              enddo
! !              close(ofid)
! !           endif
! !        endif

     if (allocated(h_A)) deallocate(h_A)
     if (allocated(h_A_T)) deallocate(h_A_T)
     if (allocated(h_A1)) deallocate(h_A1)
     if (allocated(h_A2)) deallocate(h_A2)
     if (allocated(d_A)) deallocate(d_A)
     if (allocated(d_A_T)) deallocate(d_A_T)
     if (allocated(h_QR)) deallocate(h_QR)
     if (allocated(v_QR)) deallocate(v_QR)
     call cleanup_cublas()
     call MPI_FINALIZE(err)

  end subroutine test_structs


  subroutine test_qr(nxb, nyb, nzb, vpart, save_output, rseed)
  use rnd_lib, only: init_random_seed_mpi
  use matrix_ops
  use cudafor
  use distributed_arrays
  use cublas_aux, only: setup_cublas, cleanup_cublas
  use distributed_comms !, only  : super_comm
  use communicators, only: mpi_global_comm, mpi_sub_comm
  use time_lib, only: timef
  use matrix_util, only: pprint_matrix
  implicit none
  integer, intent(IN) :: nxb,nyb,nzb
  logical, intent(IN) :: vpart, save_output
  integer, intent(IN) :: rseed

  character(*), parameter      :: subnam = "[test_qr]"
  type(super_comm)             :: cpu_comm
  integer                      :: err, failed
  double precision             :: t1, t2, dt, nrm
  double precision, parameter  :: tol = 1d-12

  type(dPart) :: core, hp1, hp2,hp3,hp4, hp5
  type(stack) :: vCore, hCore
  type(stack) :: hst1, hst2, hst3, hst4, hst5, vst1, vst2, vst3, vst4, vst5
  type(stack) :: hQ, hR, vQ, vR, hC, vC

  integer :: i,j,k,nr,rank,nranks,idx,OFFSET,sh(2)
  integer :: i1, nx, ny, nz

  double precision, allocatable     :: h_A(:,:,:), h_A_T(:,:,:), h_A1(:,:,:), h_A2(:,:,:)
  double precision, allocatable     :: X_d_A1(:), X_d_A(:,:), X_d_A2(:,:), Y_d_A(:,:), Y_d_A2(:,:)
  double precision, allocatable, target, device  :: d_A(:,:,:), d_A_T(:,:,:), d_R2(:,:)
  double precision, pointer, device :: d_AV(:), d_AV_T(:)
  double precision, allocatable, target          :: h_QR(:,:,:), v_QR(:,:,:), QR(:,:)
  double precision, pointer                      :: A1d(:)

  ! local output file name
  character(len=8) lrankstr, numstr1, numstr2
  character(len=32) lrank_outfname
  integer :: ofid

     err  = 0
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
     write(lrankstr, '(I8)') rank
     lrankstr= adjustl(lrankstr)

     nx = nxb; ny = nyb; nz = nzb*nranks
     if (rank == 0) then
        print '("=========================")'
        print '("Running dray structures test on ",I5," ranks")', nranks
        print '("Block 3D: (nxb,nyb,nzb) = (",3(I5,1X),")")', nxb, nyb, nzb
        print '("Global 3D:(nx, ny, nz)  = (",3(I5,1X),")")', nx, ny, nz
        print '("hStack: (M, N) = (",I5," rows X ",I5," columns)")', nx, ny*nz
        print '("vStack: (M, N) = (",I5," rows X ",I5," columns)")', nx*nz, ny
        print '("---")'
     endif
     call init_random_seed_mpi(rseed)
     if (save_output) then
        if (rank.eq.0) print '(/,"--- local output files")'
        write(lrank_outfname, '("test_dray.OUT_",A)') trim(lrankstr)
        print '("Rank ",A,": opening """,A,""" for local rank output")', &
              trim(lrankstr), trim(lrank_outfname)
        ofid = 105 + rank
        open (ofid,file=trim(lrank_outfname),action="write")
     endif

     if (rank.eq.0) print '(/,"--- supercomm: init_cuda()")'
     t1=timef()
     err = cpu_comm% init_cuda(VRBZ_=3)
     dt = timef() - t1
     if (err==0) then
         print '(A,F9.2,A)',"    [+] CUDA backend initialized OK: dt= ", dt,'s'
     else
         error stop "ERROR: unable to initialize CUDA backend; exiting"
     endif
     t1=timef()
     call setup_cublas()

     allocate(h_A(nxb,nyb,nzb), h_A_T(nyb,nxb,nzb), h_A1(nxb,nyb,nzb), h_A2(nxb,nyb,nzb))
     allocate(d_A(nxb,nyb,nzb), d_A_T(nyb,nxb,nzb))

!     if (rank.eq.0) print '(/,"--- initial data")'
!     dt = timef() - t1
!     OFFSET = nxb*nyb*nzb*rank
!     idx = 1
!     do k=1,nzb
!        do j=1,nyb
!           do i = 1,nxb
!              h_A(i,j,k) = OFFSET + idx
!              h_A1(i,j,k) = 10.d0*h_A(i,j,k)
!              h_A2(i,j,k) = 100.d0*h_A(i,j,k)
!              idx = idx + 1
!           end do !i
!        end do !j
!     end do !k
!     dt = timef() - t1
!     call cpu_comm%barrier()
!     print '(A,F9.2)', &
!           "    [+]["//trim(lrankstr)//"] initial data finished: dt= ",dt
!
!!!! >>> begin horizontal stacking test
!     if (rank.eq.0) print '(/,"--- testing dPart constructor")'
!     t1=timef()
!     core = dPart(cpu_comm, h_A)
!     call cpu_comm%barrier()
!     dt = timef() - t1
!     print '(A,F9.2)', &
!           "    [+]["//trim(lrankstr)//"]core = dPart(..) OK, dt= ",dt
!
!     if (rank == 0) print '(/,"--- testing hstack(dPart)")'
!     t1=timef()
!     hCore = stack(core, 'h')
!     call cpu_comm%barrier()
!     dt = timef() - t1
!     X_d_A = hCore%d_A2
!     sh = shape(X_d_A)
!     if (rank == 0) then
!        print '(A)', "    [+]["//trim(lrankstr)//"] hCore.d_A2 = "
!        call pprint_matrix(X_d_A, frmt_='(F7.0)')
!        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
!     endif
!     call cpu_comm%barrier()
!     if (rank /= 0) then
!        print '(A)', "    [+]["//trim(lrankstr)//"] hCore.d_A = "
!        call pprint_matrix(X_d_A, frmt_='(F7.0)')
!     endif
!     call cpu_comm%barrier()
!
!     t1=timef()
!     hst1 = stack(cpu_comm, h_A1, 'h')
!     call print_info(hst1)
!     call cpu_comm%barrier()
!     dt = timef() - t1
!     X_d_A2 = hst1%d_A2
!     sh = shape(X_d_A2)
!     call cpu_comm%barrier()
!     if (rank == 0) then
!        print '(A)', "    [+]["//trim(lrankstr)//"] hst1.d_A2 = "
!        call pprint_matrix(X_d_A2, frmt_='(F7.0)')
!        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
!     endif
!     call cpu_comm%barrier()
!
!     t1=timef()
!     hst2 = hst1 + hCore
!     call cpu_comm%barrier()
!     dt = timef() - t1
!     print '(A,F9.2)', "    [+]["//trim(lrankstr)//"] "// &
!           "hst2 = hst1 + hCore computed OK| dt = ",dt
!     X_d_A2 = hst2%d_A2
!     sh = shape(X_d_A2)
!     call cpu_comm%barrier()
!     if (rank == 0) then
!        print '(A)', "    [+]["//trim(lrankstr)//"] hst2.d_A2 = "
!        call pprint_matrix(X_d_A2, frmt_='(F7.0)')
!        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
!     endif
!     call cpu_comm%barrier()
!     if (rank /= 0) then
!        print '(A)', "    [+]["//trim(lrankstr)//"] hst2.d_A2 = "
!        call pprint_matrix(X_d_A2, frmt_='(F7.0)')
!     endif
!     call cpu_comm%barrier()
!
!     if (rank == 0) print '(/,"--- testing vstack(dPart)")'
!     t1=timef()
!     call hard_copy_dPart(vCore, core)
!     call vstack_dPart(vCore, VRBZ_=0)
!     call cpu_comm%barrier()
!     dt = timef() - t1
!     print '(A,F9.2)', "    [+]["//trim(lrankstr)//"] " // &
!           " vCore = vstack(core) constructor OK, dt= ",dt
!     call cpu_comm%barrier()
!     Y_d_A = vCore%d_A2
!     sh = shape(Y_d_A)
!     call cpu_comm%barrier()
!     if (rank == 0) then
!        print '(A)', "    [+]["//trim(lrankstr)//"] vCore.d_A2 = "
!        call pprint_matrix(Y_d_A, frmt_='(F7.0)')
!        print '(A,F9.2,A)', "    [+]["//trim(lrankstr)//"] dt= ",dt," s"
!     endif
!     call cpu_comm%barrier()
!     if (rank /= 0) then
!        print '(A)', "    [+]["//trim(lrankstr)//"] vCore.d_A2 = "
!        call pprint_matrix(Y_d_A, frmt_='(F7.0)')
!     endif
!     call cpu_comm%barrier()

     if (.not.vpart) then
        if (rank == 0) print '(/,"--- testing hStack QR decomposition")'
        t1 = timef()
        h_QR = generate_dist_3dcore(nxb, nyb, nzb*nranks, nxb, nyb, nzb, 0, 0, rank)
        hst2 = stack(cpu_comm, h_QR, 'h')
        call cpu_comm%barrier()
        dt = timef() - t1
        print '("    [+][rank",A,"]: generate_matrix, stack OK| dt= ",F8.3," s")', &
              trim(lrankstr), dt

        if (save_output) then
           write (ofid,'("=== Testing QR factorization of hStack ===")')
           write (ofid,'("Matrix A (before converted to hst2):")')
           sh = shape(QR)
           do i = 1,nxb
              do k = 1,nzb
                 do j = 1,nyb
                    write(ofid, '(999(ES14.7,1X))', advance='no') h_QR(i,j,k)
                 enddo
                 if (k.lt.nzb) write(ofid,'(" | ")', advance='no')
              enddo
              write(ofid,*)
           enddo
        endif

        hst2% h_A2 = hst2% d_A2
        sh = shape(hst2% h_A2)
        if (save_output) then
           write (numstr1, '(I8)') sh(1)
           write (numstr2, '(I8)') sh(2)
           write (ofid,'("local matrix A(",A,",",A,") = ")') &
                  trim(adjustl(numstr1)), trim(adjustl(numstr2))
           do i = 1,sh(1)
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') hst2% h_A2(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

        t1 = timef()
        call hst2% qr(hQ, hR) !, VRBZ_=4)
        call cpu_comm%barrier()
        dt = timef() - t1
        print '("    [+][rank",A,"]: hst2% qr(hQ, hR) computed OK| dt= ",F8.3," s")', &
              trim(lrankstr), dt

        if (save_output) then
           write(ofid,'("Matrix Q:")')
           QR = hQ% d_A2
           sh = shape(QR)
           do i = 1,sh(1)
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo

           write(ofid,'("Matrix R:")')
           QR = hR% d_A2
           sh = shape(QR)
           do i = 1,sh(1)
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

        call cpu_comm%barrier()
        !hQ% N = hQ% M
        hC = AxB_stacks(hQ, hR, VRBZ_=0)
        if (save_output) then
           write(ofid,'("Matrix hC = AxB_stacks(hQ, hR):")')
           QR = hC% d_A2
           sh = shape(QR)
           do i = 1,sh(1) ! print transposed
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

        hst3 = hst2 - hC
        nrm = hst3% nrm2()
        if (rank == 0) print '("|hA - hQ@hR|_2 = ",ES14.7)', nrm
        if (save_output) then
           write(ofid,'("Matrix (hQ@hR - hA)% d_A:")')
           QR = hst3% d_A2
           sh = shape(QR)
           do i = 1,sh(1) ! print transposed
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

        if (hR% well_stacked()) then
           call hR%transpose_slices()
           if (save_output) then
              write(ofid,'("Matrix hR.T:")')
              QR = hR% d_A2
              sh = shape(QR)
              do i = 1,sh(1)
                 do j = 1,sh(2)
                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
                 enddo
                 write(ofid,*)
              enddo
           endif
        else
           print '("[w]'//subnam//': hR is not well-stacked, ' // &
                 'applying transpose_matrix instead")'
           call hR%transpose_matrix()
           if (save_output) then
              write(ofid,'("Matrix hR.T after transpose_matrix:")')
              QR = hR% d_A2
              sh = shape(QR)
              do i = 1,sh(1)
                 do j = 1,sh(2)
                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
                 enddo
                 write(ofid,*)
              enddo
              close(ofid)
           endif
        endif
!!! <<< end horizontal stacking test
     else !< v-Stack QR test

        if (rank == 0) print '(/,"--- testing vStack QR decomposition")'
        call cpu_comm% reset_cuda_handles()

        t1 = timef()
        !v_QR = generate_dist_3dcore(nxb, nyb, nzb*nranks, nxb, nyb, nzb, 0, 0, rank)
        !if (save_output) then
        !   write (ofid,*)
        !   write (ofid,'("=== Testing QR factorization of vStack ===")')
        !   write (ofid,'("Matrix A (before converted to vst2):")')
        !   sh = shape(QR)
        !   do k = 1,nzb
        !      do i = 1,nxb
        !         do j = 1,nyb
        !            !v_QR(i,j,k) = 100*i + 10*j + k ! use this to debug transposing
        !            write(ofid, '(999(ES14.5,1X))', advance='no') v_QR(i,j,k)
        !         enddo
        !         write(ofid,*)
        !      enddo
        !      write(ofid,*)
        !   enddo
        !endif
        !vst2 = stack(cpu_comm, v_QR, 'v')
M = nxb*nzb*nranks; N = nyb        
vst2 = shallow_stack(cpu_comm, M, N, 'v')
vst2% h_A2 = generate_dist_matrix(M, N, vst2%MB, vst2%NB, &
                                  rank, 0, nranks, 1)
vst2% d_A = vst2% h_A
        call cpu_comm%barrier()
        dt = timef() - t1
        print '("    [+][rank",A,"]: generate_matrix, stack OK| dt = ",F8.3," s")', &
              trim(lrankstr), dt

        if (save_output) then
           write(ofid,'("Matrix vst2% d_A:")')
           QR = vst2% d_A2
           sh = shape(QR)
           do i = 1,sh(1) ! print transposed
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

        t1 = timef()
!-!call print_info(vst2)
!-!call vst2% save_to_ascii("tmp/vst_"//str(N)//"x"//str(N))        
        call vst2% qr(vQ, vR)
!-!call vQ% save_to_ascii("tmp/vQ_"//str(N)//"x"//str(N))
!-!call vR% save_to_ascii("tmp/vR_"//str(N)//"x"//str(N))
        call cpu_comm%barrier()
        dt = timef() - t1
        print '("    [+][rank",A,"]: vst2% qr(vQ, vR) computed OK| dt= ",F8.3," s")', &
              trim(lrankstr), dt

        if (save_output) then
           write(ofid,'("Matrix Q:")')
           QR = vQ% d_A2
           sh = shape(QR)
           do i = 1,sh(1)
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo

           write(ofid,'("Matrix R:")')
           QR = vR% d_A2
           sh = shape(QR)
           do i = 1,sh(1)
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

!vQ = stack("tmp/vQ_300x300_122", cpu_comm)        
!vR = stack("tmp/vR_300x300_122", cpu_comm)        
!-!call vR% save_to_ascii("tmp/vR_"//str(N)//"x"//str(N))        
        vC = AxB_stacks(vQ, vR, VRBZ_=0)
        if (save_output) then
           write(ofid,'("Matrix vC = AxB_stacks(vQ, vR):")')
           QR = vC% d_A2
           sh = shape(QR)
           do i = 1,sh(1) ! print transposed
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

!        call vst2% switch_stacking(VRBZ_=0)
!        if (rank == 0) then
!           print '("Stack vst2:")'
!           call print_info(vst2)
!        endif
!        if (save_output) then
!           write(ofid,'("Matrix vst2% d_A after stacking was switched:")')
!           QR = vst2% d_A2
!           sh = shape(QR)
!           do i = 1,sh(1) ! print transposed
!              do j = 1,sh(2)
!                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
!              enddo
!              write(ofid,*)
!           enddo
!        endif
!
!        call vst2% switch_stacking(VRBZ_=0)
!        if (rank == 0) then
!           print '("Stack vst2:")'
!           call print_info(vst2)
!        endif
!        if (save_output) then
!           write(ofid,'("Matrix vst2% d_A after switching back:")')
!           QR = vst2% d_A2
!           sh = shape(QR)
!           do i = 1,sh(1) ! print transposed
!              do j = 1,sh(2)
!                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
!              enddo
!              write(ofid,*)
!           enddo
!        endif

        vst3 = vst2 - vC
        nrm = vst3% nrm2()
        if (rank == 0) print '("|vA - vQ@vR|_2 = ",ES14.7)', nrm
        if (save_output) then
           write(ofid,'("Matrix (vQ@vR - vA)% d_A:")')
           QR = vst3% d_A2
           sh = shape(QR)
           do i = 1,sh(1) ! print transposed
              do j = 1,sh(2)
                 write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
              enddo
              write(ofid,*)
           enddo
        endif

        if (vR% well_stacked()) then
           call vR%transpose_slices()
           if (save_output) then
              write(ofid,'("Matrix R.T:")')
              QR = vR% d_A2
              sh = shape(QR)
              do i = 1,sh(1)
                 do j = 1,sh(2)
                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
                 enddo
                 write(ofid,*)
              enddo
              close(ofid)
           endif
        else
           print '("[w]'//subnam//': vR is not well-stacked, ' // &
                 'applying transpose_matrix instead")'
           call vR%transpose_matrix()
           if (save_output) then
              write(ofid,'("Matrix vR.T after transpose_matrix:")')
              QR = vR% d_A2
              sh = shape(QR)
              do i = 1,sh(1)
                 do j = 1,sh(2)
                    write(ofid, '(999(ES14.7,1X))', advance='no') QR(i,j)
                 enddo
                 write(ofid,*)
              enddo
              close(ofid)
           endif
        endif
     endif ! if vpart

     if (allocated(h_A)) deallocate(h_A)
     if (allocated(h_A_T)) deallocate(h_A_T)
     if (allocated(h_A1)) deallocate(h_A1)
     if (allocated(h_A2)) deallocate(h_A2)
     if (allocated(d_A)) deallocate(d_A)
     if (allocated(d_A_T)) deallocate(d_A_T)
     if (allocated(h_QR)) deallocate(h_QR)
     if (allocated(v_QR)) deallocate(v_QR)
     call cleanup_cublas()
     call MPI_FINALIZE(err)

  end subroutine test_qr


  subroutine test_stack_multiplication(M, N, K, nzA, nzB, nbA, nbB, vpart, trA, trB, rseed)
  use rnd_lib, only: init_random_seed_mpi
  use matrix_ops
  use cudafor
  use distributed_arrays
  use cublas_aux, only: setup_cublas, cleanup_cublas
  use distributed_comms !, only  : super_comm
  use communicators, only: mpi_global_comm, mpi_sub_comm
  use time_lib, only: timef
  use matrix_util, only: pprint_matrix
  implicit none
  integer, intent(IN) :: M, N, K  !< matrix dimensions: for A@B, (MxK) @ (K@N)
  integer, intent(IN) :: nzA, nzB !< z-dimensions of matrices A and B
  integer, intent(IN) :: nbA, nbB !< 3D blocking dimension for partitioning
  logical, intent(IN) :: vpart    !< test vertical partitioning?
  logical, intent(IN) :: trA, trB !< is A or B comes in a transposed form?
  integer, intent(IN) :: rseed

  character(*),parameter :: subnam='[test_stack_multiplication]'
  type(super_comm)             :: cpu_comm
  integer                      :: err, failed
  double precision             :: t1, t2, dt, nrm
  double precision, parameter  :: tol = 1d-12

  type(dPart) :: dpA, dpB
  type(stack) :: A, B, C

  integer :: nx1, ny1, nz1 !< global dimensions of A
  integer :: nx2, ny2, nz2 !< global dimensions of B
  integer :: i,j,l,nr,rank,nranks,localRank,idx
  integer :: i1, j1, l1

  ! local output file name
  character(len=8) lrankstr, numstr1, numstr2
  character(len=32) lrank_outfname
  integer :: ofid

     print '("### Testing stack multiplication ###")'
     print '("Input parameters: M="I5", N="I5", K="I5", nzA="I5", nzB="I5)', &
                                M, N, K, nzA, nzB
     print '("                  nbA="I5", nbB="I5", v-Stacks:"L1", transposed:"2(L1))', &
                                nbA, nbB, vpart, trA, trB
     ! arguments consistency check
     if (vpart) then
        ! TODO test only works for z-part
        if (.not.trA .and. .not.trB) then
           if (mod(M, nzA) /= 0) &
              error stop "[!]"//subnam//": M must be divisible by nzA"
           if (mod(K, nzB) /= 0) &
              error stop "[!]"//subnam//": K must be divisible by nzB"
           nx1 = M/nzA; ny1 = K; nz1 = nzA
           nx2 = K/nzB; ny2 = N; nz2 = nzB
        else if (.not.trA .and. trB) then
           if (mod(M, nzA) /= 0) &
              error stop "[!]"//subnam//": M must be divisible by nzA"
           if (mod(N, nzB) /= 0) &
              error stop "[!]"//subnam//": N must be divisible by nzB"
           nx1 = M/nzA; ny1 = K; nz1 = nzA
           nx2 = N/nzB; ny2 = K; nz2 = nzB
        else if (trA .and. .not.trB) then
           if (mod(K, nzA) /= 0) &
              error stop "[!]"//subnam//": K must be divisible by nzA"
           if (mod(K, nzB) /= 0) &
              error stop "[!]"//subnam//": K must be divisible by nzB"
           nx1 = K/nzA; ny1 = M; nz1 = nzA
           nx2 = K/nzB; ny2 = N; nz2 = nzB
        else
           error stop "[!]"//subnam//": A^T @ B^T vpart combo is not supported"
        endif
     else
        if (.not.trA .and. .not.trB) then
           if (mod(K, nzA) /= 0) &
              error stop "[!]"//subnam//": K must be divisible by nzA"
           if (mod(N, nzB) /= 0) &
              error stop "[!]"//subnam//": N must be divisible by nzB"
           nx1 = M; ny1 = K/nzA; nz1 = nzA
           nx2 = K; ny2 = N/nzB; nz2 = nzB
        else if (.not.trA .and. trB) then
           if (mod(K, nzA) /= 0) &
              error stop "[!]"//subnam//": K must be divisible by nzA"
           if (mod(K, nzB) /= 0) &
              error stop "[!]"//subnam//": K must be divisible by nzB"
           nx1 = M; ny1 = K/nzA; nz1 = nzA
           nx2 = N; ny2 = K/nzB; nz2 = nzB
        else if (trA .and. .not.trB) then
           if (mod(M, nzA) /= 0) &
              error stop "[!]"//subnam//": M must be divisible by nzA"
           if (mod(N, nzB) /= 0) &
              error stop "[!]"//subnam//": N must be divisible by nzB"
           nx1 = K; ny1 = M/nzA; nz1 = nzA
           nx2 = K; ny2 = N/nzB; nz2 = nzB
        else
           error stop "[!]"//subnam//": A^T @ B^T hpart combo is not supported"
        endif
     endif

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
     call init_random_seed_mpi(rseed)

     rank= cpu_comm% rank
     nranks= cpu_comm% comm_size
     write(lrankstr, '(I8)') rank
     lrankstr= adjustl(lrankstr)

     if (save_output) then
        if (rank.eq.0) print '(/,"--- local output files")'
        write(lrank_outfname, '("test_dray.OUT_",A)') trim(lrankstr)
        print '("Rank ",A,": opening """,A,""" for local rank output")', &
              trim(lrankstr), trim(lrank_outfname)
        ofid = 105 + rank
        open (ofid,file=trim(lrank_outfname),action="write")
     endif

     if (rank.eq.0) print '(/,"--- supercomm: init_cuda()")'
     t1=timef()
     err = cpu_comm% init_cuda(VRBZ_=3)
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

     if (rank.eq.0) print '(/,"--- initial data")'
     if (rank.eq.0) print '(/,"    1. construct dpA:")'
     if (nbA > -1) then
        dpA = dPart(cpu_comm, [nx1,ny1,nz1], .true., part_axis_=3, mb_=nbA, VRBZ_=3)
     else
        dpA = dPart(cpu_comm, [nx1,ny1,nz1], .true., part_axis_=3, VRBZ_=3)
     endif
     do l = 1,dpA% nzl
        l1 = iloc2glob(l, dpA%nzb, rank, nranks)
        do j = 1,dpA% nyl
           j1 = iloc2glob(j, dpA%nyb, rank, nranks)
           do i = 1,dpA% nxl
              i1 = iloc2glob(i, dpA%nxb, rank, nranks)
              idx = i + dpA%nxl*((j-1) + dpA%nyl*(l-1))
              dpA% h_A(idx) = generate_3dcore(i1, j1, l1)
           end do !i
        end do !j
     end do !l
     dpA% d_A = dpA% h_A

     if (rank.eq.0) print '(/,"    2. construct dpB:")'
     if (nbB > -1) then
        dpB = dPart(cpu_comm, [nx2,ny2,nz2], .true., part_axis_=3, mb_=nbB, VRBZ_=3)
     else
        dpB = dPart(cpu_comm, [nx2,ny2,nz2], .true., part_axis_=3, VRBZ_=3)
     endif
     do l = 1,dpB% nzl
        l1 = iloc2glob(l, dpB%nzb, rank, nranks)
        do j = 1,dpB% nyl
           j1 = iloc2glob(j, dpB%nyb, rank, nranks)
           do i = 1,dpB% nxl
              i1 = iloc2glob(i, dpB%nxb, rank, nranks)
              idx = i + dpB%nxl*((j-1) + dpB%nyl*(l-1))
              dpB% h_A(idx) = generate_3dcore(i1, j1, l1)
           end do !i
        end do !j
     end do !l
     dpB% d_A = dpB% h_A
     call cpu_comm%barrier()
     print '(A,F9.4)', &
           "    [+]["//trim(lrankstr)//"] initial data finished: dt= ",dt
     call print_info(dpA, label_='--->>>> dpA:')
     call print_info(dpB, label_='--->>>> dpB:')

     if (rank.eq.0) print '(/,"    3. make Stacks A and B:")'
     if (vpart) then
        A = stack(dpA, 'v')
        B = stack(dpB, 'v')
     else
        A = stack(dpA, 'h')
        B = stack(dpB, 'h')
     endif
     call print_info(A, label_='--->>>> stack A:')
     call print_info(B, label_='--->>>> stack B:')

     if (.not.trA .and. .not.trB) then
        if (rank.eq.0) print '(/,"    3. multiply AxB_stacks:")'
        C = AxB_stacks(A, B, VRBZ_=3)
     elseif (.not.trA .and. trB) then
        if (rank.eq.0) print '(/,"    3. multiply AxBT_stacks:")'
        if (.not.vpart .and. (A% NB /= B% NB)) then
           print '("ERROR: A% NB =/= B% NB: "I5" =/= "I5)', A% NB, B% NB
           error stop
        endif
        C = AxBT_stacks(A, B, VRBZ_=3)
     elseif (trA .and. .not.trB) then
        if (rank.eq.0) print '(/,"    3. multiply ATxB_stacks:")'
        if (vpart .and. (A% MB /= B% MB)) then
           print '("ERROR: A% MB =/= B% MB: "I5" =/= "I5)', A% MB, B% MB
           error stop
        endif
        C = ATxB_stacks(A, B, VRBZ_=3)
     else
        ! TODO
     endif
     call print_info(C, label_='--->>>> stack C:')

     ! TODO: check if multiplication was done correctly
     call cleanup_cublas()
     call MPI_FINALIZE(err)
  end subroutine test_stack_multiplication


  subroutine test_svd(M, N, vpart, num_batches, rseed)
  use rnd_lib, only: init_random_seed_mpi
  use string_lib, only: str
  use distributed_arrays
  use distributed_comms !, only  : super_comm
  use rnd_lib, only: randn
  use time_lib, only: timef
  use matrix_util, only: pprint_matrix
  use cublas_aux, only: setup_cublas, cleanup_cublas
  implicit none
  integer, intent(IN) :: M, N
  logical, intent(IN) :: vpart
  integer, intent(IN) :: rseed
  integer, intent(IN) :: num_batches
  !
  character(*), parameter :: subnam = '[test_svd]'
  type(super_comm) :: cpu_comm
  type(stack) :: A, U, VT
  double precision, allocatable, device, target :: d_S(:)
  integer     :: err, rank, nranks, i
  character   :: st
  double precision, allocatable, target :: A3(:,:,:)

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
     call init_random_seed_mpi(rseed)

     rank= cpu_comm% rank
     nranks= cpu_comm% comm_size
     if (rank.eq.0) print '(/,"--- supercomm: init_cuda()")'
     t1=timef()
     err = cpu_comm% init_cuda(VRBZ_=3)
     dt = timef() - t1
     if (err==0) then
         print '(A,F9.2,A)',"    [+] CUDA backend initialized OK: dt= ", dt,'s'
     else
         error stop "ERROR: unable to initialize CUDA backend; exiting"
     endif
     t1=timef()
     call setup_cublas()

     st = 'h'; if (vpart) st = 'v'

     A = shallow_stack(cpu_comm, M, N, st)
!     A% h_A = randn(A%ML*A%NL)
     A% h_A2 = generate_dist_matrix(A%M,  A%N, &
                                    A%MB, A%NB, &
                                    A%dev_row_idx, A%dev_col_idx, &
                                    A%numRowDevices,A%numColDevices)
     A% d_A = A% h_A
     do nr = 1,nranks
        if (nr-1.eq.rank) then
           print '("['//str(rank)//'/'//str(nranks)//']")'
           call pprint_matrix(A% h_A2)
        endif
        call cpu_comm% barrier()
     enddo

     t1=timef()
     do i = 1, num_batches
        call A% svd(U, d_S, VT, VRBZ_=3)
        if (rank.eq.0) print '("running '//str(i)//'th time...")'
        if (num_batches > 1 .and. i == 1) t1=timef()
     enddo
     dt = timef() - t1
     if (num_batches > 1) dt = dt/dble(num_batches - 1)
     if (rank.eq.0) print '("stack_svd time: ", ES14.7)', dt

     call cleanup_cublas()
     call MPI_FINALIZE(err)
  end subroutine test_svd


  subroutine test_sytrd(M, N, vpart, rseed)
  use rnd_lib, only: init_random_seed_mpi
  use distributed_arrays
  use distributed_comms !, only  : super_comm
  use rnd_lib, only: randn
  use time_lib, only: timef
  use cublas_aux, only: setup_cublas, cleanup_cublas
  use matrix_util, only: pprint_matrix
  implicit none
  integer, intent(IN) :: M, N
  logical, intent(IN) :: vpart
  integer, intent(IN) :: rseed
  !
  character(*), parameter :: subnam = '[test_sytrd]'
  type(super_comm) :: cpu_comm
  type(stack) :: A, P, T, TP1, A1, A2
  double precision, allocatable, device, target :: d_S(:)
  integer     :: err, rank, nranks
  character   :: st
  double precision :: nrm

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
     call init_random_seed_mpi(rseed)

     rank= cpu_comm% rank
     nranks= cpu_comm% comm_size
     if (rank.eq.0) print '(/,"--- supercomm: init_cuda()")'
     t1=timef()
     err = cpu_comm% init_cuda(VRBZ_=3)
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

     st = 'h'; if (vpart) st = 'v'
     A1 = shallow_stack(cpu_comm, M, M, st)
     A1% h_A = randn(A%ML*A%NL)
     A1% d_A = A1% h_A

     ! symmetrize
     A2 = A1; call A2% transpose_matrix()
     call A1% redist_to_vstack()
     call A2% redist_to_vstack()
     A = A1 + A2

     ! call the SYTRD routine
     call A% sytrd(P, T)

     TP1 = AxB_stacks(T, P, 'nt')
     A1  = AxB_stacks(P, TP1) - A

     do nr=1, nranks
        if (rank + 1.eq.nr) then
           print '(/"['//str(nr)//'/'//str(nranks)//']=== P*T*P^ - A: =====")'
           A1% h_A2 = A1% d_A2
           call pprint_matrix(A1% h_A2)
        endif
        call cpu_comm% barrier()
     enddo

     nrm = A1% nrm2()
     if (rank == 0) print '(/"|P*T*P^ - A|_2 = ",ES12.3)', nrm

     call cleanup_cublas()
     call MPI_FINALIZE(err)
  end subroutine test_sytrd


  function generate_dist_matrix(M, N, MB, NB, iB, jB, nRows, nCols) result(A)
  use cublasmp_aux, only: numroc, iloc2glob
  implicit none
  integer, intent(in) :: M, N
  integer, intent(in) :: MB, NB
  integer, intent(in) :: iB, jB
  integer, intent(in) :: nRows, nCols
  double precision, allocatable :: A(:,:)
  !
  integer :: i, j, il, jl, kl, ML, NL, M_N_min

     ML = numroc(M, MB, iB, nRows)
     NL = numroc(N, NB, jB, nCols)

     allocate(A(ML,NL))
     A = 0d0

     do jl = 1, NL
        j = iloc2glob(jl, NB, jB, nCols)
        do il = 1, ML
           i = iloc2glob(il, MB, iB, nRows)
           !A(il, jl) = 10*i + j  !< test
           A(il, jl) = 1.13d0*(mod(i*i + j*j,13)) &
                     - 1d0/(2.471 + mod(i,6) + 0.17d0*j) &
                     - 0.2d0/(4.1 + mod(i,3) + 0.17d0*mod(j,7))
        enddo
     enddo

  end function generate_dist_matrix


  function generate_dist_3dcore(M, N, L, MA, NA, LA, iB, jB, kB) result(A)
  implicit none
  integer, intent(in) :: M, N, L
  integer, intent(in) :: MA, NA, LA
  integer, intent(in) :: iB, jB, kB
  double precision, allocatable :: A(:,:,:)
  !
  integer :: i, j, k, il, jl, kl, ML, NL, LL

     ML = min(MA, M - iB*MA)
     NL = min(NA, N - jB*NA)
     LL = min(LA, L - kB*LA)

     allocate(A(ML,NL,LL))
     A = 0d0

     do kl = 1, LL
        k = kB*LA + kl
        do jl = 1, NL
           j = jB*NA + jl
           do il = 1, ML
              i = iB*MA + il
              !A(il, jl, kl) = 100*i + 10*j + k !< test
              A(il, jl, kl) = 0.01d0*(i*i + j*j + k*k) &
                            - 1d0/(2.4 + mod(i,6) + 0.1d0*j - 1.3d0*k)
              !A(il, jl, kl) = mod(j, 6) - 4*mod(i, 3) &
              !              - 3d0/(dble(2*i/7) + dble(0.2d0 - j*k/6));
           enddo
        enddo
     enddo

  end function generate_dist_3dcore


  pure elemental double precision function generate_3dcore(i, j, k)
  implicit none
  integer, intent(in) :: i, j, k
  !
     generate_3dcore = dble(i) + 0.01d0*dble(j) + 0.0001d0*dble(k)
  end function generate_3dcore


end program main

