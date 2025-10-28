!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_cusolvermp_aux.f90
!! @author Ismael Boureima, Oleg Korobkin
!! @date  March 2025
!! @brief Testing cuSOLVERMp
!!

program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: m = 10         !< m in matrix A(m x k)
integer :: n = 10         !< n in matrix B(n x k)
integer :: k = 10         !< k: contracted dimension in both A and B
integer :: trA = 0        !< is A transposed? 0=no, 1=yes
integer :: trB = 0        !< is A transposed? 0=no, 1=yes
integer :: which_test = 1 !< which test to run?
                          !  1: test_qr    - test QR factorization
                          !  2: test_sytrd - test sym.matrix -> tridiag
logical :: vgrid = .true. !< vertical partitioning?
double precision, allocatable :: rnd(:)

  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     elseif (trim(arg).eq.'-m') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) m
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
     elseif (trim(arg).eq.'-n') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) n
        if (err.ne.0) error stop "ERROR: when parsing -n <???>"
     elseif (trim(arg).eq.'v'.or.trim(arg).eq.'V') then
        vgrid = .true.
     elseif (trim(arg).eq.'h'.or.trim(arg).eq.'H') then
        vgrid = .false.
     elseif (trim(arg).eq.'--test-qr') then
        which_test = 1
     elseif (trim(arg).eq.'--test-sytrd') then
        which_test = 2
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  select case (which_test)
  case(1)
     call test_qr(m, n, vgrid)
  case(2)
     call test_sytrd(m, n)
  case default
     error stop "ERROR: invalid test number"
  endselect

contains

  subroutine help_msg()
  print '(13(/,A))', &
   "Test various cusolvermp factorizations:", &
   " -w 1: QR factorization;", &
   " -w 2: SYTR, tri-diagonal factorization of a symmetric matrix;", &
   "Usage: ./test.exe [-h] [-m <M>] [-n <N>] [v|V|h|H] [-w <testnum>]", &
   "Options:", &
   " -h           : print this help message", &
   " -m <M>       : number of rows in the matrix [10]", &
   " -n <N>       : number of rows in B and columns in C [10]", &
   " [v|V]        : partition the matrix vertically, by rows [T]", &
   " [h|H]        : partition the matrix horizontally, by cols [F]", &
   " --test-qr    : run QR decomposition test [default]", &
   " --test-sytrd : run SYTRD decomposotion test"
  end subroutine


  !> test QR factorization
  subroutine test_qr(MM, NN, vgrid)
  use cusolvermp_aux
  use cuda_comm_aux, only: get_hw_topology, cuda_comm_init, cal_comm_destroy
  use mpi_f08
  use cudafor               !, ONLY: device
  use time_lib, only: timef
  use matrix_ops, only: create_eye_matrix, zero_below_diag
  integer, intent(IN) :: MM, NN  !< matrix dimension: MM x NN
  logical, intent(IN) :: vgrid   !< is MP grid horizontal? [T]
  !
  !--- compile-time parameters
  integer, parameter :: VRBZ = 2,     & !< high-verbosity level
                        VRBZ0= 0,     & !< low-verbosity level
                        ctx_color= 11   !< GPU context color

  integer :: rank, commSize, localRank, localCommSize, ctxRank, &
             ctxCommSize,  nDevLocal, color
  integer :: err, info, sh(2), i, j, i1, j1, nr
  double precision             :: t1, dt
  type(MPI_COMM), target   :: localComm, ctxComm
  type(c_ptr)        :: gridA, descrA, gridQ, descrQ, descrR
  integer :: M, N, MB, NB, numRowDevices, numColDevices, iRow, jCol
  double precision, allocatable :: full_A(:,:)
  double precision, allocatable, device, target :: d_A(:), d_tau(:)
  double precision, pointer, device :: d_Q(:)
  double precision, device, pointer :: d_A2(:,:), d_Q2(:,:)
  double precision, allocatable :: h_A(:),h_A2(:,:),h_tau(:),h_Q(:),h_Q2(:,:)
  integer, target :: ML, NL

     ! Get Hardware topological info
     call get_hw_topology(ctx_color, rank, commSize, localComm, localRank, &
                          localCommSize, ctxComm, ctxRank, ctxCommSize, &
                          nDevLocal, color, VRBZ_=VRBZ0)
     print '("[MBIN][MPI][get_hw_topology()]: rank="I4"| size="I4"'// &
           '| localRank="I4"| localCommSize="I4)', rank, commSize, &
           localRank, localCommSize
     
     ! Create CAL-COMM handle
     call cuda_comm_init(rank, commSize, localRank, VRBZ_=VRBZ0)

     ! Create CuSOLVERMP handle
     call cusolvermp_init(localRank, VRBZ_=VRBZ0)
     call show_cusolverMpHandle

     M = MM; N = NN

     if (vgrid) then ! vertical device grid
        numRowDevices = commSize;         numColDevices = 1
        MB = ceiling(real(M) / commSize); NB = N
        iRow = rank;                      jCol = 0
     else ! horizontal device grid
        numRowDevices = 1;                numColDevices = commSize
        MB = M;                           NB = ceiling(real(N) / commSize)
        iRow = 0;                         jCol = rank
     endif

     ! Get size of the local matrix (-> ML, NL)
     call set_dist_mat_row_col(M,N, MB,NB, rank, numRowDevices,numColDevices, &
                               ML,NL, VRBZ_=VRBZ0)
     call MPI_BARRIER(MPI_COMM_WORLD)     
     print '("[MBIN]["I4"/"I4"]: M,N,MB,NB,ML,NL=",10(I8,1X))', &
           rank, commSize, M,N,MB,NB,ML,NL

     ! initialize local matrix
     if (ML*NL > 0) then
        allocate(h_A2(ML,NL))
        do j = 1,NL
           j1 = j + jCol*NB
           do i = 1,ML
              i1 = i + iRow*MB
              h_A2(i,j)= Aij(i1,j1)
           enddo
        enddo
        h_A = reshape(h_A2, (/ML*NL/))
        d_A = h_A
     endif
     call print_local_matrix(h_A2, "matrix A:", rank, commSize)
     if (ML*NL > 0) deallocate(h_A, h_A2)

     ! Create CuSOLVERMP grids
     gridA= set_cusolverMpGrid(numRowDevices, numColDevices, VRBZ_=VRBZ)

     ! Allocate d_tau
     allocate(d_tau(N))

     ! Create CuSOLVERMP matrixDescriptor
     descrA= set_cusolverMpMatrixDesc(gridA, M, N, MB, NB, numRowDevices, iRow, VRBZ_=VRBZ)

     ! Distributed QR:
     call check_status_mpi("[MBIN]: dist_qr", &
        dist_geqrf(descrA, d_A, M, N, d_tau, VRBZ_=VRBZ))

     if (ML*NL > 0) then
        h_A = d_A
        h_A2 = reshape(h_A, (/ML, NL/))
     endif
     call print_local_matrix(h_A2, "matrix QR:", rank, commSize)
     if (ML*NL > 0) deallocate(h_A, h_A2)

     ! Create a distributed unit matrix
     call create_eye_matrix(d_Q, M, N, MB, NB, iRow, jCol, ML,NL, &
                            numRowDevices, numColDevices)
     call MPI_BARRIER(MPI_COMM_WORLD)     
     if (ML*NL > 0) then
        h_A = d_Q
        h_A2 = reshape(h_A, (/ML, NL/))
     endif
     call print_local_matrix(h_A2, "matrix I:", rank, commSize)
     if (ML*NL > 0) deallocate(h_A, h_A2)

     h_tau = d_tau
     write(*,'("d_tau = ",75(F6.2,1X))') h_tau(:)

     ! Create CuSOLVERMP grids
     gridQ= set_cusolverMpGrid(numRowDevices, numColDevices, VRBZ_=VRBZ)

     ! Create CuSOLVERMP matrixDescriptor for the identity matrix
     descrQ= set_cusolverMpMatrixDesc(gridQ, M, N, MB, NB, numRowDevices, iRow, VRBZ_=VRBZ0)

     call MPI_BARRIER(MPI_COMM_WORLD)     
     call check_status_mpi("[MBIN]: dist_omrqr", &
        dist_ormqr(descrA, d_A, 0, 0, M, N, descrQ, d_Q, N, d_tau, VRBZ_=4))

     if (ML*NL > 0) then
        h_Q = d_Q
        h_Q2 = reshape(h_Q, (/ML, NL/))
     endif
     call print_local_matrix(h_Q2, "matrix Q:", rank, commSize)
     if (ML*NL > 0) deallocate(h_Q, h_Q2)

     call destroy_cusolverMpMatrixDesc(descrQ, VRBZ_=VRBZ0)
     call destroy_cusolverMpGrid(gridQ, VRBZ_=VRBZ0)

     !call check_status_mpi("[MBIN]: zero_below_diag(QR) -> matrix R", &
     d_A2(1:ML,1:N) => d_A
     call zero_below_diag(d_A2,M,N,MB,NB,iRow,jCol,ML,NL,numRowDevices,numColDevices)

     h_A = d_A
     h_A2= reshape(h_A, (/ML, NL/))
     call print_local_matrix(h_A2, "matrix R:", rank, commSize)
     deallocate(h_A2, h_A)

     call destroy_cusolverMpMatrixDesc(descrA, VRBZ_=VRBZ0)
     call destroy_cusolverMpGrid(gridA, VRBZ_=VRBZ0)

     ! Destroy CuSOLVERMP handle
     call cusolvermp_destroy(VRBZ_=VRBZ0)

     ! Destroy CAL-COMM handle
     call cal_comm_destroy(VRBZ_=VRBZ0)
     call MPI_FINALIZE(err)

  end subroutine test_qr


  !> testing reduction of a matrix to a tridiagonal form
  subroutine test_sytrd(MM, NN)
  use matrix_util, only: pprint_matrix
  use string_lib, only: str
  use cusolvermp_aux
  use cuda_comm_aux, only: get_hw_topology, cuda_comm_init, cal_comm_destroy
  use mpi_f08
  use cudafor               !, ONLY: device
  use matrix_ops, only: create_eye_matrix
  integer, intent(IN) :: MM, NN
  !
  !--- compile-time parameters
  integer, parameter :: VRBZ = 2,     & !< high-verbosity level
                        VRBZ0= 0,     & !< low-verbosity level
                        ctx_color= 11   !<  GPU context color

  integer :: rank, commSize, localRank, localCommSize, ctxRank, &
             ctxCommSize,  nDevLocal, color
  integer            :: err, info, sh(2), i, j, nr
  type(MPI_COMM), target   :: localComm, ctxComm
  type(c_ptr)        :: gridA, descrA, gridQ, descrQ, descrR
  integer :: M, N, MA, NA, numRowDevices, numColDevices
  double precision, allocatable :: full_A(:,:)
  double precision, allocatable, device, target :: d_A(:), d_d(:), d_e(:), d_tau(:), d_Q(:)
  double precision, device, pointer :: d_A2(:,:)
  double precision, allocatable :: h_A(:),h_A2(:,:),h_tau(:),h_d(:),h_e(:)

     ! Get Hardware topological info
     call get_hw_topology(ctx_color, rank, commSize, localComm, localRank, &
                          localCommSize, ctxComm, ctxRank, ctxCommSize, &
                          nDevLocal, color, VRBZ_=VRBZ0)
     print '("[MAIN][MPI][get_hw_topology()]: rank="I4"| size="I4"'// &
           '| localRank="I4"| localCommSize="I4)', rank, commSize, &
           localRank, localCommSize

     ! Create CAL-COMM handle
     call cuda_comm_init(rank, commSize, localRank, VRBZ_=VRBZ0)

     ! Create CuSOLVERMP handle
     call cusolvermp_init(localRank, VRBZ_=VRBZ0)
     call show_cusolverMpHandle

     M = MM
     MA = MM/commSize
     if (MA*commSize /= MM) error stop "ERROR: M must be divisible by # of ranks"

     allocate(full_A(M, M))
     call random(full_A)
     do i = 1,M
        full_A(i,i+1:M)= 0d0 ! zero upper triangle
        !full_A(i,1:i-1)= 0d0 ! zero lower triangle
     enddo
     h_A2 = full_A(1+rank*MA:(rank+1)*MA,1:M)
     h_A = reshape(h_A2, (/MA*M/))
     d_A = h_A
     call MPI_BARRIER(MPI_COMM_WORLD)

     do nr = 1,commSize
        if (rank.eq.nr-1) then
           print '("[MAIN][MPI][",I1,"/",I1,"]: local matrix A = ")', rank, commSize
           sh = shape(h_A2)
           do i = 1,sh(1)
              print '(75(F8.4,1X))', h_A2(i,:)
           enddo
        endif
        call MPI_BARRIER(MPI_COMM_WORLD)
     enddo
     deallocate(h_A, h_A2)

     call MPI_BARRIER(MPI_COMM_WORLD)

     ! Create CuSOLVERMP grids
     numRowDevices= commSize
     numColDevices= 1
     gridA= set_cusolverMpGrid(numRowDevices, numColDevices, VRBZ_=VRBZ)

     ! Allocate d_tau
     allocate(d_tau(M),d_d(M),d_e(M))

     ! Create CuSOLVERMP matrixDescriptor
     descrA= set_cusolverMpMatrixDesc(gridA, M, M, MA, MA, numRowDevices, rank, VRBZ_=VRBZ0)

     ! Distributed A -> P @ T @ P' factorization:
     call MPI_BARRIER(MPI_COMM_WORLD)
     call check_status_mpi("[MAIN]: dist_sytrd", &
        dist_sytrd(M, descrA, d_A, d_d, d_e, d_tau))

     h_A= d_A
     h_A2 = reshape(h_A, (/MA, M/))
     h_d= d_d
     h_e= d_e
     do nr = 1,commSize
        if (rank.eq.nr-1) then
           print '("[MAIN][MPI][",I1,"/",I1,"]: matrix T = ")', rank, commSize
           sh = shape(h_A2)
           do i = 1,sh(1)
              print '(75(F8.4,1X))', h_A2(i,:)
           enddo
           print '("[MAIN][MPI][",I1,"/",I1,"]: main diagonal = ")', rank, commSize
           print '(75(F8.4,1X))', h_d(1:M)
           print '("[MAIN][MPI][",I1,"/",I1,"]: subdiagonal = ")', rank, commSize
           print '(75(F8.4,1X))', h_e(1:M)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD)
     enddo
     deallocate(h_A, h_A2, h_d, h_e)

     gridQ= set_cusolverMpGrid(numRowDevices, numColDevices, VRBZ_=VRBZ)

     ! Create CuSOLVERMP matrixDescriptor
     descrQ= set_cusolverMpMatrixDesc(gridQ, M, M, MA, MA, numRowDevices, rank, VRBZ_=VRBZ0)

     ! Allocate d_Q
     allocate(d_Q(M*MA))
     d_Q = 0d0

     call MPI_BARRIER(MPI_COMM_WORLD)
     ! Tridiagonal matrix eigenvalue solve
     call check_status_mpi("[MAIN]: dist_stedc", &
        dist_stedc(M, d_d, d_e, descrQ, d_Q, VRBZ_=3))

     h_d= d_d
     h_A= d_Q
     h_A2 = reshape(h_A, (/MA, M/))
     do nr = 1,commSize
        if (rank.eq.nr-1) then
           print '("[MAIN][MPI][",I1,"/",I1,"]: eigenvalues = ")', rank, commSize
           print '(75(F8.4,1X))', h_d(1:M)
           print '("[MAIN][MPI][",I1,"/",I1,"]: matrix Q = ")', rank, commSize
           sh = shape(h_A2)
           do i = 1,sh(1)
              print '(75(F8.4,1X))', h_A2(i,:)
           enddo
        endif
        call MPI_BARRIER(MPI_COMM_WORLD)
     enddo
     deallocate(h_A, h_A2, h_d)
     ! Destroy CuSOLVERMP handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cusolvermp_destroy(VRBZ_=VRBZ0)

     ! Destroy CAL-COMM handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cal_comm_destroy(VRBZ_=VRBZ0)

     call MPI_FINALIZE(err)

  end subroutine test_sytrd


  subroutine check_status_mpi(msg, err)
  use mpi_f08
  implicit none
  integer, intent(in)      :: err
  character(*), intent(in) :: msg
  !
  integer :: rank, commsize

     call MPI_Comm_rank(MPI_COMM_WORLD, rank)
     call MPI_Comm_size(MPI_COMM_WORLD, commsize)
     if (err.eq.0) then
        print '("[+][",I1,"/",I1,"]",A," : OK")', rank, commsize, msg
     else
        print '("[!][",I1,"/",I1,"]",A," : FAIL")',rank, commsize, msg

     endif

  end subroutine check_status_mpi

!  typedef enum{
!    CUSOLVER_STATUS_SUCCESS=0,
!    CUSOLVER_STATUS_NOT_INITIALIZED=1,
!    CUSOLVER_STATUS_ALLOC_FAILED=2,
!    CUSOLVER_STATUS_INVALID_VALUE=3,
!    CUSOLVER_STATUS_ARCH_MISMATCH=4,
!    CUSOLVER_STATUS_MAPPING_ERROR=5,
!    CUSOLVER_STATUS_EXECUTION_FAILED=6,
!    CUSOLVER_STATUS_INTERNAL_ERROR=7,
!    CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED=8,
!    CUSOLVER_STATUS_NOT_SUPPORTED = 9,
!    CUSOLVER_STATUS_ZERO_PIVOT=10,
!    CUSOLVER_STATUS_INVALID_LICENSE=11,
!    CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED=12,
!    CUSOLVER_STATUS_IRS_PARAMS_INVALID=13,
!    CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC=14,
!    CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE=15,
!    CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER=16,
!    CUSOLVER_STATUS_IRS_INTERNAL_ERROR=20,
!    CUSOLVER_STATUS_IRS_NOT_SUPPORTED=21,
!    CUSOLVER_STATUS_IRS_OUT_OF_RANGE=22,
!    CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES=23,
!    CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED=25,
!    CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED=26,
!    CUSOLVER_STATUS_IRS_MATRIX_SINGULAR=30,
!    CUSOLVER_STATUS_INVALID_WORKSPACE=31
!} cusolverStatus_t;

  double precision function Aij(i, j)
  integer, intent(IN) :: i, j
  Aij = dble(mod(j,7) + mod(i,3))/dble(1 + mod(i,5) + mod(j,11)) - 1d0
  !Aij = dble(10*i + j)
  end function


  subroutine print_local_matrix(A2, label, rank, nranks)
  use mpi_f08
  character(*), intent(IN)     :: label
  double precision, intent(IN) :: A2(:,:)
  integer, intent(IN)          :: rank, nranks
  integer :: i, nr, sh(2)
     sh = shape(A2)
     if (sh(1).gt.20.or.sh(2).gt.20) return
     do nr = 1,nranks
        if (rank.eq.nr-1) then
           print '(A, "---------------")', label
           do i = 1,sh(1)
              print '(A,"["I2"/"I2"]: ",100(ES10.3,1X))', &
                     label, rank, nranks, A2(i,:)
           enddo
        endif
        call MPI_BARRIER(MPI_COMM_WORLD)
     enddo
  end subroutine


end program main
