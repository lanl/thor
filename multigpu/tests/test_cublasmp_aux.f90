!!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_cublasmp_aux.f90
!! @author Ismael Boureima, Oleg Korobkin
!! @date  March 2025
!! @brief Testing cuBLASMp
!!
program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: m = 10  ! m in matrix A(m x k)
integer :: n = 10  ! n in matrix B(n x k)
integer :: k = 10  ! k: contracted dimension in both A and B
integer :: m_b =-1 ! blocking dimension in M; (M - 1)/p + 1 by default
integer :: n_b =-1 ! blocking dimension in N; (N - 1)/p + 1 by default
integer :: k_b =-1 ! blocking dimension in K; (K - 1)/p + 1 by default
integer :: trA = 0 ! is A transposed? 0=no, 1=yes
integer :: trB = 0 ! is A transposed? 0=no, 1=yes
logical :: vpart = .true. ! test vertical partitioning [.true.]
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
     elseif (trim(arg).eq.'-k') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) k
        if (err.ne.0) error stop "ERROR: when parsing -k <???>"
     elseif (trim(arg).eq.'-mb') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) m_b
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
     elseif (trim(arg).eq.'-mb') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) m_b
        if (err.ne.0) error stop "ERROR: when parsing -mb <???>"
     elseif (trim(arg).eq.'-nb') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) n_b
        if (err.ne.0) error stop "ERROR: when parsing -nb <???>"
     elseif (trim(arg).eq.'-kb') then
        i=i+1; call getarg(i, arg)
        read(arg,*,iostat=err) k_b
        if (err.ne.0) error stop "ERROR: when parsing -kb <???>"
     elseif (trim(arg).eq.'-nn') then
        trA = 0; trB = 0
     elseif (trim(arg).eq.'-nt') then
        trA = 0; trB = 1
     elseif (trim(arg).eq.'-tn') then
        trA = 1; trB = 0
     elseif (trim(arg).eq.'-tt') then
        trA = 1; trB = 1
     elseif (trim(arg).eq.'v' .or. trim(arg).eq.'V') then
        vpart = .true.
     elseif (trim(arg).eq.'h' .or. trim(arg).eq.'H') then
        vpart = .false.
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  !call test_numroc(m)
  !call test_gemm(m, n, k, m_b, n_b, k_b, trA, trB, vpart)
  !call test_gemr2d(m, n)
  !call test_truncate(m, n)
  !call test_gemm_mn(m, n)
  call test_redistrib(m, n, m_b, n_b)

contains

  subroutine help_msg()
  print '(13(/,A))', &
   "Test matrix-matrix multiplication: C = A.Op @ B.Op", &
   "Usage: ./test.exe [-h] [-m <M>] [-n <N>] [-k <K>] ", &
   "Usage:            [-mb <Mb>] [-nb <Nb>] [-kb <Kb>] [-nn|nt|tn|tt]", &
   "Options:", &
   " -h           : print this help message", &
   " -m <M>       : number of rows in A and C [10]", &
   " -n <N>       : number of columns in B and C [10]", &
   " -k <K>       : number of columns in A and rows in B [10]", &
   " -mb <Mb>     : block dimension for M [int((M - 1)/nranks) + 1]", &
   " -nb <Nb>     : block dimension for N [int((N - 1)/nranks) + 1]", &
   " -kb <Kb>     : block dimension for K [int((K - 1)/nranks) + 1]", &
   " -nn|nt|tn|tt : if A and B are transposed or not [nn]"
  end subroutine


  !> Output computed values of the numroc function and iloc2glob
  subroutine test_numroc(M)
  use cublasmp_aux
  implicit none
  integer, intent(IN):: M
  !
  integer, parameter :: nranks_max = 10
  integer, parameter :: RC_SRCA = 0 !< cublasmp only works for this = 0
  integer :: nranks, MB, i, r
  integer :: nr(nranks_max), cnr(nranks_max), d(nranks_max)
  integer, allocatable :: ind(:), iglo(:)
     print '("# M = "I8", RC_SRCA = "I8)', M, RC_SRCA
     allocate(ind(max(nranks_max, M)),iglo(max(nranks_max, M)))
     do r = 1, max(nranks_max, M)
        ind(r) = r
     enddo
     do nranks = 1, nranks_max
        do MB = 1, 2*M
           nr(1:nranks)= numroc(M, MB, ind(1:nranks) - 1, nranks)
           !cnr(1:nranks)= cublasMpNumroc(M, MB, ind(1:nranks) - 1, 0, nranks)
           !d(1:nranks)= cnr(1:nranks) - nr(1:nranks)
           print '("R= "I3" MB= "I3" NUMROC="24(I3))', nranks, MB, nr(1:nranks)
           do r = 1, nranks
              print '("  - r = "I3" : ",24(I3))', r - 1, &
                     iloc2glob(ind(1:nr(r)), MB, r-1, nranks)
           enddo
           d(1:M)= iglob2rank(ind(1:M), MB, nranks)
           print '("  - global index  = ",24(I3))', ind(1:M)
           print '("  - local rank    = ",24(I3))', d(1:M)
           print '("  - local index   = ",24(I3))', &
                     iglob2loc(ind(1:M), MB, d(1:M), nranks)
        enddo
     enddo
     deallocate(ind)
  end subroutine test_numroc


  subroutine test_gemm(M, N, K, M_B, N_B, K_B, trA, trB, vpart)
  use mpi_f08
  use cuda_comm_aux
  use cublasmp_aux
  use time_lib, only: timef
  implicit none
  integer, intent(IN) :: M, N, K, M_B, N_B, K_B, trA, trB
  logical, intent(IN) :: vpart !< vertical partition?
  !
  integer, target     :: rank, nranks, localRank, localCommSize, ctxRank,&
                         ctxCommSize,  nDevLocal, color, ctx_color
  integer             :: err, info, sh(2), i, j, i1, nr, VRBZ, VRBZ0
  type(MPI_COMM), target   :: localComm, ctxComm
  type(c_ptr)         :: mp_grid, descrA, descrB, descrC
  integer :: MA,NA, MB,NB, MLA,NLA, MBA,NBA, MBB,NBB, MBC,NBC
  integer :: MLB,NLB, MLC,NLC, p,q
  double precision, allocatable, device, target :: d_A(:), d_B(:), d_C(:)
  double precision, allocatable :: h_A(:),h_A2(:,:), h_B(:),h_B2(:,:), h_C(:),h_C2(:,:)
  double precision             :: t1, t2, dt, nrm

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Get Hardware topological info
     ctx_color=11; VRBZ=2; VRBZ0=0
     call get_hw_topology(ctx_color, rank, nranks, localComm, localRank, &
                          localCommSize, ctxComm, ctxRank, ctxCommSize, &
                          nDevLocal, color, VRBZ_=VRBZ0)
     print '("[MAIN][MPI][get_hw_topology()]: rank="I4"| size="I4"'// &
           '| localRank="I4"| localCommSize="I4)', rank, nranks, &
           localRank, localCommSize
     ! Create CAL-COMM handle
     call cuda_comm_init(rank, nranks, localRank, VRBZ_=VRBZ0)

     call cublasmp_init(VRBZ_=VRBZ0)
     call show_cublasmp_handle


     !! ========================================================================
     !! vertical partitioning
     if (vpart) then
        p = nranks; q = 1
        if (trA.eq.0 .and. trB.eq.0) then ! ------------------------------------
           MA = M; NA = K
           MB = K; NB = N
           MBA = M_B; if (M_B == -1) MBA = (M - 1)/p + 1
           MBB = K_B; if (K_B == -1) MBB = (K - 1)/p + 1
           NBB = MBA  !< block size condition similar to NB = MA
           MBC = MBA  !< final matrix takes MB from matrix A,
           NBC = NBB  !<          and takes NB from matrix B
           NBA = MBB  !< O.K.: extra (undocumented) condition

        elseif (trA.eq.0 .and. trB.eq.1) then ! --------------------------------
           MA = M; NA = K
           MB = N; NB = K
           MBA = M_B; if (M_B == -1) MBA = (M - 1)/p + 1
           MBB = N_B; if (N_B == -1) MBB = (N - 1)/p + 1
           MBC = MBA; NBC = MBB
           NBA = MBA; NBB = MBA

        elseif (trA.eq.1 .and. trB.eq.0) then ! --------------------------------
           MA = K; NA = M
           MB = K; NB = N
           MBA = K_B; if (K_B == -1) MBA = (K - 1)/p + 1
           MBB = MBA
           NBA = MBA
           !NBA = M_B; if (M_B == -1) NBA = (M - 1)/p + 1
           NBB = MBA
           MBC = NBA
           NBC = MBA

        elseif (trA.eq.1 .and. trB.eq.1) then ! --------------------------------
print '("WARNING!!! THIS DOES NOT WORK IF N=/=M ")'
           MBA = K_B; if (K_B == -1) MBA = (K - 1)/p + 1
           NBB = MBA

           MBB = N_B; if (N_B == -1) MBB = (N - 1)/p + 1
           NBC = MBB

           NBA = M_B; if (M_B == -1) NBA = (M - 1)/p + 1
           MBC = NBA

           MA = K; NA = M
           MB = N; NB = K

        endif
     !=========================================================================
     else  ! horizontal partitioning a|b|...|z
        p = 1; q = nranks
        if (trA.eq.0 .and. trB.eq.0) then ! ------------------------------------
           MA = M; NA = K
           MB = K; NB = N
           NBA = K_B; if (K_B == -1) NBA = (K - 1)/q + 1
           NBB = N_B; if (N_B == -1) NBB = (N - 1)/q + 1
           MBA = NBB
           MBB = NBA
           MBC = MBA  !< final matrix takes MB from matrix A,
           NBC = NBB  !<          and takes NB from matrix B

        elseif (trA.eq.0 .and. trB.eq.1) then ! --------------------------------
           MA = M; NA = K
           MB = N; NB = K
           NBA = K_B; if (K_B == -1) NBA = (K - 1)/q + 1
           NBB = NBA
           MBA = NBA; MBB = NBB
           MBC = MBA; NBC = MBB

        elseif (trA.eq.1 .and. trB.eq.0) then ! --------------------------------
           MA = K; NA = M
           MB = K; NB = N
           NBA = M_B; if (M_B == -1) NBA = (M - 1)/q + 1
           NBB = N_B; if (N_B == -1) NBB = (N - 1)/q + 1
           MBB = NBA
           MBA = NBA
           MBC = NBA
           NBC = NBB

        elseif (trA.eq.1 .and. trB.eq.1) then ! --------------------------------
print '("WARNING!!! THIS DOES NOT WORK IF N=/=M ")'
           MA = K; NA = M
           MB = N; NB = K
           NBA = M_B; if (M_B == -1) NBA = (M - 1)/q + 1
           NBB = K_B; if (K_B == -1) NBB = (K - 1)/q + 1
           MBA = NBA
           MBB = NBB
           MBC = NBA
           NBC = MBB

        endif

     endif

     call MPI_BARRIER(MPI_COMM_WORLD)
     MLA = numroc(MA, MBA, rank, p)
     NLA = numroc(NA, NBA, rank, q)
     MLB = numroc(MB, MBB, rank, p)
     NLB = numroc(NB, NBB, rank, q)
     MLC = numroc(M,  MBC, rank, p)
     NLC = numroc(N,  NBC, rank, q)

     print '("["I2"/"I2"]: M,N,K="3(I3)" | MBA,NBA="2(I3)" |'// &
                       ' MBB,NBB="2(I3)" | MBC,NBC = "2(I3))', &
           rank,nranks,    M,N,K, MBA,NBA, MBB,NBB, MBC,NBC
     print '("["I2"/"I2"]: MLA,NLA="2(I3)" |'// &
                         ' MLB,NLB="2(I3)" | MLC,NLC = "2(I3))', &
           rank,nranks,    MLA,NLA, MLB,NLB, MLC,NLC

     ! initialize local part of matrices
     if (MLA*NLA > 0) then
        allocate(h_A2(MLA,NLA)); h_A2 = 0d0
        do j=1,NLA
           do i=1,MLA
              h_A2(i,j)= Aij(iloc2glob(i, MBA, rank, p), &
                             iloc2glob(j, NBA, rank, q))
           enddo
        enddo
        h_A = reshape(h_A2, (/MLA*NLA/))
     else
        allocate(h_A(1), h_A2(1,1))
        h_A2(1,1) = -756d0 !< poison value to indicate array is empty
     endif
     d_A = h_A
     call print_local_matrix(h_A2, "[MAIN] Matrix A:", rank, nranks)

     if (MLB*NLB > 0) then
        allocate(h_B2(MLB,NLB)); h_B2 = 0d0
        do j=1,NLB
           do i=1,MLB
              h_B2(i,j)= Bij(iloc2glob(i, MBB, rank, p), &
                             iloc2glob(j, NBB, rank, q))
           enddo
        enddo
        h_B = reshape(h_B2, (/MLB*NLB/))
     else
        allocate(h_B(1), h_B2(1,1))
        h_B2(1,1) = -756d0 !< poison value to indicate array is empty
     endif
     d_B = h_B
     call print_local_matrix(h_B2, "[MAIN] Matrix B:", rank, nranks)

     if (MLC*NLC > 0) then
        allocate(d_C(MLC*NLC))
     else
        allocate(d_C(1))
     endif
     deallocate(h_A, h_A2, h_B, h_B2)
     call MPI_BARRIER(MPI_COMM_WORLD)

     mp_grid = set_cublasMpGrid(p, q, VRBZ_=VRBZ)

     ! matrices A, B, C
     descrA = set_cublasMpMatrixDesc(mp_grid, MA, NA, MBA, NBA, p, rank, VRBZ_=VRBZ)
     descrB = set_cublasMpMatrixDesc(mp_grid, MB, NB, MBB, NBB, p, rank, VRBZ_=VRBZ)
     descrC = set_cublasMpMatrixDesc(mp_grid, M,  N,  MBC, NBC, p, rank, VRBZ_=VRBZ)

     if (trA.eq.0 .and. trB.eq.0) then ! ------------------------------------
        call check_status_mpi("[MAIN]: AxB", &
           AxB(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_=3))
     else if (trA.eq.0 .and. trB.eq.1) then
        call check_status_mpi("[MAIN]: AxBT", &
           AxBT(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_=3))
     else if (trA.eq.1 .and. trB.eq.0) then
        call check_status_mpi("[MAIN]: ATxB", &
           ATxB(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_=3))
     else
        call check_status_mpi("[MAIN]: ATxBT", &
           ATxBT(M, N, K, descrA, d_A, descrB, d_B, descrC, d_C, VRBZ_=3))
     endif

     h_C = d_C
     h_C2 = reshape(h_C, (/MLC, NLC/))
     call print_local_matrix(h_C2, "[MAIN] Matrix C:", rank, nranks)
     call check_matrix_product(h_C2, M, N, K, MBC, NBC, trA, trB, rank, p, q)
     deallocate(h_C, h_C2)

     ! Destroy CuSOLVERMP grids and matrix descriptors
     call destroy_cublasMpMatrixDesc(descrC, VRBZ_=VRBZ0)
     call destroy_cublasMpMatrixDesc(descrB, VRBZ_=VRBZ0)
     call destroy_cublasMpMatrixDesc(descrA, VRBZ_=VRBZ0)
     call destroy_cublasMpGrid(mp_grid, VRBZ_=VRBZ0)

     ! Destroy CuSOLVERMP handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cublasmp_destroy(VRBZ_=VRBZ0)

     ! Destroy CAL-COMM handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cal_comm_destroy(VRBZ_=VRBZ0)

     call MPI_FINALIZE(err)

  end subroutine test_gemm


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


  !!
  !! Copy a matrix A(M x N), distributed vertically, to a matrix B(M x N)
  !! distributed horizontally.
  !!
  subroutine test_gemr2d(MM, NN)
  use mpi_f08
  use cuda_comm_aux
  use cublasmp_aux
  implicit none
  integer, intent(IN) :: MM, NN
  !
  integer, target     :: rank, commSize, localRank, localCommSize, ctxRank,&
                         ctxCommSize,  nDevLocal, color, ctx_color
  integer             :: err, info, sh(2), i, j, i1, nr, VRBZ, VRBZ0
  type(MPI_COMM), target   :: localComm, ctxComm
  type(c_ptr)         :: mp_gridA, mp_gridB, descrA, descrB
  integer :: M,N, MA,NA, MB,NB, MMA,NNA, MMB,NNB, MBA,NBA, MBB,NBB, p,q
  double precision, allocatable, device, target :: d_A(:), d_B(:)
  double precision, allocatable :: h_A(:),h_A2(:,:), h_B(:),h_B2(:,:)

     ! Get Hardware topological info
     ctx_color=11; VRBZ=2; VRBZ0=0
     call get_hw_topology(ctx_color, rank, commSize, localComm, localRank, &
                          localCommSize, ctxComm, ctxRank, ctxCommSize, &
                          nDevLocal, color, VRBZ_=VRBZ0)
     print '("[MAIN][MPI][get_hw_topology()]: rank="I4"| size="I4"'// &
           '| localRank="I4"| localCommSize="I4)', rank, commSize, &
           localRank, localCommSize
     M = MM; N = NN
     p = commSize; q = 1

     ! Matrix A: M x N, partitioned vertically
     ! Matrix B: M x N, partitioned horizontally
     MBA = M/p
     if (p*MBA.lt.M) then
        print '("INFO: M is not divisible by p, padding with zeros")'
        MBA = MBA + 1
        M = p*MBA
     endif

     NBB = N/p
     if (p*NBB.lt.N) then
        print '("INFO: N is not divisible by p, padding with zeros")'
        NBB = NBB + 1
        N = p*NBB
     endif

     MA = M; MMA = MM
     NA = N; NNA = NN
     NBA = NBB

     MB = MA; MMB = MMA
     NB = NA; NNB = NNA
     MBB = MBA

     ! initialize local part of matrix A
     allocate(h_A2(MBA,NA)); h_A2 = 0d0
     do j = 1,NNA; do i = 1,MBA
        i1 = i + rank*MBA
        if (i1 <= MMA) h_A2(i,j)= Aij(i1,j)
     enddo; enddo
     h_A = reshape(h_A2, (/MBA*NA/))
     d_A = h_A
     call print_local_matrix(h_A2, "[test_gemr2d] Matrix A:", rank, commSize)

     ! initialize local part of matrix B
     allocate(h_B2(MB,NBB)); h_B2 = 0d0
     h_B = reshape(h_B2, (/MB*NBB/))
     d_B = h_B
     call print_local_matrix(h_B2, "[test_gemr2d] Matrix B before:", rank, commSize)

     call cuda_comm_init(rank, commSize, localRank, VRBZ_=VRBZ0)
     call cublasmp_init(VRBZ_=VRBZ0)
     call show_cublasmp_handle

     mp_gridA = set_cublasMpGrid(p, 1, VRBZ_=VRBZ)
     mp_gridB = set_cublasMpGrid(1, p, VRBZ_=VRBZ)

     ! matrix descriptors
     descrA = set_cublasMpMatrixDesc(mp_gridA, MA, NA, MBA, NBA, p, rank, VRBZ_=VRBZ)
     descrB = set_cublasMpMatrixDesc(mp_gridB, MB, NB, MBB, NBB, 1, 0, VRBZ_=VRBZ)

     call check_status_mpi("[test_gemr2d]: redistribute_rectMatrix", &
        redistribute_rectMatrix(M, N, descrA, d_A, descrB, d_B, VRBZ_=VRBZ))
     h_B = d_B
     h_B2 = reshape(h_B, [MB, NBB])
     call print_local_matrix(h_B2, "[test_gemr2d] Matrix B after:", rank, commSize)

     ! Destroy CuSOLVERMP grids and matrix descriptors
     call destroy_cublasMpMatrixDesc(descrA, VRBZ_=VRBZ0)
     call destroy_cublasMpMatrixDesc(descrB, VRBZ_=VRBZ0)
     call destroy_cublasMpGrid(mp_gridA, VRBZ_=VRBZ0)
     call destroy_cublasMpGrid(mp_gridB, VRBZ_=VRBZ0)

     ! Destroy CuSOLVERMP handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cublasmp_destroy(VRBZ_=VRBZ0)

     ! Destroy CAL-COMM handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cal_comm_destroy(VRBZ_=VRBZ0)

     call MPI_FINALIZE(err)

  end subroutine test_gemr2d


  !!
  !! Test matrix truncation: truncate a "tall" matrix A(M > N)
  !! to a square matrix
  !!
  subroutine test_truncate(MM, NN)
  use mpi_f08
  use cuda_comm_aux
  use cublasmp_aux
  implicit none
  integer, intent(IN) :: MM, NN
  !
  integer, target     :: rank, commSize, localRank, localCommSize, ctxRank,&
                         ctxCommSize,  nDevLocal, color, ctx_color
  integer             :: err, info, sh(2), i, j, i1, nr, VRBZ, VRBZ0
  type(MPI_COMM), target   :: localComm, ctxComm
  type(c_ptr)         :: mp_grid_v, mp_grid_h, descrA, descrA_trunc, descrB
  integer :: M,N, MA,NA, MMA,NNA, p,q, MB,MMB,MBA,NBA,MBB
  double precision, allocatable, device, target :: d_A(:), d_B(:)
  double precision, allocatable :: h_A(:),h_A2(:,:), h_B(:),h_B2(:,:)

     ! Get Hardware topological info
     ctx_color=11; VRBZ=2; VRBZ0=0
     call get_hw_topology(ctx_color, rank, commSize, localComm, localRank, &
                          localCommSize, ctxComm, ctxRank, ctxCommSize, &
                          nDevLocal, color, VRBZ_=VRBZ0)
     print '("[MAIN][MPI][get_hw_topology()]: rank="I4"| size="I4"'// &
           '| localRank="I4"| localCommSize="I4)', rank, commSize, &
           localRank, localCommSize
     M = MM; N = NN
     if (M < N) error stop "[!][test_truncate]: I need M > N"
     p = commSize; q = 1

     ! Matrix A: M x N, partitioned like v-stack, or vertically (=)
     ! Matrix B: N x N, truncated, partitioned like h-stack, horizontally (||)
     MBA = M/p
     if (p*MBA.lt.M) then
        print '("INFO: M is not divisible by p, padding with zeros")'
        MBA = MBA + 1
        M = p*MBA
     endif

     NBB = N/p
     if (p*NBB.lt.N) then
        print '("INFO: N is not divisible by p, padding with zeros")'
        NBB = NBB + 1
        N = p*NBB
     endif

     ! --------------
     ! A A A A      B B | B B
     ! A A A A      B B | B B
     ! A A A A      B B | B B
     ! -------  --> B B | B B
     ! A A A A
     ! A A A A
     ! A A A A
     !

     MA = M; MMA = MM
     NA = N; NNA = NN
     NBA = NBB

     MB = MA; MMB = MMA
     NB = NA; NNB = NNA
     MBB = MBA

     ! initialize local part of matrix A - partitioned like v-stack (-)
     allocate(h_A2(MBA,NA)); h_A2 = 0d0
     do j = 1,NNA; do i = 1,MBA
        i1 = i + rank*MBA
        if (i1 <= MMA) h_A2(i,j)= Aij(i1,j)
     enddo; enddo
     h_A = reshape(h_A2, [MBA*NA])
     d_A = h_A
     call print_local_matrix(h_A2, "[test_truncate] Matrix A:", rank, commSize)

     ! initialize local part of matrix B - partitioned like h-stack (.|.)
     allocate(h_B2(NB,NBB)); h_B2 = 0d0
     h_B = reshape(h_B2, [NB*NBB])
     d_B = h_B
     call print_local_matrix(h_B2, "[test_truncate] Matrix B before:", rank, commSize)

     call cuda_comm_init(rank, commSize, localRank, VRBZ_=VRBZ0)
     call cublasmp_init(VRBZ_=VRBZ0)
     call show_cublasmp_handle

     mp_grid_v= set_cublasMpGrid(p, 1, VRBZ_=VRBZ)
     mp_grid_h= set_cublasMpGrid(1, p, VRBZ_=VRBZ)

     ! matrix descriptors
     descrA = set_cublasMpMatrixDesc(mp_grid_v, MA, NA, MBA, NBA, p, rank, VRBZ_=VRBZ)
     descrB = set_cublasMpMatrixDesc(mp_grid_h, NA, NA, NBA, NBA, 1, 0, VRBZ_=VRBZ)
     descrA_trunc = set_cublasMpMatrixDesc(mp_grid_v, NA, NA, MBA, NBA, p, rank, VRBZ_=VRBZ)

     call check_status_mpi("[test_truncate]: redistribute_rectMatrix", &
        redistribute_rectMatrix(N, N, descrA_trunc, d_A, descrB, d_B, VRBZ_=VRBZ))
     h_B = d_B
     h_B2 = reshape(h_B, [NB, NBB])
     call print_local_matrix(h_B2, "[test_truncate] Matrix B after:", rank, commSize)

     ! Destroy CuSOLVERMP grids and matrix descriptors
     call destroy_cublasMpMatrixDesc(descrA_trunc, VRBZ_=VRBZ0)
     call destroy_cublasMpMatrixDesc(descrA, VRBZ_=VRBZ0)
     call destroy_cublasMpMatrixDesc(descrB, VRBZ_=VRBZ0)
     call destroy_cublasMpGrid(mp_grid_v, VRBZ_=VRBZ0)
     call destroy_cublasMpGrid(mp_grid_h, VRBZ_=VRBZ0)

     ! Destroy CuSOLVERMP handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cublasmp_destroy(VRBZ_=VRBZ0)

     ! Destroy CAL-COMM handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cal_comm_destroy(VRBZ_=VRBZ0)

     call MPI_FINALIZE(err)

  end subroutine test_truncate


  !!
  !! Test matrix redistribution
  !!
  subroutine test_redistrib(M, N, M_B, N_B)
  use mpi_f08
  use cuda_comm_aux
  use cublasmp_aux
  implicit none
  integer, intent(IN) :: M, N, M_B, N_B
  !
  character(*), parameter :: subnam = '[test_redistrib]'
  integer, target     :: rank, commSize, localRank, localCommSize, ctxRank,&
                         ctxCommSize,  nDevLocal, color, ctx_color
  integer             :: err, info, sh(2), i, j, i1, nr, VRBZ, VRBZ0
  type(MPI_COMM), target   :: localComm, ctxComm
  integer :: MA,NA, p,q, MLA,NLA, MLB,NLB
  double precision, allocatable, device, target :: d_A(:), d_B(:)
  double precision, allocatable :: h_A(:),h_A2(:,:), h_B(:),h_B2(:,:)
  type(c_ptr) :: gridA, descrA, gridB, descrB

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Get Hardware topology
     ctx_color=11; VRBZ=2; VRBZ0=0
     call get_hw_topology(ctx_color, rank, nranks, localComm, localRank, &
                          localCommSize, ctxComm, ctxRank, ctxCommSize, &
                          nDevLocal, color, VRBZ_=VRBZ0)
     print '("[MAIN][MPI][get_hw_topology()]: rank="I4"| size="I4"'// &
           '| localRank="I4"| localCommSize="I4)', rank, nranks, &
           localRank, localCommSize
     ! Create CAL-COMM handle
     call cuda_comm_init(rank, nranks, localRank, VRBZ_=VRBZ0)
     call cublasmp_init(VRBZ_=VRBZ0)
     call show_cublasmp_handle

     ! Matrix dimensions / block sizes
     p = nranks; q = 1
     MA = M_B; if (M_B == -1) MA = (M - 1)/nranks + 1
     NA = N_B; if (N_B == -1) NA = (N - 1)/nranks + 1

     ! Redistribute blocking: {M_B, 1} -> {1, N_B}
     MLA = numroc(M, MA, rank, nranks)
     NLA = numroc(N, 1,  rank, 1)
     MLB = numroc(M, 1,  rank, 1)
     NLB = numroc(N, NA, rank, nranks)

     print '("[",I2,"/",I2,"]: M,N=",2(I3)," | MB,NB=",2(I3))', &
              rank,nranks,     M,N,            MA,NA
     print '("[",I2,"/",I2,"]: MLA,NLA=",2(I3)," | MLB,NLB=",2(I3))', &
              rank,nranks,     MLA,NLA,        MLB,NLB

     ! initialize local part of matrices
     if (MLA*NLA > 0) then
        allocate(h_A2(MLA,NLA)); h_A2 = 0d0
        do j=1,NLA
           do i=1,MLA
              h_A2(i,j)=  100d0*iloc2glob(i, MA, rank, p) + &
                                iloc2glob(j, NA, rank, q)
           enddo
        enddo
        h_A = reshape(h_A2, (/MLA*NLA/))
     else
        allocate(h_A(1), h_A2(1,1))
        h_A2(1,1) = -756d0 !< poison value to indicate array is empty
     endif
     d_A = h_A
     call print_local_matrix(h_A2, "[MAIN] Matrix A:", rank, nranks, frmt_='F5.0')

     ! allocate d_B
     if (MLB*NLB > 0) then
        allocate(d_B(MLB*NLB))
     else
        allocate(d_B(1))
     endif

     ! create grids and descriptors
     gridA= set_cublasMpGrid(nranks, 1)
     gridB= set_cublasMpGrid(1, nranks)
     descrA= set_cublasMpMatrixDesc(gridA, M, N, MA, NA, nranks, rank, VRBZ_=VRBZ)
     descrB= set_cublasMpMatrixDesc(gridB, M, N, MA, NA, 1, rank, VRBZ_=VRBZ)
     call MPI_BARRIER(MPI_COMM_WORLD)

     ! redistribute
     err= redistribute_rectMatrix(M, N, descrA, d_A, descrB, d_B, VRBZ_=VRBZ)
     if (err /= 0) then
        print '("[!]'//subnam//': in redistribute_rectMatrix, err=",I5)', err
        error stop
     endif

     ! Print matrix B
     h_B = d_B
     h_B2 = reshape(h_B, [MLB, NLB])
     call print_local_matrix(h_B2, "[MAIN] Matrix B:", rank, nranks, frmt_='F5.0')

     ! Cleanup
     deallocate(d_A, d_B)
     deallocate(h_A, h_B, h_A2, h_B2)

     ! Destroy cuBLASMp grid and matrix descriptors
     call destroy_cublasMpMatrixDesc(descrA)
     call destroy_cublasMpMatrixDesc(descrB)
     call destroy_cublasMpGrid(gridA)
     call destroy_cublasMpGrid(gridB)

     ! Destroy cuBLASMp handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cublasmp_destroy(VRBZ_=VRBZ0)

     ! Destroy CAL-COMM handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cal_comm_destroy(VRBZ_=VRBZ0)

     call MPI_FINALIZE(err)

  end subroutine test_redistrib


  pure elemental double precision function Aij(i, j)
  integer, intent(IN) :: i, j
  Aij = dble(mod(j,7) + mod(i,3))/dble(1 + mod(i,5) + mod(j,11)) - 1d0
  !Aij = dble(10*i + j)
  end function


  pure elemental double precision function Bij(i, j)
  integer, intent(IN) :: i, j
  Bij = dble(mod(j,4) + mod(i,5))/dble(1 + mod(i,2) + mod(j,7)) - 1d0
  !Bij = dble(10*i + j)
  end function


  subroutine print_local_matrix(A2, label, rank, nranks, frmt_)
  use mpi_f08
  character(*),           intent(IN) :: label
  double precision,       intent(IN) :: A2(:,:)
  integer,                intent(IN) :: rank, nranks
  character(*), optional, intent(IN) :: frmt_
  !
  integer :: i, nr, sh(2)
  character(10) :: frmt
     frmt = 'ES10.3'; if (present(frmt_)) frmt = frmt_
     sh = shape(A2)
     if (sh(1).gt.20.or.sh(2).gt.20) return
     do nr = 1,nranks
        if (rank.eq.nr-1) then
           print '(A, "---------------")', label
           if (A2(1,1) /= -756d0) then
              do i = 1,sh(1)
                 print '(A,"["I2"/"I2"]: ",100('//frmt//',1X))', &
                        label, rank, nranks, A2(i,:)
              enddo
           endif
        endif
        call MPI_BARRIER(MPI_COMM_WORLD)
     enddo
  end subroutine


  subroutine check_matrix_product(C2, MM, NN, KK, &
                                      MB, NB, trA, trB, rank, p, q)
  use cublasmp_aux
  implicit none
  double precision, intent(IN) :: C2(:,:)
  integer, intent(IN)          :: MM, NN, KK, MB, NB, trA, trB
  integer, intent(IN)          :: rank, p, q
  !
  double precision :: err, Cij
  integer :: i, j, i1, j1, MB_C, sh(2), ML, NL

     err = .0d0
     sh = shape(C2)
     ML = numroc(MM, MB, rank, p)
     NL = numroc(NN, NB, rank, q)
     if (ML/=sh(1) .or. NL/=sh(2)) error stop "[!][check_matrix_product] shape"
     do j = 1,NL
        j1 = iloc2glob(j, NB, rank, q)
        do i = 1,ML
           i1 = iloc2glob(i, MB, rank, p)
           Cij = 0d0
           if     (trA.eq.0.and.trB.eq.0) then
              do k = 1,KK
                 Cij = Cij + Aij(i1,k)*Bij(k,j1)
              enddo
           elseif (trA.eq.0.and.trB.eq.1) then
              do k = 1,KK
                 Cij = Cij + Aij(i1,k)*Bij(j1,k)
              enddo
           elseif (trA.eq.1.and.trB.eq.0) then
              do k = 1,KK
                 Cij = Cij + Aij(k,i1)*Bij(k,j1)
              enddo
           else
              do k = 1,KK
                 Cij = Cij + Aij(k,i1)*Bij(j1,k)
              enddo
           endif
           err = err + (Cij - C2(i,j))**2
        enddo
     enddo
     print '("[+][check_matrix_product]["I2"/"I2"]: error = ",ES10.3)', &
            rank, max(p,q), dsqrt(err)

  end subroutine

end program main
