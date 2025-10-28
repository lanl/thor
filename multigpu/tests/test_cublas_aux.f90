!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_cublas_aux.f90
!! @author Ismael Boureima, Oleg Korobkin
!! @date  March 2025
!! @brief Testing cuBLAS
!!
program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: m = 10  ! in level-1 BLAS tests: array length
integer :: n = 10  ! (unused)
double precision, allocatable :: rnd(:)

  ! parse command-line arguments
  do i=1,command_argument_count()
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     elseif (trim(arg).eq.'-m') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) m
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
     elseif (trim(arg).eq.'-n') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) n
        if (err.ne.0) error stop "ERROR: when parsing -n <???>"
     else
        ! todo
     endif
  enddo

  call test_level1(m, n)

contains

  subroutine help_msg()
  print '(13(/,A))', &
   "Test CUDA BLAS subroutines", &
   "Usage: ./test.exe [-h] [-m <M>] [-n <N>]", &
   "Options:", &
   " -h           : print this help message", &
   " -m <M>       : in Level-1 tests: length of a vector A [10]", &
   " -n <N>       : number of columns in the matrix A [10]"
  end subroutine


  subroutine test_level1(MM, NN)
  use mpi_f08
  use cuda_comm_aux
  use cublas_aux
  use mat_lib, only: normfro
  implicit none
  integer, intent(IN) :: MM, NN
  !
  integer, target     :: rank, commSize, localRank, localCommSize, ctxRank,&
                         ctxCommSize,  nDevLocal, color, ctx_color
  integer             :: err, info, sh(2), i, j, i0, j0, i1, nr, VRBZ, VRBZ0, ival
  type(MPI_COMM), target   :: localComm, ctxComm
  type(c_ptr)         :: mp_grid, descrA, descrB, descrC
  integer :: M,N,K, MA,NA, MB,NB, MMA,NNA, MMB,NNB, MBA,NBA, MBB,NBB, MBC,NBC, p,q
  double precision, allocatable, device, target :: d_A(:), d_B(:), d_C(:), d_A2(:,:), d_B2(:,:)
  double precision, allocatable :: h_A(:),h_A2(:,:), h_B(:),h_B2(:,:), h_C(:),h_C2(:,:)
  double precision :: dval

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
     call MPI_BARRIER(MPI_COMM_WORLD)

     ! initialize vectors A and B
     allocate(h_A(M),h_B(M))
     do i = 1, M
        h_A(i)= dble(i)*0.1d0
        h_B(i)= dble(i)*10d0
     enddo
     d_A = h_A
     d_B = h_B
     deallocate(h_A, h_B)

     ! initialize vectors A and B
     allocate(h_A2(M,N),h_B2(2*M,2*N))
     h_B2= 0d0
     do j = 1, N
        do i = 1, M
           h_A2(i,j)= 1d1*dble(i) + dble(j)
        enddo
     enddo
     d_A2 = h_A2
     d_B2 = h_B2

     deallocate(h_A2, h_B2)
     call MPI_BARRIER(MPI_COMM_WORLD)

     call cuda_comm_init(rank, commSize, localRank, VRBZ_=VRBZ0)
     call setup_cublas(VRBZ0)

     ival= cublas_iamax(d_A)
     dval= cublas_maxval(d_A)
     print '("imax = ",I8,", maxval = ",ES12.5)', ival,dval

     ival= cublas_iamin(d_A)
     dval= cublas_minval(d_A)
     print '("imin = ",I8,", minval = ",ES12.5)', ival,dval

     dval= cublas_asum(d_A)
     print '("asum = ",ES12.5)', dval

     call cublas_axpy(-100d0, d_A, d_B)
     print '("|100*A - B|_2 = ",ES12.5)', cublas_nrm2(d_B)

     call cublas_copy(d_A, d_B)
     h_B = d_B
     print '("B = ",100(ES12.5,1X))', h_B
     deallocate(h_B)

     dval= cublas_dot(d_A, d_B)
     print '("(A, B) = ",ES12.5)', dval

     dval= cublas_nrm2(d_A)
     print '("nrm2^2 = ",ES12.5)', dval**2

     call cublas_scal(-1d0, d_A)
     h_A = d_A
     print '("A = ",100(ES12.5,1X))', h_A
     deallocate(h_A)

     call cublas_swap(d_A, d_B)
     h_A=d_A; h_B=d_B
     print '("A = ",100(ES12.5,1X))', h_A
     print '("B = ",100(ES12.5,1X))', h_B
     deallocate(d_A,d_B)

     i0 = max(MM/2,1); j0 = max(NN/2,1)
     call cublas_copy_2d(d_A2, d_B2, yi0_=i0, yj0_=j0)
     h_B2 = d_B2
     h_A2 = d_A2
     h_B2(i0:MM+i0-1,j0:NN+j0-1)= h_B2(i0:MM+i0-1,j0:NN+j0-1) - h_A2(1:MM,1:NN)
     print '("|B2 - A2|_2 = ",ES12.5)', normfro(h_B2)
     ! Cleanup
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cleanup_cublas(VRBZ0)

     ! Destroy CAL-COMM handle
     call MPI_BARRIER(MPI_COMM_WORLD)
     call cal_comm_destroy(VRBZ_=VRBZ0)

     call MPI_FINALIZE(err)

  end subroutine test_level1


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


  double precision function Aij(i, j)
  integer, intent(IN) :: i, j
  Aij = dble(mod(j,7) + mod(i,3))/dble(1 + mod(i,5) + mod(j,11)) - 1d0
  !Aij = dble(10*i + j)
  end function


  double precision function Bij(i, j)
  integer, intent(IN) :: i, j
  Bij = dble(mod(j,4) + mod(i,5))/dble(1 + mod(i,2) + mod(j,7)) - 1d0
  !Bij = dble(10*i + j)
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

  subroutine check_matrix_product(C2, MM, NN, KK, trA, trB, rank, nranks)
  double precision, intent(IN) :: C2(:,:)
  integer, intent(IN)          :: MM, NN, KK, trA, trB
  integer, intent(IN)          :: rank, nranks
  !
  double precision :: err, Cij
  integer :: i, j, i1, MB_C, sh(2)

     err = .0d0
     sh = shape(C2)
     MB_C= sh(1)
     do j = 1,NN; do i = 1,MB_C
        i1 = i + rank*MB_C
        if (i1.gt.MM) cycle
        Cij = 0d0
        if     (trA.eq.0.and.trB.eq.0) then
           do k = 1,KK
              Cij = Cij + Aij(i1,k)*Bij(k,j)
           enddo
        elseif (trA.eq.0.and.trB.eq.1) then
           do k = 1,KK
              Cij = Cij + Aij(i1,k)*Bij(j,k)
           enddo
        elseif (trA.eq.1.and.trB.eq.0) then
           do k = 1,KK
              Cij = Cij + Aij(k,i1)*Bij(k,j)
           enddo
        else
           do k = 1,KK
              Cij = Cij + Aij(k,i1)*Bij(j,k)
           enddo
        endif
        err = err + (Cij - C2(i,j))**2
     enddo; enddo
     print '("[+][check_matrix_product]["I2"/"I2"]: error = ",ES10.3)', &
            rank, nranks, dsqrt(err)

  end subroutine

end program main
