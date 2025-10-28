!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file mat_ops_kernels.f90
!! @author Ismael Djibrilla Boureima, Oleg Korobkin
!! @date  May 2024
!! @brief CUDA kernels for matrix mainupulation
!!
module matrix_ops_kernels
  implicit none
contains


  ! simple copy kernel
  !
  ! used as reference case representing best
  ! effictive bandwidth

  attributes(global) subroutine copy(idata, odata, nx, ny, TILE_DIM, BLOCK_ROWS)
    implicit none
    integer, value          :: nx, ny, TILE_DIM, BLOCK_ROWS
    double precision, intent(in) :: idata(nx,ny)
    double precision, intent(out) :: odata(nx,ny)
    integer :: x, y, j

    x = (blockIdx%x-1) * TILE_DIM + threadIdx%x
    y = (blockIdx%y-1) * TILE_DIM + threadIdx%y
    if (x > nx .or. y > ny) return

    do j = 0, TILE_DIM-1, BLOCK_ROWS
       if (y + j > ny) exit
       odata(x,y+j) = idata(x,y+j)
    end do
  end subroutine copy

  ! copy kernel using shared memory
  !
  ! also used as reference case, demonstrating effect of
  ! using shared memory

  attributes(global) subroutine copySharedMem(idata, odata, nx, ny, TILE_DIM, BLOCK_ROWS)
    implicit none
    integer, value        :: nx, ny, TILE_DIM, BLOCK_ROWS
    double precision, intent(in) :: idata(nx,ny)
    double precision, intent(out) :: odata(nx,ny)
    double precision, shared :: tile(TILE_DIM, TILE_DIM)
    integer :: x, y, j

    x = (blockIdx%x-1) * TILE_DIM + threadIdx%x
    y = (blockIdx%y-1) * TILE_DIM + threadIdx%y

    do j = 0, TILE_DIM-1, BLOCK_ROWS
       tile(threadIdx%x, threadIdx%y+j) = idata(x,y+j)
    end do

    call syncthreads()

    do j = 0, TILE_DIM-1, BLOCK_ROWS
       odata(x,y+j) = tile(threadIdx%x, threadIdx%y+j)
    end do
  end subroutine copySharedMem

  ! naive transpose
  !
  ! simplest transpose - doesn't use shared memory
  ! reads from global memory are coalesced but not writes

  attributes(global) subroutine transposeNaive(idata, odata, nx, ny, TILE_DIM, BLOCK_ROWS)
    implicit none
    integer, value               :: nx, ny, TILE_DIM, BLOCK_ROWS
    double precision, intent(in) :: idata(nx,ny)
    double precision, intent(out) :: odata(ny,nx)
    integer :: x, y, j

    x = (blockIdx%x-1) * TILE_DIM + threadIdx%x
    y = (blockIdx%y-1) * TILE_DIM + threadIdx%y
    if (x > nx .or. y > ny) return

    transpose_loop: do j = 0, TILE_DIM-1, BLOCK_ROWS
       if (y + j > ny) exit transpose_loop
       odata(y+j,x) = idata(x,y+j)
    enddo transpose_loop
  end subroutine transposeNaive


  attributes(global) subroutine transpose1D(A, A_T, M, N)
  implicit none
  integer, value                  :: M, N
  double precision, intent(in)    :: A(M*N)
  double precision, intent(inout) :: A_T(M*N)
  integer :: i, i_t, off
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    off = i - 1
    i_t = 1 + (mod(off,M)*N) + (off/M)
    if (i <= M*N) A_T(i_t) = A(i)
  end subroutine transpose1D


  !> transpose last two dimensions of a 3D array
  !!
  !! Given A(:) that is a flat representation of a 3D array A3(:,:,:),
  !! compute A3_T s.t. A3_T(:,i,j) = A3(:,j,i)
  !!
  attributes(global) subroutine transpose_flat3d_132(A, A_T, M, N, L)
  implicit none
  integer, value                  :: M, N, L
  double precision, intent(in)    :: A(M*N*L)
  double precision, intent(inout) :: A_T(M*N*L)
  integer :: i,j,k,ijk,ijk_t
    ijk = blockDim%x*(blockIdx%x - 1) + threadIdx%x - 1
    i = mod(ijk, M)
    j = mod((ijk/M), N)
    k = (ijk/M)/N
    ijk_t = 1 + i + M*(k + L*j)
    if (ijk_t <= M*N*L .and. ijk < M*N*L) A_T(ijk_t) = A(ijk + 1)
  end subroutine transpose_flat3d_132


  ! coalesced transpose
  !
  ! uses shared memory to achieve coalesing in both reads
  ! and writes
  !
  ! tile size causes shared memory bank conflicts

  attributes(global) subroutine transposeCoalesced(idata, odata, nx, ny, TILE_DIM, BLOCK_ROWS)
    implicit none
    integer, value          :: nx, ny, TILE_DIM, BLOCK_ROWS
    double precision, intent(in) :: idata(nx,ny)
    double precision, intent(out) :: odata(ny,nx)
    double precision, shared :: tile(TILE_DIM, TILE_DIM)
    integer :: x, y, j

    x = (blockIdx%x-1) * TILE_DIM + threadIdx%x
    y = (blockIdx%y-1) * TILE_DIM + threadIdx%y

    do j = 0, TILE_DIM-1, BLOCK_ROWS
       tile(threadIdx%x, threadIdx%y+j) = idata(x,y+j)
    end do

    call syncthreads()

    x = (blockIdx%y-1) * TILE_DIM + threadIdx%x
    y = (blockIdx%x-1) * TILE_DIM + threadIdx%y

    do j = 0, TILE_DIM-1, BLOCK_ROWS
       odata(x,y+j) = tile(threadIdx%y+j, threadIdx%x)
    end do
  end subroutine transposeCoalesced

  ! no bank-conflict transpose
  !
  ! same as transposeCoalesced except the first tile dimension is padded
  ! to avoid shared memory bank conflicts

  attributes(global) subroutine transposeNoBankConflicts(idata, odata, nx, ny, TILE_DIM, BLOCK_ROWS)
    implicit none
    integer, value               :: nx, ny, TILE_DIM, BLOCK_ROWS
    double precision, intent(in) :: idata(nx,ny)
    double precision, intent(out) :: odata(ny,nx)
    double precision, shared :: tile(TILE_DIM+1, TILE_DIM)
    integer :: x, y, j

    x = (blockIdx%x-1) * TILE_DIM + threadIdx%x
    y = (blockIdx%y-1) * TILE_DIM + threadIdx%y

    do j = 0, TILE_DIM-1, BLOCK_ROWS
       tile(threadIdx%x, threadIdx%y+j) = idata(x,y+j)
    end do

    call syncthreads()

    x = (blockIdx%y-1) * TILE_DIM + threadIdx%x
    y = (blockIdx%x-1) * TILE_DIM + threadIdx%y

    do j = 0, TILE_DIM-1, BLOCK_ROWS
       odata(x,y+j) = tile(threadIdx%y+j, threadIdx%x)
    end do
  end subroutine transposeNoBankConflicts


  !> kernel to copy 2D device arrays: Y(i1:i2,j1:j2)= X(p1:p2,q1:q2)
  attributes(global) subroutine cu_copy_arr2d(X,Y, i1,i2,j1,j2, p1,q1)
  implicit none
  double precision, intent(in)    :: X(:,:) !< part of the matrix on this rank
  double precision, intent(inout) :: Y(:,:) !< part of the matrix on this rank
  integer, value                  :: i1,i2,j1,j2,p1,q1
  !
  integer :: j, k, p2 !< global index of the matrix column

     j = (blockIdx%y - 1)*blockDim%y + threadIdx%y - 1
     p2 = p1 + i2 - i1
     if (j <= j2 - j1) Y(p1:p2,q1+j)= X(i1:i2,j1+j)

  end subroutine cu_copy_arr2d


  !> kernel for aX+Y op on 2D arrays: Y(i1:i2,j1:j2)+= a*X(p1:p2,q1:q2)
  attributes(global) subroutine cu_axpy_arr2d(alpha,X,Y, i1,i2,j1,j2, p1,q1)
  implicit none
  double precision, value        :: alpha
  double precision,   intent(IN) :: X(:,:)
  double precision,intent(INOUT) :: Y(:,:)
  integer, value                 :: i1,i2,j1,j2,p1,q1
  !
  integer :: j, p2 !< global index of the matrix column

     j = (blockIdx%y - 1)*blockDim%y + threadIdx%y - 1
     p2 = p1 + i2 - i1
     if (j <= j2 - j1) &
        Y(p1:p2,q1+j)= Y(p1:p2,q1+j) + alpha*X(i1:i2,j1+j)

  end subroutine cu_axpy_arr2d


  !> kernel to copy 3D device arrays: Y(i1:i2,j1:j2,k1:k2)= X(p1:p2,q1:q2,r1:r2)
  attributes(global) subroutine cu_copy_arr3d(X,Y, i1,i2,j1,j2,k1,k2, p1,q1,r1)
  implicit none
  double precision, intent(in)    :: X(:,:,:) !< part of the matrix on this rank
  double precision, intent(inout) :: Y(:,:,:) !< part of the matrix on this rank
  integer, value                  :: i1,i2,j1,j2,k1,k2,p1,q1,r1
  !
  integer :: j, k, p2 !< global index of the matrix column

     j = (blockIdx%y - 1)*blockDim%y + threadIdx%y - 1
     k = (blockIdx%z - 1)*blockDim%z + threadIdx%z - 1
     p2 = p1 + i2 - i1
     if (j <= j2 - j1 .and. k <= k2 - k1) Y(p1:p2,q1+j,r1+k)= X(i1:i2,j1+j,k1+k)

  end subroutine cu_copy_arr3d


  !> kernel for aX+Y op on 3D arrays: Y(i1:i2,j1:j2,k1:k2)+= a*X(p1:p2,q1:q2,r1:r2)
  attributes(global) subroutine cu_axpy_arr3d(alpha,X,Y, i1,i2,j1,j2,k1,k2, p1,q1,r1)
  implicit none
  double precision, value        :: alpha
  double precision,   intent(IN) :: X(:,:,:)
  double precision,intent(INOUT) :: Y(:,:,:)
  integer, value                 :: i1,i2,j1,j2,k1,k2,p1,q1,r1
  !
  integer :: j, k, p2 !< global index of the matrix column

     j = (blockIdx%y - 1)*blockDim%y + threadIdx%y - 1
     k = (blockIdx%z - 1)*blockDim%z + threadIdx%z - 1
     p2 = p1 + i2 - i1
     if (j <= j2 - j1 .and. k <= k2 - k1) &
        Y(p1:p2,q1+j,r1+k)= Y(p1:p2,q1+j,r1+k) + alpha*X(i1:i2,j1+j,k1+k)

  end subroutine cu_axpy_arr3d


  !> Use local index to compute global index in a distributed matrix
  attributes(device) &
  integer function cu_iloc2glob(j, MB, rank, nranks)
  integer, value :: j, MB, rank, nranks
     cu_iloc2glob = 1 + MB*(mod(rank,nranks) + ((j - 1)/MB)*nranks) &
                  + mod(j - 1, MB)
  end function


  !> kernel to zero out elements below matrix diagonal
  attributes(global) subroutine &
  cu_zero_below_diag(locA, M, N, MB, NB, ML, NL, iB, jB, nRowDevs, nColDevs)
  implicit none
  double precision, intent(inout) :: locA(:, :) !< part of the matrix on this rank
  integer, value                  :: M, N       !< global size of the matrix (rows, columns)
  integer, value                  :: MB, NB     !< blocking size of the matrix (rows, columns)
  integer, value                  :: ML, NL     !< local size of the matrix (rows, columns)
  integer, value                  :: iB, jB     !< local matrix index in the global grid
  integer, value                  :: nRowDevs, nColDevs !< number of devices in the grid

  integer :: i, j   !< global indices of the matrix element
  integer :: il, jl

     il = (blockIdx%x - 1)*blockDim%x + threadIdx%x
     jl = (blockIdx%y - 1)*blockDim%y + threadIdx%y

     i = cu_iloc2glob(il, MB, iB, nRowDevs)
     j = cu_iloc2glob(jl, NB, jB, nColDevs)

     if (il<= ML .and. jl<= NL.and. i > j) locA(il,jl) = 0d0

  end subroutine cu_zero_below_diag


  !> kernel to zero out elements above matrix diagonal
  attributes(global) subroutine &
  cu_zero_above_diag(locA, M, N, MB, NB, ML, NL, iB, jB, nRowDevs, nColDevs)
  implicit none
  double precision, intent(inout) :: locA(:, :) !< part of the matrix on this rank
  integer, value                  :: M, N       !< global size of the matrix (rows, columns)
  integer, value                  :: MB, NB     !< blocking size of the matrix (rows, columns)
  integer, value                  :: ML, NL     !< local size of the matrix (rows, columns)
  integer, value                  :: iB, jB     !< local matrix index in the global grid
  integer, value                  :: nRowDevs, nColDevs !< number of devices in the grid

  integer :: i, j   !< global indices of the matrix element
  integer :: il, jl

     il = (blockIdx%x - 1)*blockDim%x + threadIdx%x
     jl = (blockIdx%y - 1)*blockDim%y + threadIdx%y

     i = cu_iloc2glob(il, MB, iB, nRowDevs)
     j = cu_iloc2glob(jl, NB, jB, nColDevs)

     if (il<= ML .and. jl<= NL.and. i < j) locA(il,jl) = 0d0

  end subroutine cu_zero_above_diag


  !> kernel to create an identity matrix (rectangular)
  attributes(global) subroutine &
  cu_eye(locA, M, N, MB, NB, ML, NL, iB, jB, nRowDevs, nColDevs)
  implicit none
  double precision, intent(inout) :: locA(:) !< part of the matrix on this rank
  integer, value :: M, N    !< global size of the matrix (rows, columns)
  integer, value :: MB, NB  !< blocking size of the matrix (rows, columns)
  integer, value :: ML, NL  !< local size of the matrix (rows, columns)
  integer, value :: iB, jB  !< local matrix index in the global grid
  integer, value :: nRowDevs, nColDevs !< number of devices in the grid

  integer :: i, j   !< global indices of the matrix element
  integer :: il, jl

     il = (blockIdx%x - 1)*blockDim%x + threadIdx%x
     jl = (blockIdx%y - 1)*blockDim%y + threadIdx%y

     if (il<= ML .and. jl<= NL) then
        i = cu_iloc2glob(il, MB, iB, nRowDevs)
        j = cu_iloc2glob(jl, NB, jB, nColDevs)
        locA(il + ML*(jl-1)) = i.eq.j
     endif

  end subroutine cu_eye


  !> kernel to add to [sub|super]duaginal of a distributed matrix
  attributes(global) subroutine &
  cu_add_matrix_diag(locA, diag, ndiag, M, N, MB, NB, ML, NL, iB, jB, nRowDevs, nColDevs)
  implicit none
  double precision, intent(inout) :: locA(:) !< part of the matrix on this rank
  double precision,    intent(in) :: diag(:) !< part of the matrix on this rank
  integer, value :: ndiag   !< which diagonal: 0=main, +1=first super, -1=first sub etc. 
  integer, value :: M, N    !< global size of the matrix (rows, columns)
  integer, value :: MB, NB  !< blocking size of the matrix (rows, columns)
  integer, value :: ML, NL  !< local size of the matrix (rows, columns)
  integer, value :: iB, jB  !< local matrix index in the global grid
  integer, value :: nRowDevs, nColDevs !< number of devices in the grid

  integer :: i, j   !< global indices of the matrix element
  integer :: il, jl

     il = (blockIdx%x - 1)*blockDim%x + threadIdx%x
     jl = (blockIdx%y - 1)*blockDim%y + threadIdx%y

     if (il<= ML .and. jl<= NL) then
        i = cu_iloc2glob(il, MB, iB, nRowDevs)
        j = cu_iloc2glob(jl, NB, jB, nColDevs)
        if (j - i .eq. ndiag) locA(il + ML*(jl-1)) = diag(min(i,j))
     endif

  end subroutine cu_add_matrix_diag


  !> kernel that shifts a 2D array to the left by k positions
  attributes(global) subroutine &
  cu_shift_left_arr2d(A, M, N, k)
  implicit none
  double precision, intent(inout) :: A(:,:) !< part of the matrix on this rank
  integer, value                  :: M, N   !< array dimensions (rows, columns)
  integer, value                  :: k      !< shift by this many positions
  !
  integer :: i, j   !< global indices of the matrix element

     i = (blockIdx%x - 1)*blockDim%x + threadIdx%x
     if (i<= M) then
        do j=1,N-k
           A(i,j) = A(i,j+k)
        enddo
        A(i,N-k+1:N) = 0d0
     endif

  end subroutine cu_shift_left_arr2d


  !> kernel that shifts a 2D array to the right by k positions
  attributes(global) subroutine &
  cu_shift_right_arr2d(A, M, N, k)
  implicit none
  double precision, intent(inout) :: A(:, :) !< part of the matrix on this rank
  integer, value                  :: M, N    !< array dimensions (rows, columns)
  integer, value                  :: k       !< shift by this many positions
  !
  integer :: i, j   !< global indices of the matrix element

     i = (blockIdx%x - 1)*blockDim%x + threadIdx%x
     if (i<= M) then
        do j=N,k+1,-1
           A(i,j) = A(i,j-k)
        enddo
        A(i,1:k) = 0d0
     endif

  end subroutine cu_shift_right_arr2d


  !> kernel to shift a 2D array upward by k positions
  attributes(global) subroutine &
  cu_shift_up_arr2d(A, M, N, k)
  implicit none
  double precision, intent(inout) :: A(:,:) !< part of the matrix on this rank
  integer, value                  :: M, N   !< array dimensions (rows, columns)
  integer, value                  :: k      !< shift by this many positions
  !
  integer :: i, j   !< global indices of the matrix element

     j = (blockIdx%y - 1)*blockDim%y + threadIdx%y
     if (j<= N) then
        do i=1,M-k
           A(i,j) = A(i+k,j)
        enddo
        A(M-k+1:M,j) = 0d0
     endif

  end subroutine cu_shift_up_arr2d


  !> kernel to shifts a 2D array downward by k positions
  attributes(global) subroutine &
  cu_shift_down_arr2d(A, M, N, k)
  implicit none
  double precision, intent(inout) :: A(:,:) !< part of the matrix on this rank
  integer, value                  :: M, N   !< array dimensions (rows, columns)
  integer, value                  :: k      !< shift by this many positions
  !
  integer :: i, j   !< global indices of the matrix element

     j = (blockIdx%y - 1)*blockDim%y + threadIdx%y
     if (j<= N) then
        do i=M,k+1,-1
           A(i,j) = A(i-k,j)
        enddo
        A(1:k,j) = 0d0
     endif

  end subroutine cu_shift_down_arr2d


end module matrix_ops_kernels

!
!
! MAT/VECT OPS
!
!

module matrix_ops
  use cudafor
  use matrix_ops_kernels
  implicit none

contains

  subroutine transpose_device_array_vect(Arr, transp_Arr, M,N, block_threads_)
  implicit none
  integer, intent(IN)                     :: M, N
  double precision, device, intent(IN)    :: Arr(M*N)
  double precision, device, intent(INOUT) :: transp_Arr(M*N)
  integer, intent(IN), optional           :: block_threads_
  integer                                 :: block_threads
  type(dim3)                              :: grid, tBlock
  integer :: NN

     block_threads = 256; if(present(block_threads_)) block_threads = block_threads_
     NN = M*N
     tBlock = dim3(block_threads_,1,1)
     grid = dim3(ceiling(real(NN)/tBlock%x),1,1)
     call transpose1D<<<grid, tBlock>>>(Arr, transp_Arr, M,N)

  end subroutine transpose_device_array_vect


  subroutine transpose_device_flat3d_132(Arr, transp_Arr, M, N, L, block_threads_)
  implicit none
  integer, intent(in)                     :: M, N, L
  double precision, device, intent(in)    :: Arr(M*N*L)
  double precision, device, intent(inout) :: transp_Arr(M*N*L)
  integer, intent(in), optional           :: block_threads_
  !
  integer                                 :: block_threads
  type(dim3)                              :: grid, tBlock
  integer :: NN

     block_threads = 256; if(present(block_threads_)) block_threads = block_threads_
     NN = M*N*L
     tBlock = dim3(block_threads,1,1)
     grid = dim3((NN - 1)/tBlock%x + 1,1,1)
     call transpose_flat3d_132<<<grid, tBlock>>>(Arr, transp_Arr, M, N, L)

  end subroutine transpose_device_flat3d_132


  subroutine transpose_device_array_naive(Arr, transp_Arr, TILE_DIM, BLOCK_ROWS)
     implicit none
     integer, intent(IN)                   :: TILE_DIM, BLOCK_ROWS
     double precision, device, intent(IN)  :: Arr(:,:)
     double precision, device, intent(INOUT) :: transp_Arr(:,:)
     type (dim3)                          :: dimGrid, dimBlock
     integer :: i, j, istat, sh(2), nx, ny
     sh = shape(Arr)
     nx = sh(1); ny = sh(2)
     ! check parameters and calculate execution configuration
     if (mod(nx, TILE_DIM) /= 0 .or. mod(ny, TILE_DIM) /= 0) then
        write(*,*) 'nx and ny must be a multiple of TILE_DIM'
        stop
     end if
     if (mod(TILE_DIM, BLOCK_ROWS) /= 0) then
        write(*,*) 'TILE_DIM must be a multiple of BLOCK_ROWS'
        stop
     end if
     dimGrid = dim3(nx/TILE_DIM, ny/TILE_DIM, 1)
     dimBlock = dim3(TILE_DIM, BLOCK_ROWS, 1)
     call transposeNaive<<<dimGrid, dimBlock>>>(Arr, transp_Arr, nx, ny, TILE_DIM, BLOCK_ROWS)
  end subroutine transpose_device_array_naive

  subroutine transpose_device_array_coalesced(Arr, transp_Arr, TILE_DIM, BLOCK_ROWS)
     implicit none
     integer, intent(IN)                   :: TILE_DIM, BLOCK_ROWS
     double precision, device, intent(IN)  :: Arr(:,:)
     double precision, allocatable, device, intent(INOUT) :: transp_Arr(:,:)
     type (dim3)                          :: dimGrid, dimBlock
     integer :: i, j, istat, sh(2), nx, ny
     sh = shape(Arr)
     nx = sh(1); ny = sh(2)
     ! check parameters and calculate execution configuration
     if (mod(nx, TILE_DIM) /= 0 .or. mod(ny, TILE_DIM) /= 0) then
        write(*,*) 'nx and ny must be a multiple of TILE_DIM'
        stop
     end if
     if (mod(TILE_DIM, BLOCK_ROWS) /= 0) then
        write(*,*) 'TILE_DIM must be a multiple of BLOCK_ROWS'
        stop
     end if
     dimGrid = dim3(nx/TILE_DIM, ny/TILE_DIM, 1)
     dimBlock = dim3(TILE_DIM, BLOCK_ROWS, 1)
     call transposeCoalesced<<<dimGrid, dimBlock>>>(Arr, transp_Arr, nx, ny, TILE_DIM, BLOCK_ROWS)
  end subroutine transpose_device_array_coalesced

  subroutine transpose_device_array_NoBankConflict(Arr, transp_Arr, TILE_DIM, BLOCK_ROWS)
     implicit none
     integer, intent(IN)                   :: TILE_DIM, BLOCK_ROWS
     double precision, device, intent(IN)  :: Arr(:,:)
     double precision, allocatable, device, intent(INOUT) :: transp_Arr(:,:)
     type (dim3)                          :: dimGrid, dimBlock
     integer :: i, j, istat, sh(2), nx, ny
     sh = shape(Arr)
     nx = sh(1); ny = sh(2)
     ! check parameters and calculate execution configuration
     if (mod(nx, TILE_DIM) /= 0 .or. mod(ny, TILE_DIM) /= 0) then
        write(*,*) 'nx and ny must be a multiple of TILE_DIM'
        stop
     end if
     if (mod(TILE_DIM, BLOCK_ROWS) /= 0) then
        write(*,*) 'TILE_DIM must be a multiple of BLOCK_ROWS'
        stop
     end if
     dimGrid = dim3(nx/TILE_DIM, ny/TILE_DIM, 1)
     dimBlock = dim3(TILE_DIM, BLOCK_ROWS, 1)
     call transposeNoBankConflicts<<<dimGrid, dimBlock>>>(Arr, transp_Arr, nx, ny, TILE_DIM, BLOCK_ROWS)
   end subroutine transpose_device_array_NoBankConflict

  !function transpose_device_array(Arr, TILE_DIM, BLOCK_ROWS) result(transp_Arr)
  subroutine transpose_device_array(Arr, transp_Arr, TILE_DIM, BLOCK_ROWS)
     implicit none
     integer, intent(IN)                   :: TILE_DIM, BLOCK_ROWS
     double precision, device, intent(IN)  :: Arr(:,:)
     double precision, allocatable, device, intent(INOUT) :: transp_Arr(:,:)
     !double precision, allocatable         :: Hrr(:,:), transp_Hrr(:,:)
     type (dim3)                          :: dimGrid, dimBlock
     integer :: i, j, istat, sh(2), nx, ny
     sh = shape(Arr)
     nx = sh(1); ny = sh(2)
     !allocate(transp_Hrr(ny, nx))
     !Hrr = Arr
     !print*, "[+] Hrr(100,200) = ",  Hrr(100,200)
     !print*, "[+] Hrr(200,100) = ",  Hrr(200,100)
     ! check parameters and calculate execution configuration
     if (mod(nx, TILE_DIM) /= 0 .or. mod(ny, TILE_DIM) /= 0) then
        write(*,*) 'nx and ny must be a multiple of TILE_DIM'
        stop
     end if
     if (mod(TILE_DIM, BLOCK_ROWS) /= 0) then
        write(*,*) 'TILE_DIM must be a multiple of BLOCK_ROWS'
        stop
     end if
     dimGrid = dim3(nx/TILE_DIM, ny/TILE_DIM, 1)
     dimBlock = dim3(TILE_DIM, BLOCK_ROWS, 1)
     !
     ! Various transpose functions are taken from:
     !
     ! https://github.com/NVIDIA-developer-blog/code-samples/blob/master/series/cuda-fortran/transpose.cuf
     !
     !write(*,'(a25)', advance='NO') 'conflict-free transpose'
     !call transposeNoBankConflicts<<<dimGrid, dimBlock>>>(Arr, transp_Arr, nx, ny, TILE_DIM, BLOCK_ROWS)
     call transposeNaive<<<dimGrid, dimBlock>>>(Arr, transp_Arr, nx, ny, TILE_DIM, BLOCK_ROWS)
     !call transposeCoalesced<<<dimGrid, dimBlock>>>(Arr, transp_Arr, nx, ny, TILE_DIM, BLOCK_ROWS)
     !print*, "[+] shape(transp_Arr) = ", shape(transp_Arr)
     !transp_Hrr = transp_Arr
     !print*, "[+] transp_Hrr(100,200) = ",  transp_Hrr(100,200)
     !print*, "[+] transp_Hrr(200,100) = ",  transp_Hrr(200,100)

  end subroutine transpose_device_array


  !> Device array copy X->Y: Y(i1:i2,j1:j2)= X(p1:p2,q1:q2)
  subroutine copy_arr2d(X,Y, i1,i2,j1,j2, p1,q1, block_threads_)
  implicit none
  double precision, device,    intent(IN) :: X(:,:)
  double precision, device, intent(INOUT) :: Y(:,:)
  integer,                     intent(IN) :: i1,i2,j1,j2,p1,q1
  integer,          optional,  intent(IN) :: block_threads_
  !
  integer    :: block_threads, shA(2), shB(2)
  type(dim3) :: grid, tBlock
  character(*), parameter :: subnam = "[copy_arr2d]"

     if (i1 > i2 .or. j1 > j2) return
     block_threads = 256
     if(present(block_threads_)) block_threads = block_threads_
     shA = shape(X); shB = shape(Y)
     if (i1 < 1 .or. i1 > shA(1)) error stop '[!]'//subnam//' i1 out of range'
     if (j1 < 1 .or. j1 > shA(2)) error stop '[!]'//subnam//' j1 out of range'
     if (i2 < 1 .or. i2 > shA(1)) error stop '[!]'//subnam//' i2 out of range'
     if (j2 < 1 .or. j2 > shA(2)) error stop '[!]'//subnam//' j2 out of range'
     if (p1 < 1 .or. p1 > shB(1)) error stop '[!]'//subnam//' p1 out of range'
     if (q1 < 1 .or. q1 > shB(2)) error stop '[!]'//subnam//' q1 out of range'
     if (p1+i2-i1 > shB(1)) error stop '[!]'//subnam//': p1+i2-i1 > shB(1)'
     if (q1+j2-j1 > shB(2)) error stop '[!]'//subnam//': q1+j2-j1 > shB(2)'
     tBlock = dim3(1,block_threads,1)
     grid = dim3(1, (j2 - j1 - 1)/tBlock%y + 1, 1)
     call cu_copy_arr2d<<<grid,tBlock>>>(X,Y,i1,i2,j1,j2,p1,q1)

  end subroutine copy_arr2d


  !> The a*X + Y -> Y operation: Y(i1:i2,j1:j2)+= alpha*X(p1:p2,q1:q2)
  subroutine axpy_arr2d(alpha,X,Y, i1,i2,j1,j2, p1,q1, block_threads_)
  implicit none
  double precision,            intent(IN) :: alpha
  double precision, device,    intent(IN) :: X(:,:)
  double precision, device, intent(INOUT) :: Y(:,:)
  integer,                     intent(IN) :: i1,i2,j1,j2,p1,q1
  integer,          optional,  intent(IN) :: block_threads_
  !
  integer    :: block_threads, shA(2), shB(2)
  type(dim3) :: grid, tBlock
  character(*), parameter :: subnam = "[axpy_arr3d]"

     if (i1 > i2 .or. j1 > j2) return
     block_threads = 256
     if(present(block_threads_)) block_threads = block_threads_
     shA = shape(X); shB = shape(Y)
     if (i1 < 1 .or. i1 > shA(1)) error stop '[!]'//subnam//' i1 out of range'
     if (j1 < 1 .or. j1 > shA(2)) error stop '[!]'//subnam//' j1 out of range'
     if (i2 < 1 .or. i2 > shA(1)) error stop '[!]'//subnam//' i2 out of range'
     if (j2 < 1 .or. j2 > shA(2)) error stop '[!]'//subnam//' j2 out of range'
     if (p1 < 1 .or. p1 > shB(1)) error stop '[!]'//subnam//' p1 out of range'
     if (q1 < 1 .or. q1 > shB(2)) error stop '[!]'//subnam//' q1 out of range'
     if (p1+i2-i1 > shB(1)) error stop '[!]'//subnam//': p1+i2-i1 > shB(1)'
     if (q1+j2-j1 > shB(2)) error stop '[!]'//subnam//': q1+j2-j1 > shB(2)'
     tBlock = dim3(1,block_threads,1)
     grid = dim3(1, (j2 - j1 - 1)/tBlock%y + 1, 1)
     call cu_axpy_arr2d<<<grid,tBlock>>>(alpha,X,Y,i1,i2,j1,j2,p1,q1)

  end subroutine axpy_arr2d


  !> Device array copy X->Y: Y(i1:i2,j1:j2,k1:k2)= X(p1:p2,q1:q2,r1:r2)
  subroutine copy_arr3d(X,Y, i1,i2,j1,j2,k1,k2, p1,q1,r1, block_threads_)
  implicit none
  double precision, device,    intent(IN) :: X(:,:,:)
  double precision, device, intent(INOUT) :: Y(:,:,:)
  integer,                     intent(IN) :: i1,i2,j1,j2,k1,k2,p1,q1,r1
  integer,          optional,  intent(IN) :: block_threads_
  !
  integer    :: block_threads, shA(3), shB(3)
  type(dim3) :: grid, tBlock
  character(*), parameter :: subnam = "[copy_arr3d]"

     block_threads = 16
     if(present(block_threads_)) block_threads = block_threads_
     shA = shape(X); shB = shape(Y)
     if (i1 < 1 .or. i1 > shA(1)) error stop '[!]'//subnam//'i1 out of range'
     if (j1 < 1 .or. j1 > shA(2)) error stop '[!]'//subnam//'j1 out of range'
     if (k1 < 1 .or. k1 > shA(3)) error stop '[!]'//subnam//'k1 out of range'
     if (i2 < 1 .or. i2 > shA(1)) error stop '[!]'//subnam//'i2 out of range'
     if (j2 < 1 .or. j2 > shA(2)) error stop '[!]'//subnam//'j2 out of range'
     if (k2 < 1 .or. k2 > shA(3)) error stop '[!]'//subnam//'k2 out of range'
     if (p1 < 1 .or. p1 > shB(1)) error stop '[!]'//subnam//'p1 out of range'
     if (q1 < 1 .or. q1 > shB(2)) error stop '[!]'//subnam//'q1 out of range'
     if (r1 < 1 .or. r1 > shB(3)) error stop '[!]'//subnam//'r1 out of range'
     if (p1+i2-i1 > shB(1)) error stop '[!]'//subnam//': p1+i2-i1 > shB(1)'
     if (q1+j2-j1 > shB(2)) error stop '[!]'//subnam//': q1+j2-j1 > shB(2)'
     if (r1+k2-k1 > shB(3)) error stop '[!]'//subnam//': r1+k2-k1 > shB(3)'
     if (i1 > i2 .or. j1 > j2 .or. k1 > k2) return
     tBlock = dim3(1,block_threads,block_threads)
     grid = dim3(1, (j2 - j1 - 1)/tBlock%y + 1, (k2 - k1 - 1)/tBlock%z + 1)
     call cu_copy_arr3d<<<grid,tBlock>>>(X,Y,i1,i2,j1,j2,k1,k2,p1,q1,r1)

  end subroutine copy_arr3d


  !> The a*X + Y -> Y operation: Y(i1:i2,j1:j2,k1:k2)+= alpha*X(p1:p2,q1:q2,r1:r2)
  subroutine axpy_arr3d(alpha,X,Y, i1,i2,j1,j2,k1,k2, p1,q1,r1, block_threads_)
  implicit none
  double precision,            intent(IN) :: alpha
  double precision, device,    intent(IN) :: X(:,:,:)
  double precision, device, intent(INOUT) :: Y(:,:,:)
  integer,                     intent(IN) :: i1,i2,j1,j2,k1,k2,p1,q1,r1
  integer,          optional,  intent(IN) :: block_threads_
  !
  integer    :: block_threads, shA(3), shB(3)
  type(dim3) :: grid, tBlock
  character(*), parameter :: subnam = "[axpy_arr3d]"

     block_threads = 16
     if(present(block_threads_)) block_threads = block_threads_
     shA = shape(X); shB = shape(Y)
     if (i1 < 1 .or. i1 > shA(1)) error stop '[!]'//subnam//'i1 out of range'
     if (j1 < 1 .or. j1 > shA(2)) error stop '[!]'//subnam//'j1 out of range'
     if (k1 < 1 .or. k1 > shA(3)) error stop '[!]'//subnam//'k1 out of range'
     if (i2 < 1 .or. i2 > shA(1)) error stop '[!]'//subnam//'i2 out of range'
     if (j2 < 1 .or. j2 > shA(2)) error stop '[!]'//subnam//'j2 out of range'
     if (k2 < 1 .or. k2 > shA(3)) error stop '[!]'//subnam//'k2 out of range'
     if (p1 < 1 .or. p1 > shB(1)) error stop '[!]'//subnam//'p1 out of range'
     if (q1 < 1 .or. q1 > shB(2)) error stop '[!]'//subnam//'q1 out of range'
     if (r1 < 1 .or. r1 > shB(3)) error stop '[!]'//subnam//'r1 out of range'
     if (p1+i2-i1 > shB(1)) error stop '[!]'//subnam//': p1+i2-i1 > shB(1)'
     if (q1+j2-j1 > shB(2)) error stop '[!]'//subnam//': q1+j2-j1 > shB(2)'
     if (r1+k2-k1 > shB(3)) error stop '[!]'//subnam//': r1+k2-k1 > shB(3)'
     if (i1 > i2 .or. j1 > j2 .or. k1 > k2) return
     tBlock = dim3(1,block_threads,block_threads)
     grid = dim3(1, (j2 - j1 - 1)/tBlock%y + 1, (k2 - k1 - 1)/tBlock%z + 1)
     call cu_axpy_arr3d<<<grid,tBlock>>>(alpha,X,Y,i1,i2,j1,j2,k1,k2,p1,q1,r1)

  end subroutine axpy_arr3d


  !> zero out elements below matrix diagonal
  subroutine zero_below_diag(locA, M,N, MB,NB, iB,jB, ML,NL, nRowDevs, nColDevs)
  implicit none
  double precision, intent(inout), device :: locA(:, :)
  integer, intent(in) :: M, N, MB, NB, iB, jB, ML, NL, nRowDevs, nColDevs
  !
  type(dim3) :: grid, thr_block

     if (ML*NL.eq.0) return

     thr_block = dim3(32, 32, 1)
     grid      = dim3((ML - 1)/thr_block% x + 1, &
                      (NL - 1)/thr_block% y + 1, 1)
     call cu_zero_below_diag<<<grid, thr_block>>>(locA, M,N, MB,NB, ML,NL,&
                                                  iB,jB, nRowDevs,nColDevs)

  end subroutine zero_below_diag


  !> zero out elements above matrix diagonal
  subroutine zero_above_diag(locA, M,N, MB,NB, iB,jB, ML,NL, nRowDevs, nColDevs)
  implicit none
  double precision, intent(inout), device :: locA(:, :)
  integer, intent(in) :: M, N, MB, NB, iB, jB, ML, NL, nRowDevs, nColDevs
  !
  type(dim3) :: grid, thr_block

     if (ML*NL.eq.0) return

     thr_block = dim3(32, 32, 1)
     grid      = dim3((ML - 1)/thr_block% x + 1, &
                      (NL - 1)/thr_block% y + 1, 1)
     call cu_zero_above_diag<<<grid, thr_block>>>(locA, M,N, MB,NB, ML,NL,&
                                                  iB,jB, nRowDevs,nColDevs)

  end subroutine zero_above_diag


  !> initialize a distributed identity matrix
  subroutine create_eye_matrix(locA, M,N, MB,NB, iB,jB, ML,NL, nRowDevs, nColDevs)
  implicit none
  double precision, pointer, device, intent(INOUT) :: locA(:)
  integer, intent(IN) :: M, N, MB, NB, iB, jB, ML, NL, nRowDevs, nColDevs
  !
  type(dim3) :: grid, thr_block

     allocate(locA(ML*NL))

     thr_block = dim3(32, 32, 1)
     grid      = dim3((ML - 1)/thr_block% x + 1, &
                      (NL - 1)/thr_block% y + 1, 1)
     call cu_eye<<<grid, thr_block>>>(locA, M,N, MB,NB, ML,NL,&
                                      iB,jB, nRowDevs,nColDevs)
  end subroutine create_eye_matrix


  !> add an array to the n-th sub/super diagonal of a distributed matrix
  subroutine add_matrix_diag(locA, diag, ndiag, M,N,MB,NB,ML,NL,iB,jB,nRowDevs,nColDevs)
  implicit none
  double precision, pointer, device, intent(INOUT) :: locA(:) !< size ML*NL (local matrix)
  double precision, pointer, device,    intent(IN) :: diag(:) !< size min(M,L) - ndiag
  integer, intent(IN) :: ndiag !< which diagonal: 0=main, +1=1st super, -1=1st subdiagonal etc.
  integer, intent(IN) :: M, N  !< global matrix dimension (rows, columns)
  integer, intent(IN) :: MB,NB !< blocking dimensions in rows, columns
  integer, intent(IN) :: ML,NL !< local matrix: ML x NL 
  integer, intent(IN) :: iB,jB !< index of the local matrix in the global device grid
  integer, intent(IN) :: nRowDevs,nColDevs !< global device grid sizes
  !
  type(dim3) :: grid, thr_block

     thr_block = dim3(32, 32, 1)
     grid      = dim3((ML - 1)/thr_block% x + 1, &
                      (NL - 1)/thr_block% y + 1, 1)
     call cu_add_matrix_diag<<<grid, thr_block>>>(locA, diag, ndiag, &
                                                  M,N, MB,NB, ML,NL, &
                                                  iB,jB, nRowDevs,nColDevs)
  end subroutine add_matrix_diag


  !> shift local 2D array to the left by k positions: A(i,:) = A(i+k,:)
  subroutine shift_left_arr2d(A, M, N, k, block_threads_)
  implicit none
  double precision, pointer, device, intent(INOUT) :: A(:,:)
  integer, intent(IN) :: M, N, k
  integer, intent(IN), optional           :: block_threads_
  !
  character(*), parameter :: subnam = "[shift_left_arr2d]"
  integer    :: block_threads
  type(dim3) :: grid, thr_block
     
     block_threads = 256
     if(present(block_threads_)) block_threads = block_threads_
     if(k <= 0) error stop '[!]'//subnam//': k must be positive'
     thr_block = dim3(block_threads, 1, 1)
     grid      = dim3((M - 1)/thr_block% x + 1, 1, 1)
     call cu_shift_left_arr2d<<<grid, thr_block>>>(A, M, N, k)

  end subroutine shift_left_arr2d

  !> shift local 2D array to the right by k positions: A(i,:) = A(i-k,:)
  subroutine shift_right_arr2d(A, M, N, k, block_threads_)
  implicit none
  double precision, pointer, device, intent(INOUT) :: A(:,:)
  integer, intent(IN) :: M, N, k
  integer, intent(IN), optional           :: block_threads_
  !
  character(*), parameter :: subnam = "[shift_right_arr2d]"
  integer    :: block_threads
  type(dim3) :: grid, thr_block
     
     block_threads = 256
     if(present(block_threads_)) block_threads = block_threads_
     if(k <= 0) error stop '[!]'//subnam//': k must be positive'
     thr_block = dim3(block_threads, 1, 1)
     grid      = dim3((M - 1)/thr_block% x + 1, 1, 1)
     call cu_shift_right_arr2d<<<grid, thr_block>>>(A, M, N, k)

  end subroutine shift_right_arr2d


  !> shift local 2D array upward by k positions: A(i,:) = A(i+k,:)
  subroutine shift_up_arr2d(A, M, N, k, block_threads_)
  implicit none
  double precision, pointer, device, intent(INOUT) :: A(:,:)
  integer, intent(IN) :: M, N, k
  integer, intent(IN), optional           :: block_threads_
  !
  character(*), parameter :: subnam = "[shift_up_arr2d]"
  integer    :: block_threads
  type(dim3) :: grid, thr_block
     
     block_threads = 256
     if(present(block_threads_)) block_threads = block_threads_
     if(k <= 0) error stop '[!]'//subnam//': k must be positive'
     thr_block = dim3(1, block_threads, 1)
     grid      = dim3(1, (N - 1)/thr_block% x + 1, 1)
     call cu_shift_up_arr2d<<<grid, thr_block>>>(A, M, N, k)

  end subroutine shift_up_arr2d


  !> shift local 2D array downward by k positions: A(i,:) = A(i+k,:)
  subroutine shift_down_arr2d(A, M, N, k, block_threads_)
  implicit none
  double precision, pointer, device, intent(INOUT) :: A(:,:)
  integer, intent(IN) :: M, N, k
  integer, intent(IN), optional           :: block_threads_
  !
  character(*), parameter :: subnam = "[shift_down_arr2d]"
  integer    :: block_threads
  type(dim3) :: grid, thr_block
     
     block_threads = 256
     if(present(block_threads_)) block_threads = block_threads_
     if(k <= 0) error stop '[!]'//subnam//': k must be positive'
     thr_block = dim3(1, block_threads, 1)
     grid      = dim3(1, (N - 1)/thr_block% x + 1, 1)
     call cu_shift_down_arr2d<<<grid, thr_block>>>(A, M, N, k)

  end subroutine shift_down_arr2d

end module matrix_ops
