!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_mat_ops_kernels.f90
!! @author Ismael Boureima, Oleg Korobkin
!! @date  March 2025
!! @brief Testing CUDA kernels
!!
program main
use rnd_lib, only: random
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: m = 10, n = 8  ! in matrix transpose test: MxN
integer :: nx = 12, ny = 3, nz = 4 ! in copy3d test: nx,ny,nz
integer :: k = 2 ! shift step
double precision, allocatable :: rnd(:)
integer :: what_test = 1

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
     elseif (trim(arg).eq.'--test-transpose') then
        what_test = 1
     elseif (trim(arg).eq.'--test-copy-arr2d') then
        what_test = 2
     elseif (trim(arg).eq.'--test-copy-arr3d') then
        what_test = 3
     elseif (trim(arg).eq.'--test-shift-arr2d') then
        what_test = 4
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  select case (what_test)
  case (1)
     call test_transpose(m, n)
  case (2)
     call test_copy_arr2d(m, n)
  case (3)
     call test_copy_arr3d(nx,ny,nz)
  case (4)
     call test_shift_arr2d(m, n, k)
  endselect

contains

  subroutine help_msg()
  print '(13(/,A))', &
   "Test CUDA BLAS subroutines", &
   "Usage: ./test.exe [-h] [-m <M>] [-n <N>]i [k <k>]", &
   "        [--test-transpose |--test-copy-arr2d|",  &
   "        |--test-copy-arr3d|--test-shift-arr2d]", &
   "Options:", &
   " -h           : print this help message", &
   " -m <M>       : in Level-1 tests: length of a vector A [10]", &
   " -n <N>       : number of columns in the matrix A [10]"
  end subroutine


  subroutine test_transpose(M, N)
  use matrix_ops
  use matrix_util, only: pprint_matrix
  implicit none
  integer, intent(IN) :: M, N !< size of the matrix to transpose
  !
  integer, parameter  :: TILE_DIM=32
  integer, parameter  :: BLOCK_ROWS=8
  double precision, allocatable :: h_A(:,:),h_A_T(:,:)
  double precision, allocatable, device, target :: d_A(:,:), d_A_T(:,:)
  double precision, pointer, device :: d_AV(:), d_AV_T(:)

     allocate(h_A(M,N), h_A_T(N,M))

     call random_number(h_A)
     print '(/"Testing transpose_device_array_vect:")'
     print '("[+] h_A(1,N) = ",ES14.7)', h_A(1,N)
     print '("[+] h_A(M,1) = ",ES14.7)', h_A(M,1)
     call pprint_matrix(h_A)

     allocate(d_A(M,N), d_A_T(N,M))
     d_A   = h_A  ! H2D
     d_A_T = 0.d0 !*h_A
     d_AV(1:M*N) => d_A
     d_AV_T(1:M*N) => d_A_T

     !!d_A_T = transpose_device_array_naive(d_A, TILE_DIM, BLOCK_ROWS)
     !!call transpose_device_array(d_A, d_A_T, TILE_DIM, BLOCK_ROWS)
     !call transpose_device_array_naive(d_A, d_A_T, TILE_DIM, BLOCK_ROWS)
     call transpose_device_array_vect(d_AV, d_AV_T, M, N, 256)
     !!print*, "transpose computed OK"
     h_A_T = d_A_T ! D2H
     print '(/"After transpose:")'

     print '("[+] h_A_T(1,M) = ",ES14.7)', h_A_T(1,M)
     print '("[+] h_A_T(N,1) = ",ES14.7)', h_A_T(N,1)
     call pprint_matrix(h_A_T)

     deallocate(h_A,h_A_T,d_A,d_A_T)

  end subroutine test_transpose


  subroutine test_copy_arr2d(nx, ny)
  use matrix_ops
  use matrix_util, only: pprint_matrix
  implicit none
  integer, intent(IN) :: nx, ny !< size of the array
  !
  character(*), parameter :: subnam = '[test_copy_arr2d]'
  integer, parameter  :: TILE_DIM=32
  integer, parameter  :: BLOCK_ROWS=8
  double precision, allocatable :: h_A2(:,:),h_B2(:,:)
  double precision, allocatable, device :: d_A2(:,:), d_B2(:,:)

     print '(/"Testing copy_arr3d, {nx,ny} = {",2(I4),"}:")',nx,ny
     if (nx<6 .or. ny<3) error stop "[!]"//subnam//": size too small"
     allocate(h_A2(nx,ny))
     h_A2 = 0d0; h_A2(2:4,2:3) = -1d0
     print '("Array d_A2 before copy (h_A2(2:4,2:3) -> h_B2(4:6,1:2)):")'
     call pprint_matrix(h_A2,frmt_='(F3.0)')

     d_A2 = h_A2
     allocate(d_B2(nx,ny))
     d_B2 = 2d0
     call copy_arr2d(d_A2,d_B2, 2,4, 2,3, 4,1)
     h_B2 = d_B2
     print '(/"Array d_B2 after copy:")'
     call pprint_matrix(h_B2,frmt_='(F3.0)')

  end subroutine test_copy_arr2d


  subroutine test_copy_arr3d(nx, ny, nz)
  use matrix_ops
  use matrix_util, only: pprint_matrix, pprint_matrix3d
  implicit none
  integer, intent(IN) :: nx, ny, nz !< size of the array
  !
  character(*), parameter :: subnam = '[test_copy_arr3d]'
  integer, parameter  :: TILE_DIM=32
  integer, parameter  :: BLOCK_ROWS=8
  double precision, allocatable :: h_A3(:,:,:),h_B3(:,:,:)
  double precision, allocatable, device :: d_A3(:,:,:), d_B3(:,:,:)

     print '(/"Testing copy_arr3d, {nx,ny,nz} = {",3(I4),"}:")',nx,ny,nz
     if (nx<6 .or. ny<3 .or. nz<3) error stop "[!]"//subnam//": size too small"
     allocate(h_A3(nx,ny,nz))
     h_A3 = 0d0; h_A3(2:4,2:3,3) = -1d0
     print '("Array d_A3 before copy (h_A3(2:4,2:3,3) -> h_B3(4:6,1:2,2)):")'
     call pprint_matrix3d(h_A3,frmt_='(F3.0)')

     d_A3 = h_A3
     allocate(d_B3(nx,ny,nz))
     d_B3 = 2d0
     call copy_arr3d(d_A3,d_B3, 2,4, 2,3, 3,3, 4,1,2)
     h_B3 = d_B3
     print '(/"Array d_B3 after copy:")'
     call pprint_matrix3d(h_B3,frmt_='(F3.0)')

  end subroutine test_copy_arr3d


  subroutine test_shift_arr2d(m, n, k)
  use matrix_ops
  use matrix_util, only: pprint_matrix
  implicit none
  integer, intent(IN) :: m, n !< size of the array
  integer, intent(IN) :: k    !< shift step
  !
  character(*), parameter :: subnam = '[test_shift_arr2d]'
  double precision, allocatable :: h_A2(:,:)
  double precision, pointer, device :: d_A2(:,:)
  integer :: i, j, sh(2)

     print '(/"Testing shift_arr3d, {m,n} = {",2(I4),"}:")',m,n
     allocate(h_A2(m,n))
     do i=1,m
        h_A2(i,:) = 100*(i - 1)
     enddo
     do j=1,n
        h_A2(:,j) = h_A2(:,j) + j - 1
     enddo
     print '("Array d_A2 before shift:")'
     call pprint_matrix(h_A2,frmt_='(F6.0)')

     allocate(d_A2(m,n))
     d_A2 = h_A2
     call shift_left_arr2d(d_A2, m, n, k)
     print '("Array d_A2 after shift to the left:")'
     h_A2 = d_A2
     call pprint_matrix(h_A2,frmt_='(F6.0)')

     call shift_right_arr2d(d_A2, m, n, k)
     print '("Array d_A2 after shift to the right:")'
     h_A2 = d_A2
     call pprint_matrix(h_A2,frmt_='(F6.0)')

     call shift_up_arr2d(d_A2, m, n, k)
     print '("Array d_A2 after downward shift:")'
     h_A2 = d_A2
     call pprint_matrix(h_A2,frmt_='(F6.0)')

     call shift_down_arr2d(d_A2, m, n, k)
     print '("Array d_A2 after upward shift:")'
     h_A2 = d_A2
     call pprint_matrix(h_A2,frmt_='(F6.0)')

     deallocate(d_A2)

  end subroutine test_shift_arr2d


end program main
