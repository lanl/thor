!> file: compress_dataset.f90
!!
!! Reads a range of binary *.npy files in a directory
!! Compresses them into a TT format
!!
program process_dataset
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use shmor_mod, only: shmor_svd, read_npy  ! Assuming the library is called "shmor_mod"
  implicit none

  character(len=*), parameter :: npy_ext = ".npy"
  character(len=:), allocatable :: db_location
  character(len=100) :: filename
  integer :: num_files, ierr, i, j, k, l
  integer :: nx, ny, nz  ! Dimensions of the 3D arrays
  logical :: done
  real(dp), allocatable :: data(:,:,:), DS(:,:,:,:)

  ! Read command-line argument for database location
  call get_command_argument(1, db_location, ierr)
  if (ierr /= 0 .or. len_trim(db_location) == 0) then
    print *, "Usage: process_dataset <database_location>"
    stop
  end if

  ! Get the number of .npy files and their dimensions
  call count_and_get_dimensions(db_location, num_files, nx, ny, nz)

  ! Allocate the 4D array DS
  allocate(DS(nx, ny, nz, num_files))

  ! Read each file and copy the data into the 4D array DS
  do i = 0, num_files-1
    write(filename, '(A,I4.4,A)') 'data', i, npy_ext
    call read_npy(trim(db_location) // '/' // trim(filename), data)
    DS(:,:,:,i+1) = data
  end do

  ! Apply the Shmor SVD function
  call shmor_svd(DS)

contains

  subroutine count_and_get_dimensions(location, num_files, nx, ny, nz)
    character(len=*), intent(in) :: location
    integer, intent(out) :: num_files, nx, ny, nz
    character(len=100) :: command, line
    character(len=100), allocatable :: filenames(:)
    real(dp), allocatable :: temp(:,:,:)
    integer :: ios, unit, count

    ! Initialize
    num_files = 0
    nx = 0
    ny = 0
    nz = 0

    ! Get the list of .npy files in the directory
    open(unit=unit, file="list_of_files.txt", status="replace")
    write(command, '(A)') "ls -1 " // trim(location) // "/data*.npy > list_of_files.txt"
    call execute_command_line(trim(command))

    ! Count the files
    open(unit=unit, file="list_of_files.txt", status="old", action="read")
    count = 0
    do while (.true.)
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      count = count + 1
    end do
    close(unit)

    num_files = count
    if (num_files == 0) then
      print *, "No .npy files found in the specified directory."
      stop
    end if

    ! Get the dimensions of the first file
    open(unit=unit, file="list_of_files.txt", status="old", action="read")
    read(unit, '(A)', iostat=ios) line
    close(unit)
    call read_npy(trim(line), temp)
    nx = size(temp, 1)
    ny = size(temp, 2)
    nz = size(temp, 3)
    deallocate(temp)
  end subroutine count_and_get_dimensions

end program process_dataset

