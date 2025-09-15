!> file: poisson_fd.f90
!!
!! Solves Poisson equation with full grid and tensor-train methods
!!
program poisson_fd
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer :: N
  real(dp), allocatable :: u(:)
  integer :: ierr

  ! Read N from command line
  call get_command_argument(1, N, ierr)
  if (ierr /= 0 .or. N <= 0) then
    print *, "Usage: poisson_fd <number_of_points>"
    stop
  end if

  ! Call the solver function
  call full_grid_solver(N, u)

  ! Call the tensor-train solver
  call tt_solver(N, u)

  ! Output the results
  print *, "x, u"
  do i = 1, N+2
    print *, (i - 1) / real(N + 1, dp), u(i)
  end do

contains

  subroutine full_grid_solver(N, u)
    integer, intent(in) :: N
    real(dp), allocatable, intent(out) :: u(:)
    real(dp), allocatable :: x(:), f(:)
    real(dp), parameter :: pi = 3.141592653589793d0
    real(dp), parameter :: h = 1.0d0 / (N + 1)
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: b
    integer :: i

    allocate(x(N+2), u(N+2), f(N+2))
    allocate(A(N,N), b(N))

    ! Initialize x, u, and f
    do i = 1, N+2
      x(i) = (i - 1) * h
      u(i) = 0.0d0
      if (i == 1 .or. i == N+2) then
        f(i) = 0.0d0
      else
        f(i) = sin(pi * x(i))  ! Example source term
      end if
    end do

    ! Set up the matrix A and right-hand side b
    A = 0.0d0
    b = 0.0d0
    do i = 1, N
      if (i > 1) A(i,i-1) = 1.0d0 / h**2
      A(i,i) = -2.0d0 / h**2
      if (i < N) A(i,i+1) = 1.0d0 / h**2
      b(i) = f(i+1)
    end do

    ! Solve the linear system A * u = b using Gaussian elimination
    call gauss_elimination(A, b, u(2:N+1))

  end subroutine full_grid_solver

  subroutine gauss_elimination(A, b, u)
    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(inout) :: b(:)
    real(dp), intent(out) :: u(:)
    integer :: n, i, j, k
    real(dp) :: factor

    n = size(b)
    u = 0.0d0

    ! Forward elimination
    do k = 1, n-1
      do i = k+1, n
        factor = A(i,k) / A(k,k)
        A(i,k) = 0.0d0
        A(i,k+1:n) = A(i,k+1:n) - factor * A(k,k+1:n)
        b(i) = b(i) - factor * b(k)
      end do
    end do

    ! Back substitution
    u(n) = b(n) / A(n,n)
    do i = n-1, 1, -1
      u(i) = (b(i) - dot_product(A(i,i+1:n), u(i+1:n))) / A(i,i)
    end do

  end subroutine gauss_elimination

end program poisson_fd

