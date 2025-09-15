!> file: poisson2d_dg.f90
!!
!! Solve Poisson equation in 2D using finite elements
!!
program poisson_2d_dg
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: Nx = 10, Ny = 10  ! Number of elements in x and y direction
  real(dp), parameter :: Lx = 1.0_dp, Ly = 1.0_dp
  real(dp) :: hx, hy
  real(dp), allocatable :: u(:,:), f(:,:), u_exact(:,:)
  integer :: i, j

  ! Compute element sizes
  hx = Lx / Nx
  hy = Ly / Ny

  ! Allocate arrays
  allocate(u(Nx+1, Ny+1), f(Nx+1, Ny+1), u_exact(Nx+1, Ny+1))

  ! Initialize source term and exact solution (for testing)
  call initialize(f, u_exact, hx, hy)

  ! Solve the Poisson equation using DG method
  call solve_poisson_dg(Nx, Ny, hx, hy, f, u)

  ! Output the results
  print *, "Solution u:"
  do j = 0, Ny
    do i = 0, Nx
      print *, i*hx, j*hy, u(i,j)
    end do
  end do

contains

  subroutine initialize(f, u_exact, hx, hy)
    real(dp), intent(out) :: f(:,:), u_exact(:,:)
    real(dp), intent(in) :: hx, hy
    integer :: i, j

    ! Example source term and exact solution (for testing)
    do j = 0, size(f,2)-1
      do i = 0, size(f,1)-1
        f(i,j) = sin(pi*i*hx) * sin(pi*j*hy)
        u_exact(i,j) = sin(pi*i*hx) * sin(pi*j*hy) / (2.0_dp*pi*pi)
      end do
    end do
  end subroutine initialize

  subroutine solve_poisson_dg(Nx, Ny, hx, hy, f, u)
    integer, intent(in) :: Nx, Ny
    real(dp), intent(in) :: hx, hy
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(out) :: u(:,:)
    real(dp), allocatable :: A(:,:), b(:), u_flat(:)
    integer :: i, j, k, l, index, n, info

    ! Number of degrees of freedom
    n = (Nx+1) * (Ny+1)
    allocate(A(n,n), b(n), u_flat(n))

    ! Initialize matrices
    A = 0.0_dp
    b = 0.0_dp
    u_flat = 0.0_dp

    ! Assemble the system (simple example, not optimized)
    do j = 0, Ny
      do i = 0, Nx
        index = i + j*(Nx+1)
        ! Apply source term
        b(index) = hx*hy * f(i,j)
        ! Apply boundary conditions
        if (i == 0 .or. i == Nx .or. j == 0 .or. j == Ny) then
          A(index, index) = 1.0_dp
          b(index) = 0.0_dp
        else
          A(index, index) = -2.0_dp*(1.0_dp/hx**2 + 1.0_dp/hy**2)
          if (i > 0) A(index, index-1) = 1.0_dp/hx**2
          if (i < Nx) A(index, index+1) = 1.0_dp/hx**2
          if (j > 0) A(index, index-(Nx+1)) = 1.0_dp/hy**2
          if (j < Ny) A(index, index+(Nx+1)) = 1.0_dp/hy**2
        end if
      end do
    end do

    ! Solve the linear system A*u_flat = b using LAPACK
    call solve_linear_system(A, b, u_flat)

    ! Reshape the solution back to 2D
    do j = 0, Ny
      do i = 0, Nx
        index = i + j*(Nx+1)
        u(i,j) = u_flat(index)
      end do
    end do

    ! Deallocate temporary arrays
    deallocate(A, b, u_flat)
  end subroutine solve_poisson_dg

  subroutine solve_linear_system(A, b, x)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    integer :: n, info, ipiv(size(b))

    ! Copy b into x to use as the right-hand side
    x = b
    n = size(b)

    ! Solve the linear system using LAPACK
    call dgesv(n, 1, A, n, ipiv, x, n, info)
    if (info /= 0) then
      print *, "Error: LAPACK solver failed with info =", info
      stop
    end if
  end subroutine solve_linear_system

end program poisson_2d_dg

