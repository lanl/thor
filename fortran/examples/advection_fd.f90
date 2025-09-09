!>
!! Solves advection equation using finite differences
!! 1: using a full-grid solver
!! 2: using tensor-trains approach
!!
program advection_fd
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer :: Nx, num_steps
  character(len=128) arg
  real(dp) :: dt, dx, c, courant
  real(dp), allocatable :: u(:)
  integer :: i, ierr

  ! Read Nx and num_steps from command line
  if (command_argument_count() .ne. 2) then
    print *, "Usage: advection_fd <number_of_points> <number_of_time_steps>"
    stop
  end if

  Nx = 0; num_steps = 0
  call getarg(1, arg)
  read(arg,*,iostat=ierr) Nx
  call getarg(2, arg)
  read(arg,*,iostat=ierr) num_steps
  
  ! Define advection speed and time step
  c = 1.0d0
  dx = 1.0d0 / real(Nx, dp)
  courant = 0.75
  dt = dx/c * courant

  ! Call the solver function
  call advection_solver(Nx, num_steps, dt, c, u)

  ! Output the results
  print *, "# x, u"
  do i = 1, Nx
    print "(ES12.5,1X,ES14.7)", (i - 1) / real(Nx, dp), u(i)
  end do

contains

  subroutine advection_solver(Nx, num_steps, dt, c, u)
    integer, intent(in) :: Nx, num_steps
    real(dp), intent(in) :: dt, c
    real(dp), allocatable, intent(out) :: u(:)
    real(dp), allocatable :: u_new(:), x(:)
    real(dp), parameter :: pi = 3.141592653589793d0
    real(dp) :: dx, courant
    integer :: i, n

    dx = 1.0d0 / Nx
    courant = c * dt / dx
    if (courant > 1.0d0) then
      print *, "Courant number too large! Stability may be compromised."
      stop
    end if

    allocate(x(Nx), u(Nx), u_new(Nx))

    ! Initialize x and u with a Gaussian pulse as the initial condition
    do i = 1, Nx
      x(i) = (i - 1) * dx
      u(i) = exp(-100.0d0 * (x(i) - 0.5d0)**2)
    end do

    ! Time-stepping loop
    do n = 1, num_steps
      ! Update the solution using the explicit finite difference method
      do i = 2, Nx-1
        u_new(i) = u(i) - courant * (u(i) - u(i-1))
      end do
      ! Periodic boundary conditions
      u_new(1) = u(1) - courant * (u(1) - u(Nx))
      u_new(Nx) = u(Nx) - courant * (u(Nx) - u(Nx-1))

      ! Swap pointers
      u = u_new
    end do

  end subroutine advection_solver

end program advection_fd

