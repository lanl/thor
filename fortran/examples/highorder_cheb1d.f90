!> file: highorder_cheb1d.f0-
!!
!! Solve an arbitrary nth-order linear ODE 
!!
program chebyshev_spectral
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: N = 20  ! Number of Chebyshev nodes
  integer :: n, i, j
  real(dp) :: a(0:N), b
  real(dp), allocatable :: D(:,:), x(:), u(:)
  real(dp), allocatable :: A(:,:), rhs(:)

  ! Define the coefficients of the differential equation
  ! Example: 2nd order equation: u'' - 2u' + u = sin(x)
  n = 2
  a(2) = 1.0_dp  ! Coefficient of u''
  a(1) = -2.0_dp ! Coefficient of u'
  a(0) = 1.0_dp  ! Coefficient of u
  b = 0.0_dp     ! Source term

  ! Initialize arrays
  allocate(D(N,N), x(N), u(N))
  allocate(A(N,N), rhs(N))

  ! Compute Chebyshev nodes
  call chebyshev_nodes(N, x)

  ! Compute Chebyshev differentiation matrix
  call chebyshev_differentiation_matrix(N, D)

  ! Formulate the linear system
  A = 0.0_dp
  rhs = 0.0_dp
  do i = 1, N
    rhs(i) = b
    do j = 1, N
      A(i,j) = a(n) * D(i,j)
      if (n > 1) then
        A(i,j) = A(i,j) + a(n-1) * D(i,j)
      endif
      if (n > 0) then
        A(i,j) = A(i,j) + a(n-2) * D(i,j)
      endif
      A(i,j) = A(i,j) + a(0)
    end do
  end do

  ! Solve the linear system
  call solve_linear_system(A, rhs, u)

  ! Output the results
  print *, "x, u"
  do i = 1, N
    print *, x(i), u(i)
  end do

contains

  subroutine chebyshev_nodes(N, x)
    integer, intent(in) :: N
    real(dp), intent(out) :: x(N)
    integer :: k

    do k = 1, N
      x(k) = cos(pi * (k - 1) / (N - 1))
    end do
  end subroutine chebyshev_nodes

  subroutine chebyshev_differentiation_matrix(N, D)
    integer, intent(in) :: N
    real(dp), intent(out) :: D(N,N)
    integer :: i, j
    real(dp) :: c(N), x(N)

    ! Chebyshev nodes
    do i = 1, N
      x(i) = cos(pi * (i - 1) / (N - 1))
      if (i == 1 .or. i == N) then
        c(i) = 2.0_dp
      else
        c(i) = 1.0_dp
      endif
    end do

    ! Compute differentiation matrix
    do i = 1, N
      do j = 1, N
        if (i /= j) then
          D(i,j) = c(i) / c(j) * (-1.0_dp)**(i + j) / (x(i) - x(j))
        else
          if (i == 1) then
            D(i,j) = (2.0_dp * (N - 1)**2 + 1.0_dp) / 6.0_dp
          elseif (i == N) then
            D(i,j) = -(2.0_dp * (N - 1)**2 + 1.0_dp) / 6.0_dp
          else
            D(i,j) = -x(i) / (2.0_dp * (1.0_dp - x(i)**2))
          endif
        endif
      end do
    end do
  end subroutine chebyshev_differentiation_matrix

  subroutine solve_linear_system(A, b, x)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    integer :: info

    ! Solve the linear system using LAPACK
    call gesv(size(A,1), 1, A, size(A,1), ipiv, b, size(b,1), info)
    if (info /= 0) then
      print *, "Error: LAPACK solver failed."
      stop
    end if
    x = b
  end subroutine solve_linear_system

end program chebyshev_spectral

