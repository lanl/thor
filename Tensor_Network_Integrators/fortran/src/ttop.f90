!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file ttop.f90
!! @author Oleg Korobkin, korobkin@lanl.gov
!! @author Ismael Djibrilla Boureima iboureima@lanl.gov
!! @date   August 2024
!! @brief  Common operators: Laplace, Finite differencing, shirt etc.
!!
module ttop_lib
 use thor_lib
 use qtt_lib
 implicit none
 interface qtt_laplace_dd
  module procedure qtt_laplace_dd_1dim, qtt_laplace_dd_ndim
 end interface

contains

  !!
  !! Shift operator
  !!

  !> Creates a shift operator
  !!
  !!           / 0 1 0 0 0 0 0 0 \
  !!           | 0 0 1 0 0 0 0 0 |
  !!           | 0 0 0 1 0 0 0 0 |
  !! m = -1:   | 0 0 0 0 1 0 0 0 |
  !!           | 0 0 0 0 0 1 0 0 |
  !!           | 0 0 0 0 0 0 1 0 |
  !!           | 0 0 0 0 0 0 0 1 |
  !!           \ 0 0 0 0 0 0 0 0 /
  !!
  !! For an integer m, returns the non-periodic downward m-position
  !! shift matrix of size 2^d x 2^d, in the QTT format
  !! The shift is (-m)-position upward for m < 0
  !!  - m = 0 corresponds to the identity matrix
  !!  - m > 0 generates lower-diagonal shifts;
  !!          2^d corresponds to the zero matrix
  !!  - m < 0 generates upper-diagonal shifts;
  !!         -2^d corresponds to the zero matrix
  !! TODO: only implemented for m = 1 or -1
  function qtt_shift(d, m) result(s)
  use qtt_lib
  implicit none
  integer, intent(in) :: d !< dimensionality of the QTT tensor (size = 2^d)
  integer, intent(in) :: m !< shift
  type(qtt_matrix) :: s
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2)
  integer, allocatable :: r(:)

     allocate(r(0:d))
     r(0)= 1; r(1:d-1)= 2; r(d)= 1
     s = qtt_matrix(d, r_=r)
     deallocate(r)

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)

     if (d.eq.1) then
        select case(m)
        case(0)
           s% u4(1)% p(1,:,:,1)= I
        case(1)
           s% u4(1)% p(1,:,:,1)= J1
        case(-1)
           s% u4(1)% p(1,:,:,1)= J
        end select
        return
     endif

     if (m.eq.-1) then
        s% u4(1)% p(1,:,:,1)= J
        s% u4(1)% p(1,:,:,2)= J1
        do k=2,d-1
           s% u4(k)% p(1,:,:,1)= I
           s% u4(k)% p(2,:,:,1)= J
           s% u4(k)% p(2,:,:,2)= J1
        enddo
        s% u4(d)% p(1,:,:,1)= I
        s% u4(d)% p(2,:,:,1)= J

     elseif (m.eq.1) then
        s% u4(1)% p(1,:,:,1)= J1
        s% u4(1)% p(1,:,:,2)= J
        do k=2,d-1
           s% u4(k)% p(1,:,:,1)= I
           s% u4(k)% p(2,:,:,1)= J1
           s% u4(k)% p(2,:,:,2)= J
        enddo
        s% u4(d)% p(1,:,:,1)= I
        s% u4(d)% p(2,:,:,1)= J1

     else
        error stop "ERROR qtt_shift(): this value for m is not implemented"
     endif

  end function qtt_shift


  !> Finite difference operator (1st order) cf. Kazeev (2010)
  !!
  !!       /-1  1  0  0  0  0  0  0 \
  !!       | 0 -1  1  0  0  0  0  0 |
  !!       | 0  0 -1  1  0  0  0  0 |
  !! FD1 = | 0  0  0 -1  1  0  0  0 | / h
  !!       | 0  0  0  0 -1  1  0  0 |
  !!       | 0  0  0  0  0 -1  1  0 |
  !!       | 0  0  0  0  0  0 -1  1 |
  !!       \ 0  0  0  0  0  0  0 -1 /
  !!                              ^(d-2)
  !!                      [ I  0 ]       [I]
  !!     = [(J - I)  J' ]x|      |  |x|  | | / h
  !!                      [ J  J']       [J]
  !!
  !!
  !! Input:
  !!  - d : dimensionality exponent of the matrix (2^d x 2^d)
  !!  - h_: grid step (optional; default = 1)
  !!
  function qtt_fd1(d, h_) result(fd)
  use qtt_lib
  implicit none
  integer, intent(in) :: d !< dimensionality of the QTT tensor (size = 2^d)
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: fd
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2)
  integer, allocatable :: r(:)
  double precision :: h

     h= 1d0; if (present(h_)) h = h_

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)

     allocate(r(0:d))
     r(0:d)= 2; r(0)= 1; r(d)= 1
     fd = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        fd% u4(1)% p(1,:,:,1)= (J - I)/h
        return
     endif

     fd% u4(1)% p(1,:,:,1)= J - I
     fd% u4(1)% p(1,:,:,2)= J1
     do k=2,d-1
        fd% u4(k)% p(1,:,:,1)= I
        fd% u4(k)% p(2,:,:,1)= J
        fd% u4(k)% p(2,:,:,2)= J1
     enddo
     fd% u4(d)% p(1,:,:,1)= I/h
     fd% u4(d)% p(2,:,:,1)= J/h

  end function qtt_fd1


  !> Finite difference operator (2nd order)
  !!
  !!       / 0  1  0  0  0  0  0  0 \
  !!       |-1  0  1  0  0  0  0  0 |
  !!       | 0 -1  0  1  0  0  0  0 |
  !! FD2 = | 0  0 -1  0  1  0  0  0 | / (2*h)
  !!       | 0  0  0 -1  0  1  0  0 |
  !!       | 0  0  0  0 -1  0  1  0 |
  !!       | 0  0  0  0  0 -1  0  1 |
  !!       \ 0  0  0  0  0  0 -1  0 /
  !!
  !!                             ^(d-2)
  !!                  [ I  J  J']     [J'- J ]
  !!     = [I  J  J']x|    J'   | |x| |  -J' | / (2*h)
  !!                  [       J ]     [   J  ]
  !!
  !! Input:
  !!  - d : dimensionality exponent of the matrix (2^d x 2^d)
  !!  - h_: grid step (optional; default = 0.5)
  !!
  function qtt_fd2(d, h_) result(fd)
  use qtt_lib
  implicit none
  integer, intent(in) :: d
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: fd
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2)
  integer, allocatable :: r(:)
  double precision :: h

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)

     h= .5d0; if (present(h_)) h = h_

     allocate(r(0:d))
     r(0:d)= 3; r(0)= 1; r(d)= 1
     fd = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        fd% u4(1)% p(1,:,:,1)= (J1 - J)/(2d0*h)
        return
     endif


     fd% u4(1)% p(1,:,:,1)=-J1 + J
     fd% u4(1)% p(1,:,:,2)=-J
     fd% u4(1)% p(1,:,:,3)= J1
     do k=2,d-1
        fd% u4(k)% p(1,:,:,1)= I
        fd% u4(k)% p(2,:,:,1)= J1
        fd% u4(k)% p(3,:,:,1)= J
        fd% u4(k)% p(2,:,:,2)= J
        fd% u4(k)% p(3,:,:,3)= J1
     enddo
     fd% u4(d)% p(1,:,:,1)= I /(2d0*h)
     fd% u4(d)% p(2,:,:,1)= J1/(2d0*h)
     fd% u4(d)% p(3,:,:,1)= J /(2d0*h)

  end function qtt_fd2


  !> Finite difference operator (2nd order)
  !!
  !!       /-2  2  0  0  0  0  0  0 \
  !!       |-1  0  1  0  0  0  0  0 |
  !!       | 0 -1  0  1  0  0  0  0 |
  !! FD2 = | 0  0 -1  0  1  0  0  0 | / (2*h)
  !!       | 0  0  0 -1  0  1  0  0 |
  !!       | 0  0  0  0 -1  0  1  0 |
  !!       | 0  0  0  0  0 -1  0  1 |
  !!       \ 0  0  0  0  0  0 -2  2 /
  !!                                          ^(d-2)
  !!                        [ I  J  J'       ]         [ J'- J    ]
  !!                        [    J'          ]         [  -J'     ]
  !!     = [I  J  J' I1 I2]x|       J        |   |x|   |   J      | / (2*h)
  !!                        [          I1    ]         [ J - 2*I1 ]
  !!                        [             I2 ]         [-J'+ 2*I2 ]
  !!
  !! Input:
  !!  - d : dimensionality exponent of the matrix (2^d x 2^d)
  !!  - h_: grid step (optional; default = 0.5)
  !!
  function qtt_fd21(d, h_) result(fd)
  use qtt_lib
  implicit none
  integer, intent(in) :: d
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: fd
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2), I1(2,2), I2(2,2)
  integer, allocatable :: r(:)
  double precision :: h

     I1= reshape([1d0, 0d0, &
                  0d0, 0d0], [2,2])
     I2= reshape([0d0, 0d0, &
                  0d0, 1d0], [2,2])
     I = I1 + I2
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)

     h= .5d0; if (present(h_)) h = h_

     allocate(r(0:d))
     r(0:d)= 5; r(0)= 1; r(d)= 1
     fd = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        fd% u4(1)% p(1,:,:,1)= (-I1 + J -J1 + I2)/h
        return
     endif


     fd% u4(1)% p(1,:,:,1)= J - J1
     fd% u4(1)% p(1,:,:,2)=-J
     fd% u4(1)% p(1,:,:,3)= J1
     fd% u4(1)% p(1,:,:,4)= J  - 2d0*I1
     fd% u4(1)% p(1,:,:,5)=-J1 + 2d0*I2
     do k=2,d-1
        fd% u4(k)% p(1,:,:,1)= I
        fd% u4(k)% p(2,:,:,1)= J1
        fd% u4(k)% p(3,:,:,1)= J
        fd% u4(k)% p(2,:,:,2)= J
        fd% u4(k)% p(3,:,:,3)= J1
        fd% u4(k)% p(4,:,:,4)= I1
        fd% u4(k)% p(5,:,:,5)= I2
     enddo
     fd% u4(d)% p(1,:,:,1)= I /(2d0*h)
     fd% u4(d)% p(2,:,:,1)= J1/(2d0*h)
     fd% u4(d)% p(3,:,:,1)= J /(2d0*h)
     fd% u4(d)% p(4,:,:,1)= I1/(2d0*h)
     fd% u4(d)% p(5,:,:,1)= I2/(2d0*h)

  end function qtt_fd21


  !> Laplace operator with Dirichlet boundaries (Kazeev 2010, Lemma 2.1)
  !!
  !!       / 2 -1  0  0  0  0  0  0 \
  !!       |-1  2 -1  0  0  0  0  0 |
  !!       | 0 -1  2 -1  0  0  0  0 |
  !! LDD = | 0  0 -1  2 -1  0  0  0 | / h^2
  !!       | 0  0  0 -1  2 -1  0  0 |
  !!       | 0  0  0  0 -1  2 -1  0 |
  !!       | 0  0  0  0  0 -1  2 -1 |
  !!       \ 0  0  0  0  0  0 -1  2 /
  !!                            ^(d-2)
  !!                 [ I  J' J ]     [ 2*I - J - J']
  !!     = [I  J' J]x|    J    | |x| |      -J     | / h^2
  !!                 [       J']     [      -J'    ]
  !!
  !!
  function qtt_laplace_dd_1dim(d, h_) result(lapl)
  use qtt_lib
  implicit none
  integer, intent(in) :: d !< dimensionality of the QTT tensor (size = 2^d)
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: lapl
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2)
  integer, allocatable :: r(:)
  double precision :: h

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)

     h= 1d0; if (present(h_)) h = h_

     allocate(r(0:d))
     r(0:d)= 3; r(0)= 1; r(d)= 1
     lapl = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        lapl% u4(1)% p(1,:,:,1)= (2d0*I - J1 - J)/(h*h)
        return
     endif

     lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
     lapl% u4(1)% p(1,:,:,2)=-J1
     lapl% u4(1)% p(1,:,:,3)=-J
     do k=2,d-1
        lapl% u4(k)% p(1,:,:,1)= I
        lapl% u4(k)% p(2,:,:,1)= J
        lapl% u4(k)% p(3,:,:,1)= J1
        lapl% u4(k)% p(2,:,:,2)= J1
        lapl% u4(k)% p(3,:,:,3)= J
     enddo
     lapl% u4(d)% p(1,:,:,1)= I / (h*h)
     lapl% u4(d)% p(2,:,:,1)= J / (h*h)
     lapl% u4(d)% p(3,:,:,1)= J1/ (h*h)

  end function qtt_laplace_dd_1dim


  !> Laplace operator with periodic boundaries (Kazeev 2010, Lemma 2.2)
  !!
  !!       / 2 -1  0  0  0  0  0 -1 \
  !!       |-1  2 -1  0  0  0  0  0 |
  !!       | 0 -1  2 -1  0  0  0  0 |
  !! LDD = | 0  0 -1  2 -1  0  0  0 | / h^2
  !!       | 0  0  0 -1  2 -1  0  0 |
  !!       | 0  0  0  0 -1  2 -1  0 |
  !!       | 0  0  0  0  0 -1  2 -1 |
  !!       \-1  0  0  0  0  0 -1  2 /
  !!                                ^(d-2)
  !!                     [ I  J' J ]   [ 2*I - J - J']
  !!     = [I  P  P] |x| [    J    ]|x||      -J     | / h^2
  !!                     [       J']   [      -J'    ]
  !!
  !!
  function qtt_laplace_p(d, h_) result(lapl)
  use qtt_lib
  implicit none
  integer, intent(in) :: d !< dimensionality of the QTT tensor (size = 2^d)
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: lapl
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2), P(2,2)
  integer, allocatable :: r(:)
  double precision :: h

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)
     P = reshape([0d0, 1d0, &
                  1d0, 0d0], [2,2])

     h= 1d0; if (present(h_)) h = h_

     allocate(r(0:d))
     r(0)= 1; r(d)= 1; r(1:d-1)= 3
     lapl = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        lapl% u4(1)% p(1,:,:,1)= 2d0*(I - P)/(h*h)
        return
     endif

     lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
     lapl% u4(1)% p(1,:,:,2)=-J1
     lapl% u4(1)% p(1,:,:,3)=-J

     lapl% u4(d)% p(1,:,:,1)= I / (h*h)
     lapl% u4(d)% p(2,:,:,1)= P / (h*h)
     lapl% u4(d)% p(3,:,:,1)= P / (h*h)

     do k=2,d-1
        lapl% u4(k)% p(1,:,:,1)= I
        lapl% u4(k)% p(2,:,:,1)= J
        lapl% u4(k)% p(3,:,:,1)= J1
        lapl% u4(k)% p(2,:,:,2)= J1
        lapl% u4(k)% p(3,:,:,3)= J
     enddo

  end function qtt_laplace_p


  !> Laplace operator with DN boundaries (Kazeev 2010, Lemma 2.2)
  !!
  !!       / 2 -1  0  0  0  0  0  0 \
  !!       |-1  2 -1  0  0  0  0  0 |
  !!       | 0 -1  2 -1  0  0  0  0 |
  !! LDD = | 0  0 -1  2 -1  0  0  0 | / h^2
  !!       | 0  0  0 -1  2 -1  0  0 |
  !!       | 0  0  0  0 -1  2 -1  0 |
  !!       | 0  0  0  0  0 -1  2 -1 |
  !!       \ 0  0  0  0  0  0 -1  1 / / h
  !!                                 ^(d-2)
  !!                    [ I  J' J   ]     [ 2*I - J - J']
  !!     = [I  J' J I2]x|    J      | |x| |      -J     | / h^2
  !!                    [       J'  ]     [      -J'    ]
  !!                    [         I2]     [      -I2*h  ]
  !!
  function qtt_laplace_dn(d, h_) result(lapl)
  use qtt_lib
  implicit none
  integer, intent(in) :: d
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: lapl
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2), I2(2,2)
  integer, allocatable :: r(:)
  double precision :: h, h2

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)
     I2= reshape([0d0, 0d0, &
                  0d0, 1d0], [2,2])

     h= 1d0; if (present(h_)) h = h_
     h2= h*h

     allocate(r(0:d))
     r(0:d)= 4; r(0)= 1; r(d)= 1
     lapl = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        lapl% u4(1)% p(1,:,:,1)= transpose(reshape( &
           [2d0/h2, -1d0/h2, &
           -1d0/h,   1d0/h], [2,2]))
        return
     endif

     lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
     lapl% u4(1)% p(1,:,:,2)=-J1
     lapl% u4(1)% p(1,:,:,3)=-J
     lapl% u4(1)% p(1,:,:,4)=-I2
     do k=2,d-1
        lapl% u4(k)% p(1,:,:,1)= I
        lapl% u4(k)% p(2,:,:,1)= J
        lapl% u4(k)% p(3,:,:,1)= J1
        lapl% u4(k)% p(2,:,:,2)= J1
        lapl% u4(k)% p(3,:,:,3)= J
        lapl% u4(k)% p(4,:,:,4)= I2
     enddo
     lapl% u4(d)% p(1,:,:,1)= I / h2
     lapl% u4(d)% p(2,:,:,1)= J / h2
     lapl% u4(d)% p(3,:,:,1)= J1/ h2
     lapl% u4(d)% p(4,:,:,1)= I2/ h !< 1st-order derivative

  end function qtt_laplace_dn


  !> Laplace operator with ND boundaries (Kazeev 2010, Lemma 2.2)
  !!
  !!       / 1 -1  0  0  0  0  0  0 \ / h
  !!       |-1  2 -1  0  0  0  0  0 |
  !!       | 0 -1  2 -1  0  0  0  0 |
  !! LDD = | 0  0 -1  2 -1  0  0  0 | / h^2
  !!       | 0  0  0 -1  2 -1  0  0 |
  !!       | 0  0  0  0 -1  2 -1  0 |
  !!       | 0  0  0  0  0 -1  2 -1 |
  !!       \ 0  0  0  0  0  0 -1  2 /
  !!                                 ^(d-2)
  !!                    [ I  J' J   ]     [ 2*I - J - J']
  !!     = [I  J' J I1]x|    J      | |x| |      -J     | / h^2
  !!                    [       J'  ]     [      -J'    ]
  !!                    [         I1]     [      -I1*h  ]
  !!
  function qtt_laplace_nd(d, h_) result(lapl)
  use qtt_lib
  implicit none
  integer, intent(in) :: d
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: lapl
  !
  integer :: k
  double precision :: I(2,2), J(2,2), J1(2,2), I1(2,2)
  integer, allocatable :: r(:)
  double precision :: h, h2

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)
     I1= reshape([1d0, 0d0, &
                  0d0, 0d0], [2,2])

     h= 1d0; if (present(h_)) h = h_
     h2= h*h

     allocate(r(0:d))
     r(0:d)= 4; r(0)= 1; r(d)= 1
     lapl = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        lapl% u4(1)% p(1,:,:,1)= transpose(reshape( &
           [1d0/h, -1d0/h, &
           -1d0/h2, 2d0/h2], [2,2]))
        return
     endif

     lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
     lapl% u4(1)% p(1,:,:,2)=-J1
     lapl% u4(1)% p(1,:,:,3)=-J
     lapl% u4(1)% p(1,:,:,4)=-I1
     do k=2,d-1
        lapl% u4(k)% p(1,:,:,1)= I
        lapl% u4(k)% p(2,:,:,1)= J
        lapl% u4(k)% p(3,:,:,1)= J1
        lapl% u4(k)% p(2,:,:,2)= J1
        lapl% u4(k)% p(3,:,:,3)= J
        lapl% u4(k)% p(4,:,:,4)= I1
     enddo
     lapl% u4(d)% p(1,:,:,1)= I / h2
     lapl% u4(d)% p(2,:,:,1)= J / h2
     lapl% u4(d)% p(3,:,:,1)= J1/ h2
     lapl% u4(d)% p(4,:,:,1)= I1/ h !< 1st-order derivative

  end function qtt_laplace_nd


  !> Laplace operator with NN boundaries (Kazeev 2010, Lemma 2.2)
  !!
  !!       / 1 -1  0  0  0  0  0  0 \ / h
  !!       |-1  2 -1  0  0  0  0  0 |
  !!       | 0 -1  2 -1  0  0  0  0 |
  !! LDD = | 0  0 -1  2 -1  0  0  0 | / h^2
  !!       | 0  0  0 -1  2 -1  0  0 |
  !!       | 0  0  0  0 -1  2 -1  0 |
  !!       | 0  0  0  0  0 -1  2 -1 |
  !!       \ 0  0  0  0  0  0 -1  1 / / h
  !!
  !!                    [ I  J' J     I1]
  !!     = [I  J' J I2]x|    J          |
  !!                    [       J'      ]
  !!                    [         I2 -I1]
  !!                                    ^(d-4)
  !!                    [ I  J' J      ]
  !!                    |    J         |
  !!                |x| [       J'     ]
  !!                    [         I2   ]
  !!                    [            I1]
  !!
  !!        [ I    J'  J      ]     [ 2*I - J - J']
  !!        |      J          |     |             |
  !!        [          J'     ]     [      -J     ] / h^2
  !!    |x| [              I2 ] |x| [             ]
  !!        [ I1  I1  I1      ]     [      -J'    ]
  !!        [--   -   -   -I1 ]     [             ]
  !!        [ 2   2   2       ]     [      -I2    ] / h
  !!
  function qtt_laplace_nn(d, h_) result(lapl)
  use qtt_lib
  implicit none
  integer, intent(in) :: d
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: lapl
  !
  integer :: k,k1,k2
  double precision :: I(2,2),J(2,2),J1(2,2),I1(2,2),I2(2,2)
  integer, allocatable :: r(:)
  double precision :: h, h2

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)
     I1= reshape([1d0, 0d0, &
                  0d0, 0d0], [2,2])
     I2= reshape([0d0, 0d0, &
                  0d0, 1d0], [2,2])

     h= 1d0; if (present(h_)) h = h_; h2= h*h

     allocate(r(0:d))
     r(0)= 1; r(d)= 1
     r(1:d-1)= 4
     if (d.le.3) r(1:d-1)= 5
     r(2:d-2)= 5
     lapl = qtt_matrix(d, r_=r)
     deallocate(r)

     if (d.eq.1) then
        lapl% u4(1)% p(1,:,:,1)= (I - J - J1)/h
        return
     endif

     lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
     lapl% u4(1)% p(1,:,:,2)=-J1
     lapl% u4(1)% p(1,:,:,3)=-J
     lapl% u4(1)% p(1,:,:,4)=-I2

     lapl% u4(d)% p(1,:,:,1)= I / h2
     lapl% u4(d)% p(2,:,:,1)= J / h2
     lapl% u4(d)% p(3,:,:,1)= J1/ h2
     lapl% u4(d)% p(4,:,:,1)= I2/ h

     if (d.gt.3) then
        lapl% u4(2)% p(1,:,:,1)= I
        lapl% u4(2)% p(2,:,:,1)= J
        lapl% u4(2)% p(3,:,:,1)= J1
        lapl% u4(2)% p(2,:,:,2)= J1
        lapl% u4(2)% p(3,:,:,3)= J
        lapl% u4(2)% p(4,:,:,4)= I2
        lapl% u4(2)% p(1,:,:,5)= -.5d0*I1
        lapl% u4(2)% p(2,:,:,5)=  .5d0*I1
        lapl% u4(2)% p(3,:,:,5)=  .5d0*I1
        lapl% u4(2)% p(4,:,:,5)= -I1

        lapl% u4(3)% p(1,:,:,1)= I
        lapl% u4(3)% p(2,:,:,1)= J
        lapl% u4(3)% p(3,:,:,1)= J1
        lapl% u4(3)% p(5,:,:,1)= I1
        lapl% u4(3)% p(2,:,:,2)= J1
        lapl% u4(3)% p(3,:,:,3)= J
        lapl% u4(3)% p(4,:,:,4)= I2
        lapl% u4(3)% p(5,:,:,4)= -I1

        k1= 3; k2= d-2
     else
        lapl% u4(1)% p(1,:,:,5)=-I1
        lapl% u4(d)% p(5,:,:,1)= I1/ h

        k1= 2; k2= d-1
     endif

     do k=k1,k2
        lapl% u4(k)% p(1,:,:,1)= I
        lapl% u4(k)% p(2,:,:,1)= J
        lapl% u4(k)% p(3,:,:,1)= J1
        lapl% u4(k)% p(2,:,:,2)= J1
        lapl% u4(k)% p(3,:,:,3)= J
        lapl% u4(k)% p(4,:,:,4)= I2
        lapl% u4(k)% p(5,:,:,5)= I1
     enddo

  end function qtt_laplace_nn


  !> n-dimensional Laplace, Dirichlet bdry (Kazeev 2010, Cor. 2.4)
  !!
  function qtt_laplace_dd_ndim(d, h_) result(lapl)
  use qtt_lib
  implicit none
  integer, intent(in) :: d(:) !< N-dim QTT tensor {d_1, d_2, ..., d_D}
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: lapl
  !
  integer :: k, m, n, ds, nd
  double precision :: I(2,2), J(2,2), J1(2,2)
  integer, allocatable :: r(:)
  double precision :: h, h2

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = reshape([0d0, 0d0, &
                  1d0, 0d0], [2,2])
     J1= transpose(J)

     h= 1d0; if (present(h_)) h = h_; h2= h*h

     nd= size(d)
     ds= sum(d)
     select case(nd)
     case(1)
        lapl = qtt_laplace_dd_1dim(ds, h_)

     case(2)
        allocate(r(0:ds))
        r(0)= 1
        r(1)= 3
        r(2:d(1)-1)= 4

        r(d(1))= 2
        r(d(1)+1:ds-1)= 3
        r(ds)= 1

        lapl = qtt_matrix(ds, r_=r)
        deallocate(r)

        select case(d(1))
        case(1)
           lapl% u4(1)% p(1,:,:,1)= (2d0*I - J1 - J)/h2
           lapl% u4(1)% p(1,:,:,2)= I / h2

        case(2)
           lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
           lapl% u4(1)% p(1,:,:,2)=-J1
           lapl% u4(1)% p(1,:,:,3)=-J

           lapl% u4(d(1))% p(1,:,:,1)= I / h2
           lapl% u4(d(1))% p(2,:,:,1)= J / h2
           lapl% u4(d(1))% p(3,:,:,1)= J1/ h2
           lapl% u4(d(1))% p(1,:,:,2)= .5d0*I / h2
           lapl% u4(d(1))% p(2,:,:,2)=-.5d0*I / h2
           lapl% u4(d(1))% p(3,:,:,2)=-.5d0*I / h2

        case default
           lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
           lapl% u4(1)% p(1,:,:,2)=-J1
           lapl% u4(1)% p(1,:,:,3)=-J

           lapl% u4(2)% p(1,:,:,1)= I
           lapl% u4(2)% p(2,:,:,1)= J
           lapl% u4(2)% p(3,:,:,1)= J1
           lapl% u4(2)% p(2,:,:,2)= J1
           lapl% u4(2)% p(3,:,:,3)= J
           lapl% u4(2)% p(1,:,:,4)= .5d0*I
           lapl% u4(2)% p(2,:,:,4)=-.5d0*I
           lapl% u4(2)% p(3,:,:,4)=-.5d0*I

           do k=3,d(1)-1
              lapl% u4(k)% p(1,:,:,1)= I
              lapl% u4(k)% p(2,:,:,1)= J
              lapl% u4(k)% p(3,:,:,1)= J1
              lapl% u4(k)% p(2,:,:,2)= J1
              lapl% u4(k)% p(3,:,:,3)= J
              lapl% u4(k)% p(4,:,:,4)= I
           enddo

           lapl% u4(d(1))% p(1,:,:,1)= I / h2
           lapl% u4(d(1))% p(2,:,:,1)= J / h2
           lapl% u4(d(1))% p(3,:,:,1)= J1/ h2
           lapl% u4(d(1))% p(4,:,:,2)= I / h2
        end select

        select case(d(2))
        case (1)
           lapl% u4(ds)% p(1,:,:,1)= I / h2
           lapl% u4(ds)% p(2,:,:,1)= (2d0*I - J1 - J)/h2

        case (2)
           lapl% u4(d(1)+1)% p(1,:,:,1)= I
           lapl% u4(d(1)+1)% p(2,:,:,1)= 2d0*I - J1 - J
           lapl% u4(d(1)+1)% p(2,:,:,2)=-J1
           lapl% u4(d(1)+1)% p(2,:,:,3)=-J

           lapl% u4(ds)% p(1,:,:,1)= I / h2
           lapl% u4(ds)% p(2,:,:,1)= J / h2
           lapl% u4(ds)% p(3,:,:,1)= J1/ h2

        case default
           lapl% u4(d(1)+1)% p(1,:,:,1)= I
           lapl% u4(d(1)+1)% p(2,:,:,1)= 2d0*I - J1 - J
           lapl% u4(d(1)+1)% p(2,:,:,2)=-J1
           lapl% u4(d(1)+1)% p(2,:,:,3)=-J

           do k=d(1)+2,ds-1
              lapl% u4(k)% p(1,:,:,1)= I
              lapl% u4(k)% p(2,:,:,1)= J
              lapl% u4(k)% p(3,:,:,1)= J1
              lapl% u4(k)% p(2,:,:,2)= J1
              lapl% u4(k)% p(3,:,:,3)= J
           enddo

           lapl% u4(ds)% p(1,:,:,1)= I / h2
           lapl% u4(ds)% p(2,:,:,1)= J / h2
           lapl% u4(ds)% p(3,:,:,1)= J1/ h2
        end select

     case default ! nd > 2

        allocate(r(0:ds))
        r(0)= 1
        r(1)= 3
        r(2:d(1)-1)= 4

        m= d(1)
        do k=2,nd-1
           r(m)= 2
           r(m+1:m+d(k)-1)= 4
           m= m + d(k)
        enddo

        m= ds - d(nd)
        r(m)= 2
        r(m+1:ds-1)= 3
        r(ds)= 1

        lapl = qtt_matrix(ds, r_=r)
        deallocate(r)

        ! first dimension (2**d(1))
        select case(d(1))
        case(1)
           lapl% u4(1)% p(1,:,:,1)= (2d0*I - J1 - J)/h2
           lapl% u4(1)% p(1,:,:,2)= I / h2

        case(2)
           lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
           lapl% u4(1)% p(1,:,:,2)=-J1
           lapl% u4(1)% p(1,:,:,3)=-J

           lapl% u4(d(1))% p(1,:,:,1)= I / h2
           lapl% u4(d(1))% p(2,:,:,1)= J / h2
           lapl% u4(d(1))% p(3,:,:,1)= J1/ h2
           lapl% u4(d(1))% p(1,:,:,2)= .5d0*I / h2
           lapl% u4(d(1))% p(2,:,:,2)=-.5d0*I / h2
           lapl% u4(d(1))% p(3,:,:,2)=-.5d0*I / h2

        case default
           lapl% u4(1)% p(1,:,:,1)= 2d0*I - J1 - J
           lapl% u4(1)% p(1,:,:,2)=-J1
           lapl% u4(1)% p(1,:,:,3)=-J
           lapl% u4(1)% p(1,:,:,4)= I

           lapl% u4(2)% p(1,:,:,1)= I
           lapl% u4(2)% p(2,:,:,1)= J
           lapl% u4(2)% p(3,:,:,1)= J1
           lapl% u4(2)% p(2,:,:,2)= J1
           lapl% u4(2)% p(3,:,:,3)= J
           lapl% u4(2)% p(1,:,:,4)= .5d0*I
           lapl% u4(2)% p(2,:,:,4)=-.5d0*I
           lapl% u4(2)% p(3,:,:,4)=-.5d0*I

           do k=3,d(1)-1
              lapl% u4(k)% p(1,:,:,1)= I
              lapl% u4(k)% p(2,:,:,1)= J
              lapl% u4(k)% p(3,:,:,1)= J1
              lapl% u4(k)% p(2,:,:,2)= J1
              lapl% u4(k)% p(3,:,:,3)= J
              lapl% u4(k)% p(4,:,:,4)= I
           enddo

           lapl% u4(d(1))% p(1,:,:,1)= I / (h*h)
           lapl% u4(d(1))% p(2,:,:,1)= J / (h*h)
           lapl% u4(d(1))% p(3,:,:,1)= J1/ (h*h)
           lapl% u4(d(1))% p(4,:,:,2)= I / (h*h)

        end select

        ! intermediate dimensions
        m= d(1)
        do n=2,nd-1
           select case(d(n))
           case (1)
              lapl% u4(m+1)% p(1,:,:,1)= I / h2
              lapl% u4(m+1)% p(2,:,:,2)= I / h2
              lapl% u4(m+1)% p(2,:,:,1)= (2d0*I - J1 - J)/h2

           case default
              lapl% u4(m+1)% p(1,:,:,1)= I
              lapl% u4(m+1)% p(2,:,:,1)= 2d0*I - J1 - J
              lapl% u4(m+1)% p(2,:,:,2)=-J1
              lapl% u4(m+1)% p(2,:,:,3)=-J
              lapl% u4(m+1)% p(2,:,:,4)= I

              do k=m+2,m+d(n)-1
                 lapl% u4(k)% p(1,:,:,1)= I
                 lapl% u4(k)% p(2,:,:,1)= J
                 lapl% u4(k)% p(3,:,:,1)= J1
                 lapl% u4(k)% p(2,:,:,2)= J1
                 lapl% u4(k)% p(3,:,:,3)= J
                 lapl% u4(k)% p(4,:,:,4)= I
              enddo

              lapl% u4(m+d(n))% p(1,:,:,1)= I / (h*h)
              lapl% u4(m+d(n))% p(2,:,:,1)= J / (h*h)
              lapl% u4(m+d(n))% p(3,:,:,1)= J1/ (h*h)
              lapl% u4(m+d(n))% p(4,:,:,2)= I / (h*h)
           end select

           m= m + d(n)
        enddo

        ! the last dimension
        select case(d(nd))
        case (1)
           lapl% u4(ds)% p(1,:,:,1)= I / h2
           lapl% u4(ds)% p(2,:,:,1)= (2d0*I - J1 - J)/h2

        case (2)
           lapl% u4(ds-1)% p(1,:,:,1)= I
           lapl% u4(ds-1)% p(2,:,:,1)= 2d0*I - J1 - J
           lapl% u4(ds-1)% p(2,:,:,2)=-J1
           lapl% u4(ds-1)% p(2,:,:,3)=-J

           lapl% u4(ds)% p(1,:,:,1)= I / h2
           lapl% u4(ds)% p(2,:,:,1)= J / h2
           lapl% u4(ds)% p(3,:,:,1)= J1/ h2

        case default
           lapl% u4(m+1)% p(1,:,:,1)= I
           lapl% u4(m+1)% p(2,:,:,1)= 2d0*I - J1 - J
           lapl% u4(m+1)% p(2,:,:,2)=-J1
           lapl% u4(m+1)% p(2,:,:,3)=-J

           do k=m+2,ds-1
              lapl% u4(k)% p(1,:,:,1)= I
              lapl% u4(k)% p(2,:,:,1)= J
              lapl% u4(k)% p(3,:,:,1)= J1
              lapl% u4(k)% p(2,:,:,2)= J1
              lapl% u4(k)% p(3,:,:,3)= J
           enddo

           lapl% u4(ds)% p(1,:,:,1)= I / (h*h)
           lapl% u4(ds)% p(2,:,:,1)= J / (h*h)
           lapl% u4(ds)% p(3,:,:,1)= J1/ (h*h)

        end select

     end select

  end function qtt_laplace_dd_ndim


  !> n-dimensional inverse Laplace, Dirichlet bdry (Kazeev 2010, Lemma 3.3)
  !!
  function qtt_laplace_dd_inv(d, h_) result(lapl)
  use qtt_lib
  implicit none
  integer, intent(in) :: d
  double precision, intent(in), optional :: h_
  type(qtt_matrix) :: lapl
  !
  integer :: m, n
  double precision :: I(2,2), J(2,2), K(2,2), L(2,2), E(2,2), F(2,2), P(2,2)
  double precision :: J1(2,2), L1(2,2)
  integer, allocatable :: r(:)
  double precision :: h, h2, x, y, z

     I = reshape([1d0, 0d0, &
                  0d0, 1d0], [2,2])
     J = transpose(reshape([0d0, 1d0, &
                            0d0, 0d0], [2,2]))
     J1= transpose(J)
     P = J + J1
     K = reshape([-1d0, 0d0, &
                   0d0, 1d0], [2,2])
     E = I + P
     F = I - P
     L = J1- J
     L1= transpose(L)

     h= 1d0; if (present(h_)) h = h_; h2= h*h

     allocate(r(0:d))
     r(0)= 1; r(1:d-1)= 5; r(d)= 1
     lapl = qtt_matrix(d, r_=r)
     deallocate(r)

     lapl% u4(1)% p(1,:,:,1)= (I + E)/3d0
     lapl% u4(1)% p(1,:,:,2)= 2d0*E
     lapl% u4(1)% p(1,:,:,3)= 1d0/18d0*F
     lapl% u4(1)% p(1,:,:,4)= 2d0/3d0*K
     lapl% u4(1)% p(1,:,:,5)=-2d0/3d0*L1

     do n=2,d-1
        x = dble(lshift(1,n-1) + 1)/dble(lshift(1,n) + 1)
        y = dble(lshift(1,n-2))/dble(lshift(1,n) + 1)
        z = (1d0 + 1d0/dble(lshift(1,n-1)))*x

        lapl% u4(n)% p(1,:,:,1)= I
        lapl% u4(n)% p(2,:,:,1)= .25d0*(x*I + z*P)
        lapl% u4(n)% p(3,:,:,1)= x*I - z*P
        lapl% u4(n)% p(4,:,:,1)=-.5d0*x*K
        lapl% u4(n)% p(5,:,:,1)= .5d0*z*L1

        lapl% u4(n)% p(2,:,:,2)= 2d0*E

        lapl% u4(n)% p(2,:,:,3)= 2d0*y*y*F
        lapl% u4(n)% p(3,:,:,3)= 2d0*x*x*E
        lapl% u4(n)% p(4,:,:,3)= 2d0*x*y*K
        lapl% u4(n)% p(5,:,:,3)= 2d0*x*y*L1

        lapl% u4(n)% p(2,:,:,4)= 4d0*y*K
        lapl% u4(n)% p(4,:,:,4)= 2d0*x*E

        lapl% u4(n)% p(2,:,:,5)=-4d0*y*L1
        lapl% u4(n)% p(5,:,:,5)= 2d0*x*E

     enddo

     x = dble(lshift(1,d-1) + 1)/dble(lshift(1,d) + 1)
     y = dble(lshift(1,d-2))/dble(lshift(1,d) + 1)
     z = (1d0 + 1d0/dble(lshift(1,d-1)))*x

     lapl% u4(d)% p(1,:,:,1)= I
     lapl% u4(d)% p(2,:,:,1)= .25d0*(x*I + z*P)
     lapl% u4(d)% p(3,:,:,1)= x*I - z*P
     lapl% u4(d)% p(4,:,:,1)=-.5d0*x*K
     lapl% u4(d)% p(5,:,:,1)= .5d0*z*L1

  end function qtt_laplace_dd_inv

  ! DEL OPERATOR
  type(dtt_tensor) function del_dtt(tt, dir, h, ngc) result(d_tt)
    !
    ! This function applies the del opterator to tt using the centered finite difference method
    ! ARGUMENTS:
    !  tt  - tt-tensor being interpolated
    !  dir - Direction of del
    !  h   - double precision: step size
    !  ngc - integer: number of gard-cells
    !
    use thor_lib 
    implicit none
    type(dtt_tensor), intent(inout)                :: tt !, d_tt
    integer, intent(in)                     :: dir
    double precision,intent(in),optional    :: h
    integer ,intent(in),optional            :: ngc
    integer                                 :: rl,n,rr, gc, info
    double precision                        :: dh
    double precision,parameter                :: zero=0.d0

    if(present(h))   then; dh=2.0*h; else; dh=2.0; endif
    if(present(ngc)) then; gc=ngc;   else; gc=0   ; endif
    d_tt = tt * zero
    rl = tt%r(dir-1); n=tt%n(dir); rr=tt%r(dir)
    d_tt%u(dir)%p(:, gc+1 : n-gc+1, :) =  ( tt%u(dir)%p(:,gc+2:n-gc+2,:) - tt%u(dir)%p(:, gc:n-gc, :) ) / dh
  end function del_dtt

end module ttop_lib
