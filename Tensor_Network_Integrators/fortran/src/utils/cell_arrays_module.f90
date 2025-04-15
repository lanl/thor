!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file cell_array_module.f90
!! @author Oleg Korobkin
!! @date  October 2023
!! @brief Data structures and functions inspired by Matlab toolbox
!!


module cell_arrays_module
implicit none

type array2d
   double precision, allocatable:: arr(:,:)
end type array2d

type array3d
   double precision, allocatable:: arr(:,:,:)
end type array3d

type array4d
   double precision, allocatable:: arr(:,:,:,:)
end type array4d


type cell3d_array
   type(array3d), allocatable:: cells(:)
contains
   final :: delete_cell3d_array
end type cell3d_array
interface cell3d_array
   procedure :: new_cell3d_array
end interface cell3d_array

type cell4d_array
   type(array4d), allocatable:: cells(:)
contains
   final :: delete_cell4d_array
end type cell4d_array
interface cell4d_array
   procedure :: new_cell4d_array
end interface cell4d_array

contains ! --- module  ------------------------------------
  !> Generate cell of 3D arrays with random numbers
  function rand_cell3d(n, r_, rank_) result(y)
  integer, allocatable, intent(in) :: n(:)
  integer, allocatable, intent(in), optional :: r_(:) ! array of ranks
  integer, intent(in), optional :: rank_              ! one rank for all
  type(array3d), allocatable :: y(:)
  !
  integer :: d, rank, i
  integer, allocatable :: r(:)

    rank = 1; if (present(rank_)) rank = rank_
    d = size(n)
    allocate(y(d))

    if (present(r_)) then
       r = r_
    else
       allocate(r(d+1))
       r = rank; r(1) = 1; r(d+1) = 1
    endif
    do i=1,d
       allocate(y(i)%arr(r(i),n(i),r(i+1)))
       call random_number(y(i)%arr)
    enddo
    deallocate(r)

  end function rand_cell3d

  !> Generate cell of 4D arrays with random numbers
  function rand_cell4d(m, n, r_, rank_) result(y)
  integer, allocatable, intent(in) :: m(:), n(:)
  integer, allocatable, intent(in), optional :: r_(:) ! array of ranks
  integer, intent(in), optional :: rank_              ! one rank for all
  type(array4d), allocatable :: y(:)
  !
  integer :: d, rank, i
  integer, allocatable :: r(:)

    rank = 1; if (present(rank_)) rank = rank_
    d = size(n)
    allocate(y(d))

    if (present(r_)) then
       r = r_
    else
       allocate(r(d+1))
       r = rank; r(1) = 1; r(d+1) = 1
    endif
    do i=1,d
       allocate(y(i)%arr(r(i),m(i),n(i),r(i+1)))
       call random_number(y(i)%arr)
    enddo
    deallocate(r)

  end function rand_cell4d

  type(cell3d_array) function new_cell3d_array(n, r)
  integer, intent(in):: n(:), r(:)
  integer k, d

     d= size(n)
     allocate(new_cell3d_array% cells(d))
     do k=1,d
        allocate(new_cell3d_array% cells(k)% arr(r(k),n(k),r(k+1)))
     enddo

  end function new_cell3d_array

  subroutine delete_cell3d_array(self)
  type(cell3d_array), intent(inout) :: self
  integer k
     if(allocated(self% cells)) then
        do k=1,size(self% cells)
           deallocate(self% cells(k)% arr)
        enddo
        deallocate(self% cells)
     endif
  end subroutine delete_cell3d_array

  type(cell4d_array) function new_cell4d_array(n, m, r)
  integer, intent(in):: n(:), m(:), r(:)
  integer k, d

     d= size(n)
     allocate(new_cell4d_array% cells(d))
     do k=1,d
        allocate(new_cell4d_array% cells(k)% arr(r(k),n(k),m(k),r(k+1)))
     enddo

  end function new_cell4d_array

  subroutine delete_cell4d_array(self)
  type(cell4d_array), intent(inout) :: self
  integer k
     if(allocated(self% cells)) then
        do k=1,size(self% cells)
           deallocate(self% cells(k)% arr)
        enddo
        deallocate(self% cells)
     endif
  end subroutine delete_cell4d_array

  !> Cumulative sum: returns an array of cumulative sums of A
  function cumsum(A) result(B)
  integer, intent(in)  :: A(:)
  integer, allocatable :: B(:)
  integer:: k, d, A1

     d= size(A)
     allocate(B(d+1))
     B(1)= A(1)
     do k=1,d
        B(k+1)= B(k) + A(k)
     enddo

  end function cumsum


end module cell_arrays_module

