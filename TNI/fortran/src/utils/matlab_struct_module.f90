!--------------------------------------------------------------------------~*
!! Copyright (c) 2023 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file matlab_struct_module.f90
!! @author Oleg Korobkin
!! @date  October 2023
!! @brief Data structures and functions inspired by Matlab toolbox
!!


module matlab_struct_module
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

  !> Cumulative sum
  !!
  !! Syntax:
  !!  B = cumsum(A)
  !!  //B = cumsum(A,dim)             (not implemented)
  !!  //B = cumsum(___,direction)
  !!  //B = cumsum(___,nanflag)
  !!
  !! Description:
  !!
  !! `B = cumsum(A)` returns the cumulative sum of A starting at the beginning
  !! of the first array dimension in A whose size is greater than 1.
  !!
  !! If A is a vector, then B is a vector of the same size containing
  !! the cumulative sum of A.
  !!
  !! // If A is a matrix, then B is a matrix of the same size containing (not implemented)
  !! // the cumulative sum in each column of A.
  !!
  !! // If A is a multidimensional array, then B is an array of the same
  !! // size containing the cumulative sum along the first array dimension
  !! // of A whose size is greater than 1.
  !!
  !! // If A is a table or timetable, then M is a table or timetable of
  !! // the same size containing the cumulative sum in each variable
  !! // of A. (since R2023a)
  !!
  function cumsum(A) result(B)
  integer, dimension(:), intent(in):: A
  integer, allocatable, dimension(:):: B
  integer:: k, d, A1

     d= size(A)
     allocate(B(d+1))
     B(1)= A(1)
     do k=1,d
        B(k+1)= B(k) + A(k)
     enddo

  end function cumsum

  !> Utility function to produce pretty matrix output
  !!
  !! If matrix is too large, columns or rows are skipped
  !! Arguments:
  !! - A(:,:)    : matrix to print
  !! - frmt_     :[optional, default = (F8.4)] format for each value
  !! - line_len_ :[optional, default = 80] maximum length of the line
  !! - max_rows_ :[optional, default = 10] maximum number of rows to print
  !! - out_unit_ :[optional, default = stdout] direct output to out_unit_
  !! - strings_  :[optional] write output to array of strings_ instead of unit
  !!
  subroutine pprint_matrix(A,frmt_,line_len_,max_rows_,out_unit_,strings_)
  use iso_fortran_env
  implicit none
  double precision, intent(in)  :: A(:,:)
  character(len=*), intent(in), optional :: frmt_
  integer, intent(in), optional :: line_len_
  integer, intent(in), optional :: max_rows_
  integer, intent(in), optional :: out_unit_
  character(:), allocatable, intent(inout), optional :: strings_(:)
  !
  character(len=20) :: frmt  ! format of each number
  integer :: line_len        ! maximum width of the line
  integer :: max_rows        ! maximum number of rows to print
  integer :: out_unit        ! output unit
  integer :: i, j, imx, jmx, Nrows, Ncols, vlen, strpos, dots_width
  logical :: too_wide, too_long
  character(len=50) :: teststr
  character(:), allocatable :: strings(:)

     Nrows = size(A, 1)
     Ncols = size(A, 2)
     frmt = "(F8.4)"; if (present(frmt_)) frmt= trim(frmt_)
     out_unit = output_unit; if (present(out_unit_)) out_unit= out_unit_
     line_len = 80; if (present(line_len_)) line_len= line_len_
     max_rows = 10; if (present(max_rows_)) max_rows= max_rows_
     if (frmt(2:2).eq.'I' .or. frmt(2:2).eq.'i') then
        write (teststr, frmt) 1
     else
        write (teststr, frmt) -1.2d54
     endif
     vlen = len(trim(teststr))
     too_wide = line_len .lt. Ncols*(vlen + 1) - 1
     too_long = max_rows .lt. Nrows
     imx = Nrows
     if (too_long) then
        imx = max_rows - 2
     else
        max_rows = imx
     endif
     dots_width = 9
     jmx = Ncols
     if (too_wide) then
        jmx = (line_len - dots_width)/(vlen + 1) - 1
        line_len = (vlen + 1)*(jmx + 1) + dots_width
     else
        line_len = (vlen + 1)*jmx
     endif
     if (allocated(strings)) deallocate(strings)
     allocate(character(len=line_len) :: strings(max_rows))
     write (teststr, '("(",I5,"("" ""))")') line_len
     print_matrix_line: do i=1,imx
        write (strings(i), teststr)
        strpos = 1
        do j=1,jmx
           if (frmt(2:2).eq.'I' .or. frmt(2:2).eq.'i') then
              write (strings(i)(strpos:strpos+vlen), frmt) int(A(i,j))
           else
              write (strings(i)(strpos:strpos+vlen), frmt) A(i,j)
           endif
           strpos = strpos + vlen + 1
        enddo
        if (too_wide) then
           write (strings(i)(strpos:strpos+dots_width-1), '("  . . .  ")')
           strpos = strpos + dots_width
           if (strpos + vlen .gt. line_len) exit print_matrix_line
           if (frmt(2:2).eq.'I' .or. frmt(2:2).eq.'i') then
              write (strings(i)(strpos:strpos+vlen), frmt) int(A(i, Ncols))
           else
              write (strings(i)(strpos:strpos+vlen), frmt) A(i, Ncols)
           endif
        endif
     enddo print_matrix_line
     if (too_long) then
        ! print line of dots
        write (strings(imx+1), teststr)
        strpos = 3
        do j=1,((vlen + 1)*(jmx + 1) + 9)/2
           write (strings(imx+1)(strpos:strpos+1), '(". ")')
           strpos = strpos + 2
        enddo

        ! last line
        write (strings(imx+2), teststr)
        strpos = 1
        do j=1,jmx
           if (frmt(2:2).eq.'I' .or. frmt(2:2).eq.'i') then
              write (strings(imx+2)(strpos:strpos+vlen), frmt) int(A(Nrows,j))
           else
              write (strings(imx+2)(strpos:strpos+vlen), frmt) A(Nrows,j)
           endif
           strpos = strpos + vlen + 1
        enddo
        if (too_wide) then
           write (strings(imx+2)(strpos:strpos+8), '("  . . .  ")')
           strpos = strpos + 9
           if (frmt(2:2).eq.'I' .or. frmt(2:2).eq.'i') then
              write (strings(imx+2)(strpos:strpos+vlen), frmt) int(A(Nrows, Ncols))
           else
              write (strings(imx+2)(strpos:strpos+vlen), frmt) A(Nrows, Ncols)
           endif
        endif
     endif

     if (.not.present(strings_)) then
        do i=1,max_rows
           write (out_unit, "(A)") strings(i)
        enddo
        deallocate(strings)
     else
        strings_ = strings(1:max_rows)
     endif

  end subroutine pprint_matrix


  !> Utility function to produce pretty output of a 3D-matrix
  !!
  !! If matrix is too large, columns or rows are skipped
  !! Arguments:
  !! - A(:,:,:) : matrix to print
  !! - frmt_     ;[optional, default = (F8.4)] format for each value
  !! - line_len_ :[optional, default = 80] maximum length of the line
  !! - max_rows_ :[optional, default = 40] maximum number of rows to print
  !! - out_unit_ :[optional, default = stdout] direct output to out_unit_
  !! - flip23_   :[optional, default = false] flip 2nd and 3rd index
  !! - strings_  :[optional] write output to array of strings_ instead of unit
  !!
  subroutine pprint_matrix3d(A,frmt_,line_len_,max_rows_,out_unit_,flip23_,strings_)
  use iso_fortran_env
  implicit none
  double precision, intent(in) :: A(:,:,:)
  character(len=*), intent(in), optional :: frmt_
  integer, intent(in), optional :: line_len_
  integer, intent(in), optional :: max_rows_
  integer, intent(in), optional :: out_unit_
  logical, intent(in), optional :: flip23_
  character(:), allocatable, intent(out), optional :: strings_(:)
  !
  character(len=20) :: frmt  ! format of each number
  integer :: line_len        ! maximum width of the line
  integer :: max_rows        ! maximum number of rows to print
  integer :: out_unit        ! output unit [output_unit]
  logical :: flip23          ! transpose the last two indices of A(ikj)
  integer :: i, j, k, imx, jmx, Nrows, Ncols, Nblocks, vlen, strpos
  integer :: bpos1, bpos2, block_width, ll2, ll3, nb
  logical :: too_wide, too_long, too_many_blocks
  character(len=50) :: teststr
  character(:), allocatable :: strings(:), ss(:)

     frmt = "(F8.4)"; if (present(frmt_)) frmt= trim(frmt_)
     out_unit = output_unit; if (present(out_unit_)) out_unit= out_unit_
     line_len = 80; if (present(line_len_)) line_len= line_len_
     max_rows = 20; if (present(max_rows_)) max_rows= max_rows_
     flip23 = .false.; if (present(flip23_)) flip23 = flip23_
     Nrows = size(A, 1)
     if (flip23) then
        Ncols = size(A, 3); Nblocks = size(A, 2)
     else
        Ncols = size(A, 2); Nblocks = size(A, 3)
     endif
     write (teststr, frmt) -1.2d54
     vlen = len(trim(teststr))

     ! Four possible cases:
     ! 1. all values in all columns and all blocks fit in line_len_
     !    Nblocks*((Ncols*(vlen + 1) - 1) + 3) - 3 < line_len_
     !      1 2 ! 1 2 | 1 2
     !      2 3 ! 2 3 | 3 4
     !
     ! 2. single block fits into line_len_, but all blocks don't
     !    Nblocks*((Ncols*(vlen + 1) - 1) + 3) - 3 > line_len_
     !                                           => "too_many_blocks => True"
     !    Ncols*(vlen + 1) - 1 < line_len_
     !    - subcase 2a: line_len_ can fit 4 blocks
     !      1 2 ! 1 2 |     | 1 2
     !      2 3 ! 2 3 | * * | 3 4
     !
     !    - subcase 2b: line_len_ cannot fit 4 blocks
     !      1 2 ... 4 | 1 2 ... 4 |     | 1 2 ... 4|
     !      1 2 ... 4 | 1 2 ... 4 | * * | 1 2 ... 4|
     !      1 2 ... 4 | 1 2 ... 4 |     | 1 2 ... 4|
     !
     ! 3. single block doesn't fit into line_len_
     !    Ncols*(vlen + 1) - 1 > line_len_   => "too_wide => True"
     !    > this case is printed similarly to 2b
     !
     block_width = Ncols*(vlen + 1) - 1
     too_wide = block_width .gt. line_len
     too_many_blocks = Nblocks*(block_width + 3) - 3 .gt. line_len
     too_long = max_rows .lt. Nrows
     imx = Nrows
     if (too_long) then
        imx = max_rows - 2
     else
        max_rows = imx
     endif

     if (.not.too_wide .and..not.too_many_blocks) then
        ! case 1: print
        !      1 2 ! 1 2 | 1 2
        !      2 3 ! 2 3 | 3 4
        line_len= Nblocks*(block_width + 3) - 3
        allocate(character(len=line_len) :: strings(0:max_rows))
        write (teststr, '("(",I5,"("" ""))")') line_len
        do i=0,max_rows
           write (strings(i), teststr)
        enddo
        bpos1= 1
        do k=1,Nblocks
           if (flip23) then
              write (strings(0)(bpos1:), '("(:,",I3,",:):")') k
              call pprint_matrix(A(:,k,:),frmt,block_width,max_rows,out_unit,ss)
           else
              write (strings(0)(bpos1:), '("(:,:,",I3,"):")') k
              call pprint_matrix(A(:,:,k),frmt,block_width,max_rows,out_unit,ss)
           endif
           do i=1,max_rows
              bpos2= bpos1 + len(ss(i)) - 1
              strings(i)(bpos1:bpos2)= ss(i)
              if (k.lt.Nblocks) then
                 bpos2 = bpos1 + block_width + 1
                 strings(i)(bpos2:bpos2) = '|'
              endif
           enddo
           bpos1= bpos1 + block_width + 3
        enddo
     else if (3*(block_width + 3) + 11 - 3 .le. line_len) then
        ! subcase 2a: line_len_ can fit 4 blocks
        !  1 2 ! 1 2 |     | 1 2
        !  2 3 ! 2 3 | * * | 3 4
        block_width = Ncols*(vlen + 1)
        line_len = 3*(block_width + 3) + 11 - 3
        allocate(character(len=line_len) :: strings(0:max_rows))
        write (teststr, '("(",I5,"("" ""))")') line_len
        do i=0,max_rows
           write (strings(i), teststr)
        enddo
        bpos1= 1
        k= 1
        do while (k.le.Nblocks)
           if (k.eq.3) then
              k = Nblocks
              bpos1= line_len - block_width
              bpos2= bpos1 - 11
              strings(1+max_rows/2)(bpos2:bpos2+8) = '  . . .  '
           endif
           if (flip23) then
              write (strings(0)(bpos1:), '("(:,",I3,",:):")') k
              call pprint_matrix(A(:,k,:),frmt,block_width,max_rows,out_unit,ss)
           else
              write (strings(0)(bpos1:), '("(:,:,",I3,"):")') k
              call pprint_matrix(A(:,:,k),frmt,block_width,max_rows,out_unit,ss)
           endif
           do i=1,max_rows
              bpos2= bpos1 + len(ss(i)) - 1
              strings(i)(bpos1:bpos2)= ss(i)
              if (k.lt.Nblocks) then
                 bpos2 = bpos1 + block_width + 1
              else
                 bpos2 = bpos1 - 2
              endif
              strings(i)(bpos2:bpos2) = '|'
           enddo
           bpos1= bpos1 + block_width + 3
           k= k + 1
        enddo

     else
        ! - subcase 2b: each block is large, and there are many of them
        !   1 2 ... 4 | 1 2 ... 4 |     | 1 2 ... 4|
        !   1 2 ... 4 | 1 2 ... 4 | * * | 1 2 ... 4|
        !   1 2 ... 4 | 1 2 ... 4 |     | 1 2 ... 4|
        !
        block_width = 3*(vlen + 2) + 7
        line_len = 3*(block_width + 3) + 11 - 3
        allocate(character(len=line_len) :: strings(0:max_rows))
        write (teststr, '("(",I5,"("" ""))")') line_len
        do i=0,max_rows
           write (strings(i), teststr)
        enddo
        bpos1= 1; k= 1
        do while (k.le.Nblocks)
           if (k.eq.3) then
              k = Nblocks
              bpos1= line_len - block_width
              bpos2= bpos1 - 11
              strings(1+max_rows/2)(bpos2:bpos2+8) = '  . . .  '
           endif
           if (flip23) then
              write (strings(0)(bpos1:), '("(:,",I3,",:):")') k
              call pprint_matrix(A(:,k,:),frmt,block_width,max_rows,out_unit,ss)
           else
              write (strings(0)(bpos1:), '("(:,:,",I3,"):")') k
              call pprint_matrix(A(:,:,k),frmt,block_width,max_rows,out_unit,ss)
           endif
           do i=1,max_rows
              bpos2= bpos1 + len(ss(i)) - 1
              strings(i)(bpos1:bpos2)= ss(i)
              if (k.lt.Nblocks) then
                 bpos2 = bpos1 + block_width + 1
              else
                 bpos2 = bpos1 - 2
              endif
              strings(i)(bpos2:bpos2) = '|'
           enddo
           bpos1= bpos1 + block_width + 3
           k= k + 1
        enddo
     endif

     if (.not.present(strings_)) then
        do i=0,max_rows
           write (out_unit, "(A)") strings(i)
        enddo
        deallocate(strings)
     else
        strings_ = strings(1:max_rows)
     endif
     if (allocated(ss)) deallocate(ss)

  end subroutine pprint_matrix3d

end module matlab_struct_module

