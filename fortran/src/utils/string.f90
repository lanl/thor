!! @file string.f90
!! @author [OGK]
!! @date  March 2024
!! @brief library of simple string operations
!!
!!--------------------------------------------------------------------------~
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~

module string_lib
  implicit none
  interface str
     module procedure str_i, str_d, str_f
  end interface

contains

  !> splits a string to array of substrings using a delimiter
  !!
  !! Example: split('may;april;june', ss, cdel_=';')
  !! -> returns array ss(1:3) = ["may\0\0","april","june\0"]
  !!
  subroutine split(strarr, instr, cdel_)
  character(:), allocatable, intent(OUT) :: strarr(:)
  character(*), intent(in)  :: instr
  character(*), intent(in), optional:: cdel_ ! string of delimiters [' ']
  !
  character(:), allocatable :: cdel
  integer :: i, j
  integer :: tkn    !< token counter
  integer :: tklen  !< token length: needed to get maximum token length
  integer :: tknum  !< number of tokens
  integer :: ctlen  !< current token ('under cursor') length
  logical :: lcciat !< curent character 'ch' in a token
  logical :: lpciat !< previous character in a token (not a delimiter)
  character(len=1) :: ch

     ! optional parameter initialization
     if (present(cdel_)) then
        allocate(character(len=len(cdel_)) :: cdel)
        write(cdel, '(A)') cdel_
     else
        allocate(character(len=1) :: cdel)
        cdel = ' '
     endif

     ! count the number of tokens and maximum length
     tklen = 0; tknum = 0; ctlen = 0
     lpciat = .false.
     charloop: do i=1,len(instr)
        lcciat = .true.
        delimloop: do j=1,len(cdel)
           lcciat = lcciat.and.instr(i:i).ne.cdel(j:j)
        enddo delimloop
        if (lcciat) then
           if (lpciat) then
              ctlen = ctlen + 1
           else
              ctlen = 1
           endif
        else
           if (lpciat) then
              tknum = tknum + 1
              tklen = max(tklen, ctlen)
           endif
        endif
        lpciat = lcciat
     enddo charloop
     if (lpciat) then
        tknum = tknum + 1
        tklen = max(tklen, ctlen)
     endif

     ! allocate character array
     allocate(character(len=tklen) :: strarr(tknum))
     tkn = 0; ctlen = 0
     lpciat = .false.
     do i=1,len(instr)
        lcciat = .true.
        ch = instr(i:i)
        do j=1,len(cdel)
           lcciat = lcciat.and.ch.ne.cdel(j:j)
        enddo
        if (lcciat.and..not.lpciat) then
           ctlen = 1
           tkn = tkn + 1
           strarr(tkn) = repeat(' ',tklen)
           strarr(tkn)(1:1) = ch
        else
           if (lcciat) then
              ctlen = ctlen + 1
              strarr(tkn)(ctlen:ctlen) = ch
           endif
        endif
        lpciat = lcciat
     enddo

     deallocate(cdel)

  end subroutine split


  !> parses a list of integers from a string
  !!
  !! Example: ints = isplit('12 14 15', cdel_=' ')
  !! -> returns array ints(1:3) = [12, 14, 15]
  !!
  function isplit(instr, cdel_) result(iarr)
  implicit none
  character(*), intent(in):: instr
  character(*), intent(in), optional:: cdel_ ! string of delimiters [' ']
  integer, allocatable :: iarr(:)
  !
  integer :: i
  character(len=:), allocatable :: cdel, token, strarr(:)

     if (present(cdel_)) then
        allocate(character(len=len(cdel_)) :: cdel)
        write(cdel, '(A)') cdel_
     else
        allocate(character(len=3) :: cdel)
        cdel = ' ,;'
     endif

     call split(strarr, trim(instr), cdel_=cdel)
     allocate(iarr(size(strarr)))
     token = "1"
     read (token,*) iarr(1)
     do i=1,size(strarr)
        read (strarr(i),*) iarr(i)
     enddo
     deallocate(cdel, strarr)

  end function isplit


  function str_i(n, frmt_)
  integer, intent(IN) :: n
  character(*), intent(IN), optional :: frmt_
  character(:), allocatable :: str_i
  !
  character(len=20) :: str_n
  character(len=30) :: frmt
     frmt='(I20)'; if (present(frmt_)) frmt= frmt_
     write(str_n, trim(frmt)) n
     str_i= trim(adjustl(str_n))
  end function str_i


  function str_d(x, frmt_)
  double precision, intent(IN) :: x
  character(*), intent(IN), optional :: frmt_
  character(:), allocatable :: str_d
  character(len=25) :: str_x
  character(len=30) :: frmt
     frmt='(ES24.17)'; if (present(frmt_)) frmt= frmt_
     write(str_x, trim(frmt)) x
     str_d= trim(adjustl(str_x))
  end function str_d


  function str_f(x, frmt_)
  real, intent(IN) :: x
  character(*), intent(IN), optional :: frmt_
  character(:), allocatable :: str_f
  character(len=20) :: str_x
  character(len=30) :: frmt
     frmt='(ES14.7)'; if (present(frmt_)) frmt= frmt_
     write(str_x, trim(frmt)) x
     str_f= trim(adjustl(str_x))
  end function str_f

end module string_lib
