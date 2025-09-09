!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file thorio.f90
!! @author Oleg Korobkin, korobkin@lanl.gov
!! @author Ismael Djibrilla Boureima iboureima@lanl.gov
!! @date   March 2025
!! @brief  THOR input/output routines
!!
module thorio_lib

  use thor_lib, only: tt_size, dtt_tensor, dtt_full, alloc
  implicit none

  type, private :: tthead
   sequence
   character(len=8) :: txt='TT      '
   integer(kind=4)  :: ver(2)=(/ 1, 0 /)
   integer(kind=4)  :: inf(4)=(/tt_size, 0, 0, 0/)
   character(len=64):: comment
   integer(kind=4)  :: i(8)
  end type
  integer(kind=4), private :: sdv_format_ver(2)=(/ 1, 0 /)

contains

  subroutine read_array_from_file(fname, nx,ny,nz,arr)
  !! Auxilliary routine for reading array from file
  !! AGUMENTS:
  !!    fname     [str]      :: Name of the file containing the flatten array
  !! OUTPUT:
  !!    nx,ny,ny  [int]      :: X,Y,Z resolution of the array
  !!    arr       [real(8)]  :: Flattened Fortran array being read
  implicit none
  character(len=*),     intent(in)  :: fname
  integer,              intent(out) :: nx, ny,nz
  real(8), allocatable, intent(out) :: arr(:,:,:)
  !
  integer                                 :: i,j,k, info
     ! Read density array
     open (771, file = fname, status = 'old')
     read(771,*) nx, ny, nz
     allocate( arr(nx,ny,nz), stat=info)         ! Helium density
     if(info.ne.0) error stop '[!][read_array_from_file()] Problem allocating'
     do k = 1,nz
        do j = 1, ny
           do i = 1, nx
              read(771,*) arr(i,j,k)
           end do !i
        end do !j
     end do !k
     close(771)
  end subroutine read_array_from_file

  subroutine read_ttv_from_file(fname, d,sz,n,r,ps,cr)
     implicit none
     !character(:), allocatable, intent(in) :: fname
     character(len=*), intent(in)            :: fname
     integer,                   intent(out)  :: d, sz
     integer, allocatable,      intent(out)  :: n(:)
     integer, allocatable,      intent(out)  :: r(:)
     integer, allocatable,      intent(out)  :: ps(:)
     real(8), allocatable,      intent(out)  :: cr(:)
     integer                                 :: info

     ! Read density array
        open (773, file = fname, status = 'old')
            read(773,*) d, sz
            allocate( n(d), r(d+1), ps(d+1), cr(sz), stat=info)
            if(info.ne.0)then;write(*,*)'[!] tt-Vector meta';stop;endif
            read(773,*) n
            read(773,*) r
            read(773,*) ps
            read(773,*) cr
        close(773)
        print*,"[+] tt-Vector ",fname," read OK"
  end subroutine read_ttv_from_file

  subroutine read_ttm_from_file(fname, d,sz,m,n,r,ps,cr)
     implicit none
     !character(:), allocatable, intent(in) :: fname
     character(len=*), intent(in)            :: fname
     integer,                   intent(out)  :: d, sz
     integer, allocatable,      intent(out)  :: m(:)
     integer, allocatable,      intent(out)  :: n(:)
     integer, allocatable,      intent(out)  :: r(:)
     integer, allocatable,      intent(out)  :: ps(:)
     real(8), allocatable,      intent(out)  :: cr(:)
     integer                                 :: info

     ! Read density array
        open (774, file = fname, status = 'old')
            read(774,*) d, sz
            allocate( m(d), n(d), r(d+1), ps(d+1), cr(sz), stat=info)
            if(info.ne.0)then;write(*,*)'[!] tt-Matrix meta';stop;endif
            read(774,*) m
            read(774,*) n
            read(774,*) r
            read(774,*) ps
            read(774,*) cr
        close(774)
        print*,"[+] tt-Matrix ",fname," read OK"
  end subroutine read_ttm_from_file


  !> count the number of lines in a file
  integer function file_nlines(fname) result(nlines)
  character(*), intent(in) :: fname
  integer, parameter :: fd = 776

     nlines = 0
     open (fd, file=fname, status='old')
     do
        read(fd, *, end=10)
        nlines = nlines + 1
     enddo
10   close(fd)
  end function file_nlines


  !> Writes tt-tensor into a file in ASCII format
  !!
  !! Data format:
  !! - first row: tensor dimensions in all modes, separated by commas
  !! - second row: tensor ranks in all modes, separated by commas
  !! - rest of the file: tensor cores, flattened, in sequential order,
  !!   single value per line
  !!
  !! File example:
  !! >>>>>>>>>>>>>>>> BEGIN FILE >>>>>>>>>>>>>>
  !! 128,128,128,128
  !! 1,12,14,12,1
  !! -0.132E+00  \
  !!  1.333E+03  |
  !!  1.333E+03   \  total of 40080 =
  !!  1.333E+03   /   = 12*128 + 12*128*14 + 14*128*12 + 128*1 lines
  !!  ...        |
  !! -1.000E-0.1 /
  !! <<<<<<<<<<<<<<<< END FILE <<<<<<<<<<<<<<<<
  !!
  subroutine dtt_write_ascii_file(tt, fname, verb_)
  implicit none
  class(dtt_tensor), intent(in)  :: tt
  character(*), intent(in)      :: fname
  logical, optional, intent(in) :: verb_
  !
  logical :: verb
  integer :: i, sz, nlines, k, l, m
  integer, parameter :: fd = 778
  integer, allocatable :: n(:), r(:)
  double precision, pointer :: p1(:)
  character(20) :: fmtstr

     verb = .false.; if (present(verb_)) verb = verb_
     if (verb) then
        print '("writing tt-tensor to a file: ",A)', fname
        print '("tt-tensor: ")'
        call tt% say()
     endif

     open (fd, file=fname, status='replace')
     l = tt% l; m = tt% m

     ! write modes
     do k = l, m - 1
        i = tt%n(k)
        write (fmtstr, '("(i",i2,",'','')")') int(log10(dble(i)))+1
        write (fd, fmtstr, advance='no') i
     enddo
     i = tt%n(m)
     write (fmtstr, '("(i",i2,")")') int(log10(dble(i)))+1
     write (fd, fmtstr) i

     ! write ranks
     do k = l - 1, m - 1
        i = tt%r(k)
        write (fmtstr, '("(i",i2,",'','')")') int(log10(dble(i)))+1
        write (fd, fmtstr, advance='no') i
     enddo
     i = tt%r(m)
     write (fmtstr, '("(i",i2,")")') int(log10(dble(i)))+1
     write (fd, fmtstr) i

     ! write cores
     do k=l,m
        sz = tt% r(k-1) * tt% n(k) * tt% r(k)
        p1(1:sz)=> tt% u(k)% p
        do i=1,sz
           write(fd, '(ES22.14)') p1(i)
        enddo
     enddo

     close(fd)

     if (verb) then
        nlines = file_nlines(fname)
        print '(I12," lines written.")', nlines
     endif

  end subroutine dtt_write_ascii_file


  !> Reads tt-tensor from a file in ASCII format
  !!
  !! Data format: see description to dtt_write_ascii_file:
  !!
  function dtt_read_ascii_file(fname, verb_) result(tt)
  use string_lib, only: isplit
  implicit none
  character(*), intent(in)      :: fname
  logical, optional, intent(in) :: verb_
  type(dtt_tensor) :: tt
  !
  logical :: verb
  integer :: i, sz, nlines
  integer, parameter :: fd = 775
  integer, allocatable :: n(:), r(:)
  double precision, allocatable :: cr(:)
  character(:), allocatable :: line, strarr(:)

     verb = .false.; if (present(verb_)) verb = verb_
     if (verb) print '("reading tt-tensor from file: ",A)', fname

     nlines = file_nlines(fname)
     open (fd, file=fname, status='old')
     allocate (character(len=2048*8) :: line)
     read (fd, '(A)') line
     n = isplit(line, cdel_=',')
     read (fd, '(A)') line
     r = isplit(line, cdel_=',')
     sz = 0
     do i=1,size(n)
        sz = sz + r(i)*n(i)*r(i+1)
     enddo
     if (sz.ne.nlines - 2) then
        print '("ERROR: Cannot read the file ",A)', fname
        print '(A22,I12,A6)', "expected file size: ", sz + 2, " lines"
        print '(A22,I12,A6)', "actual file size: ", nlines, " lines"
        error stop
     endif
     allocate(cr(sz))
     do i=1,sz
        read(fd, *) cr(i)
     enddo
     close (fd)
     tt = dtt_tensor(cr, n, r)
     deallocate(line,n,r,cr)
     if (verb) then
        print '("lines read: ",I12)', nlines
        print '("tt-tensor: ")'
        call tt% say()
     endif

  end function dtt_read_ascii_file


  !> Writes tt-matrix into a file in ASCII format
  !!
  !! Data format: same as for tt-tensor (see dtt_write_ascii_file),
  !! but the first two rows are modes q and s of the matrix
  !!
  subroutine dttm_write_ascii_file(ttm, fname, verb_)
  use thor_lib, only: dtt_matrix
  implicit none
  class(dtt_matrix), intent(IN) :: ttm
  character(*), intent(IN)      :: fname
  logical, optional, intent(IN) :: verb_
  !
  logical :: verb
  integer :: i, sz, nlines, k, l, m
  integer, parameter :: fd = 778
  integer, allocatable :: q(:), s(:), r(:)
  double precision, pointer :: p1(:)
  character(20) :: fmtstr

     verb = .false.; if (present(verb_)) verb = verb_
     if (verb) then
        print '("writing tt-matrix to a file: ",A)', fname
        print '("tt-matrix: ")'
        call ttm% say()
     endif

     open (fd, file=fname, status='replace')
     l = ttm% l; m = ttm% m

     ! write modes
     do k = l, m - 1
        i = ttm%q(k)
        write (fmtstr, '("(i",i2,",'','')")') int(log10(dble(i)))+1
        write (fd, fmtstr, advance='no') i
     enddo
     i = ttm%q(m)
     write (fmtstr, '("(i",i2,")")') int(log10(dble(i)))+1
     write (fd, fmtstr) i

     do k = l, m - 1
        i = ttm%s(k)
        write (fmtstr, '("(i",i2,",'','')")') int(log10(dble(i)))+1
        write (fd, fmtstr, advance='no') i
     enddo
     i = ttm%s(m)
     write (fmtstr, '("(i",i2,")")') int(log10(dble(i)))+1
     write (fd, fmtstr) i

     ! write ranks
     do k = l - 1, m - 1
        i = ttm%r(k)
        write (fmtstr, '("(i",i2,",'','')")') int(log10(dble(i)))+1
        write (fd, fmtstr, advance='no') i
     enddo
     i = ttm%r(m)
     write (fmtstr, '("(i",i2,")")') int(log10(dble(i)))+1
     write (fd, fmtstr) i

     ! write cores
     do k=l,m
        sz = ttm% r(k-1) * ttm% n(k) * ttm% r(k)
        p1(1:sz)=> ttm% u(k)% p
        do i=1,sz
           write(fd, '(ES22.14)') p1(i)
        enddo
     enddo

     close(fd)

     if (verb) then
        nlines = file_nlines(fname)
        print '(I12," lines written.")', nlines
     endif

  end subroutine dttm_write_ascii_file


  !> Reads tt-matrix from a file in ASCII format
  !!
  !! Data format: see description to dtt_write_ascii_file:
  !!
  function dttm_read_ascii_file(fname, verb_) result(ttm)
  use thor_lib, only: dtt_matrix
  use string_lib, only: isplit
  implicit none
  character(*), intent(in)      :: fname
  logical, optional, intent(in) :: verb_
  type(dtt_matrix) :: ttm
  !
  logical :: verb
  integer :: i, sz, nlines
  integer, parameter :: fd = 775
  integer, allocatable :: q(:), s(:), r(:)
  double precision, allocatable :: cr(:)
  character(:), allocatable :: line, strarr(:)

     verb = .false.; if (present(verb_)) verb = verb_
     if (verb) print '("reading tt-matrix from file: ",A)', fname

     nlines = file_nlines(fname)
     open (fd, file=fname, status='old')
     allocate (character(len=2048*8) :: line)
     read (fd, '(A)') line
     q = isplit(line, cdel_=',')
     read (fd, '(A)') line
     s = isplit(line, cdel_=',')
     read (fd, '(A)') line
     r = isplit(line, cdel_=',')
     if (size(q)/=size(s).or.size(r)-1/=size(q)) then
        error stop "Incompatible size of q, s, and r in file "//fname
     endif
     sz = 0
     do i=1,size(q)
        sz = sz + r(i)*q(i)*s(i)*r(i+1)
     enddo
     if (sz.ne.nlines - 3) then
        print '("ERROR: Cannot read the file ",A)', fname
        print '(A22,I12,A6)', "expected file size: ", sz + 3, " lines"
        print '(A22,I12,A6)', "actual file size: ", nlines, " lines"
        error stop
     endif
     allocate(cr(sz))
     do i=1,sz
        read(fd, *) cr(i)
     enddo
     close (fd)
     ttm = dtt_matrix(cr, q, s, r)
     deallocate(line,q,s,r,cr)
     if (verb) then
        print '("lines read: ",I12)', nlines
        print '("tt-matrix: ")'
        call ttm% say()
     endif

  end function dttm_read_ascii_file


  !> write fortran SVD format
  !!
  !! Adopted from TT-toolbox 'ttio.f90' file
  !!
  subroutine dtt_write_sdv_file(arg,fname,info)
  implicit none
  class(dtt_tensor),intent(in) :: arg
  character(*),intent(in) :: fname
  integer,intent(out),optional :: info
  !
  character(*),parameter :: subnam='dtt_write_sdv_file'
  type(tthead) :: head
  integer :: io,u,i,j,k,b,sz
  integer(kind=4) :: l,m,n(tt_size),r(0:tt_size)
  real(kind=8),allocatable :: x(:)
  character(len=*),parameter :: frm='unformatted', acc='stream'
  integer,parameter :: un=51
  logical :: ex,op

     if(present(info))info=-11
     inquire (file=fname, exist=ex, opened=op)
     if(op)then
        write(*,*)subnam,': file is open, trying to close: ',fname
        inquire(file=fname, number=u)
        write(*,*)subnam,': establish unit: ',u
        close(unit=u,status='keep')
        write(*,*)subnam,': closed ok'
     endif

     u=un; op=.true.
     do while(op)
        inquire(unit=u,opened=op)
        if(op)then
          write(*,*)subnam,': unit ',u,' is busy, trying next '
          u=u+1
        endif
     enddo

     sz= arg% size()
     if(sz.le.0)write(*,*)subnam,': tt structure has invalid size: ',sz
     allocate(x(max(sz,1)),stat=i)
     if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif

     l=arg%l; m=arg%m
     n=arg%n; r=arg%r
     sz=0
     do b=l,m
        forall(i=1:r(b-1), j=1:n(b), k=1:r(b)) &
           x(sz + i + (j-1)*r(b-1) + (k-1)*r(b-1)*n(b))= arg% u(b)% p(i,j,k)
        sz=sz+r(b-1)*n(b)*r(b)
     end do

     open(unit=u,file=fname,form=frm,access=acc,action='write', &
          position='rewind',status='replace',err=101,iostat=io)

     head%i(1)=l
     head%i(2)=m
     write(u,err=111,iostat=io) head
     write(u,err=112,iostat=io) l,m
     write(u,err=113,iostat=io) n(l:m),r(l-1:m)
     write(u,err=114,iostat=io) x
     close(u,err=121,iostat=io)
     if(present(info))info=0
     deallocate(x,stat=i)
     if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif
     return

 101 continue
     write(*,*) subnam,': error opening file: ',io
     if(present(info))info=io
     return
 111 continue
     write(*,*) subnam,': error writing header: ',io
     if(present(info))info=io
     return
 112 continue
     write(*,*) subnam,': error writing lm: ',io
     if(present(info))info=io
     return
 113 continue
     write(*,*) subnam,': error writing nr: ',io
     if(present(info))info=io
     return
 114 continue
     write(*,*) subnam,': error writing cores: ',io
     if(present(info))info=io
     return
 121 continue
     write(*,*) subnam,': error closing file: ',io
     if(present(info))info=io
     return
  end subroutine dtt_write_sdv_file


  !> read fortran SVD format
  !!
  !! Adopted from TT-toolbox 'ttio.f90' file
  !!
  function dtt_read_sdv_file(fname,info) result(arg)
  implicit none
  character(*),intent(in)      :: fname
  integer,intent(out),optional :: info
  type(dtt_tensor) :: arg
  !
  character(*),parameter :: subnam='dtt_read_sdv_file'
  type(tthead) :: head
  integer :: io,u,i,j,k,b,sz
  integer(kind=4) :: l,m,n(tt_size),r(0:tt_size)
  real(kind=8),allocatable :: x(:)
  logical :: ex,op
  integer,parameter :: un=51, rec1=1
  character(len=*),parameter :: frm='unformatted', acc='stream'

     if(present(info))info=-11
     inquire (file=fname, exist=ex, opened=op)
     if(.not.ex)then
        write(*,*)subnam,': file not exist: ',fname
        if(present(info))info=-1
        return
     endif
     if(op)then
        write(*,*)subnam,': file is open, trying to close: ',fname
        inquire(file=fname, number=u)
        write(*,*)subnam,': establish unit: ',u
        close(unit=u,status='keep')
        write(*,*)subnam,': closed ok'
     endif

     u=un; op=.true.
     do while(op)
      inquire(unit=u,opened=op)
      if(op)then
         write(*,*)subnam,': unit ',u,' is busy, trying next '
         u= u + 1
      end if
     end do

     open(unit=u,file=fname,form=frm,access=acc, &
          action='read',position='rewind',status='old', &
          err=101,iostat=io)
     read(u,err=111,iostat=io) head

     if(head%txt(1:2).ne.'TT')then
        write(*,*)subnam,': not TT header in file: ',fname
        if(present(info))info=-2
        return
     end if
     if(head%ver(1).ne.sdv_format_ver(1))then
        write(*,*)subnam,': not correct version of TT file: ',head%ver
        if(present(info))info=-3
        return
     end if

     read(u,err=112,iostat=io) l,m
     arg%l=l; arg%m=m
     if(l.lt.0)then
        write(*,*)subnam,': read strange l,m: ',l,m
     end if

     read(u,err=113,iostat=io) n(l:m),r(l-1:m)
     arg= dtt_tensor(n, r_=r, l_=l, m_=m)

     sz= arg% size()
     if(sz.le.0)write(*,*)subnam,': tt structure has invalid size: ',sz
     allocate(x(max(sz,1)),stat=i)
     if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif

     read(u,err=114,iostat=io) x
     sz=0
     do b=l,m
        forall(i=1:r(b-1), j=1:n(b), k=1:r(b)) &
           arg%u(b)%p(i,j,k)= x(sz + i + (j-1)*r(b-1) + (k-1)*r(b-1)*n(b))
        sz= sz + r(b-1)*n(b)*r(b)
     end do

     close(u,err=121,iostat=io)
     if(present(info))info=0
     deallocate(x,stat=i)
     if(i.ne.0)then;write(*,*)subnam,': cannot allocate core array';stop;endif
     return

 101 continue
     write(*,*) subnam,': error opening file: ',io
     if(present(info))info=io
     return
 111 continue
     write(*,*) subnam,': error reading header: ',io
     if(present(info))info=io
     return
 112 continue
     write(*,*) subnam,': error reading lm: ',io
     if(present(info))info=io
     return
 113 continue
     write(*,*) subnam,': error reading nr: ',io
     if(present(info))info=io
     return
 114 continue
     write(*,*) subnam,': error writing cores: ',io
     if(present(info))info=io
     return
 121 continue
     write(*,*) subnam,': error closing file: ',io
     if(present(info))info=io
     return
  end function dtt_read_sdv_file
!
! W.I.P
!
end module thorio_lib
