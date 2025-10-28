!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file thorio.f90
!! @author Oleg Korobkin, korobkin@lanl.gov
!! @author Ismael Djibrilla Boureima iboureima@lanl.gov
!! @date   Septober 2025
!! @brief  THOR input/output routines
!!
module dthorio_lib
use distributed_comms,  only: super_comm

contains

  !> Writes a distributed THOR tensor into a file in ASCII format
  !!
  !! Data format:
  !! - <file_prefix>_<rank#>.dat: one file per rank;
  !! - first row: tensor dimensions in all modes, separated by commas;
  !! - second row: tensor ranks in all modes, separated by commas
  !! - third row: nzb-s, blocking dimensions in the z-direction, s.b.c.;
  !! - fourth row: nzl-s, local z-sizes, separated by commas;
  !! - rest of the file: localor cores, flattened, in sequential order,
  !!   single value per line
  !!
  !! File example:
  !! >>>>>>>>>>>>>>>> BEGIN FILE Q_0.dat out ot 3 files >>>>>>>>>>>>>>
  !! 120,125,128,131
  !! 1,12,14,12,1
  !! 5,5,5,5     <- blocking dimensions for each core
  !! 40,45,45,45 <- local z-size for each core on rank 0
  !! -0.132E+00
  !!  1.333E+03  |
  !!  1.333E+03   \  total of 15645 =
  !!  1.333E+03   /   = 12*40 + 12*14*45 + 14*12*45 + 1*45 lines
  !!  ...        |
  !! -1.000E-0.1 /
  !! >>>>>>>>>>>>>>>> BEGIN FILE Q_1.dat out ot 3 files >>>>>>>>>>>>>>
  !! 120,125,128,131
  !! 1,12,14,12,1
  !! 5,5,5,5
  !! 40,40,43,45 <- local z-size for each core on rank 1
  !! -0.132E+00  \
  !!  1.333E+03  |
  !!  1.333E+03   \  total of 14469 =
  !!  1.333E+03   /   = 12*40 + 12*14*40 + 14*12*43 + 1*45 lines
  !!  ...        |
  !! -1.000E-0.1 /
  !! >>>>>>>>>>>>>>>> BEGIN FILE Q_2.dat out ot 3 files >>>>>>>>>>>>>>
  !! 120,125,128,131
  !! 1,12,14,12,1
  !! 5,5,5,5
  !! 40,40,40,41 <- local z-size for each core on rank 2
  !! -0.132E+00  \
  !!  1.333E+03  |
  !!  1.333E+03   \  total of 13961 =
  !!  1.333E+03   /   = 12*40 + 12*14*40 + 14*12*40 + 1*41 lines
  !!  ...        |
  !! -1.000E-0.1 /
  !! <<<<<<<<<<<<<<<< END OF 3 FILES <<<<<<<<<<<<<<<<
  !!
  subroutine dthor_write_ascii_files(DT, fnamprefix, VRBZ_)
  use thor
  use string_lib, only: str
  implicit none
  class(dthor),      intent(IN) :: DT
  character(*),      intent(IN) :: fnamprefix
  integer, optional, intent(IN) :: VRBZ_
  !
  character(*), parameter :: subnam = '[dthor_write_ascii_files]'
  integer :: i, sz, nlines, k, d, rank, nranks, VRBZ
  integer, parameter :: fd = 778
  integer, allocatable :: n(:), r(:)
  double precision, pointer :: p1(:)
  double precision, allocatable :: h_A(:)
  character(20) :: fmtstr
  character(200) :: fname

     VRBZ = 0; if (present(VRBZ_)) VRBZ = VRBZ_
     d = DT% d
     if (d < 1) then
        print '("[w]'//subnam//': tensor is empty")'
        return
     endif

     rank = DT% core(1)% rank
     nranks = DT% core(1)% nranks
     if (nranks > 1) then
         write(fname,'(A)') fnamprefix//"_"//str(rank)//".dat"
     else
         write(fname,'(A)') fnamprefix//".dat"
     endif

     if (VRBZ) then
        print '("writing thor object to a file: ",A)', fname
     endif
     open (fd, file=fname, status='replace')

     ! write the first line
     if (DT% is_hstack()) then
        write (fd, '("# hstack ['//str(rank)//'/'//str(nranks)//']")')
     else
        write (fd, '("# vstack ['//str(rank)//'/'//str(nranks)//']")')
     endif

     ! write modes
     do k = 1, d - 1
        write (fd, '(a,",")', advance='no') str(DT%n(k))
     enddo
     write (fd, '(a)') str(DT%n(d))

     ! write ranks
     do k = 0, d - 1
        write (fd, '(a,",")', advance='no') str(DT%r(k))
     enddo
     write (fd, '(a)') str(DT%r(d))

     ! write blocking dimensions
     do k = 1, d - 1
        write (fd, '(a,",")', advance='no') str(DT%core(k)%nzb)
     enddo
     write (fd, '(a)') str(DT%core(d)%nzb)

     ! write local z-size
     do k = 1, d - 1
        write (fd, '(a,",")', advance='no') str(DT%core(k)%nzl)
     enddo
     write (fd, '(a)') str(DT%core(d)%nzl)

     ! write local cores
     do k=1,d
        sz = (DT% core(k)% ML)*(DT% core(k)% NL)
        allocate(h_A(sz))
        h_A(1:sz) = DT% core(k)% d_A(1:sz)
        do i=1,sz
           write(fd, '(ES22.14)') h_A(i)
        enddo
        deallocate(h_A)
        !write(fd, '("=== finished core '//str(k)//'")')
     enddo

     close(fd)
     print '("File '//trim(fname)//' has been written.")'

     !!DEBUG!! TODO: if (VRBZ) then
     !!DEBUG!! TODO:    nlines = file_nlines(fname)
     !!DEBUG!! TODO:    print '(I12," lines written.")', nlines
     !!DEBUG!! TODO: endif

  end subroutine dthor_write_ascii_files


  !> Reads dthor object from files in ASCII format (1 file per rank)
  !!
  !! Data format: see description to dthor_write_ascii_files
  !!
  function dthor_read_ascii_files(fnamprefix, comm, VRBZ_) result(DT)
  use thor
  use string_lib, only: split, isplit, str
  use distributed_arrays, only: file_nlines
  implicit none
  type(dthor) :: DT
  type(super_comm),  intent(IN), target :: comm
  character(*),      intent(IN) :: fnamprefix
  integer, optional, intent(IN) :: VRBZ_
  !
  character(*), parameter :: subnam = '[dthor_read_ascii_files]'
  integer :: VRBZ
  logical :: is_vstack
  integer :: i, sz, nlines, rank, nranks, d, k
  integer, parameter :: fd = 775
  integer, allocatable :: n(:), r(:), nzb(:), nzl(:)
  double precision, allocatable :: cr(:)
  character(:), allocatable :: line, strarr(:)
  character(200) :: fname, rankstr
  character :: st

     VRBZ = 0; if (present(VRBZ_)) VRBZ = VRBZ_
     rank = comm% ctx_rank
     nranks = comm% ctx_size
     if (VRBZ>0) &
        print '("reading dthor object with prefix: ",A)', trim(fnamprefix)
     if (nranks > 1) then
         write(fname,'(A)') fnamprefix//"_"//str(rank)//".dat"
     else
         write(fname,'(A)') fnamprefix//".dat"
     endif

     nlines = file_nlines(fname)
     open (fd, file=fname, status='old')
     allocate (character(len=2048*8) :: line)
     read (fd, '(A)') line
     call split(strarr, line, cdel_=' ')
     if (strarr(2).eq.'vstack') then
        st = 'v'
     else if (strarr(2).eq.'hstack') then
        st = 'h'
     else
        error stop '[!]'//subnam//': unknown stacking '//strarr(2)
     endif

     rankstr = '['//str(rank)//'/'//str(nranks)//']'
     if (strarr(3).ne.trim(rankstr)) then
        error stop '[!]'//subnam//': wrong number of ranks'
     endif

     read (fd, '(A)') line; n = isplit(line, cdel_=',')
     read (fd, '(A)') line; r = isplit(line, cdel_=',')
     read (fd, '(A)') line; nzb = isplit(line, cdel_=',')
     read (fd, '(A)') line; nzl = isplit(line, cdel_=',')
     if (VRBZ>2) then
        print '("[i]'//trim(rankstr)//subnam//': n = ",100(I10,1X))', n
        print '("[i]'//trim(rankstr)//subnam//': r = ",100(I10,1X))', r
        print '("[i]'//trim(rankstr)//subnam//': nzb = ",100(I10,1X))', nzb
        print '("[i]'//trim(rankstr)//subnam//': nzl = ",100(I10,1X))', nzl
     endif

     DT = empty_dthor(comm, n, stacking_=st, r_=r)
     d = size(n)
     if (VRBZ>2) then
        print '("[i]'//trim(rankstr)//subnam//': nzb = ",100(I10,1X))', &
              DT% core(1:d)% nzb
        print '("[i]'//trim(rankstr)//subnam//': nzl = ",100(I10,1X))', &
              DT% core(1:d)% nzl
     endif

     sz = 0
     do k = 1, d
        sz = sz + r(k)*nzl(k)*r(k+1)
     enddo
     if (sz.ne.nlines - 5) then
        print '("ERROR: Cannot read the file ",A)', fname
        print '(A22,I12,A6)', "expected file size: ", sz + 2, " lines"
        print '(A22,I12,A6)', "actual file size: ", nlines, " lines"
        error stop
     endif

     nlines = 5
     do k = 1, d
        sz = r(k)*nzl(k)*r(k+1)
        allocate(cr(sz))
        do i=1,sz
           read(fd, *) cr(i)
           nlines = nlines + 1
        enddo
        DT% core(k)% d_A = cr
        deallocate(cr)
     enddo
     close (fd)

     deallocate(n,r,nzl,nzb,strarr)
     if (VRBZ) then
        print '("[i]'//trim(rankstr)//subnam//': lines read = ",I10)', nlines
        if (rank.eq.0) then
           print '("[i]'//subnam//', dthor object: ")'
           call DT% info()
        endif
     endif

  end function dthor_read_ascii_files


end module dthorio_lib

