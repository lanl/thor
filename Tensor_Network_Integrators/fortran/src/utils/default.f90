module default_lib
 implicit none
 interface default
  module procedure default_i, default_d, default_z, default_a
 end interface
 interface readarg
  module procedure readarg_i, readarg_d, readarg_a
 end interface
contains
!
! DEFAULT
!
 integer function default_i(def,opt) result(fun)
  integer,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 double precision function default_d(def,opt) result(fun)
  double precision,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 double complex function default_z(def,opt) result(fun)
  double complex,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
 character function default_a(def,opt) result(fun)
  character,intent(in) :: def,opt
  optional :: opt
  fun=def
  if(present(opt))fun=opt
 end function
!
! READARG
!
 subroutine readarg_i(pos,arg,opt)
  integer,intent(in) :: pos
  integer,intent(out) :: arg
  integer,intent(in),optional :: opt
  character(len=128) :: str
  if(pos.le.0)then;write(*,*)'readarg_i: illegal pos:',pos;stop;endif
  call getarg(pos,str)
  if(present(opt).and. str.eq.' ')then
    arg=opt
  else
   read(str,*)arg
  end if 
 end subroutine
 subroutine readarg_d(pos,arg,opt)
  integer,intent(in) :: pos
  double precision,intent(out) :: arg
  double precision,intent(in),optional :: opt
  character(len=128) :: str
  if(pos.le.0)then;write(*,*)'readarg_d: illegal pos:',pos;stop;endif
  call getarg(pos,str)
  if(present(opt).and. str.eq.' ')then
    arg=opt
  else
   read(str,*)arg
  end if 
 end subroutine
 subroutine readarg_a(pos,arg,opt)
  integer,intent(in) :: pos
  character,intent(out) :: arg
  character,intent(in),optional :: opt
  character(len=128) :: str
  if(pos.le.0)then;write(*,*)'readarg_a: illegal pos:',pos;stop;endif
  call getarg(pos,str)
  if(present(opt).and. str.eq.' ')then
    arg=opt
  else
   read(str,'(a1)')arg
  end if 
 end subroutine
 
end module
