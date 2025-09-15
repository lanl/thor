program main
 use thor_lib
 use ttaux_lib
 use dmrgg_lib
 use time_lib
 !use quad_lib
 !use default_lib
 !use mat_lib
 use omp_lib
 use tt_array_tools_lib, only : darray_to_tt
 implicit none
 include 'mpif.h'
 type(tt_tensor) :: tt,qq, var_tt, tt2
 integer :: i,j,p,m,n,nx,r,piv,decay,info,nproc,me,adj
 integer(kind=8) :: neval
 double precision :: f,bnd,t1,t2,tcrs, einf,efro,ainf,afro, acc,val,tru,h,w,t, eps
 double precision,allocatable::par(:),mat(:,:)
 character(len=1) :: a
 character(len=32) :: aa
 logical :: rescale
 integer :: var_d
 integer, allocatable :: var_n(:), var_r(:)
 double precision,external :: dfunc_inv_tt, dfunc_tt_gt !,dfunc_ising_discr
 double precision          :: val0, inv_val0
 integer                   :: ind(4)
 double precision, allocatable  ::  crtt(:)
 integer, allocatable           ::  ps(:)
 ! Read params
 call readarg(1,a,'c')      ! type of the integral, e.g. 'c', 'd' or 'e'
 call readarg(2,m,2)        ! index of the integral, e.g 6 to calculate C_6
 call readarg(3,n,4)       ! quadrature mode size
 call readarg(4,r,4)       ! max TT rank
 call readarg(5,piv,1)      ! pivoting strategy for cross interpolation


  print*, "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  print*, "|               Test 01: Testing DMRG for approximation of 1/tt             |"
  print*, "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  print*,""


 call mpi_init(info)
 if(info.ne.0)then;write(*,*)'mpi: init fail: ',info;stop;endif
 call mpi_comm_size(MPI_COMM_WORLD,nproc,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_size fail: ',info;stop;endif
 call mpi_comm_rank(MPI_COMM_WORLD,me,info)
 if(info.ne.0)then;write(*,*)'mpi: comm_rank fail: ',info;stop;endif

!$OMP PARALLEL
  if (omp_get_thread_num().eq.0) then
   write(*,'(3x,a,i10)')'OMP thrds:', omp_get_num_threads()
  end if 
!$OMP END PARALLEL
 acc=500*epsilon(1.d0)

 allocate(par(2*n+1), stat=info)
 if(info.ne.0)then;write(*,*)'cannot allocate par:',info;stop;endif
 
 !t1=timef()
 tt%l=1;tt%m=m;tt%n=n;tt%r=r;call alloc(tt)
 tt2%l=1;tt2%m=m;tt2%n=n;tt2%r=r;call alloc(tt2)
 call ones(tt)
 call ones(tt2)

 ! integer :: var_d
 ! integer, allocatable :: var_n(:), var_r(:), 

 ! Define tt metadata
    var_d  = m         ! dimension
    allocate(var_n(var_d), var_r(var_d+1), stat=info)
    if(info.ne.0)then;write(*,*)'[!] Problem allocating var_tt meta var_n, var_r';stop;endif
    var_n=n; var_r(1)=1; var_r(2:var_d)=r; var_r(var_d+1)=1;
    par(2) = 0.d50

  ! [!] dtt_rnd(tt) needs some meta data of tt initilized
  ! Set minimum tt meta data to make the dtt_rnd function call
    var_tt%l = 1      ! index of the leftmost core
    var_tt%m = var_d  ! index of the rightmost core
    var_tt%n = var_n  ! mode sizes (storage for matrices)
    var_tt%r = var_r  ! TT ranks
  ! [!] Note that the initialion of the cores var_tt%u, neither of var_tt%q,
  !          var_tt%s, or var_tt%t, are required by dtt_rnd(var_tt)

    print*,""
    print*, "[+] m=",m, "n=",n
    print*,"[+] dim(var_tt) = ", var_d
    print*,"[+] n(var_tt)   = ", var_n


  ! Creat random tt
  allocate(crtt(16), ps(var_d+1),stat=info)
  if(info.ne.0)then;write(*,*)'[!] Problem allocating crtt';stop;endif
  crtt = (/0.8641, 0.9178 ,0.8869,0.8460,0.9757,0.7043,1.2966,1.5632,1.0847,0.8228,1.3630,1.5189, 0.0,0.0,0.0,0.0/)


  print*,""
  print*,"[1] Input Matrix var =", crtt

  eps = 1.0e-16
  call darray_to_tt(crtt,var_n,var_d,eps,var_tt)
  !call svd(var_n, crtt, var_tt, eps)

  !ps(1)=1;
  !do i=1,m
  !  ps(i+1) = ps(i) + var_r(i)*var_n(i)*var_r(i+1)
  !end do
  !call alloc(var_tt)
  !do i=1,m
  !     call dcopy(var_r(i)*var_n(i)*var_r(i+1),crtt(ps(i)),1,var_tt%u(i)%p,1)
  !end do
  
  !  call dtt_rnd(var_tt)
  ! call ones(var_tt)

 print*,""
 print*, "[2] --> in TT-format  var_tt:"
 call say(var_tt)
 call sayfull(var_tt)

 print*,""
 print*, "[3] Initial guess for 1/var tt:"
 call say(tt)
 call sayfull(tt)

 print*,""
 print*, "[4] Running DMRG on 1/tt_var:"
 t1=timef()
 !call dtt_dmrgg(tt,dfunc_ising_discr,par,maxrank=r,accuracy=acc,pivoting=piv,neval=neval,quad=qq)
 !call dtt_dmrgg_tt(tt, dfunc_inv_tt, var_tt, par, maxrank=r,accuracy=acc,pivoting=piv,neval=neval,quad=qq)
 !call dtt_dmrgg_tt(tt, dfunc_inv_tt, var_tt, par, maxrank=r,accuracy=acc,pivoting=piv)
 !call dtt_dmrgg_tt(tt, dfunc_inv_tt, var_tt, par, maxrank=r,accuracy=acc)
 call dtt_dmrgg_tt(tt, dfunc_inv_tt, var_tt, par, accuracy=acc)
 !call dtt_dmrgg_tt(tt, dfunc_tt_gt, var_tt, par, maxrank=r,accuracy=acc)
 !call dtt_dmrgg_tt(tt2, dfunc_inv_tt, tt, par, maxrank=r,accuracy=acc)
 t2=timef()
 tcrs=t2-t1
 if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '...with',neval,' evaluations completed in ',tcrs,' sec.'

 print*,""
 print*, "[5] Final guess for 1/var tt:"
 call say(tt)
 call sayfull(tt)

 !print*, "[2] tt2:"
 !call say(tt2)
 !call sayfull(tt2)

 call dealloc(tt)
 call dealloc(tt2)
 call dealloc(var_tt)
 call mpi_finalize(info)
 if(info.ne.0)then;write(*,*)'mpi: finalize fail: ',info;stop;endif
end program




!
! INVERSE 1/tt
!


double precision function dfunc_inv_tt(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used dtt_dmrg_tt to interpolate 1/var_tt
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( abs(val) .le. 1.0E-12)  then
     f = 0.d0
 else
     f = 1.d0/val
 endif
end function





!
! CONDITIONALS tt (.eq., .ge., .gt., .le., .lt.) par(2)
!




double precision function dfunc_tt_eq(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt == par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( abs(val - par(2)) .le. 1.0E-12)  then
     f = 1.d0
 else
     f = 0.d0
 endif
end function


double precision function dfunc_tt_gt(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt > par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 !if ( val .gt. par(2))  then
 if ( val .gt. 0.50)  then
     f = 1.d0
 else
     f = 0.d0
 endif
end function


double precision function dfunc_tt_ge(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt >= par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( val .ge. par(2))  then
     f = 1.d0
 else
     f = 0.d0
 endif
end function


double precision function dfunc_tt_lt(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt < par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( val .lt. par(2))  then
     f = 1.d0
 else
     f = 0.d0
 endif
end function



double precision function dfunc_tt_le(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt <= par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( val .le. par(2))  then
     f = 1.d0
 else
     f = 0.d0
 endif
end function






!
! WHERE  where( tt {.eq., .ge., .gt., .le., .lt.} par(2) )
!



double precision function dfunc_where_tt_eq(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt == par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( abs(val - par(2)) .le. 1.0E-12)  then
     f = dble(product(ind))
 else
     f = -99.d0
 endif
end function


double precision function dfunc_where_tt_gt(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt > par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( val .gt. par(2))  then
     f = dble(product(ind))
 else
     f = -99.d0
 endif
end function


double precision function dfunc_where_tt_ge(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt >= par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( val .ge. par(2))  then
     f = dble(product(ind))
 else
     f = -99.d0
 endif
end function


double precision function dfunc_where_tt_lt(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt < par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( val .lt. par(2))  then
     f = dble(product(ind))
 else
     f = -99.d0
 endif
end function



double precision function dfunc_where_tt_le(m, ind, n, var_tt, par) result(f)
 !
 ! This function is used by dtt_dmrg_tt to check indeces at which var_tt <= par(2)
 ! ARGUMENTS:
 !  m      - index of most right mode of var_tt (=dimension of var_tt)
 !  ind    - indeces of the diferent modes at which var_tt is being evaluated
 !  n      - modes of var_tt
 !  var_tt - tt-vector being interpolated
 !  par    - auxiliary array where to hold parameters useful for the function call
 !
 use thor_lib, only                        : tt_tensor, tijk
 implicit none
 type(tt_tensor), intent(inout)                :: var_tt
 integer,intent(in)                      :: m
 integer,intent(in)                      :: ind(m),n(m)
 double precision,intent(inout),optional :: par(*)
 double precision                        :: val
 val=tijk(var_tt, ind)
 !print*,"[+] val(",ind,") = ", val
 if ( val .le. par(2))  then
     f = dble(product(ind))
 else
     f = -99.d0
 endif
end function

