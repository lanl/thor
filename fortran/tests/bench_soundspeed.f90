program main
 use thor_lib
 use ttaux_lib
 !use dmrgg_lib
 use time_lib
 !use quad_lib
 !use default_lib
 !use mat_lib
 use omp_lib

 use pagosa_FLIP_Module, only         : Use_2nd_Order_Q
 use pagosa_FLIP_util_Module, only    : DoPlus
 use pagosa_soundspeed_lib, only      : Sound_Speed, Original_Sound_Speed
 use pagosa_Options_Module, only      : cutd, zeps

 use tt_array_tools_lib, only : darray_to_tt, dtt_to_array, read_array_from_file
 use interpolate_specfial_funcs, only : inv_tt, div_tt, max_tt, sqrt_tt, merge_tt
 !use interpolate_specfial_funcs, only : set_val_where_tt_cond_else, set_func_where_tt_cond_else
 use ttop_lib, only : del_tt !_dtt
 implicit none
 include 'mpif.h'
 type(tt_tensor) :: tt,qq, var_tt, tt2 !, tt3, tt4, tt5, tt6, tt7, tt8, tt9, tt10
 type(tt_tensor) :: mtt(3)
 integer :: i,j,p,m,n,nx,r,piv,decay,info,nproc,me,adj, resx,resy,resz,res_total, resMax
 integer(kind=8) :: neval
 double precision :: f,bnd,t1,t2,tcrs, einf,efro,ainf,afro, acc,val,tru,h,w,t, eps
 double precision,allocatable   :: par(:),mat(:,:)
 character(len=1)               :: a
 character(len=32)              :: aa
 character(len=2)               :: cond
 logical                        :: rescale
 integer                        :: var_d
 integer, allocatable           :: var_n(:), var_r(:)
 double precision,external      :: dfunc_inv_tt, dfunc_tt_gt !,dfunc_ising_discr
 double precision,external      :: func1, func1_default
 double precision               :: val0, inv_val0
 integer                        :: ind(4)
 double precision, allocatable  ::  crtt(:), arr(:,:,:)
 integer, allocatable           ::  ps(:)
 integer                        :: del_dir, ngc
 double precision               :: dh
 ! Read params
 call readarg(1,a,'c')      ! type of the integral, e.g. 'c', 'd' or 'e'
 call readarg(2,m,2)        ! index of the integral, e.g 6 to calculate C_6
 call readarg(3,n,4)       ! quadrature mode size
 call readarg(4,r,4)       ! max TT rank
 call readarg(5,piv,1)      ! pivoting strategy for cross interpolation


  print*, "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  print*, "|            BENCH 01: Benchmarking DMRG Interpolation functions            |"
  print*, "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  print*,""

!$OMP PARALLEL
  if (omp_get_thread_num().eq.0) then
   write(*,'(3x,a,i10)')'OMP thrds:', omp_get_num_threads()
  end if 
!$OMP END PARALLEL
 acc=500*epsilon(1.d0)

 ! Read Input Array
    !call read_array_from_file('bench01_array.dat', resx,resy,resz,arr)
    call read_array_from_file('bench02_small.dat', resx,resy,resz,arr)
    print*,""
    print*,"[1-1] Input Matrix Read OK:"
    print*,"      Input Matrix:"
    print*,"    resx =,",resx,"| resy =",resy,"| resz =",resz
    res_total = resx*resy*resz
    resMax = max(max(resx,resy),resz)
    print*,"    --> total resolution = ", res_total
    allocate(crtt(res_total), stat=info)
    crtt=reshape(arr,(/res_total/))
 ! Define tt metadata
    var_d = 3! dimension
    allocate(var_n(var_d), var_r(var_d+1), stat=info)
    if(info.ne.0)then;write(*,*)'[!] Problem allocating var_tt meta var_n, var_r';stop;endif
    var_n = (/resx,resy,resz/)
    var_r(1)=1; 
    var_r(2:var_d)= 64 !resMax; ! Max TT-rank
    var_r(var_d+1)=1;
  ! [!] dtt_rnd(tt) needs some meta data of tt initilized
  ! Set minimum tt meta data to make the dtt_rnd function call
    var_tt%l = 1      ! index of the leftmost core
    var_tt%m = var_d  ! index of the rightmost core
    var_tt%n = var_n  ! mode sizes (storage for matrices)
    var_tt%r = var_r  ! TT ranks
  ! [!] Note that the initialion of the cores var_tt%u, neither of var_tt%q,
  !          var_tt%s, or var_tt%t, are required by dtt_rnd(var_tt)
    print*,""
    print*,"[+] m=",m, "n=",n
    print*,"[+] dim(var_tt) = ", var_d
    print*,"[+] n(var_tt)   = ", var_n

    allocate(par(2*resMax+1), stat=info)
    if(info.ne.0)then;write(*,*)'cannot allocate par:',info;stop;endif
    !t1=timef()
    tt%l=1  ; tt%m=var_d  ; tt%n=var_n  ;tt%r=var_r ; call alloc(tt)
    tt2%l=1 ; tt2%m=var_d ; tt2%n=var_n ;tt2%r=var_r ; call alloc(tt2)
    !tt3%l=1 ; tt3%m=var_d ; tt3%n=var_n ;tt3%r=var_r ; call alloc(tt3)
    !tt4%l=1 ; tt4%m=var_d ; tt4%n=var_n ;tt4%r=var_r ; call alloc(tt4)
    !tt5%l=1 ; tt5%m=var_d ; tt5%n=var_n ;tt5%r=var_r ; call alloc(tt5)
    !tt6%l=1 ; tt6%m=var_d ; tt6%n=var_n ;tt6%r=var_r ; call alloc(tt6)
    !tt7%l=1 ; tt7%m=var_d ; tt7%n=var_n ;tt7%r=var_r ; call alloc(tt7)
    !tt8%l=1 ; tt8%m=var_d ; tt8%n=var_n ;tt8%r=var_r ; call alloc(tt8)
    !tt9%l=1 ; tt9%m=var_d ; tt9%n=var_n ;tt9%r=var_r ; call alloc(tt9)
    call ones(tt);  call ones(tt2); !call ones(tt3); call ones(tt4); call ones(tt5)
    !call ones(tt6); call ones(tt7); call ones(tt8); !call ones(tt9)
    t1=timef()
  ! Convert into tt Vector format
    par(2) = 0.d50
    eps = 1.0e-12
    call darray_to_tt(crtt,var_n,var_d,eps,var_tt, rmax=10)
    t2=timef()
    tcrs=t2-t1
    print*,""
    print*, "[2] --> in TT-format  var_tt:"
    if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '[',2,'] completed in ',tcrs,' sec.'
    call say(var_tt)
    call sayfull(var_tt)
  ! Show Initial Guess for 1/var_tt
    print*,""
    print*, "[3] Initial guess for 1/var tt:"
    call say(tt)
    call sayfull(tt)
  ! Approximate 1/var_tt with DMRGG
    print*,""
    print*, "[4] Running DMRGG on 1/tt_var:"
    call mpi_init(info)
    if(info.ne.0)then;write(*,*)'mpi: init fail: ',info;stop;endif
    t1=timef()
    par(2) = 1.0e-15   ! eps for DMRG approx
    !call inv_tt(var_tt, tt, par, accuracy=acc, maxrank=r)
    !call inv_tt(var_tt, tt, par, accuracy=acc)
    !call inv_tt(var_tt, tt, par)
    call inv_tt(var_tt, tt)
    t2=timef()
    tcrs=t2-t1
    if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '[',4,'] completed in ',tcrs,' sec.'
  ! Show tt, Axproximation of 1/var_tt
    print*,""
    print*, "[5] Final guess for 1/var tt:"
    call say(tt)
    call sayfull(tt)
    t1=timef()
  ! Approximation of Max(tt, 1.0):
    par(2) = 1.0e-15   ! eps for DMRG approx
    par(3) = 1.d0      ! This is how you set the max treshold
    call max_tt(tt, tt2, par)
    t2=timef()
    tcrs=t2-t1
    print*, "[6] triming Final guess for tt2 = max(1/var,1.0) tt2:"
    if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '... [',6,'] completed in ',tcrs,' sec.'
    call say(tt2)
    call sayfull(tt2)
  ! Deallocate all tt vectors
    call dealloc(tt); call dealloc(tt2) !; call dealloc(tt3); call dealloc(tt4); call dealloc(tt5); call dealloc(var_tt);
    !call dealloc(tt6); call dealloc(tt7); call dealloc(tt8); call dealloc(tt9); call dealloc(tt10)
    !call dealloc(mtt(1)); call dealloc(mtt(2)); call dealloc(mtt(3))
    print*, "[10] TT-vect tt, tt2, tt3, tt4, tt5 and var_tt deallocated OK:"
    call mpi_finalize(info)
    if(info.ne.0)then;write(*,*)'mpi: finalize fail: ',info;stop;endif
end program


    double precision function func1(m, ind, n, tt, par) result(f)
      !
      ! ARGUMENTS:
      !  m      - index of most right mode of tt (=dimension of tt)
      !  ind    - indeces of the diferent modes at which tt is being evaluated
      !  n      - modes of tt
      !  tt     - tt-vector format
      !  par    - auxiliary array where to hold parameters useful for the function call
      !
      use thor_lib! only                        : tt_tensor, tijk
      implicit none
      type(tt_tensor), intent(inout)                :: tt
      integer,intent(in)                      :: m
      integer,intent(in)                      :: ind(m),n(m)
      double precision,intent(inout),optional :: par(*)
      double precision                        :: D, eps, zeps, zero, sol
      D    = tijk(tt, ind)
      if(present(par)) then; eps=par(2); else; eps=1.0E-15; endif
      zeps = 1.d0 !par(3)
      zero = 0.d0
      !f = SQRT( MAX(C + Tmp0/ MAX(D,zeps) ,zero) )
      sol  = 1.0*D !Max(D,0.0)
      sol  = MAX(sol,0.0)
      sol  = SQRT( MAX(7.0 + 5.0/ MAX(sol,zeps) ,zero) )
      if(sol.ge.1.d0) then
         f = sol
      else
         f = 9.d0
      endif
      !f = SQRT( MAX(f ,zero) )
      !f = 9.0 !D !sol
    end function func1
    double precision function func1_default(m, ind, n, tt, par) result(f)
      !
      ! ARGUMENTS:
      !  m      - index of most right mode of tt (=dimension of tt)
      !  ind    - indeces of the diferent modes at which tt is being evaluated
      !  n      - modes of tt
      !  tt     - tt-vector format
      !  par    - auxiliary array where to hold parameters useful for the function call
      !
      use thor_lib !,! only                        : tt_tensor, tijk
      implicit none
      type(tt_tensor), intent(inout)                :: tt
      integer,intent(in)                      :: m
      integer,intent(in)                      :: ind(m),n(m)
      double precision,intent(inout),optional :: par(*)
      double precision                        :: D !, Tmp0, C, eps, zeps, zero
      !D    = tijk(tt, ind)
      !if(present(par)) then; eps=par(2); else; eps=1.0E-15; endif
      !zeps = 1.d0 !par(3)
      !zero = 0.d0
      !!f = SQRT( MAX(C + Tmp0/ MAX(D,zeps) ,zero) )
      !f = MAX(D,zeps)
      !if(f.ge.10) then
      !   f = D/f
      !else
      !   f = 0.d0
      !endif
      f = 1.d0 !0.d0 !SQRT( MAX(f ,zero) )
    end function func1_default




