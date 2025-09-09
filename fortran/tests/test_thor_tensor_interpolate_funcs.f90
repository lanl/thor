program main
character(len=*), parameter :: fname_in = 'tests/testdata/test01_small.dat'

  call test_tensor_interpolate_funcs(fname_in)
  !call test_bug(fname_in)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_tensor_interpolate_funcs(fname_input_array)
  use thor_lib !TODO: CLEANUP!!!, only: dtt_tensor,dtt_tensor_ones,dtt_tensor_ones_like,sayfull
  use ttaux_lib
  !use dmrgg_lib
  use time_lib
  !use quad_lib
  use default_lib, ONLY: readarg
  !use mat_lib
  use omp_lib
  use thorio_lib, only : read_array_from_file
  use interpolate_special_funcs, only : &
      inv_tt, div_tt, max_tt, sqrt_tt, merge_tt
  use interpolate_special_funcs, only : &
      set_val_where_tt_cond_else, set_func_where_tt_cond_else
  use ttop_lib, only : del_tt !_dtt
  implicit none
  character(*), intent(in) :: fname_input_array
  include 'mpif.h'
  type(dtt_tensor) :: tt,qq,var_tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10
  type(dtt_tensor) :: mtt(3)
  integer :: i,j,p,m,n,nx,r,piv,decay,info,nproc,me,adj
  integer :: resx,resy,resz,res_total, resMax
  integer(kind=8) :: neval
  double precision :: f,bnd,t1,t2,tcrs,einf,efro,ainf,afro
  double precision :: acc,val,tru,h,w,t, eps
  double precision,allocatable   :: par(:),mat(:,:),crtt(:),arr(:,:,:)
  character(len=1)               :: a
  character(len=32)              :: aa
  character(len=2)               :: cond
  logical                        :: rescale
  integer                        :: var_d
  integer, allocatable           :: var_n(:), var_r(:)
  double precision,external      :: dfunc_inv_tt, dfunc_tt_gt !,dfunc_ising_discr
  double precision,external      :: func1, func1_default
  double precision               :: val0, inv_val0, dh
  integer                        :: ind(4), del_dir, ngc, failed
  integer, allocatable           :: ps(:)

     print '("")'
     print '("Test 05: DMRG interpolation functions")'
     print '("")'
     failed= 0

     call mpi_init(info)
     if(info.ne.0) error stop '[MPI]: init fail'

     ! Read params (optional)
     call readarg(1,a,'c')      ! type of the integral, e.g. 'c', 'd' or 'e'
     call readarg(2,m,2)        ! index of the integral, e.g 6 to calculate C_6
     call readarg(3,n,4)        ! quadrature mode size
     call readarg(4,r,4)        ! max TT rank
     call readarg(5,piv,1)      ! pivoting strategy for cross interpolation

!     !$OMP PARALLEL
!     if (omp_get_thread_num().eq.0) &
!        write(*,'(3x,a,i10)')'OMP thrds:', omp_get_num_threads()
!     !$OMP END PARALLEL
     acc=500*epsilon(1.d0)
 
     !call read_array_from_file('bench01_array.dat', resx,resy,resz,arr)
     !call read_array_from_file('bench02_small.dat', resx,resy,resz,arr)
     print '("-- [1-1] Read input array with read_array_from_file(..):")'
     !call read_array_from_file(fname_input_array, resx,resy,resz,arr)
     call read_array_from_file('tests/testdata/test01_small.dat', resx,resy,resz,arr)
     print '("")'
     print '("   Input Matrix Read OK:")'
     print '("   Input Matrix:")'
     print '("   resx= ",I6,"| resy= ",I6,"| resz= ",I6)', resx,resy,resz
     res_total = resx*resy*resz
     resMax = max(max(resx,resy),resz)
     print '("    --> total resolution = ",I12)', res_total
     allocate(crtt(res_total), stat=info)
     crtt=reshape(arr,(/res_total/))

     print '("   Define tt matrix parameters:")'
     eps = 1.0e-12
     var_n = (/resx,resy,resz/)
     var_d = size(var_n)! dimension
     allocate(var_r(0:var_d))
     var_r = 128; var_r(0) = 1; var_r(var_d) = 1
     !print '("    - m= ",I4, "n= ",I4)',m,n
     print '("    - dim(var_tt) = ",I5)', var_d
     print '("    - n(var_tt)   = ",12(I5,1X))', var_n
     print '("    - eps         = ",ES12.5)', eps
     print '("")'

     print '("-- [2] Conversion from array to var_tt in TT-format:")'
     t1= timef()
     var_tt% l = 1
     var_tt% r = var_r
     var_tt = dtt_tensor(crtt,var_n,eps_=eps,rmax_=resx)
     t2= timef()
     tcrs= t2 - t1
     print '("   completed in ",ES12.5," sec.")', tcrs
     call var_tt% say()
     call sayfull(var_tt)
     print '("")'

     allocate(par(2*resMax+1))
     par(2) = 0.d50

     tt = dtt_tensor_ones(var_n)
     tt2= dtt_tensor_ones_like(tt)
     tt3= dtt_tensor_ones_like(tt)
     tt4= dtt_tensor_ones_like(tt)
     tt5= dtt_tensor_ones_like(tt)
     tt6= dtt_tensor_ones_like(tt)
     tt7= dtt_tensor_ones_like(tt)
     tt8= dtt_tensor_ones_like(tt)
     tt9= dtt_tensor_ones_like(tt)

     print '("-- [3] Initial guess for 1/var_tt:")'
     call tt% say()
     call sayfull(tt)
     print '("")'

     print '("-- [4] Calculating 1/var_tt using DMRG:")'
     t1= timef()
     par(2) = 1.0e-15   ! eps for DMRG approx
     !call inv_tt(var_tt, tt, par, accuracy=acc, maxrank=r)
     !call inv_tt(var_tt, tt, par, accuracy=acc)
     !call inv_tt(var_tt, tt, par)
     call inv_tt(var_tt, tt)
     t2= timef()
     tcrs= t2-t1
     print '("   completed in ",ES12.5," sec.")', tcrs
     print '("")'

     print '("-- [5] Final result for computing 1/var_tt:")'
     call tt% say()
     call sayfull(tt)
     print '("")'

     print '("-- [6] Approximating tt2 = max(tt, 1.0):")'
     t1=timef()
     par(2) = 1.0e-15   ! eps for DMRG approx
     par(3) = 1.d0      ! This is how you set the max treshold
     call max_tt(tt, tt2, par)
     t2=timef()
     tcrs=t2-t1
     print '("   completed in ",ES12.5," sec.")', tcrs
     call tt2% say()
     call sayfull(tt2)
     print '("")'

     print '("-- [7] Approximating tt3 = sqrt(tt2):")'
     t1=timef()
     par(2) = 1.0e-15   ! eps for DMRG approx
     call sqrt_tt(tt2, tt3, par)
     t2=timef()
     tcrs=t2-t1
     print '("   completed in ",ES12.5," sec.")', tcrs
     call tt3% say()
     call sayfull(tt3)
     print '("")'

!!     print '("-- [8-1] Testing div_tt: tt4 = var_tt/tt")'
!!     par(2) = 1.0e-15   ! eps for DMRG approx
!!     call div_tt(var_tt, tt, tt4, par)
!!     call say(tt4)
!!     call sayfull(tt4)
!!
!!     print '("-- [8-1] Testing div_tt: tt5 = var_tt/tt3")'
!!     par(2) = 1.0e-15   ! eps for DMRG approx
!!     call div_tt(var_tt, tt3, tt5, par)
!!     call say(tt5)
!!     call sayfull(tt5)
 
     ! Test Conditional branch set_val_where_tt_cond_else()
       print*, "[9] Testing set_val_where_tt_cond_else(): tt6 = set_val_where_tt_cond_else(tt5>=0,set 1, else set 0):"
       t1=timef()
       par(3)=0.0     ! target
       par(4)=1.d0  ! where value
       par(5)=1.0     ! else_swhich
       par(6)=0.d0 ! elsewhere value
       cond  = "GT"
       call set_val_where_tt_cond_else(tt, tt6, cond, par) !< FAILS HERE
       t2=timef()
       tcrs=t2-t1
       if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '... [',9,'] completed in ',tcrs,' sec.'
       call tt6% say()
       call sayfull(tt6)
    ! Test Conditional branch set_func_where_tt_cond_else()
      print*, "[10] Testing set_func_where_tt_cond_else(): tt7 = set_func_where_tt_cond_else(tt6 < 0.0, set func1, else set func1_default):"
      t1=timef()
      par(3)=0.0     ! target
      par(5)=1.0     ! else_swhich
      cond  = "LT"
      call set_func_where_tt_cond_else(tt6, tt7, cond, func1, func1_default, par)
      t2=timef()
      tcrs=t2-t1
      if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '... [',10,'] completed in ',tcrs,' sec.'
      call tt7% say()
      call sayfull(tt7)
    ! Test merge_tt(tt5,tt6,mask=tt7)
      print*, "[11] Testing merge_tt(): tt8 = merge(tt, tt7, mask=tt6):"
      t1=timef()
      call merge_tt(tt, tt7, tt6, tt8, par)
      t2=timef()
      tcrs=t2-t1
      if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '... [',11,'] completed in ',tcrs,' sec.'
      call tt8% say()
      call sayfull(tt8)
    ! Test del_tt(tt8,dir=1,dx=0.1)
      print*, "[12] Testing del_tt(): tt9 = del(tt8, dir=1, h=0.1):"
      t1=timef()
      del_dir = 1; dh=0.1; ngc=1
      tt9 = del_tt(tt8, del_dir, dh, ngc)
      t2=timef()
      tcrs=t2-t1
      if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '... [',12,'] completed in ',tcrs,' sec.'
      !tt9 = del_dtt(tt8, del_dir, dh, ngc)
      !call del_dtt(tt8, tt9, del_dir, dh, ngc)
      call tt9% say()
      call sayfull(tt9)
    ! Test operator overloading "*"
      print*, "[13] tt1*tt2: tt10 = tt8*tt9"
      t1=timef()
      tt10 = tt9 * tt8
      t2=timef()
      tcrs=t2-t1
      if(me.eq.0)write(*,'(a,i12,a,e12.4,a)') '... [',13,'] completed in ',tcrs,' sec.'
      call tt10% say()
      call sayfull(tt10)
    ! Test multidimentional tt mtt = tt(3)
      print*, "[14] testing multidim tt: mtt"
      mtt(1) = tt5 + tt6
      mtt(2) = tt6 + tt7
      mtt(3) = tt7 + tt8

      print*, "mtt(1) = tt5 + tt6:"
      call mtt(1)% say()
      !call sayfull(mtt(1))
      print*, "mtt(2) = tt6 + tt7:"
      call mtt(2)% say()
      !call sayfull(mtt(2))
      print*, "mtt(3) = tt7 + tt8:"
      call mtt(3)% say()
      !call sayfull(mtt(3))




     ! Deallocate all tt vectors
      deallocate(arr)
      call dealloc(mtt(1)); call dealloc(mtt(2)); call dealloc(mtt(3))
      print*, "[10] TT-vect tt, tt2, tt3, tt4, tt5 and var_tt deallocated OK:"
      call mpi_finalize(info)
      if(info.ne.0)then;write(*,*)'mpi: finalize fail: ',info;stop;endif
  end subroutine test_tensor_interpolate_funcs

end program


    !>
    !! ARGUMENTS:
    !!  m      - index of most right mode of tt (=dimension of tt)
    !!  ind    - indeces of the diferent modes at which tt is being evaluated
    !!  n      - modes of tt
    !!  tt     - tt-vector format
    !!  par    - auxiliary array where to hold parameters useful for the function call
    !!
    double precision function func1(m, ind, n, tt, par) result(f)
    use thor_lib
    implicit none
    type(dtt_tensor), intent(inout)         :: tt
    integer,intent(in)                      :: m
    integer,intent(in)                      :: ind(m),n(m)
    double precision,intent(inout),optional :: par(*)
    double precision                        :: D, eps, zeps, zero, sol

       D    = tijk(tt, ind)
       f = SQRT(7d0 + 5d0/ MAX(D, 1d0))
       if (sol.lt.1d0) f = 9d0

    end function func1


    !>
    !! ARGUMENTS:
    !!  m      - index of most right mode of tt (=dimension of tt)
    !!  ind    - indeces of the diferent modes at which tt is being evaluated
    !!  n      - modes of tt
    !!  tt     - tt-vector format
    !!  par    - auxiliary array where to hold parameters useful for the function call
    !!
    double precision function func1_default(m, ind, n, tt, par) result(f)
    use thor_lib! only                        : tt_tensor, tijk
    implicit none
    type(dtt_tensor), intent(inout)                :: tt
    integer,intent(in)                      :: m
    integer,intent(in)                      :: ind(m),n(m)
    double precision,intent(inout),optional :: par(*)
    double precision                        :: D !, Tmp0, C, eps, zeps, zero

       f = 1.d0 !0.d0 !SQRT( MAX(f ,zero) )
    
    end function func1_default




