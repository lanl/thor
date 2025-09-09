!--------------------------------------------------------------------------~*
!! Copyright (c) 2024 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file test_qtt_tensor_ops.f90
!! @author [IDB], [OGK]
!! @date  July 2024
!! @brief Testing QTT tensor operations
!!

program main
use rnd_lib, only: random, init_random_seed
implicit none
integer i
character(len=128) arg
!
integer :: err
integer :: d = 3  ! tensor order
integer :: n = 4  ! tensor dimension size
integer :: r = 2  ! tensor rank
integer :: l = 1  ! starting index
integer :: m      ! ending index
integer :: rseed = 0
integer, allocatable :: rr(:)
logical :: lrandomize_nr = .false.
logical :: m_given = .false.
logical :: lsave = .false.
double precision, allocatable :: rnd(:)
double precision :: tol = 1d-12

  m = d
  ! parse command-line arguments
  i = 1; do while (i <= command_argument_count())
     call getarg(i, arg)
     if (trim(arg).eq.'-h') then
        call help_msg
        stop
     elseif (trim(arg).eq.'-d') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) d
        if (err.ne.0) error stop "ERROR: when parsing -d <???>"
        i= i + 1
     elseif (trim(arg).eq.'-l') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) l
        if (err.ne.0) error stop "ERROR: when parsing -l <???>"
        i= i + 1
     elseif (trim(arg).eq.'-m') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) m
        if (err.ne.0) error stop "ERROR: when parsing -m <???>"
        m_given = .true.
        i= i + 1
     elseif (trim(arg).eq.'-t') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) tol
        if (err.ne.0) error stop "ERROR: when parsing -t <???>"
        i= i + 1
     elseif (trim(arg).eq.'--rseed') then
        call getarg(i+1, arg)
        read(arg,*,iostat=err) rseed
        if (err.ne.0) error stop "ERROR: when parsing --rseed <???>"
        i= i + 1
     elseif (trim(arg).eq.'--randomize-r') then
        lrandomize_nr= .true.
     elseif (trim(arg).eq.'--save'.or.trim(arg).eq.'-s') then
        lsave= .true.
     else
        error stop "ERROR: unknown command-line argument '"//arg//"'"
     endif
     i= i + 1
  enddo

  call init_random_seed(rseed)
  allocate(rr(d+1), rnd(d))
  rr(1)= 1; rr(d+1)= 1
  if (lrandomize_nr) then
     call random(rnd)
     do i = 2,d
        rr(i)= 1 + int(rnd(i)*min(2**i,2**(d+1-i)))
     enddo
  else
     rr(2:d)= 2
  endif
  if (.not.m_given) m = d
  call test_qtt_ops(d, rr, tol, lsave)
  deallocate(rr, rnd)


contains

  subroutine help_msg()
  print '(11(/,A))', &
   "Test of tensor operations", &
   "Usage: ./test.exe [-h] [options]", &
   "Options:", &
   " -h : print this help message", &
   " -d <ord>  : tensor order <ord> [3]", &
   " -l <num>  : starting index [1]", &
   " -m <num>  : ending index [d]", &
   " -t <tol>  : use this tolerance in SVD step [1e-12]", &
   " --randomize-r: generate random ranks r(d) [F]", &
   " --rseed <num> : supply random seed [none]", &
   " --save|-s : save generated tensors to /tmp [F]"
  end subroutine

  subroutine test_qtt_ops(d, r, tol, lsave)
  use time_lib
  use mat_lib, only : normfro
  use qtt_lib
  use ttop_lib
  use rnd_lib, only : random
  use matrix_util, only: interlace_3d
  use matlab_struct_module, only : array3d, pprint_matrix
  use ttamen_lib
  implicit none
  integer, intent(in) :: d, r(:)
  double precision, intent(in) :: tol
  logical, intent(in) :: lsave
  !
  integer, parameter  :: dmax4print = 12
  include 'mpif.h'
  type(qtt_tensor) :: v1, v2, v3, dv
  type(qtt_matrix) :: A1, A2, I2, A3
  integer :: failed, i, j, k, N1, fd, d2
  integer, allocatable :: dd(:), nn(:)
  double precision, parameter :: pi = 4d0*atan(1d0)
  double precision :: alpha, dt, dx, x
  double precision, allocatable:: A(:), B(:), arr(:), ar2(:,:), b2(:,:), ar3(:,:,:)

     ! initialize system timer
     call system_timer_init()
     failed = 0

     print '("-- [14-1] testing:  v1 = qtt_tensor(d):")'
     call system_timer_start
     v1 = qtt_tensor(d)
     call system_timer_stop
     call v1% say()
     dt= system_dt
     write (*,'("[+][TEST14-1]["F7.2" s]: qtt tensor created")') dt

     print '(/,"-- [14-2] testing:  v1 = qtt_tensor(arr):")'
     if (d.lt.21) then
         N1 = 2**d - 2**(max(1,d-4))
     else
         print '("2^d too big, using N1 = 100000 instead")'
         N1 = 100000
     endif
     allocate(A(N1))
     dx= 1e-3
     do i=1,N1
        x = dble(i)*dx
        A(i)= dabs(dsin(x)*dexp(-1d1*x/(dx*N1)))
     enddo
     v1 = qtt_tensor(A, eps_=tol)
     call v1% say()

     B = v1% full()
     if (lsave) then
        open (fd, file='/tmp/AB.dat', status='replace')
        write (fd, '("# testing qtt% full ")')
        write (fd, '("# size(B) = ",I10)') size(B)
        write (fd, '("# 1:x 2:A 3:B")')
        do i=1,N1
           write (fd, '(10(ES14.7,1X))') i*dx, A(i), B(i)
        enddo
        do i=N1,size(B)
           write (fd, '(10(ES14.7,1X))') i*dx, 0d0, B(i)
        enddo
        close(fd)
        print '("file /tmp/AB.dat written")'
     endif

     print '(/,"-- [14-3] testing qtt assignment v2 = v1:")'
     call system_timer_start
     v2 = v1
     call system_timer_stop
     call v2% say()
     dt= system_dt
     write (*,'("[+][TEST14-3]["F7.2" s]: qtt assignment operator")') dt

     print '(/,"-- [14-4] testing:  v1 = qtt_ones(d):")'
     call system_timer_start
     v2 = qtt_ones(d)
     call system_timer_stop
     call v2% say()
     dt= system_dt
     write (*,'("[+][TEST14-4]["F7.2" s]: qtt_ones created")') dt

     print '(/,"-- [14-5] testing:  A1 = qtt_matrix(d) - all ranks 1:")'
     call system_timer_start
     A1 = qtt_matrix(d)
     call system_timer_stop
     call A1% say()
     dt= system_dt
     write (*,'("[+][TEST14-5]["F7.2" s]: qtt_matrix created")') dt

     print '(/,"-- [14-6] testing:  A1 = qtt_matrix(d, r):")'
     call system_timer_start
     A1 = qtt_matrix(d, r_=r)
     call system_timer_stop
     call A1% say()
     dt= system_dt
     write (*,'("[+][TEST14-6]["F7.2" s]: qtt_matrix created")') dt

     print '(/,"-- [14-6b] testing:  A1 = qtt_matrix(ar2(2^d,2^d), r):")'
     allocate(ar2(2**d,2**d))
     do i=1,2**d; do j=1,2**d
        ar2(i,j)= dsin(2*pi*(i-1)/dble(2**d))*dcos(2*pi*(j-1)/dble(2**d))
     enddo; enddo
     call system_timer_start
     A1 = qtt_matrix(ar2, eps_=tol)
     call system_timer_stop
     call A1% say()
     dt= system_dt
     write (*,'("[+][TEST14-6]["F7.2" s]: qtt_matrix created")') dt

     print '(/,"-- [14-6c] testing:  A1 = dtt_tensor(ar3(2^d,2^d,2^d), r):")'
     allocate(ar3(2**d,2**d,2**d), nn(d))
     nn = 8
     do i=1,2**d; do j=1,2**d; do k=1,2**d
        ar3(i,j,k)= dsin(2*pi*(i-1)/dble(2**d)) &
                   *dcos(2*pi*(j-1)/dble(2**d)) &
                   *dcos(2*pi*(k-1)/dble(2**d))
     enddo; enddo; enddo
     call system_timer_start
     A1 = dtt_tensor(interlace_3d(ar3), nn, eps_=tol)
     call system_timer_stop
     call A1% say()
     dt= system_dt
     write (*,'("[+][TEST14-6]["F7.2" s]: qtt_matrix created")') dt

     print '(/,"-- [14-7] testing assignment:  A2 = A1:")'
     call system_timer_start
     A1 = qtt_matrix(d, r_=r)
     call system_timer_stop
     call A1% say()
     dt= system_dt
     write (*,'("[+][TEST14-7]["F7.2" s]: qtt_matrix created")') dt

     print '(/,"-- [14-8] testing Kronecker product:  kron(v1, v2)")'
     call system_timer_start
     v3 = kron(v1, v2)
     call system_timer_stop
     call v3% say()
     dt= system_dt
     write (*,'("[+][TEST14-8]["F7.2" s]: kronecker product")') dt

     print '(/,"-- [14-9] testing 2x2 identity matrix")'
     call system_timer_start
     I2 = qtt_eye(1)
     call system_timer_stop
     dt= system_dt
     call I2% say()
     if (d.le.dmax4print) then
        ar2= I2% full_matrix()
        call pprint_matrix(ar2)
     endif
     write (*,'("[+][TEST14-9]["F7.2" s]: 2x2 identity matrix")') dt

     print '(/,"-- [14-10] testing construction of the d-dim I")'
     call system_timer_start
     A2 = qtt_eye(d)
     call system_timer_stop
     dt= system_dt
     call A2% say()
     if (d.le.dmax4print) then
        ar2= A2% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-10]["F7.2" s]: construction of d-dim I")') dt

     print '(/,"-- [14-11] testing multiplication by a unit matrix")'
     call system_timer_start
     call amen_mv(v2, A2, v1, 1d-8)
     call system_timer_stop
     dt= system_dt
     print '("Vector v2 after multiplying:")'
     call v2% say()
     dv = v1 - v2
     !print '("|v - I*v|/|v| = "ES14.7)', v3% normb()/v1% normb()
     write (*,'("[+][TEST14-11]["F7.2" s]: multiplication by a unit matrix")') dt

     print '(/,"-- [14-12] qtt shift operator")'
     call system_timer_start
     A1 = qtt_shift(d, -1)
     call system_timer_stop
     dt= system_dt
     call A1% say()
     if (d.le.dmax4print) then
        ar2= A1% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-12]["F7.2" s]: qtt shift operator")') dt

     print '(/,"-- [14-13] applying the qtt shift operator")'
     call system_timer_start
     call amen_mv(v2, A1, v1, 1d-8)
     call system_timer_stop
     dt= system_dt
     if (d.le.dmax4print.and.lsave) then
        A = v1% full()
        B = v2% full()
        open (fd, file='/tmp/AB_shift.dat', status='replace')
        write (fd, '("# testing qtt% full ")')
        write (fd, '("# size(B) = ",I10)') size(B)
        write (fd, '("# 1:x 2:A 3:B")')
        do i=1,size(B)
           write (fd, '(10(ES14.7,1X))') i*dx, A(i), B(i)
        enddo
        close(fd)
        print '("file /tmp/AB_shift.dat written")'
     endif
     write (*,'("[+][TEST14-13]["F7.2" s]: applying qtt shift operator")') dt

     print '(/,"-- [14-14] 1st-order qtt finite difference operator")'
     call system_timer_start
     A3 = qtt_fd1(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-14]["F7.2" s]: 1st-order FD operator")') dt

     print '(/,"-- [14-15] applying 1st-order fd operator")'
     N1= 2**d
     deallocate(A); allocate(A(N1))
     deallocate(B); allocate(B(N1))
     dx= 1d-3
     do i=1,N1
        x = (dble(i) - 1d0)*dx
        A(i)= dsin(x)*dexp(-x/(dx*N1))
        x = (dble(i) - .5d0)*dx
        B(i)= dexp(-x/(dx*N1))*(dcos(x) - dsin(x)/(dx*N1))
     enddo
     v1= qtt_tensor(A, eps_=tol)
     v2= qtt_tensor(B, eps_=tol)
     print '("v1 = ")'; call v1% say()
     print '("v2 = ")'; call v2% say()
     call system_timer_start
     A3 = qtt_fd1(d, dx)
     call amen_mv(v3, A3, v1, 1d-8)
     if (d.le.dmax4print.and.lsave) then
        A = v3% full()
        B = v2% full()
        open (fd, file='/tmp/dvdx_fd1.dat', status='replace')
        write (fd, '("# testing dv/dx - D*v")')
        write (fd, '("# 1:x 2:dv/dx 3:D*v")')
        do i=1,size(B)
           write (fd, '(10(ES14.7,1X))') i*dx, B(i), A(i)
        enddo
        close(fd)
        print '("file /tmp/dvdx_fd1.dat written")'
     endif
     call system_timer_stop
     dt= system_dt
     dv = v3 - v2
     print '("|dv/dx - D*v|/|v| = "ES14.7)', dv% normb()/v2% normb()
     write (*,'("[+][TEST14-15]["F7.2" s]: applying 1st-order FD operator")') dt

     print '(/,"-- [14-16] 2nd-order qtt finite difference operator")'
     call system_timer_start
     A3 = qtt_fd21(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-16]["F7.2" s]: 2nd-order FD operator")') dt

     print '(/,"-- [14-17] applying FD21 operator")'
     N1= 2**d
     deallocate(A); allocate(A(N1))
     deallocate(B); allocate(B(N1))
     dx= 1d-3
     do i=1,N1
        x = (dble(i) - 1d0)*dx
        A(i)= dsin(x)*dexp(-x/(dx*N1))
        !if (i.eq.1)  x = .5d0*dx
        !if (i.eq.N1) x = x - .5d0*dx
        B(i)= dexp(-x/(dx*N1))*(dcos(x) - dsin(x)/(dx*N1))
     enddo
     v1= qtt_tensor(A, eps_=tol)
     v2= qtt_tensor(B, eps_=tol)
     print '("v1 = ")'; call v1% say()
     print '("v2 = ")'; call v2% say()
     call system_timer_start
     A3 = qtt_fd21(d, dx)
     call amen_mv(v3, A3, v1, 1d-8)
     if (d.le.dmax4print.and.lsave) then
        A = v3% full()
        B = v2% full()
        open (fd, file='/tmp/dvdx_fd21.dat', status='replace')
        write (fd, '("# testing dv/dx - D*v")')
        write (fd, '("# 1:x 2:dv/dx 3:D*v")')
        do i=1,size(B)
           write (fd, '(10(ES14.7,1X))') i*dx, B(i), A(i)
        enddo
        close(fd)
        print '("file /tmp/dvdx_fd21.dat written")'
     endif
     call system_timer_stop
     dt= system_dt
     dv = v3 - v2
     print '("|dv/dx - D*v|/|v| = "ES14.7)', dv% normb()/v2% normb()
     write (*,'("[+][TEST14-17]["F7.2" s]: applying 1st-order FD operator")') dt

     print '(/,"-- [14-18] 2nd-order qtt Laplace operator")'
     call system_timer_start
     A3 = qtt_laplace_dd(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        b2= ar2
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-18]["F7.2" s]: Laplace FD operator")') dt

     print '(/,"-- [14-19] 2nd-order qtt Laplace operator, DN")'
     call system_timer_start
     A3 = qtt_laplace_dn(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-19]["F7.2" s]: Laplace FD operator, DN")') dt

     print '(/,"-- [14-20] 2nd-order qtt Laplace operator, ND")'
     call system_timer_start
     A3 = qtt_laplace_nd(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-20]["F7.2" s]: Laplace operator, ND")') dt

     print '(/,"-- [14-21] 2nd-order qtt Laplace operator, NN")'
     call system_timer_start
     A3 = qtt_laplace_nn(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-21]["F7.2" s]: Laplace operator, NN")') dt

     print '(/,"-- [14-22] 2nd-order qtt Laplace operator, periodic")'
     call system_timer_start
     A3 = qtt_laplace_p(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        call pprint_matrix(ar2, frmt_='(I3)')
     endif
     write (*,'("[+][TEST14-22]["F7.2" s]: Laplace operator, periodic")') dt

     print '(/,"-- [14-23] 2nd-order 2D QTT Laplace operator")'
     dd= [d, d-1]
     call system_timer_start
     A3 = qtt_laplace_dd(dd)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (dd(1) + dd(2).le.dmax4print) then
        ar2= A3% full_matrix()
        open (fd, file='/tmp/A2d.dat', status='replace')
        call pprint_matrix(ar2, frmt_='(I3)', out_unit_=fd, &
                           line_len_=10000, &
                           max_rows_=10000)
        close(fd)
        print '("file /tmp/A2d.dat written")'

     endif
     write (*,'("[+][TEST14-23]["F7.2" s]: 2nd-order 2D QTT Laplace operator")') dt

     print '(/,"-- [14-24] 3D QTT Laplace operator")'
     dd= [d-1, d, d-1]
     call system_timer_start
     A3 = qtt_laplace_dd(dd)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (sum(dd).le.dmax4print.and.lsave) then
        ar2= A3% full_matrix()
        open (fd, file='/tmp/A3d.dat', status='replace')
        call pprint_matrix(ar2, frmt_='(I3)', out_unit_=fd, &
                           line_len_=10000, &
                           max_rows_=10000)
        close(fd)
        print '("file /tmp/A3d.dat written")'

     endif
     write (*,'("[+][TEST14-24]["F7.2" s]: 3D QTT Laplace operator")') dt

     print '(/,"-- [14-25] 1D QTT inverse Laplacian")'
     call system_timer_start
     A3 = qtt_laplace_dd_inv(d)
     call system_timer_stop
     dt= system_dt
     call A3% say()
     if (d.le.dmax4print) then
        ar2= A3% full_matrix()
        call pprint_matrix(ar2*(2**d+1), frmt_='(ES11.4)')
        !call pprint_matrix(matmul(ar2,b2), frmt_='(ES11.4)')
        !print '("---")'
        I2 = qtt_eye(d)
        print '("|Laplace @ InvLaplace - I|_2 = ",ES11.4)', &
               normfro(matmul(ar2,b2) - I2% full_matrix())
     endif
     write (*,'(/,"[+][TEST14-25]["F7.2" s]: 1D QTT inverse Laplacian")') dt

     if(failed.eq.0) then
        print '("[+][TEST14]: PASSED")'
     else
        print '("[-][TEST14]: FAILED with ",I2," error(s)")', failed
        error stop failed
     end if

  end subroutine test_qtt_ops

end program
