program main
    use iso_fortran_env, dp => real64
    
    implicit none
    
    character(len=20) :: func
    integer :: ncell,d,norder, tt_size, max_rank
    real(dp) :: err, total_time, time_cross, avg_rank

    character(len=100) :: filename, filename_time, filename_size, filename_time_cross, filename_avg_rank, filename_max_rank, file_id, arg
    real(dp), allocatable :: xgrid(:)
    
    func = 'f8'
    d = 100
    ncell = 2
    norder = 10
    
    !call get_command_argument(1, arg)
    !read(arg, *) norder

    write(file_id, '(i0)') norder

    filename            = 'fortran_plot_data/output_' // trim(adjustl(file_id)) // '.dat'
    filename_time       = 'fortran_plot_data/output_' // trim(adjustl(file_id)) // '_time.dat'
    filename_size       = 'fortran_plot_data/output_' // trim(adjustl(file_id)) // '_size.dat'
    filename_time_cross = 'fortran_plot_data/output_' // trim(adjustl(file_id)) // '_time_cross.dat'
    filename_avg_rank   = 'fortran_plot_data/output_' // trim(adjustl(file_id)) // '_avg_rank.dat'
    filename_max_rank   = 'fortran_plot_data/output_' // trim(adjustl(file_id)) // '_max_rank.dat'

    
    call tt_GL_fn(func,d,ncell,norder,err,total_time,tt_size,time_cross,avg_rank,max_rank)
    
    !open(unit=10, file=filename, status='replace', action='write')
    !write(10,*) err   
    !close(10)
    !open(unit=20, file=filename_time, status='replace', action='write')
    !write(20,*) total_time
    !close(20)
    !open(unit=30, file=filename_size, status='replace', action='write')
    !write(30,*) tt_size
    !close(30)
    !open(unit=40, file=filename_time_cross, status='replace', action='write')
    !write(40,*) time_cross
    !close(40)
    !open(unit=50, file=filename_avg_rank, status='replace', action='write')
    !write(50,*) avg_rank
    !close(50)
    !open(unit=60, file=filename_max_rank, status='replace', action='write')
    !write(60,*) max_rank
    !close(60)
    
    print *,'Error:',err
    print *,'Total time:',total_time
    print *,'Size of TT:',tt_size
    print *,'Time for DMRG cross:',time_cross
    print *,'Average rank:',avg_rank
    print *,'Maximum rank:',max_rank


contains

subroutine tt_GL_fn(func,d,ncell,norder,err,total_time,tt_size,time_cross,avg_rank,max_rank)
    use stdlib_quadrature, only: gauss_legendre
    use thor_lib, only: dtt_tensor, core2cell_dtt, dtt_tensor_ones
    use dmrgg_lib, only: dtt_dmrgg
    use matlab_struct_module, only: array3d
    
    implicit none

    include 'mpif.h'

    ! Variable declarations
    integer, intent(in) :: ncell,d,norder
    integer :: i,j,k,info, tt_size, max_rank
    integer, allocatable :: modes(:)
    real(dp) :: a, b, dx, I_exact, err, start_time, finish_time, total_time, tc1, tc2, time_cross, avg_rank
    real(dp), dimension(2) :: interval
    real(dp), allocatable :: w(:), tempx(:), tempw(:)
    real(dp), allocatable :: S1(:,:),S2(:,:),S3(:,:)
    type(dtt_tensor) :: ftt
    type(array3d), allocatable :: cc(:) 
    character(len=20),intent(in) :: func

    call mpi_init(info)

    call cpu_time(start_time)

    a = 0.0d0
    b = 1.0d0
    
    dx = (b - a) / real(ncell, dp)

    allocate(modes(d))
    modes(:) = ncell*norder
    
    allocate(xgrid(ncell * norder))
    allocate(w(ncell * norder))

    do i = 1, ncell
        interval(1) = a + real(i - 1, dp) * dx
        interval(2) = interval(1) + dx
        
        allocate(tempx(norder), tempw(norder))
        call gauss_legendre( tempx, tempw, interval )

        xgrid((i-1)*norder+1:i*norder) = tempx
        w((i-1)*norder+1:i*norder) = tempw

        deallocate(tempx, tempw)
    end do

    ftt = dtt_tensor_ones(modes)

    select case (trim(func))
        case ('f1')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f1, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = 0.5d0
        case ('f2')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f2, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = 1.0d0
        case ('f3')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f3, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = 1.0d0/gamma(real(d,dp)+2.0d0)
        case ('f4')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f4, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = ( sqrt(3.141592653589793d0)/2.0d0 * erf(1.0d0) )**real(d,dp)
        case ('f5')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f5, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = (1.0d0-1.0d0/exp(1.0d0))**real(d,dp)
        case ('f6')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f6, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = 1.0d0
        case ('f7')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f7, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = 1.0d0
        case ('f8')
            call cpu_time(tc1)
            call dtt_dmrgg(ftt, f8, accuracy=1.0d-20)
            call cpu_time(tc2)
            I_exact = -1.0d0/99.0d0
    end select
    
    cc = core2cell_dtt(ftt)
    
    allocate(S1(1,1))
    S1(1,1) = 1.0d0

    do k = 1,d
    
        allocate(S2( size(cc(k)%arr(:,1,1)),size(cc(k)%arr(1,1,:))  ))
        
        do i = 1,size(S2(:,1))
            do j = 1,size(S2(1,:))
                S2(i,j) = dot_product(  cc(k)%arr(i,:,j),w(:)  ) 
            end do
        end do    
    
        S3 = matmul(S1,S2)
    
        deallocate(S1)
        deallocate(S2)
    
        allocate(S1( size(S3(:,1)),size(S3(1,:))  ))
        S1 = S3
        
    end do

    call cpu_time(finish_time)

    total_time = finish_time - start_time

    time_cross = tc2 - tc1
    
    err = abs(I_exact - S3(1,1))

    tt_size = ftt%size()
    
    avg_rank = sum(ftt%r(:ftt%m))/real(size(ftt%r(:ftt%m)),dp)

    max_rank = maxval(ftt%r(:ftt%m))

    call mpi_finalize(info)

end subroutine tt_GL_fn





pure function f1(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)
    real(dp) :: c, w, sum_val

    c = 0.5d0
    w = 0.0d0

    sum_val = 0.0d0
    do i = 1,d
        x(i) = xgrid(ind(i))
        sum_val = sum_val + c*x(i)
    end do
    
    y = cos( 2.0d0 * 3.141592653589793d0 * w + 2.0d0 * 3.141592653589793d0 * sum_val )
    y = y * y

end function f1


pure function f2(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)
    real(dp) :: fx,wl, cl, s
    
    wl = 0.0d0
    cl = 1.0d0
    s  = 4.0d0 / 3.141592653589793d0
    
    do i = 1,d
        x(i) = xgrid(ind(i))
    end do
    
    y = 1.0d0
    do i = 1, d
        y = y * ( s / ( cl**2 + (x(i) - wl)**2 ) )   
    end do
    
end function f2



pure function f3(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)
    real(dp) :: c, sum_val

    c = 1.0d0

    sum_val = 0.0d0
    do i = 1,d
        x(i) = xgrid(ind(i))
        sum_val = sum_val + c*x(i)
    end do
    
    y = (1+ sum_val)**real((-d-1),dp)

end function f3


pure function f4(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)
    real(dp) :: c, sum_val

    c = 1.0d0

    sum_val = 0.0d0
    do i = 1,d
        x(i) = xgrid(ind(i))
        sum_val = sum_val + c*x(i)**2
    end do
    
    y = exp(-sum_val)

end function f4


pure function f5(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)
    real(dp) :: c, sum_val

    c = 1.0d0

    sum_val = 0.0d0
    do i = 1,d
        x(i) = xgrid(ind(i))
        sum_val = sum_val + c*x(i)
    end do
    
    y = exp(-sum_val)

end function f5



pure function f6(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)
    real(dp) :: ak
    
    ak = 1.0d0

    y = 1.0d0
    do i = 1,d
        x(i) = xgrid(ind(i))
        y = y * abs(4.0d0 * x(i) - 3.0d0 + ak)/(1.0d0+ak)
    end do
    
end function f6


pure function f7(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)

    y = 1.0d0
    do i = 1,d
        x(i) = xgrid(ind(i))
        y = y * abs(4.0d0 * x(i) - 3.0d0) - 0.25d0
    end do
    
end function f7


pure function f8(d,ind,nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    real(dp) :: sgn
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    double precision :: x(tt_size)
    real(dp) :: poly_sum
  
    poly_sum = 0.0d0
    
    do i = 1, d
     x(i) = xgrid(ind(i))
     poly_sum = poly_sum + (512.0d0 * x(i)**10 - 1280.0d0 * x(i)**8 + 1120.0d0 * x(i)**6 - &
                             400.0d0 * x(i)**4 + 50.0d0 * x(i)**2 - 1.0d0)
    end do
    
    if ((x(1) - 0.5d0) > 0.0d0) then
      sgn = 1.0d0
    else if ((x(1) - 0.5d0) < 0.0d0) then
      sgn = -1.0d0
    else
      sgn = 0.0d0
    end if
    
    y = (x(1) - 0.5d0)**10 * sgn + (1.0d0 / real(d,dp)) * poly_sum
    
end function f8
    
end program 
    
