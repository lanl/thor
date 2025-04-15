module ttstsc_data
use, intrinsic :: iso_fortran_env, only: dp => real64
real(dp), allocatable ::XPt(:),Tvec(:)
end module ttstsc_data

program main
    use ttstsc_data
    use thor_lib, only: dttm_from_cell_array,dtt_norm,dtt_matrix,dttm_add,dtt_tensor_ones,dtt_tensor,core2cell_dtt,dtt_from_cell_array,dttm_sub, &
    dttm_from_dtt_tensor,dtt_round,dtt_sub,core2cell_dttm,dtt_size,dtt_from_dttm
    use matlab_struct_module, only: cell4d_array,array3d
    use dmrgg_lib, only: dtt_dmrgg
    use ttamen_lib, only: amen_mm,amen_mv,dtt_amen_solve


    implicit none

    include 'mpif.h'
    
    integer :: Ns,i,TN,N1,info,modes(4),j,k,l
    integer, allocatable :: It(:),Ix(:)
    real(dp) :: a,b,tol,Err,start_time,finish_time
    real(dp), allocatable :: ET(:),DMA(:,:),Lapmat(:,:),At(:,:),Iden(:,:),tempten(:,:,:,:)
    type(cell4d_array) :: cell_array
    type(dtt_matrix) :: Btt,Bmaptt,Att,Amaptt,ttm1,ttm2,ttm3,Laptt,Lapmaptt,amen_mm_out,Koptt,boptt,Convtt1,Convtt2,Convtt3,Convtt,Convmaptt
    type(dtt_matrix) :: coptt,Reactmaptt,Reacttt
    type(array3d), allocatable :: temp_cell1(:),temp_cell2(:),G(:)
    type(dtt_tensor) :: kfuntt,btemptt,cfuntt,ftt,Gbc,Gint,Gtt,amen_mv_out,RHStt,ttsol,exacttt
    character(len=25) :: arg
    character(len=100) :: filename_err, filename_time, file_id


    call get_command_argument(1, arg)
    read(arg, *) Ns
    
    call mpi_init(info)

    start_time = MPI_Wtime()

    tol = 1.0d-12
    
    TN = Ns
    N1 = Ns

    a = -1.0d0
    b = 1.0d0 

    ET =  ZECHGL(Ns-1)
    Tvec = (ET+1)*0.5

    DMA = dmchgl(Ns-1,ET)
    At = 2.0d0*DMA

    XPt = (b-a)/2.0d0 * ET + (a+b)/2.0d0
    DMA = (2.0d0/(b-a))*DMA

    Lapmat = matmul(DMA,DMA)

    modes(1)  = size(Tvec)
    modes(2:) = size(Xpt)

    allocate(Iden(Ns,Ns))

    Iden = 0.0d0

    do i = 1,Ns
        Iden(i,i) = 1.0d0
    end do

    allocate(It(1:TN-1))
    It = [(i, i=2,TN)]

    allocate(Ix(1:N1-2))
    Ix = [(i, i=2,N1-1)]


    allocate(cell_array%cells(4))
    
    cell_array%cells(1)%arr = reshape((At(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1]) 
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1]) 
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1]) 
    Btt = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((At(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    Bmaptt = dttm_from_cell_array(cell_array)

    Att = Btt
    Amaptt = Bmaptt


    !!!!!!!!!!!!!!!!!! Laplace operator !!!!!!!!!!!!!!!!!!

    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Lapmat(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    ttm1 = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Lapmat(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    ttm2 = dttm_from_cell_array(cell_array)


    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Lapmat(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    ttm3 = dttm_from_cell_array(cell_array)

    Laptt = dttm_add(ttm1,dttm_add(ttm2,ttm3))

    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((Lapmat(Ix,:)),[1,size(Ix),TN,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    ttm1 = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(3)%arr = reshape((Lapmat(Ix,:)),[1,size(Ix),TN,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    ttm2 = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(4)%arr = reshape((Lapmat(Ix,:)),[1,size(Ix),TN,1])  
    ttm3 = dttm_from_cell_array(cell_array)


    Lapmaptt = dttm_add(ttm1,dttm_add(ttm2,ttm3))

    kfuntt = dtt_tensor_ones(modes)
    call dtt_dmrgg(kfuntt, kfun_dmrg, accuracy=tol)
    
    temp_cell1 = core2cell_dtt(kfuntt)
    allocate(temp_cell2(4))
    temp_cell2(1)%arr = temp_cell1(1)%arr(:,2:TN,:)
    do i = 2,4
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:N1-1,:)
    end do
    kfuntt = dtt_from_cell_array(temp_cell2)
    deallocate(temp_cell2)

   
    G = core2cell_dtt(kfuntt)
    do i=1,4
    allocate(  tempten(kfuntt%r(i-1),size(G(i)%arr(j,:,k)),size(G(i)%arr(j,:,k)),kfuntt%r(i))  )
    tempten = 0.0d0
      do j = 1,kfuntt%r(i-1)
        do k = 1,kfuntt%r(i)
            do l = 1,size(G(i)%arr(j,:,k))
            tempten(j,l,l,k) = G(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    Koptt = dttm_from_cell_array(cell_array)

    call amen_mm(amen_mm_out, Koptt, Laptt, tol)
    Laptt = amen_mm_out

    call amen_mm(amen_mm_out, Koptt, Lapmaptt, tol)
    Lapmaptt = amen_mm_out


    Att = dttm_sub(Att,Laptt)

    Amaptt = dttm_sub(Amaptt,Lapmaptt)


    !!!!!!!!!!! Convection operator  !!!!!!!!!!!


    btemptt = dtt_tensor_ones(modes)
    call dtt_dmrgg(btemptt, bfun_dmrg, accuracy=tol)
    
    temp_cell1 = core2cell_dtt(btemptt)
    allocate(temp_cell2(4))
    temp_cell2(1)%arr = temp_cell1(1)%arr(:,2:N1,:)
    do i = 2,4
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:N1-1,:)
    end do
    btemptt = dtt_from_cell_array(temp_cell2)
    deallocate(temp_cell2)


    G = core2cell_dtt(btemptt)
    do i=1,4
    allocate(  tempten(btemptt%r(i-1),size(G(i)%arr(j,:,k)),size(G(i)%arr(j,:,k)),btemptt%r(i))  )
    tempten = 0.0d0
      do j = 1,btemptt%r(i-1)
        do k = 1,btemptt%r(i)
            do l = 1,size(G(i)%arr(j,:,k))
            tempten(j,l,l,k) = G(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    boptt = dttm_from_cell_array(cell_array)


    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((DMA(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    Convtt1 = dttm_from_cell_array(cell_array)

    call amen_mm(amen_mm_out, boptt, Convtt1, tol)
    Convtt1 = amen_mm_out


    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((DMA(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    Convtt2 = dttm_from_cell_array(cell_array)

    call amen_mm(amen_mm_out, boptt, Convtt2, tol)
    Convtt2 = amen_mm_out


    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((DMA(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    Convtt3 = dttm_from_cell_array(cell_array)

    call amen_mm(amen_mm_out, boptt, Convtt3, tol)
    Convtt3 = amen_mm_out

    Convtt= dttm_add(Convtt1,dttm_add(Convtt2,Convtt3))



    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((DMA(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    Convtt1 = dttm_from_cell_array(cell_array)


    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(3)%arr = reshape((DMA(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    Convtt2 = dttm_from_cell_array(cell_array)
    

    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(4)%arr = reshape((DMA(Ix,:)),[1,size(Ix),N1,1])  
    Convtt3 = dttm_from_cell_array(cell_array)


    Convmaptt= dttm_add(Convtt1,dttm_add(Convtt2,Convtt3))

    Att = dttm_add(Att,Convtt)
    Amaptt = dttm_add(Amaptt,Convmaptt)

    !!!!!!!!!!! Reaction Term  !!!!!!!!!!!

    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    Reacttt = dttm_from_cell_array(cell_array)


    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),TN,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),N1,1])  
    Reactmaptt = dttm_from_cell_array(cell_array)



    cfuntt = dtt_tensor_ones(modes)
    call dtt_dmrgg(cfuntt, cfun_dmrg, accuracy=tol)


    temp_cell1 = core2cell_dtt(cfuntt)
    allocate(temp_cell2(4))
    temp_cell2(1)%arr = temp_cell1(1)%arr(:,2:TN,:)
    do i = 2,4
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:N1-1,:)
    end do
    cfuntt = dtt_from_cell_array(temp_cell2)
    deallocate(temp_cell2)


    G = core2cell_dtt(cfuntt)
    do i=1,4
    allocate(  tempten(cfuntt%r(i-1),size(G(i)%arr(j,:,k)),size(G(i)%arr(j,:,k)),cfuntt%r(i))  )
    tempten = 0.0d0
      do j = 1,cfuntt%r(i-1)
        do k = 1,cfuntt%r(i)
            do l = 1,size(G(i)%arr(j,:,k))
            tempten(j,l,l,k) = G(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    coptt = dttm_from_cell_array(cell_array)


    call amen_mm(amen_mm_out, coptt, Reacttt, tol)
    Reacttt = amen_mm_out

    call amen_mm(amen_mm_out, coptt, Reactmaptt, tol)
    Reactmaptt = amen_mm_out

    Att = dttm_add(Att,Reacttt)
    Amaptt = dttm_add(Amaptt,Reactmaptt)

    Att = dttm_from_dtt_tensor(dtt_round(Att,tol_=tol),Att%q(1:4),Att%s(1:4))
    Amaptt = dttm_from_dtt_tensor(dtt_round(Amaptt,tol_=tol),Amaptt%q(1:4),Amaptt%s(1:4))

    !!!!!!!!!!!!!!!!!!

    ftt = dtt_tensor_ones(modes)
    call dtt_dmrgg(ftt, rhsfn_dmrg, accuracy=tol)


    temp_cell1 = core2cell_dtt(ftt)
    allocate(temp_cell2(4))
    temp_cell2(1)%arr = temp_cell1(1)%arr(:,2:TN,:)
    do i = 2,4
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:N1-1,:)
    end do
    ftt = dtt_from_cell_array(temp_cell2)
    deallocate(temp_cell2)


    Gtt = dtt_tensor_ones(modes)
    call dtt_dmrgg(Gtt, gfn_dmrg, accuracy=tol)


    G = core2cell_dtt(Gtt)
    G(1)%arr(:,1,:) = 0.0d0
    do i = 2,4
        G(i)%arr(:,1,:) = 0.0d0
        G(i)%arr(:,N1,:) = 0.0d0
    end do
    Gint = dtt_from_cell_array(G)


    Gbc = dtt_sub(Gtt, Gint)

    call amen_mv(amen_mv_out, Amaptt, Gbc, tol=tol)

    RHStt = dtt_sub(ftt,amen_mv_out)
    
    ttsol = dtt_amen_solve(Att, RHStt, tol=tol)

    finish_time = MPI_Wtime()

    exacttt = dtt_tensor_ones(modes)
    call dtt_dmrgg(exacttt, exactfn_dmrg, accuracy=tol)
    

    temp_cell1 = core2cell_dtt(exacttt)
    allocate(temp_cell2(4))
    temp_cell2(1)%arr = temp_cell1(1)%arr(:,2:TN,:)
    do i = 2,4
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:N1-1,:)
    end do
    exacttt = dtt_from_cell_array(temp_cell2)
    deallocate(temp_cell2)

    Err = dtt_norm(dtt_round(dtt_sub(ttsol,exacttt)))/dtt_norm(dtt_round(exacttt))
    

    print *,'########'
    print *,'Error=',Err
    print *,'Compression of Att=',real(dtt_size(dtt_from_dttm(Att)),dp) / product( real( Att%q(:4) * Att%s(:4),dp ))
    print *,'Time=',finish_time-start_time
    print *,'########'


    write(file_id, '(i0)') Ns
    filename_err   = '../plot_data/Fortran/SC_err_' // trim(adjustl(file_id)) // '.dat'
    filename_time  = '../plot_data/Fortran/SC_time_' // trim(adjustl(file_id)) // '.dat'

    open(unit=10, file=filename_err, status='replace', action='write')
    write(10,*) Err   
    close(10)
    open(unit=20, file=filename_time, status='replace', action='write')
    write(20,*) finish_time-start_time
    close(20)



    call mpi_finalize(info)



contains


function ZECHGL(N) result(ET)
    implicit none
    integer, intent(in) :: N
    real(dp), allocatable :: ET(:)
    real(dp) :: PI, C, ETX
    integer :: N2, I

    ! Handle N == 0
    if (N == 0) then
        allocate(ET(1))
        ET(1) = 0.0_8
        return
    end if

    ! Allocate ET(0:N) since MATLAB indices start at 1 but Fortran at 1 too here
    allocate(ET(0:N))

    ! Initialize ET(0) to -1
    ET(0) = -1.0_8

    ! Handle N == 1
    if (N == 1) then
        return
    end if

    ! Initialize ET(N) to 1
    ET(N) = 1.0_8

    ! Handle N == 2
    if (N == 2) then
        return
    end if

    ! Calculate N2
    N2 = (N - 1) / 2

    ! Initialize ET(N2) to 0
    ET(N2) = 0.0_8

    ! Define PI
    PI = 3.14159265358979323846_8

    ! Calculate constant C
    C = PI / real(N, kind=8)

    ! Loop to calculate ET values
    do I = 1, N2
        ETX = cos(C * real(I, kind=8))
        ET(I)     = -ETX
        ET(N - I) = ETX
    end do

end function ZECHGL





function dmchgl(N, ET) result(DMA)
    implicit none
    integer, intent(in) :: N
    real(8), intent(in) :: ET(:)
    real(8) :: DMA(N+1,N+1)
    
    integer :: I, J, DN, N2, SN
    real(8) :: CN, EJ, EI
    integer :: SGN, SGNI, SGNJ

    ! Initialize DMA with zeros
    DMA = 0.0d0

    if (N == 1) then
        DMA(1,1) = 0.0d0
        return
    end if

    DN = N
    CN = (2.0d0 * DN * DN + 1.0d0) / 6.0d0

    if (mod(N, 2) == 0) then
        N2 = N / 2
        SN = 1 + 4 * N2 - 2 * N
    else
        N2 = (N - 1) / 2
        SN = 1 + 4 * N2 - 2 * N
    end if

    DMA(1,1)       = -CN
    DMA(N+1,N+1)   =  CN
    DMA(1,N+1)     = -0.5d0 * SN
    DMA(N+1,1)     =  0.5d0 * SN

    if (N == 2) return

    SGN = -1
    do J = 2, N
        EJ = ET(J)
        DMA(1,J)     = (-2.0d0 * SGN) / (1.0d0 + EJ)
        DMA(N+1,J)   = ( 2.0d0 * SGN * SN) / (1.0d0 - EJ)
        DMA(J,1)     = ( 0.5d0 * SGN) / (1.0d0 + EJ)
        DMA(J,N+1)   = -(0.5d0 * SGN * SN) / (1.0d0 - EJ)
        SGN = -SGN
    end do

    SGNI = -1
    do I = 2, N
        EI = ET(I)
        SGNJ = -1
        do J = 2, N
            if (I /= J) then
                EJ = ET(J)
                DMA(I,J) = (SGNI * SGNJ) / (EI - EJ)
            else
                DMA(I,J) = -0.5d0 * EI / (1.0d0 - EI * EI)
            end if
            SGNJ = -SGNJ
        end do
        SGNI = -SGNI
    end do

end function dmchgl



function kfun(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans

  ans = 1.0d0
end function kfun


function kfun_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    use ttstfd_data
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = Tvec(ind(1))

    do i = 2,dpt
        x_temp(i) = XPt(ind(i))
    end do

    y = kfun(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function kfun_dmrg


function bfun(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans

  ans = 1.0d0
end function bfun


function bfun_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    use ttstsc_data
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = Tvec(ind(1))

    do i = 2,dpt
        x_temp(i) = XPt(ind(i))
    end do

    y = bfun(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function bfun_dmrg


function cfun(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans

  ans = 1.0d0
end function cfun



function cfun_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    use ttstsc_data
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = Tvec(ind(1))

    do i = 2,dpt
        x_temp(i) = XPt(ind(i))
    end do

    y = cfun(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function cfun_dmrg





function rhsfn(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans
  real(dp) :: pi, sum1, theta

  pi = 4 * atan (1.0d0)
  sum1 = t + x + y + z
  theta = pi * sum1 * (2.0d0 / 11.0d0)

  ans = t / 2.0d0 + x / 2.0d0 + y / 2.0d0 + z / 2.0d0 + &
      sin(theta) + &
      pi * cos(theta) * (8.0d0 / 11.0d0) + &
      pi**2 * sin(theta) * (12.0d0 / 121.0d0) + &
      sum1**2 / 16.0d0 - 3.0d0 / 8.0d0
      
end function rhsfn



function rhsfn_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    use ttstsc_data
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = Tvec(ind(1))

    do i = 2,dpt
        x_temp(i) = XPt(ind(i))
    end do

    y = rhsfn(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function rhsfn_dmrg




function gfn(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans
  real(dp) :: pi, sum1, theta

  pi = 4 * atan (1.0d0)
  sum1 = t + x + y + z
  theta = pi * sum1 * (2.0d0 / 11.0d0)

  ans = sin(theta) + sum1**2 / 16.0d0
end function gfn



function gfn_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    use ttstsc_data
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = Tvec(ind(1))

    do i = 2,dpt
        x_temp(i) = XPt(ind(i))
    end do

    y = gfn(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function gfn_dmrg





function exactfn(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans
  real(dp) :: pi, sum1
  
  pi = 4 * atan (1.0d0)
  sum1 = t + x + y + z

  ans = sin(pi * sum1 * (2.0d0 / 1.1d1)) + (sum1**2) / 1.6d1
end function exactfn



function exactfn_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    use ttstsc_data
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = Tvec(ind(1))

    do i = 2,dpt
        x_temp(i) = XPt(ind(i))
    end do

    y = exactfn(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function exactfn_dmrg





end program 
    
