module data
use, intrinsic :: iso_fortran_env, only: dp => real64
use thor_lib, only: dtt_matrix,dtt_tensor
use matlab_struct_module, only: cell4d_array

real(dp), allocatable ::XPt(:),Tvec(:)
real(dp) :: tol
type(dtt_matrix) :: Bmaptt,Lapmaptt,Convttmap1,Convttmap2,Convttmap3,Btt,Convtt1,Convtt2,Convtt3
type(cell4d_array) :: cell_array
type(dtt_tensor) :: Gbc

end module data

program main
    use data
    use matlab_struct_module, only: cell4d_array,array3d
    use thor_lib, only: dttm_from_cell_array,dtt_matrix,dtt_norm,dttm_full,dttm_add,dtt_tensor_ones,dtt_tensor,dtt_sub,core2cell_dtt,dtt_from_cell_array,dtt_tensor_zeros,dtt_round,dtt_mul_d,dtt_add,dtt_full
    use dmrgg_lib, only: dtt_dmrgg
    use ttamen_lib, only: dtt_amen_solve
    
    implicit none

    include 'mpif.h'

    integer :: Ns,SN,N1,i,modes(4),info,maxiter,k,t
    real(dp) :: a,b,t0,t1,nrmFU0,epsk,alphas(5),local_err,res_norm,Newtoneps,Err,start_time,finish_time
    real(dp), allocatable :: temppts(:),tempDMA(:,:),DMA(:,:),Lapmat(:,:),At(:,:),Iden(:,:),It(:),Ix(:)
    type(dtt_matrix) :: ttm1,ttm2,ttm3,Laptt,opt1,opt2,Jmatk
    type(dtt_tensor) :: Gtt,Gint,U0,U,U1,du,b_fortran,U2,exacttt
    type(array3d), allocatable :: G(:),temp_cell1(:),temp_cell2(:)
    double precision,allocatable :: full_check(:,:),full_check_vec(:)
    character(len=25) :: arg
    character(len=100) :: filename_err, filename_time, file_id


    call get_command_argument(1, arg)
    read(arg, *) Ns
    
    call mpi_init(info)

    start_time = MPI_Wtime()

    Newtoneps = 1.0d-5
    
    SN = Ns-1
    N1 = SN+1

    a = 0.0d0
    b = 6.0d0
    tol = 1.0d-5
    

    temppts =  ZECHGL(SN)
    tempDMA = dmchgl(SN,temppts)

    XPt = ((b - a) / 2.0d0) * temppts + (a + b) / 2.0d0
    DMA = (2.0d0/(b-a))*tempDMA;

    Lapmat = matmul(DMA,DMA)

    t0 = 0.0d0;
    t1 = 1.0d0;

    Tvec = ((t1-t0)/2.0d0) * temppts + (t1+t0)/2.0d0

    modes(1)  = size(Tvec)
    modes(2:) = size(Xpt)

    At = (2.0d0/(t1-t0))*tempDMA

    allocate(Iden(Ns,Ns))
    Iden = 0.0d0
    do i = 1,Ns
        Iden(i,i) = 1.0d0
    end do

    allocate(It(1:Ns-1))
    It = [(i, i=2,Ns)]

    allocate(Ix(1:Ns-2))
    Ix = [(i, i=2,Ns-1)]


    allocate(cell_array%cells(4))
    
    cell_array%cells(1)%arr = reshape((At(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1]) 
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1]) 
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1]) 
    Btt = dttm_from_cell_array(cell_array)   

    cell_array%cells(1)%arr = reshape((At(It,:)),[1,size(It),Ns,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    Bmaptt = dttm_from_cell_array(cell_array)




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


    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),Ns,1]) 
    cell_array%cells(2)%arr = reshape((Lapmat(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    ttm1 = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),Ns,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(3)%arr = reshape((Lapmat(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    ttm2 = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),Ns,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(4)%arr = reshape((Lapmat(Ix,:)),[1,size(Ix),Ns,1])  
    ttm3 = dttm_from_cell_array(cell_array)


    Lapmaptt = dttm_add(ttm1,dttm_add(ttm2,ttm3))


    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((DMA(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    Convtt1 = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((DMA(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    Convtt2 = dttm_from_cell_array(cell_array)

    cell_array%cells(1)%arr = reshape((Iden(It,It)),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    cell_array%cells(4)%arr = reshape((DMA(Ix,Ix)),[1,size(Ix),size(Ix),1])  
    Convtt3 = dttm_from_cell_array(cell_array)


    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),Ns,1]) 
    cell_array%cells(2)%arr = reshape((DMA(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    Convttmap1 = dttm_from_cell_array(cell_array)


    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),Ns,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(3)%arr = reshape((DMA(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(4)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    Convttmap2 = dttm_from_cell_array(cell_array)
    

    cell_array%cells(1)%arr = reshape((Iden(It,:)),[1,size(It),Ns,1]) 
    cell_array%cells(2)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(3)%arr = reshape((Iden(Ix,:)),[1,size(Ix),Ns,1])  
    cell_array%cells(4)%arr = reshape((DMA(Ix,:)),[1,size(Ix),Ns,1])  
    Convttmap3 = dttm_from_cell_array(cell_array)

    !!!!!

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

    Gbc = dtt_round(Gbc,tol)

    U0 = dtt_tensor_zeros([N1-1,N1-2,N1-2,N1-2]) 

    call F(U0,U)
    nrmFU0 = dtt_norm(U)

    maxiter = 100 

    du = U0
    
    epsk = 1.0d-1


    !!!!

    do k = 1,maxiter

    call Jfuntt(U0,opt1) 
    call JFbc(U0,opt2) 

    Jmatk = dttm_add(opt1,opt2)

    call F(U0,U)
    b_fortran = dtt_mul_d(U,-1.0d0)
    
    du = dtt_amen_solve(Jmatk, b_fortran, tol=epsk,nswp=30)

    alphas = [1.0d0,1.0d0/2.0d0,1.0d0/4.0d0,1.0d0/8.0d0,1.0d0/16.0d0]

    do t=1,5
    U1 = dtt_mul_d(du,alphas(t))
    U1 = dtt_add(U0,U1)
    U1 = dtt_round(U1,tol_=max(epsk,tol))
    call F(U1,U2)
    call F(U0,U)
    if (dtt_norm(dtt_round(U2,tol)) < dtt_norm(dtt_round(U,tol))) then
    exit
    end if
    end do

    
    local_err= dtt_norm( dtt_round(dtt_sub(U1,U0),tol)   )/dtt_norm(dtt_round(U0,tol))

    res_norm = dtt_norm(dtt_round(U2,tol))/nrmFU0

    epsk = min(min(local_err, res_norm), epsk)

    if ((local_err < Newtoneps) .or. (res_norm < Newtoneps)) then
    exit
    end if

    U0 = U1

    end do

    finish_time = MPI_Wtime()

    exacttt = dtt_tensor_ones(modes)
    call dtt_dmrgg(exacttt, exactfn_dmrg, accuracy=tol)


    temp_cell1 = core2cell_dtt(exacttt)
    allocate(temp_cell2(4))
    temp_cell2(1)%arr = temp_cell1(1)%arr(:,2:N1,:)
    do i = 2,4
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:N1-1,:)
    end do
    exacttt = dtt_from_cell_array(temp_cell2)
    deallocate(temp_cell2)

    

    Err = dtt_norm(dtt_round(dtt_sub(U1,exacttt)))/dtt_norm(dtt_round(exacttt))

    print *,'########'
    print *,'Error=',Err
    print *,'Time=',finish_time-start_time
    print *,'########'


    write(file_id, '(i0)') Ns
    filename_err   = '../plot_data/Fortran/Burger_err_' // trim(adjustl(file_id)) // '.dat'
    filename_time  = '../plot_data/Fortran/Burger_time_' // trim(adjustl(file_id)) // '.dat'

    open(unit=10, file=filename_err, status='replace', action='write')
    write(10,*) Err   
    close(10)
    open(unit=20, file=filename_time, status='replace', action='write')
    write(20,*) finish_time-start_time
    close(20)


    call mpi_finalize(info)
    

contains 





subroutine Jfuntt(Uin,opt) 
    use data
    use thor_lib, only: dtt_tensor,dtt_matrix,dttm_sub,dttm_add,dttm_from_dtt_tensor
    use matlab_struct_module, only: array3d
    use ttamen_lib, only: amen_mv,amen_mm

    implicit none 

    type(dtt_tensor),intent(in)  :: Uin
    type(dtt_matrix),intent(out)  :: opt
    type(array3d), allocatable :: G_temp(:)
    real(dp), allocatable :: tempten(:,:,:,:)
    integer :: i,j,k,l

    type(dtt_matrix) :: diagbf,diagdbf1,diagdbf2,diagdbf3,amen_mm_out1,amen_mm_out2,amen_mm_out3
    type(dtt_tensor) :: temp1,temp2,temp3

    !!!!
    G_temp = core2cell_dtt(Uin)
    do i=1,4
    allocate(  tempten(Uin%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),Uin%r(i))  )
    tempten = 0.0d0
      do j = 1,Uin%r(i-1)
        do k = 1,Uin%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    diagbf = dttm_from_cell_array(cell_array)
    !!!!

    call amen_mv(temp1, Convtt1, Uin, tol=tol)
    call amen_mv(temp2, Convtt2, Uin, tol=tol)
    call amen_mv(temp3, Convtt3, Uin, tol=tol)


    !!!!
    G_temp = core2cell_dtt(temp1)
    do i=1,4
    allocate(  tempten(temp1%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),temp1%r(i))  )
    tempten = 0.0d0
      do j = 1,temp1%r(i-1)
        do k = 1,temp1%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    diagdbf1 = dttm_from_cell_array(cell_array)
    !!!!
    !!!!
    G_temp = core2cell_dtt(temp2)
    do i=1,4
    allocate(  tempten(temp2%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),temp2%r(i))  )
    tempten = 0.0d0
      do j = 1,temp2%r(i-1)
        do k = 1,temp2%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    diagdbf2 = dttm_from_cell_array(cell_array)
    !!!!
    !!!!
    G_temp = core2cell_dtt(temp3)
    do i=1,4
    allocate(  tempten(temp3%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),temp3%r(i))  )
    tempten = 0.0d0
      do j = 1,temp3%r(i-1)
        do k = 1,temp3%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    diagdbf3 = dttm_from_cell_array(cell_array)
    !!!!

    call amen_mm(amen_mm_out1, diagbf, Convtt1, tol)
    call amen_mm(amen_mm_out2, diagbf, Convtt2, tol)
    call amen_mm(amen_mm_out3, diagbf, Convtt3, tol)

    

    amen_mm_out1 = dttm_add(amen_mm_out1,diagdbf1)
    amen_mm_out2 = dttm_add(amen_mm_out2,diagdbf2)
    amen_mm_out3 = dttm_add(amen_mm_out3,diagdbf3)

    opt = dttm_sub(Btt,Laptt)
    opt = dttm_sub(dttm_sub(dttm_sub(opt,amen_mm_out1),amen_mm_out2),amen_mm_out3)

    opt = dttm_from_dtt_tensor(dtt_round(opt,tol_=tol),opt%q(1:4),opt%s(1:4))

end subroutine Jfuntt


subroutine JFbc(Uin,opt) 
    use data
    use thor_lib, only: dtt_tensor,dtt_matrix,dtt_mult,dttm_add
    use ttamen_lib, only: amen_mv
    use matlab_struct_module, only: array3d

    implicit none 

    type(dtt_tensor),intent(in)  :: Uin
    type(dtt_matrix),intent(out)  :: opt

    type(dtt_tensor) :: temp1,temp2,temp3
    type(dtt_matrix) :: temp_m1,temp_m2,temp_m3
    type(array3d), allocatable :: G_temp(:)
    real(dp), allocatable :: tempten(:,:,:,:)
    integer :: i,j,k,l

    call amen_mv(temp1, Convttmap1, Gbc, tol=tol)
    call amen_mv(temp2, Convttmap2, Gbc, tol=tol)
    call amen_mv(temp3, Convttmap3, Gbc, tol=tol)

    temp1 = dtt_mult(Uin,temp1)
    temp2 = dtt_mult(Uin,temp2)
    temp3 = dtt_mult(Uin,temp3)

    !!!!
    G_temp = core2cell_dtt(temp1)
    do i=1,4
    allocate(  tempten(temp1%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),temp1%r(i))  )
    tempten = 0.0d0
      do j = 1,temp1%r(i-1)
        do k = 1,temp1%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    temp_m1 = dttm_from_cell_array(cell_array)
    !!!!
    !!!!
    G_temp = core2cell_dtt(temp2)
    do i=1,4
    allocate(  tempten(temp2%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),temp2%r(i))  )
    tempten = 0.0d0
      do j = 1,temp2%r(i-1)
        do k = 1,temp2%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    temp_m2 = dttm_from_cell_array(cell_array)
    !!!!
    !!!!
    G_temp = core2cell_dtt(temp3)
    do i=1,4
    allocate(  tempten(temp3%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),temp3%r(i))  )
    tempten = 0.0d0
      do j = 1,temp3%r(i-1)
        do k = 1,temp3%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    temp_m3 = dttm_from_cell_array(cell_array)
    !!!!

    opt = dttm_add(dttm_add (temp_m1,temp_m2),temp_m3)
    
    
end subroutine JFbc



subroutine F(Uin,Uout)
    use data
    use thor_lib, only: dtt_tensor,dtt_add,dtt_sub,dtt_round,dtt_mult
    use ttamen_lib, only: amen_mv

    implicit none

    type(dtt_tensor) :: temp1,temp2,temp3,temp4,temp5

    type(dtt_tensor),intent(in)  :: Uin
    type(dtt_tensor),intent(out) :: Uout

    call amen_mv(temp1, Convtt1, Uin, tol=tol)
    call amen_mv(temp2, Convtt2, Uin, tol=tol)
    call amen_mv(temp3, Convtt3, Uin, tol=tol)

    temp1 = dtt_round(dtt_mult(Uin,temp1),tol)
    temp2 = dtt_round(dtt_mult(Uin,temp2),tol)
    temp3 = dtt_round(dtt_mult(Uin,temp3),tol)

    call amen_mv(temp4, Btt, Uin, tol=tol)
    call amen_mv(temp5, Laptt, Uin, tol=tol)

    call Fbc(Uin,Uout)

    Uout = dtt_add(dtt_add(dtt_add(dtt_add(dtt_sub(temp4,temp5),temp1),temp2),temp3),Uout)

    Uout = dtt_round(Uout,tol)
    
end subroutine F


subroutine Fbc(Uin,Uout)
  use data
  use thor_lib, only: dtt_tensor,core2cell_dtt,dttm_from_cell_array,dtt_norm,dtt_round,dttm_sub,dttm_add,dttm_from_dtt_tensor
  use matlab_struct_module, only: array3d
  use ttamen_lib, only: amen_mm,amen_mv
  implicit none
  
  type(dtt_tensor),intent(in)  :: Uin
  type(dtt_tensor),intent(out) :: Uout

  integer :: i,j,k,l
  type(array3d), allocatable :: G_temp(:)
  real(dp), allocatable :: tempten(:,:,:,:)
  type(dtt_matrix) :: diagbf,amen_mm_out1,amen_mm_out2,amen_mm_out3,Amap


    G_temp = core2cell_dtt(Uin)
    do i=1,4
    allocate(  tempten(Uin%r(i-1),size(G_temp(i)%arr(j,:,k)),size(G_temp(i)%arr(j,:,k)),Uin%r(i))  )
    tempten = 0.0d0
      do j = 1,Uin%r(i-1)
        do k = 1,Uin%r(i)
            do l = 1,size(G_temp(i)%arr(j,:,k))
            tempten(j,l,l,k) = G_temp(i)%arr(j,l,k)
            end do
        end do
      end do
      cell_array%cells(i)%arr = tempten
      deallocate(tempten)
    end do
    diagbf = dttm_from_cell_array(cell_array)


  call amen_mm(amen_mm_out1, diagbf, Convttmap1, tol)
  call amen_mm(amen_mm_out2, diagbf, Convttmap2, tol)
  call amen_mm(amen_mm_out3, diagbf, Convttmap3, tol)

  Amap = dttm_add(dttm_add(dttm_add(dttm_sub(Bmaptt,Lapmaptt), amen_mm_out1),amen_mm_out2),amen_mm_out3)

  Amap = dttm_from_dtt_tensor(dtt_round(Amap,tol_=tol),Amap%q(1:4),Amap%s(1:4))

  call amen_mv(Uout, Amap, Gbc, tol=tol)

  Uout = dtt_round(Uout,tol)
    
end subroutine Fbc







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



function gfn(t, x, y, z) result(ans)
  use data
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans
  real(dp) :: pi 
  
  pi = 4 * atan (1.0d0)

  ans = ( pi * exp(t * (-3.289868133696453d0)) * &
            sin( (pi * (x + y + z)) / 3.0d0 ) * (2.0d0 / 3.0d0) ) / &
          ( exp(t * (-3.289868133696453d0)) * &
            cos( (pi * (x + y + z)) / 3.0d0 ) + 5.0d0 )
            
end function gfn


function gfn_dmrg(dpt, ind, nn) result(y)
    use data
    use thor_lib, only: tt_size
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
  use data
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans
  real(dp) :: pi,exparg,ang
  
  pi = 4 * atan (1.0d0)

  exparg = exp(t * (-3.289868133696453d0))
  
  ang    = (pi * (x + y + z)) / 3.0d0

  ans = (pi * exparg * sin(ang) * (2.0d0/3.0d0)) / (exparg * cos(ang) + 5.0d0)
            
end function exactfn




function exactfn_dmrg(dpt, ind, nn) result(y)
    use data
    use thor_lib, only: tt_size
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
