program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use thor_lib, only: dtt_tensor, dtt_matrix,dtt_tensor_ones,core2cell_dtt,dttm_from_cell_array,dttm_add,dtt_round,dtt_matrix_zeros,&
    dtt_from_cell_array,dtt_sub,dtt_norm,dttm_from_dtt_tensor,dtt_from_dttm,dtt_size
    use dmrgg_lib, only: dtt_dmrgg
    use ttamen_lib, only: amen_mv,dtt_amen_solve
    use matlab_struct_module, only: array3d,cell4d_array
    implicit none

    include 'mpif.h'
    
    integer :: d,jjns,Ex,Nx,i,j,k,info,nt1,nt2,j1,j2
    integer, dimension(5) :: NEs
    integer, allocatable :: modes(:)
    real(dp) :: Lx, Ly, Lz, hx, tt_tol,Err,norm1,norm2,start_time,finish_time
    real(dp), allocatable :: X_grid(:)
    real(dp), allocatable :: diag_CC(:,:,:,:), AG(:,:,:,:), MG(:,:,:,:)
    real(dp), allocatable :: CC(:,:)
    type(dtt_tensor) :: att,LLtt,F_newtt,uqtt,uexacttt,diff_tt
    type(array3d), allocatable :: G(:),temp_cell1(:),temp_cell2(:)
    real(dp), allocatable :: NNL(:,:), NNR(:,:),ML(:,:)
    real(dp), allocatable :: AA1(:,:), AA2(:,:), MM1(:,:), MM2(:,:)
    real(dp) :: A1(2,2), A2(2,2), M1(2,2), M2(2,2), B1(2,2), B2(2,2)
    real(dp) :: Bs(2,2)
    integer, allocatable :: Ix(:)
    type(cell4d_array) :: cell_array
    type(dtt_matrix) :: ttm1,ttm2,ttm3,ttm_f,MMtt,Agttcur,Agtt

    d = 3
    NEs(1) = 17 
    NEs(2) = 33 
    NEs(3) = 65 
    NEs(4) = 129 
    NEs(5) = 257 
    
    Lx = 1.0d0
    Ly = 1.0d0
    Lz = 1.0d0

    call mpi_init(info)
    
    do jjns = 1,size(NEs)

    start_time = MPI_Wtime()
    
    Ex = NEs(jjns)
    hx = Lx / real(Ex, dp)
    Nx = Ex + 1
    
    
    tt_tol = 0.01d0 * hx**2
    
    allocate(X_grid(Nx))
    do i = 1, Nx
        X_grid(i) = real(i - 1, dp) * hx
    end do
    
    allocate(modes(d))
    modes(:) = size(X_grid)

    att = dtt_tensor_ones(modes)
    call dtt_dmrgg(att, afn_dmrg, accuracy=tt_tol)

    call build_interpolation_matrix(NNL, Nx, Ex)

    NNR = transpose(NNL) 
    
    call nonlinear_mat(X_grid(1), X_grid(2), A1, A2, M1, M2, B1, B2)

    
    AA1 = kron_eye_block(A1, Ex)  
    AA2 = kron_eye_block(A2, Ex) 
    MM1 = kron_eye_block(M1, Ex)
    MM2 = kron_eye_block(M2, Ex)

    Bs = (hx / 6.0_dp) * reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2,2])    

    allocate(ML(Ex+1, Ex+1))
    ML = matmul(NNL, matmul(kron_eye_block(Bs, Ex), NNR))

    G = core2cell_dtt(att)
    nt1 = att%r(1)
    nt2 = att%r(2)

    allocate(Ix(Ex - 1))
    Ix = [(i, i = 2, Ex)]


    allocate(CC(d,2*Nx))
    allocate(diag_CC(d, 2, 2*Nx-2, 2*Nx-2))
    allocate(AG(d, 2, Nx, Nx))
    allocate(MG(d, 2, Nx, Nx))
    allocate(cell_array%cells(d))

    
    diag_CC = 0.0d0

    Agtt = dtt_matrix_zeros([size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
    do j1 = 1,nt1
        do j2=1,nt2
        
        do i = 1, Nx
        CC(1,2*i - 1) =G(1)%arr(1, i, j1)
        CC(1,2*i   ) = G(1)%arr(1, i, j1)

        CC(2,2*i - 1) =G(2)%arr(j1, i, j2)
        CC(2,2*i   ) = G(2)%arr(j1, i, j2)

        CC(3,2*i - 1) =G(3)%arr(j2, i, 1)
        CC(3,2*i   ) = G(3)%arr(j2, i, 1)
        end do
        
        do j=1,d
            do i = 1, 2*Nx-2
            diag_CC(j,1,i,i) = CC(j, i)
            end do

            do i = 3, 2*Nx
            diag_CC(j,2,i-2,i-2) = CC(j, i)
            end do
            
        end do


        do j =1,d
            AG(j,1,:,:) =  matmul(NNL,matmul(matmul(AA1,diag_CC(j,1,:,:)),NNR))
            AG(j,2,:,:) =  matmul(NNL,matmul(matmul(AA2,diag_CC(j,2,:,:)),NNR))
            MG(j,1,:,:) =  matmul(NNL,matmul(matmul(MM1,diag_CC(j,1,:,:)),NNR))
            MG(j,2,:,:) =  matmul(NNL,matmul(matmul(MM2,diag_CC(j,2,:,:)),NNR))
        end do

        Agttcur  = dtt_matrix_zeros([size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
        do i = 1,2
            do j = 1,2
              do k = 1,2
                  cell_array%cells(1)%arr = reshape(AG(1,i,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(2)%arr = reshape(MG(2,j,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(3)%arr = reshape(MG(3,k,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  ttm1 = dttm_from_cell_array(cell_array)
                  
                  cell_array%cells(1)%arr = reshape(MG(1,i,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(2)%arr = reshape(AG(2,j,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(3)%arr = reshape(MG(3,k,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  ttm2 = dttm_from_cell_array(cell_array)

                  cell_array%cells(1)%arr = reshape(MG(1,i,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(2)%arr = reshape(MG(2,j,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(3)%arr = reshape(AG(3,k,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  ttm3 = dttm_from_cell_array(cell_array)
                  
                  ttm_f = dttm_add(dttm_add(ttm1,ttm2), ttm3)

                  ttm_f = dttm_from_dtt_tensor(dtt_round(ttm_f,tol_=tt_tol),[size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
                  
                  Agttcur = dttm_add(Agttcur,ttm_f)

                  Agttcur = dttm_from_dtt_tensor(dtt_round(Agttcur,tol_=tt_tol),[size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
              
              end do 
            end do
        end do
        
        Agtt = dttm_add(Agtt,Agttcur)
        Agtt = dttm_from_dtt_tensor(dtt_round(Agtt,tol_=tt_tol),[size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
        
        end do
    end do

    LLtt = dtt_tensor_ones(modes)
    call dtt_dmrgg(LLtt, f_dmrg, accuracy=tt_tol)

    cell_array%cells(1)%arr = reshape(ML(Ix,:),[1,Nx,size(Ix),1])
    cell_array%cells(2)%arr = reshape(ML(Ix,:),[1,Nx,size(Ix),1])
    cell_array%cells(3)%arr = reshape(ML(Ix,:),[1,Nx,size(Ix),1])
    MMtt = dttm_from_cell_array(cell_array)

    F_newtt = dtt_tensor_ones([size(Ix),size(Ix),size(Ix)])
    call amen_mv(F_newtt, MMtt, LLtt, tol=tt_tol)

    uqtt = dtt_amen_solve(Agtt, F_newtt, tol=tt_tol)

    finish_time = MPI_Wtime()

    uexacttt = dtt_tensor_ones(modes)
    call dtt_dmrgg(uexacttt, uexactfn_dmrg, accuracy=tt_tol)
    temp_cell1 = core2cell_dtt(uexacttt)
    allocate(temp_cell2(d))
    do i = 1,d
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:Nx-1,:)
    end do
    uexacttt = dtt_from_cell_array(temp_cell2)

    diff_tt = dtt_sub(uqtt,uexacttt)

    
    norm1 = dtt_norm(dtt_round(diff_tt,tol_=tt_tol))
    norm2 = dtt_norm(dtt_round(uexacttt,tol_=tt_tol))
    Err = norm1/norm2
    print *,'Error=',Err
    print *,'Compression=',real(dtt_size(dtt_from_dttm(Agtt)),dp)/(real(size(Ix),dp)**2.0d0)**real(d,dp)
    print *,'Time=',finish_time-start_time
    print *,'########'
        
    deallocate(temp_cell2)
    deallocate(cell_array%cells)
    deallocate(CC)
    deallocate(AG)
    deallocate(MG)
    deallocate(diag_CC)
    deallocate(Ix)
    deallocate(ML)
    deallocate(modes)
    deallocate(X_grid)
    
    end do

    call mpi_finalize(info)


contains




function kron_eye_block(mat, Ex) result(global_mat)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  real(dp), intent(in) :: mat(:,:)
  integer, intent(in) :: Ex
  real(dp), allocatable :: global_mat(:,:)

  integer :: i, m, n
  m = size(mat,1)
  n = size(mat,2)

  allocate(global_mat(Ex*m, Ex*n))
  global_mat = 0.0d0

  do i = 0, Ex-1
    global_mat(1+i*m:i*m+m, 1+i*n:i*n+n) = mat
  end do
end function kron_eye_block



function phi1(x, x2, hx) result(val)
  real(dp), intent(in) :: x, x2, hx
  real(dp) :: val
  val = (x2 - x) / hx
end function phi1


function phi2(x, x1, hx) result(val)
  real(dp), intent(in) :: x, x1, hx
  real(dp) :: val
  val = (x - x1) / hx
end function phi2



subroutine nonlinear_mat(x1, x2, A1, A2, M1, M2, B1, B2)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  real(dp), intent(in) :: x1, x2
  real(dp), intent(out) :: A1(2,2), A2(2,2), M1(2,2), M2(2,2), B1(2,2), B2(2,2)

  real(dp) :: hx, x_md
  real(dp) :: DPhi1, DPhi2
  real(dp) :: p1, p2, pm
  real(dp) :: p11, p12, p1m
  real(dp) :: p21, p22, p2m

  hx = x2 - x1
  x_md = 0.5_dp * (x1 + x2)

  DPhi1 = -1.0_dp / hx
  DPhi2 = 1.0_dp / hx

  p11 = phi1(x1, x2, hx)
  p1m = phi1(x_md, x2, hx)
  p12 = phi1(x2, x2, hx)

  p21 = phi2(x1, x1, hx)
  p2m = phi2(x_md, x1, hx)
  p22 = phi2(x2, x1, hx)

  ! Mass matrices M1, M2
  M1(1,1) = (hx/6.0_dp) * (p11**3 + 4*p1m**3 + p12**3)
  M1(1,2) = (hx/6.0_dp) * (p11**2*p21 + 4*p1m**2*p2m + p12**2*p22)
  M1(2,1) = (hx/6.0_dp) * (p11*p21*p11 + 4*p1m*p2m*p1m + p12*p22*p12)
  M1(2,2) = (hx/6.0_dp) * (p11*p21**2 + 4*p1m*p2m**2 + p12*p22**2)

  M2(1,1) = (hx/6.0_dp) * (p21*p11**2 + 4*p2m*p1m**2 + p22*p12**2)
  M2(1,2) = (hx/6.0_dp) * (p21*p11*p21 + 4*p2m*p1m*p2m + p22*p12*p22)
  M2(2,1) = (hx/6.0_dp) * (p21**2*p11 + 4*p2m**2*p1m + p22**2*p12)
  M2(2,2) = (hx/6.0_dp) * (p21**3 + 4*p2m**3 + p22**3)

  ! Stiffness matrices A1, A2
  A1(1,1) = (hx/6.0_dp) * (p11*DPhi1*DPhi1 + 4*p1m*DPhi1*DPhi1 + p12*DPhi1*DPhi1)
  A1(1,2) = (hx/6.0_dp) * (p11*DPhi1*DPhi2 + 4*p1m*DPhi1*DPhi2 + p12*DPhi1*DPhi2)
  A1(2,1) = (hx/6.0_dp) * (p11*DPhi2*DPhi1 + 4*p1m*DPhi2*DPhi1 + p12*DPhi2*DPhi1)
  A1(2,2) = (hx/6.0_dp) * (p11*DPhi2*DPhi2 + 4*p1m*DPhi2*DPhi2 + p12*DPhi2*DPhi2)

  A2(1,1) = (hx/6.0_dp) * (p21*DPhi1*DPhi1 + 4*p2m*DPhi1*DPhi1 + p22*DPhi1*DPhi1)
  A2(1,2) = (hx/6.0_dp) * (p21*DPhi1*DPhi2 + 4*p2m*DPhi1*DPhi2 + p22*DPhi1*DPhi2)
  A2(2,1) = (hx/6.0_dp) * (p21*DPhi2*DPhi1 + 4*p2m*DPhi2*DPhi1 + p22*DPhi2*DPhi1)
  A2(2,2) = (hx/6.0_dp) * (p21*DPhi2*DPhi2 + 4*p2m*DPhi2*DPhi2 + p22*DPhi2*DPhi2)

  ! Nonlinear matrices B1, B2
  B1(1,1) = (hx/6.0_dp) * (p11*DPhi1*p11 + 4*p1m*DPhi1*p1m + p12*DPhi1*p12)
  B1(1,2) = (hx/6.0_dp) * (p11*DPhi2*p11 + 4*p1m*DPhi2*p1m + p12*DPhi2*p12)
  B1(2,1) = (hx/6.0_dp) * (p11*DPhi1*p21 + 4*p1m*DPhi1*p2m + p12*DPhi1*p22)
  B1(2,2) = (hx/6.0_dp) * (p11*DPhi2*p21 + 4*p1m*DPhi2*p2m + p12*DPhi2*p22)

  B2(1,1) = (hx/6.0_dp) * (p21*DPhi1*p11 + 4*p2m*DPhi1*p1m + p22*DPhi1*p12)
  B2(1,2) = (hx/6.0_dp) * (p21*DPhi2*p11 + 4*p2m*DPhi2*p1m + p22*DPhi2*p12)
  B2(2,1) = (hx/6.0_dp) * (p21*DPhi1*p21 + 4*p2m*DPhi1*p2m + p22*DPhi1*p22)
  B2(2,2) = (hx/6.0_dp) * (p21*DPhi2*p21 + 4*p2m*DPhi2*p2m + p22*DPhi2*p22)

end subroutine nonlinear_mat



subroutine build_interpolation_matrix(NNL, Nx, Ex)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: Nx, Ex
  real(dp), allocatable, intent(out) :: NNL(:,:)
  integer :: i

  allocate(NNL(Nx, 2*Ex))
  NNL = 0.0d0

  ! Interior rows
  do i = 2, Ex
    NNL(i, 2*(i-1)) = 1.0d0  
    NNL(i, 2*i-1)   = 1.0d0  
  end do

  ! Boundary rows
  NNL(1,1) = 1.0d0
  NNL(Nx, 2*Ex) = 1.0d0
end subroutine build_interpolation_matrix


pure function uexactfn(x, y, z) result(u)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  real(dp), intent(in) :: x, y, z
  real(dp) :: pi,u

  pi = 4 * atan (1.0d0)
  
  u = sin(pi*x) * sin(pi*y) * sin(pi*z)
  
end function uexactfn


pure function uexactfn_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)

    do i = 1,d
        x(i) = X_grid(ind(i))
    end do

    y = uexactfn(x(1),x(2),x(3))
    
end function uexactfn_dmrg


pure function afn(x, y, z) result(a)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  real(dp), intent(in) :: x, y, z
  real(dp) :: pi,a
  
  pi = 4 * atan (1.0d0)
  
  a = 1.0d0 + cos(pi*(x + y)) * cos(pi*z)
  
end function afn


pure function afn_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)

    do i = 1,d
        x(i) = X_grid(ind(i))
    end do

    y = afn(x(1),x(2),x(3))
    
end function afn_dmrg


pure function f(x, y, z) result(fval)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  real(dp), intent(in) :: x, y, z
  real(dp) :: fval, u, a
  real(dp) :: grad_dot, pi

  pi = 4 * atan (1.0d0)

  ! u and a
  u = uexactfn(x,y,z) 
  a = afn(x,y,z) 

  ! grad(a) · grad(u)
  grad_dot = -pi**2 * sin(pi*(x + y)) * cos(pi*z) * ( &
               cos(pi*x) * sin(pi*y) + sin(pi*x) * cos(pi*y) ) * sin(pi*z) &
             - pi**2 * cos(pi*(x + y)) * sin(pi*z) * sin(pi*x) * sin(pi*y) * cos(pi*z)

  ! Final expression
  fval = -grad_dot + 3.0d0 * pi**2 * a * u
end function f


pure function f_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x(tt_size)

    do i = 1,d
        x(i) = X_grid(ind(i))
    end do

    y = f(x(1),x(2),x(3))
    
end function f_dmrg


end program 
    
