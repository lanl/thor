program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use thor_lib, only: dttm_from_cell_array,dtt_matrix,dtt_norm,dtt_tensor_ones,dtt_tensor,core2cell_dtt,dtt_matrix_zeros,&
    dttm_add,dtt_round,dttm_from_dtt_tensor,dttm_kron,dtt_from_cell_array,dtt_sub,dtt_size,dtt_from_dttm
    use dmrgg_lib, only: dtt_dmrgg
    use ttamen_lib, only: amen_mv,dtt_amen_solve
    use matlab_struct_module, only: cell4d_array,array3d

    implicit none

    include 'mpif.h'
    
    integer :: d,jjns,Ex,Nx,Et,Nt,i,info,nt1,nt2,j1,j2,j,k
    integer :: NE
    real(dp) :: Lx,Ly,Lz,Lt,hx,ht,tt_tol,Err,norm1,norm2,start_time,finish_time
    real(dp), allocatable :: X(:),T(:),NNL(:,:),NNR(:,:),NNLt(:,:),NNRt(:,:),ML(:,:)
    real(dp), allocatable :: CC(:,:),diag_CC(:,:,:,:), AG(:,:,:,:), MG(:,:,:,:),AB(:,:,:,:)
    real(dp) :: A1(2,2), A2(2,2), M1(2,2), M2(2,2), B1(2,2), B2(2,2),Bs(2,2)
    real(dp) :: A1t(2,2), A2t(2,2), M1t(2,2), M2t(2,2), B1t(2,2), B2t(2,2),Bt(2,2)
    real(dp) :: MT(2,2)
    real(dp), allocatable :: AA1(:,:), AA2(:,:),MM1(:,:), MM2(:,:),BB1(:,:), BB2(:,:)
    real(dp), allocatable :: ZT(:,:),MLt(:,:),MX1(:,:),MX2(:,:)
    integer, allocatable :: Ix(:),It(:)
    type(cell4d_array) :: cell_array
    type(dtt_matrix) :: AT,Agtt,ttm1,ttm2,ttm3,ttm_f,Agttcur,AD,ACxtt,ACytt,ACztt,ACspace,AC,AR,Att_global,MMxtt,MMtt
    type(dtt_tensor) :: att,btt,ctt,LLtt,F_newtt,utt,uexacttt,diff_tt
    integer, allocatable :: modes(:)
    type(array3d), allocatable :: G(:),temp_cell1(:),temp_cell2(:)
    character(len=25) :: arg

    d = 3
    
    call get_command_argument(1, arg)
    read(arg, *) NE
   

    Lx = 1.0d0
    Ly = 1.0d0
    Lz = 1.0d0
    Lt = 1.0d0

    call mpi_init(info)

    !do jjns = 1,size(NEs)

    start_time = MPI_Wtime()

    Ex = NE 
    hx = Lx / real(Ex, dp)
    Nx = Ex + 1
    Et = Ex-1
    ht = Lt / real(Et, dp)
    Nt = Et + 1

    tt_tol = 0.01 * hx**2


    allocate(X(Nx))
    allocate(T(Nt))
    do i = 1, Nx
        X(i) = real(i - 1, dp) * hx
    end do
    do i = 1, Nt
        T(i) = real(i - 1, dp) * ht
    end do

    allocate(modes(d+1))
    modes(1)  = size(T)
    modes(2:) = size(X)


    call build_interpolation_matrix(NNL , Nx, Ex)
    call build_interpolation_matrix(NNLt, Nt, Et)

    NNR = transpose(NNL) 
    NNRt = transpose(NNLt) 

    call nonlinear_mat(X(1), X(2), A1, A2, M1, M2, B1, B2)

    AA1 = kron_eye_block(A1, Ex)  
    AA2 = kron_eye_block(A2, Ex)
    MM1 = kron_eye_block(M1, Ex)  
    MM2 = kron_eye_block(M2, Ex) 
    BB1 = kron_eye_block(B1, Ex)  
    BB2 = kron_eye_block(B2, Ex) 

    Bs = (hx / 6.0_dp) * reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2,2])    

    allocate(ML(Ex+1, Ex+1))
    ML = matmul(NNL, matmul(kron_eye_block(Bs, Ex), NNR))


    call nonlinear_mat(T(1), T(2), A1t, A2t, M1t, M2t, B1t, B2t)

    Bt = (ht / 6.0_dp) * reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2,2])    

    call Mat_Time(T(1:2), MT)

    allocate(ZT(Et+1, Et+1))
    allocate(MLt(Et+1, Et+1))
    allocate(MX1(Et+1, Et+1))
    allocate(MX2(Et+1, Et+1))
    
    ZT  = matmul(NNLt, matmul(kron_eye_block(MT, Et), NNRt))
    MLt = matmul(NNLt, matmul(kron_eye_block(Bt, Et), NNRt))
    MX1 = matmul(NNLt, matmul(kron_eye_block(M1t, Et), NNRt))
    MX2 = matmul(NNLt, matmul(kron_eye_block(M2t, Et), NNRt))


    allocate(Ix(Ex - 1))
    Ix = [(i, i = 2, Ex)]

    allocate(It(Et))
    It = [(i, i = 2, Et+1)]


    allocate(cell_array%cells(4))
    cell_array%cells(1)%arr = reshape(ZT(It,It),[1,size(It),size(It),1]) 
    cell_array%cells(2)%arr = reshape(ML(Ix,Ix),[1,size(Ix),size(Ix),1]) 
    cell_array%cells(3)%arr = reshape(ML(Ix,Ix),[1,size(Ix),size(Ix),1]) 
    cell_array%cells(4)%arr = reshape(ML(Ix,Ix),[1,size(Ix),size(Ix),1]) 
    AT = dttm_from_cell_array(cell_array)
    deallocate(cell_array%cells)


    att = dtt_tensor_ones(modes(2:))
    call dtt_dmrgg(att, afn_dmrg, accuracy=tt_tol)
    
    G = core2cell_dtt(att)
    nt1 = att%r(1)
    nt2 = att%r(2)

    allocate(CC(d,2*Nx))
    allocate(diag_CC(d, 2, 2*Nx-2, 2*Nx-2))
    allocate(AG(d, 2, Nx, Nx))
    allocate(MG(d, 2, Nx, Nx))
    allocate(cell_array%cells(d))
    

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
    deallocate(cell_array%cells)
    deallocate(AG)



    allocate(cell_array%cells(1))
    cell_array%cells(1)%arr = reshape(MX1(It,It),[1,size(It),size(It),1]) 
    AD = dttm_kron(dttm_from_cell_array(cell_array),Agtt)
    cell_array%cells(1)%arr = reshape(MX2(It,It),[1,size(It),size(It),1]) 
    AD = dttm_add(AD , dttm_kron(dttm_from_cell_array(cell_array),Agtt))
    AD = dttm_from_dtt_tensor(dtt_round(AD,tol_=tt_tol),[size(Ix),size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix),size(Ix)])
    deallocate(cell_array%cells)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Build Convection operator - x term !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    btt = dtt_tensor_ones(modes(2:))
    call dtt_dmrgg(btt, bfn_1_dmrg, accuracy=tt_tol)

    G = core2cell_dtt(btt)
    nt1 = btt%r(1)
    nt2 = btt%r(2)


    allocate(AB(d, 2, Nx, Nx))
    allocate(cell_array%cells(d))

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
            AB(j,1,:,:) =  matmul(NNL,matmul(matmul(BB1,diag_CC(j,1,:,:)),NNR))
            AB(j,2,:,:) =  matmul(NNL,matmul(matmul(BB2,diag_CC(j,2,:,:)),NNR))
            MG(j,1,:,:) =  matmul(NNL,matmul(matmul(MM1,diag_CC(j,1,:,:)),NNR))
            MG(j,2,:,:) =  matmul(NNL,matmul(matmul(MM2,diag_CC(j,2,:,:)),NNR))
        end do

        Agttcur  = dtt_matrix_zeros([size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
        do i = 1,2
            do j = 1,2
              do k = 1,2
                  cell_array%cells(1)%arr = reshape(AB(1,i,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(2)%arr = reshape(MG(2,j,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(3)%arr = reshape(MG(3,k,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  ttm_f = dttm_from_cell_array(cell_array)
                  
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

    ACxtt = Agtt


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Build Convection operator - y term !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    btt = dtt_tensor_ones(modes(2:))
    call dtt_dmrgg(btt, bfn_2_dmrg, accuracy=tt_tol)

    G = core2cell_dtt(btt)
    nt1 = btt%r(1)
    nt2 = btt%r(2)

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
            AB(j,1,:,:) =  matmul(NNL,matmul(matmul(BB1,diag_CC(j,1,:,:)),NNR))
            AB(j,2,:,:) =  matmul(NNL,matmul(matmul(BB2,diag_CC(j,2,:,:)),NNR))
            MG(j,1,:,:) =  matmul(NNL,matmul(matmul(MM1,diag_CC(j,1,:,:)),NNR))
            MG(j,2,:,:) =  matmul(NNL,matmul(matmul(MM2,diag_CC(j,2,:,:)),NNR))
        end do

        Agttcur  = dtt_matrix_zeros([size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
        do i = 1,2
            do j = 1,2
              do k = 1,2
                  cell_array%cells(1)%arr = reshape(MG(1,i,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(2)%arr = reshape(AB(2,j,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(3)%arr = reshape(MG(3,k,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  ttm_f = dttm_from_cell_array(cell_array)
                  
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

    ACytt = Agtt

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Build Convection operator - z term !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    btt = dtt_tensor_ones(modes(2:))
    call dtt_dmrgg(btt, bfn_3_dmrg, accuracy=tt_tol)

    G = core2cell_dtt(btt)
    nt1 = btt%r(1)
    nt2 = btt%r(2)

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
            AB(j,1,:,:) =  matmul(NNL,matmul(matmul(BB1,diag_CC(j,1,:,:)),NNR))
            AB(j,2,:,:) =  matmul(NNL,matmul(matmul(BB2,diag_CC(j,2,:,:)),NNR))
            MG(j,1,:,:) =  matmul(NNL,matmul(matmul(MM1,diag_CC(j,1,:,:)),NNR))
            MG(j,2,:,:) =  matmul(NNL,matmul(matmul(MM2,diag_CC(j,2,:,:)),NNR))
        end do

        Agttcur  = dtt_matrix_zeros([size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
        do i = 1,2
            do j = 1,2
              do k = 1,2
                  cell_array%cells(1)%arr = reshape(MG(1,i,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(2)%arr = reshape(MG(2,j,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(3)%arr = reshape(AB(3,k,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  ttm_f = dttm_from_cell_array(cell_array)
                  
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
    deallocate(cell_array%cells)
    deallocate(AB)

    ACztt = Agtt

    ACspace = dttm_add(ACxtt,dttm_add( ACytt , ACztt))

    

    allocate(cell_array%cells(1))
    cell_array%cells(1)%arr = reshape(MX1(It,It),[1,size(It),size(It),1]) 
    AC = dttm_kron(dttm_from_cell_array(cell_array),ACspace)
    cell_array%cells(1)%arr = reshape(MX2(It,It),[1,size(It),size(It),1]) 
    AC = dttm_add(AC , dttm_kron(dttm_from_cell_array(cell_array),ACspace))
    AC = dttm_from_dtt_tensor(dtt_round(AC,tol_=tt_tol),[size(Ix),size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix),size(Ix)])
    deallocate(cell_array%cells)



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Build Reaction operator !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ctt = dtt_tensor_ones(modes(2:))
    call dtt_dmrgg(ctt, cfn_dmrg, accuracy=tt_tol)

    G = core2cell_dtt(ctt)
    nt1 = btt%r(1)
    nt2 = btt%r(2)

    allocate(cell_array%cells(d))

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
            MG(j,1,:,:) =  matmul(NNL,matmul(matmul(MM1,diag_CC(j,1,:,:)),NNR))
            MG(j,2,:,:) =  matmul(NNL,matmul(matmul(MM2,diag_CC(j,2,:,:)),NNR))
        end do

        Agttcur  = dtt_matrix_zeros([size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix)])
        do i = 1,2
            do j = 1,2
              do k = 1,2
                  cell_array%cells(1)%arr = reshape(MG(1,i,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(2)%arr = reshape(MG(2,j,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  cell_array%cells(3)%arr = reshape(MG(3,k,Ix, Ix),[1,size(Ix),size(Ix),1]) 
                  ttm_f = dttm_from_cell_array(cell_array)
                  
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
    deallocate(cell_array%cells)
    deallocate(CC)
    deallocate(diag_CC)
    deallocate(MG)


    allocate(cell_array%cells(1))
    cell_array%cells(1)%arr = reshape(MX1(It,It),[1,size(It),size(It),1]) 
    AR = dttm_kron(dttm_from_cell_array(cell_array),Agtt)
    cell_array%cells(1)%arr = reshape(MX2(It,It),[1,size(It),size(It),1]) 
    AR = dttm_add(AR , dttm_kron(dttm_from_cell_array(cell_array),Agtt))
    AR = dttm_from_dtt_tensor(dtt_round(AR,tol_=tt_tol),[size(Ix),size(Ix),size(Ix),size(Ix)],[size(Ix),size(Ix),size(Ix),size(Ix)])
    deallocate(cell_array%cells)


    !!!!!!!!!!!!!!!! Build global !!!!!!!!!!!!!!!!

    Att_global = dttm_add(AT,dttm_add(AC,dttm_add(AD, AR)))

    LLtt = dtt_tensor_ones(modes)
    call dtt_dmrgg(LLtt, f_dmrg, accuracy=tt_tol)


    allocate(cell_array%cells(d))
    cell_array%cells(1)%arr = reshape(ML(Ix,:),[1,Nx,size(Ix),1])
    cell_array%cells(2)%arr = reshape(ML(Ix,:),[1,Nx,size(Ix),1])
    cell_array%cells(3)%arr = reshape(ML(Ix,:),[1,Nx,size(Ix),1])
    MMxtt = dttm_from_cell_array(cell_array)
    deallocate(cell_array%cells)




    allocate(cell_array%cells(1))
    cell_array%cells(1)%arr = reshape(MLt(It,:),[1,Nt,size(It),1]) 
    MMtt = dttm_kron(dttm_from_cell_array(cell_array),MMxtt)
    deallocate(cell_array%cells)


    F_newtt = dtt_tensor_ones([size(Ix),size(Ix),size(Ix)])
    call amen_mv(F_newtt, MMtt, LLtt, tol=tt_tol)

    utt = dtt_amen_solve(Att_global, F_newtt, tol=tt_tol)

    finish_time = MPI_Wtime()

    uexacttt = dtt_tensor_ones(modes)
    call dtt_dmrgg(uexacttt, uexactfn_dmrg, accuracy=tt_tol)
    temp_cell1 = core2cell_dtt(uexacttt)
    allocate(temp_cell2(d+1))
    temp_cell2(1)%arr = temp_cell1(1)%arr(:,2:Nt,:)
    do i = 2,d+1
        temp_cell2(i)%arr = temp_cell1(i)%arr(:,2:Nx-1,:)
    end do
    uexacttt = dtt_from_cell_array(temp_cell2)


    diff_tt = dtt_sub(utt,uexacttt)    
    norm1 = dtt_norm(dtt_round(diff_tt,tol_=tt_tol))
    norm2 = dtt_norm(dtt_round(uexacttt,tol_=tt_tol))
    Err = norm1/norm2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print *,'Error=',Err
    print *,'Compression=',real(dtt_size(dtt_from_dttm(Att_global)),dp) / product( real( Att_global%q(:d+1) * Att_global%s(:d+1),dp ))
    print *,'Time=',finish_time-start_time
    print *,'########'
    
    deallocate(temp_cell2)
    deallocate(Ix)
    deallocate(It)
    deallocate(MX2)
    deallocate(MX1)
    deallocate(MLt)
    deallocate(ZT)
    deallocate(ML)
    deallocate(X)
    deallocate(T)
    deallocate(modes)

    !end do

    call mpi_finalize(info)
    

contains


function Phit1(t, t1, t2, ht) result(fy)
  implicit none
  real(8), intent(in) :: t, t1, t2, ht
  real(8) :: fy
  fy = (t2 - t) / ht
end function Phit1


function Phit2(t, t1, t2, ht) result(fy)
  implicit none
  real(8), intent(in) :: t, t1, t2, ht
  real(8) :: fy
  fy = (t - t1) / ht
end function Phit2



subroutine Mat_Time(T_Loc, MT)
  implicit none
  real(dp), intent(in)  :: T_Loc(2)         
  real(dp), intent(out) :: MT(2,2)          
  real(dp) :: t1, t2, ht, t_md
  real(dp) :: DPhi1, DPhi2
  integer :: i, j


  t1 = T_Loc(1)
  t2 = T_Loc(2)
  ht = t2 - t1
  t_md = 0.5d0 * (t1 + t2)

  DPhi1 = -1.0d0 / ht
  DPhi2 =  1.0d0 / ht

  MT(1,1) = (ht / 6.0d0) * (DPhi1 * Phit1(t1, t1, t2, ht) + &
                            4.0d0 * DPhi1 * Phit1(t_md, t1, t2, ht) + &
                            DPhi1 * Phit1(t2, t1, t2, ht))

  MT(1,2) = (ht / 6.0d0) * (DPhi2 * Phit1(t1, t1, t2, ht) + &
                            4.0d0 * DPhi2 * Phit1(t_md, t1, t2, ht) + &
                            DPhi2 * Phit1(t2, t1, t2, ht))

  MT(2,1) = (ht / 6.0d0) * (DPhi1 * Phit2(t1, t1, t2, ht) + &
                            4.0d0 * DPhi1 * Phit2(t_md, t1, t2, ht) + &
                            DPhi1 * Phit2(t2, t1, t2, ht))

  MT(2,2) = (ht / 6.0d0) * (DPhi2 * Phit2(t1, t1, t2, ht) + &
                            4.0d0 * DPhi2 * Phit2(t_md, t1, t2, ht) + &
                            DPhi2 * Phit2(t2, t1, t2, ht))

end subroutine Mat_Time




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


subroutine build_interpolation_matrix(NNL_xt, N, E)
  implicit none
  integer, intent(in) :: N, E
  real(dp), allocatable, intent(out) :: NNL_xt(:,:)
  integer :: i

  allocate(NNL_xt(N, 2*E))
  NNL_xt = 0.0d0

  do i = 2, E
    NNL_xt(i, 2*(i-1)) = 1.0d0  
    NNL_xt(i, 2*i-1)   = 1.0d0  
  end do

  NNL_xt(1,1) = 1.0d0
  NNL_xt(N, 2*E) = 1.0d0
end subroutine build_interpolation_matrix



function uexactfn(t, x, y, z) result(u)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: u
  real(dp) :: pi

  pi = 4 * atan (1.0d0)

  u = sin(pi*x) * sin(pi*y) * sin(pi*z) * sin(pi*t) + &
      sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) * sin(2*pi*t) + &
      sin(3*pi*x) * sin(3*pi*y) * sin(3*pi*z) * sin(3*pi*t)
end function




function uexactfn_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = T(ind(1))

    do i = 2,dpt
        x_temp(i) = X(ind(i))
    end do

    y = uexactfn(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function uexactfn_dmrg





function afn(x, y, z) result(a)
  implicit none
  real(dp), intent(in) :: x, y, z
  real(dp) :: a
  real(dp) :: pi

  pi = 4 * atan (1.0d0)

  a = cos(pi*x) * cos(pi*y) * cos(pi*z) + 1.0d0
end function



function afn_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    do i = 1,d
        x_temp(i) = X(ind(i))
    end do

    y = afn(x_temp(1),x_temp(2),x_temp(3))
    
end function afn_dmrg


function bfn_1_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    do i = 1,d
        x_temp(i) = X(ind(i))
    end do

    y = x_temp(1)
    
end function bfn_1_dmrg


function bfn_2_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    do i = 1,d
        x_temp(i) = X(ind(i))
    end do

    y = x_temp(2)
    
end function bfn_2_dmrg


function bfn_3_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    do i = 1,d
        x_temp(i) = X(ind(i))
    end do

    y = x_temp(3)
    
end function bfn_3_dmrg



function cfn(x, y, z) result(c)
  implicit none
  real(dp), intent(in) :: x, y, z
  real(dp) :: c

  c = exp(-x - y - z)
end function



function cfn_dmrg(d, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: d, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    do i = 1,d
        x_temp(i) = X(ind(i))
    end do

    y = cfn(x_temp(1),x_temp(2),x_temp(3))
    
end function cfn_dmrg




function dudt(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans
  real(dp), parameter :: pi = 4 * atan (1.0d0)

  ans =  pi   * sin(pi*x)   * sin(pi*y)   * sin(pi*z)   * cos(pi*t) + &
         2*pi * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) * cos(2*pi*t) + &
         3*pi * sin(3*pi*x) * sin(3*pi*y) * sin(3*pi*z) * cos(3*pi*t)
end function dudt




function Dfn(t, x, y, z) result(diff_term)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: diff_term
  real(dp), parameter :: pi = 4 * atan (1.0d0)
  real(dp) :: a, ax, ay, az, ux,uy,uz,d2udx2,d2udy2,d2udz2

  a =   afn(x, y, z)
  ax = -pi*sin(pi*x)*cos(pi*y)*cos(pi*z)
  ay = -pi*cos(pi*x)*sin(pi*y)*cos(pi*z)
  az = -pi*cos(pi*x)*cos(pi*y)*sin(pi*z)

  ux = pi*cos(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*t) + &
       2*pi*cos(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)*sin(2*pi*t) + &
       3*pi*cos(3*pi*x)*sin(3*pi*y)*sin(3*pi*z)*sin(3*pi*t)

  uy = pi*sin(pi*x)*cos(pi*y)*sin(pi*z)*sin(pi*t) + &
       2*pi*sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*z)*sin(2*pi*t) + &
       3*pi*sin(3*pi*x)*cos(3*pi*y)*sin(3*pi*z)*sin(3*pi*t)

  uz = pi*sin(pi*x)*sin(pi*y)*cos(pi*z)*sin(pi*t) + &
       2*pi*sin(2*pi*x)*sin(2*pi*y)*cos(2*pi*z)*sin(2*pi*t) + &
       3*pi*sin(3*pi*x)*sin(3*pi*y)*cos(3*pi*z)*sin(3*pi*t)


  d2udx2 =  -pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*t) &
            -4*pi**2*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)*sin(2*pi*t) &
            -9*pi**2*sin(3*pi*x)*sin(3*pi*y)*sin(3*pi*z)*sin(3*pi*t)

  d2udy2 =  -pi**2 * sin(pi*x) * sin(pi*y) * sin(pi*z) * sin(pi*t) &
            -4*pi**2 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) * sin(2*pi*t) &
            -9*pi**2 * sin(3*pi*x) * sin(3*pi*y) * sin(3*pi*z) * sin(3*pi*t)

  d2udz2 =  -pi**2 * sin(pi*x) * sin(pi*y) * sin(pi*z) * sin(pi*t) &
            -4*pi**2 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) * sin(2*pi*t) &
            -9*pi**2 * sin(3*pi*x) * sin(3*pi*y) * sin(3*pi*z) * sin(3*pi*t)

  diff_term = -((ax * ux + a * d2udx2) + &
                (ay * uy + a * d2udy2) + &
                (az * uz + a * d2udz2))

end function Dfn




function Rfn(t, x, y, z) result(rhs)
  implicit none
  real(8), intent(in) :: t, x, y, z
  real(8) :: rhs

  rhs = cfn(x, y, z) * uexactfn(t, x, y, z)
end function Rfn



function Cvecfn(t, x, y, z) result(convection_term)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: convection_term
  real(dp), parameter :: pi = 4 * atan (1.0d0)
  real(dp) :: ux, uy, uz

  ux = pi   * cos(pi*x)   * sin(pi*y)   * sin(pi*z)   * sin(pi*t) + &
       2*pi * cos(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) * sin(2*pi*t) + &
       3*pi * cos(3*pi*x) * sin(3*pi*y) * sin(3*pi*z) * sin(3*pi*t)

  uy = pi   * sin(pi*x)   * cos(pi*y)   * sin(pi*z)   * sin(pi*t) + &
       2*pi * sin(2*pi*x) * cos(2*pi*y) * sin(2*pi*z) * sin(2*pi*t) + &
       3*pi * sin(3*pi*x) * cos(3*pi*y) * sin(3*pi*z) * sin(3*pi*t)

  uz = pi   * sin(pi*x)   * sin(pi*y)   * cos(pi*z)   * sin(pi*t) + &
       2*pi * sin(2*pi*x) * sin(2*pi*y) * cos(2*pi*z) * sin(2*pi*t) + &
       3*pi * sin(3*pi*x) * sin(3*pi*y) * cos(3*pi*z) * sin(3*pi*t)

  convection_term = x * ux + y * uy + z * uz
end function Cvecfn


function f(t, x, y, z) result(ans)
  implicit none
  real(dp), intent(in) :: t, x, y, z
  real(dp) :: ans

  ans = dudt(t,x,y,z) + Dfn(t,x,y,z) + Cvecfn(t,x,y,z) + Rfn(t,x,y,z)
end function f



function f_dmrg(dpt, ind, nn) result(y)
    use thor_lib, only: tt_size
    implicit none
    real(dp) :: y
    integer, intent(IN) :: dpt, ind(1:tt_size), nn(1:tt_size)
    integer :: i
    real(dp) :: x_temp(tt_size)

    x_temp(1) = T(ind(1))

    do i = 2,dpt
        x_temp(i) = X(ind(i))
    end do

    y = f(x_temp(1),x_temp(2),x_temp(3),x_temp(4))
    
end function f_dmrg





end program 
