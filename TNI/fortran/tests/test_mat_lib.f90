!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                     LANL 2024/01/19
!!                                           T-DIV
!!                Ismael Djibrilla Boureima [IDB], Oleg Korobkin [OK], Ryan Barron [RCB]
!!
!!_________________________________________________________________________________________
!!   SCOPE: testing thor utils/mat.f90 (module mat_lib)
!!_________________________________________________________________________________________
!! Fun facts: 1. svd routine in mat_lib returns u in a transposed state,
!!               relative to the Matlab function
!!            2. svd determines eigenvectors up to a sign (plus or minus)
program main
  use mat_lib, only: d_svd, matinv_lu_d, matinv, norm => normfro, d_chol, &
                     diag => matlab_diag, eye => matlab_eye, qr => matlab_qr, &
                     svdgram
  use lr_lib, only: d2_lual                   
  use matrix_util, only: unravel_index, ravel_multi_index
  use matlab_struct_module, only : pprint_matrix, pprint_matrix3d, array3d
  use thor_lib, only: dtt_tensor, dtt_matrix, pprint, dtt_tensor_ones, core2cell, &
                     operator(*), dtt_random_ortho
  implicit none
  double precision, allocatable :: A(:,:), A_inv(:,:), g(:), col(:,:,:)
  double precision, allocatable :: R2(:,:), R2_g(:,:), R2_inv(:,:)
  double precision, pointer :: u(:,:), v(:,:), s(:)
  double precision, allocatable :: ug(:,:), vg(:,:), sg(:), C(:,:), R(:,:)
  double precision, allocatable :: Ucg(:,:), Uc(:,:)
  double precision, parameter :: err_tol = 1d-15
  double precision :: err
  integer :: i
  character(:), allocatable :: strs(:)
  type(dtt_tensor) :: ta, tx, ty
  type(dtt_matrix) :: mA
  type(array3d), allocatable :: tz_c(:), ty0_c(:)
  integer :: failed
  integer, allocatable :: mind(:), ishp(:)


  print '("")'
  print '("testing thor utils/mat.f90 (module mat_lib)")'
  print '("")'
  failed = 0

  allocate(A(50, 45))
  call random_number(A)
  A = A - 0.5d0
  call pprint_matrix(A,strings_=strs,line_len_=44,max_rows_=5,frmt_="(ES9.2)")
  do i=1,size(strs)
     if (i.eq.3) then
        print '("A = |",A, "|")', strs(i)
     else
        print '("    |",A, "|")', strs(i)
     endif
  enddo
  deallocate(strs)
  call pprint_matrix(reshape(A, [size(A), 1]))
  call pprint_matrix(reshape(A, [1, size(A)]))
  !call pprint_matrix(R2, strings_=strs, &
  !                         line_len_=44,max_rows_=5,frmt_="(ES9.2)")
  !do i=1,size(strs)
  !   if (i.eq.3) then
  !      print '("A = |",A, "|")', strs(i)
  !   else
  !      print '("    |",A, "|")', strs(i)
  !   endif
  !enddo
  deallocate(A)
  print '("")'

  ! in Matlab:
  ! >> A = [1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12. 13. 14. 15. 16. 17. 18. 19. 20. 21. 22. 23. 24.]
  ! >> reshape(A, 3,8)
  !ans =
  !     1     4     7    10    13    16    19    22
  !     2     5     8    11    14    17    20    23
  !     3     6     9    12    15    18    21    24
  A = reshape([1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, &
               9.0d0, 1.0d1, 1.1d1, 1.2d1, 1.3d1, 1.4d1, 1.5d1, 1.6d1, &
               1.7d1, 1.8d1, 1.9d1, 2.0d1, 2.1d1, 2.2d1, 2.3d1, 2.4d1 ], [3, 8])
  print '("A(3,8) = ")'
  call pprint_matrix (A, frmt_='(F3.0)')

  ! >> reshape(A, 3,4,2)
  !
  ! ans(:,:,1) =
  !      1     4     7    10
  !      2     5     8    11
  !      3     6     9    12
  !
  ! ans(:,:,2) =
  !     13    16    19    22
  !     14    17    20    23
  !     15    18    21    24
  print '("A(3,4,2) = ")'
  call pprint_matrix3d (reshape(A, [3,4,2]), &
                        flip23_=.true., frmt_='(ES14.7)', line_len_=80)
  deallocate(A)

  allocate(A(30, 40*50))
  call random_number(A)
  A = A - 0.5
  print '("A(12, 100, 50) = ")'
  call pprint_matrix3d (reshape(A, [12,100,50]), &
                        flip23_=.true., frmt_='(ES14.7)', line_len_=120)
  deallocate(A)

  A = transpose(reshape([0.8147d0, 0.9134d0, 0.2785d0, &
                         0.9058d0, 0.6324d0, 0.5469d0, &
                         0.1270d0, 0.0975d0, 0.9575d0], [3, 3]))
  print '("A = ")'
  call pprint_matrix (A)

  R2 = matmul(transpose(A), A)   ! R2 = A'*A
  R2_g = transpose(reshape([ &
    1.5004d0, 1.3293d0, 0.8439d0, &
    1.3293d0, 1.2436d0, 0.6936d0, &
    0.8439d0, 0.6936d0, 1.2935d0] &
  , [3,3]))
  print '("R2 error: ", ES12.5)', sum((R2 - R2_g)**2)
  print '("Frobenius norm of the R2 error: ", ES12.5)', norm(R2 - R2_g)

  R2_inv = R2
  print '("testing call matinv_lu_d(R2, R2_inv):")'
  call matinv_lu_d(R2, R2_inv)
  call pprint_matrix (matmul(R2_inv,R2))
  print '("|I - R2_inv*R2|_2 = ",ES14.7)', &
        norm(eye(3) - matmul(R2, R2_inv))
  print '("testing call matinv(R2, R2_inv, alg=`t`):")'
  call matinv(R2, R2_inv, alg='t')
  call pprint_matrix (matmul(R2_inv,R2))
  print '("|I - R2_inv*R2|_2 = ",ES14.7)', &
        norm(eye(3) - matmul(R2, R2_inv))
  print '("testing call matinv(R2, R2_inv, alg=`s`):")'
  call matinv(R2, R2_inv, alg='s')
  call pprint_matrix (matmul(R2_inv,R2))
  print '("|I - R2_inv*R2|_2 = ",ES14.7)', &
        norm(eye(3) - matmul(R2, R2_inv))

  call d_svd(R2,u,s,v) ! [u,s,v]=svd(R2, 'econ');
  v = transpose(v)     ! need to transpose v for compativility with matlab
  u = matmul(A, v)     ! u = A*v
  s = sum(u**2, 1)     ! s = sum(u.^2, 1);
  s = sqrt(s)          ! s = sqrt(s.') "s.'" is a non-conjugate transpose (simple transpose)
  print '("After s = sqrt(...):")'
  print '("          s = [ ",3(F6.4,1X),"]")', s
  print '("Expected: s = [ 1.8168 0.8389 0.1815 ]")'

  call d_svd(A, u, s, v)
  print '("A = ")'
  call pprint_matrix (A)
  print '("u = ")'
  call pprint_matrix (u)
  print '("v = ")'
  call pprint_matrix (v)
  print '("s = ",1000(F7.4,1X))', s(:)



  ug = transpose(reshape([ &
   -0.6612d0,   -0.4121d0,   -0.6269d0, &
   -0.6742d0,   -0.0400d0,    0.7374d0, &
   -0.3290d0,    0.9103d0,   -0.2514d0],&
   [3, 3]))
  vg = reshape([ &
   -0.6556d0,   -0.3056d0,    0.6905d0, &
   -0.5848d0,   -0.3730d0,   -0.7204d0, &
   -0.4777d0,    0.8761d0,   -0.0659d0],&
   [3, 3])
  sg = [ 1.8168d0, 0.8389d0, 0.1815d0]

  print '("|u(1) - ug(1)|_2 = ", ES12.5)', sum((u(:,1) + ug(1,:))**2)
  print '("|u(2) - ug(2)|_2 = ", ES12.5)', sum((u(:,2) + ug(2,:))**2)
  print '("|u(3) - ug(3)|_2 = ", ES12.5)', sum((u(:,3) - ug(3,:))**2)
  print '("|s - sg|_2 = ", ES12.5)', sum((s - sg)**2)
  print '("|v(1) - vg(1)|_2 = ", ES12.5)', sum((v(1,:) + vg(1,:))**2)
  print '("|v(2) - vg(2)|_2 = ", ES12.5)', sum((v(2,:) + vg(2,:))**2)
  print '("|v(3) - vg(3)|_2 = ", ES12.5)', sum((v(3,:) - vg(3,:))**2)
  print '("ug = ")'
  call pprint_matrix (ug)
  print '("vg = ")'
  call pprint_matrix (vg)
  print '("|A - u*s*v_g|_2 = ",ES14.7)', &
        norm(A - matmul(ug, matmul(diag(sg), vg)))
  print '("|A - u*s*v|_2 = ",ES14.7)', &
        norm(A - matmul(u, matmul(diag(s), v)))

  print '("-----")'
  print '("testing matrix inversion")'
  call random_number(A)
  print '("A = ")'
  call pprint_matrix(A)
  A_inv = A
  call matinv(A, A_inv, 't')
  print '("A^{-1} = ")'
  call pprint_matrix(A_inv)
  print '("A^{-1} * A = ")'
  call pprint_matrix(matmul(A_inv, A))
  print '("|I - A*A^{-1}|_2 = ",ES14.7)', &
        norm(eye(3) - matmul(A, A_inv))

  print '("-----")'
  print '("testing QR")'
  A = transpose(reshape( &
      [0.8147d0, 0.6324d0, 0.9575d0, 0.9572d0, &
       0.9058d0, 0.0975d0, 0.9649d0, 0.4854d0, &
       0.1270d0, 0.2785d0, 0.1576d0, 0.8003d0, &
       0.9134d0, 0.5469d0, 0.9706d0, 0.1419d0],&
      [4,4]))
  call qr(C, R, A)
  print '("A = ")'
  call pprint_matrix(A)
  print '("Q = ")'
  call pprint_matrix(C)
  print '("R = ")'
  call pprint_matrix(R)
  print '("|A - Q*R|_2 = ",ES14.7)', &
        norm(A - matmul(C, R))
  print '("Q is orthogonal: |I - Q''*Q|_2 = ",ES14.7)', &
        norm(eye(4) - matmul(transpose(C), C))

  print '("-----")'
  print '("testing QR with rectangular matrix 4x3:")'
  A = A(:,1:3)
  call qr(C, R, A)
  print '("Q = ")'
  call pprint_matrix(C)
  print '("R = ")'
  call pprint_matrix(R)
  print '("|A - Q*R|_2 = ",ES14.7)', &
        norm(A - matmul(C, R))
  print '("Q is orthogonal: |I - Q''*Q|_2 = ",ES14.7)', &
        norm(eye(3) - matmul(transpose(C), C))

  print '("-----")'
  print '("testing QR with rectangular matrix 3x4:")'
  A = reshape(A, [3, 4])
  call qr(C, R, A)
  print '("Q = ")'
  call pprint_matrix(C)
  print '("R = ")'
  call pprint_matrix(R)
  print '("A - Q*R = ")'
  call pprint_matrix(A - matmul(C,R))
  print '("|A - Q*R|_2 = ",ES14.7)', &
        norm(A - matmul(C, R))
  print '("Q is orthogonal: |I - Q''*Q|_2 = ",ES14.7)', &
        norm(eye(3) - matmul(transpose(C), C))

  print '("-----")'
  print '("testing Cholesky factorization:")'
  A = transpose(reshape( &
      [2.3346d0, 1.1384d0, 2.5606d0, 1.4507d0, &
       1.1384d0, 0.7860d0, 1.2743d0, 0.9531d0, &
       2.5606d0, 1.2743d0, 2.8147d0, 1.6487d0, &
       1.4507d0, 0.9531d0, 1.6487d0, 1.8123d0],&
      [4,4]))
  print '("A = ")'
  call pprint_matrix(A)
  Uc = d_chol(A)
  print '("Uc = ")'
  call pprint_matrix(Uc)
  Ucg = transpose(reshape( &
        [1.5279d0,    0.7451d0,  1.6759d0,  0.9494d0, &
              0d0,    0.4805d0,  0.0534d0,  0.5113d0, &
              0d0,         0d0,  0.0580d0,  0.5216d0, &
              0d0,         0d0,       0d0,  0.6143d0],&
        [4,4]))
  print '("Error = ",ES14.7)', norm(Uc - Ucg)/sqrt(10d0)
  print '("|A - U''*U|_2 = ",ES14.7)', norm(A - matmul(transpose(Uc),Uc))

  print '("-----")'
  print '("testing the call to d2_lual with zero input")'
  allocate(g(1)); g(1) = 0d0
  allocate(col(1,2,1)); col = 0d0
  call d2_lual(2,1,g,col)
  err= sum(col**2)
  if (dabs(err) < 1d-15) then
     print '("[+] call to d2_lual with zero input: PASS")'
  else
     print '("[-] call to d2_lual with zero input: FAIL")'
     failed= failed + 1
  end if
  deallocate(g,col)

  print '("-----")'
  print '("testing svdgram")'
  A = transpose(reshape( &
      [0.8147d0, 0.9134d0, 0.2785d0, &
       0.9058d0, 0.6324d0, 0.5469d0, &
       0.1270d0, 0.0975d0, 0.9575d0],&
      [3, 3]))
  print '("|A - u*s*v_g|_2 = ",ES14.7)', &
        norm(A - matmul(ug, matmul(diag(sg), vg)))
  call svdgram(u,s,v,A,1d-15)
  print '("u = ")'
  call pprint_matrix(u)
  print '("ug = ")'
  call pprint_matrix(ug)
  print '("sv = ")'
  call pprint_matrix(v)
  print '("s*vg = ")'
  call pprint_matrix(matmul(diag(s),vg))
  print '("s = ")'
  call pprint_matrix(diag(s))
  print '("sg = ")'
  call pprint_matrix (diag(sg))
  print '("|A - u*sv|_2 = ",ES14.7)', &
        norm(A - matmul(u, v))
  print '("|A - u*s*v_g|_2 = ",ES14.7)', &
        norm(A - matmul(ug, matmul(diag(sg), vg)))

  deallocate(C, R)

  print '("-----")'
  print '("testing dtt_random_ortho")'
  ta= dtt_random_ortho([3,3,3,3], r_=[1,2,2,2,1])
  print '("ta = ")'
  call pprint(ta)
  C = reshape(ta% u(1)% p, [3,2])
  print '("Core 1 of ta is orthogonal: |I - C''*C|_2 = ",ES14.7)', &
        norm(eye(2) - matmul(transpose(C),C))
  C = reshape(ta% u(2)% p, [6,2])
  print '("Core 2 of ta is orthogonal: |I - C''*C|_2 = ",ES14.7)', &
        norm(eye(2) - matmul(transpose(C),C))
  C = reshape(ta% u(3)% p, [6,2])
  print '("Core 3 of ta is orthogonal: |I - C''*C|_2 = ",ES14.7)', &
        norm(eye(2) - matmul(transpose(C),C))
  C = reshape(ta% u(4)% p, [6,1])
  print '("Core 4 of ta is orthogonal: |I - C''*C|_2 = ",ES14.7)', &
        norm(eye(1) - matmul(transpose(C),C))
  print '("testing simple form, with right-to-left ortho:")'
  ta= dtt_random_ortho(3, 4, r_=2, lr_=.false.)
  print '("ta = ")'
  call pprint(ta)
  C = reshape(ta% u(1)% p, [1,6])
  print '("Core 1^T of ta is orthogonal: |I - C*C''|_2 = ",ES14.7)', &
        norm(eye(1) - matmul(C,transpose(C)))
  C = reshape(ta% u(2)% p, [2,6])
  print '("Core 2^T of ta is orthogonal: |I - C*C''|_2 = ",ES14.7)', &
        norm(eye(2) - matmul(C,transpose(C)))
  C = reshape(ta% u(3)% p, [2,6])
  print '("Core 3^T of ta is orthogonal: |I - C*C''|_2 = ",ES14.7)', &
        norm(eye(2) - matmul(C,transpose(C)))
  C = reshape(ta% u(4)% p, [2,3])
  print '("Core 4^T of ta is orthogonal: |I - C*C''|_2 = ",ES14.7)', &
        norm(eye(2) - matmul(C,transpose(C)))
  C = reshape(ta% u(2)% p, [2,6])
  print '("ta%u(2), reshaped to [2,6] = ")'
  call pprint_matrix(C)

  print '("-----")'
  print '("testing index ravel/unravel:")'
  ishp = [3, 5]; i = 13
  mind = unravel_index(i, ishp)
  print '(" - shape: ["2(I2,1X)"]")', ishp
  print '(" - unraveling index ij="I3": i,j="2(I2,1X))', i, mind
  print '(" - raveling multiindex i,j="2(I2,1X)": ij="I3)', mind, ravel_multi_index(mind, ishp)
  ishp = [4, 11, 10, 13, 5]; i = 0
  mind = unravel_index(i, ishp)
  print '(" - shape: ["5(I2,1X)"]")', ishp
  print '(" - unraveling index ijk="I4": i,j..="5(I2,1X))', i, mind
  print '(" - raveling multiindex i,j..="5(I2,1X)": ijk="I4)', mind, ravel_multi_index(mind, ishp)

end program
