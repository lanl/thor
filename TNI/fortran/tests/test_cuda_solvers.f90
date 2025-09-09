program main

use cuda_solvers , only: dgesvd 
  implicit none
  real(8), allocatable :: A_mat(:,:), u(:,:),vt(:,:),s(:)
  integer :: m,n,mn, lda, info
  m = 3; n=2
  mn = min(m,n); lda=m
  allocate(A_mat(m,n),stat=info)
  if(info.ne.0)then;print*,"[!][test_cusolver] Unable to allocate A_mat";return;endif
  A_mat = reshape((/1.0, 4.0, 2.0, 2.0, 5.0, 1.0/), (/m,n/))
  allocate(u(lda,m),  vt(lda,n))
  allocate(s(1:min(m,n)))

  call dgesvd(m,n,mn, A_mat, u, s, vt)
  print*," U = ", u
  print*," S = ", s
  print*," V = ", vt
  deallocate(a_mat,u,s,vt)

end program main
