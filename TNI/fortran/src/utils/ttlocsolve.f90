module ttlocsolve_lib
  use trans_lib
  implicit none

  !!!!!
  !! Standard matrix routines (transposed, slices), and
  !! ALS-related local operations: matvecs, gmres, precs, etc.
  !!!!!

contains


  ! 1324 - permutation. Necessary for tt-matrices, etc.
  subroutine dperm1324(n1,n2,n3,n4, A, B)
    integer, intent(in) :: n1,n2,n3,n4
    double precision, intent (in) :: A(n1,n2,n3,n4)
    double precision, intent (inout) :: B(n1,n3,n2,n4)
    integer i2,i3,i4

    do i4=1,n4
       do i3=1,n3
          do i2=1,n2
             call dcopy(n1, A(1,i2,i3,i4), 1, B(1,i3,i2,i4), 1)
          end do
       end do
    end do
  end subroutine

  subroutine zperm1324(n1,n2,n3,n4, A, B)
    integer, intent(in) :: n1,n2,n3,n4
    double complex, intent (in) :: A(n1,n2,n3,n4)
    double complex, intent (inout) :: B(n1,n3,n2,n4)
    integer i2,i3,i4

    do i4=1,n4
       do i3=1,n3
          do i2=1,n2
             call zcopy(n1, A(1,i2,i3,i4), 1, B(1,i3,i2,i4), 1)
          end do
       end do
    end do
  end subroutine



  !!!!!
  !! Local ALS matvec using given phi1, A, phi2
  !!!!!
  subroutine d2d_mv(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y, work1, work2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    double precision, intent(in) :: phi1(*), A(*), phi2(*), x(*)
    double precision, intent(out) :: y(*)
    double precision, intent(inout),pointer,optional :: work1(:), work2(:)
    character(len=*),parameter :: subnam='d2d_mv'
    double precision, pointer :: res1(:), res2(:)
    integer :: info
    logical :: loc1,loc2

    loc1=.true.; if(present(work1))then;if(associated(work1))loc1=.false.;endif
    loc2=.true.; if(present(work2))then;if(associated(work2))loc2=.false.;endif
    if(loc1)then
     allocate(res1(rx1*max(m*ra2, ra1*n)*ry2),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate res1';stop;endif
    else
     res1=>work1
    endif
    if(loc2)then
     allocate(res2(rx1*max(m*ra2, ra1*n)*ry2),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate res2';stop;endif
    else
     res2=>work2
    endif

    !   phi2: rx2, ra2, ry2
    !  x: rx1,m,rx2
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2
    call trans2d(rx1, m*ra2*ry2, res1, res2)
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res2, m*ra2, 0d0, res1, ra1*n)
    !     res2: ra1,n,ry2,rx1
    call trans2d(ra1*n*ry2, rx1, res1, res2)
    !    phi1: ry1, rx1, ra1
    call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res2, rx1*ra1, 0d0, y, ry1)

    if(loc1)deallocate(res1)
    if(loc2)deallocate(res2)
  end subroutine

  subroutine z2d_mv(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y, work1, work2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    double complex, intent(in) :: phi1(*), A(*), phi2(*), x(*)
    double complex, intent(inout) :: y(*)
    double complex, intent(inout),pointer,optional :: work1(:), work2(:)
    character(len=*),parameter :: subnam='d2d_mv'
    double complex, pointer :: res1(:), res2(:)
    integer :: info
    logical :: loc1,loc2

    loc1=.true.; if(present(work1))then;if(associated(work1))loc1=.false.;endif
    loc2=.true.; if(present(work2))then;if(associated(work2))loc2=.false.;endif
    if(loc1)then
     allocate(res1(rx1*max(m*ra2, ra1*n)*ry2),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate res1';stop;endif
    else
     res1=>work1
    endif
    if(loc2)then
     allocate(res2(rx1*max(m*ra2, ra1*n)*ry2),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate res2';stop;endif
    else
     res2=>work2
    endif

    !   phi2: rx2, ra2, ry2
    !  x: rx1,m,rx2
    call zgemm('N', 'N', rx1*m, ra2*ry2, rx2, (1d0,0d0), x, rx1*m, phi2, rx2, (0d0,0d0), res1, rx1*m)
    !    res1: rx1,m,ra2,ry2
    call trans2z(rx1, m*ra2*ry2, res1, res2)
    call zgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, (1d0,0d0), A, ra1*n, res2, m*ra2, (0d0,0d0), res1, ra1*n)
    !     res2: ra1,n,ry2,rx1
    call trans2z(ra1*n*ry2, rx1, res1, res2)
    !    phi1: ry1, rx1, ra1
    call zgemm('N', 'N', ry1, n*ry2, rx1*ra1, (1d0,0d0), phi1, ry1, res2, rx1*ra1, (0d0,0d0), y, ry1)

    if(loc1)deallocate(res1)
    if(loc2)deallocate(res2)
  end subroutine


!!! Generate full local matrix
  subroutine d2d_fullmat(rx1, m, rx2, ra1, ra2, phi1, A, phi2, B)
    integer, intent(in) :: rx1, m, rx2, ra1, ra2
    double precision, intent(in) :: phi1(*), A(*), phi2(*)
    double precision, intent(inout) :: B(*)
    double precision, allocatable :: res1(:), res2(:), phi2t(:)

    allocate(res1(rx1*m*rx1*m*max(ra2, rx2*rx2)), res2(rx1*m*rx1*m*max(ra2, rx2*rx2)), phi2t(ra2*rx2*rx2))

    ! phi1: ry1,rx1,ra1
    call dgemm('N', 'N', rx1*rx1, m*m*ra2, ra1, 1d0, phi1, rx1*rx1, A, ra1, 0d0, res1, rx1*rx1)
    ! res1: ry1,rx1,n,m,ra2
    call dperm1324(rx1, rx1, m, m*ra2, res1, res2)
    ! phi2: rx2,ra2,ry2
    call trans2d(rx2, ra2*rx2, phi2, phi2t)
    call dgemm('N', 'N', rx1*m*rx1*m, rx2*rx2, ra2, 1d0, res2, rx1*m*rx1*m, phi2t, ra2, 0d0, res1, rx1*m*rx1*m);
    call dperm1324(rx1*m, rx1*m, rx2, rx2, res1, B)
    deallocate(res1, res2, phi2t)
  end subroutine

  subroutine z2d_fullmat(rx1, m, rx2, ra1, ra2, phi1, A, phi2, B)
    integer, intent(in) :: rx1, m, rx2, ra1, ra2
    double complex, intent(in) :: phi1(*), A(*), phi2(*)
    double complex, intent(inout) :: B(*)
    double complex, allocatable :: res1(:), res2(:), phi2t(:)

    allocate(res1(rx1*m*rx1*m*max(ra2, rx2*rx2)), res2(rx1*m*rx1*m*max(ra2, rx2*rx2)), phi2t(ra2*rx2*rx2))

    ! phi1: ry1,rx1,ra1
    call zgemm('N', 'N', rx1*rx1, m*m*ra2, ra1, (1d0,0d0), phi1, rx1*rx1, A, ra1, (0d0,0d0), res1, rx1*rx1)
    ! res1: ry1,rx1,n,m,ra2
    call zperm1324(rx1, rx1, m, m*ra2, res1, res2)
    ! phi2: rx2,ra2,ry2
    call trans2z(rx2, ra2*rx2, phi2, phi2t)
    call zgemm('N', 'N', rx1*m*rx1*m, rx2*rx2, ra2, (1d0,0d0), res2, rx1*m*rx1*m, phi2t, ra2, (0d0,0d0), res1, rx1*m*rx1*m);
    call zperm1324(rx1*m, rx1*m, rx2, rx2, res1, B)
    deallocate(res1, res2, phi2t)
  end subroutine


  subroutine d2d_jac_gen(ptype, rx1, n, rx2, ra1, ra2, Phi1, A, Phi2, jacs)
    ! sizes of work1, work2:
    !   center:
    !       rx*ra, rx1*n*n*rx2, rx2*ra2*rx2
    !   left:
    !       rx*ra, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2, rx2*ra2*rx2
    !   right:
    !       rx*ra, rx1*n*n*rx2*rx2, rx2*ra2*rx2
    ! sizes of jacs:
    !   center:
    !       rx1*n*n*ra2,  rx1*n*n*rx2
    !   left:
    !       rx1*rx1*n*n*ra2,  rx1*rx1*n*n*rx2
    !   right:
    !       n*n*rx2*rx2*rx1, rx1*n*n*ra2
    character, intent(in):: ptype
    integer, intent(in) :: rx1, n, rx2, ra1, ra2
    double precision, intent(in) :: Phi1(*), A(*), Phi2(*)
    double precision, intent(inout),pointer :: jacs(:) !, work1(*), work2(*)
    double precision, allocatable :: work1(:), work2(:)
    integer, allocatable :: ipiv(:)
    integer :: i, info


    if ((ptype=='c').or.(ptype=='C')) then

       i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2)
       allocate(work1(i))
       allocate(work2(i))
       allocate(ipiv(n))

       call trans2d(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do
       do i=0,ra1-1
          call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do

       call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
       call dgemm('N', 'T', rx1*n*n, rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2, 0d0, work1, rx1*n*n)
       call trans2d(rx1, n*n*rx2, work1, work2)
       ! inversion
       do i=1,n*n
          work1(i)=0d0
       end do
       do i=1,n
          work1(i+(i-1)*n)=1d0
       end do
       do i=0,rx2*rx1-1
          call dcopy(n*n, work1, 1, jacs(i*n*n+1:(i+1)*n*n), 1)
          call dgesv(n, n, work2(i*n*n+1), n, ipiv, jacs(i*n*n+1), n, info)
       end do
    end if

    if ((ptype=='l').or.(ptype=='L')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx1*n
       allocate(ipiv(i))

       call trans2d(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do

       call dgemm('N', 'N', rx1*rx1, n*n*ra2, ra1, 1d0, Phi1, rx1*rx1, A, ra1, 0d0, work1, rx1*rx1)
       call dperm1324(rx1, rx1, n, n*ra2, work1, jacs)
       call dgemm('N', 'T', rx1*n*rx1*n, rx2, ra2, 1d0, jacs, rx1*n*rx1*n, work2, rx2, 0d0, work1, rx1*n*rx1*n)
       ! inversion
       do i=1,rx1*n*rx1*n
          work2(i)=0d0
       end do
       do i=1,rx1*n
          work2(i+(i-1)*rx1*n)=1d0
       end do
       do i=0,rx2-1
          call dcopy(rx1*n*rx1*n, work2, 1, jacs(i*rx1*n*rx1*n+1:(i+1)*rx1*n*rx1*n), 1)
          call dgesv(rx1*n, rx1*n, work1(i*rx1*n*rx1*n+1), rx1*n, ipiv, jacs(i*rx1*n*rx1*n+1), rx1*n, info)
          if(info.ne.0)then;write(*,*)': dgesv problem ', info;endif
       end do
    end if

    if ((ptype=='r').or.(ptype=='R')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx2*n
       allocate(ipiv(i))

       call trans2d(rx2*ra2, rx2, Phi2, work2)
       do i=0,ra1-1
          call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do

       call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
       call dgemm('N', 'T', rx1*n*n, rx2*rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2*rx2, 0d0, work1, rx1*n*n)
       call dperm1324(rx1*n, n, rx2, rx2, work1, work2)
       call trans2d(rx1, n*rx2*n*rx2, work2, work1)
       ! inversion
       do i=1,rx2*n*rx2*n
          work2(i)=0d0
       end do
       do i=1,rx2*n
          work2(i+(i-1)*rx2*n)=1d0
       end do
       do i=0,rx1-1
          call dcopy(rx2*n*rx2*n, work2, 1, jacs(i*rx2*n*rx2*n+1:(i+1)*rx2*n*rx2*n), 1)
          call dgesv(rx2*n, rx2*n, work1(i*rx2*n*rx2*n+1), rx2*n, ipiv, jacs(i*rx2*n*rx2*n+1), rx2*n, info)
       end do
    end if

    deallocate(work1,work2, ipiv)
  end subroutine

  subroutine z2d_jac_gen(ptype, rx1, n, rx2, ra1, ra2, Phi1, A, Phi2, jacs)
    ! sizes of work1, work2:
    !   center:
    !       rx*ra, rx1*n*n*rx2, rx2*ra2*rx2
    !   left:
    !       rx*ra, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2, rx2*ra2*rx2
    !   right:
    !       rx*ra, rx1*n*n*rx2*rx2, rx2*ra2*rx2
    ! sizes of jacs:
    !   center:
    !       rx1*n*n*ra2,  rx1*n*n*rx2
    !   left:
    !       rx1*rx1*n*n*ra2,  rx1*rx1*n*n*rx2
    !   right:
    !       n*n*rx2*rx2*rx1, rx1*n*n*ra2
    character, intent(in):: ptype
    integer, intent(in) :: rx1, n, rx2, ra1, ra2
    double complex, intent(in) :: Phi1(*), A(*), Phi2(*)
    double complex, intent(inout),pointer :: jacs(:) !, work1(*), work2(*)
    double complex, allocatable :: work1(:), work2(:)
    integer, allocatable :: ipiv(:)
    integer :: i, info


    if ((ptype=='c').or.(ptype=='C')) then

       i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2)
       allocate(work1(i))
       allocate(work2(i))
       allocate(ipiv(n))

       call trans2z(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          call zcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do
       do i=0,ra1-1
          call zcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do

       call zgemm('N', 'N', rx1, n*n*ra2, ra1, (1d0,0d0), work1, rx1, A, ra1, (0d0,0d0), jacs, rx1)
       call zgemm('N', 'T', rx1*n*n, rx2, ra2, (1d0,0d0), jacs, rx1*n*n, work2, rx2, (0d0,0d0), work1, rx1*n*n)
       call trans2z(rx1, n*n*rx2, work1, work2)
       ! inversion
       do i=1,n*n
          work1(i)=(0d0,0d0)
       end do
       do i=1,n
          work1(i+(i-1)*n)=(1d0,0d0)
       end do
       do i=0,rx2*rx1-1
          call zcopy(n*n, work1, 1, jacs(i*n*n+1), 1)
          call zgesv(n, n, work2(i*n*n+1), n, ipiv, jacs(i*n*n+1), n, info)
       end do
    end if

    if ((ptype=='l').or.(ptype=='L')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx1*n
       allocate(ipiv(i))

       call trans2z(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          call zcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do

       call zgemm('N', 'N', rx1*rx1, n*n*ra2, ra1, (1d0,0d0), Phi1, rx1*rx1, A, ra1, (0d0,0d0), work1, rx1*rx1)
       call zperm1324(rx1, rx1, n, n*ra2, work1, jacs)
       call zgemm('N', 'T', rx1*n*rx1*n, rx2, ra2, (1d0,0d0), jacs, rx1*n*rx1*n, work2, rx2, (0d0,0d0), work1, rx1*n*rx1*n)
       ! inversion
       do i=1,rx1*n*rx1*n
          work2(i)=(0d0,0d0)
       end do
       do i=1,rx1*n
          work2(i+(i-1)*rx1*n)=(1d0,0d0)
       end do
       do i=0,rx2-1
          call zcopy(rx1*n*rx1*n, work2, 1, jacs(i*rx1*n*rx1*n+1), 1)
          call zgesv(rx1*n, rx1*n, work1(i*rx1*n*rx1*n+1), rx1*n, ipiv, jacs(i*rx1*n*rx1*n+1), rx1*n, info)
          if(info.ne.0)then;write(*,*)': dgesv problem ', info;endif
       end do
    end if

    if ((ptype=='r').or.(ptype=='R')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx2*n
       allocate(ipiv(i))

       call trans2z(rx2*ra2, rx2, Phi2, work2)
       do i=0,ra1-1
          call zcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do

       call zgemm('N', 'N', rx1, n*n*ra2, ra1, (1d0,0d0), work1, rx1, A, ra1, (0d0,0d0), jacs, rx1)
       call zgemm('N', 'T', rx1*n*n, rx2*rx2, ra2, (1d0,0d0), jacs, rx1*n*n, work2, rx2*rx2, (0d0,0d0), work1, rx1*n*n)
       call zperm1324(rx1*n, n, rx2, rx2, work1, work2)
       call trans2z(rx1, n*rx2*n*rx2, work2, work1)
       ! inversion
       do i=1,rx2*n*rx2*n
          work2(i)=(0d0,0d0)
       end do
       do i=1,rx2*n
          work2(i+(i-1)*rx2*n)=(1d0,0d0)
       end do
       do i=0,rx1-1
          call zcopy(rx2*n*rx2*n, work2, 1, jacs(i*rx2*n*rx2*n+1), 1)
          call zgesv(rx2*n, rx2*n, work1(i*rx2*n*rx2*n+1), rx2*n, ipiv, jacs(i*rx2*n*rx2*n+1), rx2*n, info)
       end do
    end if

    deallocate(work1,work2, ipiv)
  end subroutine


  subroutine d2d_jac_apply(ptype, rx1, n, rx2, jacs, x, y, res1)
    ! sizes of work1: rx1*n*rx2
    character, intent(in) :: ptype
    integer, intent(in) :: rx1, n, rx2
    double precision, intent(in) :: jacs(*), x(*)
    double precision, intent(inout) :: y(*)
    double precision, intent(inout),pointer,optional :: res1(:)
    double precision, pointer :: work1(:)
    integer i, info
    logical :: loc1

    loc1=.true.; if(present(res1))then;if(associated(res1))loc1=.false.;endif
    if(loc1)then
     allocate(work1(rx1*n*rx2),stat=info)
     if(info.ne.0)then;write(*,*)': cannot allocate work1';stop;endif
    else
     work1=>res1
    endif

    if ((ptype=='c').or.(ptype=='C')) then
       ! jacs is n,n,rx2,rx1
       call trans2d(rx1, n*rx2, x, work1)
       do i=0,(rx2*rx1-1)
          call dgemv('N', n, n, 1d0, jacs(i*n*n+1), n, work1(i*n+1), 1, 0d0, y(i*n+1), 1)
       end do
       call trans2d(n*rx2, rx1, y, work1)
       call dcopy(rx1*n*rx2, work1, 1, y, 1)
    end if

    if ((ptype=='l').or.(ptype=='L')) then
       ! jacs is rx1*n,rx1*n, rx2
       call dcopy(rx1*n*rx2, x, 1, work1, 1)
       do i=0,rx2-1
          call dgemv('N', rx1*n, rx1*n, 1d0, jacs(i*rx1*n*rx1*n+1), rx1*n, work1(i*rx1*n+1), 1, 0d0, y(i*rx1*n+1), 1)
       end do
    end if

    if ((ptype=='r').or.(ptype=='R')) then
       ! jacs is n*rx2, n*rx2, rx1
       call trans2d(rx1, n*rx2, x, work1)
       do i=0,rx1-1
          call dgemv('N', n*rx2, n*rx2, 1d0, jacs(i*n*rx2*n*rx2+1), n*rx2, work1(i*n*rx2+1), 1, 0d0, y(i*n*rx2+1), 1)
       end do
       call trans2d(n*rx2, rx1, y, work1)
       call dcopy(rx1*n*rx2, work1, 1, y, 1)
    end if

    if(loc1)deallocate(work1)
  end subroutine

  subroutine z2d_jac_apply(ptype, rx1, n, rx2, jacs, x, y, res1)
    ! sizes of work1: rx1*n*rx2
    character, intent(in) :: ptype
    integer, intent(in) :: rx1, n, rx2
    double complex, intent(in) :: jacs(*), x(*)
    double complex, intent(inout) :: y(*)
    double complex, intent(inout),pointer,optional :: res1(:)
    double complex, pointer :: work1(:)
    integer i, info
    logical :: loc1

    loc1=.true.; if(present(res1))then;if(associated(res1))loc1=.false.;endif
    if(loc1)then
     allocate(work1(rx1*n*rx2),stat=info)
     if(info.ne.0)then;write(*,*)': cannot allocate work1';stop;endif
    else
     work1=>res1
    endif


    if ((ptype=='c').or.(ptype=='C')) then
       ! jacs is n,n,rx2,rx1
       call trans2z(rx1, n*rx2, x, work1)
       do i=0,(rx2*rx1-1)
          call zgemv('N', n, n, (1d0,0d0), jacs(i*n*n+1), n, work1(i*n+1), 1, (0d0,0d0), y(i*n+1), 1)
       end do
       call trans2z(n*rx2, rx1, y, work1)
       call zcopy(rx1*n*rx2, work1, 1, y, 1)
    end if

    if ((ptype=='l').or.(ptype=='L')) then
       ! jacs is rx1*n,rx1*n, rx2
       call zcopy(rx1*n*rx2, x, 1, work1, 1)
       do i=0,rx2-1
          call zgemv('N', rx1*n, rx1*n, (1d0,0d0), jacs(i*rx1*n*rx1*n+1), rx1*n, work1(i*rx1*n+1), 1, (0d0,0d0), y(i*rx1*n+1), 1)
       end do
    end if

    if ((ptype=='r').or.(ptype=='R')) then
       ! jacs is n*rx2, n*rx2, rx1
       call trans2z(rx1, n*rx2, x, work1)
       do i=0,rx1-1
          call zgemv('N', n*rx2, n*rx2, (1d0,0d0), jacs(i*n*rx2*n*rx2+1), n*rx2, work1(i*n*rx2+1), 1, (0d0,0d0), y(i*n*rx2+1), 1)
       end do
       call trans2z(n*rx2, rx1, y, work1)
       call zcopy(rx1*n*rx2, work1, 1, y, 1)
    end if

    if(loc1)deallocate(work1)
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!
!!!!! GMRES !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
  subroutine d2d_gmresr(Phi1, A, Phi2, rhs, rx1, n, rx2, ra1, ra2, nrestart, tol, niters, ptype, jacs, sol, verb)
    ! right preconditioned - for residual tolerance
    ! This one with Householder tranforms
    double precision,intent(in):: Phi1(*), A(*), Phi2(*), rhs(*)
    double precision,intent(in), pointer :: jacs(:)
    integer, intent(in) :: rx1,n,rx2, ra1,ra2, nrestart, niters, verb
    double precision, intent(in) :: tol
    character, intent(in) :: ptype
    double precision, intent(inout) :: sol(*)
    character(len=*),parameter :: subnam='d2d_gmresr'

    integer :: i,j,it, sz, nrestart0,info
    double precision, pointer :: U(:,:), w(:), R(:,:), JJ(:,:), tau(:), res1(:), res2(:)
    double precision  nrmr, curres, nrmrhs, dbeta, dalpha

    double precision dnrm2
    double precision ddot


    sz = rx1*n*max(ra2,ra1)*rx2;
    allocate(res1(sz), res2(sz),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot alloc res12';stop;endif

    sz = rx1*n*rx2;

    nrestart0 = min(sz, nrestart)

    allocate(U(sz, nrestart0+1), w(sz), tau(nrestart0+1), JJ(2,nrestart0+1), R(nrestart0+1, nrestart0+1))

    do j=1,sz
       sol(j)=0d0
    end do

    do it=0,niters-1
       ! r0
       if (.not.((ptype=='n').or.(ptype=='N'))) then
          call d2d_jac_apply(ptype, rx1, n, rx2, jacs, sol, w, res1=res1)
          call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2)
       end if
       if ((ptype=='n').or.(ptype=='N')) then
          call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, sol, w, work1=res1,work2=res2)
       end if
       call daxpy(sz,-1d0,rhs,1,w,1)

       nrmr = dnrm2(sz,w,1);
       if (it==0) then
          nrmrhs = nrmr
       end if
       if (nrmr==0d0) then
          exit
       end if
       ! initial HHT
       dbeta = nrmr;
       if (w(1)<0d0) then
          dbeta = -dbeta
       end if
       w(1) = w(1)+dbeta
       tau(1) = -dbeta
       nrmr = dnrm2(sz,w,1)
       dbeta = 1d0/nrmr
       call dscal(sz,dbeta,w,1)
       call dcopy(sz,w,1,U,1)

       do j=1,nrestart0
          ! HHT on last U
          call dcopy(sz,U(1:sz,j),1,w,1)
          dbeta = -2d0*U(j,j)
          call dscal(sz,dbeta,w,1)
          w(j) = w(j) + 1d0
          do i=j-1,1,-1
             dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
             call daxpy(sz,dbeta,U(1,i),1,w,1);
          end do
          dbeta = dnrm2(sz,w,1);
          dbeta = 1d0/dbeta;
          call dscal(sz,dbeta,w,1);

          ! precvec, matvec
          if (.not.((ptype=='n').or.(ptype=='N'))) then
             call d2d_jac_apply(ptype, rx1, n, rx2, jacs, w, w, res1=res1);
             call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2);
          end if
          if ((ptype=='n').or.(ptype=='N')) then
             call d2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2);
          end if

          ! Orthog w to j projectors
          do i=1,j
             dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
             call daxpy(sz,dbeta,U(1,i),1,w,1);
          end do

          ! new P_{j+1}
          if (j<sz) then
             do i=1,j
                U(i,j+1)=0d0;
             end do
             i = sz-j;
             call dcopy(i, w(j+1:j+i), 1, U(j+1:j+i, j+1), 1);
             dalpha = dnrm2(i, U(j+1:j+i, j+1), 1);
             if (.not.(dalpha==0d0)) then
                if (w(j+1)<0d0) then
                   dalpha = -dalpha;
                end if
                U(j+1, j+1) = U(j+1,j+1) + dalpha;
                dbeta = dnrm2(i, U(j+1:j+i, j+1), 1);
                dbeta = 1d0/dbeta;
                call dscal(i,dbeta,U(j+1, j+1),1);

                w(j+1) = -dalpha;
                do i=j+2,sz
                   w(i)=0d0;
                end do
             end if
          end if

          ! Givens rotators to the top of w
          do i=1,j-1
             dbeta = w(i);
             w(i) = JJ(1,i)*w(i) + JJ(2,i)*w(i+1);
             w(i+1) = -JJ(2,i)*dbeta + JJ(1,i)*w(i+1);
          end do

          ! New rotator
          if (j<sz) then
             dalpha = dnrm2(2, w(j:j+1), 1)
             JJ(1,j) = w(j)/dalpha;
             JJ(2,j) = w(j+1)/dalpha;
             tau(j+1) = -JJ(2,j)*tau(j);
             tau(j) = JJ(1,j)*tau(j);
             w(j) = dalpha;
             w(j+1) = 0d0;
          end if

          call dcopy(j, w, 1, R(1:j,j), 1);

          ! residual
          curres = dabs(tau(j+1))/nrmrhs;

          if (curres<tol) then
             exit
          end if
       end do ! inner

       if (j>nrestart0) then
          j = nrestart0
       end if

       call dtrsv('u','n','n',j,R,nrestart0+1,tau,1);

       ! Correction
       call dcopy(sz, U(1:sz,j), 1, w, 1);
       dbeta = -2d0*U(j,j)*tau(j);
       call dscal(sz, dbeta, w, 1);
       w(j) = w(j) + tau(j);
       do i=j-1,1,-1
          w(i) = w(i)+tau(i);
          dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
          call daxpy(sz,dbeta,U(1,i),1,w,1);
       end do
       dalpha=-1d0;
       call daxpy(sz,dalpha,w,1,sol,1);
       if (curres<tol) then
        exit;
       end if
    end do ! iters

    if (verb>0) then
       write(*,"(A,I0,A,I0,A,ES10.3)") 'gmres conducted[', it, ',', j, '] iters to relres ', curres
    end if

    deallocate(U,w,tau,JJ,R,res1,res2)
  end subroutine

  ! Matlab's MKL has broken zdotc. I don't know why
  subroutine zdot(n, a, b, res)
    integer,intent(in) :: n
    double complex,intent(in) :: a(*), b(*)
    double complex,intent(inout) :: res
    integer i
    res=(0d0,0d0)
    do i=1,n; res = res+dconjg(a(i))*b(i); end do
  end subroutine
  subroutine z2d_gmresr(Phi1, A, Phi2, rhs, rx1, n, rx2, ra1, ra2, nrestart, tol, niters, ptype, jacs, sol, verb)
    ! right preconditioned - for residual tolerance
    ! This one with Householder tranforms
    double complex,intent(in):: Phi1(*), A(*), Phi2(*), rhs(*)
    double complex,intent(in), pointer :: jacs(:)
    integer, intent(in) :: rx1,n,rx2, ra1,ra2, nrestart, niters, verb
    double precision, intent(in) :: tol
    character, intent(in) :: ptype
    double complex, intent(inout) :: sol(*)
    character(len=*),parameter :: subnam='d2d_gmresr'

    integer :: i,j,it, sz, nrestart0,info
    double complex, allocatable :: U(:,:), w(:), R(:,:), JJ(:,:), tau(:)
    double complex, pointer :: res1(:), res2(:)
    double precision  nrmr, curres, nrmrhs
    double complex  zbeta, zalpha

    double precision,external :: dznrm2

    sz = rx1*n*max(ra2,ra1)*rx2;
    allocate(res1(sz), res2(sz),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot alloc res12';stop;endif

    sz = rx1*n*rx2;

    nrestart0 = min(sz, nrestart)

    allocate(U(sz, nrestart0+1), w(sz), tau(nrestart0+1), JJ(2,nrestart0+1), R(nrestart0+1, nrestart0+1))

    do j=1,sz
       sol(j)=(0d0,0d0)
    end do

    do it=0,niters-1
       ! r0
       if (.not.((ptype=='n').or.(ptype=='N'))) then
          call z2d_jac_apply(ptype, rx1, n, rx2, jacs, sol, w, res1=res1)
          call z2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2)
       end if
       if ((ptype=='n').or.(ptype=='N')) then
          call z2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, sol, w, work1=res1,work2=res2)
       end if
       call zaxpy(sz,(-1d0,0d0),rhs,1,w,1)

       nrmr = dznrm2(sz,w,1);
       if (it==0) then
          nrmrhs = nrmr
       end if
       if (nrmr==0d0) then
          exit
       end if
       ! initial HHT
       zbeta = dcmplx(nrmr,0d0)
       nrmr = abs(w(1))
       if (.not.(nrmr==0d0))zbeta = w(1)*zbeta/dcmplx(nrmr,0d0)
       w(1) = w(1)+zbeta
       tau(1) = -zbeta
       nrmr = dznrm2(sz,w,1)
       zbeta = (1d0,0d0)/dcmplx(nrmr,0d0)
       call zscal(sz,zbeta,w,1)
       call zcopy(sz,w,1,U(1,1),1)

       do j=1,nrestart0
          ! HHT on last U
          call zcopy(sz,U(1,j),1,w,1)
          zbeta = (-2d0,0d0)*dconjg(U(j,j))
          call zscal(sz,zbeta,w,1)
          w(j) = w(j) + (1d0,0d0)
          do i=j-1,1,-1
             call zdot(sz, U(1,i), w, zbeta); zbeta = zbeta*(-2d0,0d0);
             call zaxpy(sz,zbeta,U(1,i),1,w,1);
          end do
          nrmr = dznrm2(sz,w,1);
          zbeta = (1d0,0d0)/dcmplx(nrmr,0d0);
          call zscal(sz,zbeta,w,1);

          ! precvec, matvec
          if (.not.((ptype=='n').or.(ptype=='N'))) then
             call z2d_jac_apply(ptype, rx1, n, rx2, jacs, w, w, res1=res1);
             call z2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2);
          end if
          if ((ptype=='n').or.(ptype=='N')) then
             call z2d_mv(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, work1=res1,work2=res2);
          end if

          ! Orthog w to j projectors
          do i=1,j
             call zdot(sz, U(1,i), w, zbeta); zbeta = zbeta*(-2d0,0d0);
             call zaxpy(sz,zbeta,U(1,i),1,w,1);
          end do

          ! new P_{j+1}
          if (j<sz) then
             do i=1,j
                U(i,j+1)=(0d0,0d0);
             end do
             i = sz-j;
             call zcopy(i, w(j+1), 1, U(j+1, j+1), 1);
             nrmr = dznrm2(i, U(j+1, j+1), 1);
             if (.not.(nrmr==0d0)) then
                zalpha = dcmplx(nrmr,0d0)
                nrmr = abs(w(j+1))
                if (.not.(nrmr==0d0))zalpha = w(j+1)*zalpha/dcmplx(nrmr,0d0)
                U(j+1, j+1) = U(j+1,j+1) + zalpha;
                nrmr = dznrm2(i, U(j+1, j+1), 1);
                zbeta = (1d0,0d0)/dcmplx(nrmr,0d0);
                call zscal(i,zbeta,U(j+1, j+1),1);

                w(j+1) = -zalpha;
                do i=j+2,sz
                   w(i)=(0d0,0d0);
                end do
             end if
          end if

          ! Givens rotators to the top of w
          do i=1,j-1
             zbeta = w(i);
             w(i) = dconjg(JJ(1,i))*w(i) + dconjg(JJ(2,i))*w(i+1);
             w(i+1) = -JJ(2,i)*zbeta + JJ(1,i)*w(i+1);
          end do

          ! New rotator
          if (j<sz) then
             nrmr = dznrm2(2, w(j), 1)
             zalpha = dcmplx(nrmr,0d0)
             JJ(1,j) = w(j)/zalpha;
             JJ(2,j) = w(j+1)/zalpha;
             tau(j+1) = -JJ(2,j)*tau(j);
             tau(j) = dconjg(JJ(1,j))*tau(j);
             w(j) = zalpha;
             w(j+1) = (0d0, 0d0);
          end if

          call zcopy(j, w, 1, R(1,j), 1);

          ! residual
          curres = abs(tau(j+1))/nrmrhs;

          if (curres<tol) then
             exit
          end if
       end do ! inner

       if (j>nrestart0) then
          j = nrestart0
       end if

       call ztrsv('u','n','n',j,R,nrestart0+1,tau,1);

       ! Correction
       call zcopy(sz, U(1,j), 1, w, 1);
       zbeta = (-2d0,0d0)*dconjg(U(j,j))*tau(j);
       call zscal(sz, zbeta, w, 1);
       w(j) = w(j) + tau(j);
       do i=j-1,1,-1
          w(i) = w(i)+tau(i);
!           zbeta = (-2d0,0d0)*zdotc(sz,U(1,i),1,w,1);
          call zdot(sz,U(1,i), w, zbeta); zbeta = zbeta*(-2d0,0d0);
          call zaxpy(sz,zbeta,U(1,i),1,w,1);
       end do
       zalpha=(-1d0,0d0);
       call zaxpy(sz,zalpha,w,1,sol,1);
       if (curres<tol) then
        exit;
       end if
    end do ! iters

    if (verb>0) then
       write(*,"(A,I0,A,I0,A,ES10.3)") 'gmres conducted[', it, ',', j, '] iters to relres ', curres
    end if

    deallocate(U,w,tau,JJ,R,res1,res2)
  end subroutine


!!!!!!!!!!!!!!
!!!! Unfinished Choletskiy !!!!!!!!
!!!!!!!!!!!!!
  subroutine duchol(trans, m, n, dA, numvec, du, lam, real_numvec)
    character, intent(in) :: trans
    integer, intent(in) :: m,n,numvec
    double precision, intent(in) :: dA(m,n)
    double precision, intent(inout) :: du(n,numvec), lam(*)
    integer, intent(inout) :: real_numvec

    double precision,allocatable :: ddiag(:), du2(:,:), work(:), dx(:), dv(:,:)
    integer i,j,cnt, maxi, lwork
    double precision maxa, maxx

    ! functions
    double precision ddot
    integer idamax

    lwork = 3*numvec
    allocate(dx(numvec), dv(numvec, numvec), work(lwork), du2(n,numvec), ddiag(n))

    do j=1,numvec
       lam(j)=0d0
    end do

    !/* Initial diagonal */
    if ((trans=='n').or.(trans=='N')) then
       do i=1,n
          ddiag(i) = ddot(m, dA(1,i), 1, dA(1,i), 1)
       end do
    else
       do i=1,n
          ddiag(i) = ddot(m, dA(i,1), n, dA(i,1), n)
       end do
    end if

    !/* Main cycle */
    do cnt=1,numvec
       !/* maxel */
       maxi = idamax(n, ddiag, 1)
       maxa = ddiag(maxi)

       !/* residual */
       if ((trans=='n').or.(trans=='N')) then
      call dgemv('T', m,n, 1d0, dA, m, dA(1,maxi), 1, 0d0, du(1,cnt), 1)
       else
      call dgemv('N', n,m, 1d0, dA, n, dA(maxi,1), n, 0d0, du(1,cnt), 1)
       end if

       do j=1,cnt-1
          call daxpy(n, -lam(j)*du(maxi,j), du(1,j), 1, du(1,cnt),1)
       end do
       maxa = dsqrt(dabs(du(maxi,cnt)));
       if (maxa<1d-300) then ! // We've got the exact zero matrix
          exit;
       end if
       call dscal(n, 1d0/maxa, du(1,cnt),1)
       do i=1,n
          ddiag(i) = ddiag(i)-(du(i,cnt)*du(i,cnt));
       end do
       !/* reorth */
       call dgeqrf(n, cnt, du, n, du2, work, lwork, i);
       do j=1,cnt
      dx(j) = du(j,cnt)
       end do
       call dorgqr(n, cnt, cnt, du, n, du2, work, lwork, i);

       !/* new eig */
       do i=1,cnt
      do j=1,cnt
             dv(j,i)=dx(j)*dx(i)
      end do
      dv(i,i)=dv(i,i) + lam(i)
       end do
       j=cnt+1;
       call dsyev('V', 'U', cnt, dv, numvec, lam, work, lwork, i);

       !/* update u */
       call dgemm('N','N',n,cnt,cnt,1d0,du,n,dv,numvec,0d0,du2,n);
       call dcopy(n*cnt, du2, 1, du, 1)
    end do ! cnt

    if (cnt>numvec) then
       cnt=numvec
    end if

    deallocate(ddiag, dx, dv, work, du2)
    real_numvec = cnt
  end subroutine

end module
