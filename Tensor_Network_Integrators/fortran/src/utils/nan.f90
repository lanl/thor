module nan_lib
 implicit none
 interface nan
  module procedure d0_nan,d1_nan,d2_nan,d3_nan, z0_nan,z1_nan,z2_nan,z3_nan
 end interface

contains
 logical pure function d0_nan(x) result (isnan)
  double precision,intent(in) :: x
  isnan=.not.((x.lt.0.d0).or.(x.ge.0.d0))
 end function
 logical pure function d1_nan(x) result (isnan)
  double precision,intent(in) :: x(:)
  integer :: n,i
  n=size(x)
  do i=1,n
   isnan=.not.((x(i).lt.0.d0).or.(x(i).ge.0.d0))
   if(isnan)return
  end do
 end function
 logical pure function d2_nan(x) result (isnan)
  double precision,intent(in) :: x(:,:)
  integer :: m,n,i,j
  m=size(x,1);n=size(x,2)
  do j=1,n
   do i=1,m
    isnan=.not.((x(i,j).lt.0.d0).or.(x(i,j).ge.0.d0))
    if(isnan)return
   end do 
  end do
 end function
 logical pure function d3_nan(x) result (isnan)
  double precision,intent(in) :: x(:,:,:)
  integer :: m,n,p,i,j,k
  m=size(x,1);n=size(x,2);p=size(x,3)
  do k=1,p
   do j=1,n
    do i=1,m
     isnan=.not.((x(i,j,k).lt.0.d0).or.(x(i,j,k).ge.0.d0))
     if(isnan)return
    end do 
   end do 
  end do
 end function
 
 logical pure function z0_nan(x) result (isnan)
  double complex,intent(in) :: x
  isnan=.not.( ((real(x).lt.0.d0).or.(real(x).ge.0.d0)) .and. ((imag(x).lt.0.d0).or.(imag(x).ge.0.d0)) )
 end function
 logical pure function z1_nan(x) result (isnan)
  double complex,intent(in) :: x(:)
  integer :: n,i
  n=size(x)
  do i=1,n
   isnan=.not.( ((real(x(i)).lt.0.d0).or.(real(x(i)).ge.0.d0)) .and. ((imag(x(i)).lt.0.d0).or.(imag(x(i)).ge.0.d0)) )
   if(isnan)return
  end do
 end function
 logical pure function z2_nan(x) result (isnan)
  double complex,intent(in) :: x(:,:)
  integer :: m,n,i,j
  m=size(x,1);n=size(x,2)
  do j=1,n
   do i=1,m
    isnan=.not.( ((real(x(i,j)).lt.0.d0).or.(real(x(i,j)).ge.0.d0)) .and. ((imag(x(i,j)).lt.0.d0).or.(imag(x(i,j)).ge.0.d0)) )
    if(isnan)return
   end do 
  end do
 end function
 logical pure function z3_nan(x) result (isnan)
  double complex,intent(in) :: x(:,:,:)
  integer :: m,n,p,i,j,k
  m=size(x,1);n=size(x,2);p=size(x,3)
  do k=1,p
   do j=1,n
    do i=1,m
     isnan=.not.( ((real(x(i,j,k)).lt.0.d0).or.(real(x(i,j,k)).ge.0.d0)) .and. ((imag(x(i,j,k)).lt.0.d0).or.(imag(x(i,j,k)).ge.0.d0)) )
     if(isnan)return
    end do 
   end do 
  end do
 end function
 
end module
