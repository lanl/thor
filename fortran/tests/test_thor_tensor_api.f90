!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                     LANL 04/21/2023
!!                                           T-3
!!                                Ismael Djibrilla Boureima [IDB]
!!
!!_________________________________________________________________________________________
!!   SCOP:  Test program 01:  This is to illustrate how to use the thor tensor object
!!          defined below, allong with usefule auxiliary procedures:
!!_________________________________________________________________________________________
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
character(len=*), parameter :: fname_in = 'tests/testdata/test01_small.dat'

  call test_tt_tensor_api(fname_in)

contains

  subroutine test_tt_tensor_api(fname_input_array)
  use thor_lib, only  : dtt_tensor, operator(-), pprint
  use thorio_lib, only : read_array_from_file
  use time_lib, only           : timef
  implicit none
  character(*), intent(in) :: fname_input_array
  type(dtt_tensor) :: v, arr_tt, crr_tt, d_tt
  integer l, failed, nx, ny, nz, res_total, d, info
  integer, allocatable :: n(:), r(:)
  double precision           :: t1,t2, dt, err
  double precision, parameter :: tol = 1d-12
  double precision, external  :: dnrm2
  double precision, allocatable :: arr(:,:,:), crr(:)

     print '("")'
     print '("Test 03: Reading a tensor object from a file")'
     print '("")'
     failed= 0

     ! Define tt metadata
     l = 1
     n = [3, 3, 3]
     r = [1, 2, 2, 1]

     print '("-- [3-1] testing:  v1 = dtt_tensor (empty):")'
     v = dtt_tensor(n, r_=r)
     call v% say
     print*,"[+] Exec time = ", dt,'s'
     err= dnrm2(size(v% full()), v% full(), 1)
     print '("|v|_2 = ",ES12.5)', err
     if (err.lt.tol) then
        print '("[+][TEST03-1]: empty tensor v1 created: PASS")'
     else
        print '("[-][TEST03-1]: tensor v1 not created: FAIL")'
        failed= failed + 1
     end if

     !!! {WIP} not pretty at all: produces SEGFAULT
     call pprint(v, label_="dtt_tensor v:")

     print '("-- [3-2] testing: read 3D array from file:")'
     call read_array_from_file(fname_input_array, nx,ny,nz, arr)
     !dt = timef() - t1
     print*,"[+] Exec time = ", dt,'s'
     if (nx.ne.128.or.ny.ne.128.or.nz.ne.128) failed= failed + 1

     if(failed.eq.0) then
        print '("[+][TEST03]: PASSED")'
     else
        print '("[-][TEST03]: FAILED with ",I2," error(s)")', failed
        error stop failed
     end if
     print '("-- [1-3] testing: Contruct TT from array read above:")'
     !Test constructor with 3D array
     !t1 = timef()
     arr_tt = dtt_tensor(arr, eps_=tol)
     !dt = timef() - t1
     print*,"arr_tt:"
     call arr_tt% say
     print*,"[+] Exec time = ", dt,'s'
     ! Test constructor with MultiD array
     res_total = nx*ny*nz
     n = (/nx,ny,nz/)
     allocate(crr(res_total), stat=info)
     crr=reshape(arr,(/res_total/))
     crr_tt = dtt_tensor(crr, n, eps_=tol)
     print*,"crr_tt:"
     call arr_tt% say

     !t1 = timef()
     d_tt = arr_tt - crr_tt
     !dt = timef() - t1
     print*,"d_tt = arr_tt - crr_tt:"
     call arr_tt% say
     err = sum(d_tt% full())
     print*,"[+] Exec time = ", dt,'s'
     if(err .le. tol) then
        print '("[+][TEST01][1-3]: Contruct TT from array PASSED")'
     else
        print '("[-][TEST01][1-3]: FAILED to accuratly construct TT from array with ",I2," error(s)")', failed
        error stop failed
     end if
     print*,"[+] err=", err
  end subroutine test_tt_tensor_api

end program
