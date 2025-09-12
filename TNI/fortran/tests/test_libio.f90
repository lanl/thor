program test_thor_io_lib
use mat_lib, only: norm => normfro
use thor_lib, only: dtt_tensor, sayfull, pprint, dtt_tensor_rand, operator(-)
use thorio_lib, only: read_array_from_file, dtt_read_sdv_file, dtt_write_sdv_file, &
      dtt_read_ascii_file, dtt_write_ascii_file
implicit none
integer                :: info, nx, ny, nz !, sz, He, Vac, i,j,k
double precision, allocatable :: arr(:,:,:), arr1(:)
integer, allocatable :: nn(:)
double precision, allocatable :: nrm(:)
type(dtt_tensor) :: tt, tt2, dtt
character(:), allocatable :: sdv_fname
double precision, parameter :: err_tol = 1d-14
double precision :: err, x
integer :: failed

   print '(/,"[test_thor_io] Testing THOR I/O functions",/)'
   failed= 0

   print '(/,"--- read_array_from_file(file test01_small.dat)")'
   call read_array_from_file('tests/testdata/test01_small.dat', nx,ny,nz, arr)
   
   nn = [nx, ny, nz]
   tt = dtt_tensor(reshape(arr, [product(nn)]), nn, eps_=1d-12)

   call tt% say
   call sayfull(tt)
   call pprint(tt)

   arr1 = tt% full() 
   err= dabs(norm(arr1 - reshape(arr,[product(nn)])) / norm(arr))
   print '("    Norm of the difference between the arrays: ", ES14.7)', err
   deallocate(arr, arr1)
   if (err.lt.err_tol) then
      print '("    PASS")'
   else
      print '("    FAIL")'
      failed= failed + 1
   endif

   print '(/,"--- Write an ASCII file, read it back and compare:")'
   tt = dtt_tensor_rand([12,10,11,2], r_=[1,2,4,3,1])
   call dtt_write_ascii_file(tt, '/tmp/tt_rand.dat', verb_=.true.)

   tt2 = dtt_read_ascii_file('/tmp/tt_rand.dat',verb_=.true.)
   dtt = tt - tt2
   print '("read it back and compute the difference: dtt = tt - tt2")'
   call dtt% say
   err= dabs(dtt% normb())
   print '("difference: |tt - tt2|_2:   ",ES14.7)', err
   if (err.lt.err_tol) then
      print '("    PASS")'
   else
      print '("    FAIL")'
      failed= failed + 1
   endif

   print '(/,"--- Write an svd file, read it back and compare:")'
   call dtt_write_sdv_file(tt, '/tmp/tt_rand.sdv')
   tt2 = dtt_read_sdv_file('/tmp/tt_rand.sdv')
   dtt = tt - tt2
   print '("tensor dtt = tt - tt2:")'
   call dtt% say
   err= dabs(dtt% normb())
   print '("|tt - tt2|_2: ",ES14.7)', err
   if (err.lt.err_tol) then
      print '("    PASS")'
   else
      print '("    FAIL")'
      failed= failed + 1
   endif

   sdv_fname = 'tests/testdata/b.sdv'
   print '(/,"--- Open an SVD file ",A,":")', sdv_fname
   tt = dtt_read_sdv_file(sdv_fname)
   call tt% say
   print '("    PASS")'

   sdv_fname = 'tests/testdata/A.sdv'
   print '(/,"--- Open an SVD file ",A,":")', sdv_fname
   tt = dtt_read_sdv_file(sdv_fname)
   call tt% say
   print '("    PASS")'

   if(failed.eq.0) then
      print '(/,"[+][test_thor_io]: PASSED")'
   else
      print '(/,"[-][test_thor_io]: FAILED with ",I2," error(s)")', failed
      error stop failed
   end if

end program test_thor_io_lib
