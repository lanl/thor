!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                     LANL 04/21/2023
!!                                           T-3
!!                                Ismael Djibrilla Boureima [IDB]
!!
!!_________________________________________________________________________________________
!!   SCOP:  Test program 02:  This is to illustrate how to use the thor tensor object
!!                            Testing the mpi/nccl communicator
!!_________________________________________________________________________________________
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main

  call test_thor_mpi_comm()
  call MPI_FINALIZE(err)

contains

  subroutine test_thor_mpi_comm()
  use distributed_comms !, only  : super_comm
  use communicators, ONLY: mpi_global_comm, mpi_sub_comm
  use time_lib, only           : timef
  implicit none
  integer                      :: err, failed
  type(super_comm)             :: cpu_comm
  type(mpi_global_comm)        :: gcomm
  type(mpi_sub_comm)           :: lcomm
  double precision             :: t1,t2, dt
  double precision, parameter  :: tol = 1d-12
     print*, '("")'
     print*, '("[MPI-COMM][Test][01]: INIT")'
     print*, '("")'
     failed= 0
     err  = -1
     ! Define tt metadata
     print*, '("-- [1-1] testing: Class intenciation:")'
     t1=timef()
     cpu_comm = super_comm()
     dt = timef() - t1
     !print*,"[+]","[",TRIM(ADJUSTL(cpu_comm%host_name)),"][MPI-COMM] INSTENCIATED OK: Exec time = ", dt,'s'
     print*,"   Global rank: rank                                     = ",cpu_comm%rank
     print*,"   Size of global communicator: size                     = ",cpu_comm%comm_size
     print*,"   Local rank: local_rank                                = ",cpu_comm%local_rank
     print*,"   Size of local communicator: local_comm_size           = ",cpu_comm%local_comm_size 
     print*,"   Cuda device in use: cuda_device                       = ",cpu_comm%cuda_device
     print*,"   Flag for curent proc in contex: in_gpu_ctx           = ",cpu_comm%in_gpu_ctx
     print*,"   Flag to check wether comm is initialized: initialized = ",cpu_comm%mpi_comm_init
     print '("-- [1-2] testing: COMM Init:")'
     t1=timef()
     err = cpu_comm % init_mpi()
     dt = timef() - t1
     if (err==0) then
        ! print*,"[+]","[",TRIM(ADJUSTL(cpu_comm%host_name)),"][MPI-COMM] INITIALIZED OK: Exec time = ", dt,'s'
        print*,"OK"
     else
         !print*,"[-]","[",TRIM(ADJUSTL(cpu_comm%host_name)),"][MPI-COMM] UNABLE TO INITIALIZE MPI COMM: Exec time = ", dt,'s'
         print*,"NOT OK"
     endif
      print*,"   Global rank: rank                                     = ",cpu_comm%rank
     print*,"   Size of global communicator: size                     = ",cpu_comm%comm_size
     print*,"   Local rank: local_rank                                = ",cpu_comm%local_rank
     print*,"   Size of local communicator: local_comm_size           = ",cpu_comm%local_comm_size
     print*,"   Cuda device in use: cuda_device                       = ",cpu_comm%cuda_device
     print*,"   Flag for curent proc in contex: in_gpu_ctx           = ",cpu_comm%in_gpu_ctx
     print*,"   Flag to check wether comm is initialized: initialized = ",cpu_comm%mpi_comm_init
     t1=timef()
     err = cpu_comm % init_cuda()
     if (err==0) then
         !print*,"[+]","[",TRIM(ADJUSTL(cpu_comm%host_name)),"][MPI-COMM] CUDA-BACKEND INITIALIZED OK: Exec time = ", dt,'s'
         print*,"OK"
     else
         !print*,"[-]","[",TRIM(ADJUSTL(cpu_comm%host_name)),"][MPI-COMM] UNABLE TO INITIALIZE CUDA_BACKEND: Exec time = ", dt,'s'
         print*,"NOT OK"
     endif
      print*,"   Global rank: rank                                     = ",cpu_comm%rank
     print*,"   Size of global communicator: size                     = ",cpu_comm%comm_size
     print*,"   Local rank: local_rank                                = ",cpu_comm%local_rank
     print*,"   Size of local communicator: local_comm_size           = ",cpu_comm%local_comm_size
     print*,"   Cuda device in use: cuda_device                       = ",cpu_comm%cuda_device
     print*,"   Flag for curent proc in contex: in_gpu_ctx           = ",cpu_comm%in_gpu_ctx
     print*,"   Flag to check wether comm is initialized: initialized = ",cpu_comm%mpi_comm_init

     err = cpu_comm%init_nccl_comm()
  end subroutine test_thor_mpi_comm

end program
