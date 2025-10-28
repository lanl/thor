program main
use cuda_comm_aux
use string_lib, only: str
use mpi_f08
implicit none
integer, target          :: rank, commSize, localRank, localCommSize, ctxRank,&
                            ctxCommSize,  nDevLocal, color
integer                  :: ctx_color, VRBZ, err
type(MPI_COMM), target   :: localComm, ctxComm
character(len=128) :: prefix

  ctx_color=11; VRBZ=3
  call get_hw_topology(ctx_color, rank, commSize, localComm, localRank, &
                         localCommSize, ctxComm, ctxRank, ctxCommSize, nDevLocal, &
                         color, VRBZ_=VRBZ)
  prefix = '[L'//str(localRank)//'/G'//str(rank)//']'
  print '("'//trim(prefix)//': rank='//str(rank)// &
          '| size='//str(commSize)// &
          '| localRank='//str(localRank) // &
          '| localCommSize='//str(localCommSize)//'")'
  call MPI_BARRIER(MPI_COMM_WORLD)

  call cuda_comm_init(rank, commSize, localRank, VRBZ)
  print '("'//trim(prefix)//' CUDA communicator(s) initialized OK")'
  call MPI_BARRIER(MPI_COMM_WORLD)

  print '("'//trim(prefix)//' CAL communicator initialized: ",L1)', is_cal_initialized
  print '("'//trim(prefix)//' NCCL communicator initialized: ",L1)', is_nccl_initialized
 
  call cal_comm_destroy
  print '("'//trim(prefix)//' CAL communicator destroyed OK")'
  call MPI_BARRIER(MPI_COMM_WORLD)

  call MPI_FINALIZE(err)
end program main
