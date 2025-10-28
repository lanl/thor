module communicators
 implicit none
  USE mpi_f08
 type,public             :: mpi_comm               ! modern mpi communicator
   integer               :: rank                  = -1  ! mpi rank
   integer               :: comm_size             = -2  ! size of communicator
   integer               :: comm                  = -3
   integer               :: local_rank            = -4  ! local mpi rank
   integer               :: local_comm_size       = -5  ! size of local communicator
   integer               :: local_comm            = -6  ! Local comunicator ref
   integer               :: cuda_device           = -7  ! Cuda device
   integer               :: n_local_cuda_device   = -8  ! Local number of cuda GPUs
   integer               :: n_global_cuda_device  = -9  ! Global number of cuda GPUs
   integer               :: ctx_comm              = -10 ! Inter-node CTX communicator
   integer               :: ctx_local_comm        = -11 ! Intra-node CTX communicator
   integer               :: ctx_rank              = -12 ! rank on inter-nodr ctx communicator
   integer               :: ctx_local_rank        = -13 ! rank on intra-nodr ctx communicator
   integer               :: ctx_size              = -14 ! rank on inter-nodr ctx communicator
   integer               :: local_ctx_size        = -15 ! rank on intra-nodr ctx communicator
   logical               :: in_cuda_ctx           = .false.
   logical               :: mpi_comm_init         = .false.
   logical               :: cuda_backend_init     = .false.
   logical               :: local_nccl_comm_init  = .false.
   logical               :: global_nccl_comm_init = .false.
   type(ncclUniqueId)    :: global_nccl_id
   type(ncclUniqueId)    :: local_nccl_id
   type(ncclComm)        :: global_nccl_comm
   type(ncclComm)        :: local_nccl_comm
 contains
  procedure       :: init_mpi         => init_mpi_comm
  procedure       :: init_cuda        => init_cuda_back_end
  procedure       :: init_nccl_comm   => buld_global_nccl_comm
  procedure       :: get_rank         => get_comm_rank
  procedure       :: get_size         => get_comm_size
  procedure       :: barrier          => comm_barrier
  procedure       :: finalize         => finalize_comm
  procedure       :: get_mpi_op       => get_mpi_op
  procedure       :: translate_mpi_op => get_mpi_op
  !procedure       :: send             => global_comm_send
  !procedure       :: recv             => global_comm_recv
  !procedure       :: bcast            => global_comm_bcast
  !procedure       :: all_reduce       => global_comm_all_reduce
  final           :: dealloc_comm
 end type mpi_comm
end module communicators
