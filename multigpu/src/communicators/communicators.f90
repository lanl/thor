module communicators
  USE mpi_f08 , ONLY        : MPI_COMM, MPI_INIT, MPI_COMM_FREE, MPI_FINALIZE, MPI_Initialized, MPI_COMM_SPLIT, MPI_COMM_SPLIT_TYPE, MPI_COMM_WORLD,&
                             MPI_COMM_TYPE_SHARED, MPI_INFO_NULL, MPI_BARRIER
  USE comm_allreduce, ONLY : ALLREDUCE_ => allreduce
  USE comm_bcast, ONLY     : BCAST_     => bcast
  USE comm_recv, ONLY      : RECV_      => recv
  !USE comm_send, ONLY      : SEND_      => send
  USE nccl, only          : ncclUniqueId, ncclResult
  implicit none
 type,public               :: mpi_global_comm       ! global mpi communicator
   type(MPI_COMM)          :: comm
   integer                 :: rank                  = -1  ! mpi rank
   integer                 :: comm_size             = -2  ! size of communicator
   logical                 :: is_initialized        = .false.
   character*16            :: hostname

 contains
  procedure       :: sync             => sync_global_comm
  procedure       :: barrier          => sync_global_comm
  procedure       :: split            => split_global_comm
  procedure       :: finalize         => finalize_mpi_global_comm
  !final           ::                     close_mpi_global_comm
  procedure, pass(this)  :: mpi_global_comm_send_scal_i2,mpi_global_comm_recv_scal_i2,&
                            mpi_global_comm_bcast_scal_i2,mpi_global_comm_allreduce_scal_i2, &
                            mpi_global_comm_send_scal_i4,mpi_global_comm_recv_scal_i4,&
                            mpi_global_comm_bcast_scal_i4,mpi_global_comm_allreduce_scal_i4, &
                            mpi_global_comm_send_scal_i8,mpi_global_comm_recv_scal_i8,&
                            mpi_global_comm_bcast_scal_i8,mpi_global_comm_allreduce_scal_i8, &
                            mpi_global_comm_send_scal_r4,mpi_global_comm_recv_scal_r4,&
                            mpi_global_comm_bcast_scal_r4,mpi_global_comm_allreduce_scal_r4, &
                            mpi_global_comm_send_scal_r8,mpi_global_comm_recv_scal_r8,&
                            mpi_global_comm_bcast_scal_r8,mpi_global_comm_allreduce_scal_r8, &
                            mpi_global_comm_send_vect_i2,mpi_global_comm_recv_vect_i2,&
                            mpi_global_comm_bcast_vect_i2,mpi_global_comm_allreduce_vect_i2, &
                            mpi_global_comm_send_vect_i4,mpi_global_comm_recv_vect_i4,&
                            mpi_global_comm_bcast_vect_i4,mpi_global_comm_allreduce_vect_i4, &
                            mpi_global_comm_send_vect_i8,mpi_global_comm_recv_vect_i8,&
                            mpi_global_comm_bcast_vect_i8,mpi_global_comm_allreduce_vect_i8, &
                            mpi_global_comm_send_vect_r4,mpi_global_comm_recv_vect_r4,&
                            mpi_global_comm_bcast_vect_r4,mpi_global_comm_allreduce_vect_r4, &
                            mpi_global_comm_send_vect_r8,mpi_global_comm_recv_vect_r8,&
                            mpi_global_comm_bcast_vect_r8,mpi_global_comm_allreduce_vect_r8, &
                            mpi_global_comm_send_mat_i2,mpi_global_comm_recv_mat_i2,&
                            mpi_global_comm_bcast_mat_i2,mpi_global_comm_allreduce_mat_i2, &
                            mpi_global_comm_send_mat_i4,mpi_global_comm_recv_mat_i4,&
                            mpi_global_comm_bcast_mat_i4,mpi_global_comm_allreduce_mat_i4, &
                            mpi_global_comm_send_mat_i8,mpi_global_comm_recv_mat_i8,&
                            mpi_global_comm_bcast_mat_i8,mpi_global_comm_allreduce_mat_i8, &
                            mpi_global_comm_send_mat_r4,mpi_global_comm_recv_mat_r4,&
                            mpi_global_comm_bcast_mat_r4,mpi_global_comm_allreduce_mat_r4, &
                            mpi_global_comm_send_mat_r8,mpi_global_comm_recv_mat_r8,&
                            mpi_global_comm_bcast_mat_r8,mpi_global_comm_allreduce_mat_r8, &
                            mpi_global_comm_send_vol_i2,mpi_global_comm_recv_vol_i2,&
                            mpi_global_comm_bcast_vol_i2,mpi_global_comm_allreduce_vol_i2, &
                            mpi_global_comm_send_vol_i4,mpi_global_comm_recv_vol_i4,&
                            mpi_global_comm_bcast_vol_i4,mpi_global_comm_allreduce_vol_i4, &
                            mpi_global_comm_send_vol_i8,mpi_global_comm_recv_vol_i8,&
                            mpi_global_comm_bcast_vol_i8,mpi_global_comm_allreduce_vol_i8, &
                            mpi_global_comm_send_vol_r4,mpi_global_comm_recv_vol_r4,&
                            mpi_global_comm_bcast_vol_r4,mpi_global_comm_allreduce_vol_r4, &
                            mpi_global_comm_send_vol_r8,mpi_global_comm_recv_vol_r8,&
                            mpi_global_comm_bcast_vol_r8,mpi_global_comm_allreduce_vol_r8
  generic                :: send => mpi_global_comm_send_scal_i2, &
                                      mpi_global_comm_send_scal_i4, &
                                      mpi_global_comm_send_scal_i8, &
                                      mpi_global_comm_send_scal_r4, &
                                      mpi_global_comm_send_scal_r8, &
                                      mpi_global_comm_send_vect_i2, &
                                      mpi_global_comm_send_vect_i4, &
                                      mpi_global_comm_send_vect_i8, &
                                      mpi_global_comm_send_vect_r4, &
                                      mpi_global_comm_send_vect_r8, &
                                      mpi_global_comm_send_mat_i2, &
                                      mpi_global_comm_send_mat_i4, &
                                      mpi_global_comm_send_mat_i8, &
                                      mpi_global_comm_send_mat_r4, &
                                      mpi_global_comm_send_mat_r8, &
                                      mpi_global_comm_send_vol_i2, &
                                      mpi_global_comm_send_vol_i4, &
                                      mpi_global_comm_send_vol_i8, &
                                      mpi_global_comm_send_vol_r4, &
                                      mpi_global_comm_send_vol_r8
  generic                :: recv => mpi_global_comm_recv_scal_i2, &
                                       mpi_global_comm_recv_scal_i4, &
                                       mpi_global_comm_recv_scal_i8, &
                                       mpi_global_comm_recv_scal_r4, &
                                       mpi_global_comm_recv_scal_r8, &
                                       mpi_global_comm_recv_vect_i2, &
                                       mpi_global_comm_recv_vect_i4, &
                                       mpi_global_comm_recv_vect_i8, &
                                       mpi_global_comm_recv_vect_r4, &
                                       mpi_global_comm_recv_vect_r8, &
                                       mpi_global_comm_recv_mat_i2, &
                                       mpi_global_comm_recv_mat_i4, &
                                       mpi_global_comm_recv_mat_i8, &
                                       mpi_global_comm_recv_mat_r4, &
                                       mpi_global_comm_recv_mat_r8, &
                                       mpi_global_comm_recv_vol_i2, &
                                       mpi_global_comm_recv_vol_i4, &
                                       mpi_global_comm_recv_vol_i8, &
                                       mpi_global_comm_recv_vol_r4, &
                                       mpi_global_comm_recv_vol_r8
  generic                :: bcast => mpi_global_comm_bcast_scal_i2, &
                                        mpi_global_comm_bcast_scal_i4, &
                                        mpi_global_comm_bcast_scal_i8, &
                                        mpi_global_comm_bcast_scal_r4, &
                                        mpi_global_comm_bcast_scal_r8, &
                                        mpi_global_comm_bcast_vect_i2, &
                                        mpi_global_comm_bcast_vect_i4, &
                                        mpi_global_comm_bcast_vect_i8, &
                                        mpi_global_comm_bcast_vect_r4, &
                                        mpi_global_comm_bcast_vect_r8, &
                                        mpi_global_comm_bcast_mat_i2, &
                                        mpi_global_comm_bcast_mat_i4, &
                                        mpi_global_comm_bcast_mat_i8, &
                                        mpi_global_comm_bcast_mat_r4, &
                                        mpi_global_comm_bcast_mat_r8, &
                                        mpi_global_comm_bcast_vol_i2, &
                                        mpi_global_comm_bcast_vol_i4, &
                                        mpi_global_comm_bcast_vol_i8, &
                                        mpi_global_comm_bcast_vol_r4, &
                                        mpi_global_comm_bcast_vol_r8
  generic                :: allreduce => mpi_global_comm_allreduce_scal_i2, &
                                            mpi_global_comm_allreduce_scal_i4, &
                                            mpi_global_comm_allreduce_scal_i8, &
                                            mpi_global_comm_allreduce_scal_r4, &
                                            mpi_global_comm_allreduce_scal_r8, &
                                            mpi_global_comm_allreduce_vect_i2, &
                                            mpi_global_comm_allreduce_vect_i4, &
                                            mpi_global_comm_allreduce_vect_i8, &
                                            mpi_global_comm_allreduce_vect_r4, &
                                            mpi_global_comm_allreduce_vect_r8, &
                                            mpi_global_comm_allreduce_mat_i2, &
                                            mpi_global_comm_allreduce_mat_i4, &
                                            mpi_global_comm_allreduce_mat_i8, &
                                            mpi_global_comm_allreduce_mat_r4, &
                                            mpi_global_comm_allreduce_mat_r8, &
                                            mpi_global_comm_allreduce_vol_i2, &
                                            mpi_global_comm_allreduce_vol_i4, &
                                            mpi_global_comm_allreduce_vol_i8, &
                                            mpi_global_comm_allreduce_vol_r4, &
                                            mpi_global_comm_allreduce_vol_r8
 end type mpi_global_comm
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type,public               :: mpi_sub_comm                    ! modern mpi communicator
   type(MPI_COMM)          :: comm
   integer                 :: rank                  = -1  ! mpi rank
   integer                 :: comm_size             = -2  ! size of communicator
   integer                 :: color                 = -4  ! Sub_comm group color
   logical                 :: is_initialized        = .false.
   character*16            :: hostname
 contains
  procedure       :: sync             => sync_sub_comm
  procedure       :: barrier          => sync_sub_comm
  procedure       :: split            => split_sub_comm
  procedure       :: finalize         => finalize_mpi_sub_comm
  !final           ::                     close_mpi_sub_comm
  procedure, pass(this)  :: mpi_sub_comm_send_scal_i2,mpi_sub_comm_recv_scal_i2,&
                            mpi_sub_comm_bcast_scal_i2,mpi_sub_comm_allreduce_scal_i2, &
                            mpi_sub_comm_send_scal_i4,mpi_sub_comm_recv_scal_i4,&
                            mpi_sub_comm_bcast_scal_i4,mpi_sub_comm_allreduce_scal_i4, &
                            mpi_sub_comm_send_scal_i8,mpi_sub_comm_recv_scal_i8,&
                            mpi_sub_comm_bcast_scal_i8,mpi_sub_comm_allreduce_scal_i8, &
                            mpi_sub_comm_send_scal_r4,mpi_sub_comm_recv_scal_r4,&
                            mpi_sub_comm_bcast_scal_r4,mpi_sub_comm_allreduce_scal_r4, &
                            mpi_sub_comm_send_scal_r8,mpi_sub_comm_recv_scal_r8,&
                            mpi_sub_comm_bcast_scal_r8,mpi_sub_comm_allreduce_scal_r8, &
                            mpi_sub_comm_send_vect_i2,mpi_sub_comm_recv_vect_i2,&
                            mpi_sub_comm_bcast_vect_i2,mpi_sub_comm_allreduce_vect_i2, &
                            mpi_sub_comm_send_vect_i4,mpi_sub_comm_recv_vect_i4,&
                            mpi_sub_comm_bcast_vect_i4,mpi_sub_comm_allreduce_vect_i4, &
                            mpi_sub_comm_send_vect_i8,mpi_sub_comm_recv_vect_i8,&
                            mpi_sub_comm_bcast_vect_i8,mpi_sub_comm_allreduce_vect_i8, &
                            mpi_sub_comm_send_vect_r4,mpi_sub_comm_recv_vect_r4,&
                            mpi_sub_comm_bcast_vect_r4,mpi_sub_comm_allreduce_vect_r4, &
                            mpi_sub_comm_send_vect_r8,mpi_sub_comm_recv_vect_r8,&
                            mpi_sub_comm_bcast_vect_r8,mpi_sub_comm_allreduce_vect_r8, &
                            mpi_sub_comm_send_mat_i2,mpi_sub_comm_recv_mat_i2,&
                            mpi_sub_comm_bcast_mat_i2,mpi_sub_comm_allreduce_mat_i2, &
                            mpi_sub_comm_send_mat_i4,mpi_sub_comm_recv_mat_i4,&
                            mpi_sub_comm_bcast_mat_i4,mpi_sub_comm_allreduce_mat_i4, &
                            mpi_sub_comm_send_mat_i8,mpi_sub_comm_recv_mat_i8,&
                            mpi_sub_comm_bcast_mat_i8,mpi_sub_comm_allreduce_mat_i8, &
                            mpi_sub_comm_send_mat_r4,mpi_sub_comm_recv_mat_r4,&
                            mpi_sub_comm_bcast_mat_r4,mpi_sub_comm_allreduce_mat_r4, &
                            mpi_sub_comm_send_mat_r8,mpi_sub_comm_recv_mat_r8,&
                            mpi_sub_comm_bcast_mat_r8,mpi_sub_comm_allreduce_mat_r8, &
                            mpi_sub_comm_send_vol_i2,mpi_sub_comm_recv_vol_i2,&
                            mpi_sub_comm_bcast_vol_i2,mpi_sub_comm_allreduce_vol_i2, &
                            mpi_sub_comm_send_vol_i4,mpi_sub_comm_recv_vol_i4,&
                            mpi_sub_comm_bcast_vol_i4,mpi_sub_comm_allreduce_vol_i4, &
                            mpi_sub_comm_send_vol_i8,mpi_sub_comm_recv_vol_i8,&
                            mpi_sub_comm_bcast_vol_i8,mpi_sub_comm_allreduce_vol_i8, &
                            mpi_sub_comm_send_vol_r4,mpi_sub_comm_recv_vol_r4,&
                            mpi_sub_comm_bcast_vol_r4,mpi_sub_comm_allreduce_vol_r4, &
                            mpi_sub_comm_send_vol_r8,mpi_sub_comm_recv_vol_r8,&
                            mpi_sub_comm_bcast_vol_r8,mpi_sub_comm_allreduce_vol_r8
  generic                :: send => mpi_sub_comm_send_scal_i2, &
                                      mpi_sub_comm_send_scal_i4, &
                                      mpi_sub_comm_send_scal_i8, &
                                      mpi_sub_comm_send_scal_r4, &
                                      mpi_sub_comm_send_scal_r8, &
                                      mpi_sub_comm_send_vect_i2, &
                                      mpi_sub_comm_send_vect_i4, &
                                      mpi_sub_comm_send_vect_i8, &
                                      mpi_sub_comm_send_vect_r4, &
                                      mpi_sub_comm_send_vect_r8, &
                                      mpi_sub_comm_send_mat_i2, &
                                      mpi_sub_comm_send_mat_i4, &
                                      mpi_sub_comm_send_mat_i8, &
                                      mpi_sub_comm_send_mat_r4, &
                                      mpi_sub_comm_send_mat_r8, &
                                      mpi_sub_comm_send_vol_i2, &
                                      mpi_sub_comm_send_vol_i4, &
                                      mpi_sub_comm_send_vol_i8, &
                                      mpi_sub_comm_send_vol_r4, &
                                      mpi_sub_comm_send_vol_r8
  generic                :: recv => mpi_sub_comm_recv_scal_i2, &
                                       mpi_sub_comm_recv_scal_i4, &
                                       mpi_sub_comm_recv_scal_i8, &
                                       mpi_sub_comm_recv_scal_r4, &
                                       mpi_sub_comm_recv_scal_r8, &
                                       mpi_sub_comm_recv_vect_i2, &
                                       mpi_sub_comm_recv_vect_i4, &
                                       mpi_sub_comm_recv_vect_i8, &
                                       mpi_sub_comm_recv_vect_r4, &
                                       mpi_sub_comm_recv_vect_r8, &
                                       mpi_sub_comm_recv_mat_i2, &
                                       mpi_sub_comm_recv_mat_i4, &
                                       mpi_sub_comm_recv_mat_i8, &
                                       mpi_sub_comm_recv_mat_r4, &
                                       mpi_sub_comm_recv_mat_r8, &
                                       mpi_sub_comm_recv_vol_i2, &
                                       mpi_sub_comm_recv_vol_i4, &
                                       mpi_sub_comm_recv_vol_i8, &
                                       mpi_sub_comm_recv_vol_r4, &
                                       mpi_sub_comm_recv_vol_r8
  generic                :: bcast => mpi_sub_comm_bcast_scal_i2, &
                                        mpi_sub_comm_bcast_scal_i4, &
                                        mpi_sub_comm_bcast_scal_i8, &
                                        mpi_sub_comm_bcast_scal_r4, &
                                        mpi_sub_comm_bcast_scal_r8, &
                                        mpi_sub_comm_bcast_vect_i2, &
                                        mpi_sub_comm_bcast_vect_i4, &
                                        mpi_sub_comm_bcast_vect_i8, &
                                        mpi_sub_comm_bcast_vect_r4, &
                                        mpi_sub_comm_bcast_vect_r8, &
                                        mpi_sub_comm_bcast_mat_i2, &
                                        mpi_sub_comm_bcast_mat_i4, &
                                        mpi_sub_comm_bcast_mat_i8, &
                                        mpi_sub_comm_bcast_mat_r4, &
                                        mpi_sub_comm_bcast_mat_r8, &
                                        mpi_sub_comm_bcast_vol_i2, &
                                        mpi_sub_comm_bcast_vol_i4, &
                                        mpi_sub_comm_bcast_vol_i8, &
                                        mpi_sub_comm_bcast_vol_r4, &
                                        mpi_sub_comm_bcast_vol_r8
  generic                :: allreduce => mpi_sub_comm_allreduce_scal_i2, &
                                            mpi_sub_comm_allreduce_scal_i4, &
                                            mpi_sub_comm_allreduce_scal_i8, &
                                            mpi_sub_comm_allreduce_scal_r4, &
                                            mpi_sub_comm_allreduce_scal_r8, &
                                            mpi_sub_comm_allreduce_vect_i2, &
                                            mpi_sub_comm_allreduce_vect_i4, &
                                            mpi_sub_comm_allreduce_vect_i8, &
                                            mpi_sub_comm_allreduce_vect_r4, &
                                            mpi_sub_comm_allreduce_vect_r8, &
                                            mpi_sub_comm_allreduce_mat_i2, &
                                            mpi_sub_comm_allreduce_mat_i4, &
                                            mpi_sub_comm_allreduce_mat_i8, &
                                            mpi_sub_comm_allreduce_mat_r4, &
                                            mpi_sub_comm_allreduce_mat_r8, &
                                            mpi_sub_comm_allreduce_vol_i2, &
                                            mpi_sub_comm_allreduce_vol_i4, &
                                            mpi_sub_comm_allreduce_vol_i8, &
                                            mpi_sub_comm_allreduce_vol_r4, &
                                            mpi_sub_comm_allreduce_vol_r8
 end type mpi_sub_comm
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![           NCCL         ]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type,public               :: nccl_sub_comm               ! modern NCCL communicator
   integer                 :: rank                  = -1  ! rank
   integer                 :: comm_size             = -2  ! size of communicator
   logical                 :: is_initialized        = .false.
   type(ncclUniqueId)      :: nccl_id
   type(ncclResult)        :: nccl_result
 contains
 end type nccl_sub_comm
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  [MODULE INTERFACES]
!
 !! CONSTRUCTOR\DESTRUCTOR
   interface mpi_global_comm ! global comm constructor interface
     procedure :: empty_mpi_global_comm
   end interface mpi_global_comm
   interface mpi_sub_comm ! sub_communicator constructor interface
     procedure :: local_mpi_comm, &                              !This should give the direct local sub_communicator
                  mpi_sub_comm_from_split, mpi_sub_comm_from_global_comm_split, mpi_sub_comm_from_sub_comm_split
   end interface mpi_sub_comm
 !! SEND
    interface send
      procedure mpi_global_comm_send_scal_i2, mpi_sub_comm_send_scal_i2, &
                        mpi_global_comm_send_scal_i4, mpi_sub_comm_send_scal_i4, &
                        mpi_global_comm_send_scal_i8, mpi_sub_comm_send_scal_i8, &
                        mpi_global_comm_send_scal_r4, mpi_sub_comm_send_scal_r4, &
                        mpi_global_comm_send_scal_r8, mpi_sub_comm_send_scal_r8, &
                        mpi_global_comm_send_vect_i2, mpi_sub_comm_send_vect_i2, &
                        mpi_global_comm_send_vect_i4, mpi_sub_comm_send_vect_i4, &
                        mpi_global_comm_send_vect_i8, mpi_sub_comm_send_vect_i8, &
                        mpi_global_comm_send_vect_r4, mpi_sub_comm_send_vect_r4, &
                        mpi_global_comm_send_vect_r8, mpi_sub_comm_send_vect_r8, &
                        mpi_global_comm_send_mat_i2, mpi_sub_comm_send_mat_i2, &
                        mpi_global_comm_send_mat_i4, mpi_sub_comm_send_mat_i4, &
                        mpi_global_comm_send_mat_i8, mpi_sub_comm_send_mat_i8, &
                        mpi_global_comm_send_mat_r4, mpi_sub_comm_send_mat_r4, &
                        mpi_global_comm_send_mat_r8, mpi_sub_comm_send_mat_r8, &
                        mpi_global_comm_send_vol_i2, mpi_sub_comm_send_vol_i2, &
                        mpi_global_comm_send_vol_i4, mpi_sub_comm_send_vol_i4, &
                        mpi_global_comm_send_vol_i8, mpi_sub_comm_send_vol_i8, &
                        mpi_global_comm_send_vol_r4, mpi_sub_comm_send_vol_r4, &
                        mpi_global_comm_send_vol_r8, mpi_sub_comm_send_vol_r8
    end interface send
 !! RECV
    interface recv
      procedure mpi_global_comm_recv_scal_i2, mpi_sub_comm_recv_scal_i2, &
                        mpi_global_comm_recv_scal_i4, mpi_sub_comm_recv_scal_i4, &
                        mpi_global_comm_recv_scal_i8, mpi_sub_comm_recv_scal_i8, &
                        mpi_global_comm_recv_scal_r4, mpi_sub_comm_recv_scal_r4, &
                        mpi_global_comm_recv_scal_r8, mpi_sub_comm_recv_scal_r8, &
                        mpi_global_comm_recv_vect_i2, mpi_sub_comm_recv_vect_i2, &
                        mpi_global_comm_recv_vect_i4, mpi_sub_comm_recv_vect_i4, &
                        mpi_global_comm_recv_vect_i8, mpi_sub_comm_recv_vect_i8, &
                        mpi_global_comm_recv_vect_r4, mpi_sub_comm_recv_vect_r4, &
                        mpi_global_comm_recv_vect_r8, mpi_sub_comm_recv_vect_r8, &
                        mpi_global_comm_recv_mat_i2, mpi_sub_comm_recv_mat_i2, &
                        mpi_global_comm_recv_mat_i4, mpi_sub_comm_recv_mat_i4, &
                        mpi_global_comm_recv_mat_i8, mpi_sub_comm_recv_mat_i8, &
                        mpi_global_comm_recv_mat_r4, mpi_sub_comm_recv_mat_r4, &
                        mpi_global_comm_recv_mat_r8, mpi_sub_comm_recv_mat_r8, &
                        mpi_global_comm_recv_vol_i2, mpi_sub_comm_recv_vol_i2, &
                        mpi_global_comm_recv_vol_i4, mpi_sub_comm_recv_vol_i4, &
                        mpi_global_comm_recv_vol_i8, mpi_sub_comm_recv_vol_i8, &
                        mpi_global_comm_recv_vol_r4, mpi_sub_comm_recv_vol_r4, &
                        mpi_global_comm_recv_vol_r8, mpi_sub_comm_recv_vol_r8
    end interface recv
 !! BCAST
    interface bcast
      procedure mpi_global_comm_bcast_scal_i2, mpi_sub_comm_bcast_scal_i2, &
                        mpi_global_comm_bcast_scal_i4, mpi_sub_comm_bcast_scal_i4, &
                        mpi_global_comm_bcast_scal_i8, mpi_sub_comm_bcast_scal_i8, &
                        mpi_global_comm_bcast_scal_r4, mpi_sub_comm_bcast_scal_r4, &
                        mpi_global_comm_bcast_scal_r8, mpi_sub_comm_bcast_scal_r8, &
                        mpi_global_comm_bcast_vect_i2, mpi_sub_comm_bcast_vect_i2, &
                        mpi_global_comm_bcast_vect_i4, mpi_sub_comm_bcast_vect_i4, &
                        mpi_global_comm_bcast_vect_i8, mpi_sub_comm_bcast_vect_i8, &
                        mpi_global_comm_bcast_vect_r4, mpi_sub_comm_bcast_vect_r4, &
                        mpi_global_comm_bcast_vect_r8, mpi_sub_comm_bcast_vect_r8, &
                        mpi_global_comm_bcast_mat_i2, mpi_sub_comm_bcast_mat_i2, &
                        mpi_global_comm_bcast_mat_i4, mpi_sub_comm_bcast_mat_i4, &
                        mpi_global_comm_bcast_mat_i8, mpi_sub_comm_bcast_mat_i8, &
                        mpi_global_comm_bcast_mat_r4, mpi_sub_comm_bcast_mat_r4, &
                        mpi_global_comm_bcast_mat_r8, mpi_sub_comm_bcast_mat_r8, &
                        mpi_global_comm_bcast_vol_i2, mpi_sub_comm_bcast_vol_i2, &
                        mpi_global_comm_bcast_vol_i4, mpi_sub_comm_bcast_vol_i4, &
                        mpi_global_comm_bcast_vol_i8, mpi_sub_comm_bcast_vol_i8, &
                        mpi_global_comm_bcast_vol_r4, mpi_sub_comm_bcast_vol_r4, &
                        mpi_global_comm_bcast_vol_r8, mpi_sub_comm_bcast_vol_r8
    end interface bcast
 !! ALLREDUCE
    interface allreduce
      procedure mpi_global_comm_allreduce_scal_i2, mpi_sub_comm_allreduce_scal_i2, &
                        mpi_global_comm_allreduce_scal_i4, mpi_sub_comm_allreduce_scal_i4, &
                        mpi_global_comm_allreduce_scal_i8, mpi_sub_comm_allreduce_scal_i8, &
                        mpi_global_comm_allreduce_scal_r4, mpi_sub_comm_allreduce_scal_r4, &
                        mpi_global_comm_allreduce_scal_r8, mpi_sub_comm_allreduce_scal_r8, &
                        mpi_global_comm_allreduce_vect_i2, mpi_sub_comm_allreduce_vect_i2, &
                        mpi_global_comm_allreduce_vect_i4, mpi_sub_comm_allreduce_vect_i4, &
                        mpi_global_comm_allreduce_vect_i8, mpi_sub_comm_allreduce_vect_i8, &
                        mpi_global_comm_allreduce_vect_r4, mpi_sub_comm_allreduce_vect_r4, &
                        mpi_global_comm_allreduce_vect_r8, mpi_sub_comm_allreduce_vect_r8, &
                        mpi_global_comm_allreduce_mat_i2, mpi_sub_comm_allreduce_mat_i2, &
                        mpi_global_comm_allreduce_mat_i4, mpi_sub_comm_allreduce_mat_i4, &
                        mpi_global_comm_allreduce_mat_i8, mpi_sub_comm_allreduce_mat_i8, &
                        mpi_global_comm_allreduce_mat_r4, mpi_sub_comm_allreduce_mat_r4, &
                        mpi_global_comm_allreduce_mat_r8, mpi_sub_comm_allreduce_mat_r8, &
                        mpi_global_comm_allreduce_vol_i2, mpi_sub_comm_allreduce_vol_i2, &
                        mpi_global_comm_allreduce_vol_i4, mpi_sub_comm_allreduce_vol_i4, &
                        mpi_global_comm_allreduce_vol_i8, mpi_sub_comm_allreduce_vol_i8, &
                        mpi_global_comm_allreduce_vol_r4, mpi_sub_comm_allreduce_vol_r4, &
                        mpi_global_comm_allreduce_vol_r8, mpi_sub_comm_allreduce_vol_r8
    end interface allreduce

 contains
!
!  [CONSTRUCTORS]
   !!!!              GLOBAL
   function empty_mpi_global_comm() result(this)
   type(mpi_global_comm)          :: this
   !
   integer                        :: err, stat
   character(len=*),parameter     :: subname='[empty_mpi_global_comm]:'
   INTEGER*4                      :: hostnm
      stat = hostnm( this%hostname )
      if(stat .ne. 0) error stop '[!]'//subname//": Unable to obtain hostname"
      if(.not. mpi_is_initialized()) then
         print '("[w]'//subname//': MPI_COMM not initialized")'
         print '("[w]'//subname//': REINITIALIZING MPI_COMM ...")'
         CALL MPI_INIT(err)
         if(err.ne.0) error stop '[!]'//subname//': Unable to initialize MPI_COMM: '
         print '("[i]'//subname//': MPI_COMM initialized OK")'
      endif
      err = init_mpi_global_comm(this)
   end function empty_mpi_global_comm
   !!!!              SUBCOMM
   function local_mpi_comm() result(this)
   type(mpi_sub_comm)             :: this
   !
   integer                        :: info, stat
   character(len=*),parameter     :: subname='[local_mpi_comm]:'
   INTEGER*4                      :: hostnm
      stat = hostnm( this%hostname )
      if(stat .ne. 0) error stop '[!]'//subname//": Unable to obtain hostname"
      if(.not. mpi_is_initialized()) error stop '[!]'//subname//": MPI is not initialized"
      call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, this%comm, info)
      call reconf_comm_rank(this%comm, this%rank, subname)
      call reconf_comm_size(this%comm, this%comm_size, subname)
      this%is_initialized = .true.
   end function local_mpi_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! TODO: not sure how this is compiled with optional argument in the middle
   function mpi_sub_comm_from_split(color, parent_comm, key) result(this)
       integer, intent(IN)                    :: color, key
       type(MPI_COMM), intent(IN), optional   :: parent_comm
       type(mpi_sub_comm)                     :: this
       type(MPI_COMM)                         :: parent
       integer                                :: info, stat
       character(len=*),parameter             :: subname='[!][instantiate_mpi_sub_comm]:'
       INTEGER*4                              :: hostnm
       stat = hostnm( this%hostname )
       if(stat .ne. 0) then; write(*,*) subname, "[ERR]: Unable to obtain hostname"; stop;endif
       if(.not. mpi_is_initialized()) then; write(*,*) subname, "[ERR]: MPI is  not Initialized"; stop;endif
       if(present(parent_comm)) then
           parent = parent_comm
       else
           parent = MPI_COMM_WORLD
       endif
       call MPI_COMM_SPLIT(parent, color, key, this%comm, info)
       if(info.ne.0)then;write(*,*) subname,": Unable to SPLIT MPI_COMM(",parent,"using color", color, ": info=",info;stop;endif
       call reconf_comm_rank(this%comm, this%rank, subname)
       call reconf_comm_size(this%comm, this%comm_size, subname)
       this%is_initialized = .true.
   end function mpi_sub_comm_from_split
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_sub_comm_from_global_comm_split(color, parent_comm, key_) result(this)
   type(mpi_sub_comm)                 :: this
   integer,                intent(IN) :: color
   class(mpi_global_comm), intent(IN) :: parent_comm
   integer,      optional, intent(IN) :: key_
   !
   integer :: info, key, stat
   character(*),parameter :: subname='[mpi_sub_comm_from_global_comm_split]'
   INTEGER*4                          :: hostnm
      stat = hostnm( this%hostname )
      if(stat.ne.0) error stop '[!]'//subname//": Unable to obtain hostname"
      if(.not. mpi_is_initialized()) error stop '[!]'//subname//": MPI is not initialized"
      if(.not. parent_comm%is_initialized) error stop '[!]'//subname//": parent comm is not initialized"
      key = parent_comm%rank; if(present(key_)) key=key_
      call MPI_COMM_SPLIT(parent_comm%comm, color, key, this%comm, info)
      if(info.ne.0)then;write(*,*) subname,": Unable to SPLIT MPI_COMM(",parent_comm%comm,"using color", color, ": info=",info;stop;endif
      call reconf_comm_rank(this%comm, this%rank, subname)
      call reconf_comm_size(this%comm, this%comm_size, subname)
      this%color = color
      this%is_initialized = .true.
   end function mpi_sub_comm_from_global_comm_split
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_sub_comm_from_sub_comm_split(color, parent_comm, key) result(this)
       integer, intent(IN)                           :: color
       Class(mpi_sub_comm), intent(IN)               :: parent_comm
       integer, intent(IN), optional                 :: key
       type(mpi_sub_comm)                            :: this
       integer                                       :: info, key_, stat
       character(len=*),parameter                    :: subname='[!][mpi_sub_comm_from_sub_comm_split]:'
       INTEGER*4                                     :: hostnm
       stat = hostnm( this%hostname )
       if(stat .ne. 0) then; write(*,*) subname, "[ERR]: Unable to obtain hostname"; stop;endif

       if(.not. mpi_is_initialized()) then; write(*,*) subname, "[MPI][ERR] MPI is not initialized"; stop; endif
       if( present(key)) then
         key_ = key
       else
         key_ = parent_comm%rank
       end if
       call MPI_COMM_SPLIT(parent_comm%comm, color, key_, this%comm, info)
       if(info.ne.0)then;write(*,*) subname,": Unable to SPLIT MPI_COMM(",parent_comm%comm,"using color", color, ": info=",info;stop;endif
       call reconf_comm_rank(this%comm, this%rank, subname)
       call reconf_comm_size(this%comm, this%comm_size, subname)
       this%color = color
       this%is_initialized = .true.
   end function mpi_sub_comm_from_sub_comm_split
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  [DESTRUCTORS]
!
   subroutine close_mpi_global_comm(this)
       type(mpi_global_comm), intent(inout)  :: this
       integer                               :: err
       character(len=*),parameter            :: subname='[!][close_mpi_global_comm]:'
       if (this%is_initialized) then
         call MPI_COMM_FREE(this%comm, err)
         !call MPI_FINALIZE(err)
         if(err.ne.0)then;write(*,*)subname,': Problem finalizing MPI COMM: ',err;stop;endif
         this%comm_size             = -2  ! size of communicator
         this%is_initialized        = .false.
         write(*,*)subname,':[GLOBAL_COMM] FINALIZED'
       else
         write(*,*)subname,':[GLOBAL_COMM] ALREADY FINALIZED'
       end if
       !call MPI_FINALIZE(err)
   end subroutine close_mpi_global_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine close_mpi_sub_comm(this)
       type(mpi_sub_comm), intent(inout)     :: this
       integer                               :: err
       character(len=*),parameter            :: subname='[!][close_mpi_sub_comm]:'
       if (this%is_initialized) then
         call MPI_COMM_FREE(this%comm, err)
         if(err.ne.0)then;write(*,*)subname,': Problem finalizing MPI COMM: ',err;stop;endif
         this%comm_size             = -2  ! size of communicator
         this%is_initialized        = .false.
         write(*,*)subname,':[SUB_COMM] FINALIZED'
       else
         write(*,*)subname,':[SUB_COMM] ALREADY FINALIZED'
       end if
   end subroutine close_mpi_sub_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine finalize_mpi_global_comm(this)
       Class(mpi_global_comm), intent(inout) :: this
       integer                               :: err
       character(len=*),parameter            :: subname='[!][finalize_mpi_global_comm]:'
       if (this%is_initialized) then
         call MPI_COMM_FREE(this%comm, err)
         !call MPI_FINALIZE(err)
         if(err.ne.0)then;write(*,*)subname,': Problem finalizing MPI COMM: ',err;stop;endif
         this%comm_size             = -2  ! size of communicator
         this%is_initialized        = .false.
       !  write(*,*)subname,':[GLOBAL_COMM] FINALIZED'
       else
         write(*,*)subname,':[GLOBAL_COMM] ALREADY FINALIZED'
       end if
       !call MPI_FINALIZE(err)
   end subroutine finalize_mpi_global_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine finalize_mpi_sub_comm(this)
       Class(mpi_sub_comm), intent(inout)    :: this
       integer                               :: err
       character(len=*),parameter            :: subname='[!][finalize_mpi_sub_comm]:'
       if (this%is_initialized) then
         call MPI_COMM_FREE(this%comm, err)
         if(err.ne.0)then;write(*,*)subname,': Problem finalizing MPI COMM: ',err;stop;endif
         this%comm_size             = -2  ! size of communicator
         this%is_initialized        = .false.
       !  write(*,*)subname,':[SUB_COMM] FINALIZED'
       else
       !  write(*,*)subname,':[SUB_COMM] ALREADY FINALIZED'
       end if
   end subroutine finalize_mpi_sub_comm

!
!  [INIT]
!
   function mpi_is_initialized() result (is_initialized)
       logical                                :: is_initialized
       call MPI_Initialized(is_initialized)
   end function mpi_is_initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function init_mpi_global_comm(this, VRBZ_) result (err)
   class(mpi_global_comm), intent(INOUT) :: this
   integer, optional,         intent(IN) :: VRBZ_
   !
   integer :: err, info, VRBZ
   character(*),parameter :: subname='[init_mpi_global_comm]:'
      VRBZ=0; if(present(VRBZ_)) VRBZ =VRBZ_
      if (.not. mpi_is_initialized()) error stop "[!]"//subname//": MPI_COMM is not initialized"
      if (this%is_initialized) then
         print '("[w]'//subname//': MPI-GLOBAL-COMM CLASS already initialized")'
      else
         ! Global communicator
         call MPI_COMM_DUP(MPI_COMM_WORLD, this%comm, info)
         if (info.ne.0) error stop '[!]'//subname//': Unable to duplicate global mpi'
         ! Get Global rank
         call MPI_COMM_RANK(MPI_COMM_WORLD, this%rank, info)
         if (info.ne.0) error stop '[!]'//subname//': Unable to determine global mpi rank'
         ! Get Size of Global communicator
         call MPI_COMM_SIZE(MPI_COMM_WORLD, this%comm_size, info)
         if (info.ne.0) error stop '[!]'//subname//': Unable to determine mpi comm size'
         this%is_initialized = .true.
         print '("[i]'//subname//': MPI-GLOBAL-COMM initialized successfully")'
      end if
   end function init_mpi_global_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function init_mpi_sub_comm(this) result (err)
   class(mpi_sub_comm), intent(inout) :: this
   !
   integer :: err, info
   character(*),parameter :: subname='[init_mpi_sub_comm]:'
      if (.not. mpi_is_initialized()) error stop '[!]'//subname//": MPI is not initialized"
      if (this%is_initialized) then
         print '("[w]'//subname//'MPI sub-communicator is already initialized")'
      else ! Comm check (this%comm .ne. -3)
         call reconf_comm_rank(this%comm, this%rank, subname)
         call reconf_comm_size(this%comm, this%comm_size, subname)
         this%is_initialized = .true.
         print '("[i]'//subname//': MPI subcomm initialized successfully")'
      end if ! Comm check
   end function init_mpi_sub_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![           CUDA         ]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  [RECONF]
!
    subroutine reconf_comm_rank(comm, rank, subname, VRBZ)
        type(MPI_COMM), intent(IN)              :: comm
        integer, intent(INOUT)                  :: rank
        character(len=*), intent(IN), optional  :: subname
        logical, intent(IN), optional           :: VRBZ
        integer                                 :: info
        call MPI_COMM_RANK(comm, rank, info)
        if(info.ne.0) then
          if(present(VRBZ)) then
            if(present(subname)) then
                write(*,*)subname,": [ERR] Unable to determine MPI SUB-COMM",comm,"'s rank info:", info; stop
            else
                write(*,*) "[ERR] Unable to determine MPI SUB-COMM",comm,"'s rank info:", info; stop
            end if
          !else
          !  write(*,*)"[ERR] Unable to determine MPI SUB-COMM",comm,"'s rank info:", info; stop
          end if
        else
          if(present(VRBZ)) then
            if(present(subname)) then
                write(*,*)subname,": [WARN] MPI SUB-COMM",comm,"'s rank info UPDATED"
            else
                write(*,*)": [WARN] MPI SUB-COMM",comm,"'s rank info UPDATED"
            end if
          !else
          !  write(*,*)"[WARN] MPI SUB-COMM",comm,"'s rank info UPDATED"
          end if
        endif
    end subroutine reconf_comm_rank
    subroutine reconf_comm_size(comm, comm_size, subname, VRBZ)
        type(MPI_COMM), intent(IN)              :: comm
        integer, intent(INOUT)                  :: comm_size
        character(len=*), intent(IN), optional  :: subname
        logical, intent(IN), optional           :: VRBZ
        integer                                 :: info
        !character(len=*)                        :: subname_
        call MPI_COMM_SIZE(comm, comm_size, info)
        if(info.ne.0) then
            if(present(VRBZ)) then
                if(present(subname)) then
                    write(*,*) subname,": [ERR] Unable to determine MPI SUB-COMM",comm,"'s size info:", info; stop
                else
                    write(*,*) ": [ERR] Unable to determine MPI SUB-COMM",comm,"'s size info:", info; stop
                endif
            !else
            !     write(*,*) "[ERR] Unable to determine MPI SUB-COMM",comm,"'s size info:", info; stop
            end if
        else
            if(present(VRBZ)) then
              if(present(subname)) then
                write(*,*) subname,": [WARN] MPI SUB-COMM",comm,"'s size info UPDATED"
              else
                write(*,*) subname,": [WARN] MPI SUB-COMM",comm,"'s size info UPDATED"
              endif
            !else
            !    write(*,*) "[WARN] MPI SUB-COMM",comm,"'s size info UPDATED"
            end if
        endif
    end subroutine reconf_comm_size

!
!  [BARRIER/SYNC]
!
   subroutine sync_global_comm(this)
       Class(mpi_global_comm), intent(inout)  :: this
       integer                               :: err
       character(len=*),parameter            :: subname='[!][sync_global_comm]:'
       if(.not. mpi_is_initialized()) then; write(*,*) subname, "[MPI][ERR] MPI is not initialized"; stop; endif
       call MPI_BARRIER( this%comm, err)
       if(err.ne.0)then;write(*,*) subname,": Unable to SYNCHRONIZE GLOBAL COMMUNICATOR: err=",err;stop;endif
   end subroutine sync_global_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine sync_sub_comm(this)
       Class(mpi_sub_comm), intent(inout)  :: this
       integer                               :: err
       character(len=*),parameter            :: subname='[!][sync_sub_comm]:'
       if(.not. mpi_is_initialized()) then; write(*,*) subname, "[MPI][ERR] MPI is not initialized"; stop; endif
       call MPI_BARRIER( this%comm, err)
       if(err.ne.0)then;write(*,*) subname,": Unable to SYNCHRONIZE SUB-COMMUNICATOR (", this%comm,"): err=",err;stop;endif
   end subroutine sync_sub_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!m

!
! SPLIT
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function split_global_comm(this, color, key) result(new_comm)
       Class(mpi_global_comm), intent(IN)            :: this
       integer, intent(IN)                           :: color
       type(mpi_sub_comm)                            :: new_comm
       integer, intent(IN), optional                 :: key
       integer                                       :: info, key_
       character(len=*),parameter                    :: subname='[!][split_global_comm]:'
       if(.not. mpi_is_initialized()) then; write(*,*) subname, "[MPI][ERR] MPI is not initialized"; stop; endif
       if( present(key)) then
         key_ = key
       else
         key_ = this%rank
       end if
       call MPI_COMM_SPLIT(this%comm, color, key_, new_comm%comm, info)
       if(info.ne.0)then;write(*,*) subname,": Unable to SPLIT MPI_COMM(",this%comm,"using color", color,", with key = ", key_, ": info=",info;stop;endif
       call reconf_comm_size(new_comm%comm, new_comm%comm_size, subname)
       new_comm%is_initialized = .true.
   end function split_global_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function split_sub_comm(this, color, key) result(new_comm)
       Class(mpi_sub_comm), intent(IN)            :: this
       integer, intent(IN)                           :: color
       type(mpi_sub_comm)                            :: new_comm
       integer, intent(IN), optional                 :: key
       integer                                       :: info, key_
       character(len=*),parameter                    :: subname='[!][split_sub_comm]:'
       if(.not. mpi_is_initialized()) then; write(*,*) subname, "[MPI][ERR] MPI is not initialized"; stop; endif
       if( present(key)) then
         key_ = key
       else
         key_ = this%rank
       end if
       call MPI_COMM_SPLIT(this%comm, color, key_, new_comm%comm, info)
       if(info.ne.0)then;write(*,*) subname,": Unable to SPLIT MPI_COMM(",this%comm,"using color", color,", with key = ", key_, ": info=",info;stop;endif
       call reconf_comm_size(new_comm%comm, new_comm%comm_size, subname)
       new_comm%is_initialized = .true.
   end function split_sub_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  [SEND]
!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Send   Scalars [0D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_send_scal_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_scal_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_scal_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_scal_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_scal_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_scal_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_scal_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_scal_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_scal_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_scal_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Send   Vectors [1D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_send_vect_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vect_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vect_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vect_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vect_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vect_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vect_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vect_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vect_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vect_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Send   Matrices [2D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_send_mat_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_mat_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_mat_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_mat_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_mat_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_mat_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_mat_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_mat_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_mat_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_mat_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Send   Volumes [3D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_send_vol_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vol_i2(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vol_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vol_i4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vol_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vol_i8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vol_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vol_r4(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_send_vol_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_send_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_send_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_send_vol_r8(this, a, dest, tag, err)
        USE comm_send, ONLY      : send
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: dest
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_send_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call send(a, dest, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_send_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  [RECV]
!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Recv   Scalars [0D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_recv_scal_i2(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_scal_i2(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_scal_i4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_scal_i4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_scal_i8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_scal_i8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_scal_r4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_scal_r4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_scal_r8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_scal_r8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Recv   Vectors [1D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_recv_vect_i2(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vect_i2(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vect_i4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vect_i4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vect_i8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vect_i8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vect_r4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vect_r4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vect_r8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vect_r8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Recv   Matrices [2D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_recv_mat_i2(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_mat_i2(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_mat_i4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_mat_i4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_mat_i8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_mat_i8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_mat_r4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_mat_r4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_mat_r8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_mat_r8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ Recv   Volumes [3D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_recv_vol_i2(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vol_i2(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vol_i4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vol_i4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vol_i8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vol_i8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vol_r4(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vol_r4(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_recv_vol_r8(this, a, root, tag, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_recv_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_recv_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_recv_vol_r8(this, a, root, tag, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(IN),    optional     :: tag
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_recv_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          tag_ = 123456789; if(present(tag)) tag_ = tag
          call RECV_(a, root, this%comm, tag_, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_recv_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  [BCAST]
!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ bcast   Scalars [0D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_bcast_scal_i2(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_scal_i2(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_scal_i4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_scal_i4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_scal_i8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_scal_i8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_scal_r4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_scal_r4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_scal_r8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_scal_r8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ bcast   Vectors [1D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_bcast_vect_i2(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vect_i2(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vect_i4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vect_i4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vect_i8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vect_i8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vect_r4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vect_r4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vect_r8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vect_r8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ bcast   Matrices [2D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_bcast_mat_i2(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_mat_i2(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_mat_i4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_mat_i4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_mat_i8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_mat_i8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_mat_r4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_mat_r4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_mat_r8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_mat_r8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! [ bcast   Volumes [3D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_global_comm_bcast_vol_i2(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vol_i2(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vol_i4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vol_i4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vol_i8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vol_i8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vol_r4(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vol_r4(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_global_comm_bcast_vol_r8(this, a, root, err)
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_global_comm_bcast_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_global_comm_bcast_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpi_sub_comm_bcast_vol_r8(this, a, root, err)
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)        :: a(:,:,:)
        integer, intent(IN)                  :: root
        integer, intent(INOUT), optional     :: err
        integer                              :: err_, tag_
        character(len=*),parameter           :: subname='[!][mpi_sub_comm_bcast_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
              call BCAST_(a, root, this%comm, err_)
          if(present(err)) err = err_
  end  subroutine mpi_sub_comm_bcast_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  [ALLREDUCE]
!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! [ aredux   Scalars [0D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_global_comm_allreduce_scal_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a
        integer(2)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_scal_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a
        integer(2)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_scal_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_scal_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_scal_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a
        integer(4)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_scal_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a
        integer(4)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_scal_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_scal_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_scal_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a
        integer(8)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_scal_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a
        integer(8)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_scal_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_scal_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_scal_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a
        real(4)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_scal_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a
        real(4)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_scal_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_scal_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_scal_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a
        real(8)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_scal_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a
        real(8)                      :: res
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_scal_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_scal_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! [ aredux   Vectors [1D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_global_comm_allreduce_vect_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a(:)
        integer(2)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vect_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a(:)
        integer(2)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vect_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vect_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vect_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a(:)
        integer(4)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vect_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a(:)
        integer(4)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vect_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vect_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vect_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a(:)
        integer(8)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vect_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a(:)
        integer(8)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vect_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vect_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vect_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a(:)
        real(4)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vect_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a(:)
        real(4)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vect_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vect_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vect_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a(:)
        real(8)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vect_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a(:)
        real(8)                      :: res(size(a,1))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vect_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vect_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! [ aredux   Matrices [2D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_global_comm_allreduce_mat_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a(:,:)
        integer(2)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_mat_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a(:,:)
        integer(2)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_mat_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_mat_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_mat_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a(:,:)
        integer(4)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_mat_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a(:,:)
        integer(4)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_mat_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_mat_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_mat_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a(:,:)
        integer(8)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_mat_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a(:,:)
        integer(8)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_mat_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_mat_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_mat_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a(:,:)
        real(4)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_mat_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a(:,:)
        real(4)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_mat_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_mat_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_mat_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a(:,:)
        real(8)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_mat_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a(:,:)
        real(8)                      :: res(size(a,1),size(a,2))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_mat_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_mat_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!! [ aredux   Volumes [3D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_global_comm_allreduce_vol_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a(:,:,:)
        integer(2)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vol_i2(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(2), intent(IN)          :: a(:,:,:)
        integer(2)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vol_i2]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vol_i2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vol_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a(:,:,:)
        integer(4)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vol_i4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(4), intent(IN)          :: a(:,:,:)
        integer(4)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vol_i4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vol_i4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vol_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a(:,:,:)
        integer(8)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vol_i8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        integer(8), intent(IN)          :: a(:,:,:)
        integer(8)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vol_i8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vol_i8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vol_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a(:,:,:)
        real(4)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vol_r4(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(4), intent(IN)          :: a(:,:,:)
        real(4)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vol_r4]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vol_r4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_global_comm_allreduce_vol_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_global_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a(:,:,:)
        real(8)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_global_comm_allreduce_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_global_comm_allreduce_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mpi_sub_comm_allreduce_vol_r8(this, a, oper) result (res)
        USE comm_types, Only : get_mpi_op
        Class(mpi_sub_comm), intent(in)  :: this
        !type(cls), intent(in)  :: this
        real(8), intent(IN)          :: a(:,:,:)
        real(8)                      :: res(size(a,1),size(a,2),size(a,3))
        character(len=*), intent(IN), optional  :: oper
        character(len=*),parameter             :: subname='[!][mpi_sub_comm_allreduce_vol_r8]:'
          if (.not. this%is_initialized) then; write(*,*) subname, "MPI communicator Not Initialized"; stop; end if
          if(present(oper)) then
            res = ALLREDUCE_(a, oper, this%comm)
          else
            res = ALLREDUCE_(a, "sum", this%comm)
          endif
  end  function mpi_sub_comm_allreduce_vol_r8
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module communicators
