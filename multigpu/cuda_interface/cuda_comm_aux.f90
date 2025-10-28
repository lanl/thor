!!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/
!!
!! @file cuda_comm_aux.f90
!! @author Ismael Djibrilla Boureima, Oleg Korobkin
!! @date  November 2023
!! @brief Fortran interface to C wrappers for NCCL and CAL (see cuda_cal_aux.c)
!!
module cuda_comm_aux
  use iso_c_binding
  use nccl

  !> ---------------------------------------------------------------
  !! Public module globals
  type(c_ptr)           :: cal_handle     !< opaque pointer to CAL communicator
  type(c_ptr)           :: nccl_handle    !< opaque pointer to NCCL communicator
  type(c_ptr)           :: local_stream   !< opaque poitner to CUDA stream
  type(ncclUniqueId)    :: global_nccl_id !< custom Fortran type to NCCL ID
  type(ncclComm)        :: nccl_comm      !< custom Fortran type to NCCL communicator

  !> ----------------------------------------------------------------
  !! Private globals of the module
  !!
  !! The following two module variables are made private to prohibit 
  !! external modification. 
  !! To inquire the status of CAL and NCCL initialization, use the public 
  !! functions "is_cal_initialized" and "is_nccl_initialized", resp.  
  logical, private      :: cal_initialized = .false.
  logical, private      :: nccl_initialized = .false.

  !! Prohibit calling C interface function directly
  !private :: set_localStream
  private :: init_cal_comm_c,     &
             cal_stream_sync_c,   &
             cal_get_comm_size_c, &
             cal_get_rank_c,      &
             show_cal_comm_c,     &
             destroy_cal_comm_c,  &
             get_hw_topology_c,   &
             set_localstream_c,   &
             show_localstream_c,  &
             calaux_stop_if_error

  interface
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int set_localStream(cudaStream_t *localStream, int VRBZ)
   integer(c_int) function set_localStream_c(localStream, VRBZ) bind(C,name="set_localStream")
     use iso_c_binding
     implicit none
     type(c_ptr)           :: localStream
     integer(c_int), value :: VRBZ
   end function

   ! C-PROTOTYPE: show_localStream(cudaStream_t *localStream);
   integer(c_int) function show_localStream_c(localStream) bind(C,name="show_localStream")
     use iso_c_binding
     implicit none
     type(c_ptr) :: localStream
   end function

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int init_cal_comm(cal_comm_t *cal_comm, int *rank, int *commSize, int *localRank, int VRBZ)
   integer(c_int) function init_cal_comm_c(cal_comm, rank, commSize, localRank, VRBZ) bind(C,name="init_cal_comm")
     use iso_c_binding
     implicit none
     type(c_ptr)           :: cal_comm
     type(c_ptr), value    :: rank, commSize, localRank
     integer(c_int), value :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-Prototype: calError_t cal_stream_sync(cal_comm_t comm, cudaStream_t stream);
   integer(c_int) function cal_stream_sync_c(comm, stream) bind(C,name="cal_stream_sync")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: comm
     type(c_ptr),value :: stream
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-Prototype: calError_t cal_get_comm_size( cal_comm_t comm, int* size );
   integer(c_int) function cal_get_comm_size_c(comm, siz) bind(C,name="cal_get_comm_size")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: comm
     type(c_ptr),value :: siz
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-Prototype: calError_t cal_get_rank( cal_comm_t comm, int* rank );
   integer(c_int) function cal_get_rank_c(comm, rank) bind(C,name="cal_get_rank")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: comm
     type(c_ptr),value :: rank
   end function

   ! C-Prototype: int show_cal_comm(cal_comm_t *cal_comm);
   integer(c_int) function show_cal_comm_c(cal_comm) bind(C,name="show_cal_comm")
     use iso_c_binding
     implicit none
     type(c_ptr)        :: cal_comm
   end function

   ! C-PROTOTYPE: int destroy_cal_comm(cal_comm_t *cal_comm,  int VRBZ)
   integer(c_int) function destroy_cal_comm_c(cal_comm, VRBZ) bind(C,name="destroy_cal_comm")
     use iso_c_binding
     implicit none
     type(c_ptr)           :: cal_comm
     integer(c_int), value :: VRBZ
   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C-PROTOTYPE: int get_hw_topology(int *rank, int *commSize, MPI_Comm *localComm, int *localRank, int *localCommSize,
   !             MPI_Comm *ctxComm, int *ctxRank, int *ctxCommSize, int *nDevLocal, int ctx_color, int *color, int VRBZ);
   integer(c_int) function get_hw_topology_c(rank, commSize, localComm, localRank, localCommSize, &
                           ctxComm, ctxRank, ctxCommSize, nDevLocal, ctx_color, color, VRBZ) bind(C,name="get_hw_topo")
     use iso_c_binding
     implicit none
     type(c_ptr), value  :: rank, commSize, localComm, localRank, localCommSize, &
                            ctxComm,ctxRank, ctxCommSize, nDevLocal, color
     integer(c_int), value :: ctx_color, VRBZ
   end function
  end interface

contains

  subroutine calaux_stop_if_error(err, msg)
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: error_unit
  integer(c_int), intent(in) :: err
  character(*), intent(in)   :: msg
     if (err.ne.0) then
        write(error_unit,'(A,I8)') msg, err
        error stop
     endif
  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Wrappers with native Fortran types
  !!!!!!
  subroutine get_hw_topology(ctx_color, rank, comm_Size, localComm, &
                   localRank, localCommSize, ctxComm, ctxRank, ctxCommSize, &
                   nDevLocal, color, VRBZ_)
  use iso_c_binding
  use mpi_f08
  implicit none
  integer, intent(IN) :: ctx_color
  integer, intent(INOUT), target:: rank, comm_Size, localRank, localCommSize, &
                                   ctxRank, ctxCommSize, nDevLocal, color
  type(MPI_COMM), intent(INOUT),target :: localComm, ctxComm
  integer, optional   :: VRBZ_
  !
  integer(c_int) :: VRBZ, err

     err= get_hw_topology_c(C_LOC(rank), C_LOC(comm_Size), C_LOC(localComm), &
                            C_LOC(localRank), C_LOC(localCommSize), &
                            C_LOC(ctxComm), C_LOC(ctxRank), C_LOC(ctxCommSize), &
                            C_LOC(nDevLocal), int(ctx_color,c_int), C_LOC(color), VRBZ)
     call calaux_stop_if_error(err, "[!][get_hw_topology], err=")

  end subroutine get_hw_topology


  subroutine cuda_comm_init(rank, commSize, localRank, VRBZ_)
  use iso_c_binding
  use mpi_f08, only: MPI_COMM, MPI_BCAST, MPI_INT, MPI_BYTE, MPI_COMM_WORLD
  use nccl
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  integer, intent(INOUT), target:: rank, commSize, localRank
  integer, optional   :: VRBZ_
  !
  character(*), parameter :: subnam = '[cal_comm_init]'
  integer(c_int)          :: VRBZ, err
  type(ncclResult)        :: nccl_err
  type(MPI_COMM)          :: my_mpi_comm
  integer                 :: ierr

     VRBZ=0; if (present(VRBZ_)) VRBZ=int(VRBZ_,c_int)
     if (cal_initialized) then
        if (VRBZ > 1) write (error_unit, '("WARNING: CAL already initialized")')
     else
        err= init_cal_comm_c(cal_handle, C_LOC(rank), C_LOC(commSize), C_LOC(localRank), VRBZ)
        call calaux_stop_if_error(err, "[!]"//subnam//", init_cal_comm_c, err=")

        err= set_localStream_c(local_stream, VRBZ)
        call calaux_stop_if_error(err, "[!]"//subnam//", set_localStream_c, err=")

        cal_initialized= .true.
     endif

     if (nccl_initialized) then
        if (VRBZ > 1) write (error_unit, '("WARNING: NCCL already initialized")')
     else 
        if (rank.eq.0) then
           nccl_err = ncclGetUniqueId(global_nccl_id)
           call calaux_stop_if_error(nccl_err% member, &
                                    "[!]"//subnam//", ncclGetUniqueId, err=")
        endif
        call MPI_Bcast(global_nccl_id, int(sizeof(global_nccl_id)), &
                       MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
        call calaux_stop_if_error(int(err,c_int), "[!]"//subnam//", MPI_Bcast, err=")

        nccl_err= ncclCommInitRank(nccl_comm, commSize, global_nccl_id, rank);
        call calaux_stop_if_error(nccl_err% member, &
                                "[!]"//subnam//", ncclCommInitRank, err=")
        nccl_handle = nccl_comm% member
        nccl_initialized= .true.
     endif 

  end subroutine cuda_comm_init


  logical function is_cal_initialized()
  is_cal_initialized= cal_initialized
  end function


  logical function is_nccl_initialized()
  is_nccl_initialized= nccl_initialized
  end function


  subroutine show_cal_comm()
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: error_unit
  !
  integer(c_int) :: err

     if (.not.cal_initialized) then
        write (error_unit, '("[!][show_cal_comm] ERROR: CAL is not initialized")')
        return
     endif
     err= show_cal_comm_c(cal_handle)
     call calaux_stop_if_error(err, "[!][show_cal_comm], err=")

  end subroutine show_cal_comm


  subroutine cal_comm_destroy(VRBZ_)
  use iso_c_binding
  implicit none
  integer, optional   :: VRBZ_
  !
  integer(c_int) :: VRBZ, err

     VRBZ=0; if (present(VRBZ_)) VRBZ=int(VRBZ_,c_int)
     err= destroy_cal_comm_c(cal_handle, VRBZ)
     call calaux_stop_if_error(err, "[!][destroy_cal_comm_f], err=")
     cal_initialized= .false.

  end subroutine cal_comm_destroy

end module cuda_comm_aux

