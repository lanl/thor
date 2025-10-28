module distributed_comms
use mpi_f08, only: MPI_COMM, &
                   MPI_Initialized, &
                   MPI_COMM_SPLIT, &
                   MPI_COMM_SPLIT_TYPE, &
                   MPI_COMM_WORLD, &
                   MPI_COMM_TYPE_SHARED,&
                   MPI_INFO_NULL, &
                   MPI_FLOAT, &
                   MPI_INT, &
                   MPI_DOUBLE, &
                   MPI_BYTE
use communicators, only: mpi_global_comm, mpi_sub_comm
use cudafor
use nccl
use cublasmp_aux
use cusolvermp_aux
use string_lib, only: str
implicit none
!include 'mpif.h'
type,public              :: super_comm                  ! modern mpi communicator
   type(mpi_global_comm) :: comm                        ! Global MPI communicator
   type(mpi_sub_comm)    :: local_comm                  ! Local MPI communicator
   type(mpi_sub_comm)    :: masters_comm                ! Local Masters communicator
   type(mpi_sub_comm)    :: ctx_comm                    ! Inter-node CTX communicator
   type(mpi_sub_comm)    :: ctx_local_comm              ! Intra-node CTX communicator
   integer               :: host_number           = -1  ! Node number (host is a reserved fortran name)
   integer               :: rank                  = -1  ! mpi rank
   integer               :: comm_size             = -2  ! size of communicator
   integer               :: local_rank            = -4  ! local mpi rank
   integer               :: local_comm_size       = -5  ! size of local communicator
   integer               :: cuda_device           = -7  ! Cuda device
   integer               :: n_local_cuda_device   = -8  ! Local number of cuda GPUs
   integer               :: n_global_cuda_device  = -9  ! Global number of cuda GPUs
   integer               :: ctx_rank              = -12 ! rank on inter-nodr ctx communicator
   integer               :: ctx_local_rank        = -13 ! rank on intra-nodr ctx communicator
   integer               :: ctx_size              = -14 ! rank on inter-nodr ctx communicator
   integer               :: ctx_local_size        = -15 ! rank on intra-nodr ctx communicator
   logical               :: in_gpu_ctx            = .false.
   logical               :: GMASTER               = .false.
   logical               :: LMASTER               = .false.
   logical               :: mpi_comm_init         = .false.
   logical               :: cuda_backend_init     = .false.
   type(c_ptr)           :: cal_comm
   type(MPI_COMM)        :: cal_local_comm, cal_ctx_comm
   logical               :: local_cal_comm_init   = .false.
   logical               :: global_cal_comm_init  = .false.
   type(c_ptr)           :: nccl_handle
   logical               :: nccl_comm_init        = .false.
   logical               :: local_nccl_comm_init  = .false.
   type(ncclUniqueId)    :: nccl_id
   type(ncclUniqueId)    :: local_nccl_id
   type(ncclComm)        :: nccl_comm
   type(ncclComm)        :: local_nccl_comm
   type(c_ptr)           :: localStream
   type(c_ptr)           :: cuBlasMP_handle, cuSolverMP_handle
 contains
   procedure :: init_mpi         => init_mpi_comm
   procedure :: init_cuda        => init_cuda_backend
   procedure :: init_nccl_comm   => build_global_nccl_comm
   procedure :: get_rank         => get_comm_rank
   procedure :: get_size         => get_comm_size
   procedure :: barrier          => comm_barrier
   procedure :: finalize         => finalize_comm
   procedure :: get_mpi_op       => get_mpi_op
   procedure :: translate_mpi_op => get_mpi_op
   procedure :: reset_cuda_handles
   !procedure       :: send             => global_comm_send
   !procedure       :: recv             => global_comm_recv
   !procedure       :: bcast            => global_comm_bcast
   !procedure       :: all_reduce       => global_comm_all_reduce
   !final           :: dealloc_comm
end type super_comm

!
!  [MODULE INTERFACES]
!
 !! CONSTRUCTOR\DESTRUCTOR
 interface super_comm ! distributed comm constructor interface
    procedure :: empty_super_comm
 end interface super_comm
 !! GLOBAL SEND
 interface global_comm_send
     module procedure global_comm_send_i, global_comm_send_i1, global_comm_send_i2, global_comm_send_i3,&
                      global_comm_send_f, global_comm_send_f1, global_comm_send_f2, global_comm_send_f3,&
                      global_comm_send_d, global_comm_send_d1, global_comm_send_d2, global_comm_send_d3
 end interface global_comm_send
 !! GLOBAL RECV
 interface global_comm_recv
     module procedure global_comm_recv_i, global_comm_recv_i1, global_comm_recv_i2, global_comm_recv_i3,&
                      global_comm_recv_f, global_comm_recv_f1, global_comm_recv_f2, global_comm_recv_f3,&
                      global_comm_recv_d, global_comm_recv_d1, global_comm_recv_d2, global_comm_recv_d3
 end interface global_comm_recv
 !! GLOBAL BCAST
 interface global_comm_bcast
     module procedure global_comm_bcast_i, global_comm_bcast_i1, global_comm_bcast_i2, global_comm_bcast_i3,&
                      global_comm_bcast_f, global_comm_bcast_f1, global_comm_bcast_f2, global_comm_bcast_f3,&
                      global_comm_bcast_d, global_comm_bcast_d1, global_comm_bcast_d2, global_comm_bcast_d3,&
                      global_comm_bcast_nccl_uid
 end interface global_comm_bcast
!! GLOBAL ALLREDUCE
 interface global_comm_all_reduce
     module procedure global_comm_all_reduce_i, global_comm_all_reduce_i1, global_comm_all_reduce_i2, global_comm_all_reduce_i3,&
                      global_comm_all_reduce_f, global_comm_all_reduce_f1, global_comm_all_reduce_f2, global_comm_all_reduce_f3,&
                      global_comm_all_reduce_d, global_comm_all_reduce_d1, global_comm_all_reduce_d2, global_comm_all_reduce_d3
 end interface global_comm_all_reduce


contains
!
!  [CONSTRUCTORS]
   !!!!              GLOBAL
   function empty_super_comm() result(this)
       type(super_comm)               :: this
       integer                        :: err
       character(len=*),parameter     :: subname='[!][intenciate_super_comm]:'
       !err = init_mpi_comm(this)
   end function empty_super_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  [INIT]
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!![           MPI         ]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function init_mpi_comm(this) result (err)
   class(super_comm), intent(inout) :: this
   !
   integer :: err, n_local_cuda_device, color
   character(*),parameter :: subname='[init_mpi_comm]:'
      if (.not. this%mpi_comm_init) then
         call MPI_INIT(err)
         if(err.ne.0) error stop '[!]'//subname//': Unable to initialize MPI COMM'
         this%comm = mpi_global_comm() !MPI_COMM_WORLD ! Global communicator
         this%rank = this%comm%rank  ! Get Global rank
         this%comm_size = this%comm%comm_size ! Get Size of Global communicator
         this%local_comm = mpi_sub_comm() ! Local communicator
         this%local_rank = this%local_comm%rank  ! Get local rank
         this%local_comm_size = this%local_comm%comm_size ! Get size of local comm
         if (this%local_rank == 0) then
            color  = 101
            this%LMASTER = .true.
         else
            color  = 100
            this%LMASTER = .false.
         end if
         this%masters_comm = mpi_sub_comm(color, this%comm)
         if (this%LMASTER) then
            this%host_number = this%masters_comm%rank
         else
            this%host_number = -1
         end if
         this%mpi_comm_init = .true.
         print '("[i]'//subname//': MPI communicator initialized")'
      else
         print '("[w]'//subname//': MPI communicator already initialized")'
      end if
   end function init_mpi_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!![           CUDA         ]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function init_cuda_backend(this, VRBZ_) result (err)
   use cuda_comm_aux, only  : get_hw_topology, cuda_comm_init, local_stream, &
                              is_cal_initialized, cal_handle, &
                              is_nccl_initialized, nccl_comm, nccl_handle, global_nccl_id
   use cublasmp_aux, only   : cublasmp_init, cublasmp_is_initialized, cublasmp_handle
   use cusolvermp_aux, only : cusolvermp_init, cusolvermp_is_initialized, cusolvermp_handle
   use string_lib, only     : str
   class(super_comm), intent(INOUT) :: this
   integer, intent(IN), optional :: VRBZ_
   !
   integer :: ierr, info, lcolor, ctx_color, VRBZ, localDeviceId
   integer :: n_local_cuda_device, n_global_cuda_device
   character(len=*),parameter :: subnam='[init_cuda_backend]:'
   integer :: rank, commSize, localRank, localCommSize, ctxRank, &
              ctxCommSize,  nDevLocal, color
   type(MPI_COMM) :: localComm, ctxComm

      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
      if (.not.this% mpi_comm_init) then
         print '("[!] MPI COMM is not initiliazed")'
         info = this% init_mpi()
         if(info /= 0) error stop '[!]'//subnam//': in init_mpi(), info='//str(info)
         print '("[+] MPI Initialized OK")'
      endif
      info = cudaGetDeviceCount(n_local_cuda_device)
      if(info /= 0) error stop '[!]'//subnam//': in cudaGetDeviceCount, info='//str(info)
      this% n_local_cuda_device =  n_local_cuda_device
      if (this%local_rank < n_local_cuda_device) then
          info = cudaSetDevice(this% local_rank)
          if(info /= 0) error stop '[!]'//subnam//': in cudaSetDevice, rank=' // &
                                   str(this%local_rank)//', info='//str(info)
          this%cuda_device = this%local_rank
          this%in_gpu_ctx  = .true.
          color            = 11  ! TODO: remove hardcoded values
          lcolor           = 33
      else
          color            = 22
          lcolor           = 66
      end if
      err = this%barrier(); if(err.ne.0) error stop '[!]'//subnam//': in sync [1/4]'
      if (this% lmaster) &
         n_global_cuda_device= this% masters_comm% allreduce(n_local_cuda_device, "sum")

      err = this%barrier(); if(err.ne.0) error stop '[!]'//subnam//': in sync [2/4]'
      call this%local_comm%bcast(n_global_cuda_device, 0, ierr)
      err = this%barrier(); if(err.ne.0) error stop '[!]'//subnam//': in sync [3/4]'
      this% n_global_cuda_device = n_global_cuda_device

      this%ctx_comm = mpi_sub_comm(color, this%comm)
      this%ctx_size = this%ctx_comm%comm_size
      this%ctx_rank = this%ctx_comm%rank

      err= this%barrier(); if(err.ne.0) error stop '[!]'//subnam//': in sync [4/4]'
      this%ctx_local_comm = mpi_sub_comm(lcolor, this%local_comm)
      this%ctx_local_size = this%ctx_local_comm%comm_size
      this%ctx_local_rank = this%ctx_local_comm%rank

      ctx_color=11
      localDeviceId=this% cuda_device
      call get_hw_topology(ctx_color, rank, commSize, localComm, &
                           localRank, localCommSize, ctxComm, ctxRank, ctxCommSize, &
                           nDevLocal, color, VRBZ_=VRBZ)
      call this% comm% barrier()
      this% cal_local_comm = localComm
      this% cal_ctx_comm = ctxComm
      if (this% in_gpu_ctx) then
         if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: ' &
                             //subnam//': in GPU context")'
         call cuda_comm_init(rank, commSize, localRank, VRBZ_=VRBZ)
         this% localStream = local_stream
         this% global_cal_comm_init = is_cal_initialized()
         if (this% global_cal_comm_init) then
            this% cal_comm = cal_handle
            if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: ' &
                             //subnam//': CAL communicator initialized")'
         else
            if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: ' &
                             //subnam//': CAL communicator NOT initialized")'
         endif
         this% nccl_comm_init = is_nccl_initialized()
         if (this% nccl_comm_init) then
            this% nccl_comm = nccl_comm
            this% nccl_handle = nccl_handle
            this% nccl_id = global_nccl_id
            !! TODO: local_nccl_comm?
            if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: ' &
                             //subnam//': NCCL communicator initialized")'
         else
            if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: ' &
                             //subnam//': NCCL communicator NOT initialized")'
         endif
         if (.not.cublasmp_is_initialized) then
            call cublasmp_init(VRBZ_=VRBZ)
            this% cuBlasMP_handle = cublasmp_handle
            if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: ' &
                             //subnam//': cuBLASMp library initialized")'
         endif
         if (.not.cusolvermp_is_initialized) then
            call cusolvermp_init(localDeviceId, VRBZ_=VRBZ)
            this% cuSolverMp_handle = cusolvermp_handle
            if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: ' &
                             //subnam//': cuSolverMp library initialized")'
         endif
      else
         if (VRBZ > 2) print '("[i]['//str(rank)//'/'//str(commSize)//']: '//subnam &
                             //' not in GPU context")'
      endif
      call this% comm% barrier()
      this%cuda_backend_init = .true.
   end function init_cuda_backend
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!![           NCCL         ]!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function build_global_nccl_comm(this, VRBZ_) result (err)
   class(super_comm), intent(INOUT) :: this
   integer, optional                :: VRBZ_
   !
   integer            :: info, VRBZ
   type(ncclUniqueId) :: nccl_id
   type(ncclResult)   :: nccl_result
   character(len=*),parameter :: subname='[!][build_global_nccl_comm]:'

      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
      if (.not. this%cuda_backend_init) then
         print '("[!] CUDA BACKEND is not initilized")'
         info = this% init_cuda()
         if (info /= 0) error stop subname//': in init_cuda(), info='//str(info)
         print '("[+] CUDA BACKEND Initialized OK")'
      endif
      if (this%in_gpu_ctx) then
         if (this%rank .eq. 0) then
            nccl_result = ncclGetUniqueId(this% nccl_id)
         end if
         err = this%barrier()
         if(err.ne.0) error stop subname//': in SYNC, info='//str(info)
         call global_comm_bcast(this, this% nccl_id, root=0)

         if (VRBZ>0) print '(A,I2)',"[+]["//str(this%rank)//"]: nccl_id = ", nccl_id
      end if
      this% nccl_comm_init = .true.
   end function build_global_nccl_comm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  [GETTERS]
!
   function get_comm_rank(this) result (rank)
       Class(super_comm), intent(in) :: this
       integer                     :: rank
       character(len=*),parameter  :: subname='[!][get_comm_rank]:'
       if (this%mpi_comm_init) then
           rank = this%rank
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end function get_comm_rank
   function get_comm_size(this) result (comm_size)
       Class(super_comm), intent(in) :: this
       integer                     :: comm_size
       character(len=*),parameter  :: subname='[!][get_comm_size]:'
       if (this%mpi_comm_init) then
           comm_size = this%comm_size
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end function get_comm_size
   function get_mpi_op(this, oper) result (op)
        Use mpi_f08, Only : MPI_OP, MPI_OP_NULL, MPI_SUM, MPI_MIN, MPI_MAX, MPI_PROD, MPI_LAND, MPI_BAND, &
                            MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MINLOC, MPI_MAXLOC, MPI_REPLACE

        Class(super_comm), intent(in) :: this
        character(len=*), intent(in)  :: oper
        TYPE(MPI_OP)                  :: op
        SELECT CASE (oper)
            CASE ('null')
                op = MPI_OP_NULL
            CASE ('sum')
                op = MPI_SUM
            CASE ('min')
                op = MPI_MIN
            CASE ('max')
                op = MPI_MAX
            CASE ('prod')
                op = MPI_PROD
            CASE ('land')
                op = MPI_LAND
            CASE ('band')
                op = MPI_BAND
            CASE ('lor')
                op = MPI_LOR
            CASE ('bor')
                op = MPI_BOR
            CASE ('lxor')
                op = MPI_LXOR
            CASE ('bxor')
                op = MPI_BXOR
            CASE ('arg_min')
                op = MPI_MINLOC
            CASE ('arg_max')
                op = MPI_MAXLOC
            CASE ('min_loc')
                op = MPI_MINLOC
            CASE ('max_loc')
                op = MPI_MAXLOC
            CASE ('replace')
                op = MPI_REPLACE
            CASE DEFAULT
                op = MPI_SUM
        END SELECT
   end function get_mpi_op

!
!  [SEND]
!
   subroutine global_comm_send_i(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       !type(super_comm), intent(inout)      :: this
       integer, intent(in)              :: a, dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_i]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, 1, MPI_INT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_i
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_f(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       real, intent(in)                 :: a
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_f]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, 1, MPI_FLOAT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_f
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_d(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       double precision, intent(in)     :: a
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_d]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, 1, MPI_DOUBLE, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_d
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Vector
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_i1(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       integer, intent(in)              :: a(:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_i1]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1), MPI_INT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_i1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_f1(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       real, intent(in)                 :: a(:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_f1]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1), MPI_FLOAT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_f1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_d1(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       double precision, intent(in)     :: a(:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_d1]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1), MPI_DOUBLE, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_d1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Mat
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_i2(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       integer, intent(in)              :: a(:,:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_i2]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1)*size(a,2), MPI_INT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_f2(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       real, intent(in)                 :: a(:,:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_f2]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1)*size(a,2), MPI_FLOAT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_f2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_d2(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       double precision, intent(in)     :: a(:,:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_d2]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1)*size(a,2), MPI_DOUBLE, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_d2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Vol
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_i3(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       integer, intent(in)              :: a(:,:,:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_i3]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_INT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_i3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_f3(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       real, intent(in)                 :: a(:,:,:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_f3]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_FLOAT, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_f3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_send_d3(this, a, dest, tag, err)
       Class(super_comm), intent(inout)      :: this
       double precision, intent(in)     :: a(:,:,:)
       integer, intent(in)              :: dest
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_send_d3]:'
       tag_ = 123456789
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_DOUBLE, dest, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_send_d3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  [RECIVE]
!
   subroutine global_comm_recv_i(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       integer, intent(in)              :: a, root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_i]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, 1, MPI_INT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_i
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_f(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       real, intent(in)                 :: a
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_f]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, 1, MPI_FLOAT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_f
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_d(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       double precision, intent(in)     :: a
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_d]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, 1, MPI_DOUBLE, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_d
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Vector
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_i1(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       integer, intent(in)              :: a(:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_i1]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1), MPI_INT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_i1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_f1(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       real,    intent(in)              :: a(:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_f1]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1), MPI_FLOAT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_f1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_d1(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       double precision,    intent(in)  :: a(:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_d1]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1), MPI_DOUBLE, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_d1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Mat
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_i2(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       integer, intent(in)              :: a(:,:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_i2]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1)*size(a,2), MPI_INT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_f2(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       real,    intent(in)              :: a(:,:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_f2]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1)*size(a,2), MPI_FLOAT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_f2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_d2(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       double precision,    intent(in)  :: a(:,:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_d2]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1)*size(a,2), MPI_DOUBLE, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_d2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Vol
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_i3(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       integer, intent(in)              :: a(:,:,:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_i3]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1)*size(a,2)*size(a,3), MPI_INT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_i3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_f3(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       real,    intent(in)              :: a(:,:,:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_f3]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1)*size(a,2)*size(a,3), MPI_FLOAT, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_f3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_recv_d3(this, a, root, tag, err)
       Class(super_comm), intent(in)      :: this
       double precision,    intent(in)  :: a(:,:,:)
       integer, intent(in)              :: root
       integer, intent(in),    optional :: tag
       integer, intent(inout), optional :: err
       integer                          :: err_, tag_
       character(len=*),parameter  :: subname='[!][global_comm_recv_d3]:'
       tag_ = 87654321
       if (this%mpi_comm_init) then
           if(present(tag)) tag_ = tag
           call MPI_RECV(a, size(a,1)*size(a,2)*size(a,3), MPI_DOUBLE, root, tag_, this%comm, err_)
           if(present(err)) err = err_
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end subroutine global_comm_recv_d3
!
!  [BCAST]
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Scalar
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_i(this, a, root)
       Class(super_comm), intent(in)          :: this
       integer, intent(inout)                 :: a
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_i]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, 1, MPI_INT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_i
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_f(this, a, root)
       Class(super_comm), intent(in)          :: this
       real, intent(inout)                    :: a
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_f]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, 1, MPI_FLOAT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_f
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_d(this, a, root)
       Class(super_comm), intent(in)          :: this
       double precision, intent(inout)        :: a
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_d]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, 1, MPI_DOUBLE, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_d
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_nccl_uid(this, a, root)
       Class(super_comm), intent(in)          :: this
       type(ncclUniqueId), intent(inout)      :: a
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_nccl_uid]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, sizeof(a), MPI_BYTE, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_nccl_uid
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Vector
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_i1(this, a, root)
       Class(super_comm), intent(in)          :: this
       integer, intent(inout)                 :: a(:)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_i1]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1), MPI_INT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_i1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_f1(this, a, root)
       Class(super_comm), intent(in)          :: this
       real, intent(inout)                    :: a(:)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_f1]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1), MPI_FLOAT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_f1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_d1(this, a, root)
       Class(super_comm), intent(in)          :: this
       double precision, intent(inout)        :: a(:)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_d1]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1), MPI_DOUBLE, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_d1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Mat
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_i2(this, a, root)
       Class(super_comm), intent(in)          :: this
       integer, intent(inout)                 :: a(:,:)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_i2]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1)*size(a,2), MPI_INT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_f2(this, a, root)
       Class(super_comm), intent(in)          :: this
       real, intent(inout)                    :: a(:,:)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_f2]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1)*size(a,2), MPI_FLOAT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_f2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_d2(this, a, root)
       Class(super_comm), intent(in)          :: this
       double precision, intent(inout)        :: a(:, :)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_d2]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1)*size(a,2), MPI_DOUBLE, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_d2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Vol
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_i3(this, a, root)
       Class(super_comm), intent(in)          :: this
       integer, intent(inout)                 :: a(:,:,:)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_i3]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1)*size(a,2)*size(a,3), MPI_INT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_i3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_f3(this, a, root)
       Class(super_comm), intent(in)          :: this
       real, intent(inout)                    :: a(:, :, :)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_f3]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1)*size(a,2)*size(a,3), MPI_FLOAT, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_f3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine global_comm_bcast_d3(this, a, root)
       Class(super_comm), intent(in)          :: this
       double precision, intent(inout)        :: a(:, :, :)
       integer, intent(in), optional          :: root
       integer                                :: root_, err_
       character(len=*),parameter  :: subname='[!][global_comm_bcast_d3]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(root)) then
           root_ = root
       else
           root_ = 0
       endif
       call MPI_Bcast(a, size(a,1)*size(a,2)*size(a,3), MPI_DOUBLE, root_, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
   end subroutine global_comm_bcast_d3

!
!  [ALLREDUCE]
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Scalar
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_i(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       integer, intent(in)                    :: a
       character(len=*), intent(in), optional :: oper
       integer                                :: res
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_i]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, 1, MPI_INT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_i
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_f(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       real, intent(in)                       :: a
       character(len=*), intent(in), optional :: oper
       real                                   :: res
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_f]:'
       if (.not. this%mpi_comm_init) then
            write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, 1, MPI_FLOAT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_f
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_d(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       double precision, intent(in)           :: a
       character(len=*), intent(in), optional :: oper
       double precision                       :: res
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_d]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, 1, MPI_DOUBLE, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_d
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Vect
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_i1(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       integer, intent(in)                    :: a(:)
       character(len=*), intent(in), optional :: oper
       integer                                :: res(size(a,1))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_i1]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1), MPI_INT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_i1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_f1(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       real, intent(in)                       :: a(:)
       character(len=*), intent(in), optional :: oper
       real                                   :: res(size(a,1))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_f1]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1), MPI_FLOAT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_f1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_d1(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       double precision, intent(in)           :: a(:)
       character(len=*), intent(in), optional :: oper
       double precision                       :: res(size(a,1))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_d1]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1), MPI_DOUBLE, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_d1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Mat
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_i2(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       integer, intent(in)                    :: a(:,:)
       character(len=*), intent(in), optional :: oper
       integer                                :: res(size(a,1), size(a,2))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_i2]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1)*size(a,2), MPI_INT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_f2(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       real, intent(in)                       :: a(:,:)
       character(len=*), intent(in), optional :: oper
       real                                   :: res(size(a,1), size(a,2))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_f2]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1)*size(a,2), MPI_FLOAT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_f2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_d2(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       double precision, intent(in)           :: a(:,:)
       character(len=*), intent(in), optional :: oper
       double precision                       :: res(size(a,1), size(a,2))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_d2]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1)*size(a,2), MPI_DOUBLE, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_d2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Vol
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_i3(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       integer, intent(in)                    :: a(:,:,:)
       character(len=*), intent(in), optional :: oper
       integer                                :: res(size(a,1), size(a,2), size(a,3))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_i3]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1)*size(a,2)*size(a,3), MPI_INT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_i3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_f3(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       real, intent(in)                       :: a(:,:,:)
       character(len=*), intent(in), optional :: oper
       real                                   :: res(size(a,1), size(a,2), size(a,3))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_f3]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1)*size(a,2)*size(a,3), MPI_FLOAT, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_f3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function global_comm_all_reduce_d3(this, a, oper) result (res)
       Use mpi_f08, Only : MPI_OP
       Class(super_comm), intent(in)          :: this
       double precision, intent(in)           :: a(:,:,:)
       character(len=*), intent(in), optional :: oper
       double precision                       :: res(size(a,1), size(a,2), size(a,3))
       TYPE(MPI_OP)                           :: op
       integer                                :: err_, info
       character(len=*),parameter  :: subname='[!][global_comm_all_reduce_d3]:'
       if (.not. this%mpi_comm_init) then
           write(*,*) subname,': MPI COMM is not INITIALIZED: '; stop
       endif
       if(present(oper)) then
           op = this%get_mpi_op(oper)
       else
           op = this%get_mpi_op("sum")
       endif
       call MPI_Allreduce(a, res, size(a,1)*size(a,2)*size(a,3), MPI_DOUBLE, op, this%comm, err_)
       if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_Allreduce: ',err_;stop;endif
   end function global_comm_all_reduce_d3

!
!  [SYNCHRONIZATION]
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function comm_barrier(this) result (err)
       Class(super_comm), intent(in) :: this
       integer                     :: err
       character(len=*),parameter  :: subname='[!][comm_barrier]:'
       if (this%mpi_comm_init) then
           call MPI_BARRIER(this%comm, err);
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end function
!
!  [finalize]
!
   function finalize_comm(this) result (err)
       Class(super_comm), intent(inout) :: this
       integer                        :: err
       character(len=*),parameter     :: subname='[!][finalize_comm]:'
       if (this%mpi_comm_init) then
           call MPI_FINALIZE(err)
           if(err.ne.0)then;write(*,*)subname,': Problem finalizing MPI COMM: ',err;stop;endif
           this%rank = -1
           this%comm_size = -1
           this%mpi_comm_init = .false.
       else
           write(*,*) subname, "MPI communicator INACTIVE"
           stop
       end if
   end function finalize_comm

   subroutine dealloc_comm(this)
       type(super_comm), intent(inout)  :: this
       integer                        :: err
       character(len=*),parameter     :: subname='[!][dealloc_comm]:'
       ! API: integer(c_int) function destroy_cusolverMpHandle(cusolverMpHandle, VRBZ)
       ! API: integer(c_int) function destroy_cublasMpHandle(handle, VRBZ)
       if (this%mpi_comm_init) then
           call MPI_FINALIZE(err)
           if(err.ne.0)then;write(*,*)subname,': Problem finalizing MPI COMM: ',err;stop;endif
       end if
   end subroutine dealloc_comm

   subroutine reset_cuda_handles(this, VRBZ_)
   class(super_comm), intent(inout)  :: this
   integer, optional                :: VRBZ_
   !
   integer :: err, VRBZ
   character(len=*),parameter     :: subname='[reset_cuda_handles]:'

      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
      call cusolvermp_destroy(VRBZ_=VRBZ)
      call cublasmp_destroy(VRBZ_=VRBZ)
      call cublasmp_init(VRBZ_=VRBZ)
      call cusolvermp_init(this% cuda_device, VRBZ_=VRBZ)

   end subroutine reset_cuda_handles

end module distributed_comms
