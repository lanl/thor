module comm_allreduce
  USE mpi_f08, only : MPI_Comm, MPI_Initialized, MPI_COMM_WORLD, MPI_OP
  implicit none
    interface allreduce
    module procedure      mpi_allreduce_scal_i2,mpi_allreduce_scal_i4,mpi_allreduce_scal_i8,mpi_allreduce_scal_r4,mpi_allreduce_scal_r8, &
                          mpi_allreduce_vect_i2,mpi_allreduce_vect_i4,mpi_allreduce_vect_i8,mpi_allreduce_vect_r4,mpi_allreduce_vect_r8, &
                          mpi_allreduce_mat_i2,mpi_allreduce_mat_i4,mpi_allreduce_mat_i8,mpi_allreduce_mat_r4,mpi_allreduce_mat_r8, &
                          mpi_allreduce_vol_i2,mpi_allreduce_vol_i4,mpi_allreduce_vol_i8,mpi_allreduce_vol_r4,mpi_allreduce_vol_r8
    end interface allreduce
  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Scalars [0D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_allreduce_scal_i2(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER2
     USE comm_types, Only : get_mpi_op
     integer(2), intent(IN)        :: a
     integer(2)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_scal_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, 1, MPI_INTEGER2, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_scal_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_scal_i4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER4
     USE comm_types, Only : get_mpi_op
     integer(4), intent(IN)        :: a
     integer(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_scal_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, 1, MPI_INTEGER4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_scal_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_scal_i8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER8
     USE comm_types, Only : get_mpi_op
     integer(8), intent(IN)        :: a
     integer(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_scal_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, 1, MPI_INTEGER8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_scal_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_scal_r4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL4
     USE comm_types, Only : get_mpi_op
     real(4), intent(IN)        :: a
     real(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_scal_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, 1, MPI_REAL4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_scal_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_scal_r8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL8
     USE comm_types, Only : get_mpi_op
     real(8), intent(IN)        :: a
     real(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_scal_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, 1, MPI_REAL8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_scal_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Vectors [1D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_allreduce_vect_i2(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER2
     USE comm_types, Only : get_mpi_op
     integer(2), intent(IN)        :: a(:)
     integer(2)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vect_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1), MPI_INTEGER2, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vect_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vect_i4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER4
     USE comm_types, Only : get_mpi_op
     integer(4), intent(IN)        :: a(:)
     integer(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vect_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1), MPI_INTEGER4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vect_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vect_i8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER8
     USE comm_types, Only : get_mpi_op
     integer(8), intent(IN)        :: a(:)
     integer(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vect_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1), MPI_INTEGER8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vect_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vect_r4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL4
     USE comm_types, Only : get_mpi_op
     real(4), intent(IN)        :: a(:)
     real(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vect_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1), MPI_REAL4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vect_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vect_r8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL8
     USE comm_types, Only : get_mpi_op
     real(8), intent(IN)        :: a(:)
     real(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vect_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1), MPI_REAL8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vect_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Matrices [2D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_allreduce_mat_i2(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER2
     USE comm_types, Only : get_mpi_op
     integer(2), intent(IN)        :: a(:,:)
     integer(2)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_mat_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2), MPI_INTEGER2, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_mat_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_mat_i4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER4
     USE comm_types, Only : get_mpi_op
     integer(4), intent(IN)        :: a(:,:)
     integer(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_mat_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2), MPI_INTEGER4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_mat_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_mat_i8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER8
     USE comm_types, Only : get_mpi_op
     integer(8), intent(IN)        :: a(:,:)
     integer(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_mat_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2), MPI_INTEGER8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_mat_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_mat_r4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL4
     USE comm_types, Only : get_mpi_op
     real(4), intent(IN)        :: a(:,:)
     real(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_mat_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2), MPI_REAL4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_mat_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_mat_r8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL8
     USE comm_types, Only : get_mpi_op
     real(8), intent(IN)        :: a(:,:)
     real(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_mat_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2), MPI_REAL8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_mat_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Volumes [3D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mpi_allreduce_vol_i2(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER2
     USE comm_types, Only : get_mpi_op
     integer(2), intent(IN)        :: a(:,:,:)
     integer(2)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vol_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER2, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vol_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vol_i4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER4
     USE comm_types, Only : get_mpi_op
     integer(4), intent(IN)        :: a(:,:,:)
     integer(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vol_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vol_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vol_i8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_INTEGER8
     USE comm_types, Only : get_mpi_op
     integer(8), intent(IN)        :: a(:,:,:)
     integer(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vol_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vol_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vol_r4(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL4
     USE comm_types, Only : get_mpi_op
     real(4), intent(IN)        :: a(:,:,:)
     real(4)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vol_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2)*size(a,3), MPI_REAL4, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vol_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpi_allreduce_vol_r8(a, oper, comm) result (res)
     USE mpi_f08, only:  MPI_ALLREDUCE, MPI_REAL8
     USE comm_types, Only : get_mpi_op
     real(8), intent(IN)        :: a(:,:,:)
     real(8)                    :: res
     character(len=*), intent(in), optional :: oper
     !integer, intent(IN),  optional        :: err
     TYPE(MPI_Comm), intent(IN), optional  :: comm
     TYPE(MPI_OP)                          :: op
     TYPE(MPI_Comm)                        :: comm_
     integer                               :: err_
     logical                               :: is_initialized
     character(len=*),parameter            :: subname='[!][mpi_allreduce_vol_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then; comm_=comm; else; comm_ = MPI_COMM_WORLD; end if
         if(present(oper)) then; op=get_mpi_op(oper); else; op = get_mpi_op("sum"); end if
         call MPI_ALLREDUCE(a, res, size(a,1)*size(a,2)*size(a,3), MPI_REAL8, op, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_allreduce: ',err_;stop;endif
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
       !if(present(err)) err=err_
   end function mpi_allreduce_vol_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
end module comm_allreduce
