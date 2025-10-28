module comm_send
  USE mpi_f08, only : MPI_Comm, MPI_Initialized, MPI_COMM_WORLD
  implicit none
    interface send
    module procedure      mpi_send_scal_i2,mpi_send_scal_i4,mpi_send_scal_i8,mpi_send_scal_r4,mpi_send_scal_r8, &
                          mpi_send_vect_i2,mpi_send_vect_i4,mpi_send_vect_i8,mpi_send_vect_r4,mpi_send_vect_r8, &
                          mpi_send_mat_i2,mpi_send_mat_i4,mpi_send_mat_i8,mpi_send_mat_r4,mpi_send_mat_r8, &
                          mpi_send_vol_i2,mpi_send_vol_i4,mpi_send_vol_i8,mpi_send_vol_r4,mpi_send_vol_r8
    end interface send
  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Scalars [0D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_send_scal_i2(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER2
     integer(2), intent(IN)       :: a
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_scal_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, 1, MPI_INTEGER2, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_scal_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_scal_i4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER4
     integer(4), intent(IN)       :: a
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_scal_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, 1, MPI_INTEGER4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_scal_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_scal_i8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER8
     integer(8), intent(IN)       :: a
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_scal_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, 1, MPI_INTEGER8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_scal_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_scal_r4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL4
     real(4), intent(IN)       :: a
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_scal_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, 1, MPI_REAL4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_scal_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_scal_r8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL8
     real(8), intent(IN)       :: a
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_scal_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, 1, MPI_REAL8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_scal_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Vectors [1D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_send_vect_i2(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER2
     integer(2), intent(IN)       :: a(:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vect_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1), MPI_INTEGER2, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vect_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vect_i4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER4
     integer(4), intent(IN)       :: a(:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vect_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1), MPI_INTEGER4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vect_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vect_i8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER8
     integer(8), intent(IN)       :: a(:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vect_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1), MPI_INTEGER8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vect_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vect_r4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL4
     real(4), intent(IN)       :: a(:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vect_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1), MPI_REAL4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vect_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vect_r8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL8
     real(8), intent(IN)       :: a(:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vect_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1), MPI_REAL8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vect_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Matrices [2D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_send_mat_i2(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER2
     integer(2), intent(IN)       :: a(:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_mat_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2), MPI_INTEGER2, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_mat_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_mat_i4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER4
     integer(4), intent(IN)       :: a(:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_mat_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2), MPI_INTEGER4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_mat_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_mat_i8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER8
     integer(8), intent(IN)       :: a(:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_mat_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2), MPI_INTEGER8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_mat_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_mat_r4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL4
     real(4), intent(IN)       :: a(:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_mat_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2), MPI_REAL4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_mat_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_mat_r8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL8
     real(8), intent(IN)       :: a(:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_mat_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2), MPI_REAL8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_mat_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Volumes [3D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_send_vol_i2(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER2
     integer(2), intent(IN)       :: a(:,:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vol_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER2, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vol_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vol_i4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER4
     integer(4), intent(IN)       :: a(:,:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vol_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vol_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vol_i8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_INTEGER8
     integer(8), intent(IN)       :: a(:,:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vol_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vol_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vol_r4(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL4
     real(4), intent(IN)       :: a(:,:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vol_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_REAL4, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vol_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_send_vol_r8(a, dest, comm, tag, err)
     USE mpi_f08, only:  MPI_SEND, MPI_REAL8
     real(8), intent(IN)       :: a(:,:,:)
     integer, intent(IN)                  :: dest
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN),    optional     :: tag
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: err_, tag_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_send_vol_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         tag_ = 123456789; if(present(tag)) tag_ = tag
         call MPI_SEND(a, size(a,1)*size(a,2)*size(a,3), MPI_REAL8, dest, tag_, comm_, err_)
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_send_vol_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
end module comm_send
