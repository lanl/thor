module comm_bcast
  USE mpi_f08, only : MPI_Comm, MPI_Initialized, MPI_COMM_WORLD
  implicit none
    interface bcast
    module procedure      mpi_bcast_scal_i2,mpi_bcast_scal_i4,mpi_bcast_scal_i8,mpi_bcast_scal_r4,mpi_bcast_scal_r8, &
                          mpi_bcast_vect_i2,mpi_bcast_vect_i4,mpi_bcast_vect_i8,mpi_bcast_vect_r4,mpi_bcast_vect_r8, &
                          mpi_bcast_mat_i2,mpi_bcast_mat_i4,mpi_bcast_mat_i8,mpi_bcast_mat_r4,mpi_bcast_mat_r8, &
                          mpi_bcast_vol_i2,mpi_bcast_vol_i4,mpi_bcast_vol_i8,mpi_bcast_vol_r4,mpi_bcast_vol_r8
    end interface bcast
  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Scalars [0D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_bcast_scal_i2(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER2
     integer(2)    :: a
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_scal_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, 1, MPI_INTEGER2, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_scal_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_scal_i4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER4
     integer(4)    :: a
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_scal_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, 1, MPI_INTEGER4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_scal_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_scal_i8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER8
     integer(8)    :: a
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_scal_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, 1, MPI_INTEGER8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_scal_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_scal_r4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL4
     real(4)    :: a
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_scal_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, 1, MPI_REAL4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_scal_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_scal_r8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL8
     real(8)    :: a
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_scal_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, 1, MPI_REAL8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_scal_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Vectors [1D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_bcast_vect_i2(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER2
     integer(2)    :: a(:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vect_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1), MPI_INTEGER2, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vect_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vect_i4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER4
     integer(4)    :: a(:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vect_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1), MPI_INTEGER4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vect_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vect_i8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER8
     integer(8)    :: a(:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vect_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1), MPI_INTEGER8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vect_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vect_r4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL4
     real(4)    :: a(:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vect_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1), MPI_REAL4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vect_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vect_r8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL8
     real(8)    :: a(:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vect_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1), MPI_REAL8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vect_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Matrices [2D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_bcast_mat_i2(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER2
     integer(2)    :: a(:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_mat_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2), MPI_INTEGER2, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_mat_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_mat_i4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER4
     integer(4)    :: a(:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_mat_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2), MPI_INTEGER4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_mat_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_mat_i8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER8
     integer(8)    :: a(:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_mat_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2), MPI_INTEGER8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_mat_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_mat_r4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL4
     real(4)    :: a(:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_mat_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2), MPI_REAL4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_mat_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_mat_r8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL8
     real(8)    :: a(:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_mat_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2), MPI_REAL8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_mat_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! [    Volumes [3D]   ]  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_bcast_vol_i2(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER2
     integer(2)    :: a(:,:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vol_i2]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER2, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vol_i2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vol_i4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER4
     integer(4)    :: a(:,:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vol_i4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vol_i4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vol_i8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_INTEGER8
     integer(8)    :: a(:,:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vol_i8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2)*size(a,3), MPI_INTEGER8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vol_i8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vol_r4(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL4
     real(4)    :: a(:,:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vol_r4]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2)*size(a,3), MPI_REAL4, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vol_r4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpi_bcast_vol_r8(a, root, comm, err)
     USE mpi_f08, only:  MPI_BCAST, MPI_REAL8
     real(8)    :: a(:,:,:)
     TYPE(MPI_Comm), intent(IN), optional :: comm
     integer, intent(IN), optional        :: root 
     integer, intent(INOUT), optional     :: err
     TYPE(MPI_Comm)                       :: comm_
     integer                              :: root_, err_
     logical                              :: is_initialized
     character(len=*),parameter           :: subname='[!][mpi_bcast_vol_r8]:'
     call MPI_Initialized(is_initialized)
       if (is_initialized) then
         if(present(comm)) then
           comm_ = comm
         else
           comm_ = MPI_COMM_WORLD 
         end if
         if(present(root)) then
           root_ = root
         else
           root_ = 0
         endif
         call MPI_BCAST(a, size(a,1)*size(a,2)*size(a,3), MPI_REAL8, root_, comm_, err_)
         if(err_ .ne. 0)then;write(*,*) subname,': Unable to compute MPI_bcast: ',err_;stop;endif
         if(present(err)) err = err_
       else
         write(*,*) subname, "MPI communicator INACTIVE"; stop
       end if
   end subroutine mpi_bcast_vol_r8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
end module comm_bcast
