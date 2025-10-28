!#define VERBOSE4

module distributed_arrays
use iso_c_binding, only: c_ptr
use cusolvermp_aux
use cublasmp_aux
use cudafor
use matrix_ops
use distributed_comms, only      : super_comm, mpi_sub_comm
implicit none
type,public :: dPart ! Distributed 3D array
   ! COMM Info
   type(super_comm), pointer    :: comm => null()
   ! Global dimension
   integer :: nx = -1       !< Number of elements of global array along x (# rows if nz=1)
   integer :: ny = -1       !< Number of elements of global array along y (# cols if nz=1)
   integer :: nz = -1       !< Number of elements of global array along z
   ! Blocking dimension
   integer :: nxb  = -1     !< BLock size used to distribute along X
   integer :: nyb  = -1     !< BLock size used to distribute along y
   integer :: nzb  = -1     !< BLock size used to distribute along z
   ! Local size
   integer :: nxl  = -1     !< Number of elements of local partion along x
   integer :: nyl  = -1     !< Number of elements of local partion along y
   integer :: nzl  = -1     !< Number of elements of local partion along z
   ! Rank, rank size, and partition axis
   integer :: part_axis = 3 !< Axis normal to partition plane {1,2,3} -> {x,y,z}, -1 -> undefined
   integer :: rank          !< local GPU-capable rank
   integer :: nranks        !< total number of GPU-capable ranks
   ! Local data
   double precision,pointer,contiguous,device :: d_A(:) !=> null()! local data on device
   double precision,pointer,contiguous        :: h_A(:) => null() ! local data on the host
   ! Status (TODO: get rid of these)
   logical :: is_allocated = .false.! Allocation state
   logical :: is_init      = .false.! Flag to check whether dRay was initialized

 contains
!   procedure, pass(self) :: dPart_assign
!   generic   :: assignment(=) => dPart_assign
   procedure :: allocate_arrays => dPart_allocate_arrays
   final     :: dealloc_dPart
end type dPart


!> The "Stack": 3D partition flattened to a matrix suitable for algebraic operations
type,public, extends(dPart)    :: stack
   ! COMM Info
   double precision,pointer,contiguous,device :: d_A2(:,:)           ! 2D on device
   double precision,pointer,contiguous        :: h_A2(:,:) => null() ! 2D on host
   integer :: M, N          !< # of rows, columns of the global matrix
   integer :: MB, NB        !< vertical/horizonal row/col blocking size
   integer :: ML, NL        !< local # or rows, columns
   integer :: numRowDevices !< Number of devices allong the rows of the device grid
   integer :: numColDevices !< Number of devices allong the columns of the device grid
   integer :: dev_row_idx   !< Rank GPU/Device Row   index
   integer :: dev_col_idx   !< Rank GPU/Device Column index
 contains
   procedure, pass(self) :: stack_assign
   generic :: assignment(=) => stack_assign
   procedure :: is_hstack => stack_is_hstack ! Slices of dPart are stacked horizontally
   procedure :: is_vstack => stack_is_vstack ! Slices of dPart are stacked vertically
   procedure :: well_stacked => stack_well_stacked
   procedure :: pprint => stack_pprint
   procedure :: cusolvermp_grid_handle   => stack_cusolvermp_grid_handle
   procedure :: cublasmp_grid_handle     => stack_cublasmp_grid_handle
   procedure :: cusolvermp_matrix_handle => stack_cusolvermp_matrix_handle
   procedure :: cublasmp_matrix_handle   => stack_cublasmp_matrix_handle
   procedure :: qr => stack_qr
   procedure :: transpose_slices => stack_transpose_slices
   procedure :: transpose_matrix => stack_transpose_matrix
   procedure :: switch_stacking
   procedure :: to_hstack => stack_to_hstack
   procedure :: to_vstack => stack_to_vstack
   procedure :: redist_to_vstack => stack_redist_to_vstack
   procedure :: redist_to_hstack => stack_redist_to_hstack
   procedure :: nrm2 => stack_nrm2
   procedure :: truncate_matrix => stack_truncate_matrix
   procedure :: allocate_arrays => stack_allocate_arrays
   procedure :: add_to_diag => stack_add_to_diag
   procedure :: svd => stack_svd
   procedure :: sytrd => stack_sytrd
   procedure :: save_to_ascii => stack_save_to_ascii
   final     :: dealloc_stack
end type stack
 !/////////////////////////////////////////////////////////////////////////////////
 !                               OPERATORS / OVERLOADING
 !/////////////////////////////////////////////////////////////////////////////////
  interface operator (+)
       module procedure add_stacks
  end interface
  interface operator (-)
       module procedure subtract_stacks
  end interface
  interface assignment (=)
       module procedure dPart_assign_dPart, &
                        dPart_assign_dble, &
                        dPart_assign_real, &
                        dPart_assign_int
  end interface
  ! interface operator (*)
  !     module procedure mul_dParts, mul_arr_by_dPart, mul_dPart_by_array
  ! end interface
  ! interface operator (/)
  !     module procedure div_dParts, div_arr_by_dPart, div_dPart_by_array
  ! end interface
  ! interface operator (**)
  !     module procedure dPart_to_pow
  ! end interface
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !/////////////////////////////////////////////////////////////////////////////////
 !                              CONSTRUCTOR\DESTRUCTO INTERFACE
 !/////////////////////////////////////////////////////////////////////////////////
  interface dPart ! constructor interface
     procedure :: constr_dPart
     procedure :: dPart_from_dPart  ! Deep copy
     procedure :: constr_dPart_simplified
     procedure :: constr_dPart_from_array3
  end interface dPart
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface stack
     procedure :: empty_stack
     procedure :: shallow_stack
     procedure :: stack_from_stack
     procedure :: stack_from_dPart
     procedure :: stack_from_flat_array
     procedure :: stack_from_array3
     procedure :: stack_read_from_ascii
  end interface stack
 !/////////////////////////////////////////////////////////////////////////////////
 !                              OPERATORS / OVERLOADING
 !/////////////////////////////////////////////////////////////////////////////////
 ! interface operator (+)
 !     module procedure dtt_add, ztt_add, dtt_plus_d, d_plus_dtt, dttm_add
 ! end interface
 ! interface operator (-)
 !     module procedure dtt_sub, dtt_minus_d, d_minus_dtt, minus_dtt, &
 !                      ztt_sub
 ! end interface
 ! interface operator (*)
 !     module procedure dtt_mul_d, ztt_mul_z, d_mul_dtt, z_mul_ztt, &
 !                      dtt_mult, ztt_mult
 ! end interface
 ! interface operator (/)
 !     module procedure dtt_div_d
 ! end interface
 ! interface operator (**)
 !     module procedure dtt_pwr_i
 ! end interface


 !/////////////////////////////////////////////////////////////////////////////////
 !                              UTILITY FUNCTIONs
 !/////////////////////////////////////////////////////////////////////////////////
  interface similar
     module procedure similar_dParts, similar_stacks
  end interface
  interface print_info
     module procedure print_info_dPart, print_info_stack
  end interface
  interface zeros_stack
     module procedure zeros_stack_2d, zeros_stack_3d
  end interface
  interface eye_stack
     module procedure eye_stack_2d, eye_stack_3d
  end interface
 contains




!
! UTILITY FUNCTIONS
!
   logical function stack_is_vstack(st)
   class(stack), intent(in) :: st
      stack_is_vstack = (st% N == st% ny) .and. &
                       ((st% M >  st% nx) .or.  &
                        (st% M == st% nx .and. st% numColDevices == 1))
   end function stack_is_vstack


   logical function stack_is_hstack(st)
   class(stack), intent(in) :: st
      stack_is_hstack = (st% M == st% nx) .and. &
                       ((st% N >  st% ny) .or.  &
                        (st% N == st% ny .and. st% numRowDevices == 1))
   end function stack_is_hstack


   !> Check if stack is well stacked
   logical function stack_well_stacked(st)
   class(stack), intent(in) :: st
      stack_well_stacked = (st% is_hstack() .and. st% N == st% ny * st% nz) &
                       .or.(st% is_vstack() .and. st% M == st% nx * st% nz)
   end function stack_well_stacked


   !> Pretty-print stack (mainly for debugging)
   subroutine stack_pprint(this, label_, line_len_, max_rows_)
   use matrix_util, only: pprint_matrix
   use string_lib, only: str
   class(stack),           intent(IN) :: this
   character(*), optional, intent(IN) :: label_
   integer,      optional, intent(in) :: line_len_
   integer,      optional, intent(in) :: max_rows_
   !
   integer :: nr, rank, nranks, line_len, max_rows
   character(200) :: label
   double precision, allocatable :: h_A2(:,:)

      rank = this% rank
      nranks = this% nranks
      label = '['//str(rank)//'/'//str(nranks)//']: '
      if (present(label_)) label = trim(label)//trim(label_)
      max_rows = 12; if (present(max_rows_)) max_rows = max_rows_
      line_len = 80; if (present(line_len_)) line_len = line_len_
      allocate(h_A2(this% ML, this% NL))
      h_A2 = this% d_A2
      do nr = 1, nranks
         if (rank + 1 .eq. nr) then
            print '(/,"'//trim(label)//'")'
            call pprint_matrix(h_A2, max_rows_=max_rows, line_len_=line_len)
         endif
         call this% comm% barrier()
      enddo
      deallocate(h_A2)
   end subroutine stack_pprint


   !> Check if dPart A and B are similar
   logical function similar_dParts(A, B)
   type(dPart), intent(in) :: A, B
   character(*),parameter  :: subnam='[dPart][similar_dParts]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      similar_dParts = (B%nxb==A%nxb).and.(B%nyb==A%nyb).and.(B%nzb==A%nzb) &
                 .and. (B%nxl==A%nxl).and.(B%nyl==A%nyl).and.(B%nzl==A%nzl) &
                 .and.  (B%nx==A%nx) .and. (B%ny==A%ny) .and. (B%nz==A%nz)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function similar_dParts


   !> Print dPart metadata
   subroutine print_info_dPart(D, label_)
   type(dPart), intent(IN)            :: D
   character(*), intent(IN), optional :: label_
   !
   character(*),parameter   :: subnam='[print_info_dPart]'
   integer :: r, nr
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      r = D% rank; nr = D% nranks
      if(present(label_)) print '("["I3"/"I3"]",A)',r,nr,label_
      print '("["I3"/"I3"] {nx, ny, nz}  = {"I7","I7","I7"}")',r,nr, D%nx, D%ny, D%nz
      print '("["I3"/"I3"] {nxb,nyb,nzb} = {"I7","I7","I7"}")',r,nr,D%nxb, D%nyb, D%nzb
      print '("["I3"/"I3"] {nxl,nyl,nzl} = {"I7","I7","I7"}")',r,nr,D%nxl, D%nyl, D%nzl
      print '("["I3"/"I3"] part_axis = "I5)', r,nr, D%part_axis
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine print_info_dPart


   !> Check if stacks A and B are similar
   logical function similar_stacks(A, B, VRBZ_)
   class(stack),      intent(IN) :: A, B
   integer, optional, intent(IN) :: VRBZ_
   !
   integer :: VRBZ
   character(*),parameter      :: subnam='[stack][similar_stacks]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
      similar_stacks = &
                 (B%nxb==A%nxb) .and. (B%nyb==A%nyb) .and. (B%nzb==A%nzb) &
           .and. (B%nxl==A%nxl) .and. (B%nyl==A%nyl) .and. (B%nzl==A%nzl) &
           .and.  (B%nx==A%nx)  .and.  (B%ny==A%ny)  .and.  (B%nz==A%nz)  &
           .and. (B%M==A%M).and.(B%N==A%N) .and. (B%MB==A%MB).and.(B%NB==A%NB) &
           .and. (B%rank==A%rank) .and. (B%nranks==A%nranks) &
           .and. (B%numRowDevices==A%numRowDevices) &
           .and. (B%numColDevices==A%numColDevices) &
           .and. (B%dev_row_idx==A%dev_row_idx) &
           .and. (B%dev_col_idx==A%dev_col_idx)
      if (.not.similar_stacks .and. VRBZ > 0) then
         print '("(B%nxb==A%nxb)                     : "L1)', (B%nxb==A%nxb)
         print '("(B%nyb==A%nyb)                     : "L1)', (B%nyb==A%nyb)
         print '("(B%nzb==A%nzb)                     : "L1)', (B%nzb==A%nzb)
         print '("(B%nxl==A%nxl)                     : "L1)', (B%nxl==A%nxl)
         print '("(B%nyl==A%nyl)                     : "L1)', (B%nyl==A%nyl)
         print '("(B%nzl==A%nzl)                     : "L1)', (B%nzl==A%nzl)
         print '("(B%nx==A%nx)                       : "L1)', (B%nx==A%nx)
         print '("(B%ny==A%ny)                       : "L1)', (B%ny==A%ny)
         print '("(B%nz==A%nz)                       : "L1)', (B%nz==A%nz)
         print '("(B%M==A%M)                         : "L1)', (B%M==A%M)
         print '("(B%N==A%N)                         : "L1)', (B%N==A%N)
         print '("(B%MB==A%MB)                       : "L1)', (B%MB==A%MB)
         print '("(B%NB==A%NB)                       : "L1)', (B%NB==A%NB)
         print '("(B%rank==A%rank)                   : "L1)', (B%rank==A%rank)
         print '("(B%nranks==A%nranks)               : "L1)', (B%nranks==A%nranks)
         print '("(B%numRowDevices==A%numRowDevices) : "L1)', (B%numRowDevices==A%numRowDevices)
         print '("(B%numColDevices==A%numColDevices) : "L1)', (B%numColDevices==A%numColDevices)
         print '("(B%dev_row_idx==A%dev_row_idx)     : "L1)', (B%dev_row_idx==A%dev_row_idx)
         print '("(B%dev_col_idx==A%dev_col_idx)     : "L1)', (B%dev_col_idx==A%dev_col_idx)
      endif
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function similar_stacks


   !> Print stack metadata
   subroutine print_info_stack(A, label_)
   type(stack), intent(IN)            :: A
   character(*), intent(IN), optional :: label_
   !
   character(*),parameter   :: subnam='[print_info_stack]'
   integer :: r, nr, sh(2)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      r = A% rank; nr = A% nranks
      if(present(label_)) print '("["I3"/"I3"]",A)',r,nr,label_
      print '("["I3"/"I3"] {nx,ny,nz} = {"I7","I7","I7"}")',r,nr, A%nx, A%ny, A%nz
      print '("["I3"/"I3"] {nxb,nyb,nzb} = {"I7","I7","I7"}")',r,nr, A%nxb, A%nyb, A%nzb
      print '("["I3"/"I3"] {nxl,nyl,nzl} = {"I7","I7","I7"}")',r,nr, A%nxl, A%nyl, A%nzl
      print '("["I3"/"I3"] {M,N} = {"I7","I7"}")',r,nr, A%M, A%N
      print '("["I3"/"I3"] {MB,NB} = {"I7","I7"}")',r,nr, A%MB, A%NB
      print '("["I3"/"I3"] {ML,NL} = {"I7","I7"}")',r,nr, A%ML, A%NL
      print '("["I3"/"I3"] {rank, nranks} = {"I7","I7"}")',r,nr, A%rank, A%nranks
      print '("["I3"/"I3"] {numRowDevices, numColDevices} = {"I7","I7"}")', &
            r,nr, A%numRowDevices, A%numColDevices
      print '("["I3"/"I3"] {dev_row_idx,dev_col_idx} = {"I7","I7"}")', &
            r,nr, A%dev_row_idx, A%dev_col_idx
      print '("["I3"/"I3"] {is_hstack(), is_vstack()} = {"L1","L1"}")', &
            r,nr, A%is_hstack(), A%is_vstack()
      print '("["I3"/"I3"] size(A% d_A) = "I12"}")', r,nr, size(A%d_A)
      print '("["I3"/"I3"] shape(A% d_A2) = (",2(I7),")}")', r,nr, shape(A%d_A2)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine print_info_stack


!
! CONSTRUCTOR \ DESTRUCTOR
!
   !> Main dPart constructor subroutine
   subroutine construct_dPart(self, comm, n, nb, lalloc, part_axis_, arr_, VRBZ_)
   use matrix_util, only: d1copy
   use cublasmp_aux, only: numroc
   class(dPart),            intent(INOUT) :: self
   type(super_comm), target,   intent(IN) :: comm
   integer,                    intent(IN) :: n(3)       !< global dimensions (order: (nx,ny,nz)
   integer,                    intent(IN) :: nb(3)      !< blocking sizes used to distribute
   logical,                    intent(IN) :: lalloc     !< allocate arrays?
   double precision, optional, intent(IN) :: arr_(*)    !< array data (on device!)
   integer,          optional, intent(IN) :: part_axis_ !< partition axis (1:x, 2:y, 3:z)
   integer,          optional, intent(IN) :: VRBZ_
   !
   integer :: sz, VRBZ
   character(*),parameter :: subnam='[construct_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      if (self% is_allocated) call dealloc_dPart_base(self, VRBZ_)
      self%comm => comm
      self%rank   = comm%ctx_rank
      self%nranks = comm%ctx_size
      self%nx  = n(1);  self%ny  = n(2); self%nz = n(3)
      self%nxb =nb(1);  self%nyb =nb(2); self%nzb=nb(3)

      self%part_axis = 3
      if (present(part_axis_)) self%part_axis = part_axis_
      select case (self%part_axis)
      case(1)
         self% nxl = numroc(n(1), nb(1), self% rank, self% nranks)
         self% nyl = self% ny
         self% nzl = self% nz
      case(2)
         self% nxl = self% nx
         self% nyl = numroc(n(2), nb(2), self% rank, self% nranks)
         self% nzl = self% nz
      case(3)
         self% nxl = self% nx
         self% nyl = self% ny
         self% nzl = numroc(n(3), nb(3), self% rank, self% nranks)
      case default
         error stop '[!]'//subnam//": invalid part_axis"
      end select
      if (VRBZ>0) call print_info(self)

      if (lalloc) then
         call self% allocate_arrays(VRBZ_)
         if (present(arr_)) then
            sz= self% nxl * self% nyl * self% nzl
            if (VRBZ>0) print '("[+]'//subnam//': copying arr_("I7")")', sz
            call d1copy(arr_, self% h_A, sz)
            self% d_A = self% h_A
            if (VRBZ>0 .and. sz >= 4) &
               print '("[+]'//subnam//' self% h_A(1:4) = "(4(F6.2,1X)))', self% h_A(1:4)
         endif
         call comm%barrier()
      else
         if (present(arr_)) error stop '[!]'//subnam//': not allocated'
      endif
      self%is_allocated = lalloc
      self%is_init      = .true. ! TODO: this seems useless
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine construct_dPart


   !> Main dPart constructor function
   function constr_dPart(comm, n, nb, lalloc, part_axis_, arr_, VRBZ_) result(self)
   type(dPart) :: self
   type(super_comm), target,   intent(IN) :: comm
   integer,                    intent(IN) :: n(3)       !< global dimensions (order: (nx,ny,nz)
   integer,                    intent(IN) :: nb(3)      !< blocking sizes used to distribute
   logical,                    intent(IN) :: lalloc     !< allocate arrays?
   integer,          optional, intent(IN) :: part_axis_ !< partition axis (1:x, 2:y, 3:z)
   double precision, optional, intent(IN) :: arr_(*)    !< array data (on device!)
   integer,          optional, intent(IN) :: VRBZ_
   !
   character(*),parameter :: subnam='[constr_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call construct_dPart(self, comm, n, nb, lalloc, part_axis_, arr_, VRBZ_)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function constr_dPart


   !> copy constructor
   function dPart_from_dPart(other, VRBZ_) result(self)
   type(dPart):: self
   type(dPart), intent(IN)       :: other
   integer, intent(IN), optional :: VRBZ_
   !
   character(*),parameter      :: subnam='[dPart_from_dPart]'
   integer :: VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      call construct_dPart(self, other% comm, &
           [other% nx, other% ny, other% nz], &
           [other% nxb, other% nyb, other% nzb], .true., &
           part_axis_=other% part_axis, VRBZ_=VRBZ)
      self% d_A = other% d_A
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function dPart_from_dPart


   !> Simplified dPart constructor: assumes nb := [ng/nranks] as a blocking size
   subroutine dPart_simplified(self, comm, ng, lalloc, part_axis_, arr_, mb_, VRBZ_)
   use matrix_util, only: d1copy
   class(dPart),            intent(INOUT) :: self
   type(super_comm), target,   intent(IN) :: comm
   integer,                    intent(IN) :: ng(3)
   logical,                    intent(IN) :: lalloc
   integer,          optional, intent(IN) :: part_axis_
   double precision, optional, intent(IN) :: arr_(*)
   integer,          optional, intent(IN) :: mb_   !< block size used to distribute
   integer,          optional, intent(IN) :: VRBZ_
   !
   integer :: nb(3), part_axis, nranks, VRBZ
   character(*),parameter      :: subnam='[dPart_simplified]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      if (VRBZ>0) then
         print '("[i]'//subnam//': ng="3(I7,1X)", lalloc="L1)', ng,lalloc
         if (present(part_axis_)) &
            print '("[i]'//subnam//': part_axis_= "I1)', part_axis_
         if (present(arr_)) print '("[i]'//subnam//': arr_ is present")'
         if (present(mb_))  print '("[i]'//subnam//': mb_ = "I7)', mb_
      endif
      nb = ng
      nranks = comm% ctx_size
      part_axis = 3; if(present(part_axis_)) part_axis = part_axis_
      if (present(mb_)) then
         nb(part_axis) = mb_
      else
         nb(part_axis) = (ng(part_axis) - 1)/nranks + 1
      endif
      call construct_dPart(self, comm, ng, nb, lalloc, part_axis_, arr_, VRBZ_)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dPart_simplified


   !> Function wrapper to dPart_simplified
   function constr_dPart_simplified(comm, ng, lalloc, part_axis_, arr_, mb_, VRBZ_) result(self)
   type(dPart):: self
   type(super_comm), target,   intent(IN) :: comm
   integer,                    intent(IN) :: ng(3)
   logical,                    intent(IN) :: lalloc
   integer,          optional, intent(IN) :: part_axis_
   double precision, optional, intent(IN) :: arr_(*)
   integer,          optional, intent(IN) :: mb_   !< block size used to distribute
   integer,          optional, intent(IN) :: VRBZ_
   !
   character(*),parameter      :: subnam='[constr_dPart_simplified]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call dPart_simplified(self, comm, ng, lalloc, part_axis_, arr_, mb_, VRBZ_)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function constr_dPart_simplified


   !> Construct dPart from local arrays
   !! Global size along the partition axis is computed by global reduction operation
   function constr_dPart_from_array3(comm, arr, part_axis_, VRBZ_) result(self)
   type(dPart):: self
   type(super_comm), target, intent(IN):: comm
   double precision,         intent(IN):: arr(:,:,:)
   integer,        optional, intent(IN):: part_axis_
   integer,        optional, intent(IN):: VRBZ_
   !
   integer  :: n(3), ng(3), nb(3), VRBZ, part_axis
   character(*),parameter      :: subnam='[constr_dPart_from_array3]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_

      n = shape(arr)
      ng = n; nb = n
      part_axis = 3; if(present(part_axis_)) part_axis = part_axis_
      ng(part_axis) = comm% comm% allreduce(n(part_axis), 'sum')
      nb(part_axis) = comm% comm% allreduce(n(part_axis), 'max')
      call dPart_simplified(self, comm, ng, .true., part_axis_=part_axis, &
         arr_=arr, mb_=nb(part_axis), VRBZ_=VRBZ)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function constr_dPart_from_array3

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![    stack    ]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> Main stack constructor
   subroutine construct_stack(self, comm, ng, nb, stacking, &
                              part_axis_, arr_, morn_, VRBZ_)
   type(stack),             intent(INOUT):: self
   type(super_comm), target,   intent(IN):: comm         !< supercommunicator
   integer,                    intent(IN):: ng(3), nb(3) !< global dimensions and block sizes
   character,                  intent(IN):: stacking     !< 'h' by default (horizontal)
   integer,          optional, intent(IN):: part_axis_   !< partition axis (default=3, z-axis)
   double precision, optional, intent(IN):: arr_(*)      !< array data (on host!)
   integer,          optional, intent(IN):: morn_        !< stack size in the direction of unfolding
   integer,          optional, intent(IN):: VRBZ_
   !
   character(*),parameter :: subnam='[construct_stack]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call construct_dPart(self, comm, ng, nb, .not.present(morn_), &
                           part_axis_, arr_, VRBZ_)
      select case(stacking)
      case ('h','H')
         call hstack_dPart(self, morn_, VRBZ_)
      case ('v','V')
         call vstack_dPart(self, morn_, VRBZ_)
      end select
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine construct_stack


   function empty_stack(comm, ng, nb, stacking, part_axis_, arr_, morn_, VRBZ_) result(self)
   type(stack):: self
   type(super_comm), target,   intent(IN) :: comm         !< supercommunicator
   integer,                    intent(IN) :: ng(3), nb(3) !< global dimensions and block sizes
   character,                  intent(IN) :: stacking     !< 'h' by default (horizontal)
   integer,   optional,        intent(IN) :: part_axis_   !< partition axis (default=3, z-axis)
   double precision, optional, intent(IN) :: arr_(*)      !< array data (on host!)
   integer,   optional,        intent(IN) :: morn_        !< size in the direction of unfolding
   integer,   optional,        intent(IN) :: VRBZ_
   !
   character(*),parameter    :: subnam='[empty_stack]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call construct_stack(self,comm,ng,nb,stacking,part_axis_,arr_,morn_,VRBZ_)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function empty_stack


   !> Creates a "shallow" stack with nz = n_dev, nzb = 1, and nzl = 1
   function shallow_stack(comm, M, N, stacking, MB_, part_axis_, val_, VRBZ_) result(self)
   type(stack):: self
   type(super_comm), target,   intent(IN):: comm
   integer,                    intent(IN):: M, N     !< size of the matrix (M rows, N columns)
   character,                  intent(IN):: stacking !< 'h'/'H' or 'v'/'V' ['h']
   integer,          optional, intent(IN):: MB_      !< block size (default: determined by nranks)
   integer,          optional, intent(IN):: part_axis_ !< partition axis (default=3, z-axis)
   double precision, optional, intent(IN):: val_     !< which value to assign [default: nothing]
   integer,          optional, intent(IN):: VRBZ_
   !
   character(*),parameter :: subnam='[shallow_stack]'
   integer :: nz, M_B, N_B, MorN, MorN_B, nb(3), ng(3), VRBZ, part_axis
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      nz = comm% ctx_size
      select case(stacking)
      case ('h','H')
         M_B = M
         N_B = (N - 1)/nz + 1
         if (present(MB_)) then
            N_B = MB_
            nz  = (N - 1)/N_B + 1
         endif
         MorN = N
      case ('v','V')
         M_B = (M - 1)/nz + 1
         N_B = N
         if (present(MB_)) then
            M_B = MB_
            nz = (M - 1)/M_B + 1
         endif
         MorN = M
      case default
         error stop '[!]'//subnam//': bad stacking parameter'
      end select
      part_axis = 3; if(present(part_axis_)) part_axis = part_axis_
      select case(part_axis)
      case(1)
         ng = [nz, M_B, N_B]
         nb = [1,  M_B, N_B]
      case(2)
         ng = [M_B, nz, N_B]
         nb = [M_B, 1,  N_B]
      case(3)
         ng = [M_B, N_B, nz]
         nb = [M_B, N_B, 1]
      case default
         error stop '[!]'//subnam//': bad part_axis_ parameter'
      end select
      call construct_stack(self, comm, ng, nb, stacking,  &
                           part_axis_=part_axis, morn_=MorN, &
                           VRBZ_=VRBZ)
      if (present(val_)) then
         self% h_A = val_
         self% d_A = self% h_A
      endif
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function shallow_stack


   !> Zero stack constructor
   function zeros_stack_3d(comm, ng, nb, stacking, part_axis_, morn_, VRBZ_) result(self)
   type(stack), self
   type(super_comm), target, intent(IN):: comm         !< supercommunicator
   integer,                  intent(IN):: ng(3), nb(3) !< global dimensions and block sizes
   character,                intent(IN):: stacking     !< 'h' by default (horizontal)
   integer,        optional, intent(IN):: part_axis_   !< partition axis (default=3, z-axis)
   integer,        optional, intent(IN):: morn_        !< stack blocking dimension
   integer,        optional, intent(IN):: VRBZ_
   !
   character(*),parameter      :: subnam='[zeros_stack_3d]'
   integer :: part_axis, VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      part_axis = 3;  if(present(part_axis_)) part_axis = part_axis_
      if (present(morn_)) then
         call construct_stack(self, comm, ng, nb, stacking, &
                              part_axis_=part_axis, morn_=morn_, VRBZ_=VRBZ)
         if (VRBZ>2) print '("[i]'//subnam//': calling construct_stack with morn_")'
      else
         call construct_stack(self, comm, ng, nb, stacking, &
                              part_axis_=part_axis, VRBZ_=VRBZ)
         if (VRBZ>2) print '("[i]'//subnam//': calling construct_stack w/o morn_")'
      endif
      self% h_A = 0d0
      self% d_A = self% h_A
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function zeros_stack_3d


   !> Shallow stack filled with zeros
   function zeros_stack_2d(comm, M, N, stacking, MB_, part_axis_, VRBZ_) result(self)
   type(stack):: self
   type(super_comm), target,   intent(IN):: comm
   integer,                    intent(IN):: M, N     !< size of the matrix (M rows, N columns)
   character,                  intent(IN):: stacking !< 'h'/'H' or 'v'/'V' ['h']
   integer,          optional, intent(IN):: MB_      !< block size (default: determined by nranks)
   integer,          optional, intent(IN):: part_axis_ !< partition axis (default=3, z-axis)
   integer,          optional, intent(IN):: VRBZ_
   !
   character(*),parameter      :: subnam='[zeros_stack_2d]'
   integer :: MB, part_axis, VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      part_axis = 3;  if(present(part_axis_)) part_axis = part_axis_
      if (present(MB_)) then
         MB = mb_
         self= shallow_stack(comm, M, N, stacking, MB_=MB, &
                             part_axis_=part_axis, val_=0d0, VRBZ_=VRBZ)
      else
         self= shallow_stack(comm, M, N, stacking, &
                             part_axis_=part_axis, val_=0d0, VRBZ_=VRBZ)
      endif
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function zeros_stack_2d


   !> Copy constructor
   function stack_from_stack(other) result(self)
   type(stack):: self
   type(stack), target, intent(IN) :: other
   !
   character(*), parameter :: subnam='[stack_from_stack]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (self% is_allocated) call dealloc_stack(self)
      call hard_copy_stack(self, other)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_from_stack


   !> Constructor from dPart
   function stack_from_dPart(other, stacking, VRBZ_) result(self)
   type(stack):: self
   type(dPart), target, intent(IN) :: other
   character,           intent(IN) :: stacking
   integer,   optional, intent(IN) :: VRBZ_
   !
   character(*), parameter :: subnam='[stack_from_dPart]'
   integer :: VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
      call hard_copy_dPart(self, other)
      select case (stacking)
      case ('h', 'H')
         call hstack_dPart(self, VRBZ_=VRBZ)
      case ('v', 'V')
         call vstack_dPart(self, VRBZ_=VRBZ)
      case default
         error stop '[!]'//subnam//': invalid stacking parameter'
      end selecT
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_from_dPart


   !> Stack constructor from a flat array with a character argument
   function stack_from_flat_array(comm, ng, arr, stacking, part_axis_, VRBZ_) result(self)
   type(stack):: self
   type(super_comm), target, intent(IN):: comm
   integer,                  intent(IN):: ng(3)
   double precision,         intent(IN):: arr(*)
   character,                intent(IN):: stacking
   integer,        optional, intent(IN):: part_axis_
   integer,        optional, intent(IN):: VRBZ_
   !
   character(*),parameter :: subnam='[stack_from_flat_array]'
   integer :: VRBZ, part_axis
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      part_axis = 3; if(present(part_axis_)) part_axis = part_axis_
      call dPart_simplified(self, comm, ng, .true., part_axis_=part_axis, &
                            arr_=arr, VRBZ_=VRBZ)
      select case (stacking)
      case ('h','H')
         call hstack_dPart(self, VRBZ_=VRBZ)
      case ('v','V')
         call vstack_dPart(self, VRBZ_=VRBZ)
      case default
         error stop '[!]'//subnam//': stacking ican be aither "H" or "V"'
      end select

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_from_flat_array


   !> Stack constructor from a 3D-array
   function stack_from_array3(comm, arr, stacking, part_axis_, VRBZ_) result(self)
   type(stack):: self
   type(super_comm), target, intent(IN):: comm
   double precision,         intent(IN):: arr(:,:,:)
   character,                intent(IN):: stacking
   integer,        optional, intent(IN):: part_axis_
   integer,        optional, intent(IN):: VRBZ_
   !
   character(*),parameter :: subnam='[stack_from_array3]'
   integer :: VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      self = constr_dPart_from_array3(comm, arr, part_axis_, VRBZ_)
      select case (stacking)
      case ('h', 'H')
         call hstack_dPart(self, VRBZ_=VRBZ)
      case ('v', 'V')
         call vstack_dPart(self, VRBZ_=VRBZ)
      case default
         error stop '[!]'//subnam//': stacking can be either "H" or "V"'
      end select
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_from_array3


!
!  [COPY]
!
   !> Soft copy for dPart
   subroutine soft_copy_dPart(B, A, VRBZ_)
   class(dPart),   intent(INOUT) ::  B
   class(dPart),      intent(IN) ::  A
   integer, optional, intent(IN) :: VRBZ_
   !
   character(*), parameter :: subnam = '[soft_copy_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if(.not. associated(A%comm)) error stop '[!]'//subnam//'A%comm not associated'
      if(associated(A%comm)) B%comm => A%comm
      B%nx=A%nx;   B%ny=A%ny;   B%nz=A%nz;
      B%nxb=A%nxb; B%nyb=A%nyb; B%nzb=A%nzb;
      B%nxl=A%nxl; B%nyl=A%nyl; B%nzl=A%nzl;
      B%part_axis = A%part_axis
      B%rank = A%rank
      B%nranks = A%nranks
      if(associated(A%d_A))  B%d_A  => A%d_A
      if(associated(A%h_A))  B%h_A  => A%h_A
      B%is_init=A%is_init
      B%is_allocated=A%is_allocated
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine soft_copy_dPart


   !> Soft copy for a Stack structure
   subroutine soft_copy_stack(B, A, VRBZ_)
   class(stack),   intent(INOUT) ::  B
   class(stack),      intent(IN) ::  A
   integer, optional, intent(IN) :: VRBZ_
   !
   character(*), parameter :: subnam = '[soft_copy_stack]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call soft_copy_dPart(B, A, VRBZ_)
      B%M = A%M; B%MB = A%MB; B%ML = A%ML
      B%N = A%N; B%NB = A%NB; B%NL = A%NL
      B%numRowDevices = A%numRowDevices
      B%numColDevices = A%numColDevices
      B%dev_row_idx   = A%dev_row_idx
      B%dev_col_idx   = A%dev_col_idx
      if (associated(A%d_A2)) B%d_A2 => A%d_A2
      if (associated(A%h_A2)) B%h_A2 => A%h_A2
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine soft_copy_stack


   !> Hard copy for dPart structure
   subroutine hard_copy_dPart(B, A)
   class(dPart),      intent(INOUT) :: B
   class(dPart), target, intent(IN) :: A
   !
   integer :: sz
   character(*), parameter :: subnam = '[hard_copy_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call soft_copy_dPart(B, A)
      if(associated(A%d_A)) then
         sz= size(A% d_A)
         allocate(B% d_A(sz))
         B% d_A= A% d_A
      end if
      if(associated(A%h_A)) then
         sz= size(A% h_A)
         allocate(B% h_A(sz))
         B% h_A= A% h_A
      end if
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine hard_copy_dPart


   !> Hard copy for the Stack structure
   subroutine hard_copy_stack(B, A)
   type(stack),       intent(INOUT) :: B
   class(stack), target, intent(IN) :: A
   !
   integer :: sz, sh(2)
   character(*), parameter :: subnam = '[hard_copy_stack]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call soft_copy_stack(B, A)
      !! Hard copy Pointers
      if(associated(A%d_A)) then
         sz= size(A% d_A)
         allocate(B% d_A(sz))
         B% d_A= A% d_A
         if(associated(A%d_A2)) then
             sh = shape(A%d_A2)
             B%d_A2(1:sh(1),1:sh(2)) => B% d_A
         end if
      end if
      if(associated(A%h_A)) then
         sz= size(A% h_A)
         allocate(B% h_A(sz))
         B% h_A= A% h_A
         if(associated(A%h_A2)) then
             sh = shape(A%h_A2)
             B%h_A2(1:sh(1),1:sh(2)) => B% h_A
         end if
      end if
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine hard_copy_stack

!
!  [ASSIGNEMENT]     A = B
!
   !> Assignment operator function for dPart type
   subroutine dPart_assign_dPart(self, other)
   class(dPart), intent(INOUT) :: self
   type(dPart),     intent(IN) :: other
   !
   character(*), parameter :: subnam = '[dPart_assign_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (self% is_allocated) call dealloc_dPart(self)
      call hard_copy_dPart(self, other)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dPart_assign_dPart


   !> Assignment operator function for dPart type
   subroutine dPart_assign_dble(self, x)
   class(dPart),  intent(INOUT) :: self
   double precision, intent(IN) :: x
   !
   character(*), parameter :: subnam = '[dPart_assign_dble]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (.not.self% is_allocated) error stop '[!]'c//subnam//': not allocated'
      self% d_A = x
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dPart_assign_dble


   !> Assignment operator function for dPart type
   subroutine dPart_assign_real(self, r)
   class(dPart),  intent(INOUT) :: self
   real, intent(IN) :: r
   !
   character(*), parameter :: subnam = '[dPart_assign_real]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (.not.self% is_allocated) error stop '[!]'c//subnam//': not allocated'
      self% d_A = dble(r)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dPart_assign_real


   !> Assignment operator function for dPart type
   subroutine dPart_assign_int(self, n)
   class(dPart),  intent(INOUT) :: self
   integer, intent(IN) :: n
   !
   character(*), parameter :: subnam = '[dPart_assign_int]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (.not.self% is_allocated) error stop '[!]'c//subnam//': not allocated'
      self% d_A = dble(n)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dPart_assign_int


   !> Assignment operator function for Stack type
   subroutine stack_assign(self, other)
   class(stack), intent(INOUT) :: self
   type(stack), intent(in), target :: other
   !
   character(*), parameter :: subnam = '[stack_assign]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (self% is_allocated) call dealloc_stack(self) ! TODO: Write a proper dealloc func
      call hard_copy_stack(self, other)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_assign

!
!  [ALLOCATION / DEALLOCATION]
!
   !> Allocation data arrays in dPart
   subroutine dPart_allocate_arrays(self, VRBZ_)
   class(dPart), intent(INOUT) :: self
   integer, intent(in), optional :: VRBZ_
   !
   character(*), parameter :: subnam = '[dPart_allocate_arrays]'
   integer :: sz, VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      sz= max(self%nxl * self%nyl * self%nzl, 1)
      if (associated(self% d_A)) deallocate(self% d_A)
      if (associated(self% h_A)) deallocate(self% h_A)
      allocate(self%d_A(sz))
      allocate(self%h_A(sz))
      call self% comm% barrier
      self% is_allocated = .true.
      if (VRBZ>0) &
         print '("[i]'//subnam// &
            ' allocated self% d_A("I7" = "I7"*"I7"*"I7")")', &
            sz, self%nxl, self%nyl, self%nzl
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dPart_allocate_arrays


   !> dPart destructor: base subroutine
   subroutine dealloc_dPart_base(self, VRBZ_)
   class(dPart), intent(INOUT) :: self
   integer, intent(in), optional :: VRBZ_
   !
   character(*), parameter :: subnam = '[dealloc_dPart_base]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      !! Nullify Pointers
      if(associated(self%d_A)) then
         deallocate(self%d_A)
         nullify(self%d_A)
      endif
      if(associated(self%h_A)) then
         deallocate(self%h_A)
         nullify(self%h_A)
      endif
      self% is_allocated = .false.
      self% is_init = .false.
      if(associated(self% comm)) nullify(self%comm)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dealloc_dPart_base


   !> dPart destructor
   subroutine dealloc_dPart(self)
   type(dPart), intent(INOUT) :: self
   !
   character(*), parameter :: subnam = '[dealloc_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call dealloc_dPart_base(self)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dealloc_dPart


   !> create cublasmp grid handle
   function stack_cublasmp_grid_handle(self, VRBZ_) result(grid_handle)
   type(c_ptr) :: grid_handle
   class(stack),      intent(IN) :: self
   integer, optional, intent(IN) :: VRBZ_
   !
   integer :: VRBZ
   character(*),parameter :: subnam='[stack_cublasmp_grid_handle]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      if (VRBZ>0) write(*,'(A,5(I7,1X))') '[i]'//subnam//': set_cublasMpGrid, ', &
                                self%numRowDevices, self%numColDevices
      grid_handle= set_cublasMpGrid(self%numRowDevices, self%numColDevices, VRBZ_)
      if (VRBZ>0) call show_cublasMpGrid(grid_handle)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_cublasmp_grid_handle


   !> create cusolvermp grid handle
   function stack_cusolvermp_grid_handle(self, VRBZ_) result(grid_handle)
   type(c_ptr) :: grid_handle
   class(stack),      intent(IN) :: self
   integer, optional, intent(IN) :: VRBZ_
   !
   integer                       :: VRBZ
   character(*),parameter    :: subnam='[stack_cusolvermp_grid_handle]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      if (VRBZ>0) write(*,'(A,5(I7,1X))') '[i]'//subnam//': set_cusolverMpGrid, ', &
                             self%numRowDevices, self%numColDevices
      grid_handle= set_cusolverMpGrid(self%numRowDevices, self%numColDevices, VRBZ_)
      if (VRBZ>0) call show_cusolverMpGrid(grid_handle)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_cusolvermp_grid_handle


   !> create cublasmp matrix handle
   function stack_cublasmp_matrix_handle(self, grid_handle, VRBZ_) result(matrix_handle)
   type(c_ptr) :: matrix_handle
   class(stack),      intent(IN) :: self
   type(c_ptr),       intent(IN) :: grid_handle
   integer, optional, intent(IN) :: VRBZ_
   !
   integer :: VRBZ
   character(*),parameter :: subnam='[stack_cublasmp_matrix_handle]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      if (VRBZ.ne.0) write(*,'(A,19(I7,1X))') &
                    '[i]'//subnam//': set_cublasMpMatrixDesc, ', &
                     self%M, self%N, self%MB, self%NB, &
                     self%numRowDevices, self%dev_row_idx
      matrix_handle= set_cublasMpMatrixDesc(grid_handle, &
                     self%M, self%N, self%MB, self%NB, &
                     self%numRowDevices, self%dev_row_idx, VRBZ_=VRBZ)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_cublasmp_matrix_handle


   !> create cusolvermp matrix handle
   function stack_cusolvermp_matrix_handle(self, grid_handle, VRBZ_) result(matrix_handle)
   type(c_ptr) :: matrix_handle
   class(stack),      intent(IN) :: self
   type(c_ptr),       intent(IN) :: grid_handle
   integer, optional, intent(IN) :: VRBZ_
   !
   integer :: VRBZ
   character(*),parameter :: subnam='[stack_cusolvermp_matrix_handle]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      if (VRBZ.ne.0) write(*,'(A,9(I7,1X))') &
                    '[i]'//subnam//': set_cusolverMpMatrixDesc, ', &
                     self%M, self%N, self%MB, self%NB, &
                     self%numRowDevices, self%dev_row_idx
      matrix_handle= set_cusolverMpMatrixDesc(grid_handle, &
                     self%M, self%N, self%MB, self%NB, &
                     self%numRowDevices, self%dev_row_idx, VRBZ_=VRBZ)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_cusolvermp_matrix_handle


   !> Allocating data arrays in Stack
   subroutine stack_allocate_arrays(self, VRBZ_)
   class(stack),   intent(INOUT) :: self
   integer, optional, intent(IN) :: VRBZ_
   !
   character(*), parameter :: subnam = '[stack_allocate_arrays]'
   integer :: sz, VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      sz= self%ML * self%NL
      if (VRBZ>0) &
         print '("[i]'//subnam//': allocating ("I7"x"I7")")', self%ML, self%NL
      if (associated(self% d_A)) deallocate(self% d_A)
      if (associated(self% h_A)) deallocate(self% h_A)
      if (sz > 0) then
         allocate(self%d_A(sz))
         allocate(self%h_A(sz))
         self% d_A2(1:self%ML,1:self%NL)=> self% d_A
         self% h_A2(1:self%ML,1:self%NL)=> self% h_A
         if (VRBZ>0) &
            print '("[i]'//subnam//': allocating d_A2("I7"x"I7") array")', self%ML, self%NL
      else
         allocate(self%d_A(1))
         allocate(self%h_A(1))
         self% d_A2(1:1,1:1)=> self% d_A
         self% h_A2(1:1,1:1)=> self% h_A
         if (VRBZ>0) &
            print '("[i]'//subnam//': local size = 0, allocating d_A2(1x1)")'
      endif
      call self% comm% barrier()
      self% is_allocated = .true.
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_allocate_arrays


   !> Stack deallocation helper function
   subroutine dealloc_stack_base(self, VRBZ_)
   class(stack), intent(INOUT)   :: self
   integer, intent(in), optional :: VRBZ_
   !
   integer                       :: VRBZ
   character(*),parameter    :: subnam='[stack_dealloc_stack_base]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
      call dealloc_dPart_base(self)
      nullify(self%d_A2)
      nullify(self%h_A2)
      nullify(self%comm)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine dealloc_stack_base


   !> Stack destructor
   subroutine dealloc_stack(self)
   type(stack), intent(INOUT) :: self
#  ifdef VERBOSE4
      print '("[+][dealloc_stack] entry")'
#  endif
      call dealloc_stack_base(self)
#  ifdef VERBOSE4
      print '("[+][dealloc_stack] exit")'
#  endif
   end subroutine dealloc_stack


!
!  [UNFOLDING]
!

   !> hStack constructor helper function: complete construction of the hStack
   subroutine hstack_dPart(H, N_, VRBZ_)
   type(stack),    intent(INOUT) :: H
   integer, optional, intent(IN) :: N_  !< use this for arbitrary horizontal size
   integer, optional, intent(IN) :: VRBZ_
   !
   integer :: sh(2),VRBZ,nx,ny,nz,nxb,nyb,nzb,M_B,N_B
   integer :: n_dev_row,n_dev_col,i_row,i_col
   character(*),parameter :: subnam='[hstack_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ=VRBZ_
      if(VRBZ>0) call print_info(H, label_=subnam//": H upon entry:")

      nx = H% nx; nxb = H% nxb
      ny = H% ny; nyb = H% nyb
      nz = H% nz; nzb = H% nzb
      select case(H% part_axis)
      case(1) ! x-part
         M_B = nxb;  n_dev_row = H% nranks; i_row = H% rank
         N_B = ny;   n_dev_col = 1;         i_col = 0
      case(2) ! y-part
         M_B = nx;   n_dev_row = 1;         i_row = 0
         N_B = nyb;  n_dev_col = H% nranks; i_col = H% rank
         if (ny /= nyb*H% nranks) error stop '[!]'//subnam//": in hStack, ny /= nyb*nranks"
      case(3) ! z-part
         M_B = nx;      n_dev_row = 1;         i_row = 0
         N_B = ny*nzb;  n_dev_col = H% nranks; i_col = H% rank
      end select

      H% numRowDevices = n_dev_row;    H% dev_row_idx = i_row
      H% numColDevices = n_dev_col;    H% dev_col_idx = i_col

      H% M = nx
      H% N = ny*nz; if (present(N_)) H% N = N_
      H% MB = M_B
      H% NB = N_B
      H% ML = numroc(H% M, H% MB, H% rank, H% numRowDevices)
      H% NL = numroc(H% N, H% NB, H% rank, H% numColDevices)

      call stack_allocate_arrays(H, VRBZ_)
      H% d_A2(1:H%ML,1:H%NL) => H% d_A
      H% h_A2(1:H%ML,1:H%NL) => H% h_A

      H% is_init = .true.
      if(VRBZ>0) call print_info(H, label_=subnam//": H upon exit:")
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine hstack_dPart


   !> vStack constructor helper function: complete construction of the vStack
   subroutine vstack_dPart(V, M_, VRBZ_)
   type(stack),    intent(INOUT) :: V
   integer, optional, intent(IN) :: M_   !< use this for arbitrary vertical size
   integer, optional, intent(IN) :: VRBZ_
   !
   integer :: sh(2),VRBZ,nx,ny,nz,nxb,nyb,nzb,M_B,N_B
   integer :: n_dev_row,n_dev_col,i_row,i_col
   character(*),parameter :: subnam='[vstack_dPart]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if(present(VRBZ_)) VRBZ=VRBZ_
      if(VRBZ>0) call print_info(V, label_=subnam//": V upon entry:")

      nx = V% nx; nxb = V% nxb
      ny = V% ny; nyb = V% nyb
      nz = V% nz; nzb = V% nzb
      select case(V% part_axis)
      case(1) ! x-part
         M_B = nxb;  n_dev_row = V% nranks; i_row = V% rank
         N_B = ny;   n_dev_col = 1;         i_col = 0
         if (nx /= nxb*V% nranks) error stop '[!]'//subnam//": in vStack, nx /= nxb*nranks"
      case(2) ! y-part
         M_B = nx;   n_dev_row = 1;         i_row = 0
         N_B = nyb;  n_dev_col = V% nranks; i_col = V% rank
      case(3) ! z-part
         M_B = nx*nzb;  n_dev_row = V% nranks; i_row = V% rank
         N_B = ny;      n_dev_col = 1;         i_col = 0
      end select

      V% numRowDevices = n_dev_row;    V% dev_row_idx = i_row
      V% numColDevices = n_dev_col;    V% dev_col_idx = i_col

      V% M = nx*nz; if (present(M_)) V% M = M_
      V% N = ny
      V% MB = M_B
      V% NB = N_B
      V% ML = numroc(V% M, V% MB, V% rank, V% numRowDevices)
      V% NL = numroc(V% N, V% NB, V% rank, V% numColDevices)

      call stack_allocate_arrays(V, VRBZ_)
      V% d_A2(1:V%ML,1:V%NL) => V% d_A
      V% h_A2(1:V%ML,1:V%NL) => V% h_A

      V% is_init = .true.
      if(VRBZ>0) call print_info(V, label_=subnam//": V upon exit:")
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine vstack_dPart


!
!  [TRANSPOSE]
!
   !> Transposes a stack, slice-by-slice
   !! TODO: valid only for z-part
   subroutine stack_transpose_slices(A, VRBZ_, CUDA_BLOCK_)
   class(stack),   intent(INOUT) :: A
   integer, optional, intent(IN) :: VRBZ_, CUDA_BLOCK_
   !
   character(*),parameter :: subnam='[stack_transpose_slices]'
   integer  :: VRBZ, CUDA_BLOCK
   integer  :: slice, NN, i0, i1, nrow, ncol
   double precision, pointer, contiguous, device :: B(:)
   logical  :: A_is_hstack
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (.not.A%is_init) error stop subnam//': A is not initialized'
      if (.not.associated(A%d_A)) error stop '[!]'//subnam//': A%d_A is not associated'
      if (.not.A% well_stacked()) error stop '[!]'//subnam//': A is not well stacked'

      VRBZ       = 0;  if(present(VRBZ_)) VRBZ=VRBZ_
      CUDA_BLOCK = 64; if(present(CUDA_BLOCK_)) CUDA_BLOCK=CUDA_BLOCK_
      if (VRBZ > 0) then
         print "('Before transpose_slice:')"
         print "('nxb=',i8,'| nyb=',i8,'|nzb=',i8)", A%nxb, A%nyb, A%nzb
         print "('nx=',i8,'| ny=',i8,'|nz=',i8)", A%nx, A%ny, A%nz
         print "('MB, NB, M, N = ',4(I4,1X))", A%MB, A%NB, A%M, A%N
         print "('CUDA_BLOCK =',i8)", CUDA_BLOCK
      endif
      ! it is important to determine the unfolding before swapping dimensions
      A_is_hstack = A% is_hstack()
      NN = A%nx;  A%nx =A%ny;  A%ny =NN; ! nx <-> ny
      NN = A%nxb; A%nxb=A%nyb; A%nyb=NN; ! nxb <-> nyb
      NN = A%nxl; A%nxl=A%nyl; A%nyl=NN; ! nxl <-> nyl
      if(A_is_hstack) then
         A%MB = A%nx
         A%NB = A%ny * A%nzb
         A%M  = A%nx
         A%N  = A%ny * A%nz
      else
         A%MB = A%nx * A%nzb
         A%NB = A%ny
         A%M  = A%nx * A%nz
         A%N  = A%ny
      endif

      A%ML = numroc(A%M, A%MB, A%dev_row_idx, A%numRowDevices)
      A%NL = numroc(A%N, A%NB, A%dev_col_idx, A%numColDevices)
      nrow=A%nxb; ncol=A%nyb

      if (A%ML * A%NL > 0) then
         allocate(B(A%ML*A%NL))
         B = A%d_A
         A%d_A2(1:A%ML, 1:A%NL) => A%d_A
         if (VRBZ > 0) then
            print '("After transpose_slice:")'
            print '("nxb, nyb, nzb   = ",3(I7,1X))', A%nxb, A%nyb, A%nzb
            print '("A%nx, A%ny = ",3(I7,1X))', A%nx, A%ny
            print '("MB, NB, M, N = ",4(I7,1X))', A%MB, A%NB, A%M, A%N
         endif
   !      CUDA_BLOCK = min(CUDA_BLOCK,nrow) ![TODO] Make sure this fix works for all scenarios
         NN = A%nxl * A%nyl
         do slice=0,A%nzl-1
             i0 = slice*NN + 1
             i1 =  i0 + NN - 1
             call transpose_device_array_vect(B(i0:i1), A%d_A(i0:i1), &
                  A%nyl, A%nxl, CUDA_BLOCK)
         end do
         deallocate(B)
      endif
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_transpose_slices


   !> Transposes a shallow stack (matrix), swaps the unfolding
   !! TODO: valid only for z-part
   subroutine stack_transpose_matrix(A, VRBZ_, CUDA_BLOCK_)
   class(stack), intent(INOUT) :: A
   integer,intent(in),optional :: VRBZ_, CUDA_BLOCK_
   !
   character(*),parameter :: subnam='[stack_transpose_matrix]'
   integer  :: VRBZ, CUDA_BLOCK
   integer  :: NN, i0, i1, nrow, ncol
   double precision, pointer, contiguous, device :: B(:)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif

      if (.not.A%is_init) error stop subnam//': A is not initialized'
      if (.not.associated(A%d_A)) error stop '[!]'//subnam//': A%d_A is not associated'
      if (A% nzb /= 1) error stop '[!]'//subnam//': A is not a shallow stack/matrix'

      VRBZ       = 0;  if(present(VRBZ_)) VRBZ=VRBZ_
      CUDA_BLOCK = 64; if(present(CUDA_BLOCK_)) CUDA_BLOCK=CUDA_BLOCK_
      if (VRBZ > 0) then
         print "('Before transpose_slice:')"
         print "('nxb=',i8,'| nyb=',i8,'|nzb=',i8)", A%nxb, A%nyb, A%nzb
         print "('nx=',i8,'| ny=',i8,'|nz=',i8)", A%nx, A%ny, A%nz
         print "('MB, NB, M, N = ',4(I4,1X))", A%MB, A%NB, A%M, A%N
         print "('CUDA_BLOCK =',i8)", CUDA_BLOCK
      endif
      ! it is important to determine the unfolding before swapping dimensions
      NN = A%nx;  A%nx =A%ny;  A%ny =NN; ! nx <-> ny
      NN = A%nxb; A%nxb=A%nyb; A%nyb=NN; ! nxb <-> nyb
      NN = A%nxl; A%nxl=A%nyl; A%nyl=NN; ! nxl <-> nyl
      NN = A%M;   A%M =A%N;    A%N  =NN; ! M <-> N
      NN = A%MB;  A%MB=A%NB;   A%NB =NN; ! MB <-> NB
      NN = A%ML;  A%ML=A%NL;   A%NL =NN; ! ML <-> NL

      nrow = A% numRowDevices
      ncol = A% numColDevices
      A% numRowDevices = ncol
      A% numColDevices = nrow

      nrow = A% dev_row_idx
      ncol = A% dev_col_idx
      A% dev_row_idx = ncol
      A% dev_col_idx = nrow

      if (A%ML * A%NL > 0) then
         allocate(B(A%ML*A%NL))
         B = A%d_A
         A%d_A2(1:A%ML, 1:A%NL) => A%d_A
         if (VRBZ > 0) then
            print '("After transpose_slice:")'
            print '("nxb, nyb, nzb   = ",3(I7,1X))', A%nxb, A%nyb, A%nzb
            print '("A%nx, A%ny = ",3(I7,1X))', A%nx, A%ny
            print '("MB, NB, M, N = ",4(I7,1X))', A%MB, A%NB, A%M, A%N
         endif
         call transpose_device_array_vect(B, A%d_A, A%NL, A%ML, CUDA_BLOCK)
         deallocate(B)
      endif
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_transpose_matrix


!
!  [ADD/SUBTRACT]
!
   !> Stacks addition operation
   function add_stacks(A, B) result(C)
   type(stack),  intent(in) :: A, B
   type(stack)              :: C
   !
   integer :: sz, VRBZ, info
   type(c_ptr) :: grid, descrB, descrC
   character(*),parameter  :: subnam='[add_stacks]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0
      if (.not.associated(A%d_A)) error stop '[!]'//subnam//': A% d_A is not associated'
      if (.not.associated(B%d_A)) error stop '[!]'//subnam//': B% d_A is not associated'
      if (A%M.ne.B%M .or. A%N.ne.B%N) error stop '[!]'//subnam//': A and B have different sizes'
      if (A% is_hstack().and.B% is_hstack()) then
         if (A%NL .ne. B%NL) error stop '[!]'//subnam//': A and B have different NL'
      elseif (A% is_vstack().and.B% is_vstack()) then
         if (A%ML .ne. B%ML) error stop '[!]'//subnam//': A and B have different ML'
      else
         error stop '[!]'//subnam//': A and B must be both either hstack or vstack'
      endif
      C = A
      grid = C% cublasmp_grid_handle(VRBZ_=VRBZ)
      descrB = B% cublasmp_matrix_handle(grid, VRBZ_=VRBZ)
      info= CaddA(B%M, B%N, descrB, B%d_A, descrB, C%d_A, VRBZ_=VRBZ)
      call stop_if_error(info, '[!]'//subnam//': err=')
      call destroy_cublasMpMatrixDesc(descrB, VRBZ_=VRBZ)
      call destroy_cublasMpGrid(grid, VRBZ_=VRBZ)

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function add_stacks


   !> Stacks subtraction
   function subtract_stacks(A, B) result(C) !
   type(stack),  intent(in)    :: A, B
   type(stack)                 :: C
   !
   integer :: sz, VRBZ, info
   type(c_ptr) :: grid, descrB, descrC
   character(*),parameter  :: subnam='[subtract_stacks]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0
      if (.not.associated(A%d_A)) error stop '[!]'//subnam//': A% d_A is not associated'
      if (.not.associated(B%d_A)) error stop '[!]'//subnam//': B% d_A is not associated'
      if (A%M.ne.B%M .or. A%N.ne.B%N) error stop '[!]'//subnam//': A and B have different sizes'
      if (A% is_hstack().and.B% is_hstack()) then
         if (A%NL .ne. B%NL) error stop '[!]'//subnam//': A and B have different NL'
      elseif (A% is_vstack().and.B% is_vstack()) then
         if (A%ML .ne. B%ML) error stop '[!]'//subnam//': A and B have different ML'
      else
         error stop '[!]'//subnam//': A and B must be both either hstack or vstack'
      endif
      C = A
      grid = C% cublasmp_grid_handle(VRBZ_=VRBZ)
      descrB = B% cublasmp_matrix_handle(grid, VRBZ_=VRBZ)
      info= CsubA(B%M, B%N, descrB, B%d_A, descrB, C%d_A, VRBZ_=VRBZ)
      call stop_if_error(info, '[!]'//subnam//': err=')
      call destroy_cublasMpMatrixDesc(descrB, VRBZ_=VRBZ)
      call destroy_cublasMpGrid(grid, VRBZ_=VRBZ)

      !! TODO: only if allocated
      !sz= size(A% h_A)
      !allocate(C% h_A(sz))
      !C% h_A= A% h_A - B% h_A

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function subtract_stacks


   !> Frobenius norm
   double precision function stack_nrm2(A, VRBZ_) result(nrm)
   use cublas_aux, only: cublas_nrm2
   class(stack), intent(IN) :: A
   integer, intent(IN), optional :: VRBZ_
   !
   character(*),parameter  :: subnam='[stack][stack_nrm2]'
   double precision :: nrm_local
   integer :: VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 1; if (present(VRBZ_)) VRBZ= VRBZ_
      nrm_local = cublas_nrm2(A%d_A)**2
      nrm  = dsqrt(A%comm%comm%allreduce(nrm_local))
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_nrm2


!
!  [MULT]
!
   !> Multiplication operator
   function AxB_stacks(A, B, ops_, VRBZ_) result(C)
   type(stack)                            :: C
   type(stack),             intent(INOUT) :: A, B
   character(len=2), optional, intent(IN) :: ops_ !< 'nn', 'nt', 'tn', 'tt'
   integer,          optional, intent(IN) :: VRBZ_
   !
   integer :: VRBZ, info, rank, nranks, rank_row, rank_col
   integer :: n_dev_row, n_dev_col, M, N, K, ng_C(3), nb_C(3), MorN
   integer :: MBA,MBB,MBC,NBA,NBB,NBC
   character(len=2) :: ops
   character :: vh
   character(*),parameter :: subnam='[AxB_stacks]'
   type(c_ptr) :: grid, descrA, descrB, descrC
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (.not.associated(A%d_A)) error stop '[!]'//subnam//': A% d_A is not associated'
      if (.not.associated(B%d_A)) error stop '[!]'//subnam//': B% d_A is not associated'
      VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
      ops='nn'; if (present(ops_)) ops=ops_

      ! Create matrix descriptors suitable for multiplication
      rank = A% rank; nranks = A% nranks
      vh = 'h'; if (A% is_vstack) vh = 'v'
      if (A%is_vstack() .and. B%is_vstack()) then
         vh = 'v'
      elseif (A%is_hstack() .and. B%is_hstack()) then
         vh = 'h'
      else
         error stop '[!]'//subnam// &
            ': A and B should both be either vstacks or hstacks'
      endif

      ! Initialize global matrix dimensions & do compatibility checks
      select case (ops)
      case ('nn')
         M = A%M; N = B%N; K = A%N; MorN = M
         if (B%M /= K) error stop '[!]'//subnam//': A%N =/= B%M'
      case ('nt')
         M = A%M; N = B%M; K = A%N; MorN = M
         if (B%N /= K) error stop '[!]'//subnam//': A%N =/= B%N'
      case ('tn')
         M = A%N; N = B%N; K = A%M; MorN = N
         if (B%M /= K) error stop '[!]'//subnam//': A%M =/= B%M'
      case ('tt')
         M = A%N; N = B%M; K = A%M; MorN = N
         if (B%N /= K) error stop '[!]'//subnam//': A%M =/= B%N'
      case default
         error stop '[!]'//subnam//': unknown operation ops_='//ops
      end select

      MBA = A% MB;  NBA = A% NB
      MBB = B% MB;  NBB = B% NB

      if (vh.eq.'v') then
         !! TODO: this only works for z-part!!!
         n_dev_row = A% numRowDevices
         n_dev_col = 1
         rank_row  = rank
         rank_col  = 0

         select case (ops)
         case ('nn')
            ng_C = [A% nx,  B% ny,  A% nz]
            nb_C = [A% nxb, B% nyb, A% nzb]
            NBA = MBB
            NBB = MBA
            MBC = MBA;  NBC = NBB
         case ('nt')
            ng_C = [A% nx, B% nx*B% nz,A% nz]
            nb_C = [A% nx, B% nx*B% nz,A% nzb]
            NBA = MBA
            NBB = MBA
            MBC = MBA;  NBC = MBB
         case ('tn')
            ng_C = [A% NB, B% N,(A%N - 1)/A% NB + 1]
            nb_C = [A% NB, B% N, 1]
            NBA = MBA
            NBB = MBA
            MBC = NBA;  NBC = NBB
         end select

      else !! vh.eq.'h'
         ! Create matrix descriptors suitable for multiplication
         n_dev_row = 1
         n_dev_col = A% numColDevices
         rank_row  = 0
         rank_col  = rank

         select case (ops)
         case ('nn')
            ng_C = [A% nx,  B% ny,  B% nz]
            nb_C = [A% nxb, B% nyb, B% nzb]
            MBA = NBB
            MBB = NBA
            MBC = MBA;  NBC = NBB
         case ('nt')
            MBB = (B% M - 1)/nranks + 1
            ng_C = [A% M, MBB, min(B%M, nranks)]
            nb_C = [A% M, MBB, 1]
            MBA = NBA
            MBC = MBA;  NBC = MBB
         case ('tn')
            ng_C = [A% ny*A% nz, B% ny, B% nz]
            nb_C = [A% ny*A% nz, B% ny, B% nzb]
            MBA = NBA
            MBB = NBA
            MBC = NBA;  NBC = NBB
         end select
      endif

      ! Create stack C
      C = zeros_stack(A%comm, ng_C, nb_C, vh, morn_=MorN, VRBZ_=VRBZ)
      grid = A% cublasMp_grid_handle(VRBZ_=VRBZ)

      ! Update matrix descriptors for multiplication
      if (VRBZ>2) print '("[i]["I3"/"I3"]'//subnam &
         //': initializing descrA:")', rank, nranks
      A% MB = MBA; A% NB = NBA
      descrA = A% cublasmp_matrix_handle(grid, VRBZ_=VRBZ)

      if (VRBZ>2) print '("[i]["I3"/"I3"]'//subnam &
         //': initializing descrB:")', rank, nranks
      B% MB = MBB; B% NB = NBB
      descrB = B% cublasmp_matrix_handle(grid, VRBZ_=VRBZ)

      if (VRBZ>2) print '("[i]["I3"/"I3"]'//subnam &
         //': initializing descrC:")', rank, nranks
      C% MB = MBC; C% NB = NBC
      descrC = C% cublasmp_matrix_handle(grid, VRBZ_=VRBZ)
      select case(ops)
      case ('nn')
         info = AxB(M,  N, K, descrA, A%d_A, descrB, B%d_A, descrC, C%d_A, VRBZ_=VRBZ)
      case ('nt')
         info = AxBT(M, N, K, descrA, A%d_A, descrB, B%d_A, descrC, C%d_A, VRBZ_=VRBZ)
      case ('tn')
         info = ATxB(M, N, K, descrA, A%d_A, descrB, B%d_A, descrC, C%d_A, VRBZ_=VRBZ)
      case ('tt')
         info = ATxBT(M,N, K, descrA, A%d_A, descrB, B%d_A, descrC, C%d_A, VRBZ_=VRBZ)
      end select
      call stop_if_error(info, '[!]'//subnam//': err=')

      call destroy_cublasMpMatrixDesc(descrA, VRBZ_=VRBZ)
      call destroy_cublasMpMatrixDesc(descrB, VRBZ_=VRBZ)
      call destroy_cublasMpMatrixDesc(descrC, VRBZ_=VRBZ)
      call destroy_cublasMpGrid(grid, VRBZ_=VRBZ)

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function AxB_stacks


   !> Multiplying stacks: C = A @ B^T
   function AxBT_stacks(A, B, VRBZ_) result(C) !
   type(stack)                   :: C
   type(stack),    intent(INOUT) :: A, B
   integer, optional, intent(IN) :: VRBZ_
   !
   integer   :: VRBZ
   character(*),parameter :: subnam='[stack][AxBT_stacks]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
      C = AxB_stacks(A, B, ops_='nt', VRBZ_=VRBZ)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function AxBT_stacks


   !> Multiplying stacks: C = A^T @ B
   function ATxB_stacks(A, B, VRBZ_) result(C) !
   type(stack)                   :: C
   type(stack),    intent(INOUT) :: A, B
   integer, optional, intent(IN) :: VRBZ_
   !
   integer   :: VRBZ
   character(*),parameter :: subnam='[stack][ATxB_stacks]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
      C = AxB_stacks(A, B, ops_='tn', VRBZ_=VRBZ)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function ATxB_stacks

!
!  [EYE]
!
   !> Stack that represents a unit diagonal matrix (could be rectangular)
   function eye_stack_3d(comm, ng, nb, stacking, part_axis_, VRBZ_) result(eye)
   type(stack):: eye
   type(super_comm), target, intent(IN):: comm !< supercommunicator, used for rank & size
   integer,                  intent(IN):: ng(3), nb(3) !< global dimensions and block sizes
   character,      optional, intent(IN):: stacking   !< 'h' by default (horizontal)
   integer,        optional, intent(IN):: part_axis_  !< partition axis (default=3, z-axis)
   integer,        optional, intent(IN):: VRBZ_
   !
   character(*),parameter :: subnam='[eye_stack_3d]'
   integer :: part_axis, VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if (present(VRBZ_)) VRBZ= VRBZ_
      part_axis= 3; if (present(part_axis_)) part_axis= part_axis_
      eye = zeros_stack(comm, ng, nb, stacking, part_axis_=part_axis, VRBZ_=VRBZ)

      call create_eye_matrix(eye% d_A, eye%M, eye%N, eye%MB, eye%NB, &
                             eye% dev_row_idx, eye% dev_col_idx, &
                             eye% ML, eye% NL, &
                             eye% numRowDevices, eye% numColDevices)
      eye% d_A2(1:eye%ML,1:eye%NL)=> eye% d_A

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function eye_stack_3d


   !> Stack that represents a unit diagonal matrix (could be rectangular)
   function eye_stack_2d(comm, M, N, stacking, MB_, part_axis_, VRBZ_) result(eye)
   type(stack):: eye
   type(super_comm), target,   intent(IN):: comm
   integer,                    intent(IN):: M, N     !< size of the matrix (M rows, N columns)
   character,                  intent(IN):: stacking !< 'h'/'H' or 'v'/'V' ['h']
   integer,          optional, intent(IN):: MB_      !< block size (default: determined by nranks)
   integer,          optional, intent(IN):: part_axis_ !< partition axis (default=3, z-axis)
   integer,          optional, intent(IN):: VRBZ_
   !
   character(*),parameter :: subnam='[eye_stack_2d]'
   integer :: part_axis, VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if (present(VRBZ_)) VRBZ= VRBZ_
      part_axis= 3; if (present(part_axis_)) part_axis= part_axis_
      eye = zeros_stack(comm, M, N, stacking, MB_, part_axis_, VRBZ_)

      call create_eye_matrix(eye% d_A, eye%M, eye%N, eye%MB, eye%NB, &
                             eye% dev_row_idx, eye% dev_col_idx, &
                             eye% ML, eye% NL, &
                             eye% numRowDevices, eye% numColDevices)
      eye% d_A2(1:eye%ML,1:eye%NL)=> eye% d_A

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function eye_stack_2d


   !> Unit diagonal matrix with dimensions like in the other stack
   function eye_like_stack(other, VRBZ_) result(eye)
   type(stack), intent(in) :: other
   type(stack):: eye
   integer, intent(in), optional :: VRBZ_
   !
   character(*),parameter    :: subnam='[eye_like_stack]'
   integer :: ng(3), nb(3), VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
!print '("[i]'//subnam//': VRBZ_ present? ",L1)', present(VRBZ_)
!if(present(VRBZ_)) print '("[i]'//subnam//': VRBZ_=",I1)', VRBZ_
      VRBZ=0; if (present(VRBZ_)) VRBZ= VRBZ_
      if (.not.other% well_stacked()) error stop '[!]'//subnam//': other not well stacked'
      ng = [other% nx, other% ny, other% nz]
      nb = [other% nxb, other% nyb, other% nzb]
      if (other% is_hstack()) then
         eye = eye_stack(other%comm, ng, nb, 'h', part_axis_=other%part_axis, VRBZ_=VRBZ)
      else if (other% is_vstack()) then
         eye = eye_stack(other%comm, ng, nb, 'v', part_axis_=other%part_axis, VRBZ_=VRBZ)
      endif
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function eye_like_stack


   !> Add array d_d to the (n-th sub/super) diagonal of the stack
   subroutine stack_add_to_diag(this, d_diag, ndiag, VRBZ_)
   class(stack),   intent(INOUT):: this
   double precision, device, pointer, intent(IN):: d_diag(:) !< 1D array of the diagonal
   integer,           intent(IN):: ndiag !< diagonal number: 0 - main, -1 - subdiagonal etc.
   integer, optional, intent(IN):: VRBZ_
   !
   character(*),parameter :: subnam='[stack_add_to_diag]'
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      call add_matrix_diag(this% d_A, d_diag, ndiag, &
                           this%M, this%N, &
                           this%MB,this%NB, &
                           this%ML,this%NL, &
                           this%dev_row_idx, this%dev_col_idx, &
                           this%numRowDevices, this%numColDevices)

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_add_to_diag
!
!  [QR]
!
   !> Flips stacking between vertical <-> horizontal
   subroutine switch_stacking(A, VRBZ_, CUDA_BLOCK_)
   class(stack), intent(inout)  :: A
   integer, intent(in), optional :: VRBZ_, CUDA_BLOCK_
   !
   double precision, pointer, device :: d_A_new(:)
   character(*),parameter :: subnam='[switch_stacking]'
   integer :: sz, VRBZ, CUDA_BLOCK
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif

      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
      CUDA_BLOCK=256; if(present(CUDA_BLOCK_)) CUDA_BLOCK=CUDA_BLOCK_
      if (.not.A% well_stacked()) error stop '[!]'//subnam//': A is not well stacked'

      sz= size(A% d_A)
      allocate(d_A_new(sz))
      if (A% is_hstack()) then
         call transpose_device_flat3d_132(A%d_A, d_A_new, A%nxl, A%nyl, A%nzl, CUDA_BLOCK)
         A% M         = A% nx * A% nz
         A% N         = A% ny
         A% MB        = A% nx * A% nzb
         A% NB        = A% MB
         A% numRowDevices =  A% nranks
         A% numColDevices =  1
         A% dev_row_idx   = A% rank
         A% dev_col_idx   = 0
      else if (A% is_vstack()) then
         call transpose_device_flat3d_132(A%d_A, d_A_new, A%nxl, A%nzl, A%nyl, CUDA_BLOCK)
         A% M         = A% nx
         A% N         = A% ny * A% nz
         A% NB        = A% ny * A% nzb
         A% MB        = A% NB
         A% numRowDevices =  1
         A% numColDevices =  A% nranks
         A% dev_row_idx   = 0
         A% dev_col_idx   = A% rank
      else
         error stop "[!]"//subnam//": stack is neither vStack nor hStack"
      endif

      deallocate(A% d_A)
      A%d_A(1:sz) => d_A_new
      A% ML = numroc(A% M, A% MB, A% dev_row_idx, A% numRowDevices)
      A% NL = numroc(A% N, A% NB, A% dev_col_idx, A% numColDevices)
      A%d_A2(1:A%ML, 1:A%NL) => A%d_A

      if (VRBZ>0) then
         print "('nxb= ',i8,'| nyb= ',i8,'|nzb= ',i8)", A%nxb, A%nyb, A%nzb
         print "('nx=',i8,'| ny=',i8,'|nz=',i8)", A%nx, A%ny, A%nz
         print "('MB=',i8,'| NB=',i8,'|M=',i8,'|N=',i8,'|ML= ',i8,'|NL=',i8)", &
                A%MB, A%NB, A%M, A%N, A%ML, A%NL
      endif

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine switch_stacking


   !> Converts stack to hStack
   subroutine stack_to_hstack(A, VRBZ_) !
   class(stack), intent(inout)  :: A
   integer, intent(in), optional :: VRBZ_
   !
   character(*),parameter :: subnam='[stack][to_hstack]'
   integer :: VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
      if (A% is_vstack()) call A% switch_stacking(VRBZ_=VRBZ)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_to_hstack


   !> Converts stack to hStack
   subroutine stack_to_vstack(A, VRBZ_) !
   class(stack), intent(inout)  :: A
   integer, intent(in), optional :: VRBZ_
   !
   character(*),parameter :: subnam='[stack][to_vstack]'
   integer :: VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
      if (A% is_hstack()) call A% switch_stacking(VRBZ_=VRBZ)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_to_vstack


   !> Redistributes stack to shallow (nz=nranks) vStack
   subroutine stack_redist_to_vstack(A, VRBZ_) !
   class(stack), intent(inout)  :: A
   integer, intent(in), optional :: VRBZ_
   !
   character(*),parameter :: subnam='[stack_redist_to_vstack]'
   double precision, device, pointer :: d_B(:)
   type(c_ptr) :: grid, descrA, grid_v, descrA_v
   integer :: VRBZ, M,N, MB,NB, ML,NL, err
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_

      ! check if A is already a vstack; if so, do nothing
      if (A% is_vstack()) return

      ! record the previous grid and matrix description
      grid   = A% cublasMp_grid_handle(VRBZ_=VRBZ)
      descrA = A% cublasmp_matrix_handle(grid, VRBZ_=VRBZ)

      ! compute the new distribution
      M = A% M;                                  N = A% N
      MB = (M - 1)/A% nranks + 1;                NB = N
      ML = numroc(M, MB, A% rank, A% nranks);    NL = N

      ! reset parameters of the stack A
      A% numRowDevices = A% nranks
      A% numColDevices = 1
      A% dev_row_idx = A% rank
      A% dev_col_idx = 0
      A% MB = MB; A% NB = NB
      A% ML = ML; A% NL = NL
      allocate(d_B(ML*NL))

      ! grid and matrix handles for the new vertical grid
      grid_v   = A% cublasMp_grid_handle(VRBZ_=VRBZ)
      descrA_v = A% cublasmp_matrix_handle(grid_v, VRBZ_=VRBZ)

      ! redistribute
      err= redistribute_rectMatrix(M, N, descrA, A% d_A, descrA_v, d_B, VRBZ_)
      if (err /= 0) then
         print '("[!]'//subnam//': in redistribute_rectMatrix, err=",I7)', err
         error stop
      endif
      deallocate(A% d_A)
      A% d_A(1:ML*NL)=> d_B
      A% d_A2(1:ML,1:NL)=> d_B
      nullify(d_B)

      ! update dpart
      A% nx = MB;                 A% nxb = MB
      A% ny = NB;                 A% nyb = NB
      A% nz = min(A% nranks, MB); A% nzb = 1

      A% nxl = MB
      A% nyl = NB
      A% nzl = numroc(A% nz, 1, A% rank, A% nranks)

      ! destroy the descriptors
      call destroy_cublasMpMatrixDesc(descrA, VRBZ_=VRBZ)
      call destroy_cublasMpMatrixDesc(descrA_v, VRBZ_=VRBZ)
      call destroy_cublasMpGrid(grid, VRBZ_=VRBZ)
      call destroy_cublasMpGrid(grid_v, VRBZ_=VRBZ)

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_redist_to_vstack


   !> Redistributes stack to shallow (nz=nranks) hStack
   subroutine stack_redist_to_hstack(A, VRBZ_) !
   class(stack), intent(inout)  :: A
   integer, intent(in), optional :: VRBZ_
   !
   character(*),parameter :: subnam='[stack_redist_to_hstack]'
   double precision, device, pointer :: d_B(:)
   type(c_ptr) :: grid, descrA, grid_h, descrA_h
   integer :: VRBZ, M,N, MB,NB, ML,NL, err
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_

      ! check if A is already an hstack; if so, do nothing
      if (A% is_hstack()) return

      ! record the previous grid and matrix description
      grid   = A% cublasMp_grid_handle(VRBZ_=VRBZ)
      descrA = A% cublasmp_matrix_handle(grid, VRBZ_=VRBZ)

      ! compute the new distribution
      M = A% M;  N = A% N
      MB = M;    NB = (N - 1)/A% nranks + 1
 ;    ML = M;    NL = numroc(N, NB, A% rank, A% nranks)

      ! reset parameters of the stack A
      A% numRowDevices = 1
      A% numColDevices = A% nranks
      A% dev_row_idx = 0
      A% dev_col_idx = A% rank
      A% MB = MB; A% NB = NB
      A% ML = ML; A% NL = NL
      allocate(d_B(ML*NL))

      ! grid and matrix handles for the new horizontal grid
      grid_h   = A% cublasMp_grid_handle(VRBZ_=VRBZ)
      descrA_h = A% cublasmp_matrix_handle(grid_h, VRBZ_=VRBZ)

      ! redistribute
      err= redistribute_rectMatrix(M, N, descrA, A% d_A, descrA_h, d_B, VRBZ_)
      if (err /= 0) then
         print '("[!]'//subnam//': in redistribute_rectMatrix, err=",I7)', err
         error stop
      endif
      deallocate(A% d_A)
      A% d_A(1:ML*NL)=> d_B
      A% d_A2(1:ML,1:NL)=> d_B
      nullify(d_B)

      ! update dpart
      A% nx = MB;                 A% nxb = MB
      A% ny = NB;                 A% nyb = NB
      A% nz = min(A% nranks, NB); A% nzb = 1

      A% nxl = MB
      A% nyl = NB
      A% nzl = numroc(A% nz, 1, A% rank, A% nranks)

      ! destroy the descriptors
      call destroy_cublasMpMatrixDesc(descrA, VRBZ_=VRBZ)
      call destroy_cublasMpMatrixDesc(descrA_h, VRBZ_=VRBZ)
      call destroy_cublasMpGrid(grid, VRBZ_=VRBZ)
      call destroy_cublasMpGrid(grid_h, VRBZ_=VRBZ)

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_redist_to_hstack


   !> Performs QR factorization of the stack: A -> Q @ R
   subroutine stack_qr(A, Q, R, VRBZ_) !
   use string_lib, only: str
   use matrix_util, only: pprint_matrix
   class(stack), intent(in)  :: A
   type(stack), intent(out)  :: Q, R
   integer, intent(in), optional :: VRBZ_
   !
   integer :: VRBZ, sz, info, i, ng(3), nb(3), M,N,ML,NL,rank,nranks,nr,fd
   type(c_ptr) :: grid, descrR, descrQ
   double precision, allocatable, device, target :: d_tau(:)
   character(*),parameter :: subnam='[stack_qr]'
   character :: st
   double precision,allocatable:: h_e(:), h_tau(:)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (.not.associated(A%d_A)) error stop '[!]'//subnam//': A% d_A is not associated'
      if (.not.A% well_stacked()) error stop '[!]'//subnam//': A is not well stacked'
      R = A
      VRBZ = 0; if (present(VRBZ_)) VRBZ=VRBZ_

      !if (A% is_vstack()) allocate(d_tau(R%nx))
      !if (A% is_hstack()) allocate(d_tau(R%nx*R%nzb))
      M  = R%M;  N  = R%N
      allocate(d_tau(min(N, M)))
      grid = R% cusolvermp_grid_handle(VRBZ_)
      descrR = R% cusolvermp_matrix_handle(grid, VRBZ_)
      info = dist_geqrf(descrR, R%d_A, M, N, d_tau, VRBZ_=3)
      call stop_if_error(info,'[!]'//subnam//' in dist_geqrf: err = ')
!-!call R% save_to_ascii("tmp/RR_"//str(N)//"x"//str(N))

      st = 'h'; if (A% is_vstack()) st = 'v'
      ng = [A%nx, A%ny, A%nz] ! TODO: debug
      nb = [A%nxb, A%nyb, A%nzb]
      Q = eye_stack(A%comm, ng, nb, st, part_axis_=A%part_axis, VRBZ_=VRBZ)

      ! restore Q from Householder reflectors
      ! see https://docs.nvidia.com/hpc-sdk/archive/23.5/cusolvermp/functions.html#cusolvermpormqr
      descrQ = Q% cusolvermp_matrix_handle(grid, VRBZ_)
      info = dist_ormqr(descrR, R%d_A, 0, 0, M, N, &
                        descrQ, Q%d_A, min(M, N), d_tau, VRBZ_=VRBZ)
      call stop_if_error(info, '[!]'//subnam//' in dist_ormqr: err = ')
!-!call Q% save_to_ascii("tmp/QQ_"//str(N)//"x"//str(N))

      rank = R% rank
      nranks = R% nranks
!-!if (nranks > 1) then
!-!   fd = 728; h_tau = d_tau
!-!   open (fd, file="tmp/tau_"//str(N)//"x"//str(N)//"_"//str(rank)//".dat", status='replace')
!-!   write (fd, '(4100(ES22.14,1X))') h_tau 
!-!   close (fd)
!-!else
!-!   fd = 729; h_tau = d_tau
!-!   open (fd, file="tmp/tau_"//str(N)//"x"//str(N)//".dat", status='replace')
!-!   write (fd, '(4100(ES22.14,1X))') h_tau 
!-!   close (fd)
!-!endif
      ML = R%ML; NL = R%NL
      if (R% is_hstack()) then
         call zero_below_diag(R%d_A2, M,N, R%M,R%NB, 0,rank, ML,NL, 1,nranks)
      else
         call zero_below_diag(R%d_A2, M,N, R%MB,R%N, rank,0, ML,NL, nranks,1)
      endif
      deallocate(d_tau)

      ! Trim Q or R
      if (R% M > R% N) then
         ! Q Q Q   R R R     Q Q Q   R R R
         ! Q Q Q x 0 R R ->  Q Q Q x 0 R R
         ! Q Q Q   0 0 R     Q Q Q   0 0 R
         ! Q Q Q   0 0 0     Q Q Q   . . .
         call R% truncate_matrix(R% N, R% N) !, VRBZ_=4)
      endif
      if (Q% M < Q% N) then
         ! Q Q Q 0   R R R R     Q Q Q .   R R R R
         ! Q Q Q 0 x 0 R R R ->  Q Q Q . x 0 R R R
         ! Q Q Q 0   0 0 R R     Q Q Q .   0 0 R R
         call Q% truncate_matrix(Q% M, Q% M) !, VRBZ_=4)
      endif

      call destroy_cusolverMpMatrixDesc(descrR, VRBZ_)
      call destroy_cusolverMpMatrixDesc(descrQ, VRBZ_)
      call destroy_cusolverMpGrid(grid, VRBZ_)

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_qr


   !> Truncates a vStack or hStack unfolding matrix without multi-rank rearrangement
   !!
   !! Procedure:
   !!  - assert new_M <= M, new_N <= N
   !!  - set M to new_M, N to new_N
   !!  - deallocate truncated parts
   !!  - blocking dimensions MB, NB, as well as 3D box size nxb, nyb, nzb remain unchanged
   !!
   subroutine stack_truncate_matrix(A, new_M, new_N, VRBZ_) !
   class(stack), intent(INOUT) :: A
   integer, intent(IN)  :: new_M, new_N
   integer, intent(IN), optional :: VRBZ_
   !
   integer :: rank
   character(*),parameter :: subnam='[stack_truncate_matrix]'
   double precision,pointer,contiguous,device  :: d_A2(:,:)
   integer  :: VRBZ, ML, NL, new_ML, new_NL
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      if (.not.associated(A%d_A)) error stop '[!]'//subnam//': A% d_A is not associated'
      !if (.not.A% well_stacked())   error stop '[!]'//subnam//': A is not well stacked'
      if (new_M > A%M .or. new_N > A%N) error stop '[!]'//subnam//': new_M > M or new_N > N!'
      VRBZ = 0; if (present(VRBZ_)) VRBZ=VRBZ_
      rank = A% comm% rank
      call set_dist_mat_row_col(A%M, A%N, A%MB, A%NB, rank, &
                          A%numRowDevices, A%numColDevices, &
                          ML, NL, VRBZ)
      call set_dist_mat_row_col(new_M, new_N, A%MB, A%NB, rank, &
                          A%numRowDevices, A%numColDevices, &
                          new_ML, new_NL, VRBZ)
      if (new_ML*new_NL .eq. 0) then
         if (associated(A% d_A)) deallocate(A% d_A)
         allocate(A% d_A(1))
         A% d_A2(1:1,1:1) => A% d_A
      else
         if (ML /= new_ML .or. NL /= new_NL) then
            allocate(d_A2(new_ML,new_NL))
            call copy_arr2d(A% d_A2, d_A2, 1,new_ML, 1,new_NL, 1,1)
            deallocate(A% d_A)
            A% d_A(1:new_ML*new_NL)=> d_A2
            A% d_A2(1:new_ML,1:new_NL)=> d_A2
            d_A2=> null()
         endif
      endif
      if (A% is_hstack()) then ! TODO: only works for z-part
         A%ny = min(A% ny, new_N)
         A% M = new_M;   A% N  = new_N
         A%nx = new_M;   A% nz = (new_N - 1)/A% ny + 1
      else
         A%nx = min(A% nx, new_M)
         A% M = new_M;   A% N  = new_N
         A%ny = new_N;   A% nz = (new_M - 1)/A% nx + 1
      endif
      A% ML = new_ML !numroc(new_M, A% MB, rank, A% numRowDevices)
      A% NL = new_NL !numroc(new_N, A% NB, rank, A% numColDevices)
      A% nxl = A% nx
      A% nyl = A% ny
      A% nzl = numroc(A% nz, A% nzb, rank, A% nranks)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_truncate_matrix


   !> stack two 2D shallow stacks (matrices) together
   !!
   !! Equivalent to using matlab expression: C = [A,B]
   !!  - st = 'h': stack horizontally;
   !!  - st = 'v': stack vertically;
   !!
   function stack_matrices(A, B, st, VRBZ_) result(C)
   implicit none
   type(stack) :: C
   class(stack),   intent(INOUT) :: A, B
   character,         intent(IN) :: st
   integer, optional, intent(IN) :: VRBZ_
   !
   character(*), parameter  :: subnam = '[stack_matrices]'
   integer :: i, d, m1, m2, n1, n2, VRBZ
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif

      VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
      select case (st)
      case ('h')
         if (A% is_hstack()) call A% redist_to_vstack(VRBZ_)
         if (B% is_hstack()) call B% redist_to_vstack(VRBZ_)
         ! local horisontal matrix join
         if (A%N /= B%N) error stop '[!]'//subnam//': A%N /= B%N'
         C= shallow_stack(A%comm, A%M, A%N + B%N, 'v', MB_=A%MB, VRBZ_=VRBZ)
         call copy_arr2d(A% d_A2, C% d_A2, 1,A%ML, 1,A%NL, 1,1)
         call copy_arr2d(B% d_A2, C% d_A2, 1,B%ML, 1,B%NL, 1,A%NL+1)
      case ('v')
         if (A% is_vstack()) call A% redist_to_hstack(VRBZ_)
         if (B% is_vstack()) call B% redist_to_hstack(VRBZ_)
         ! local vertical matrix join
         if (A%M /= B%M) error stop '[!]'//subnam//': A%M /= B%M'
         C= shallow_stack(A%comm, A%M + B%M, A%N, 'h', MB_=A%NB, VRBZ_=VRBZ)
         call copy_arr2d(A% d_A2, C% d_A2, 1,A%ML, 1,A%NL, 1,1)
         call copy_arr2d(B% d_A2, C% d_A2, 1,B%ML, 1,B%NL, A%ML+1,1)
      case default
         error stop '[!]'//subnam//': bad value of st'
      end select

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end function stack_matrices


   !> Compute SVD factorization of a matrix: A = U * d * VT
   subroutine stack_svd(A, U, d_S, VT, VRBZ_)
   use string_lib, only: str
   use matrix_util, only: pprint_matrix
   implicit none
   class(stack), intent(IN) :: A
   type(stack), intent(OUT) :: U, VT
   double precision, allocatable, device, target, intent(OUT) :: d_S(:)
   integer, intent(in), optional :: VRBZ_
   !
   character(*),parameter :: subnam='[stack_svd]'
   integer :: VRBZ, sz, info, rank,nranks, M,N,ML,NL,MR,MRB
   type(c_ptr) :: grid, descrA, descrU, descrV, descrQ, descrR
   type(c_ptr) :: grid_MxN, matrix_MxN, grid_2Nx2N, matrix_2Nx2N
   type(c_ptr) :: matrix_NxN
   type(stack) :: Q, R, Z1, Z2, R1, R2, Q2, W
   character :: st
   double precision,allocatable,device,target:: d_d(:), d_e(:)
   double precision,allocatable,device,target:: d_tau(:), d_tau2(:)
   double precision,pointer,device:: d_tau2p(:), eig(:)
   double precision,allocatable:: h_e(:), h_tau(:)
   character(5) :: rankstr ! [0/1]
   integer :: fd
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if (present(VRBZ_)) VRBZ=VRBZ_
      st = 'h'; if (A% is_vstack()) st = 'v'
      rank = A% rank; nranks = A% nranks

      ! [Q, R] = qr(A,'econ');% <---- cusolverMpGeqrf: A -> Q * R
      if (VRBZ>2) call print_info(A,label_="=====> input matrix A:")      
      R = A
      M = R%M;  N = R%N
      allocate(d_tau(min(N, M)))
      grid_MxN = R% cusolvermp_grid_handle(VRBZ_)
      matrix_MxN = R% cusolvermp_matrix_handle(grid_MxN, VRBZ_)
      call stop_if_error( &
           dist_geqrf(matrix_MxN, R%d_A, M, N, d_tau, VRBZ_), &
          '[!]'//subnam//' in dist_geqrf: err = ')

      ! save output to Q (size MxN): will be needed later to recover U
      Q = R

      ! truncate R to make it square
      if (M > N) then
         ! Q Q Q   R R R     Q Q Q   R R R
         ! Q Q Q x 0 R R ->  Q Q Q x 0 R R
         ! Q Q Q   0 0 R     Q Q Q   0 0 R
         ! Q Q Q   0 0 0     Q Q Q   . . .
         call R% truncate_matrix(N, N) !, VRBZ_=4)
      endif

      ! zero R below diagonal
      if (R% is_hstack()) then
         call zero_below_diag(R%d_A2, N,N, N,R%NB, 0,rank, R%ML,R%NL, 1,nranks)
      else
         call zero_below_diag(R%d_A2, N,N, R%MB,N, rank,0, R%ML,R%NL, nranks,1)
      endif

      ! transpose matrix R
      call R% transpose_matrix() !< hstack
      if (VRBZ>2) then
         call print_info(R, label_='----------> R <------------')
         call R% pprint(label_="R transposed")
      endif
!call Q% save_to_ascii("tmp/Q_800x200")
!call R% save_to_ascii("tmp/R_800x200")
!h_tau = d_tau
!if (nranks > 1) then
!   fd = 728
!   open (fd, file="tmp/tau_"//str(rank)//".dat", status='replace')
!   write (fd, '(3200(ES22.14,1X))') h_tau 
!   close (fd)
!else
!   fd = 728
!   open (fd, file="tmp/tau.dat", status='replace')
!   write (fd, '(3200(ES22.14,1X))') h_tau 
!   close (fd)
!endif
!stop
!R = stack("tmp/R_800x200_221", A%comm)
!Q = stack("tmp/Q_800x200_221", A%comm)
      ! R2 = zeros(2*n);
      ! R2(n+1:end,1:n) = R';
      ! R2(1:n,n+1:end) = R   % <---- matrix/vector copy operations
      if (R% is_hstack()) then
         Z1 = shallow_stack(R%comm,  N,N,'h',MB_=min(N,R%NB), VRBZ_=VRBZ)
         Z1 = 0d0                                !< h: N  x N
         R1 = stack_matrices(Z1,  R, 'v', VRBZ_) !< h: 2N x N
         Z2 = shallow_stack(R%comm,2*N,N,'h',mb_=min(R%NB,N), VRBZ_=VRBZ)
         Z2 = 0d0                                !< h: 2N x N
         R2 = stack_matrices(R1, Z2, 'h', VRBZ_) !< v: 2N x 2N
         MR = 2*N; MRB = R2% MB
      else
         !!! TODO: implement h-stack here
         Z1 = shallow_stack(R%comm, N,  N,'v',mb_=min(R%MB,M),VRBZ_=VRBZ)
         Z1 = 0d0                                !< v:  N x N
         R1 = stack_matrices(R,  Z1, 'h', VRBZ_) !< v:  N x 2N
         Z2 = shallow_stack(R%comm, N,2*N,'v',mb_=min(R%MB,M),VRBZ_=VRBZ)
         Z2 = 0d0                                !< v:  N x 2N
         R2 = stack_matrices(Z2, R1, 'v', VRBZ_) !< h: 2N x 2N
         MR = R2%N; MRB = R2% NB
      endif

      !!-- Perform the tridiagonal factorization
      !    Q2 = R2
      !    call dsytrd('L', 2*N, Q2, 2*N, d, e, tau2, WRK, lwork, info)
      !    write (*,'("dsytrd exit info = ", I10)') info
      Q2 = R2
      grid_2Nx2N= set_cusolverMpGrid(R2%numRowDevices, R2%numColDevices, VRBZ_=VRBZ)
      matrix_2Nx2N= set_cusolverMpMatrixDesc(grid_2Nx2N, MR,MR, MRB,MRB, &
               R2%numRowDevices, R2%rank, VRBZ_=VRBZ)
      allocate(d_tau2(MR+1),d_d(MR),d_e(MR))
      d_tau2 = 0d0
      d_tau2p(1:MR)=> d_tau2(2:MR+1)
      call stop_if_error( &
         dist_sytrd(MR, matrix_2Nx2N, Q2%d_A, d_d, d_e, d_tau2p), &
        '[!]'//subnam//': in dist_sytrd, err = ')

      !!-- Output Q2 after SYTRD
      !    write (*,'("Q2 =")')
      !    do i = 1,2*N
      !       write (*,'(75(F10.4,1X))') Q2(i,:)
      !    enddo
      if (VRBZ>2) then
         call print_info(Q2,label_="------ Q2 after sytrd")
         call Q2% pprint(label_="------ Q2 after sytrd")
      endif

      !!-- Output diagonal, subdiagonal, and tau after SYTRD
      !    write (*,'(/"d =",505(ES14.5,1X))') d
      !    write (*,'(/"e =",505(F10.4,1X))') e
      !    write (*,'(/"tau =",505(F10.4,1X))') tau2
      if (VRBZ>2) then
         h_e = d_d
         print '(/"d =",505(ES14.5,1X))', h_e
         h_e = d_e
         print '(/"e =",505(ES14.5,1X))', h_e
         h_e = d_tau2
         print '(/"tau2 =",505(ES14.5,1X))', h_e
      endif


      !!-- Shift Q2 one column to the right
      !    do i = 2*N, 2,-1
      !       Q2(:,i)= Q2(:,i-1)
      !       tau2(i)= tau2(i-1)
      !    enddo
      !    tau2(1) = 0d0
      !    write (*,'("Q2 =")')
      !    do i = 1,2*N
      !       write (*,'(75(F10.4,1X))') Q2(i,:)
      !    enddo
      call shift_right_arr2d(Q2% d_A2, Q2% ML, Q2% NL, 1)
      if (VRBZ>2) call Q2% pprint(label_="---> Q2 after the right shift:")

      !!-- Compute eigenvalues
      !    call dstedc('I',2*N,d,e,Z2,2*N,WRK,lwork,IWRK,lwork,info)
      !    write (*,'(/"d =",505(ES14.5,1X))') d
      !    write (*,*) "Z2 ="
      !    do i = 1,2*N
      !       write (*,'(75(F10.4,1X))') Z2(i,:)
      !    enddo
      !    write (*,*)
      Z2 = Q2; Z2 = 0d0
      call stop_if_error( &
         dist_stedc(MR, d_d, d_e, matrix_2Nx2N, Z2% d_A), &
        '[!]'//subnam//': in dist_stedc, err = ')
      if (VRBZ>2) then
         h_e = d_d
         print '(/"eig =",505(ES14.5,1X))', h_e
         call Z2% pprint(label_="---> Z2 after STEDC:")
      endif

      !R2 = AxB_stacks(Z2, Z2, ops_='nt') !, VRBZ_=3)
      !call R2% pprint(label_="---> Z2^*Z2:")

      !!-- Apply ormq to recover Z2: Z2 = Q2*Z2
      !    call dormqr('L','N',2*N,2*N,2*N,Q2,2*N,tau2,Z2,2*N,WRK,lwork,info)
      call stop_if_error( &
         dist_ormqr(matrix_2Nx2N, Q2% d_A, 0, 0, MR, MR, &
                    matrix_2Nx2N, Z2% d_A, MR,  d_tau2), &
        '[!]'//subnam//': in dist_ormqr(Q2,Z2), err = ')
      if (VRBZ>2) then
         call Z2% pprint(label_="---> Z2 after ORMQR:", line_len_=200)
      endif

      !! ----------
      !! Assuming R = W @ D @ V', we can read off W and V from decomposition of R2:
      !!
      !!  R2 = | 0   R | = 1/2 | W  W |*|-D  0|*| W'  V'|
      !!       | R'  0 |       | V -V | | 0  D| | W' -V'|
      !!
      !!  | W  W | -> R1 = | W |
      !!  | V -V |         | v |
      Z1 = shallow_stack(R%comm,2*N,N,'v',MB_=Z2%MB, VRBZ_=VRBZ)
      call copy_arr2d(Z2%d_A2, Z1%d_A2, 1,Z2%ML,N+1,2*N,1,1)
      call Z1% redist_to_hstack() !< h: 2N x N
      VT = shallow_stack(R%comm,N,N,'h',val_=0d0, VRBZ_=VRBZ)
      call axpy_arr2d(dsqrt(2d0), Z1%d_A2, VT%d_A2, N+1,2*N, 1,VT%NL, 1,1)

      !!-- Recover U
      !    U = 0d0
      !    U(1:N,1:N) = W(1:N,1:N)
      !    call dormqr('L','N',M,N,N,QR,M,tau,U,M,WRK,lwork,info)
      U = shallow_stack(R%comm,M,N,'h',MB_=VT%NB,val_=0d0,VRBZ_=VRBZ)
      call axpy_arr2d(dsqrt(2d0), Z1%d_A2, U%d_A2, 1,N, 1,U%NL, 1,1)
      call U% redist_to_vstack()
      call stop_if_error( &
         dist_ormqr(matrix_MxN, Q%d_A, 0, 0, M, N, &
                    matrix_MxN, U%d_A, N, d_tau, VRBZ_=VRBZ), &
         "[!]"//subnam//": in dist_ormqr(Q,U)")
      if (VRBZ>2) call U% pprint(label_="--> matrix U after ORMQR:")
      call VT% transpose_matrix()

      !! Output eigenvalues, U, and VT
      if (VRBZ>2) then
         h_e = d_d
         print '(/"eig =",505(ES14.5,1X))', h_e(N+1:2*N)
         Z1 = shallow_stack(R%comm, N,N,'v',val_=0d0,VRBZ_=VRBZ)
         eig(1:N)=> d_d(N+1:2*N)
         call Z1% add_to_diag(eig, 0)
         W = AxB_stacks(Z1, VT) ! V is actually V transposed

         call U%  pprint(label_=">>> U: ")
         call VT% pprint(label_=">>> V^T: ")
         Z1 = AxB_stacks(U, W) - A
         !call Z1% pprint(label_="--> matrix U * D * V' - A:")
         print '(/"||U * D * V^ - A||_2 = ",ES14.7)', Z1% nrm2()
      endif

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_svd


   !> Compute tridiagonal factorization A -> P @ T @ P'
   subroutine stack_sytrd(A, P, T, VRBZ_)
   use string_lib, only: str
   use matrix_util, only: pprint_matrix
   implicit none
   class(stack), intent(IN) :: A
   type(stack), intent(OUT) :: P, T
   integer, intent(in), optional :: VRBZ_
   !
   character(*),parameter :: subnam='[stack_sytrd]'
   integer :: VRBZ, sz, rank,nranks, M, MA
   type(c_ptr) :: grid, descrA
   type(stack) :: A1, II, TP1
   character :: st
   integer :: nr, ng(3), nb(3), nRowDevices, nColDevices
   double precision,allocatable,device,target:: d_d(:),d_e(:),d_tau(:)
   double precision,pointer,device:: d_tau1(:)
   double precision,allocatable:: h_e(:)
#  ifdef VERBOSE4
      print '("[+]'//subnam//' entry")'
#  endif
      VRBZ = 0; if (present(VRBZ_)) VRBZ=VRBZ_
      st = 'h'; if (A% is_vstack()) st = 'v'
      M  = A% M
      MA = A% MB
      if (M.ne.A% N) error stop '[!]'//subnam//': A is not a square matrix'
      if (st.ne.'v') error stop '[!]'//subnam//': not implemented for hstack'

      ! create CuSOLVERMP grids
      nRowDevices= A% nranks
      nColDevices= 1
      grid = set_cusolverMpGrid(nRowDevices, nColDevices, VRBZ_=VRBZ)

      ! create CuSOLVERMP matrixDescriptor
      descrA= set_cusolverMpMatrixDesc(grid, M, M, MA, MA, &
              nRowDevices, A%rank, VRBZ_=VRBZ)

      ! distributed A -> P @ T @ P' factorization:
      allocate(d_tau(M+1), d_d(M), d_e(M))
      d_tau1(1:M)=> d_tau(2:M+1)
      d_tau = 0d0
      T = A
      call stop_if_error(dist_sytrd( &
           M, descrA, T% d_A, d_d, d_e, d_tau1), &
          '[!]'//subnam//': in dist_sytrd: err = ')

      ! Lower left corner (below 1st subdiagonal) contain N-1 Householder reflectors:
      !  || d_1   *    *   ...  *  ||
      !  !!  1   d_2   *   ...  *  ||
      !  || v_1   1   d_3  ...  *  ||
      !  || v_2  u_2   1   ...  *  ||
      !  || . . . . . . . . , , ,  ||
      !  || v_n  u_n  w_n  ... d_n ||
      !
      ! To recover the orthogonal matrix P, we need to shift these columns to the right,
      ! by one, and apply the function "ormqr", similar to the one applied after dgeqrt
      ! function. It also needs factors tau_k (also present in Householder reflectors):
      !
      ! H_k = II - tau_k v_k v_k^T = II - 2 v_k v_k^T / (v_k^T @ v_k)
      !
      ! shift matrix T (TODO: generalize to T% shift_right)
      call shift_right_arr2d(T% d_A2, T% ML, T% NL, 1)

      ! recover the matrix P
      ng = [A%nx,  A%ny,  A%nz ]
      nb = [A%nxb, A%nyb, A%nzb]
      P  = eye_stack(A%comm, ng, nb, st, VRBZ_=VRBZ)
      call stop_if_error(dist_ormqr( &
           descrA, T% d_A, 0, 0, M, M, &
           descrA, P% d_A, M, d_tau, VRBZ_=VRBZ), &
          '[!]'//subnam//': in dist_ormqr: err = ')

      ! The function sytrd returns 1D arrays containing diagonal (d_d) and subdiagonal
      ! (d_e) of the symmetric tridiagonal matrix T. We recover it here in full form:
      T = 0d0
      call T% add_to_diag(d_d, 0)
      call T% add_to_diag(d_e, 1)
      call T% add_to_diag(d_e,-1)

#  ifdef VERBOSE4
      print '("[+]'//subnam//' exit")'
#  endif
   end subroutine stack_sytrd


  !> Writes a distributed array (dray) to a file in ASCII format
  !!
  !! Data format:
  !! - <file_prefix>_<rank#>.dat: one file per rank;
  !! - first row: "# [h|v]stack [<rank>/<nranks>]"
  !! - second row: 3D dimensions: nx,ny,nz
  !! - third row: blockiing dimensions: nxb,nyb,nzb
  !! - fourth row: local dimensions: nxl,nyl,nzl
  !! - fifth row: matrix dimensions: M, N, MB, NB, ML, NL
  !!
  !! File example:
  !! >>>>>>>>>>>>>>>> BEGIN FILE Q_0.dat out ot 3 files >>>>>>>>>>>>>>
  !! # hstack [0/3]
  !! 3,2,5
  !! 3,2,2
  !! 3,2,2
  !! 3,10,3,4,3,4
  !! 0.0260   0.0477  -0.0625  -0.0094   0.0160
  !! 0.0546  -0.0274   0.0049  -0.0650   0.0034
  !! 0.0212   0.0432   0.0296   0.0016  -0.0323
  !! >>>>>>>>>>>>>>>> BEGIN FILE Q_1.dat out ot 3 files >>>>>>>>>>>>>>
  !! # hstack [1/3]
  !! 3,2,5
  !! 3,2,2
  !! 3,2,2
  !! 3,10,3,4,3,4
  !! 0.0260   0.0477  -0.0625  -0.0094   0.0160
  !! 0.0546  -0.0274   0.0049  -0.0650   0.0034
  !! 0.0212   0.0432   0.0296   0.0016  -0.0323
  !! >>>>>>>>>>>>>>>> BEGIN FILE Q_2.dat out ot 3 files >>>>>>>>>>>>>>
  !! # hstack [2/3]
  !! 3,2,5
  !! 3,2,2
  !! 3,2,1
  !! 3,10,3,4,3,2
  !! 0.0260   0.0477  -0.0625
  !! 0.0546  -0.0274   0.0049
  !! 0.0212   0.0432   0.0296
  !! <<<<<<<<<<<<<<<< END OF 3 FILES <<<<<<<<<<<<<<<<
  subroutine stack_save_to_ascii(this, fnamprefix, VRBZ_)
  use string_lib, only: str
  implicit none
  class(stack),      intent(IN) :: this
  character(*),      intent(IN) :: fnamprefix
  integer, optional, intent(IN) :: VRBZ_
  !
  character(*), parameter :: subnam = '[stack_save_to_ascii]'
  integer :: i, j, rank, nranks, VRBZ
  integer, parameter :: fd = 773
  double precision, allocatable :: h_A2(:,:)
  character(200) :: fname

     VRBZ = 0; if (present(VRBZ_)) VRBZ = VRBZ_
     if (.not.this% is_allocated) then
        print '("[w]'//subnam//': stack is not allocated")'
        return
     endif

     rank = this% rank
     nranks = this% nranks
     if (nranks > 1) then
         write(fname,'(A)') fnamprefix//"_"//str(rank)//".dat"
     else
         write(fname,'(A)') fnamprefix//".dat"
     endif

     if (VRBZ) then
        print '("writing stack to a file: ",A)', fname
     endif
     open (fd, file=fname, status='replace')

     ! write the first line
     if (this% is_hstack()) then
        write (fd, '("# hstack ['//str(rank)//'/'//str(nranks)//']")')
     else
        write (fd, '("# vstack ['//str(rank)//'/'//str(nranks)//']")')
     endif

     ! write header
     write (fd,'("'//str(this%nx)//',' &
                   //str(this%ny)//',' &
                   //str(this%nz)//'")')
     write (fd,'("'//str(this%nxb)//',' &
                   //str(this%nyb)//',' &
                   //str(this%nzb)//'")')
     write (fd,'("'//str(this%nxl)//',' &
                   //str(this%nyl)//',' &
                   //str(this%nzl)//'")')
     write (fd,'("'//str(this%M) //','//str(this%N) //',' &
                   //str(this%MB)//','//str(this%NB)//',' &
                   //str(this%ML)//','//str(this%NL)//'")')

     ! write local chunk
     if (this% ML * this% NL > 0) then
        h_A2 = this% d_A2
        do i=1,this% ML
           do j=1,this% NL
              write(fd, '(ES22.14,1X)', advance='no') h_A2(i,j)
           enddo
           write(fd, '("")') ! EOL
        enddo
        deallocate(h_A2)
     endif

     close(fd)
     print '("File '//trim(fname)//' has been written.")'

  end subroutine stack_save_to_ascii


  !> Reads a stack from files in ASCII format (1 file per rank)
  !!
  !! Data format: see description to stack_save_to_ascii
  !!
  function stack_read_from_ascii(fnamprefix, comm, VRBZ_) result(this)
  use string_lib, only: split, isplit, str
  implicit none
  type(stack) :: this
  type(super_comm),  intent(IN), target :: comm
  character(*),      intent(IN) :: fnamprefix
  integer, optional, intent(IN) :: VRBZ_
  !
  character(*), parameter :: subnam = '[stack_read_from_ascii]'
  integer :: VRBZ
  integer :: i, j, nlines, rank, nranks
  integer, parameter :: fd = 775
  integer, parameter :: max_line_len = 2048*8
  integer, allocatable :: nn(:),nnb(:),nnl(:),mm(:)
  integer :: M,N,MB,NB,ML,NL
  double precision, allocatable :: h_A2(:,:)
  character(:), allocatable :: strarr(:)
  character(200) :: fname, rankstr, hline
  character :: st

     VRBZ = 0; if (present(VRBZ_)) VRBZ = VRBZ_
     rank = comm% ctx_rank
     nranks = comm% ctx_size
     if (VRBZ>0) &
        print '("reading dthor object with prefix: ",A)', trim(fnamprefix)
     if (nranks > 1) then
         write(fname,'(A)') fnamprefix//"_"//str(rank)//".dat"
     else
         write(fname,'(A)') fnamprefix//".dat"
     endif

     nlines = file_nlines(fname)
     open (fd, file=fname, status='old')
     read (fd, '(A)') hline
     call split(strarr, hline, cdel_=' ')
     if (strarr(2).eq.'vstack') then
        st = 'v'
     else if (strarr(2).eq.'hstack') then
        st = 'h'
     else
        error stop '[!]'//subnam//': unknown stacking '//strarr(2)
     endif

     rankstr = '['//str(rank)//'/'//str(nranks)//']'
     if (strarr(3).ne.trim(rankstr)) then
        error stop '[!]'//subnam//': wrong number of ranks'
     endif

     read (fd, '(A)') hline; nn = isplit(hline, cdel_=',')
     read (fd, '(A)') hline; nnb = isplit(hline, cdel_=',')
     read (fd, '(A)') hline; nnl = isplit(hline, cdel_=',')
     read (fd, '(A)') hline; mm = isplit(hline, cdel_=',')
     M  = mm(1);  N = mm(2)
     MB = mm(3); NB = mm(4)
     ML = mm(5); NL = mm(6)

     if (VRBZ>2) then
        print '("[i]'//trim(rankstr)//subnam//': nn = ",3(I10,1X))', nn
        print '("[i]'//trim(rankstr)//subnam//': nb = ",3(I10,1X))', nnb
        print '("[i]'//trim(rankstr)//subnam//': nl = ",3(I10,1X))', nnl
        print '("[i]'//trim(rankstr)//subnam//': M, N   = ",I10,1X,I10)', M, N
        print '("[i]'//trim(rankstr)//subnam//': MB, NB = ",I10,1X,I10)', MB, NB
        print '("[i]'//trim(rankstr)//subnam//': ML, NL = ",I10,1X,I10)', ML, NL
     endif

     !!! TODO: only implemented for z-partitioning
     !!! TODO: place bunch of assertions to check consistency of the input files
     
     if (st.eq.'h') then
        if (N ==  nn(2)*nn(3)) then ! well-stacked
           this = empty_stack(comm, nn, nnb, 'h', VRBZ_=VRBZ)
        else
           this = empty_stack(comm, nn, nnb, 'h', morn_=N, VRBZ_=VRBZ)
        endif
     else
        if (M ==  nn(1)*nn(3)) then ! well-stacked
           this = empty_stack(comm, nn, nnb, 'v', VRBZ_=VRBZ)
        else
           this = empty_stack(comm, nn, nnb, 'v', morn_=M, VRBZ_=VRBZ)
        endif
     endif

     if (ML*NL > 0) then
        allocate(h_A2(ML, NL))
        do i=1,this% ML
           read(fd,*) h_A2(i,1:this%NL)
        enddo
        this% d_A2(1:this%ML,1:this%NL)= h_A2(1:this%ML,1:this%NL)
        deallocate(h_A2)
     endif
     close (fd)

     deallocate(nn, nnb, nnl, mm, strarr)
     if (VRBZ > 0) then
        print '("[i]'//trim(rankstr)//subnam//': lines read = ",I10)', nlines
        if (rank.eq.0) then
           print '("[i]'//subnam//', stack: ")'
           call print_info(this)
        endif
     endif

  end function stack_read_from_ascii


  !> count the number of lines in a file
  integer function file_nlines(fname) result(nlines)
  character(*), intent(in) :: fname
  integer, parameter :: fd = 776

     nlines = 0
     open (fd, file=fname, status='old')
     do
        read(fd, *, end=10)
        nlines = nlines + 1
     enddo
10   close(fd)
  end function file_nlines


   subroutine stop_if_error(err, msg)
   use, intrinsic :: iso_fortran_env, only: error_unit
   integer, intent(in) :: err
   character(*), intent(in) :: msg
      if (err.ne.0) then
         write(error_unit,'(A,I8)') msg, err
         error stop
      endif
   end subroutine

#undef VERBOSE4
end module distributed_arrays
