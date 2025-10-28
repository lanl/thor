!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file dist_thor.f90
!! @author Ismael Djibrilla, Oleg Korobkin
!! @date  August 2025
!! @brief Distributed tensor-train format
!!
!#define VERBOSE4

module thor
use distributed_comms,  only: super_comm
use distributed_arrays, only: stack, dealloc_stack, empty_stack, &
                              AxB_stacks, AxBT_stacks, ATxB_stacks
implicit none

  !//////////////////////////////////////////////////
  !//
  !//  dthor: double precision distributed tensor in
  !//         the tensor-train format
  !//
  !//////////////////////////////////////////////////
  type, public :: dthor
     integer          :: d      ! tensor rank (number of dimensions)
     integer, pointer :: n(:)   ! mode sizes
     integer, pointer :: r(:)   ! tensor ranks
     type(stack), allocatable  :: core(:)
     logical          :: is_allocated
  contains
     procedure, pass(self) :: dthor_assign
     generic :: assignment(=) => dthor_assign
     procedure :: round_rand_orth => dthor_round_rand_orth
     procedure :: info => dthor_print_info
     procedure :: pprint_cores => dthor_pprint_local_cores
     procedure :: is_hstack => dthor_is_hstack
     procedure :: is_vstack => dthor_is_vstack
     procedure :: to_hstack => dthor_to_hstack
     procedure :: to_vstack => dthor_to_vstack
     final     :: dealloc_dthor
  end type

  interface dthor ! constructor interface
     procedure :: empty_dthor
     procedure :: dthor_from_cell_array
  end interface dthor
  ! TODO: add arithmetic operations (addition, subtraction,
  !       Hadamard multiplication)
  !! OPERATORS / OVERLOADING
  interface operator (+)
     module procedure dthor_add
  end interface
  interface operator (-)
     module procedure dthor_sub
  end interface
  interface operator (*)
     module procedure dthor_mul_d, d_mul_dthor, dthor_mul_i, i_mul_dthor
  end interface
  interface operator (/)
     module procedure dthor_div_d, dthor_div_i
  end interface

contains

  !> creates empty tensor from arrays of dimensions and ranks (optional)
  function empty_dthor(comm, n, r_, stacking_, VRBZ_) result(this)
  use distributed_arrays, only: construct_stack
  use string_lib, only: str
  type(dthor):: this
  type(super_comm),    intent(IN), target :: comm
  integer,             intent(IN) :: n(:)
  integer,   optional, intent(IN) :: r_(0:size(n))
  character, optional, intent(IN) :: stacking_
  integer,   optional, intent(IN) :: VRBZ_
  !
  character(*),parameter :: subnam='[empty_dthor]'
  integer   :: d, k, sh(3), ng(3), nb(3), nranks, VRBZ
  character :: stacking
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     stacking = 'h'; if (present(stacking_)) stacking = stacking_
     VRBZ=0; if (present(VRBZ_)) VRBZ = VRBZ_
     nranks = comm% ctx_size
     d= size(n)
     allocate(this% n(d), this% r(0:d), this% core(d))
     this% d = d
     this% n(1:d)= n(1:d)
     this% r(0:d)= 1; if (present(r_)) this% r(0:d)= r_(0:d)
     do k=1,d
        ng = [this% r(k-1), this% r(k), this% n(k)]
        nb = [this% r(k-1), this% r(k),(this% n(k) - 1)/nranks + 1]
        this% core(k) = stack(comm, ng, nb, stacking, VRBZ_=VRBZ)
        if (VRBZ>1) print '("[i]'//subnam//': this% core('//str(k)//') allocated")'
     enddo
     this% is_allocated= .true.
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end function empty_dthor


  !> constructor from a cell array
  function dthor_from_cell_array(comm, cores, stacking, VRBZ_) result(this)
  use cell_arrays_module, only: cell3d_array
  type(dthor):: this
  type(super_comm),   intent(IN) :: comm
  type(cell3d_array), intent(IN) :: cores
  character,          intent(in) :: stacking
  integer, optional,  intent(in) :: VRBZ_
  !
  integer:: d, k, info, sh(3), VRBZ
  character(*),parameter   :: subnam='[dthor_from_cell_array]:'
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
     d = size(cores% cells)
     this% r(0)= 1
     allocate(this%core(d))
     do k=1,d
        sh= shape(cores% cells(k)% arr)
        this% n(k)= sh(2)
        this% r(k)= sh(3)
        this% core(k) = stack(comm, cores% cells(k)% arr, 'v', VRBZ_=VRBZ)
        if (stacking.eq.'h') call this% core(k)% to_hstack(VRBZ_)
     enddo
     this% is_allocated= .true.
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end function dthor_from_cell_array


  !> destructor
  subroutine dealloc_dthor(self)
  type(dthor), intent(INOUT) :: self
  !
  character(*),parameter :: subnam='[dealloc_dthor]'
  integer :: k
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     if(allocated(self%core)) then
        do k=1, self% d
           call dealloc_stack(self% core(k))
        enddo
        deallocate(self%core)
        deallocate(self% n, self% r)
        nullify(self% n)
        nullify(self% r)
     endif
     self% is_allocated= .false.
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end subroutine dealloc_dthor


  !> the dthor assignment operator, A = B
  subroutine dthor_assign(self, other)
  class(dthor), intent(INOUT) :: self
  type(dthor),     intent(IN) :: other
  !
  integer:: k, d
  character(*),parameter :: subnam='[dthor_assign]'
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     d= other% d
     if (self% d.ne.d) then
        if (self% is_allocated) call dealloc_dthor(self)
        self% d = d
        allocate(self% n(d), self% r(0:d))
        allocate(self% core(d))
     endif
     if (.not.allocated(self% n)) allocate(self% n(d))
     if (.not.allocated(self% r)) allocate(self% r(0:d))
     self% n(1:d)= other% n(1:d)
     self% r(0:d)= other% r(0:d)
     if (.not.allocated(self% core)) allocate(self% core(d))
     do k=1,d
        self% core(k)= other% core(k)
     enddo
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end subroutine dthor_assign


  !> print tensor info
  subroutine dthor_print_info(this, on_all_ranks_)
  implicit none
  class(dthor),      intent(IN) :: this
  logical, optional, intent(IN) :: on_all_ranks_
  !
  character(len=1) dstr
  integer :: d
  logical :: on_all_ranks
  character(*), parameter :: subnam = '[dthor_print_info]'

     on_all_ranks = .true.
     if (present(on_all_ranks_)) on_all_ranks = on_all_ranks_

     d = this%d
     if (d < 1) then
        write (*,'("'//subnam//'": empty dthor object")')
        return
     endif
     if (on_all_ranks.and.this% core(1)% rank > 0) return
     write(*,'(A,I4,A)') 'dthor[1:', d,']:' 
     write (dstr, '(I1)') int(log10(dble(maxval(this%n)))) + 2
     if(all(this%n.eq.this%n(1)))then
        write(*,'(a,1x,i'//dstr//',a)') 'n: ',this%n(1),' for all modes'
     else
        write(*,'(a,1x,2048i'//dstr//')') 'n: ',this%n(1:d)
     end if
     if (maxval(this%r).gt.maxval(this%n)) &
       write (dstr, '(I1)') int(log10(dble(maxval(this%r(1:d))))) + 2
     write(*,'("r: ",i'//dstr//')',advance='no') this%r(0)
     if(all(this%r(2:d-1).eq.this%r(2)).and.d.gt.5) then
       write(*,'(1x,i'//dstr//'" ... "i'//dstr//')',advance='no') this%r(1),this%r(2)
     else
       write(*,'(2048i'//dstr//')',advance='no') this%r(1:d-1)
     end if
     write(*,'(i'//dstr//')') this%r(d)
  endsubroutine dthor_print_info


  !> pretty-print tensor cores
  subroutine dthor_pprint_local_cores(this, on_all_ranks_, frmt_)
  use matrix_util, only: pprint_matrix3d
  implicit none
  class(dthor),           intent(IN) :: this
  logical, optional,      intent(IN) :: on_all_ranks_
  character(*), optional, intent(IN) :: frmt_
  !
  character(len=20) :: frmt
  integer :: k, d, r1, r2, nn, nzl
  integer :: rank, nranks, r, rmax
  logical :: on_all_ranks
  double precision, allocatable, target :: h_A(:)
  double precision, pointer :: h_A3(:,:,:)
  character(*), parameter :: subnam = '[dthor_pprint_local_cores]'

     on_all_ranks = .true.
     if (present(on_all_ranks_)) on_all_ranks = on_all_ranks_
     frmt = "(F8.4)"; if (present(frmt_)) frmt= trim(frmt_)
     d = this% d
     if (d < 1) then
        print '("'//subnam//': empty dthor, nothing to print")'
        return
     endif
     rank   = this% core(1)% rank
     nranks = this% core(1)% nranks
     rmax   = nranks
     if (.not.on_all_ranks) rmax = 1
     do r = 1,rmax
        if(.not.on_all_ranks .or. r-1.eq.rank) then
           print *
           do k = 1,d
              r1 = this% core(k)% nx
              r2 = this% core(k)% ny
              nn = this% core(k)% nz
              nzl = this% core(k)% nzl
              print '("[",I4,"/",I4,"]== core #",I4,": ' // &
                    '[r_{k-1},n_k,r_k] = [",3(I6,1X)"], '    // &
                    'nzl = ",I6,", is_vstack = ",L1)', &
                    rank, nranks, k, r1, nn, r2, nzl, this% core(k)% is_vstack()
              h_A = this% core(k)% d_A
              if (this% core(k)% is_vstack()) then
                 h_A3(1:r1, 1:nzl, 1:r2)=> h_A
                 call pprint_matrix3d(h_A3, frmt_=frmt, line_len_=120, flip23_=.true.)
              else
                 h_A3(1:r1, 1:r2, 1:nzl)=> h_A
                 call pprint_matrix3d(h_A3, frmt_=frmt, line_len_=120)
              endif
           enddo
        endif ! not on_all_ranks or rank == 1
        if (on_all_ranks) call this%core(1)%comm%barrier()
     enddo
  end subroutine dthor_pprint_local_cores


  !> Create a dtt_tensor with specified dims and ranks with random cores
  !!
  !! Parameters:
  !! * comm      : supercommunicator
  !! * n(:)      : tensor modes
  !! * r_(:)     : tensor ranks [1]
  !! * stacking_ : core stacking, 'h'|'v' ['h']
  !! * VRBZ_     : verbosity flag [0]
  !!
  function random_dthor(comm, n, r_, stacking_, VRBZ_) result(this)
  use string_lib, only: str
  use rnd_lib, only: randn
  type(dthor):: this
  type(super_comm),    intent(IN), target :: comm
  integer,             intent(IN) :: n(:)
  integer,   optional, intent(IN) :: r_(0:size(n))
  character, optional, intent(IN) :: stacking_
  integer,   optional, intent(IN) :: VRBZ_
  !
  integer:: d, k, info, sh(3), VRBZ, ML, NL
  character(*),parameter          :: subnam='[random_dthor]'
  double precision, allocatable   :: local_core(:)
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     VRBZ = 0; if(present(VRBZ_)) VRBZ = VRBZ_
     d = size(n)
     this = empty_dthor(comm, n, r_, stacking_, VRBZ_)
     do k=1,d
        ML = max(1, this% core(k)% ML)
        NL = max(1, this% core(k)% NL)
        local_core = randn(ML*NL)
        this% core(k)% d_A(1:ML*NL) = local_core(1:ML*NL)
     enddo
     if (d > 0) deallocate(local_core)
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end function random_dthor



  !> Partial contraction subroutine for the randomized rounding algorithm
  !!
  !! From: Al Daas et al. (2023), "RANDOMIZED ALGORITHMS FOR ROUNDING IN THE
  !!       TENSOR-TRAIN FORMAT", doi:10.1137/21M1451191
  !! Algorithm 3.1 TT-rounding: Randomize-then-orthogonalize.
  subroutine partial_contractions_rl(W, X, Y, VRBZ_)
  use string_lib, only: str
  type(stack),       intent(IN) :: W(:)
  type(dthor),       intent(IN) :: X, Y
  integer, optional, intent(IN) :: VRBZ_
  !
  type(stack) :: Zn
  integer     :: d, info, n, VRBZ
  character(*),parameter    :: subnam='[partial_contractions_rl]'
double precision, allocatable :: h_A(:)
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
     d = size(W)
     if(X% d /= d) error stop '[!]'//subnam//': X% '//str(X%d)//' /= '//str(d)
     if(Y% d /= d) error stop '[!]'//subnam//': Y% '//str(Y%d)//' /= '//str(d)
     call X%core(d)% to_hstack(VRBZ_)
     call Y%core(d)% to_hstack(VRBZ_)
     W(d-1) = AxBT_stacks(X%core(d), Y%core(d), VRBZ_)
     do n=d-1, 2, -1
        call X%core(n)% to_vstack(VRBZ_)
        call W(n)% redist_to_vstack(VRBZ_)
        Zn = AxB_stacks(X%core(n), W(n), VRBZ_=VRBZ)
        if (VRBZ>1) print '("[i]'//subnam//': computed Zn for iteration '//str(n)//'")'
        call Y%core(n)% to_hstack(VRBZ_)
        call Zn% to_hstack(VRBZ_)
        W(n-1) = AxBT_stacks(Zn, Y%core(n), VRBZ_)
        if (VRBZ>1) print '("[i]'//subnam//': finished contracting core '//str(n)//'")'
     end do
     call W(1)% redist_to_vstack(VRBZ_)
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end subroutine partial_contractions_rl


  !> TT rounding function using randomization
  !!
  !! From: Al Daas et al. (2023), "RANDOMIZED ALGORITHMS FOR ROUNDING IN THE
  !!       TENSOR-TRAIN FORMAT", doi:10.1137/21M1451191
  !! Algorithm 3.1 TT-rounding: Randomize-then-orthogonalize.
  function dthor_round_rand_orth(Y, rn, RTT_, VRBZ_) result(X)
use string_lib, only: str
use matrix_util, only: pprint_matrix
  type(dthor) :: X
  class(dthor),          intent(IN) :: Y         !< input tensor
  integer,               intent(IN) :: rn(0:Y%d) !< array of target ranks
  type(dthor), optional, intent(IN) :: RTT_      !< initial guess [random]
  integer,     optional, intent(IN) :: VRBZ_
  !
  type(stack), allocatable :: W(:)
  type(stack) :: Zn, Yn, Mn, DUMMY
  type(dthor) :: RTT
  integer     :: d, info, n, VRBZ
  character(*),parameter     :: subnam='[dthor_round_rand_orth]'
double precision, allocatable :: h_A(:), h_A2(:,:)
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
     d = Y%d

     if (present(RTT_)) then
        RTT = RTT_
        if (RTT% d /= d) &
           error stop '[!]'//subnam//': wrong dimension in RTT_'
        if (any(RTT%n(1:d) /= Y%n(1:d))) &
           error stop '[!]'//subnam//': incompatible mode size n(:) in RTT_'
        if (any(RTT%r(0:d) /= rn(0:d))) &
           error stop '[!]'//subnam//': incompatible ranks r(:) in RTT_'
     else
        RTT =  random_dthor(Y%core(1)%comm, Y%n, r_=rn, VRBZ_=VRBZ)
     endif
     allocate(W(d))
     call partial_contractions_rl(W, Y, RTT, VRBZ_)
     X = empty_dthor(Y%core(1)%comm, Y%n, r_=rn)
     X%core(1) = Y%core(1)

     do n = 1, d-1
        call X% core(n)% to_vstack(VRBZ_)
        Zn = X% core(n)
        Yn = AxB_stacks(Zn, W(n), VRBZ_=VRBZ)  !< vStack
        call Yn% QR(X% core(n), DUMMY, VRBZ_)  !< X% core(n) is vStack
        Mn = ATxB_stacks(X%core(n), Zn, VRBZ_) !< vStack
        call Mn% redist_to_hstack(VRBZ_)       !< hStack
        call Y%core(n+1)% to_hstack(VRBZ_)     !< hStack
        X%core(n+1) = AxB_stacks(Mn, Y%core(n+1), VRBZ_=VRBZ)  !< hStack
     end do
     call X% core(d)% to_vstack(VRBZ_)
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end function dthor_round_rand_orth


  !> AXPY operation on two dthor objects: Z = alpha*X + Y
  subroutine dthor_axpy_baseclass(Z, alpha, X, Y, VRBZ_)
  use matrix_ops, only: copy_arr3d, axpy_arr3d
  implicit none
  type(dthor),         intent(INOUT) :: Z
  double precision,       intent(IN) :: alpha
  class(dthor),   target, intent(IN) :: X, Y
  integer,      optional, intent(IN) :: VRBZ_
  !
  character(*),parameter     :: subnam='[dthor_add_baseclass]'
  integer :: k,d,nn,ar1,br1,ar2,br2,VRBZ
  integer,allocatable :: r(:),n(:)
  double precision, device, pointer :: d_X3(:,:,:),d_Y3(:,:,:),d_Z3(:,:,:)
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif

     VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
     if (X%d /= Y%d) error stop '[!]'//subnam//': tensor rank mismatch'
     if (X%d < 1) return
     ! TODO: implemented for vStack only as of now (O.K.)
     if (.not.X% is_vstack()) error stop '[!]'//subnam//': X is not vStack'
     if (.not.Y% is_vstack()) error stop '[!]'//subnam//': Y is not vStack'
     d = X%d
     if (any(X%n(1:d) /= Y%n(1:d))) error stop '[!]'//subnam//': dimension mismatch'
     allocate(r(0:d),n(d))
     r(0:d) = X%r(0:d) + Y%r(0:d)
     r(0) = 1; r(d) = 1
     n(1:d) = X%n(1:d)
     Z = empty_dthor(X%core(1)%comm, n, r_=r, stacking_='v', VRBZ_=VRBZ)
     deallocate(n, r)

     nn = Z% core(1)% nzl
     ar2 = X% r(1)
     br2 = Y% r(1)
     d_Z3(1:1,1:nn,1:ar2+br2)=> Z% core(1)% d_A
     d_X3(1:1,1:nn,1:ar2)=> X% core(1)% d_A
     d_Y3(1:1,1:nn,1:br2)=> Y% core(1)% d_A
     if (dabs(alpha-1d0) < 1d-15) then
        call copy_arr3d(d_X3,d_Z3, 1,1, 1,nn, 1,ar2, 1,1,1)
     else
        d_Z3 = 0d0
        call axpy_arr3d(alpha,d_X3,d_Z3, 1,1, 1,nn, 1,ar2, 1,1,1)
     endif
     call copy_arr3d(d_Y3,d_Z3, 1,1, 1,nn, 1,br2, 1,1,ar2+1)

     do k=2,d-1
        nn = Z% core(k)% nzl
        ar1= X% r(k-1)
        ar2= X% r(k)
        br1= Y% r(k-1)
        br2= Y% r(k)

        d_Z3(1:ar1+br1, 1:nn, 1:ar2+br2)=>  Z% core(k)% d_A
        d_X3(1:ar1,1:nn,1:ar2)=> X% core(k)% d_A
        d_Y3(1:br1,1:nn,1:br2)=> Y% core(k)% d_A

        d_Z3= 0d0
        call copy_arr3d(d_X3,d_Z3, 1,ar1, 1,nn, 1,ar2, 1,1,1)
        call copy_arr3d(d_Y3,d_Z3, 1,br1, 1,nn, 1,br2, ar1+1,1,ar2+1)
     enddo

     nn = Z% core(d)% nzl
     ar1= X% r(d-1)
     br1= Y% r(d-1)

     d_Z3(1:ar1+br1,1:nn,1:1)=> Z% core(d)% d_A
     d_X3(1:ar1,1:nn,1:1)=> X% core(d)% d_A
     d_Y3(1:br1,1:nn,1:1)=> Y% core(d)% d_A
     !if (dabs(alpha-1d0) < 1d-15) then
     !   call copy_arr3d(d_X3,d_Z3, 1,ar1, 1,nn, 1,1, 1,1,1)
     !else
     !   d_Z3 = 0d0
     !   call axpy_arr3d(alpha,d_X3,d_Z3, 1,ar1, 1,nn, 1,1, 1,1,1)
     !endif
     call copy_arr3d(d_X3,d_Z3, 1,ar1, 1,nn, 1,1, 1,1,1)
     call copy_arr3d(d_Y3,d_Z3, 1,br1, 1,nn, 1,1, ar1+1,1,1)

# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end subroutine dthor_axpy_baseclass


  type(dthor) function dthor_add(A, B) result(C)
  type(dthor),intent(in),target :: A, B
# ifdef VERBOSE4
     print '("[+][dthor_add] entry")'
# endif

     call dthor_axpy_baseclass(C, 1d0, A, B)

# ifdef VERBOSE4
     print '("[+][dthor_add] exit")'
# endif
  end function dthor_add


  type(dthor) function dthor_sub(A, B) result(C)
  type(dthor),intent(in),target :: A, B
# ifdef VERBOSE4
     print '("[+][dthor_sub] entry")'
# endif

     call dthor_axpy_baseclass(C, -1d0, B, A)

# ifdef VERBOSE4
     print '("[+][dthor_sub] exit")'
# endif
  end function dthor_sub


  type(dthor) function dthor_mul_d(X,a) result(Y)
  use cublas_aux, only: cublas_scal
  class(dthor),     intent(IN) :: X
  double precision, intent(IN) :: a
# ifdef VERBOSE4
     print '("[+][dthor_mul_d] entry")'
# endif
     Y = X
     !call cublas_scal(a, Y% core(Y%d)% d_A)
     call cublas_scal(a, Y% core(1)% d_A)
# ifdef VERBOSE4
     print '("[+][dthor_mul_d] exit")'
# endif
  end function dthor_mul_d


  type(dthor) function d_mul_dthor(a,X) result(Y)
  double precision, intent(IN) :: a
  class(dthor),     intent(IN) :: X
# ifdef VERBOSE4
     print '("[+][d_mul_dthor] entry")'
# endif
     Y = dthor_mul_d(X,a)
# ifdef VERBOSE4
     print '("[+][d_mul_dthor] exit")'
# endif
  end function d_mul_dthor


  type(dthor) function dthor_mul_i(X,n) result(Y)
  class(dthor), intent(IN) :: X
  integer,      intent(IN) :: n
# ifdef VERBOSE4
     print '("[+][dthor_mul_i] entry")'
# endif
     Y = dthor_mul_d(X,dble(n))
# ifdef VERBOSE4
     print '("[+][dthor_mul_i] exit")'
# endif
  end function dthor_mul_i


  type(dthor) function i_mul_dthor(n,X) result(Y)
  integer,      intent(IN) :: n
  class(dthor), intent(IN) :: X
# ifdef VERBOSE4
     print '("[+][i_mul_dthor] entry")'
# endif
     Y = dthor_mul_d(X,dble(n))
# ifdef VERBOSE4
     print '("[+][i_mul_dthor] exit")'
# endif
  end function i_mul_dthor


  type(dthor) function dthor_div_d(X,a) result(Y)
  class(dthor),     intent(IN) :: X
  double precision, intent(IN) :: a
# ifdef VERBOSE4
     print '("[+][dthor_div_d] entry")'
# endif
     Y = dthor_mul_d(X,1d0/a)
# ifdef VERBOSE4
     print '("[+][dthor_div_d] exit")'
# endif
  end function dthor_div_d


  type(dthor) function dthor_div_i(X,n) result(Y)
  class(dthor), intent(IN) :: X
  integer,      intent(IN) :: n
# ifdef VERBOSE4
     print '("[+][dthor_div_i] entry")'
# endif
     Y = dthor_mul_d(X,1d0/dble(n))
# ifdef VERBOSE4
     print '("[+][dthor_div_i] exit")'
# endif
  end function dthor_div_i


  !> DThor object is h-stacked if its every core is an h-stack
  logical function dthor_is_hstack(this)
  class(dthor), intent(IN) :: this
  !
  character(*),parameter :: subnam='[dthor_is_hstack]'
  integer :: k
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     dthor_is_hstack = .true.
     do k = 1, this% d
        dthor_is_hstack = dthor_is_hstack .and. this% core(k)% is_hstack()
     enddo
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end function dthor_is_hstack


  !> DThor object is v-stacked if its every core is a v-stack
  logical function dthor_is_vstack(this)
  class(dthor), intent(IN) :: this
  !
  character(*),parameter :: subnam='[dthor_is_vstack]'
  integer :: k
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     dthor_is_vstack = .true.
     do k = 1, this% d
        dthor_is_vstack = dthor_is_vstack .and. this% core(k)% is_vstack()
     enddo
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end function dthor_is_vstack


  !> Converts every core to hStack
  subroutine dthor_to_hstack(this, VRBZ_) !
  class(dthor),   intent(INOUT) :: this
  integer, optional, intent(IN) :: VRBZ_
  !
  character(*),parameter :: subnam='[dthor_to_hstack]'
  integer :: VRBZ, k
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
     do k = 1, this% d
        call this% core(k)% to_hstack(VRBZ_=VRBZ)
     enddo
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end subroutine dthor_to_hstack


  !> Converts every core to vStack
  subroutine dthor_to_vstack(this, VRBZ_) !
  class(dthor),   intent(INOUT) :: this
  integer, optional, intent(IN) :: VRBZ_
  !
  character(*),parameter :: subnam='[dthor_to_vstack]'
  integer :: VRBZ, k
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     VRBZ=0; if(present(VRBZ_)) VRBZ=VRBZ_
     do k = 1, this% d
        call this% core(k)% to_vstack(VRBZ_=VRBZ)
     enddo
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end subroutine dthor_to_vstack


  !> Full tensor recovery from the distributed TT representation
  !!
  !! Returns local part of the N-dimensional array in the flattened form,
  !! assumed split along the last dimension N_d
  !!
  function dthor_full(TT, VRBZ_) result(X)
  double precision, allocatable :: X(:)  !< flattened local part of X
  class(dthor),      intent(IN) :: TT    !< tensor in TT format
  integer, optional, intent(IN) :: VRBZ_
  !
  character(*), parameter  :: subnam='[dthor_full]'
  type(stack) :: W(2)
  integer     :: d, info, n, VRBZ
# ifdef VERBOSE4
     print '("[+]'//subnam//' entry")'
# endif
     VRBZ=0; if (present(VRBZ_)) VRBZ=VRBZ_
     d = TT%d
     ! 1. Check memory requirements
     ! 2. Convert TT to hStack
     ! 3. Multiply
     ! 4. Redistribute
# ifdef VERBOSE4
     print '("[+]'//subnam//' exit")'
# endif
  end function dthor_full

end module thor
