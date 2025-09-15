module time_lib
!use omp_lib
implicit none
real,save:: system_dt

integer,private,parameter :: max_counters = 200     !< counters stack depth
integer,private,save      :: counters(max_counters) !< counters stack
integer,private,save      :: sp                     !< counters stack pointer
integer,private,save      :: cr                     !< counts rate
integer,private,save      :: cm                     !< counter max value
real,private,save         :: rate                   !< counts rate (real prec.)

contains

 double precision function timef( )
 implicit none
 !include 'mpif.h'
 integer :: c,r,m
 real*4 :: t(2),time
 real*8 :: mpitime
 real*8,external :: mpi_wtime
    timef=1.d0*mpi_wtime()
 end function


 subroutine system_timer_init()
 call system_clock(count_rate=cr)
 call system_clock(count_max=cm)
    rate = real(cr)
    sp = 0 !< install stack pointer to 0 (empty stack)
 end subroutine


 subroutine system_timer_start()
    system_dt= 0.0
    if (sp.lt.max_counters) sp= sp + 1
    call system_clock(counters(sp))
 end subroutine


 subroutine system_timer_stop()
 integer :: c2
    call system_clock(c2)
    system_dt= real(c2 - counters(sp))/rate
    if (system_dt.lt.0.0) system_dt = system_dt + cm/rate
    if (sp.gt.0) sp= sp - 1
 end subroutine

end module 
