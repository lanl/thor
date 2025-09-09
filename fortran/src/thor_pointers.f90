!--------------------------------------------------------------------------~*
!! Copyright (c) 2025 Triad National Security, LLC
!! All rights reserved.
!!--------------------------------------------------------------------------~*/

!!
!! @file thor_pointers.f90
!! @author Ismael Djibrilla Boureima [IDB] (T-3)
!! @date   04/21/2023
!! @brief  Library of pointers for host 1D, 2D, 3D, and 4D array
!!
module thor_pointers
 implicit none
 
 !
 ! HOST/CPU ARRAY POINTERS
 !
 !             INTEGERS
   type,public :: pointi
     integer,dimension(:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointi2
     integer,dimension(:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointi3
    integer,dimension(:,:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointi4
    integer,dimension(:,:,:,:),pointer,contiguous :: p=>null()
   end type
 !             SINGLE PRECISION
   type,public :: point
     real,dimension(:),pointer,contiguous :: p=>null()
   end type
   type,public :: point2
     real,dimension(:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: point3
     real,dimension(:,:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: point4
     real,dimension(:,:,:,:),pointer,contiguous :: p=>null()
   end type
 !             DOUBLE PRECISION
   type,public :: pointd
     double precision,dimension(:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointd2
     double precision,dimension(:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointd3
     double precision,dimension(:,:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointd4
     double precision,dimension(:,:,:,:),pointer,contiguous :: p=>null()
   end type
!             COMPLEX 
   type,public :: pointz
     double complex,dimension(:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointz2
     double complex,dimension(:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointz3
     double complex,dimension(:,:,:),pointer,contiguous :: p=>null()
   end type
   type,public :: pointz4
     double complex,dimension(:,:,:,:),pointer,contiguous :: p=>null()
   end type

end module thor_pointers
