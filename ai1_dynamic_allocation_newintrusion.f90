subroutine dynamic_allocation_newintrusion()

  !*****************************************************************************
  !     This file is part of the AirIntrusion Package.                         *
  !     This subroutine is an helper subroutine to dynamically allocate        *
  !     arrays.                                                                *
  !                                                                            *
  !     Author: Massimo Martina, PhD Student (MFF UK, Prague)                  *
  !             massimo.martina@matfyz.cuni.cz                                 *
  !                                                                            *
  !     20 June 2025                                                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use ai1_intrus_mod
  use par_mod, only: dp

  implicit none
  !VARIABLES
  !*********
  integer :: dsize,aux_index
  integer, allocatable, dimension(:) :: temp_int
  real, allocatable, dimension(:) :: temp_real
  real, allocatable, dimension (:,:) :: temp_real_2d
  real(kind=dp),allocatable, dimension (:,:) :: temp_real_2d_kind
  integer, allocatable, dimension (:,:) :: temp_int_2d

  !DYNAMIC ALLOCATION OF NEEDED ARRAY FOR THE INTRUSION IDENTIFICATION PROCESS
  !***************************************************************************
  !NOTE_1: the array must be intially defined as allocatable and then allocate
  !        as an array with 1 element to allow the dynamic allocation.
  !NOTE_2: the array "free_entering_time" must be defined for each parcels since
  !        we do not know a priori which parcel will enter the tropopause/stratosphere;
  !        hence, it must be allocated considering the max number of particles before entering
  !        this subroutine.

  !Check the last element of tl_particle_id, if it is different from the
  !initial value (-999) it means that the end of the array has been reached
  !(also for the other arrays)
!  write(*,*) "size tl_particle_id", size(tl_particle_id),"allocated", allocated(tl_particle_id)
!  write(*,*) "tl_particle_id", tl_particle_id
!  write(*,*) "*********************************************************************************"

!  stop
  if(tl_particle_id(size(tl_particle_id)) /= -999) then
     write(*,*) "I am doing dynamic allocation tl - dynamic_allocation_newintrusion"
     !Increase the size of the array by 100 spaces       
     dsize = size(tl_particle_id) + 1000
     aux_index = size(tl_particle_id) + 1 !index for the first new element of the arrays

     !TROPOPAUSE LAYER
!     allocate(temp_int(dsize))
!     temp_int(:dsize) = tl_count_flag
!     call move_alloc(temp_int, tl_count_flag)

     allocate(temp_int(dsize))
     temp_int(:dsize) = tl_particle_id
     call move_alloc(temp_int, tl_particle_id)

     allocate(temp_real(dsize))
     temp_real(:dsize) = tl_transition_time
     call move_alloc(temp_real, tl_transition_time)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = tl_entering_time
     call move_alloc(temp_real_2d, tl_entering_time)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = tl_residence_time
     call move_alloc(temp_real_2d, tl_residence_time)

     allocate(temp_int_2d(dsize,colnum_tl))
     temp_int_2d(:dsize,:) = tl_flux_flag
     call move_alloc(temp_int_2d, tl_flux_flag)

     allocate(temp_real_2d_kind(dsize,colnum_tl))
     temp_real_2d_kind(:dsize,:) = tl_xcross
     call move_alloc(temp_real_2d_kind, tl_xcross)

     allocate(temp_real_2d_kind(dsize,colnum_tl))
     temp_real_2d_kind(:dsize,:) = tl_ycross   
     call move_alloc(temp_real_2d_kind, tl_ycross)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = tl_zcross
     call move_alloc(temp_real_2d, tl_zcross)

     !THERMODYNAMIC Variables
     !TL
     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = pvtraj_tl
     call move_alloc(temp_real_2d,pvtraj_tl)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = wtraj_tl
     call move_alloc(temp_real_2d,wtraj_tl)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = tketraj_tl
     call move_alloc(temp_real_2d,tketraj_tl)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = pptraj_tl
     call move_alloc(temp_real_2d,pptraj_tl)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = ritraj_tl
     call move_alloc(temp_real_2d,ritraj_tl)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = eitraj_tl
     call move_alloc(temp_real_2d,eitraj_tl)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = bvfsqtraj_tl
     call move_alloc(temp_real_2d,bvfsqtraj_tl)

     allocate(temp_real_2d(dsize,colnum_tl))
     temp_real_2d(:dsize,:) = topo_cross_tl
     call move_alloc(temp_real_2d,topo_cross_tl)

     !Initialization of the values
     tl_particle_id(aux_index:dsize)           = -999
     tl_entering_time(aux_index:dsize,:)       = -999
     tl_transition_time(aux_index:dsize)       = -999
     tl_residence_time(aux_index:dsize,:)      = -999
     tl_flux_flag(aux_index:dsize,:)           = -999
     tl_xcross(aux_index:dsize,:)              = -999
     tl_ycross(aux_index:dsize,:)              = -999
     tl_zcross(aux_index:dsize,:)              = -999

     pvtraj_tl(aux_index:dsize,:)     = -999
     wtraj_tl(aux_index:dsize,:)      = -999
     tketraj_tl(aux_index:dsize,:)    = -999
     pptraj_tl(aux_index:dsize,:)     = -999
     ritraj_tl(aux_index:dsize,:)     = -999
     eitraj_tl(aux_index:dsize,:)     = -999
     bvfsqtraj_tl(aux_index:dsize,:) = -999
     topo_cross_tl(aux_index:dsize,:) = -999


  end if

  !Check the last element of strato_particle_id, if it is different from the
  !initial value (-999) it means that the end of the array has been reached
  !(also for the other arrays)
  write(*,*) strato_particle_id(size(strato_particle_id))
  if(strato_particle_id(size(strato_particle_id)) /= -999) then
     write(*,*) "I am doing dynamic allocation STRATO - dynamic_allocation_newintrusion"
     !Increase the size of the array by 100 spaces       
     dsize = size(strato_particle_id) + 100
     aux_index = size(strato_particle_id) + 1 !index for the first new element of the arrays

     !STRATOSPHERE
     allocate(temp_int(dsize))
     temp_int(:dsize) = strato_particle_id
     call move_alloc(temp_int, strato_particle_id)

     allocate(temp_real(dsize))
     temp_real(:dsize) = strato_transition_time
     call move_alloc(temp_real, strato_transition_time)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = strato_entering_time
     call move_alloc(temp_real_2d, strato_entering_time)

     allocate(temp_int_2d(dsize,colnum_st))
     temp_int_2d(:dsize,:) = strato_flux_flag
     call move_alloc(temp_int_2d, strato_flux_flag)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = strato_residence_time
     call move_alloc(temp_real_2d, strato_residence_time)

     allocate(temp_real_2d_kind(dsize,colnum_st))
     temp_real_2d_kind(:dsize,:) = strato_xcross
     call move_alloc(temp_real_2d_kind, strato_xcross)

     allocate(temp_real_2d_kind(dsize,colnum_st))
     temp_real_2d_kind(:dsize,:) = strato_ycross
     call move_alloc(temp_real_2d_kind, strato_ycross)
 
     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = strato_zcross
     call move_alloc(temp_real_2d, strato_zcross)

     !THERMODYNAMIC Variables
     !STRATO
     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = pvtraj_st
     call move_alloc(temp_real_2d,pvtraj_st)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = wtraj_st
     call move_alloc(temp_real_2d,wtraj_st)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = tketraj_st
     call move_alloc(temp_real_2d,tketraj_st)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = pptraj_st
     call move_alloc(temp_real_2d,pptraj_st)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = ritraj_st
     call move_alloc(temp_real_2d,ritraj_st)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = eitraj_st
     call move_alloc(temp_real_2d,eitraj_st)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = bvfsqtraj_st
     call move_alloc(temp_real_2d,bvfsqtraj_st)

     allocate(temp_real_2d(dsize,colnum_st))
     temp_real_2d(:dsize,:) = topo_cross_st
     call move_alloc(temp_real_2d,topo_cross_st)

     !Initialization of the values
     strato_particle_id(aux_index:dsize)       = -999
     strato_entering_time(aux_index:dsize,:)   = -999
     strato_transition_time(aux_index:dsize)   = -999
     strato_residence_time(aux_index:dsize,:)  = -999
     strato_flux_flag(aux_index:dsize,:)       = -999
     strato_xcross(aux_index:dsize,:)          = -999
     strato_ycross(aux_index:dsize,:)          = -999
     strato_zcross(aux_index:dsize,:)          = -999

     pvtraj_st(aux_index:dsize,:)     = -999
     wtraj_st(aux_index:dsize,:)      = -999
     tketraj_st(aux_index:dsize,:)    = -999
     pptraj_st(aux_index:dsize,:)     = -999
     ritraj_st(aux_index:dsize,:)     = -999
     eitraj_st(aux_index:dsize,:)     = -999
     bvfsqtraj_st(aux_index:dsize,:) = -999
     topo_cross_st(aux_index:dsize,:) = -999

  end if

end subroutine dynamic_allocation_newintrusion
