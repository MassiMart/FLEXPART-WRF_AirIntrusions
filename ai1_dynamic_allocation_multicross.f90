subroutine dynamic_allocation_multicross(multicross,intrusion_type)

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
  integer :: dsize,multicross,aux_index2
  real, allocatable, dimension (:,:) :: temp_real_2d
  real(kind=dp), allocatable,dimension(:,:) :: temp_real_2d_kind
  integer, allocatable, dimension (:,:) :: temp_int_2d
  integer :: intrusion_type

  if(intrusion_type == 1) then
     dsize = size(tl_particle_id)

     colnum_tl = multicross + 5
     aux_index2 = multicross

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
     tl_entering_time(:,aux_index2:colnum_tl)       = -999
     tl_residence_time(:,aux_index2:colnum_tl)      = -999
     tl_flux_flag(:,aux_index2:colnum_tl)           = -999
     tl_xcross(:,aux_index2:colnum_tl)              = -999
     tl_ycross(:,aux_index2:colnum_tl)              = -999
     tl_zcross(:,aux_index2:colnum_tl)              = -999

     pvtraj_tl(:,aux_index2:colnum_tl)     = -999
     wtraj_tl(:,aux_index2:colnum_tl)      = -999
     tketraj_tl(:,aux_index2:colnum_tl)    = -999
     pptraj_tl(:,aux_index2:colnum_tl)     = -999
     ritraj_tl(:,aux_index2:colnum_tl)     = -999
     eitraj_tl(:,aux_index2:colnum_tl)     = -999
     bvfsqtraj_tl(:,aux_index2:colnum_tl)  = -999
     topo_cross_tl(:,aux_index2:colnum_tl) = -999

  else if(intrusion_type == 2) then
     dsize = size(strato_particle_id)

     colnum_st = multicross + 5
     aux_index2 = multicross

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

     strato_entering_time(:,aux_index2:colnum_st)   = -999
     strato_residence_time(:,aux_index2:colnum_st)  = -999
     strato_flux_flag(:,aux_index2:colnum_st)       = -999
     strato_xcross(:,aux_index2:colnum_st)          = -999
     strato_ycross(:,aux_index2:colnum_st)          = -999
     strato_zcross(:,aux_index2:colnum_st)          = -999

     pvtraj_st(:,aux_index2:colnum_st)     = -999
     wtraj_st(:,aux_index2:colnum_st)      = -999
     tketraj_st(:,aux_index2:colnum_st)    = -999
     pptraj_st(:,aux_index2:colnum_st)     = -999
     ritraj_st(:,aux_index2:colnum_st)     = -999
     eitraj_st(:,aux_index2:colnum_st)     = -999
     bvfsqtraj_st(:,aux_index2:colnum_st) = -999
     topo_cross_st(:,aux_index2:colnum_st) = -999

  end if
  
end subroutine dynamic_allocation_multicross

