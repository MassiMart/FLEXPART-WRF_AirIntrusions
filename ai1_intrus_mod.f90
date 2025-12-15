!***********************************************************************
!* Copyright 2025                                                      *
!*                                                                     *
!* This file is part of FLEXPART-WRF package "AirIntrusions".      *
!*                                                                     *
!***********************************************************************
!*******************************************************************************
!   Definition of global variables used for the identification of those        *
!   air parcels that enter the tropopause/strato from the Planetary Boundary   *
!   Layer (PBL)/Free Atmosphere (FA) (or the other way round).                 *
!                                                                              *
!   Author: Massimo Martina PhD Student (MFF UK, Prague)                       *
!           massimo.martina@matfyz.cuni.cz                                     *
!                                                                              *
!   19 June 2025                                                               *
!                                                                              *
!  Last modifications: 19 June 2025                                            *
!                                                                              *
!*******************************************************************************
module ai1_intrus_mod


  use par_mod, only: maxnests,nxmax,nymax,nuvzmax,dp,nxmaxn,nymaxn,nzmax

  implicit none

  !****************************************************************
  !-----------Varibles for the Intrusions Subroutines-------------*
  !****************************************************************
  !PARAMETERS
  integer, parameter :: maxresidenceclass = 10
  integer, parameter :: unit_res = 60, unit_mass = 65 !, unit_deb = 67, unit_deb_flux=68, unitdebug2=69
  integer, parameter :: expected_crossing = 1 !maximum expected number of intrusion into the TL/Strato for one parcel

  !readinput.f90
  integer :: nresidenceclass, lresidence(maxresidenceclass),tropo_method, source_type
  integer :: cape_option
  !Calcpar.f90 & ai1_boundaries_computation.f90
  real :: cold_tropopause(0:nxmax-1,0:nymax-1,1,2) 
  real :: cold_tropopausen(0:nxmax-1,0:nymax-1,1,2,maxnests) 

  !ai1_thermodyn_computation.f90
  !These variables are 2D arrays:
  ! - first dim = it represents the number of parcels entering the Tl/Strato
  ! - second dim  = it represents the number of intrusions into the Tl/Strato for each parcel
  !Tropopause Layer
  real,allocatable,dimension(:,:) :: pvtraj_tl,wtraj_tl,tketraj_tl,pptraj_tl,ritraj_tl,eitraj_tl,&
                                       bvfsqtraj_tl,topo_cross_tl
  !Stratosphere
  real,allocatable,dimension(:,:) :: pvtraj_st,wtraj_st,tketraj_st,pptraj_st,ritraj_st,eitraj_st,&
                                       bvfsqtraj_st,topo_cross_st

  !ai1_up_intrusion_identifier
  integer,allocatable, dimension (:)       :: npoint_id
  integer,allocatable, dimension (:)       :: tl_particle_id,strato_particle_id  !MM_11-10-2024
  real,allocatable,dimension (:)           :: free_atm_entering_time,&
                                              tl_transition_time, strato_transition_time
  integer,allocatable,dimension (:)        :: tl_count_flag, strato_count_flag
  real,allocatable,dimension (:,:)         :: tl_entering_time     , strato_entering_time,&
                                              tl_residence_time    , strato_residence_time
  integer,allocatable,dimension (:,:)      :: tl_flux_flag         , strato_flux_flag
  real(kind=dp),allocatable,dimension (:,:):: tl_xcross            , strato_xcross, &
                                              tl_ycross            , strato_ycross
  real, allocatable,dimension(:,:)         :: tl_zcross            , strato_zcross
  integer :: colnum_tl,colnum_st
  !ai1_readwind.f90
  real :: capeh(0:nxmax-1,0:nymax-1,nuvzmax,2)
  !ai1_verttransform.f90
  real :: cape(0:nxmax-1,0:nymax-1,nzmax,2)
  !ai1_readwinf_nests.f90
  real :: capehn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
  !ai1_verttransform_nests.f90
  real :: capen(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)

  !ai1_calccrossing.f90
  real,allocatable, dimension (:,:,:,:,:)      :: flux_mass !MM_11-03-2025
  integer,allocatable, dimension (:,:,:,:)     :: flux_crossing !MM_10-10-2024
  integer,allocatable, dimension (:,:,:,:)     :: sum_entering_time !MM_10-10-2024
  integer,allocatable, dimension (:,:,:,:)     :: sum_transition_time !MM_10-10-2024
  real,allocatable,dimension(:,:,:,:)          :: sum_pvtraj,sum_wtraj,sum_tketraj,&
                                                  sum_pptraj,sum_ritraj,sum_eitraj,&
                                                  sum_bvfsqtraj

end module ai1_intrus_mod
