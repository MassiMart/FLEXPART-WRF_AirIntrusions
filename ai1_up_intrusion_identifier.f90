
subroutine up_intrusion_identifier(jtime,jpart,xold,yold,zold,ll_index,mm_index)
  !*****************************************************************************
  !     This file is part of the AirIntrusion Package                          *
  !     This subroutine performs the identification of those parcels entering  *
  !     the tropopause/stratosphere. It also computed some thermodynamic       *
  !     variables at the particle position before it enter the tropopause/     *
  !     stratosphere.                                                          *
  !                                                                            *
  !     Additional required subroutines:                                       *
  !       1) ai1_boundaries_computation.f90 for the computation of the         *
  !          tropopause boundaries. These can be computed considering          *
  !          the Lapse Rate tropopause at its bottom and the Cold Point        *
  !          tropopause at its top (tropo_method == 1), or using the           *
  !          climatological definition by Fueglistaler et al.(2009)            * 
  !          (tropo_method == 2).                                              *
  !       2) ai1_thermodyn_computation.f90 for the computation of thermodynamic*
  !          variables (such as Richardson Number, ...).                       *
  !
  !     Author: Massimo Martina, PhD Student (MFF UK, Prague)                  *
  !             massimo.martina@matfyz.cuni.cz                                 *
  !                                                                            *
  !     09 January 2025                                                        *
  !     Last modifications:                                                    *
  !     20 June 2025 substantial reorganization of the code; introduction of   *
  !                  dynamic allocation for certain arrays.                    *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! jtime            (IN) Actual temporal position of calculation              *
  ! jpart            (IN) Index of the particle considered                     *
  ! xold,yold,zold   (IN) "Memorized" old positions of the particle            *
  ! tropo_method     (IN) Index for choosing the tropopause boundaries         *
  !                       computation method                                   *
  ! ll_index         (IN) Index used to access those array associated to       *
  !                       intrusions in the tropopause                         *
  ! mm_index         (IN) Index used to access those array associated to       *
  !                       intrusions in the stratosphere                       *             
  !*****************************************************************************

  use par_mod
  use com_mod
  use ai1_intrus_mod

  implicit none
  !VARIABLES
  !******************************
  !Input
  integer           :: jtime, jpart, newtime
  real(kind=dp)     :: xold, yold
  real              :: zold
 !Internal variables 
  real(kind=dp)     :: jul_rel, jul_cross, jul_pbl                               !MM_14-06-2024
  integer           :: jjjjmmdd_rel, ihmmss_rel, jjjjmmdd_cross, ihmmss_cross,&  !MM_14-06-2024
                       jjjjmmdd_pbl, ihmmss_pbl
  character         :: adate_rel*8,atime_rel*6,adate_cross*8,atime_cross*6,&     !MM_14-06-2024
                       adate_pbl*8,atime_pbl*6
  character(len=20) :: a_rel, a_cross, a_pbl                                     !MM_14-06-2024
  integer           :: il, ind                                                   !MM_18-07-2024
  real(kind=dp)     :: jul_traj,jul_term_rel,jul_term                            !MM_18-07-2024
  integer           :: jjjjmmdd_traj, ihmmss_traj                                !MM_18-07-2024
  character         :: adate_traj*8,atime_traj*6                                 !MM_18-07-2024
  character(len=20) :: time_traj
  integer           :: ll_index, mm_index, ll_index_saved, mm_index_saved
  integer           :: kk, count_index                                     
  integer           :: dsize,aux_index
  real              :: tri_old,tri,cold_tropo_old,cold_tropo,pp_old,pp_new,hmixi,hmixi_old !BUG: PREVIOUS VERSIONS DEFINE THESE VARIABLES ONLY IN BOUNDARIES_COMPUTATION.F90  
  real              :: aux_tropo
  integer, allocatable, dimension(:) :: temp_int
  real, allocatable, dimension(:) :: temp_real
  real, allocatable, dimension (:,:) :: temp_real_2d
  integer, allocatable, dimension (:,:) :: temp_int_2d
  integer :: intrusion_type 

  integer :: jjjjmmdd,ihmmss
  real(kind=dp) :: jul
  character :: adate*8,atime*6
  !END VARIABLES
  !*****************************

  !DYNAMIC ALLOCATION OF NEEDED ARRAY FOR THE INTRUSION IDENTIFICATION PROCESS
  !***************************************************************************
!  if(colnum_tl .gt. expected_crossing) then
!    write(*,*) "colnum_tl = ", colnum_tl
!  end if

!  if(colnum_st .gt. expected_crossing) then
!    write(*,*) "colnum_st = ", colnum_st
!  end if
!  write(*,*) "I AM IN UP INTRUSION IDENTIFIER"
!  write(*,*) "size tl_particle_id", size(tl_particle_id),"allocated", allocated(tl_particle_id)
!  write(*,*) "tl_particle_id", tl_particle_id
!  write(*,*) "*********************************************************************************"
!  if(jpart == 4075 ) then
!    write(*,*) "npoint_id(jpart)",npoint_id(jpart)
!    write(*,*) "size(tl_particle_id)", size(tl_particle_id)
!    write(*,*) tl_particle_id    
!  end if
!  write(*,*) "I AM IN UP_INTRUSION_IDENTIFIER"
!  write(*,*) tl_particle_id(size(tl_particle_id))
!  write(*,*) strato_particle_id(size(strato_particle_id))
  if(tl_particle_id(size(tl_particle_id)) /= -999 .or.&
     strato_particle_id(size(strato_particle_id)) /= -999) then
     call dynamic_allocation_newintrusion()
  end if

  !*********************************************************************************!
  !                                                                                 !
  !----------METHOD 1: TROPOPAUSE DEFINED USING LRT & CPT TROPOPAUSE                !
  !----------METHOD 2: TROPOPAUSE DEFINED AS IN FUEGLISTALER ET AL. (2009)          !
  !*********************************************************************************! 
  !BOUNDARIES COMPUTATION (PBL,TROPOPAUSE)
  !***************************************
  !Old position
!  call boundaries_computation(jtime,jpart,xold,yold,zold,&
!                              tri_old,cold_tropo_old,&
!                              pp_old,hmixi_old)
  !New position
!  newtime = jtime + lsynctime
!  call boundaries_computation(newtime,jpart,xtra1(jpart),ytra1(jpart),&
!                              ztra1(jpart),tri,cold_tropo,&
!                              pp_new,hmixi)


  !DEBUG !HERE tri_old,... will be alway equal to 0!!!
!  if(jpart == 107928) then ! .and. jtime>40860 ) then
!    write(*,*) "xold",xold,"yold",yold,"zold",zold,"ztra1",ztra1(jpart)
!    write(*,*) "xtra1",xtra1(jpart),"ytra1",ytra1(jpart)
!    write(*,*) "tri_old",tri_old,"cold_tropo_old",cold_tropo_old
!    write(*,*) "tri",tri,"cold_tropo",cold_tropo
!    write(*,*) "itra1(jpart)",itra1(jpart)
!    write(*,*) "jpart", jpart
!  end if   

  !Considering particles that are still in the domain
  if(itra1(jpart) > 0) then !considering particles that are still in the domain
     !BOUNDARIES COMPUTATION (PBL,TROPOPAUSE)
     !***************************************
    ! write(*,*) "I AM BEFORE BOUNDARIES COMPUTATION"
     !Old position
     call boundaries_computation(jtime,jpart,xold,yold,zold,&
                                 tri_old,cold_tropo_old,&
                                 pp_old,hmixi_old)
     !New position
     newtime = jtime + lsynctime
     call boundaries_computation(newtime,jpart,xtra1(jpart),ytra1(jpart),&
                                 ztra1(jpart),tri,cold_tropo,&
                                 pp_new,hmixi)

    !write(*,*) "I AM AFTER BOUNDARIES COMPUTATION"

     !CHECK ON THE TROPOPAUSE BOUNDARIES
!     write(*,*) " I am in intrusion_identifier in the if(itra1(jpart) > 0"
     if(tri_old .gt. cold_tropo_old) then
        write(*,*) "Lapse rate tropopause greater than Cold Point Tropopause"
        write(*,*) "tri_old", tri_old, "cold_tropo_old", cold_tropo_old
        !inversion of the values
        aux_tropo = tri_old
        tri_old   = cold_tropo_old
        cold_tropo_old = aux_tropo
     end if
      !CHECK ON THE TROPOPAUSE BOUNDARIES
     if(tri .gt. cold_tropo) then
        write(*,*) "Lapse rate tropopause greater than Cold Point Tropopause"
        write(*,*) "tri", tri, "cold_tropo", cold_tropo
        !inversion of the values
        aux_tropo = tri
        tri   = cold_tropo
        cold_tropo = aux_tropo
     end if
       
     !INTRUSION IDENTIFICATION - SOURCES LOCATED IN THE PBL
     !*****************************************************
     !PARTICLE ENTER THE FREE ATMOSPHERE (CROSSING THE PBL TOP) 
     if(zold .le. hmixi_old .AND. ztra1(jpart) .gt. hmixi) then
          !Time at which the particle cross the PBL TOP
         free_atm_entering_time(jpart) = jtime + lsynctime/2. 
     end if !End if particle cross PBL TOP

      !PARTICLE ENTER THE TROPOPAUSE LAYER (TL)
     if(tropo_method == 1 .AND. (zold .lt. tri_old .AND. (ztra1(jpart) .ge. tri .AND. ztra1(jpart) .le. cold_tropo)) .OR.& !METHOD 1
        tropo_method == 1 .AND. (zold .gt. cold_tropo_old .AND. (ztra1(jpart) .le. cold_tropo .AND. ztra1(jpart) .ge. tri)) .OR.& !PAR entering TL from strato   
        tropo_method == 2 .AND. (pp_old .gt. 15000 .AND. (pp_new .ge. 7000 .AND. pp_new .le. 15000)) .OR.&               !METHOD 2
        tropo_method == 2 .AND. (pp_old .lt. 7000 .AND. (pp_new .ge. 7000 .AND. pp_new .le. 15000))) then                

        intrusion_type = 1
        !Counting the # of time the particle enters the TTL
     !   write(*,*) "I AM IN PARTICLE ENTER THE TROPOPAUSE LAYER"
      !  write(*,*) "jpart", jpart
        tl_count_flag(jpart)        = tl_count_flag(jpart) + 1  
        if(tl_count_flag(jpart) .gt. colnum_tl) then
           write(*,*) "Increasing second dimension of the arrays involved in the identification process"
           write(*,*) "Previous values of colnum_tl = ", colnum_tl
           write(*,*) "tl_count_flag(jpart)", tl_count_flag(jpart)
           call dynamic_allocation_multicross(tl_count_flag(jpart),intrusion_type)
           write(*,*) "New value of colnum_tl =", colnum_tl
        end if

      !  write(*,*) "tl_count_flag(jpart)", tl_count_flag(jpart)
        !Assign array index to parcel that enters in the TL for the first time
        if(tl_count_flag(jpart) == 1) then
           ll_index = ll_index + 1
           !Saving the index of the last accessed position of the array
           ll_index_saved = ll_index
           tl_particle_id(ll_index) = npoint_id(jpart)
           !tl_particle_count(ll_index) = tl_count_flag(jpart)
        !Searching for the array index associated to particles that are entering the TL
        !a second (or more) time
        else
            kk = 0
            !Do cycle over all the element of particle_id
            do 
              kk = kk + 1 
              !Check that the particle was already present in tl_particle_id, if not stop.                                                 
              if(kk .gt. size(tl_particle_id)) then
                  write(*,*) "kk index is greater than tl_particle_id"
                  write(*,*) "kk: ", kk, "jpart", jpart
                  write(*,*) "npoint_id(jpart)", npoint_id(jpart)
                  write(*,*) "size(tl_particle_id)", size(tl_particle_id)
                  write(*,*) tl_particle_id
                  write(*,*) "I AM IN PARTICLE ENTER THE TTL"
                  stop
              end if
               if(tl_particle_id(kk) == npoint_id(jpart)) then
                 !Saving the index of the last accessed position of the array
                 ll_index_saved = ll_index
                 !Index pointing to the position of the particle in the array
                 ll_index = kk
                 !tl_particle_count(ll_index) = tl_count_flag(jpart)
                 exit
              end if
            end do
        end if
        !Time at which the particle enters the TL
        tl_entering_time(ll_index,tl_count_flag(jpart)) = jtime + lsynctime/2. 
        !Flux towards the TL from the troposphere
        tl_flux_flag(ll_index,tl_count_flag(jpart)) = +1 
         
        !For the first intrusion the transition time is computed                
        if(tl_count_flag(jpart) == 1) then 
                  
           tl_transition_time(ll_index) = tl_entering_time(ll_index,tl_count_flag(jpart)) -&
                                          free_atm_entering_time(jpart)
           !Check on the transition time, if it is less than the internal time step, stop.  
           if(tl_transition_time(ll_index) < 180) then
              write(*,*) "transition time is less than internal time step: ",&
                          tl_transition_time(ll_index)
              write(*,*) "tl_entering_time", tl_entering_time(ll_index,tl_count_flag(jpart))
              write(*,*) "free_atm_entering_time", free_atm_entering_time(jpart)
              write(*,*) "ll_index", ll_index, "jpart", jpart
              stop
           end if
        end if
        !Save position of the parcel before entering the TL
        tl_xcross(ll_index,tl_count_flag(jpart)) = (xold)
        tl_ycross(ll_index,tl_count_flag(jpart)) = (yold)
        tl_zcross(ll_index,tl_count_flag(jpart)) = (zold)

        !COMPUTATION OF THE THERMODYNAMIC VARIABLES
        call thermodyn_computation(jtime,jpart,xold,yold,zold,&
                                   ll_index,tl_count_flag(jpart),&
                                   intrusion_type )
            
        !Restoring the array index to the last accessed position
        ll_index = ll_index_saved


     end if !END PARTICLE ENTER THE TL

     !PARTICLE ENTER THE STRATOSPHERE
     if(tropo_method == 1 .AND. (zold .le. cold_tropo_old .AND. ztra1(jpart) .gt. cold_tropo) .OR. & !METHOD 1
        tropo_method == 2 .AND. (pp_old .ge. 7000 .AND. pp_new .lt. 7000)) then                       !METHOD 2

        intrusion_type = 2
        !Counting the # of time the particle enters the strato
        strato_count_flag(jpart) = strato_count_flag(jpart) + 1

        if(strato_count_flag(jpart) .gt. colnum_st) then
           write(*,*) "Increasing second dimension of the arrays involved in the identification process"
           write(*,*) "Previous values of colnum_st = ", colnum_st
           call dynamic_allocation_multicross(strato_count_flag(jpart),intrusion_type)
           write(*,*) "New value of colnum_st =", colnum_st
        end if
 
        !Assign array index to parcel that enters in the strato for the first time
        if(strato_count_flag(jpart) == 1) then
           mm_index = mm_index + 1
           strato_particle_id(mm_index) = npoint_id(jpart)
           !strato_particle_count(mm_index) = strato_count_flag(jpart)
           !Saving the index of the last accessed position of the array
           mm_index_saved = mm_index
        !Searching for the array index associated to particles that are entering the TL
        !a second (or more) time
        else
            kk = 0
            !Do cycle over all the element of particle_id
            do
              kk = kk + 1
              !Check that the particle was already present in tl_particle_id, if not stop.
              if(kk .gt. size(strato_particle_id)) then
                 write(*,*) "kk index is greater than strato_particle_id"
                 write(*,*) "kk: ", kk
                 write(*,*) "strato_particle_id(jpart): ", strato_particle_id(kk)
                 write(*,*) "I AM IN PARTICLE ENTER THE STRATOSPHERE"
                 stop
              end if

              if(strato_particle_id(kk) == npoint_id(jpart)) then
                 !Saving the index of the last accessed position of the array
                 mm_index_saved = mm_index
                 !Index pointing to the position of the particle in the array
                 mm_index = kk
                 !strato_particle_count(mm_index) = strato_count_flag(jpart)
                 exit
              end if
            end do
        end if

        !Time at which the particle enters the srato
        strato_entering_time(mm_index,strato_count_flag(jpart)) = jtime + lsynctime/2. 
!        !Flux towards the strato from the troposphere
        strato_flux_flag(mm_index,strato_count_flag(jpart)) = +1
                  
        !For the first intrusion the transition time is computed
        if(strato_count_flag(jpart) == 1) then 
               strato_transition_time(mm_index) = strato_entering_time(mm_index,strato_count_flag(jpart)) -&
                                                  free_atm_entering_time(jpart)
            !Check on the transition time, if it is less than the internal time step, stop.
            if(strato_transition_time(mm_index) < 180) then
               write(*,*) "strato transition time is less then internal time step: ", &
               strato_transition_time(mm_index)
               stop
            end if
        end if
             
        !If the particle was for a crtain time in the ttl, its residence time there must be computed.
        if(tl_count_flag(jpart) > 0) then 
           kk = 0
           !Loop over all the previous elements of particle_id(l)
           do 
             kk = kk + 1
             !Check that the particle was already present in tl_particle_id, if not stop.
             if(kk .gt. size(tl_particle_id)) then
                write(*,*) "kk index is greater than tl_particle_id"
                write(*,*) "kk: ", kk, "size tl_particle_id", size(tl_particle_id)
                write(*,*) "tl_particle_id(kk): ", tl_particle_id(kk)
                write(*,*) "npoint_id", npoint_id(jpart), "jpart",jpart
                write(*,*) "tl_count_flag", tl_count_flag(jpart)
                write(*,*) "tl_residence_time", tl_residence_time(kk,tl_count_flag(jpart))
                write(*,*) "strato_count_flag",strato_count_flag(jpart)
                write(*,*) "zold",zold,"tri_old",tri_old,"cold_tropo_old",cold_tropo_old
                write(*,*) "ztra1", ztra1(jpart),"tri",tri,"cold_tropo",cold_tropo
                write(*,*) "I AM IN PARTICLE ENTER THE STRATOSPHERE - COMPUTING TL RESIDENCE TIME"
                stop
             end if
             !Finding array index location of the particle
             if(tl_particle_id(kk) == npoint_id(jpart) .AND. tl_residence_time(kk,tl_count_flag(jpart))<0) then
                tl_residence_time(kk,tl_count_flag(jpart)) = jtime+lsynctime/2.-&
                                                             tl_entering_time(kk,tl_count_flag(jpart))
                !Check. If the residence time is less than the internal time step, stop.
                if(tl_residence_time(kk,tl_count_flag(jpart)) < 180) then
                   write(*,*) "tl_residence time is less than internal time step: ", &
                               tl_residence_time(kk,tl_count_flag(jpart))
                   stop
                end if
                exit
             end if
           end do
        end if
 
        !Save position of the parcel before entering the Strato
        strato_xcross(mm_index,strato_count_flag(jpart)) = (xold)
        strato_ycross(mm_index,strato_count_flag(jpart)) = (yold)
        strato_zcross(mm_index,strato_count_flag(jpart)) = (zold)

        !COMPUTING THERMODYNAMIC VARIABLES
        call thermodyn_computation(jtime,jpart,xold,yold,zold,&    
                                   mm_index,strato_count_flag(jpart),&
                                   intrusion_type)

        !Restoring the array index to the last accessed position
        mm_index = mm_index_saved 
     end if !END PARTICLE ENTER THE STRATOSPHERE


     !PARTICLE MOVE FROM TL OR STRATOSPHERE TO TROPOSPHERE
     if(tropo_method == 1 .AND. (((zold .le. cold_tropo_old .AND. zold .ge. tri_old) .AND. ztra1(jpart) .lt. tri) .OR. & !METHOD 1
                                  (zold .gt. cold_tropo_old .AND. ztra1(jpart) .le. cold_tropo)) .OR. &  !MM_22-102025 ERROR ztra1 must be compered to cold_tropo, not tri!                          
        tropo_method == 2 .AND. (((pp_old .ge. 7000 .AND. pp_old .le. 15000) .AND. pp_new .gt. 15000) .OR.&              !METHOD 2
                                  (pp_old .lt. 7000 .AND. pp_new .gt. 7000))) then  !MM_22-102025 ERROR ztra1 must be compered to cold_tropo (70 hPa), not tri (150 hPa)!
        !Particle move from TL to Tropo
        if(tropo_method == 1 .AND. (zold .le. cold_tropo_old .AND. zold .ge. tri_old) .OR. &!METHOD 1
           tropo_method == 2 .AND. (pp_old .ge. 7000 .AND. pp_old .le. 15000)) then        !METHOD 2
           !Since the particle is exiting the ttl, it is possible to compute the residence time
           !The if statement consider only those particles entering the TL from the troposphere 
           !It does not consider the intrusion from the stratosphere for which tl_count_flag == strato_count_flag
           if(tl_count_flag(jpart) .ge.1) then
           !NOW THAT THE INTRUSION INTO THE TL FROM THE STRATO ARE CONSIDERED THE ABOVE IF STATEMENT IS NOT NEED IT ANYMORE
              kk = 0
              !Loop over all the previous elements of particle_id(l)
              do
                kk = kk + 1
                !Check that the particle was already present in tl_particle_id, if not stop
                if(kk .gt. size(tl_particle_id)) then
                   write(*,*) "kk index is greater than tl_particle_id"
                   write(*,*) "kk: ", kk
                   write(*,*) "tl_particle_id(jpart): ", tl_particle_id(jpart)
                   write(*,*) "I AM IN PARTICLE MOVING FROM TTL TO TROPOSPHERE"
                   stop
                end if
                   
                !Finding array index location of the particle
                if(tl_particle_id(kk) == npoint_id(jpart)) then
                   tl_residence_time(kk,tl_count_flag(jpart)) = jtime + lsynctime/2.-&
                                                                tl_entering_time(kk,tl_count_flag(jpart))
                   !Check. If residence time less than internal time step, stop.
                   if(tl_residence_time(kk,tl_count_flag(jpart)) < 180) then
                      write(*,*) "tl_residence time is less than internal time step: ", &
                                  tl_residence_time(kk,tl_count_flag(jpart))
                      stop
                   end if
                   exit
                end if
              end do
                     
!                elseif(tl_count_flag(jpart) == 0) then !DEBUGGING SEGMENTATION FAULT PROBLEM - 19/02/2025
!                       write(*,*) "I AM IN PARTICLE MOVE FROM TL/STRATOSPHERE TO TROPOSPHERE"
!                       write(*,*) "jpart", jpart, "tl_count_flag(jpart)", tl_count_flag(jpart)
!                       write(*,*) "strato_count_flag(jpart)", strato_count_flag(jpart)
!                       write(*,*) "zold", zold, "tri_old", tri, "cold_tropo_old", cold_tropo_old
!                       write(*,*) "ztra1(jpart)", ztra1(jpart), "tri", tri
!                       stop
           end if
            
        !Particle move from Strato to Tropo
        elseif(tropo_method == 1 .AND. (zold .gt. cold_tropo_old) .OR. &!METHOD 1
               tropo_method == 2 .AND. (pp_old .lt. 7000)) then        !METHOD 2
               !Since the particle is exiting the stratosphere, it is possible to compute the residence time
               kk = 0
               !Loop over all the previous elements of particle_id(l)
               do
                 kk = kk + 1
                 !Check that the particle was already present in tl_particle_id, if not stop.
                 if(kk .gt. size(strato_particle_id)) then
                    write(*,*) "kk index is greater than strato_particle_id"
                    write(*,*) "kk: ", kk
                    write(*,*) "strato_particle_id(jpart): ", strato_particle_id(jpart)
                    write(*,*) "I AM IN PARTICLE MOVING FROM STRATO TO TROPOSPHERE"
                    stop
                 end if
                 !Finding array index location of the particle
                 if(strato_particle_id(kk) == npoint_id(jpart)) then
                    strato_residence_time(kk,strato_count_flag(jpart)) = jtime + lsynctime/2.-&
                                                                         strato_entering_time(kk,strato_count_flag(jpart))
                    !Check. If residence time less than internal time step, stop.
                    if(strato_residence_time(kk,strato_count_flag(jpart)) < 180) then
                       write(*,*) "strato residence time is less than internal time step: ", &
                                   strato_residence_time(kk,strato_count_flag(jpart))
                       stop
                    end if
                    exit
                 end if
               end do
        end if  
     end if !END PARTICLE MOVE FROM TL/STRATO TO TROPO
     
     !PARTICLE REMAINS IN THE TL UNTIL THE END OF SIMULATION 
     !if(tropo_method == 1 .AND. ((jtime == ideltas .AND. (zold .ge. tri_old .AND. zold .le. cold_tropo_old) .AND. &           !METHOD 1 
     !                            (ztra1(jpart) .ge. tri .AND. ztra1(jpart) .le. cold_tropo) .AND. &
     !                            (tl_count_flag(jpart) .ge. 1 .AND. tl_count_flag(jpart) /= strato_count_flag(jpart)))) .OR. &  !METHOD 2
     !   tropo_method == 2 .AND. ((jtime == ideltas .AND. (pp_old .ge. 7000 .AND. pp_old .le. 15000) .AND. &
     !                            (pp_new .ge. 7000 .AND. pp_new .le. 15000) .AND. &
     !                            (tl_count_flag(jpart) .ge. 1 .AND. tl_count_flag(jpart) /= strato_count_flag(jpart))))) then 

     !if(tropo_method == 1 .AND. ((newtime == ideltas .AND. (zold .ge. tri_old .AND. zold .le. cold_tropo_old) .AND. &           !METHOD 1 
     !                            (ztra1(jpart) .ge. tri .AND. ztra1(jpart) .le. cold_tropo))) .OR. & !.AND. &
        !                         (tl_count_flag(jpart) .ge. 1 .AND. tl_count_flag(jpart) /= strato_count_flag(jpart)))) .OR. &  !METHOD 2
     if(tropo_method == 1 .AND. ((newtime == ideltas .AND. &           !METHOD 1 
                                 (ztra1(jpart) .ge. tri .AND. ztra1(jpart) .le. cold_tropo))) .OR. & !.AND. & 
       tropo_method == 2 .AND. ((newtime == ideltas .AND. 
                                (pp_new .le. 15000 .AND. pp_new .ge. 7000)))) then 

         kk = 0
         !Loop over all the previous elements of particle_id(l)
         do
           kk = kk + 1 
           !Check that the particle was already present in tl_particle_id, if not stop.
           if(kk .gt. size(tl_particle_id)) then
              write(*,*) "I AM IN PARTICLE REMAINS IN THE TTL UNTIL END SIMULATION"
              write(*,*) "kk index is greater than tl_particle_id"
              write(*,*) "kk: ", kk
              write(*,*) "tl_particle_id(jpart): ", tl_particle_id(jpart)
              write(*,*) "itime", jtime, "idelta", ideltas
              write(*,*) "jpart", jpart
              write(*,*) "tl_count_flag(jpart)", tl_count_flag(jpart)
              write(*,*) "strato_count_flag(jpart)", strato_count_flag(jpart)
              stop
           end if
           !Finding array index location of the particle
           if(tl_particle_id(kk) == npoint_id(jpart)) then
              tl_residence_time(kk,tl_count_flag(jpart)) = ideltas-&
                                                           tl_entering_time(kk,tl_count_flag(jpart))
              exit
           end if
         end do
    
    end if !END PARTICLE REMAINS IN THE TL UNTIL THE END OF SIMULATION
     
     !PARTICLE REMAINS IN THE STRATO UNTIL THE END OF SIMULATION
   !  if(tropo_method == 1 .AND. ((jtime == ideltas .AND. (zold .gt. cold_tropo_old) .AND. (ztra1(jpart) .gt. cold_tropo))) .OR. & !METHOD 1
   !     tropo_method == 2 .AND. ((jtime == ideltas .AND. (pp_old .lt. 7000) .AND. (pp_new .lt. 7000)))) then                      !METHOD 2
    !  if(tropo_method == 1 .AND. ((newtime == ideltas .AND. (zold .gt. cold_tropo_old) .AND. (ztra1(jpart) .gt. cold_tropo))) .OR. & !METHOD 1
      if(tropo_method == 1 .AND. ((newtime == ideltas .AND. (ztra1(jpart) .gt. cold_tropo))) .OR. & !METHOD 1
         tropo_method == 2 .AND. ((newtime == ideltas .AND. (pp_new .lt. 7000))) then                      !METHOD 2
         write(*,*) "I AM IN UP INTRUSION IDENTIFIER - END OF SIMULATION"
         kk = 0
         !Loop over all the previous elements of particle_id(l)
         do
          kk = kk + 1
          !Check that the particle was already present in tl_particle_id, if not stop.
          if(kk .gt. size(strato_particle_id)) then
             write(*,*) "kk index is greater than strato_particle_id"
             write(*,*) "kk: ", kk
             write(*,*) "strato_particle_id(jpart): ", strato_particle_id(jpart)
             write(*,*) "I AM IN PARTICLE REMAINS IN THE STRATO UNTIL END SIMULATION"
             stop
          end if
          !Finding array index location of the particle
          if(strato_particle_id(kk) == npoint_id(jpart)) then
             strato_residence_time(kk,strato_count_flag(jpart)) = ideltas -& 
                                                                  strato_entering_time(kk,strato_count_flag(jpart))
             exit
           end if
         end do
     end if !END PARTICLE REMAINS IN THE STRATO UNTIL THE END OF THE SIMULATION
     
  !PARTICLE REMAINS IN THE TTL/STRATOSPHERE UNTIL IT IS TERMINATED
  elseif(itra1(jpart) < 0) then

          !Old position
     call boundaries_computation(jtime,jpart,xold,yold,zold,&
                                 tri_old,cold_tropo_old,&
                                 pp_old,hmixi_old)

         !TL  
          if((tropo_method == 1 .AND. (zold .ge. tri_old .AND. zold .le. cold_tropo_old) .AND. &
              tl_count_flag(jpart) .ge. 1) .OR. &
             (tropo_method == 2 .AND. ((pp_old .ge. 7000 .AND. pp_old .le. 15000) &                                                   !METHOD 2
                               .AND.  (tl_count_flag(jpart) .ge. 1)))) then
              kk = 0
              !Loop over all the previous elements of particle_id(l)
              do      
                kk = kk + 1
                !Check that the particle was already present in tl_particle_id, if not stop.
                if(kk .gt. size(tl_particle_id)) then
                  write(*,*) "I AM IN PARTICLE REMAINS IN THE TTL UNTIL THEY ARE TERMINATED"
                  write(*,*) "kk index is greater than tl_particle_id"
                  write(*,*) "kk: ", kk
                  write(*,*) "tl_particle_id(kk): "     , tl_particle_id(kk)
                  write(*,*) "tl_count_flag(jpart)"        , tl_count_flag(jpart)
                  write(*,*) "strato_count_flag(jpart)", strato_count_flag(jpart)
                  write(*,*) "zold", zold, "cold_tropo_old", cold_tropo_old
                  write(*,*) "jpart", jpart
                  write(*,*) "--------------------------------------------------------------"
                  stop
                end if

                !Finding array index location of the particle
                if(tl_particle_id(kk) == npoint_id(jpart)) then
                   tl_residence_time(kk,tl_count_flag(jpart)) = jtime + lsynctime/2.-&
                                                                tl_entering_time(kk,tl_count_flag(jpart))

  !                  if(tl_count_flag(jpart) > 30) then
  !                     write(*,*) "tl_count_flag", tl_count_flag(jpart)
   !                 end if

                   exit
                end if
                 
              end do

          end if

          !STRATO
          if((tropo_method == 1 .AND. (zold .gt. cold_tropo_old) .AND. strato_count_flag(jpart) .ge. 1) .OR. &
             (tropo_method == 2 .AND. ((pp_old .lt. 7000) .AND. strato_count_flag(jpart) .ge. 1)))then             !METHOD 2
              kk = 0
              !Loop over all the previous elements of particle_id(l)
              do
                kk = kk + 1
                if(kk .gt. size(strato_particle_id)) then
                  write(*,*) "I AM IN PARTICLE REMAINS IN THE TTL UNTIL THEY ARE TERMINATED"
                  write(*,*) "kk index is greater than strato_particle_id"
                  write(*,*) "kk: ", kk, "size strato_particle_id", size(strato_particle_id)
                  write(*,*) "strato_particle_id(jpart): ", strato_particle_id(kk)
                  write(*,*) "tl_count_flag(jpart)", tl_count_flag(jpart)
                  write(*,*) "strato_count_flag(jpart)", strato_count_flag(jpart)
                  write(*,*) "zold", zold, "cold_tropo_old", cold_tropo_old,"ztra1",ztra1(jpart),"cold_tropo",cold_tropo
                  write(*,*) "itra1(jj)", itra1(jpart)
                  write(*,*) "jpart", jpart
                  write(*,*) "--------------------------------------------------------------"
                  stop

                end if
                if(strato_particle_id(kk) == npoint_id(jpart)) then
                   strato_residence_time(kk,strato_count_flag(jpart)) = jtime + lsynctime/2.-& 
                                                                  strato_entering_time(kk,strato_count_flag(jpart))
                   if(strato_count_flag(jpart) > 30) then
                      write(*,*) "strato_count_flag", strato_count_flag(jpart)
                   end if
                   exit
                end if
              end do
          end if       
  end if !END PARTICLE REMAINS IN THE TL/STRATO UNTIL IT IS TERMINATIED 

end subroutine up_intrusion_identifier



