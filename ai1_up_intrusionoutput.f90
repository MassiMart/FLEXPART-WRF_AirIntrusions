
subroutine up_intrusionoutput(jtime,nage)
  !*****************************************************************************
  !                                                                            *
  !     This subroutine handles the output of the intrusion in the tropopause/ *
  !     stratosphere identified by up_intrusion_identifier.f90.                *
  !     The results are written out for each hour and diveded by the           *
  !     considered residence classes                                           *
  !                                                                            *
  !     Author: Massimo Martina, PhD Student (MFF UK, Prague)                  *
  !             massimo.martina@matfyz.cuni.cz                                 *
  !                                                                            *
  !     16 January 2025                                                        *
  !     Last Modifications:                                                    *
  !     21 Juner 2025: addition of comment to improve readibility;             *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use flux_mod
  use ai1_intrus_mod

  implicit none
  !VARIABLES
  !******************************
  integer       :: tt_res,jj,k,jtime,nage,ks
  character (45):: filename_1, filename_2, filename_deb
  character :: adate*8,atime*6
  integer :: itime, time_step,jjjjmmdd,ihmmss
  real(kind=dp) :: jul
  integer       :: layer_identifier
  write(*,*) "Writing out the results from up_intrusion_identifier"
  write(*,*) "colnum_tl", colnum_tl
  write(*,*) "colnum_st", colnum_st
  write(*,*) "size(tl_particle_id)", size(tl_particle_id)
  write(*,*) "size(strato_particle_id)", size(strato_particle_id)
  if(source_type == 1) then
       !Saving ttl fluxes
       time_step = 20*lsynctime !with lsynctime = 180 s, time_step = 1h
       do itime=0,ideltas,time_step !temporal loop !MM_19-08-2025 changing itime=0,... to itime=bdate,...
                                                   !MM_02/09/2025 returning to itime = 0 as in timemanager_serial.f90
       ! Determine current calendar date, needed for the file name
       !**********************************************************

!          write(*,*) "I AM IN INTRUSIONOUTPUT - ITIME", itime

          jul=bdate+real(itime,kind=dp)/86400._dp
          call caldate(jul,jjjjmmdd,ihmmss)
          write(adate,'(i8.8)') jjjjmmdd
          write(atime,'(i6.6)') ihmmss

          do tt_res=1,nresidenceclass-1 !do cycle over the residence classes
             if(tt_res == 1) then
                filename_1 = "ttl_residence_class1"
                filename_2 = "ttl_class1_"
                filename_deb = "debug_ttl_class1_"
             elseif(tt_res == 2) then
                filename_1 = "ttl_residence_class2"
                filename_2 = "ttl_class2_"
                filename_deb = "debug_ttl_class2_"
             elseif(tt_res == 3) then
                filename_1 = "ttl_residence_class3"
                filename_2 = "ttl_class3_"
                filename_deb = "debug_ttl_class3_"
             elseif(tt_res == 4) then
                filename_1 = "ttl_residence_class4"
                filename_2 = "ttl_class4_"
                filename_deb = "debug_ttl_class4_"
             end if

             filename_1 = TRIM(filename_1)
             filename_2 = TRIM(filename_2)
             filename_deb = TRIM(filename_deb)

             open(unit_res,file=path(1)(1:length(1))//trim(filename_2)//adate//atime//".txt",form='formatted')
             !DEBUG PROCESS - 22-09-2025
 !            open(unit_deb,file=path(1)(1:length(1))//trim(filename_deb)//adate//atime//".txt",form='formatted')
             
             layer_identifier = 0 !MM_13-10-2025
             !do ks=1,nspec
             do jj=1,size(tl_particle_id)   !do cycle over the particles
                do k=1,colnum_tl      !do cycle over the multicrossing
                   !selecting particles entering the tl in the considered hour itime - itime+time_step
                   if(tl_entering_time(jj,k) .ge. itime .AND. tl_entering_time(jj,k) .lt. itime+time_step) then
                      !DEBUG - 02/09/2025
                    !  if(itime == ideltas) then
                    !     write(*,*) "tl_particle_id(jj)", tl_particle_id(jj),"k",k,"tl_entering_time(jj,k)",&
                    !                 tl_entering_time(jj,k)
                    !     write(*,*) "itime+time_step", itime+time_step
                    !  end if
                      if(tl_residence_time(jj,k) .gt. lresidence(tt_res) .AND. &
                         tl_residence_time(jj,k) .le. lresidence(tt_res+1)) then
                         if(tl_flux_flag(jj,k) /= 0) then

                            !Flux computation & writing out intrusion data in a txt file
                            call calccrossing(layer_identifier,jj,tl_xcross(jj,k),tl_ycross(jj,k),tl_entering_time(jj,k),&
                                              tl_transition_time(jj),k,pvtraj_tl(jj,k),wtraj_tl(jj,k),tketraj_tl(jj,k),&
                                              pptraj_tl(jj,k),ritraj_tl(jj,k),eitraj_tl(jj,k),bvfsqtraj_tl(jj,k))
                         end if !end if tl_flux_flag
                      end if !end if tl_residence_time
                   end if !end if tl_entering_time
                end do !end do over multicrossing
             end do !end do over particles
             !end do !end do over species
             !close the txt file
             close(unit_res)

             !FLUX OUTPUT RESULT
             call crossingoutput(filename_1,adate,atime)
          end do !end do over residence classes
       end do !end temporal loop

       !Saving strato fluxes
       do itime=0,ideltas,time_step !temporal loop
          ! Determine current calendar date, needed for the file name
          !**********************************************************

          jul=bdate+real(itime,kind=dp)/86400._dp
          call caldate(jul,jjjjmmdd,ihmmss)
          write(adate,'(i8.8)') jjjjmmdd
          write(atime,'(i6.6)') ihmmss

          do tt_res=1,nresidenceclass-1 !do cycle over the residence classes
             if(tt_res == 1) then
                filename_1 = "strato_residence_class1"
                filename_2 = "strato_class1_"
             elseif(tt_res == 2) then
                filename_1 = "strato_residence_class2"
                filename_2 = "strato_class2_"
             elseif(tt_res == 3) then
                filename_1 = "strato_residence_class3"
                filename_2 = "strato_class3_"
             elseif(tt_res == 4) then
                filename_1 = "strato_residence_class4"
                filename_2 = "strato_class4_"
             end if

             filename_1 = TRIM(filename_1)
             filename_2 = TRIM(filename_2)

             open(unit_res,file=path(1)(1:length(1))//trim(filename_2)//adate//atime//".txt",form='formatted')
             layer_identifier = 1 !MM_13-10-2025
             !do ks=1,nspec
             do jj=1,size(strato_particle_id)   !do cycle over the particles
                do k=1,colnum_st      !do cycle over the multicrossing
                   if(strato_entering_time(jj,k) .ge. itime .AND. strato_entering_time(jj,k) .lt. itime+time_step) then
                      if(strato_residence_time(jj,k) .gt. lresidence(tt_res) .AND. &
                         strato_residence_time(jj,k) .le. lresidence(tt_res+1)) then
                         if(strato_flux_flag(jj,k) /= 0) then

                            !FLUX COMPUTATION
                            !considering only first species
                            call calccrossing(layer_identifier,jj,strato_xcross(jj,k),strato_ycross(jj,k),&
                                              strato_entering_time(jj,k), strato_transition_time(jj),&
                                              k,pvtraj_st(jj,k),wtraj_st(jj,k),tketraj_st(jj,k),&
                                              pptraj_st(jj,k),ritraj_st(jj,k),eitraj_st(jj,k),bvfsqtraj_st(jj,k))
                         end if !end if strato_flux_flag
                      end if !end if strato_residence_time
                   end if !end if strato_entering time
                end do !end do multicrossing
              end do !end do over particles
              !end do !end do over species
              close(unit_res)
              call crossingoutput(filename_1,adate,atime)
          end do !end do over residence classes
       end do !end temporal loop
!!!!!!!!!!!!!!!!!!!!-------------THE FOLLOWING PART MUST BE REVISED--------------!!!!!!!!!!!!!!!!!!!
  else if(source_type == 0) then
          !Saving strato fluxes
         !  do tt_res=1,nresidenceclass-1 !do cycle over the residence classes
         !     if(tt_res == 1) then
         !        filename_1 = "intrus_residence_class1"
         !        filename_2 = "intrus_class1.txt"
         !     elseif(tt_res == 2) then
         !            filename_1 = "intrus_residence_class2"
         !            filename_2 = "intrus_class2.txt"
         !     elseif(tt_res == 3) then
         !            filename_1 = "intrus_residence_class3"
         !            filename_2 = "intrus_class3.txt"
         !     elseif(tt_res == 4) then
         !            filename_1 = "intrus_residence_class4"
         !            filename_2 = "intrus_class4.txt"
         !     end if

         !     filename_1 = TRIM(filename_1)
         !     filename_2 = TRIM(filename_2)

         !     open(unit_res,file=path(1)(1:length(1))//filename_2,form='formatted')
      !  !      do ks=1,nspec
         !     do jj=1,size(intrusion_particle_id)   !do cycle over the particles
         !        if(intrusion_type(jj) == 3) then !strong intrusion toward the PBL
         !           if(intrusion_residence_time(jj) .gt. lresidence(tt_res) .AND. &
         !              intrusion_residence_time(jj) .le. lresidence(tt_res+1)) then
         !              if(intrusion_flux_flag(jj) == 0) then
         !                 write(*,*) "problem in computing calccrossing, intrusion_flux_flag = 0 "
         !                 stop
         !              end if
         !              write(unit_res,*) intrusion_particle_id(jj), intrusion_count_flag(jj), abl_count_flag(jj),&
!                                         intrusion_exiting_time(jj), abl_entering_time(jj),&
!                                         intrusion_residence_time(jj), intrusion_transition_time(jj)
!                       call calccrossing(ks,jj,intrusion_xcross(jj), intrusion_ycross(jj),&
!                                         intrusion_exiting_time(jj), intrusion_transition_time(jj))
!                    end if
!                 end if
!              end do
!              write(*,*) "size(intrusion_particle_id)", size(intrusion_particle_id)
!              close(unit_res)
!              call crossingoutput(jtime,filename_1)
!           end do
!
!           !Saving strato fluxes
!           do tt_res=1,nresidenceclass-1 !do cycle over the residence classes
!              if(tt_res == 1) then
!                 filename_1 = "pblint_residence_class1"
!                 filename_2 = "pblint_class1.txt"
!              elseif(tt_res == 2) then
!                     filename_1 = "pblint_residence_class2"
!                     filename_2 = "pblint_class2.txt"
!              elseif(tt_res == 3) then
!                     filename_1 = "pblint_residence_class3"
!                     filename_2 = "pblint_class3.txt"
!              elseif(tt_res == 4) then
!                     filename_1 = "pblint_residence_class4"
!                     filename_2 = "pblint_class4.txt"
!              end if
!
!              filename_1 = TRIM(filename_1)
!              filename_2 = TRIM(filename_2)
!
!              open(unit_res,file=path(1)(1:length(1))//filename_2,form='formatted')
!              do jj=1,size(intrusion_particle_id)   !do cycle over the particles
!                 if(intrusion_type(jj) == 3) then !strong intrusion toward the PBL
!                    if(intrusion_residence_time(jj) .gt. lresidence(tt_res) .AND. &
!                       intrusion_residence_time(jj) .le. lresidence(tt_res+1)) then
!                       if(intrusion_flux_flag(jj) == 0) then
!                          write(*,*) "problem in computing calccrossing, intrusion_flux_flag = 0 "
!                          stop
!                       end if
!                       write(unit_res,*) intrusion_particle_id(jj), intrusion_count_flag(jj), abl_count_flag(jj), &
!                                         intrusion_exiting_time(jj), abl_entering_time(jj),&
!                                         intrusion_residence_time(jj), intrusion_transition_time(jj)
!                       call calccrossing(nage,jj,abl_xcross(jj), abl_ycross(jj),&
!                                         abl_entering_time(jj), intrusion_transition_time(jj))
!                    end if
!                 end if
!              end do
!              write(*,*) "size(intrusion_particle_id)", size(intrusion_particle_id)
!              close(unit_res)
!              call crossingoutput(jtime,filename_1)
!           end do
  end if

end subroutine up_intrusionoutput
