
subroutine calccrossing(layer_identifier,jpart,xcross,ycross,entering_time,transition_time,&
                        count_flag,pv_par,w_par,tke_par,pp_par,ri_par,ei_par,bvfsq_par)
  !*****************************************************************************
  !                                                                            *
  !     This subroutine is entirely based on the "calcfluxes.f90" subroutine,  *
  !     so all the credits goes to the Flexpart developpers.                   *
  !     This subroutine computes the "fluxes" of particles toward the TL/      *
  !     Stratosphere (when the sources are in the PBL) or of particles toward  *
  !     the PBL (when the sources are in the Stratosphere).                    *
  !     The fluxes are splitted in two parts:                                  *
  !     (A) flux that consider the first time the particle enter               *
  !         in the considered layer;                                           *
  !     (B) flux that takes into account all the later entrance of particles   *
  !         in the considered layer  (multiplecrossing)                        *
  !                                                                            *
  !     Author: Massimo Martina, PhD Student (MFF UK, Prague)                  *
  !             massimo.martina@matfyz.cuni.cz                                 *
  !                                                                            *
  !     16 January 2025                                                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nage                  Age class of the particle considered                 *
  ! jpart                 Index of the particle considered                     *
  ! xcross,ycross         Horizontal position of the particle when it enters   *
  !                       the considered layer                                 * 
  ! entering_time         Time at which the particle enters the considerd      *
  !                       layer                                                *
  ! transition_time       Time needed by the particle for moving from the      *
  !                       the source layer to the considered one               *
  ! count_flag            Number of time the particles enters the considered   *
  !                       layer (useful to distinguish among the two fluxes    *
  !*****************************************************************************

!  use outg_mod
  use par_mod
  use com_mod
  use ai1_intrus_mod
  implicit none

  integer :: jpart,ixave,jyave,kp
  integer :: k,ix,jy,i,x4i,y4i
  real(kind=dp) :: xcross,ycross
  real          :: xmean,ymean, ixave_1, jyave_1
  real:: entering_time,transition_time
  integer:: count_flag, layer_identifier,i_count
  real :: pv_par,w_par,tke_par,pp_par,ri_par,ei_par,bvfsq_par
  real :: xlon,ylat,xtmp,ytmp, x_lon1,y_lat1,x_lon2,y_lat2,tmpx,tmpy,tmplon,tmplat
  real :: xl2,yl2,xtmp4,ytmp4,xl4,yl4,xtmp1,ytmp1
!  write(*,*) "I am in calccrossing"
  ! Determine average positions
  !****************************

  if ((ioutputforeachrelease.eq.1).and.(mdomainfill.eq.0)) then
     kp=npoint(jpart)
  else
     kp=1
  endif

  xmean = xcross
  ymean = ycross

!  write(*,*) "xmean", xmean
!  write(*,*) "ymean", ymean

  ixave=int((xmean*dx+xoutshift)/dxout)
  jyave=int((ymean*dy+youtshift)/dyout)

  !Counting the intrusions into the TL/STRATO
  !*************************   
  if ((ixave.ge.0).and.(jyave.ge.0).and.(ixave.le.numxgrid-1).and. &
      (jyave.le.numygrid-1)) then

      if(outgrid_option == 1) then
        !USEFUL FOR FLUX COMPUTATION           
        xtmp = xmet0 + xcross*dx
        ytmp = ymet0 + ycross*dy
        
        call xymeter_to_ll_wrf( xtmp, ytmp, xlon, ylat )

       !DEBUG PROCESS-22-09-2025
        call xyindex_to_ll_wrf(0,xmean,ymean,x_lon1,y_lat1)
       
!       tmpx=out_xm0+(ixave)*dxout
!       tmpy=out_ym0+(jyave)*dyout
!       call xymeter_to_ll_wrf_out(tmpx,tmpy,tmplon,tmplat)

        xl2=outlon0+(ixave)*dxoutl !long 
        yl2=outlat0+(jyave)*dyoutl !lat  

       ! call xyindex_to_ll_wrf_out(0,ixave_1,jyave_1,x_lon2,y_lat2)

        xtmp4 = out_xm0 + ixave*dxout 
        ytmp4 = out_ym0 + jyave*dyout

        call xymeter_to_ll_wrf(xtmp4,ytmp4,xl4,yl4)
!        call xymeter_to_xyindex_wrf(xtmp4,ytmp4,x4i,y4i)

      else
        xtmp = out_xm0 + ixave*dx
        ytmp = out_ym0 + jyave*dy

        call xymeter_to_ll_wrf(xtmp,ytmp,xlon,ylat)

        xtmp1 = xmet0 + xcross*dx
        ytmp1 = ymet0 + ycross*dy

        call xymeter_to_ll_wrf( xtmp1, ytmp1, x_lon1, y_lat1)

      end if     
    
      if(layer_identifier == 0) then

          write(unit_res,*)  tl_particle_id(jpart), count_flag, tl_residence_time(jpart,count_flag),&
                             tl_transition_time(jpart),xlon,ylat,&
                             tl_entering_time(jpart,count_flag), tl_zcross(jpart,count_flag),&
                             topo_cross_tl(jpart,count_flag),tketraj_tl(jpart,count_flag),&
                             ritraj_tl(jpart,count_flag),pptraj_tl(jpart,count_flag)
      else if(layer_identifier == 1) then
          write(unit_res,*)  strato_particle_id(jpart), count_flag, strato_residence_time(jpart,count_flag),&
                             strato_transition_time(jpart),xlon,ylat,&
                             strato_entering_time(jpart,count_flag), strato_zcross(jpart,count_flag),&
                             topo_cross_st(jpart,count_flag),tketraj_st(jpart,count_flag),&
                             ritraj_st(jpart,count_flag),pptraj_st(jpart,count_flag)
      else
          write(*,*) "layer_identifier is not 0 (tropopause) or 1 (stratosphere)"
          write(*,*) "layer_identifier: ", layer_identifier 
          stop
      end if
         
      if(count_flag == 1) then !first transition
        !if(flux_flag /= 1) then
        !   write(*,*) "The flux must be from below"
        !   stop
        !end if
        do k=1,nspec !MM_11-03-2025. Addition of mass fluxes for each species
           flux_mass(1,ixave,jyave,k,kp)     = flux_mass(1,ixave,jyave,k,kp)+xmass1(jpart,k)
        end do
        !Number of intrusions into the tl/strato
        flux_crossing(1,ixave,jyave,kp) = flux_crossing(1,ixave,jyave,kp)+1

        !Total values of some thermodynamical variables
        sum_entering_time(1,ixave,jyave,kp)     = sum_entering_time(1,ixave,jyave,kp) +&
                                                       entering_time
        sum_transition_time(1,ixave,jyave,kp)   = sum_transition_time(1,ixave,jyave,kp) +&
                                                       transition_time

        sum_pvtraj(1,ixave,jyave,kp) = sum_pvtraj(1,ixave,jyave,kp) + pv_par
        sum_wtraj(1,ixave,jyave,kp) = sum_wtraj(1,ixave,jyave,kp) + w_par
        sum_tketraj(1,ixave,jyave,kp) =  sum_tketraj(1,ixave,jyave,kp) + tke_par
        sum_pptraj(1,ixave,jyave,kp) = sum_pptraj(1,ixave,jyave,kp) + pp_par
        sum_ritraj(1,ixave,jyave,kp) = sum_ritraj(1,ixave,jyave,kp) + ri_par
        sum_eitraj(1,ixave,jyave,kp) = sum_eitraj(1,ixave,jyave,kp) + ei_par
        sum_bvfsqtraj(1,ixave,jyave,kp) = sum_bvfsqtraj(1,ixave,jyave,kp) + bvfsq_par

!        write(*,*) "sum_transition_time", sum_transition_time(1,ixave,jyave,kp,nage)
        if(sum_transition_time(1,ixave,jyave,kp) < 0 ) then
           write(*,*) "sum transition time is less than 0: ", sum_transition_time(1,ixave,jyave,kp)
           stop
        end if
 
!        write(*,*) "sum_pvtraj", sum_pvtraj(1,ixave,jyave,kp)
        if(sum_pvtraj(1,ixave,jyave,kp) > 0 .and. sum_pvtraj(1,ixave,jyave,kp)<1e-34) then
           write(*,*) "sum pvtraj is less than 0: ",sum_pvtraj(1,ixave,jyave,kp)
           stop
        end if

     
      else !flux considering particles that enter in the tl/strato more than once

!        write(*,*) "I AM IN CALCCROSSING.F90"
!        write(*,*) "count_flag", count_flag, "jpart", jpart
        do k=1,nspec
           flux_mass(2,ixave,jyave,k,kp)  = flux_mass(2,ixave,jyave,k,kp) + xmass1(jpart,k)
        end do
        !Number of intrusions into the tl/strato
        flux_crossing(2,ixave,jyave,kp) = flux_crossing(2,ixave,jyave,kp) +&
                                                 1
        !Total values of some thermodynamical variables
        sum_entering_time(2,ixave,jyave,kp)     = sum_entering_time(2,ixave,jyave,kp) +&
                                                       entering_time

        sum_transition_time(2,ixave,jyave,kp)   = 0
  
        sum_pvtraj(2,ixave,jyave,kp) = sum_pvtraj(2,ixave,jyave,kp) + pv_par
        sum_wtraj(2,ixave,jyave,kp) = sum_wtraj(2,ixave,jyave,kp) + w_par
        sum_tketraj(2,ixave,jyave,kp) =  sum_tketraj(2,ixave,jyave,kp) + tke_par
        sum_pptraj(2,ixave,jyave,kp) = sum_pptraj(2,ixave,jyave,kp) + pp_par
        sum_ritraj(2,ixave,jyave,kp) = sum_ritraj(2,ixave,jyave,kp) + ri_par
        sum_eitraj(2,ixave,jyave,kp) = sum_eitraj(2,ixave,jyave,kp) + ei_par
        sum_bvfsqtraj(2,ixave,jyave,kp) = sum_bvfsqtraj(2,ixave,jyave,kp) + bvfsq_par

      end if

      do i=1,2 
        if(sum_transition_time(i,ixave,jyave,kp) < 0  .OR. sum_entering_time(i,ixave,jyave,kp) < 0) then
           write(*,*) "sum transition time or sum entering time is less than 0"
           write(*,*) "sum transition time: ", sum_transition_time(i,ixave,jyave,kp)
           write(*,*) "sum entering time: ", sum_entering_time(i,ixave,jyave,kp)
           stop
        end if
      end do
      
  endif

end subroutine calccrossing

