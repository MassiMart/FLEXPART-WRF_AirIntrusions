
subroutine boundaries_computation(itime,jpart,xpar,ypar,zpar,&
                                  tripar,coldpar,&
                                  pppar,hpblpar)
  !                                    
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the Planetary Boundary Layer height, the tropopause     *
  !     boundaries and the pressure at the particle position for two           *
  !     contiguous time steps.                                                 *
  !     These vaules are essential for the subroutine intrusion_identifier.f90 *
  !     whose aim is to identify those particles (1) entering the tropopause   *
  !     layer if the sources are located in the PBL or (2) entering the PBL    *
  !     if the sources are located in the stratosphere.                        *
  !     The tropopause boundaries can be computed considering two different    *
  !     definition:                                                            *
  !     (A) the tropopause is defined computing its bottom with the Lapse Rate *
  !         method and its top with the Cold Point method.                     *
  !     (B) the tropopause is defined considering the climatological definition*
  !         (using pressure surfaces) proposed in the article by               *
  !         Fueglistaler et al. (2009). In this option there is not a real     *
  !         computation of the tropopause boundaries since these are defined   *
  !         as: top = 70 hPa and bottom = 150 hPa. In this section the pressure*
  !         at the particle position is computed in order to compare it to     *
  !         the tropopause boundaries and hence identify its crossing          *
  !                                                                            *
  !     Author: Massimo Martina, PhD Student (MFF UK, Prague)                  *
  !             massimo.martina@matfyz.cuni.cz                                 *
  !                                                                            *
  !     07 January 2025                                                        *
  !     Last modifications:                                                    *
  !     19 June 2025 Reorganizing the code in a clearer and more flexible way, *
  !                  avoiding the repetition of the code for computing         *
  !                  the values for jtime and jtime+lsynctime.                 *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! itime            (IN) Actual temporal position of calculation              *
  ! jpart            (IN) Index of the particle considered                     *
  ! xpar,ypar,zpar   (IN) Positions of the particle                            *
  ! tropo_method     (IN) Index for choosing the tropopause boundaries         *
  !                       computation method                                   *
  ! tripar           (OUT) Lapse rate tropopause height                        *
  ! coldpar          (OUT) Cold Point Tropopause height                        *
  ! hpblpar          (OUT) PBL height
  !*****************************************************************************


  use par_mod
  use com_mod
  use ai1_intrus_mod

  implicit none
  !VARIABLES
  !******************************
  !Input
  integer       :: itime, jpart
  real(kind=dp) :: xpar, ypar
  real          :: zpar
  !Internal variables 
  integer       :: indexh, k, m,ngrid !MM_07-03-2025 addition of "ngrid". It was forgotten but the code works anyway, why?
  integer       :: ix,jy,ixp,jyp      !MM_07-03-2025 addition. They were forgotten but the code worked anyway, why?
  integer       :: indz, indzp        !MM_10-03-2025 addition. They were forgotten but the code worked anyway, why?
  real          :: xtn, ytn
  real          :: hpblpar,hm(2)                         
  real          :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2                  
  real          :: dz1,dz2,dz                                                 
  integer       :: il, ind                                                 
  real          :: pp_par(2), pp_parprof(2),pppar
  real          :: tr(2),coldpar,cold(2),tripar
  !END VARIABLES
  !*****************************
 
  !Some variables needed for temporal interpolation
  !*************************************************
  dt1=float(itime-memtime(1))
  dt2=float(memtime(2)-itime)
  dtt=1./(dt1+dt2)

  !Determine if nested grid are used
  !*********************************
  ! If partoutput_use_nested=0, set ngrid=0, and use the outermost grid
  ! for calculating pbl height and tropopause boundaries at the particle position
  ! Otherwise, determine the nest we are in
   
  ngrid=0
  if(partoutput_use_nested .gt. 0) then
     do k=numbnests,1,-1
        if((xpar.gt.xln(k)).and. &
           (xpar.lt.xrn(k)).and. &
           (ypar.gt.yln(k)).and. &
           (ypar.lt.yrn(k))) then
           
            ngrid=k
            goto 36
        endif
     enddo
36   continue
  endif

  !Computation of variables needed for Bilinear Interpolation
  if(ngrid .le. 0) then
     ix=int(xpar)
     jy=int(ypar)
     ddy=ypar-float(jy)
     ddx=xpar-float(ix)
  else
     xtn=(xpar-xln(ngrid))*xresoln(ngrid)
     ytn=(ypar-yln(ngrid))*yresoln(ngrid)
     ix=int(xtn)
     jy=int(ytn)
     ddy=ytn-float(jy)
     ddx=xtn-float(ix)
  endif
  ixp=ix+1
  jyp=jy+1
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy

  !PBL HEIGHT POSITION - BILINEAR INTERPALATION
  !***********************************************************
  do m=1,2
     indexh=memind(m)
     if(ngrid .le. 0) then
          !PBL height 
        hm(m)=p1*hmix(ix ,jy ,1,indexh) &
             +p2*hmix(ixp,jy ,1,indexh) &
             +p3*hmix(ix ,jyp,1,indexh) &
             +p4*hmix(ixp,jyp,1,indexh)

     else
            
        hm(m)=p1*hmixn(ix ,jy ,1,indexh,ngrid) &
             +p2*hmixn(ixp,jy ,1,indexh,ngrid) &
             +p3*hmixn(ix ,jyp,1,indexh,ngrid) &
             +p4*hmixn(ixp,jyp,1,indexh,ngrid)
     endif

  enddo

  !PBL HEIGHT POSITION - TEMPORAL INTERPALATION
  !***********************************************************
  hpblpar=(hm(1)*dt2+hm(2)*dt1)*dtt

  !TROPOPAUSE DEFINED BY LAPSE RATE AND COLD POINT METHODS
  !*******************************************************
  if(tropo_method == 1) then
     !Tropopause position - Bilinear Interpolation
     do m=1,2
        indexh=memind(m)

        if(ngrid .le. 0) then
           !Lapse Rate Tropopause
           tr(m)=p1*tropopause(ix ,jy ,1,indexh) &
                +p2*tropopause(ixp,jy ,1,indexh) &
                +p3*tropopause(ix ,jyp,1,indexh) &
                +p4*tropopause(ixp,jyp,1,indexh)
           !Cold Tropopause
           cold(m)=p1*cold_tropopause(ix ,jy ,1,indexh) &
                  +p2*cold_tropopause(ixp,jy ,1,indexh) &
                  +p3*cold_tropopause(ix ,jyp,1,indexh) &
                  +p4*cold_tropopause(ixp,jyp,1,indexh)

        else
           tr(m)=p1*tropopausen(ix ,jy ,1,indexh,ngrid) &
                +p2*tropopausen(ixp,jy ,1,indexh,ngrid) &
                +p3*tropopausen(ix ,jyp,1,indexh,ngrid) &
                +p4*tropopausen(ixp,jyp,1,indexh,ngrid)

           cold(m)=p1*cold_tropopausen(ix ,jy ,1,indexh,ngrid) &
                  +p2*cold_tropopausen(ixp,jy ,1,indexh,ngrid) &
                  +p3*cold_tropopausen(ix ,jyp,1,indexh,ngrid) &
                  +p4*cold_tropopausen(ixp,jyp,1,indexh,ngrid) 
        endif

     enddo
     !Tropopause height position - Tempoaral Interpolation
     tripar=(tr(1)*dt2+tr(2)*dt1)*dtt
     coldpar=(cold(1)*dt2+cold(2)*dt1)*dtt 

  end if

  !TROPOPAUSE DEFINED USING THE CLIAMATOLOGICAL DEFINITION - FUEGLISTALER ET AL.(2009)
  !***********************************************************************************
  if(tropo_method == 0) then
     !Position - Pressure computation
     do il=2,nz
        if(height(il).gt.zpar) then
           indz=il-1
           indzp=il
           goto 75
        endif
     enddo
75   continue

     dz1=zpar-height(indz)
     dz2=height(indzp)-zpar
     dz=1./(dz1+dz2)
          
     do ind=indz,indzp
        do m=1,2
           indexh=memind(m)
           !Bilinear Interpolation
           if(ngrid .le. 0) then
              pp_par(m)=p1*pph(ix,jy,ind,indexh) &
                       +p2*pph(ixp,jy,ind,indexh) &
                       +p3*pph(ix,jyp,ind,indexh) &
                       +p4*pph(ixp,jyp,ind,indexh)

           else
              pp_par(m)=p1*pphn(ix,jy,ind,indexh,ngrid) &
                       +p2*pphn(ixp,jy,ind,indexh,ngrid) &
                       +p3*pphn(ix,jyp,ind,indexh,ngrid) &
                       +p4*pphn(ixp,jyp,ind,indexh,ngrid)
           end if
        end do
        !Temporal interpolation
        pp_parprof(ind-indz+1)=(pp_par(1)*dt2+pp_par(2)*dt1)*dtt
     end do
     !Vertical Interpolation
     pppar = (pp_parprof(1)*dz2+pp_parprof(2)*dz1)*dz
     
     !Check on the pp value     
     if(pppar < 0) then
        write(*,*) "xpar", xpar, "ypar", ypar, "zpar", zpar, "itime", itime,"jpart", jpart
        write(*,*) "p1", p1, "p2", p2, "p3", p3,"p4", p4
        write(*,*) "dt1", dt1, "dt2", dt2, "dtt", dtt     !MM_24/09/2024 DEBUG
        write(*,*) "pp_par(1)", pp_par(1), "pp_par(2)", pp_par(2)                 !MM_24/09/2024 DEBUG
        write(*,*) "pp_parprof(1)", pp_parprof(1), "pp_parprof(2)", pp_parprof(2) !MM_24/09/2024 DEBUG
        write(*,*) "dz1", dz1, "dz2", dz2, "dz", dz       !MM_24/09/2024 DEBUG
        write(*,*) "pppar", pppar
        stop                                                                      !MM_24/09/2024 DEBUG
     end if
  end if
end subroutine boundaries_computation



