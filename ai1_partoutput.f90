!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!*                                                                     *
!* This file is part of FLEXPART WRF                                   *
!*                                                                     *
!* FLEXPART is free software: you can redistribute it and/or modify    *
!* it under the terms of the GNU General Public License as published by*
!* the Free Software Foundation, either version 3 of the License, or   *
!* (at your option) any later version.                                 *
!*                                                                     *
!* FLEXPART is distributed in the hope that it will be useful,         *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of      *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
!* GNU General Public License for more details.                        *
!*                                                                     *
!* You should have received a copy of the GNU General Public License   *
!* along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!***********************************************************************

      subroutine partoutput(itime,zold)
!                             i
!*******************************************************************************
!                                                                              *
!     Note:  This is the FLEXPART_WRF version of subroutine partoutput.        *
!                                                                              *
!     Dump all particle positions                                              *
!                                                                              *
!     Author: A. Stohl                                                         *
!     12 March 1999                                                            *
!                                                                              *
!     Dec 2005, J. Fast - Output files can be either binary or ascii.          *
!                   Topo,pv,qv,... at particle positions are calculated        *
!                   using nested fields when (partoutput_use_nested .gt. 0)    *
!                   Particle xy coords can be either lat-lon or grid-meters.   *
!                   Changed names of "*lon0*" & "*lat0*" variables             *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************
!*******************************************************************************
!             MODIFICATIONS TO USE THE AirIntrusion Package                    *
!                                                                              *
! Author: Massimo Martina (PhD student - MFF KFA Charles University)           *
!         massimo.martina@matfyz.cuni.cz                                       *
!                                                                              *
! All the changes are introduced by "!AirIntrusion Modification - date"        *
!                                                                              *
! To know more about this package, please visit its repository at:             *
! https://github.com/MassiMart/FLEXPART-WRF_AirIntrusions.git                  *
!                                                                              *
! 05 August 2025                                                               *
!*******************************************************************************

  use par_mod
  use com_mod
!AirIntrusions Modification - 06-12-2024
  use ai1_intrus_mod
!End AirIntrusion Modification
  implicit none
!      include 'includepar'
!      include 'includecom'

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  integer :: ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  integer :: numpart_out
  real :: xlon,ylat,xtmp,ytmp
  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: topo,hm(2),hmixi,pv1(2),pvprof(2),pvi,qv1(2),qvprof(2),qvi
  real :: tt1(2),ttprof(2),tti,rho1(2),rhoprof(2),rhoi
!AirIntrusions Modification - 06-12-2024
  real :: ww1(2), wwprof(2), wpar
  real :: pp_par(2), pp_parprof(2),pp_old,pp_new                      
  real :: ttv(2),ttvprof(2),uup(2),uupprof(2),vvp(2),vvpprof(2)       
  real :: thetaprof(2),thetam,bvfsqpar,ripar,pppar
  real :: tkep(2),tkepprof(2),tkepar                 
  real :: cold(2), cold_tropo  
  real :: capep(2),capepprof(2),capepar
!End AirIntrusion Modification
  real :: tr(2),tri
  character :: adate*8,atime*6
  integer :: k, ngrid
  real :: xtn, ytn
    
!AirIntrusions Modification - xx-05/06-2024
!******************************
  real, dimension (maxpart) :: zold    !MM_28-05-2024
  integer :: jj                        !MM_08-06-2024
  real(kind=dp) :: jul_rel             !MM_12-06-2024
  integer :: jjjjmmdd_rel, ihmmss_rel  !MM_12-06-2024
  character :: adate_rel*8,atime_rel*6 !MM_12-06-2024
! Determine current calendar date, needed for the file name
!**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp    ! this is the current day

      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss

! Some variables needed for temporal interpolation
!*************************************************

  dt1=float(itime-memtime(1))
  dt2=float(memtime(2)-itime)
  dtt=1./(dt1+dt2)

! Open output file and write the output
!**************************************

  if (ipout.eq.1) then
      if (iouttype.eq.0 .or. iouttype .eq.2) &
          open(unitpartout,file=path(1)(1:length(1))//'partposit_'//adate// &
          atime,form='unformatted')
      if (iouttype.eq.1) &
          open(unitpartout,file=path(1)(1:length(1))//'partposit_'//adate// &
          atime,form='formatted')
  else
        if (iouttype.eq.0 .or. iouttype .eq.2) &
            open(unitpartout,file=path(1)(1:length(1))//'partposit_end', &
                 form='unformatted')
        if (iouttype.eq.1) &
            open(unitpartout,file=path(1)(1:length(1))//'partposit_end', &
                 form='formatted')
  endif

! Write current time to file
!***************************

  numpart_out = 0
  do i=1,numpart
     if (itra1(i).eq.itime) numpart_out = numpart_out + 1
  enddo

  if (iouttype.eq.0 .or. iouttype .eq.2) write(unitpartout)   itime, &
      numpart_out, outgrid_option
  if (iouttype.eq.1) write(unitpartout,*) itime, &
      numpart_out, outgrid_option 

  do i=1,numpart

     ! Take only valid particles
     !**************************
     if (itra1(i).eq.itime) then
         xlon=xmet0+xtra1(i)*dx
         ylat=ymet0+ytra1(i)*dy
         !*********************************************************************************
         ! Interpolate several variables (PV, specific humidity, etc.) to particle position
         !*********************************************************************************
         ! If partoutput_use_nested=0, set ngrid=0, and use the outermost grid
         ! for calculating topo, pv, qv, ... at the particle position
         ! Otherwise, determine the nest we are in
         ngrid=0
         if (partoutput_use_nested .gt. 0) then
             do k=numbnests,1,-1
                if ((xtra1(i).gt.xln(k)).and. &
                    (xtra1(i).lt.xrn(k)).and. &
                    (ytra1(i).gt.yln(k)).and. &
                    (ytra1(i).lt.yrn(k))) then
                    ngrid=k
                    goto 26
                endif
             enddo
26           continue
         endif

         if (ngrid .le. 0) then
             ix=int(xtra1(i))
             jy=int(ytra1(i))
             ddy=ytra1(i)-float(jy)
             ddx=xtra1(i)-float(ix)
         else
             xtn=(xtra1(i)-xln(ngrid))*xresoln(ngrid)
             ytn=(ytra1(i)-yln(ngrid))*yresoln(ngrid)
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

         ! Topography
         !***********
         if (ngrid .le. 0) then
             topo=p1*oro(ix ,jy) &
                 + p2*oro(ixp,jy) &
                 + p3*oro(ix ,jyp) &
                 + p4*oro(ixp,jyp)
         else
             topo=p1*oron(ix ,jy ,ngrid) &
                 + p2*oron(ixp,jy ,ngrid) &
                 + p3*oron(ix ,jyp,ngrid) &
                 + p4*oron(ixp,jyp,ngrid)
         endif

         ! Potential vorticity, specific humidity, temperature, and density
         !*****************************************************************

         do il=2,nz
            if (height(il).gt.ztra1(i)) then
                indz=il-1
                indzp=il
                goto 56
            endif
         enddo
56       continue

         dz1=ztra1(i)-height(indz)
         dz2=height(indzp)-ztra1(i)
         dz=1./(dz1+dz2)

         do ind=indz,indzp
            do m=1,2
               indexh=memind(m)
    
               if (ngrid .le. 0) then
                  ! Potential vorticity
                  pv1(m)=p1*pv(ix ,jy ,ind,indexh) &
                        +p2*pv(ixp,jy ,ind,indexh) &
                        +p3*pv(ix ,jyp,ind,indexh) &
                        +p4*pv(ixp,jyp,ind,indexh)
                  ! Specific humidity
                  qv1(m)=p1*qv(ix ,jy ,ind,indexh) &
                        +p2*qv(ixp,jy ,ind,indexh) &
                        +p3*qv(ix ,jyp,ind,indexh) &
                        +p4*qv(ixp,jyp,ind,indexh)
                  ! Temperature
                  tt1(m)=p1*tt(ix ,jy ,ind,indexh) &
                        +p2*tt(ixp,jy ,ind,indexh) &
                        +p3*tt(ix ,jyp,ind,indexh) &
                        +p4*tt(ixp,jyp,ind,indexh)
                  ! Density
                  rho1(m)=p1*rho(ix ,jy ,ind,indexh) &
                         +p2*rho(ixp,jy ,ind,indexh) &
                         +p3*rho(ix ,jyp,ind,indexh) &
                         +p4*rho(ixp,jyp,ind,indexh)
!AirIntrusions Modification - 08-11-2024
                  ! Pressure
                  pp_par(m) = p1*pph(ix,jy,ind,indexh) &
                            + p2*pph(ixp,jy,ind,indexh) &
                            + p3*pph(ix,jyp,ind,indexh) &
                            + p4*pph(ixp,jyp,ind,indexh)
                  ! Vertical velocity
                  ww1(m) = p1*ww(ix,jy,ind,indexh) &
                         + p2*ww(ixp,jy,ind,indexh) &
                         + p3*ww(ix,jyp,ind,indexh) &
                         + p4*ww(ixp,jyp,ind,indexh)
                  ! Virtual temperature for Bulk Richardson number 
                  ttv(m) = p1*tt(ix,jy,ind,indexh)*(1+0.608*qv(ix,jy,ind,indexh))&  
                         + p2*tt(ixp,jy,ind,indexh)*(1+0.608*qv(ixp,jy,ind,indexh)) &
                         + p3*tt(ix,jyp,ind,indexh)*(1+0.608*qv(ix,jyp,ind,indexh)) &
                         + p4*tt(ixp,jyp,ind,indexh)*(1+0.608*qv(ixp,jyp,ind,indexh))
                  ! U-wind velocity for Bulk Richardson number
                  uup(m) = p1*uu(ix,jy,ind,indexh)  &
                         + p2*uu(ixp,jy,ind,indexh) &
                         + p3*uu(ix,jyp,ind,indexh) &
                         + p4*uu(ixp,jyp,ind,indexh)
                  ! V-wind velocity for Bulk Richardson number
                  vvp(m) = p1*vv(ix,jy,ind,indexh)  &
                         + p2*vv(ixp,jy,ind,indexh) &
                         + p3*vv(ix,jyp,ind,indexh) &
                         + p4*vv(ixp,jyp,ind,indexh)
                  ! TKE
                  tkep(m)= p1*tke(ix,jy,ind,indexh)  &
                         + p2*tke(ixp,jy,ind,indexh) &
                         + p3*tke(ix,jyp,ind,indexh) &
                         + p4*tke(ixp,jyp,ind,indexh)

                  capep(m)=p1*cape(ix,jy,ind,indexh)  &
                          + p2*cape(ixp,jy,ind,indexh) &
                          + p3*cape(ix,jyp,ind,indexh) &
                          + p4*cape(ixp,jyp,ind,indexh)
!End AirIntrusions Modification
               else
                  pv1(m)=p1*pvn(ix ,jy ,ind,indexh,ngrid) &
                        +p2*pvn(ixp,jy ,ind,indexh,ngrid) &
                        +p3*pvn(ix ,jyp,ind,indexh,ngrid) &
                        +p4*pvn(ixp,jyp,ind,indexh,ngrid)
                  qv1(m)=p1*qvn(ix ,jy ,ind,indexh,ngrid) &
                        +p2*qvn(ixp,jy ,ind,indexh,ngrid) &
                        +p3*qvn(ix ,jyp,ind,indexh,ngrid) &
                        +p4*qvn(ixp,jyp,ind,indexh,ngrid)
                  tt1(m)=p1*ttn(ix ,jy ,ind,indexh,ngrid) &
                        +p2*ttn(ixp,jy ,ind,indexh,ngrid) &
                        +p3*ttn(ix ,jyp,ind,indexh,ngrid) &
                        +p4*ttn(ixp,jyp,ind,indexh,ngrid)
                  rho1(m)=p1*rhon(ix ,jy ,ind,indexh,ngrid) &
                         +p2*rhon(ixp,jy ,ind,indexh,ngrid) &
                         +p3*rhon(ix ,jyp,ind,indexh,ngrid) &
                         +p4*rhon(ixp,jyp,ind,indexh,ngrid)
!AirIntrusions Modifications - 08-11-2024
                  pp_par(m) = p1*pphn(ix,jy,ind,indexh,ngrid) &
                            + p2*pphn(ixp,jy,ind,indexh,ngrid) &
                            + p3*pphn(ix,jyp,ind,indexh,ngrid) &
                            + p4*pphn(ixp,jyp,ind,indexh,ngrid)
                  ww1(m) = p1*wwn(ix,jy,ind,indexh,ngrid) &
                         + p2*wwn(ixp,jy,ind,indexh,ngrid) &
                         + p3*wwn(ix,jyp,ind,indexh,ngrid) &
                         + p4*wwn(ixp,jyp,ind,indexh,ngrid)
                  !Bulk Richardson Number - added on 28/08/2024
                  ttv(m) = p1*ttn(ix,jy,ind,indexh,ngrid)*(1+0.608*qvn(ix,jy,ind,indexh,ngrid))   &
                         + p2*ttn(ixp,jy,ind,indexh,ngrid)*(1+0.608*qvn(ixp,jy,ind,indexh,ngrid)) &
                         + p3*ttn(ix,jyp,ind,indexh,ngrid)*(1+0.608*qvn(ix,jyp,ind,indexh,ngrid)) &
                         + p4*ttn(ixp,jyp,ind,indexh,ngrid)*(1+0.608*qvn(ixp,jyp,ind,indexh,ngrid))
                  uup(m) = p1*uun(ix,jy,ind,indexh,ngrid)  &
                         + p2*uun(ixp,jy,ind,indexh,ngrid) &
                         + p3*uun(ix,jyp,ind,indexh,ngrid) &
                         + p4*uun(ixp,jyp,ind,indexh,ngrid)
                  vvp(m) = p1*vvn(ix,jy,ind,indexh,ngrid)  &
                         + p2*vvn(ixp,jy,ind,indexh,ngrid) &
                         + p3*vvn(ix,jyp,ind,indexh,ngrid) &
                         + p4*vvn(ixp,jyp,ind,indexh,ngrid)
                  tkep(m)= p1*tken(ix,jy,ind,indexh,ngrid)  &
                         + p2*tken(ixp,jy,ind,indexh,ngrid) &
                         + p3*tken(ix,jyp,ind,indexh,ngrid) &
                         + p4*tken(ixp,jyp,ind,indexh,ngrid)

                  capep(m)=p1*capen(ix,jy,ind,indexh,ngrid)  &
                          + p2*capen(ixp,jy,ind,indexh,ngrid) &
                          + p3*capen(ix,jyp,ind,indexh,ngrid) &
                          + p4*capen(ixp,jyp,ind,indexh,ngrid)
!End AirIntrusions Modification
               endif
            enddo !end loop over two time step
            pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
            qvprof(ind-indz+1)=(qv1(1)*dt2+qv1(2)*dt1)*dtt
            ttprof(ind-indz+1)=(tt1(1)*dt2+tt1(2)*dt1)*dtt
            rhoprof(ind-indz+1)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
!AirIntrusions Modifications - 08-11-2024    
            pp_parprof(ind-indz+1)=(pp_par(1)*dt2+pp_par(2)*dt1)*dtt
            wwprof(ind-indz+1)=(ww1(1)*dt2+ww1(2)*dt1)*dtt
            ttvprof(ind-indz+1)=(ttv(1)*dt2+ttv(2)*dt1)*dtt 
            thetaprof(ind-indz+1)=ttvprof(ind-indz+1)*(1000./pp_parprof(ind-indz+1))**(r_air/cpa) !Potential temp 
            uupprof(ind-indz+1)=(uup(1)*dt2+uup(2)*dt1)*dtt
            vvpprof(ind-indz+1)=(vvp(1)*dt2+vvp(2)*dt1)*dtt 
            tkepprof(ind-indz+1)=(tkep(1)*dt2+tkep(2)*dt1)*dtt
            capepprof(ind-indz+1)=(capep(1)*dt2+capep(2)*dt1)*dtt
!End AirIntrusions Modification
         enddo !end do over vertical axis
         pvi=(dz1*pvprof(2)+dz2*pvprof(1))*dz
         qvi=(dz1*qvprof(2)+dz2*qvprof(1))*dz
         tti=(dz1*ttprof(2)+dz2*ttprof(1))*dz
         rhoi=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz
!AirIntrusions Modifications - 08-11-2024      
         pp_new = (pp_parprof(1)*dz2+pp_parprof(2)*dz1)*dz
         wpar  = (wwprof(1)*dz2+wwprof(2)*dz1)*dz
         ripar = ga/thetaprof(1) * (thetaprof(2)-thetaprof(1))*(height(indzp)-height(indz))/ &
                 ((uupprof(2)-uupprof(1))**2 + (vvpprof(2)-vvpprof(1))**2)
         thetam = 0.5*(thetaprof(1)+thetaprof(2))
         bvfsqpar = ga/thetam*(thetaprof(2)-thetaprof(1))/(height(indzp)-height(indz)) !Brunt-Vaisala frequency
         pppar    = 0.5*(pp_parprof(1)+pp_parprof(2)) !pressure reference for Ri and Brunt-Vaisala Frequency
         tkepar   = (dz1*tkepprof(2) + dz2*tkepprof(1))*dz !MM_06032025
         capepar  = (dz1*capepprof(2) + dz2*capepprof(1))*dz
!End AirIntrusions Modification

         ! Tropopause and PBL height
         !**************************

         do m=1,2
            indexh=memind(m)
            if (ngrid .le. 0) then
                ! Tropopause
                tr(m)=p1*tropopause(ix ,jy ,1,indexh) &
                     + p2*tropopause(ixp,jy ,1,indexh) &
                     + p3*tropopause(ix ,jyp,1,indexh) &
                     + p4*tropopause(ixp,jyp,1,indexh)
                ! PBL height
                hm(m)=p1*hmix(ix ,jy ,1,indexh) &
                     + p2*hmix(ixp,jy ,1,indexh) &
                     + p3*hmix(ix ,jyp,1,indexh) &
                     + p4*hmix(ixp,jyp,1,indexh)
!AirIntrusions Modifications -             
                cold(m) = p1*cold_tropopause(ix ,jy ,1,indexh) &
                        + p2*cold_tropopause(ixp,jy ,1,indexh) &
                        + p3*cold_tropopause(ix ,jyp,1,indexh) &
                        + p4*cold_tropopause(ixp,jyp,1,indexh)
!End AirIntrusions Modification
            else
                tr(m)=p1*tropopausen(ix ,jy ,1,indexh,ngrid) &
                     + p2*tropopausen(ixp,jy ,1,indexh,ngrid) &
                     + p3*tropopausen(ix ,jyp,1,indexh,ngrid) &
                     + p4*tropopausen(ixp,jyp,1,indexh,ngrid)
                hm(m)=p1*hmixn(ix ,jy ,1,indexh,ngrid) &
                     + p2*hmixn(ixp,jy ,1,indexh,ngrid) &
                     + p3*hmixn(ix ,jyp,1,indexh,ngrid) &
                     + p4*hmixn(ixp,jyp,1,indexh,ngrid)
!AirIntrusions Modifications -  
                !cold tropopause 
!End AirIntrusions Modifications

            endif

         enddo

         hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
         tri=(tr(1)*dt2+tr(2)*dt1)*dtt
!AirIntrusions Modifications -
         cold_tropo = (cold(1)*dt2+cold(2)*dt1)*dtt
!End AirIntrusions Modifications

         ! Write the output
         !*****************
         if (outgrid_option .eq. 1) then
             xtmp = xlon
             ytmp = ylat 
             ! WRITE(*,*) "I AM IN partoutput: xlon, ylat", xlon, ylat !MM_05-03-2024
             call xymeter_to_ll_wrf( xtmp, ytmp, xlon, ylat )

         endif
      
         if(iouttype.eq.0 .or. iouttype .eq.2)  &   
            write(unitpartout) npoint(i),npoint_id(i),xlon,ylat,ztra1(i), &
                               itramem(i),&!,aw(i),
                               topo,pvi,qvi,rhoi,pppar,wpar,ripar,bvfsqpar,tkepar,capepar,hmixi,tri,tti, &
                               cold_tropo,&
                               (xmass1(i,j),j=1,nspec)
         if(iouttype.eq.1)    &                                                           
            write(unitpartout,101) npoint(i),npoint_id(i),itramem(i),xlon,ylat, & !MM_VERSION
                                   ztra1(i),&!,aw(i),
                                   topo,pvi,qvi,rhoi,pppar,wpar,ripar,bvfsqpar,tkepar,capepar,hmixi,tri,tti, &
                                   cold_tropo,&
                                   (xmass1(i,j),j=1,nspec)
             !ORIGINAL VERSION
             ! write(unitpartout,101) npoint(i),itramem(i),xlon,ylat, &
             ! ztra1(i),topo,pvi,qvi,rhoi,hmixi,tri,tti, &
             ! (xmass1(i,j),j=1,nspec)

     endif !end if take only valid particles
  enddo !end loop over the particles
 
  if(iouttype.eq.0 .or. iouttype .eq.2)  &
     write(unitpartout) -99999,-9999.9,-9999.9,-9999.9,-99999, &
                        -9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9, &
                       (-9999.9,j=1,nspec)
  if(iouttype.eq.1)  &
     write(unitpartout,101) -99999,-99999,-9999.9,-9999.9,-9999.9, &
                            -9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9, &
                           (-9999.9,j=1,nspec)

101   format( 2i10, 1p, 21e14.6 )
901   format(A8,1x,A6, i6,1x,i6,1x,i6,1x,i6)

  close(unitpartout)

end subroutine partoutput

