
subroutine thermodyn_computation(jtime,jpart,xold,yold,zold,&
                                  par_index,crossing_index,&
                                  intrusion_type)
  !                                    i     i    i    i    i           i
  !*****************************************************************************
  !                                                                            *
  !     Computation of additional thermodynamical variables in order           *
  !     to better depicts the trajectories of those particles entering         *
  !     the tropopause layer/stratosphere (in case of sources in the PBL)      *
  !     or the PBL (in case of sources located in the stratosphere.            *
  !     The additional computed variables are:                                 *
  !     - Potential Vorticity                                                  *
  !     - Bulk Richardson Number                                               *
  !     - Brunt-Vaisala Frequency                                              *
  !     - Turbulence Index (Ellrod et al., 1992)                               *
  !     - Pressure                                                             *
  !     - Vertical velocity                                                    *
  !                                                                            *
  !     Author: Massimo Martina, PhD Student (MFF UK, Prague)                  *
  !             massimo.martina@matfyz.cuni.cz                                 *
  !                                                                            *
  !     07 January 2025                                                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! jtime                 Actual temporal position of calculation              *
  ! jpart                 Index of the particle considered                     *
  ! xold,yold,zold        "Memorized" old positions of the particle            *
  ! par_index             row index representing one parcel                    *
  ! cross_index           coloumn index representing the number of time the    *
  !                       parcel entered the TL/Strato                         *
  ! intrusion_type        1 = intrusion into the TL; 2 = into the Strato       *
  !*****************************************************************************

  use par_mod
  use com_mod
  use ai1_intrus_mod
  implicit none
  !VARIABLES
  !******************************
  !Input
  integer           :: jtime, jpart
  real(kind=dp)     :: xold, yold
  real              :: zold
 !Internal variables 
  integer           :: indexh, k, m, ngrid !MM_07-03-2025, addition of "ngrid". It was forgotten but the code works anyway, why?
  integer           :: ix,jy,ixp,jyp       !MM_07-03-2025 addition. They were forgotten but the code worked anyway, why?
  integer           :: indz, indzp         !MM_10-03-2025 addition. They were forgotten but the code worked anyway, why?
  real              :: xtn, ytn
  real              :: hmixi, hm(2)   
  real              :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  real              :: dz1,dz2,dz                               
  integer           :: il, ind
  real              :: pv1(2), pvprof(2)!,pvtraj 
  real(kind=dp)     :: jul_traj,jul_term_rel,jul_term
  integer           :: jjjjmmdd_traj, ihmmss_traj          
  character         :: adate_traj*8,atime_traj*6         
  character(len=20) :: time_traj
  real              :: ww1(2), wwprof(2)!, wtraj
  real              :: pp_par(2), pp_parprof(2),pp_new 
  real              :: ttv(2),ttvprof(2),uup(2),uupprof(2),vvp(2),vvpprof(2)
  real              :: thetaprof(2),thetam!,bvfsqtraj,ritraj,pptraj          
  real              :: dz1_old,dz2_old,dz_old                               
  real              :: tkep(2),tkepprof(2)!,tketraj
  real              :: uutraj,vvtraj,ulg,vlg,xpar,ypar,xlg,ylg,vws,dst,dsh,def,cvg!,ei
  !logical           :: fract_hours                                      
  integer           :: par_index,crossing_index
  integer           :: intrusion_type
  !END VARIABLES
  !*****************************
  
!  write(*,*) "I am in thermo_dyn_computation sabroutine"  
!  write(*,*) "xold", xold,"yold",yold,"zold",zold
!  stop

  !TRAJECTORY CHARACTERIZATION - COMPUTATION OF THERMODYNAMIC VARIABLES
  !**********************************************************************
  !Variables for temporal interpolation        
  dt1=float((jtime)-memtime(1)) 
  dt2=float(memtime(2)-(jtime))
  dtt=1./(dt1+dt2)

  !Determine if nested grid are used
  !************************************************
  ! If partoutput_use_nested=0, set ngrid=0, and use the outermost grid
  ! for calculating topo, pv, qv, ... at the particle position
  ! Otherwise, determine the nest we are in

  ngrid=0
  if(partoutput_use_nested .gt. 0) then
     do k=numbnests,1,-1
        if ((xold.gt.xln(k)).and. &
            (xold.lt.xrn(k)).and. &
            (yold.gt.yln(k)).and. &
            (yold.lt.yrn(k))) then
             ngrid=k
             goto 26
        endif
     enddo
26   continue
  endif

  if(ngrid .le. 0) then
     ix=int(xold)
     jy=int(yold)
     ddy=yold-float(jy)
     ddx=xold-float(ix)
  else
     xtn=(xold-xln(ngrid))*xresoln(ngrid)
     ytn=(yold-yln(ngrid))*yresoln(ngrid)
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

  !Computation of the topography below the particle position
  if (ngrid .le. 0) then
      if(intrusion_type == 1) then
         topo_cross_tl(par_index,crossing_index) =p1*oro(ix ,jy) &
                                                 + p2*oro(ixp,jy) &
                                                 + p3*oro(ix ,jyp) &
                                                 + p4*oro(ixp,jyp)
      else if(intrusion_type == 2) then
         topo_cross_st(par_index,crossing_index) =p1*oro(ix ,jy) &
                                                 + p2*oro(ixp,jy) &
                                                 + p3*oro(ix ,jyp) &
                                                 + p4*oro(ixp,jyp)
      else
         write(*,*) "Intrusion_type must be equal to 1 'TL' or 2 'Strato'"
         write(*,*) "intrusion_type = ", intrusion_type
         stop
      end if
  else
      !       topo=p1*oron(ix ,jy ,ngrid) &
      !          + p2*oron(ixp,jy ,ngrid) &
      !          + p3*oron(ix ,jyp,ngrid) &
      !          + p4*oron(ixp,jyp,ngrid)
  endif

  !Computation of the pressure at the particle's position;
  !It is used to determine the potential temperature and
  !to obtain the pressure reference for the Bulk Richardson Number 
  !and Brunt-Vaisala Frequency
  !*******************************************************
  do il=2,nz
     if(height(il).gt.ztra1(jpart)) then
        indz=il-1
        indzp=il
        goto 65
     endif
  enddo
65 continue

   dz1=ztra1(jpart)-height(indz)
   dz2=height(indzp)-ztra1(jpart)
   dz=1./(dz1+dz2)
   do ind=indz,indzp
      do m=1,2
         indexh=memind(m)
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
      pp_parprof(ind-indz+1)=(pp_par(1)*dt2+pp_par(2)*dt1)*dtt
   end do
!   pp_new = (pp_parprof(1)*dz2+pp_parprof(2)*dz1)*dz

   !Computation of thermodynamic variables
 
!   if(mod(jtime+lsynctime,900) == 0) then !temporal if, the value are computed every 15 minutes
      
      jul_traj=bdate+real(jtime,kind=dp)/86400._dp        !this is the particle time
      call caldate(jul_traj,jjjjmmdd_traj,ihmmss_traj)
      write(adate_traj,'(i8.8)') jjjjmmdd_traj
      write(atime_traj,'(i6.6)') ihmmss_traj
      write(time_traj,*) adate_traj(1:4),"/",adate_traj(5:6),"/",adate_traj(7:8)," ",&
                         atime_traj(1:2), ":", atime_traj(3:4),":",atime_traj(5:6)

      !The values are separeted in different file in order to reduce their size.
      !The parameter "fract_hours" controls the splitting of the results in different
      !files, each of them rapresenting a "x" number of hours.
      !if(fract_hours) then
      !close(unitquantities)
      !open(unitquantities,file=path(1)(1:length(1))//'trajcharact_startdate_'//adate_traj// &
      !     atime_traj//'.txt',form='formatted')
      !end if
          
         !WRITE(*,*) "I AM ABOUT TO COMPUTE THE DYNAMICA/THERMODYNAMIC QUANTITIES"
         !Vertical variables needed for vertical interpolation
         do il=2,nz
            if(height(il).gt.zold) then
               indz=il-1
               indzp=il
               goto 56
            endif
         enddo
56       continue

         dz1=zold-height(indz)
         dz2=height(indzp)-zold
         dz=1./(dz1+dz2)
             
         !----BILINEAR INTERPOLATION----
         do ind=indz,indzp
            do m=1,2
               indexh=memind(m)
               if(ngrid .le. 0) then
                  !Potential Vorticity
                  pv1(m)=p1*pv(ix,jy,ind,indexh) &
                         +p2*pv(ixp,jy,ind,indexh) &
                         +p3*pv(ix,jyp,ind,indexh) &
                         +p4*pv(ixp,jyp,ind,indexh)
                  !Vertical Velocity
                  ww1(m)=p1*ww(ix,jy,ind,indexh) &
                        +p2*ww(ixp,jy,ind,indexh) &
                        +p3*ww(ix,jyp,ind,indexh) &
                        +p4*ww(ixp,jyp,ind,indexh)
                  !Bulk Richardson Number variables - added on 28/08/2024
                  ttv(m)=p1*tt(ix,jy,ind,indexh)*(1+0.608*qv(ix,jy,ind,indexh))   &
                        +p2*tt(ixp,jy,ind,indexh)*(1+0.608*qv(ixp,jy,ind,indexh)) &
                        +p3*tt(ix,jyp,ind,indexh)*(1+0.608*qv(ix,jyp,ind,indexh)) &
                        +p4*tt(ixp,jyp,ind,indexh)*(1+0.608*qv(ixp,jyp,ind,indexh))

                  uup(m)=p1*uu(ix,jy,ind,indexh)  &
                        +p2*uu(ixp,jy,ind,indexh) &
                        +p3*uu(ix,jyp,ind,indexh) &
                        +p4*uu(ixp,jyp,ind,indexh)

                  vvp(m)=p1*vv(ix,jy,ind,indexh)  &
                        +p2*vv(ixp,jy,ind,indexh) &
                        +p3*vv(ix,jyp,ind,indexh) &
                        +p4*vv(ixp,jyp,ind,indexh)
                  !Turbulent Kinetic Energy
                  tkep(m)=p1*tke(ix,jy,ind,indexh)  &
                         +p2*tke(ixp,jy,ind,indexh) &
                         +p3*tke(ix,jyp,ind,indexh) &
                         +p4*tke(ixp,jyp,ind,indexh)
                       
               else
                  pv1(m)=p1*pvn(ix,jy,ind,indexh,ngrid) &
                        +p2*pvn(ixp,jy,ind,indexh,ngrid) &
                        +p3*pvn(ix,jyp,ind,indexh,ngrid) &
                        +p4*pvn(ixp,jyp,ind,indexh,ngrid)

                  ww1(m)=p1*wwn(ix,jy,ind,indexh,ngrid) &
                        +p2*wwn(ixp,jy,ind,indexh,ngrid) &
                        +p3*wwn(ix,jyp,ind,indexh,ngrid) &
                        +p4*wwn(ixp,jyp,ind,indexh,ngrid)

                  !Bulk Richardson Number variables - added on 28/08/2024
                  ttv(m)=p1*ttn(ix,jy,ind,indexh,ngrid)*(1+0.608*qvn(ix,jy,ind,indexh,ngrid))   &
                        +p2*ttn(ixp,jy,ind,indexh,ngrid)*(1+0.608*qvn(ixp,jy,ind,indexh,ngrid)) &
                        +p3*ttn(ix,jyp,ind,indexh,ngrid)*(1+0.608*qvn(ix,jyp,ind,indexh,ngrid)) &
                        +p4*ttn(ixp,jyp,ind,indexh,ngrid)*(1+0.608*qvn(ixp,jyp,ind,indexh,ngrid))

                  uup(m)=p1*uun(ix,jy,ind,indexh,ngrid)  &
                        +p2*uun(ixp,jy,ind,indexh,ngrid) &
                        +p3*uun(ix,jyp,ind,indexh,ngrid) &
                        +p4*uun(ixp,jyp,ind,indexh,ngrid)

                  vvp(m)=p1*vvn(ix,jy,ind,indexh,ngrid)  &
                        +p2*vvn(ixp,jy,ind,indexh,ngrid) &
                        +p3*vvn(ix,jyp,ind,indexh,ngrid) &
                        +p4*vvn(ixp,jyp,ind,indexh,ngrid)
                  !Turbulent Kinetic Energy  
                  tkep(m)=p1*tken(ix,jy,ind,indexh,ngrid)  &
                         +p2*tken(ixp,jy,ind,indexh,ngrid) &
                         +p3*tken(ix,jyp,ind,indexh,ngrid) &
                         +p4*tken(ixp,jyp,ind,indexh,ngrid)

               end if 
            end do
            !----TEMPORAL INTERPOLATION----
            !Potential Vorticity
            pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
            !Vertical Velocity
            wwprof(ind-indz+1)=(ww1(1)*dt2+ww1(2)*dt1)*dtt
            !Bulk Richardson number variables
            ttvprof(ind-indz+1)=(ttv(1)*dt2+ttv(2)*dt1)*dtt !added on 28/08/2024
            thetaprof(ind-indz+1)=ttvprof(ind-indz+1)*(1000./pp_parprof(ind-indz+1))**(r_air/cpa) !Potential temp added on 28/08/2024
            uupprof(ind-indz+1)=(uup(1)*dt2+uup(2)*dt1)*dtt !added on 28/08/2024
            vvpprof(ind-indz+1)=(vvp(1)*dt2+vvp(2)*dt1)*dtt !added on 28/08/2024
            !Turbulent Kinetic Energy
            tkepprof(ind-indz+1)=(tkep(1)*dt2+tkep(2)*dt1)*dtt !added on 28/08/2024
         end do
         !----VERTICAL INTERPOLATION----
         if(intrusion_type == 1) then
            !Potential Vorticity
            pvtraj_tl(par_index,crossing_index) = (pvprof(1)*dz2+pvprof(2)*dz1)*dz
            !if(pvtraj < -200) then
                !write(*,*) "pvtraj is less than -200, I am going to stop the simulation"
                !stop
            !end if
            !Vertical velocity
            wtraj_tl(par_index,crossing_index)  = (wwprof(1)*dz2+wwprof(2)*dz1)*dz
            !The following line, until pptraj, added on 28/08/2024
            !Turbulent Kinetic Energy
            tketraj_tl(par_index,crossing_index) = (tkepprof(1)*dz2+tkepprof(2)*dz1)*dz
            !Bulk Richardson Number
            ritraj_tl(par_index,crossing_index) = &
                   ga/thetaprof(1) * (thetaprof(2)-thetaprof(1))*(height(indzp)-height(indz))/ &
                  ((uupprof(2)-uupprof(1))**2 + (vvpprof(2)-vvpprof(1))**2) !Bulk Richardson Number computed considering two vertical levels
                                                                            !and associated to the particle within this layer
            !Potential Temperature
            thetam = 0.5*(thetaprof(1)+thetaprof(2))
            !Brunt Vaisala Frequency
            bvfsqtraj_tl(par_index,crossing_index) = ga/thetam*(thetaprof(2)-thetaprof(1))/(height(indzp)-height(indz))
            !Pressure reference fro Bulk Richardson Number and Brunt-Vaisala Frequency
            pptraj_tl(par_index,crossing_index)    = 0.5*(pp_parprof(1)+pp_parprof(2)) 
            if(pptraj_tl(par_index,crossing_index) > 100000) then
               write(*,*) "pp_parprof(1)", pp_parprof(1)
               write(*,*) "pp_parprof(2)", pp_parprof(2)
               write(*,*) "pp_new", pp_new
               write(*,*) "pp_par(1)", pp_par(1)
               write(*,*) "pp_par(2)", pp_par(2)
               stop
            end if
            !Ellrod Index
            uutraj = (uupprof(2)*dz1 + uupprof(1)*dz2)*dz
            vvtraj = (vvpprof(2)*dz1 + vvpprof(1)*dz2)*dz
            ulg = (uu(ix,jy,indz,indexh)*dz2 + uu(ix,jy,indzp,indexh)*dz2)*dz !vertical interpolation to particle height of u at the left grid point
            vlg = (vv(ix,jy,indz,indexh)*dz2 + vv(ix,jy,indzp,indexh)*dz2)*dz
            xpar=xmet0+xold*dx !IT IS NOT IN LON DEGREES --> CHECK ELLROND COMPUTATION
            ypar=ymet0+yold*dy !IT IS NOT IN LAT DEGREES --> CHECK ELLROND COMPUTATION
            xlg =xmet0+ix*dx
            ylg =ymet0+jy*dy
            vws = ((uupprof(2)-uupprof(1))**2 + (vvpprof(2)-vvpprof(1))**2)/(height(indzp)-height(indz)) !vertical wind shear
            dst = (uutraj-ulg)/(xpar-xlg) - (vvtraj-vlg)/(ypar-ylg)
            dsh = (vvtraj-vlg)/(xpar-xlg) + (uutraj-ulg)/(ypar-ylg)
            cvg = -((uutraj-ulg)/(xpar-xlg) + (vvtraj-vlg)/(ypar-ylg))
            def = sqrt(dst**2 + dsh**2)
            eitraj_tl(par_index,crossing_index) = vws*(def+cvg) !Ellrod Index
     
         else if(intrusion_type == 2) then
            !Potential Vorticity
            pvtraj_st(par_index,crossing_index) = (pvprof(1)*dz2+pvprof(2)*dz1)*dz
            !if(pvtraj < -200) then
                !write(*,*) "pvtraj is less than -200, I am going to stop the simulation"
                !stop
            !end if
            !Vertical velocity
            wtraj_st(par_index,crossing_index)  = (wwprof(1)*dz2+wwprof(2)*dz1)*dz
            !The following line, until pptraj, added on 28/08/2024
            !Turbulent Kinetic Energy
            tketraj_st(par_index,crossing_index) = (tkepprof(1)*dz2+tkepprof(2)*dz1)*dz
            !Bulk Richardson Number
            ritraj_st(par_index,crossing_index) = &
                   ga/thetaprof(1) * (thetaprof(2)-thetaprof(1))*(height(indzp)-height(indz))/ &
                  ((uupprof(2)-uupprof(1))**2 + (vvpprof(2)-vvpprof(1))**2) !Bulk Richardson Number computed considering two vertical levels
                                                                            !and associated to the particle within this layer
            !Potential Temperature
            thetam = 0.5*(thetaprof(1)+thetaprof(2))
            !Brunt Vaisala Frequency
            bvfsqtraj_st(par_index,crossing_index) = ga/thetam*(thetaprof(2)-thetaprof(1))/(height(indzp)-height(indz))
            !Pressure reference fro Bulk Richardson Number and Brunt-Vaisala Frequency
            pptraj_st(par_index,crossing_index)    = 0.5*(pp_parprof(1)+pp_parprof(2))
            if(pptraj_st(par_index,crossing_index) > 100000) then
               write(*,*) "pp_parprof(1)", pp_parprof(1)
               write(*,*) "pp_parprof(2)", pp_parprof(2)
               write(*,*) "pp_new", pp_new
               write(*,*) "pp_par(1)", pp_par(1)
               write(*,*) "pp_par(2)", pp_par(2)
               stop
            end if
           
           !Ellrod Index
            uutraj = (uupprof(2)*dz1 + uupprof(1)*dz2)*dz
            vvtraj = (vvpprof(2)*dz1 + vvpprof(1)*dz2)*dz
            ulg = (uu(ix,jy,indz,indexh)*dz2 + uu(ix,jy,indzp,indexh)*dz2)*dz !vertical interpolation to particle height of u at the left grid point
            vlg = (vv(ix,jy,indz,indexh)*dz2 + vv(ix,jy,indzp,indexh)*dz2)*dz
            xpar=xmet0+xold*dx !IT IS NOT IN LON DEGREES --> CHECK ELLROND COMPUTATION
            ypar=ymet0+yold*dy !IT IS NOT IN LAT DEGREES --> CHECK ELLROND COMPUTATION
            xlg =xmet0+ix*dx
            ylg =ymet0+jy*dy
            vws = ((uupprof(2)-uupprof(1))**2 + (vvpprof(2)-vvpprof(1))**2)/(height(indzp)-height(indz)) !vertical wind shear
            dst = (uutraj-ulg)/(xpar-xlg) - (vvtraj-vlg)/(ypar-ylg)
            dsh = (vvtraj-vlg)/(xpar-xlg) + (uutraj-ulg)/(ypar-ylg)
            cvg = -((uutraj-ulg)/(xpar-xlg) + (vvtraj-vlg)/(ypar-ylg))
            def = sqrt(dst**2 + dsh**2)
            eitraj_st(par_index,crossing_index) = vws*(def+cvg) !Ellrod Index

         else
            write(*,*) "Intrusion_type must be equal to 1 'TL' or 2 'Strato'"
            write(*,*) "intrusion_type = ", intrusion_type
            stop
         end if

end subroutine thermodyn_computation



