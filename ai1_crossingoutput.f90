
subroutine crossingoutput(filename,adate,atime)
  !*****************************************************************************
  !                                                                            *
  !     This subroutine is entirely based on the "fluxoutput.f90" subroutine,  *
  !     so all the credits goes to the Flexpart developpers.                   *
  !     This subroutine outputs the fluxes computed by "calccrossing.f90".     *
  !     The output file is written in either sparse matrix or grid dump        *
  !     format, whichever is more efficient.                                   *
  !                                                                            *
  !     Author: Massimo Martina, PhD Student (MFF UK, Prague)                  *
  !             massimo.martina@matfyz.cuni.cz                                 *
  !                                                                            *
  !     16 January 2025                                                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! (Input)                                                                    *
  ! filename        name of the output file                                    *
  ! adate           reference date for the output                              *
  ! atime           reference time for the output                              *
  !                                                                            *
  ! (Internal)                                                                 *
  ! ncellse         number of cells with non-zero values for eastward fluxes   *
  ! sparsee         .true. if in sparse matrix format, else .false.            *
  !                                                                            *
  !*****************************************************************************

  use outg_mod
  use par_mod
  use com_mod
  use ai1_intrus_mod

  implicit none

  real(kind=dp) :: jul
  character (25) :: filename,kspec
  integer :: ix,jy,kz,k,nage,jjjjmmdd,ihmmss,kp,i
  integer :: ncell_firstcrossing, ncell1 
  logical :: sparse_firstcrossing,sparse1
  character :: adate*8,atime*6

  real :: mean_entering_time, mean_transition_time
  real :: mean_pv,mean_w,mean_tke,mean_pp,mean_ri,mean_ei,mean_bvfsq,tot_mass
  real(kind=dp) :: jul_cross              
  integer :: jjjjmmdd_cross, ihmmss_cross
  character :: adate_cross*8,atime_cross*6
  character(len=22) :: a_cross   


  !**************************************************************
  ! Check, whether output of full grid or sparse matrix format is
  ! more efficient in terms of storage space. 
  !**************************************************************

  ncell_firstcrossing = 0
  ncell1 = 0

  do kp=1,maxpointspec_act
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
            if (flux_crossing(1,ix,jy,kp).gt.0) then
!                write(*,*) "I AM IN CROSSINGOUTPUT, CHECKING SPARSE METHOD"
!                write(*,*) "flux_crossing(1,...", flux_crossing(1,ix,jy,kp)
                ncell_firstcrossing= ncell_firstcrossing+1
            end if
            if (flux_crossing(2,ix,jy,kp) .gt. 0) then
!                 write(*,*) "I AM IN CROSSINGOUTPUT, CHECKING SPARSE METHOD"
!                 write(*,*) "flux_crossing(2,...", flux_crossing(2,ix,jy,kp)
                ncell1= ncell1+1
            end if
        end do
      end do
  end do

  ! Output in sparse matrix format more efficient, if less than
  ! 2/5 of all cells contains concentrations>0
  !************************************************************
  if (4*ncell_firstcrossing.lt.numxgrid*numygrid*numzgrid) then
      sparse_firstcrossing=.true.
      write(*,*) "The sparse method for the first crossing fluxes is acrivated"
  else
      sparse_firstcrossing=.false.
      write(*,*) "The sparse method for the first crossing fluxes is NOT acrivated"
      endif
  if (4*ncell1.lt.numxgrid*numygrid*numzgrid) then
      sparse1=.true.
      write(*,*) "The sparse method for the +1 fluxes is acrivated"
  else
      sparse1=.false.
      write(*,*) "The sparse method for the +1 fluxes is NOT acrivated"
  endif

  !ACTIVATE ONLY THE SPARSE METHOD - MM - 11/10/2024
  !write(*,*) "THE SPARSE METHOD IS FORCED TO BE ACTIVE"
  !sparse_firstcrossing = .true.
  !sparse1 = .true.

  !Open file - Intrusions data
  open(unit_res,file=path(1)(1:length(1))//trim(filename)//"_"//adate//atime,form='unformatted')  

  !Open fiile - Mass data
  open(unit_mass,file=path(1)(1:length(1))//trim(filename)//"_mass_"//adate//atime,form='unformatted')
  do k =1,nspec
  do kp=1,maxpointspec_act
      if (sparse_firstcrossing) then
        write(unit_mass) 1
!        if(k == 1) then
           write(unit_res) 1
!        end if
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux_mass(1,ix,jy,k,kp).gt.0) then

                 write(unit_mass) ix+jy*numxgrid, &
                                     flux_mass(1,ix,jy,k,kp)/area(ix,jy)/outstep, area(ix,jy), outstep

              end if

                 if (flux_crossing(1,ix,jy,kp).gt.0) then
                     mean_entering_time = sum_entering_time(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     mean_transition_time = sum_transition_time(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
               
                     mean_pv = sum_pvtraj(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     mean_w = sum_wtraj(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     mean_tke = sum_tketraj(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     mean_pp = sum_pptraj(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     mean_ri = sum_ritraj(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     mean_ei = sum_eitraj(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     mean_bvfsq = sum_bvfsqtraj(1,ix,jy,kp)/flux_crossing(1,ix,jy,kp)
                     tot_mass = flux_mass(1,ix,jy,k,kp) 

                     jul_cross=bdate+real(mean_entering_time,kind=dp)/86400._dp        !this is the particle crossing time

                     call caldate(jul_cross,jjjjmmdd_cross,ihmmss_cross)
                     write(adate_cross,'(i8.8)') jjjjmmdd_cross
                     write(atime_cross,'(i6.6)') ihmmss_cross
                     write(a_cross,*) adate_cross(1:4),"/",adate_cross(5:6),"/",adate_cross(7:8)," ",&
                                      atime_cross(1:2), ":", atime_cross(3:4),":",atime_cross(5:6)

                     write(unit_res) ix+jy*numxgrid, &
                                     flux_crossing(1,ix,jy,kp),&
                                     mean_entering_time,&
                                     a_cross,&
                                     mean_transition_time,&
                                     mean_pv,mean_w,mean_tke,&
                                     mean_pp,mean_ri,mean_ei,&
                                     mean_bvfsq,tot_mass
                   end if
            end do
          end do
        write(unit_mass) -999,999.
        write(unit_res) -999,999.
      else
        write(unit_mass) 2
        write(unit_res) 2
          do ix=0,numxgrid-1
!            write(unit_resclass)  (flux_crossing(1,ix,jy,kp,nage),jy=0,numygrid-1)
!            write(unit_resclass) (sum_entering_time(1,ix,jy,kp,nage),jy=0,numygrid-1)
!            write(unit_resclass) (sum_transition_time(1,ix,jy,kp,nage),jy=0,numygrid-1)
          end do
      endif

      if (sparse1) then
        write(unit_mass) 1
        write(unit_res) 1
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux_mass(2,ix,jy,k,kp).gt.0) then

                 write(unit_mass) ix+jy*numxgrid,&
                                      flux_mass(2,ix,jy,k,kp)/area(ix,jy)/outstep
              end if
             
                 if (flux_crossing(2,ix,jy,kp).gt.0) then

                     mean_entering_time = sum_entering_time(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     mean_pv = sum_pvtraj(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     mean_w = sum_wtraj(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     mean_tke = sum_tketraj(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     mean_pp = sum_pptraj(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     mean_ri = sum_ritraj(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     mean_ei = sum_eitraj(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     mean_bvfsq = sum_bvfsqtraj(2,ix,jy,kp)/flux_crossing(2,ix,jy,kp)
                     tot_mass = flux_mass(2,ix,jy,k,kp)

                     jul_cross=bdate+real(mean_entering_time,kind=dp)/86400._dp        !this is the particle crossing time

                     call caldate(jul_cross,jjjjmmdd_cross,ihmmss_cross)
                     write(adate_cross,'(i8.8)') jjjjmmdd_cross
                     write(atime_cross,'(i6.6)') ihmmss_cross
                     write(a_cross,*) adate_cross(1:4),"/",adate_cross(5:6),"/",adate_cross(7:8)," ",&
                                      atime_cross(1:2), ":", atime_cross(3:4),":",atime_cross(5:6)

                     write(unit_res) ix+jy*numxgrid,&
                                     flux_crossing(2,ix,jy,kp),&
                                      mean_entering_time,&
                                      a_cross,&
                                      sum_transition_time(2,ix,jy,kp),&
                                      mean_pv,mean_w,mean_tke,&
                                      mean_pp,mean_ri,mean_ei,&
                                      mean_bvfsq,tot_mass

                  end if
            end do
          end do
        write(unit_mass) -999,999.
        write(unit_res) -999,999.
      else
        write(unit_mass) 2
        write(unit_res) 2
        do ix=0,numxgrid-1
      !      write(unit_resclass) (flux_crossing(2,ix,jy,k,kp),jy=0,numygrid-1) !BUG: flux_crossing(1,....) instead of flux_crossing(2,...)
      !      write(unit_resclass) (sum_entering_time(2,ix,jy,k,kp),jy=0,numygrid-1) !BUG: sum_entering_time(1,....) instead of sum_entering_time(2,...)
      !      write(unit_resclass) (sum_transition_time(2,ix,jy,k,kp),jy=0,numygrid-1) !BUG: sum_transition_time(1,....) instead of sum_transition_time(2,...)
        end do
      endif
  end do !end do over release point kp
  end do !end loop over species 
  close(unit_res)
  close(unit_mass)
  
  ! Reinitialization of grid
  !*************************
  do k=1,nspec
  do kp=1,maxpointspec_act
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
              do i=1,2
                flux_mass(i,ix,jy,k,kp)    = 0
                if(k == 1) then
                   flux_crossing(i,ix,jy,kp)=0
                   sum_entering_time(i,ix,jy,kp)=0
                   sum_transition_time(i,ix,jy,kp)=0
                   sum_pvtraj(i,ix,jy,kp)=0
                   sum_wtraj(i,ix,jy,kp)=0
                   sum_tketraj(i,ix,jy,kp)=0
                   sum_pptraj(i,ix,jy,kp)=0
                   sum_ritraj(i,ix,jy,kp)=0
                   sum_eitraj(i,ix,jy,kp)=0
                   sum_bvfsqtraj(i,ix,jy,kp)=0
                end if
              end do
      end do
    end do
  end do
  end do

end subroutine crossingoutput
