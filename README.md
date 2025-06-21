# FLEXPART-WRF - Description of the new "AirIntrusions" package
The AirIntrusions package is a modified version of some of the files forming the FLEXPART-WRF model (Brioude et al., 2013) developed in the article Martina et al. (yyyy), "Typhoon, orography and gravity waves - a recipe for rapid transport from the boundary layer to the strato.",xxxx.  
The original files have been modified only where strictly necessary in order to facilitate the implementation of the changes for potentially new users. Therefore, the changes have been organized in additional subroutines that must be added to the native code.  
In order to use the package, you need to download it and store its file in the same directory where the original FLEXPART-WRF code is hosted. All the files belonging to the AirIntrusions package are named as "ai1_*.f90", where "ai1" stands for "AirIntrusions package version 1".  
Then, you need to use "ai1_make.mom" to compile the modified version of the program (NOTE: you need to remove from the directory all those ".o" files connected to the original version of files that have been modified, such as "timemanager_serial.o" since it has been modified in some part in this new package.  
The aim of this package is, on the one hand, to provide the used code in the article Martina et al. (yyyy), on the other hand, to faciliate the use of the FLEXPART-WRF model to investigate the potential transport of tracers from the surface(stratosphere) to the stratosphere(surface).
**Herafter, a summary of the main changes is reported.**
*****************************************************************************
## NEW MODULE
   ### intrus_mod.f90
  It is a module that contains all the new introduced variables used in the AirIntrusions Package.
## MODIFIED ORIGINAL SUBROUTINES
   ### ai1_calcpar.f90 and ai1_calcpar_nest.f90
  In these subroutine the computation of the Cold Point Tropopause was added to the original version, which computed only the Lapse Rate Tropopause. 
   ### ai1_readinput.f90
   The module "ai1_intrus_mod" has been included in the subroutine.  
   Several new variables have been added among those variables that must be read from the input file "flexwrf.input", hereafter they are listed:  
   - 5
   ### releaseparticles_reg.f90  
   ### releaseparticles_irrreg.f90  
   ### makefile.mom
   
## NEW SUBROUTINES
  ### 1 INTRUSION SUBROUTINES
   ***Aim:*** The main goal of this package of subroutines is to identify those air parcels entering the tropopause layer/stratosphere from the Planetary Bounary Layer (PBL)/Free Atmosphere (FA) or the PBL/FA from the             tropopause/stratosphere.For each identified particles, the computation of the following varibables is performed:
   - the **transition time** (namely, the time rquired by the particles to reach the tropopause/stratosphere from the Planetary Boundary Layer (PBL);  
   - the **residence time** in the tropopause/stratosphere;  
   - other **thermodynamic variables** (e.g., Brunt-Vaisala frequency) to characterize the transport into the tropopause/stratosphere.
     
   ***Principal subroutines:***
   - 1 *ai1_boundaries_computation.f90*
   - 2 *ai1_thermodyn_computation.f90*
   - 3 *ai1_up_intrusion_identifier.f90*

  ***Rationale:*** In order to identify those air parcels able to penetrate into the tropopause/stratosphere from the PBL/FA(or the other way round), their vertical positions are compared to the tropopause/pbl boundaries at                     each internal time step.
  ### ai1_boundaries_computation.f90
  This subroutine computes the PBL height and the tropopause boundaries. The tropopause can be defined according to two different method:  
  - 1 The bottom and top boundaries are defined considering the Lapse Rate Tropopause and Cold Point Tropopause respectively.
  - 2 The two vertical boundaries are defined as pressure surface following the results reported in Fueglistaler et al. (2009); bottom bounday fix to 150 hPa, while the top one to 70 hPa.

  ### ai1_thermodyn_computation.f90
  This subroutine computes several thermodynamic variables at the particle position; in our study it was used to calculate the values  at the particle position before it enters the tropopause/stratosphere (or the PBL/FA).  
  The computed variables are:
  - Potential Vorticity;
  - Vertical Wind Speed;
  - Turbulent Kinetic Energy;
  - Bulk Richardson Number (Ri);
  - Brunt-Vaisala Frequency (N);
  - Pressure Reference fro (Ri and N);
  - Ellrond Index.
  Furthermore, it is computed also the tropography height at the particle position.

  ### ai1_up_intrusion_identifier.f90
  ### 2 CONVECTIVE SUBROUTINES

  ### 3 GENERAL SUBROUTINES
   ***Aim:*** These subroutines are useful, for example, for handling the allocation of the additional variables required to implement the identification process.  
   ***Principal subroutines:***  
   - 1 *ai1_static_allocation.f90*  
   - 2 *ai1_dynamic_allocation.f90*

   ### ai1_static_allocation.f90  
   It handles the static allocation of variable.  
   ### ai1_dynamic_allocation.f90  
   It manages the dynamic allocation of variable.
