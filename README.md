# FLEXPART-WRF - Description of the new "AirIntrusions" package
The AirIntrusions package is a modified version of some of the files forming the FLEXPART-WRF model (Brioude et al., 2013) developed in the article Martina et al. (2026), "The impact of orography on the troposphere-to-stratosphere transport during a typhoon event in the tropics.".  
The original files have been modified only where strictly necessary in order to facilitate the implementation of the changes for potentially new users. Therefore, the changes have been organized in additional subroutines that must be added to the native code.  
In order to use the package, you need to download it and store its files in the same directory where the original FLEXPART-WRF code is hosted. All the files belonging to the AirIntrusions package are named as "ai1_*.f90", where "ai1" stands for "AirIntrusions package version 1".  
Then, you need to use "ai1_make.mom" to compile the modified version of the program (NOTE: you need to remove from the directory all those ".o" files connected to the original version of files that have been modified, such as "timemanager_serial.o" since it has been modified in some part in this new package).  
The aim of this package is, on the one hand, to provide the used code in the article Martina et al. (2026), on the other hand, to faciliate the use of the FLEXPART-WRF model to investigate the potential transport of tracers from the surface to the stratosphere.

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
   
   - *"mquasilag"*, it is used to assign a source ID to each particle in order to associate them to their source. 0 = at each particle is assing an ID corresponding to one source; 1 = no source ID is assigned. NOTE: It is not a new variable of the AirIntrusions pacakge. It was already implemented in FLEXPART-WRF but the read command was not present in the subroutine.
   - *"cape_option"*, it can be 1 or 0. 1 means to use the CAPE values, previously computed, present in the WRF files; 0 means to not consider the CAPE values.
   - *"nresidenceclass"*, it sets the number of residence classes in which divide the results.
   - *"lresidence"*, it is an array containing the boundaries of each residence class (in seconds).
   - *"source type"*, this flag is used to set the location of the surces.  1= soruces in the PBL, 0 sources in the stratosphere (only the option 1 has been tested).
   - *"tropo_method"*, this flag is used to activate one of the two method available to compute the tropopause boundaries. 1=Lapse Rate/Cold Point Tropopause, 0=Pressure surfaces (Fueaflistaler et al. 2009).

**NOTE:** It is fundamental to add this command variables in the input file "flexwrf.input", otherwise you'll get an error during the run.
   ### ai1_releaseparticles_reg.f90 &  ai1_releaseparticles_irreg.f90  
   The module "ai1_intrus_mod" has been included in the subroutine.  
   The variable "npoint_id(jpart)" has been added in order to assign a unique ID to each particle when *"mquasilag"* == 0. This ID is fundamental in order to identify those parcels enetring the tropopause/stratopshere or the surface.
   The variable "tl_count_flag(ipart)" and "strato_count_flag(ipart) have been added and initialized to zero. These variables are used to count the number of times an air parcel enters the tropopause and stratosphere region, respectively. The variable "free_atm_antering_time(ipart) is initialized here to -999; it is used for computing the transition time from the PBL to the tropopause/stratosphere. 
   
   ### ai1_timemenager_serial.f90  
   The module "ai1_intrus_mod" has been included in the subroutine.  
   In the variables declaration, addition of "integer :: ll_index, mm_index". These are array indeces for accessing the "tl_*" and "strato_*" array used in the identification process (see ai1_up_intrusion_identifier.f90).  
   Addition of a section where the *static_allocation_**.f90* subroutin is called to allocate several different arrays.  
   Then, "ll_index" and "mm_index" are initialized to 0 before entering the do loop over the time.  
   Addition of the call to the subroutine *"ai1_up_intrusion_identifier.f90"* and to the subroutine *"ai1_intrusionoutput.f90"*  
   Addition of a deallocation statement to dischard all the new allocated arrays.

   ### ai1_calcfluxes.f90
   No main modifications, just a bug fixing in the definiton of the variables "xold" and "yold". They must be declared as REAL(kind=dp) and not just REAL. 

   ### ai1_outgrid_init_reg.f90 & ai1_outgrid_init_irreg.f90
   The module "ai1_intrus_mod" has been included in the subroutine. 
   The allocation of several different new arrays has been added to these subroutines for handling the intrusions information.

   ### ai1_partoutput.f90
   The module "ai1_intrus_mod" has been included in the subroutine.
   Addition of some variables for computing extra information to characterize the intrusions, such as cold point tropopause, kinetic energy, ellornd index, richardson number, etc.

   ### ai1_readpartpositions.f90
   This subroutine is used to re-start a simulation from a certain time step (warm start), hence FLEXPART-WRF needs to read the particle positions from its output files.
   No main modifications, just a bug fixing in the reading of the header file.

   ### ai1_readwind.f90 & ai1_readwind_nest.f90
   The module "ai1_intrus_mod" has been included in the subroutine.
   Addition of the possibility to read the Turbulent Kinetic Energy and CAPE diretly from the Weather and Research Forecasting (WRF) model.

   ### ai1_verttransform.f90 & ai1_verttransform_nests.f90
   The module "ai1_intrus_mod" has been included in the subroutine.
   Addition of the vertical interpolation of the CAPE variable.
   
   ### makefile.mom  
   In the definition of the "OBJECTS" variable, all the "ai1_*.o" object files have been included in order to produce the modified version of FLEXPART-WRF.  
   
## NEW SUBROUTINES - INTRUSION SUBROUTINES
***Aim:*** The main goal of this package of subroutines is to identify those air parcels entering the tropopause layer/stratosphere from the Planetary Bounary Layer (PBL)/Free Atmosphere (FA) or the PBL/FA from the             tropopause/stratosphere (the latter is still in development).For each identified particles, the computation of the following variables is performed:
   - the **transition time** (namely, the time required by the particles to reach the tropopause/stratosphere from the Planetary Boundary Layer (PBL);  
   - the **residence time** in the tropopause/stratosphere;  
   - other **thermodynamic variables** (e.g., Brunt-Vaisala frequency) to characterize the transport into the tropopause/stratosphere.
     
   ***Principal subroutines:***
   1) *ai1_boundaries_computation.f90*;
   2) *ai1_thermodyn_computation.f90*;
   3) *ai1_up_intrusion_identifier.f90*.
      
  ***Rationale:*** In order to identify those air parcels able to penetrate into the tropopause/stratosphere from the PBL/FA(or the other way round), their vertical positions are compared to the tropopause/pbl boundaries at                     each internal time step.
  ### ai1_boundaries_computation.f90
  This subroutine computes the PBL height and the tropopause boundaries. The tropopause can be defined according to two different method:  
   1) The bottom and top boundaries are defined considering the Lapse Rate Tropopause and Cold Point Tropopause, respectively.
   2) The two vertical boundaries are defined as pressure surface following the results reported in Fueglistaler et al. (2009); bottom bounday fix to 150 hPa, while the top one to 70 hPa.

  ### ai1_thermodyn_computation.f90
  This subroutine computes several thermodynamic variables at the particle position before it enters the tropopause/stratosphere.  
  The computed variables are:
  - Potential Vorticity;
  - Vertical Wind Speed (it needs to be checked);
  - Turbulent Kinetic Energy (TKE);
  - Bulk Richardson Number (Ri);
  - Brunt-Vaisala Frequency (N);
  - Pressure Reference for (Ri and N);
  - Ellrond Index (it needs to be checked).
  Furthermore, it is also computed the tropography height at the particle position.

  ### ai1_up_intrusion_identifier.f90  
  This subroutine performs the identification of those parcels entering the tropopause/stratosphere.  
  An important aspect is that the used array in the identification process are dynamically allocated (increasing their size of 100 spaces everytime it is reached their end) in order to reduce the computational cost.  
  This subroutine required two additional subroutines: *ai1_boundaries_computation.f90* and *ai1_thermodyn_computation.f90*.  
  The subroutine handles two types of intrusion:

   1) intrusion into the tropopause Layer (TL) (defined by ai1_boundaries_computation.f90);
   2) intrusion into the stratosphere.

  **Brief explanation of how the intrusion into the TL are identified**  
  First, the time ("free_atm_entering_time(jpart)") at which the parcel cross the Planetary Boundary Layer (PBL) top is computed and updated if the parcel oscillates around the PBL top before moving upward in the free atmosphere.  
  Second, when a particle enters the TL, several variables are associted to it:

  - *tl_count_flag*, it counts the number of times the parcel enters the TL;
  - *tl_particle_id*, it saves the ID associated to the parcel;
  - *tl_entering_time*, it saves the entering time of the parcel into the TL and it will be useful to compute both its residence time in this layer and the transition time from the PBL;
  - *tl_flux_flag*, it is used to select only those parcel entering the TL when the results are written out. If the parcel enters the TL, *tl_flux_flag* is equal tp +1. Otherwise, it is equal to -999. (see *ai1_intrusionoutput.f90*);
  - *tl_transition_time*, it computes the time required by the parcel to reach the TL from the PBL top and it is calculated as *tl_entering_time* - *free_atm_entring_time*.
    
Then, the subroutine *ai1_thermodyn_computation.f90* is called in order to compute the thermodynamic variables associated to this parcel before its entrance in the TL.  
Third, the computation of the parcel residence time is performed considering different scenarios:
 1) the parcel continues moving upward and it enters the stratosphere;
 2) the parcel returns in the troposphere from the TL; 
 4) the parcel remains in the TL until the end of the simulation;
 5) the parcel remains in the TL until it is terminated because it leaves the domain.
    
 **Brief explanation of how the intrusion into the stratosphere are identified**  
 The rationale is very similar to the one described for the TL, hence here only the differences are reported.  
 The computation of the parcel residence time is done considering three different scenarios:  

   1) the parcel returns in the troposphere from the stratosphere;
   2) the parcel remains in the stratosphere until the end of the simulation;
   3) the parcel remains in the stratosphere untile it is terminated.

### ai1_intrusionout.f90  
The intrusion data are saved at the end of the simulation using this subroutine. 
The air parcel intrusions into the TL/stratosphere are counted and divided into two fluxes: "first entrance into the layer" and "secondary entries into the layer"; the subroutine "*ai1_calccrossing.f90*" handles the computation and saves the data into ".txt" files.
The Tl and stratopshere intrusion data are saved hourly in different files:
- .txt files: "ttl_class*.txt" and "strato_class*.txt", the data are saved not on a regular grid;
- binary files: "ttl_residence_class*" and "strato_residence_class*", the data are saved hourly on a regular grid.

### ai1_calccrossing.g90
Here the "tl_class*.txt" and "strato_class*.txt" files are filled in with the data.
Then, the number of intrusions are calculated on a regular output grid; the computation is based on the original subroutine "calclfux.f90" of the FLEXPART-WRF model. Furthermore, the thermodynamical variables are also sum up for each grid point of the regular grid in order to get an avaregd value as output (done by the subroutine "*ai1_crossingoutput.f90*").

### ai1_crossingoutput.f90
The intrusions data are computed and saved on the binary files. For each grid point, the variables are:
- Number of intrusions;
- Mean entering time in the layer;
- Average transition time from the PBL to the layer;
- Average potential vorticity;
- Average vertical wind speed (it needs to be checked);
- Average Turbulent Kinetic Energy;
- Average pressure;
- Average Richardson number;
- Average Ellrond index (it needs to be checked);
- Average Brunt-Vaisala frequency;
- Total mass of the substance considered in the layer.

### ai1_dynamic_allocation_multicross_.f90  
This subroutine is an helper subroutine to dynamically allocate arrays. Specifically, it extends their second dimension that is associated to multiple intrusions into the TL/stratosphere.

### ai1_dynamic_allocation_newentry.f90  
This subroutine is an helper subroutine to dynamically allocate arrays. Specifically, it extends their first dimension that is associated to new intrusions into the TL/stratosphere.

### Acknowledgements. 
All authors were supported by GAČR project "Unravelling Subgrid-Scale Orography Effects on Composition in the Free Atmosphere (SCOPE)", 25-17683S in the final year of the study. Massimo Martina, Petr Šácha were supported by the JUNIOR STAR project n. 23-04921M “Unravelling climate impacts of atmospheric internal gravity waves” as well as Anahí Villalba Pradas in initial years of the study. Petr Šácha was also partly supported by the Charles University Research Centre program no. UNCE/24/SCI/005 and Massimo Martina by the GAUK project n°442325.
