# FLEXPART-WRF - Molave paper
In this repository you can find a modified version of the FLEXPART-WRF model (Brioude et al., 2013) used in the article Martina et al. (2026), "Typhoon, orography and gravity waves - a recipe for rapid transport from the boundary layer to the strato.",xxxx.  
The original code hase been modified only where strictly necessary in order to facilitate the implementation of the changes for potentially new users. Therefore, the changes have been organized in additional subroutines added to the native code.
**Herafter, a summary of the main changes is reported.**
*****************************************************************************
## 1 INTRUSION SUBROUTINES
   ***Aim:*** The main goal of this package of subroutines is to identify those air parcels entering the tropopause layer/stratosphere from the PBL/Free Atmosphere or the PBL/Free Atmosphere from the tropopause/stratosphere.        For each identified particles, the computation of the following varibables is performed:
   - the **transition time** (namely, the time rquired by the particles to reach the tropopause/stratosphere from the Planetary Boundary Layer (PBL);  
    - the **residence time** in the tropopause/stratosphere;  
   - other **thermodynamic variables** (e.g., Brunt-Vaisala frequency) to characterize the transport into the tropopause/stratosphere.
     
   ***Principal subroutines:***
   - 1 boundaries_computation.f90
   - 2....
   - 3...
  ### boundaries_computation.f90
  
     

