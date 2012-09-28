In order to use EDKS you need to load the right environment variables.
In jokull execute:
  source source /home/geomod/dev/edks/setup_jokull.bash
In FRAM:
  module purge
  module load openmpi/gcc44
  source /home/ortega/dev/edks/setup_fram.bash

To compute the Kernels you need to prepare 2 files: the velocity model and the edks
configuration file.

 - Velocity model file: this file has to be named XXXXXX.model and has the following 
   structure
#######
NL UF
Vs_1 Vp_1 rho_1 H_1
Vs_2 Vp_3 rho_3 H_3
:
:
Vs_NL Vp_NL rho_NL 0.0
####
 where :
   - NL is the number of layers (including the half-space below)
   - UF is a unit conversion factor to multiply all the values of Vp Vs rho and H to 
     get them to be in the units described below. (make sure it has a decimal point).
   - Vs : S wave velocity (m/s)
   - Vp : P wave velocity (m/s)
   - rho : density  (kg/m^3) 
   - H : layer thinkness  (m)

  IMPORTANT: make sure there are NO empty lines after the last layer indicating the 
             properties of the half space (the one with 0.0 thickness)


###################

 - EDKS configuration file: this file hast to be named XXXXX.edks_config (NOTE: here 
                            the XXXXX must be the same prefix used in the XXXXX.model
                            file )
   Here is an example of the file:
#######
BIN_tab4E        ${EDKS_HOME}/bin/tab4E
BIN_build_edks   ${EDKS_HOME}/bin/build_edks
ModPrefix        XXXXXXX
dR               2.5
Rmax             2200
depths           2.600  2.977  3.361  3.745  4.129  4.512  4.896  5.280  .... MaxDepth
#######
Here ModPrefix is the same prefix used in the velocity model file.
dR is the spacing of the receivers
Rmax  is the horizontal distance to the farthest receiver.
depth a list of depth in which the elementary sources are computed. IMPORTANT: HERE YOU
CANT PUT A DEPTH = 0 (ZERO) NOR AT THE BOUNDARIES OF THE LAYERS. THE CLOSEST YOU GET TO 
A LAYER, THE LARGER THE COMPUTE TIME AND ALSO INCREASES THE CHANCE OF HAVING NUMERICAL
ISSUES.
  
to run the kernel calculator using 7 parallel processes ( 1 is master process)
type:

mpirun -n 8 MPI_EDKS.py XXXXXX.edks_config

it will display in the screen the files of the depths currently being processed. 

