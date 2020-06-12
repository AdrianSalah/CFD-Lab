# CFD Lab code of Group C


**How to run this code**

*  Build the program and run it with `mpirun -np x ./sim` with `x` being an integer which specifies the number of processes
*  In the parameter file `iproc` and `jproc` can be set, which specifiy the subdivision of the whole domain in x and y-direction
*  If the number of specified subdomains exceeds the number of the processes available the respective parameters `iproc` and `jproc` are adjusted respectively
*  VTK-files for each process are stored consecutively named: `proc0`, `proc1`, `proc2`, ... following the current timestep 
*  If you want to run the program for more processes than cores avialable on your system, use the option: `mpirun -np x ./sim --oversubscribe`

**This repository contains:**

* directory "geometry" which contains .pgm files which are used to read in the geometry for the different scencarios
* directory "parameters" which contains .dat files which are used to read in the parameters for the different scencarios
* an executable "sim"
* the headers
* the files with the respective method implementations 
* example files of created visualisations: visualization examples


## Software Requirements

* VTK 7 or higher
* GCC 9 (optional) 


## Changes made to the Code Skeleton

*  used static local variables in sor.cpp
*  gpp9 set to True in CMakeLists.txt file

