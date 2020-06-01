# CFD Lab code of Group C


**To run the code for different scenarios**
*  PLease compile code and then type: `./sim x` in the command line. To run different scenarios please change `x` for an integer that refers to a certain scenario:
* `x = 1` **lid driven cavity**
* `x = 2` **plane shear flow**
* `x = 3` **karman vortex street**
* `x = 4` **flow over step**
* `x = 5` **natrual convection**
* `x = 6` **fluid trap**
* `x = 7` **rayleigh benard convection**

A folder with the respective scenario name is created to store the vtk-files for visualization. 
If no scenario is specified it will run the code for the lid-driven cavity. 
If you pass a parameter `x` which i not equal to 1-7 the programm will terminate with an error.
A check is done if the provided geometry file is solvable. Please note one obstacle cell (marked as 0 in the pgm-file) is only allowed to have at maximum two adjacent fluid neighbour cells (marked as 4 in the pgm-file).
The scenario name, the structure of the geometry, the progress and the temperature is displayed in the terminal. 
If you encounter strange looking visualisation results for the non-heat-transfer cases please change between channels in Paraview, it's due to zero temperature.

**Representation of the cells in the PGM-file**
*  fluid cell: marked 4
*  no-slip cell: marked 0
*  free-slip cell: marked 1
*  inflow cell: marked 3
*  outflw cell: marked 2

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

