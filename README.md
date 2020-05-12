# CFD Lab code of Group C


This repository contains:

* an example input file (`cavity100.dat`), with `100` stemming from the `itermax`
* the headers
* the files with the respective method implementations 
* example files of created visualisations


## Software Requirements

* VTK 7 or higher
* GCC 9 (optional) 


## Changes made to the Code Skeleton

*  used static local variables in sor.cpp
*  gpp9 set to True in CMakeLists.txt file

## Prevention of segmentation fault while reading in parameters from cavity100.dat file 

* use absolute path: 
std::string filename{"/home/someName/Documents/CFD-Lab/cfdlabcodeskeleton/cavity100.dat"};
* use relative path: 
std::string filename{"../cavity100.dat"};
* or modify CMakeList.txt file by adding following code in the end of the file:

file(GLOB CAVITY
         "${CMAKE_CURRENT_SOURCE_DIR}/*.dat"
      )
      file(COPY ${CAVITY} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
      
