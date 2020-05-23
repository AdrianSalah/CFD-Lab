#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

void boundaryvalues(int imax, int jmax, Grid& grid) {
    
    // VELOCITY - Declaration and Initialisation

    static matrix<double> u;
    static matrix<double> v;

    grid.velocity(u, velocity_type::U);
    grid.velocity(v, velocity_type::V);

    static matrix<double> p;
    grid.pressure(p);

    /* ---- Beginning of noslip boundary value implementation [Adrian] */
/*
    for(int i = 0; i < grid.imaxb(); i++){
        for(int j = 0; j < grid.jmaxb(); j++){
            Cell cell = grid.cell(i, j);
            if(cell._cellType == NOSLIP){
                //check for edge cells
                //B_N
                if(cell._nbNorth->_cellType == FLUID and cell._nbSouth->_cellType == NOSLIP and cell._nbEast->_cellType == NOSLIP and cell._nbWest->_cellType == NOSLIP){
                    v.at(i).at(j) = 0;
                    u.at(i-1).at(j) = -u.at(i-1).at(j+1);
                    u.at(i).at(j) = -u.at(i).at(j+1);
                    p.at(i).at(j) = p.at(i).at(j+1);
                }
                //B_W
                else if(cell._nbSouth->_cellType == FLUID and cell._nbNorth->_cellType == NOSLIP and cell._nbWest->_cellType == NOSLIP and cell._nbEast->_cellType == NOSLIP) {
                    u.at(i - 1).at(j) = 0;
                    v.at(i).at(j - 1) = -v.at(i - 1).at(j - 1);
                    v.at(i).at(j) = -v.at(i - 1).at(j);
                    p.at(i).at(j) = p.at(i - 1).at(j);
                }
                //B_E ...
                //B_S ...


                //check for border cells
                //B_NE
                if(cell._nbNorth->_cellType == FLUID and  cell._nbEast->_cellType == FLUID and cell._nbEast->_cellType == NOSLIP and cell._nbWest->_cellType == NOSLIP){
                    u.at(i).at(j) = 0;
                    v.at(i).at(j) = 0;
                    u.at(i-1).at(j) = -u.at(i-1).at(j+1);
                    v.at(i).at(j-1) = -v.at(i+1).at(j-1);
                    p.at(i).at(j) = (p.at(i).at(j+1) + p.at(i+1).at(j)) / 2;
                }
                //B_NW ...
                //B_SE ...
                //B_SW ...
            }



        }
    }
}


*/

    /* ------ implementation of boundary values ws1 ------ */

    /*


    // -----Boundary conditions initializaion----- //
    
    // ---BOTTOM--- //
    for (int i = 0; i < grid.imaxb(); i++) {
        u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
        v_velocity.at(i).at(0) = 0;
    }

    // ---LEFT--- //
    for (int j = 0; j < grid.jmaxb(); j++) {
        u_velocity.at(0).at(j) = 0;
        v_velocity.at(0).at(j) = -v_velocity.at(1).at(j);
    }

    // ---RIGHT--- //
    for (int j = 0; j < grid.jmaxb(); j++) {
        u_velocity.at(imax + 1).at(j) = 0;
        v_velocity.at(imax + 1).at(j) = -v_velocity.at(imax).at(j);
    }

    // Note: this is related to internal most right cells of the physical domain
    for (int j = 1; j <= jmax; j++) {
        u_velocity.at(imax).at(j) = 0;
    }

    // ---TOP--- //
    // Neuman boundary conditions u_wall = 1  ---> multiplication factor x2 from interpolation formula
    // u_wall may be later provided in the input txt read file

    for (int i = 0; i < grid.imaxb(); i++) {
        v_velocity.at(i).at(jmax + 1) = 0;
        u_velocity.at(i).at(jmax + 1) = 2 - u_velocity.at(i).at(jmax);
    }
    // Note: this is related to internal most upper cells of the physical domain
    for (int i = 1; i <= imax; i++)
        v_velocity.at(i).at(jmax) = 0;

    grid.set_velocity(u_velocity, velocity_type::U);
    grid.set_velocity(v_velocity, velocity_type::V);



    // PRESSURE - Declaration and Initialisation
    // Neuman boundary conditions

    static matrix<double> pres;
    grid.pressure(pres);

    for (int j = 1; j <= jmax; j++)
    {
        pres.at(0).at(j) = pres.at(1).at(j);                // LEFT
        pres.at(imax + 1).at(j) = pres.at(imax).at(j);      // RIGHT
    }

    for (int i = 1; i <= imax; i++)
    {
        pres.at(i).at(0) = pres.at(i).at(1);                // BOTTOM
        pres.at(i).at(jmax + 1) = pres.at(i).at(jmax);      // TOP
    }

    grid.set_pressure(pres);

 // TEMPERATURE - Declaration and Initialisation
 // Neuman boundary conditions

    static matrix<double> temp;
    grid.temperature(temp);

    for (int j = 1; j <= jmax; j++)
    {
        temp.at(0).at(j) = temp.at(1).at(j);                // LEFT
        temp.at(imax + 1).at(j) = temp.at(imax).at(j);      // RIGHT
    }

    for (int i = 1; i <= imax; i++)
    {
        temp.at(i).at(0) = temp.at(i).at(1);                // BOTTOM
        temp.at(i).at(jmax + 1) = temp.at(i).at(jmax);      // TOP
    }

    grid.set_temperature(temp);
    /*
}
