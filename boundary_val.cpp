#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

void boundaryvalues(int imax, int jmax, Grid& grid, matrix<double> & u_velocity, matrix<double>& v_velocity, matrix<double>& pres) {
    
    // VELOCITY - Declaration and Initialisation


    // -----Boundary conditions initializaion----- //
    
    // ---BOTTOM--- //
    for (int i = 0; i < imax + 2; i++) {
        u_velocity.at(i+2).at(1) = -u_velocity.at(i+2).at(2);
        v_velocity.at(i+1).at(2) = 0;
    }

    // ---LEFT--- //
    for (int j = 0; j < jmax + 2; j++) {
        u_velocity.at(2).at(j+1) = 0;
        v_velocity.at(1).at(j+2) = -v_velocity.at(2).at(j+2);
    }

    // ---RIGHT--- //
    for (int j = 0; j < jmax + 2; j++) {
        u_velocity.at(imax + 3).at(j+1) = 0;
        v_velocity.at(imax + 2).at(j+2) = -v_velocity.at(imax+1).at(j+2);
    }

    // Note: this is related to internal most right cells of the physical domain
    for (int j = 1; j <= jmax; j++) {
        u_velocity.at(imax+2).at(j+1) = 0;
    }

    // ---TOP--- //
    // Neuman boundary conditions u_wall = 1  ---> multiplication factor x2 from interpolation formula
    // u_wall may be later provided in the input txt read file

    for (int i = 0; i < imax + 2; i++) {
        v_velocity.at(i+1).at(jmax + 3) = 0;
        u_velocity.at(i+2).at(jmax + 2) = 2 - u_velocity.at(i+2).at(jmax+1);
    }
    // Note: this is related to internal most upper cells of the physical domain
    for (int i = 1; i <= imax; i++)
        v_velocity.at(i+1).at(jmax+2) = 0;



    // PRESSURE - Declaration and Initialisation
    // Neuman boundary conditions


    for (int j = 1; j <= jmax; j++)
    {
        pres.at(1).at(j+1) = pres.at(2).at(j+1);                // LEFT
        pres.at(imax + 2).at(j+1) = pres.at(imax+1).at(j+1);      // RIGHT
    }

    for (int i = 1; i <= imax; i++)
    {
        pres.at(i+1).at(1) = pres.at(i+1).at(2);                // BOTTOM
        pres.at(i+1).at(jmax + 2) = pres.at(i+1).at(jmax+1);      // TOP
    }

}
