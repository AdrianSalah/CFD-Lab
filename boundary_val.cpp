#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

void boundaryvalues(int imax, int jmax, Grid& grid) {
    
    // VELOCITY - Declaration and Initialisation

    matrix<double> u_velocity;
    matrix<double> v_velocity;

    grid.velocity(u_velocity, velocity_type::U);
    grid.velocity(v_velocity, velocity_type::V);


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

    matrix<double> pres;
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
}
