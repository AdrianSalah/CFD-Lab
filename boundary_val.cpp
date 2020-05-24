#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

void boundaryvalues(int imax, int jmax, Grid& grid, double& v_inflow, double& u_inflow, matrix<double>& F,
    matrix<double>& G, double& TD, double& kappa, double& heat_flux) {
    
    // VELOCITY - Declaration and Initialisation

    static matrix<double> u_velocity;
    static matrix<double> v_velocity;

    grid.velocity(u_velocity, velocity_type::U);
    grid.velocity(v_velocity, velocity_type::V);

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

    /* ------ implementation of boundary values ws2 ------ */

 

    static matrix<double> pres;
    grid.pressure(pres);


    static matrix<double> temp;
    grid.temperature(temp);

    // ----- Boundary conditions NO SLIP inner cells ----- //
    
    for (int i = 1; i < grid.imaxb()-1; i++) {
        for (int j = 1; j < grid.jmaxb() - 1; j++) {
            // NO slip boundary confitions
            if (grid.cell(i, j)._cellType == NOSLIP) {
                //B_NE
                if (grid.cell(i, j)._nbNorth->_cellType > 1 && grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i + 1).at(j - 1);
                    pres.at(i).at(j) = (pres.at(i).at(j + 1) + pres.at(i + 1).at(j)) / 2;
                    F[i][j] = u_velocity.at(i).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                    //B_NW
                else if (grid.cell(i, j)._nbNorth->_cellType > 1 && grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0;
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j + 1);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i - 1).at(j - 1);
                    pres.at(i).at(j) = (pres.at(i).at(j + 1) + pres.at(i - 1).at(j)) / 2;
                    F[i][j] = u_velocity.at(i).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                    //B_SW
                else if (grid.cell(i, j)._nbSouth->_cellType > 1 && grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0;
                    v_velocity.at(i).at(j - 1) = 0;
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j + 1);
                    v_velocity.at(i).at(j) = -v_velocity.at(i - 1).at(j);
                    pres.at(i).at(j) = (pres.at(i).at(j - 1) + pres.at(i - 1).at(j)) / 2;
                    F[i][j] = u_velocity.at(i).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                    //B_SE
                else if (grid.cell(i, j)._nbSouth->_cellType > 1 && grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j - 1) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j - 1);
                    v_velocity.at(i).at(j) = -v_velocity.at(i + 1).at(j);
                    pres.at(i).at(j) = (pres.at(i).at(j - 1) + pres.at(i + 1).at(j)) / 2;
                    F[i][j] = u_velocity.at(i).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                    //B_N
                else if (grid.cell(i, j)._nbNorth->_cellType > 1) {
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1);
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j + 1);
                    pres.at(i).at(j) = pres.at(i).at(j + 1);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                    //B_E
                else if (grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i + 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i + 1).at(j - 1);
                    pres.at(i).at(j) = pres.at(i + 1).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                    //B_S
                else if (grid.cell(i, j)._nbSouth->_cellType > 1) {
                    v_velocity.at(i).at(j - 1) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j - 1);
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j - 1);
                    pres.at(i).at(j) = pres.at(i).at(j - 1);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                    //B_W
                else if (grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i - 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i - 1).at(j - 1);
                    pres.at(i).at(j) = pres.at(i - 1).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
            }
        }
    }

    /* ---- BC outer cells ---- */

    // bottom and top
    for(int i = 0; i <= imax; i++){
        //noslip bottom
        if(grid.cell(i,0)._cellType == NOSLIP and grid.cell(i,0)._nbNorth->_cellType == FLUID) {
            u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
            v_velocity.at(i).at(0) = 0;
        }
        //noslip top
        if(grid.cell(i,jmax)._cellType == NOSLIP and grid.cell(i, jmax)._nbSouth->_cellType == FLUID) {
            u_velocity.at(i).at(jmax) = -u_velocity.at(i).at(jmax - 1);
            v_velocity.at(i).at(jmax) = 0;
        }
    }

    // left and right
    for(int j = 0; j <= jmax; jmax++){
        //noslip left
        if(grid.cell(0,j)._cellType == NOSLIP and grid.cell(0,j)._nbEast->_cellType == FLUID) {
            u_velocity.at(0).at(j) = 0;
            v_velocity.at(0).at(j) = -v_velocity.at(1).at(j);
        }

        //noslip right
        if(grid.cell(imax,j)._cellType == NOSLIP and grid.cell(imax, j)._nbWest->_cellType == FLUID) {
            u_velocity.at(imax + 1).at(j) = 0;
            v_velocity.at(imax + 1).at(j) = -v_velocity.at(imax).at(j);
        }

        //inflow left
        if(grid.cell(0,j)._cellType == INFLOW and grid.cell(0,j)._nbEast->_cellType == FLUID){
                v_velocity.at(0).at(j) = v_inflow;
                u_velocity.at(0).at(j) = u_inflow;
        }

        //inflow from right
        //...

        //outflow to left
        //...

        //outflow to right
        if(grid.cell(imax, j)._cellType == OUTFLOW and grid.cell(imax, j)._nbWest->_cellType == FLUID){
            v_velocity.at(imax-1).at(j) = v_velocity.at(imax).at(j); //neumann BC also for v_velocity!
            u_velocity.at(imax-1).at(j) = u_velocity.at(imax).at(j);
        }

    }



    grid.set_velocity(u_velocity, velocity_type::U);
    grid.set_velocity(v_velocity, velocity_type::V);

    grid.set_pressure(pres);

    grid.set_temperature(temp);
}


void spec_boundary_val(double &u_inflow, double &v_inflow, double& TD, double &kappa, double &heat_flux, double val_u_inflow, double val_v_inflow, double val_TD, double val_kappa, double val_heat_flux) {
    u_inflow = val_u_inflow;
    v_inflow = val_v_inflow;

    TD = val_TD;
    kappa = val_kappa;
    heat_flux = val_heat_flux;

}

// Here I tried to set boundary values for outer cells in spec_boundary_val function -> this will be discarded [Adrian]

    /*
    double v_inflow, u_inflow;

    matrix<double> u_velocity, v_velocity;

    grid.velocity(u_velocity, velocity_type::U);
    grid.velocity(v_velocity, velocity_type::V);


    switch (szenarioNumber) {
        case 1: printf("set outer BC for Plane shear flow");
            // Noslip BC
            for(int i = 0; i < grid.imaxb(); i++){
                //bottom
                u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
                v_velocity.at(i).at(0) = 0;
                //top
                u_velocity.at(i).at(jmax) = -u_velocity.at(i).at(jmax-1);
                v_velocity.at(i).at(jmax) = 0;
            }

            // Inflow BC (from left)
            v_inflow = 0.0;
            //To-Do:
            double delta_y = y_length/(jmax-2);

            for(int j = 1; j < grid.jmaxb()-1; j++){
                v_velocity.at(0).at(j) = v_inflow;
                u_velocity.at(0).at(j) = -0.5*Re*4/x_length * j*delta_y * (j*delta_y - y_length);
            }

            // Outflow BC ??
            //To-Do: check wheter that's correct
            for(int j = 1; j < grid.jmaxb()-1; j++){
                //neumann for v
                u_velocity.at(imax).at(j) = u_velocity.at(imax+1).at(j);
            }

            break;

        case 2: printf("set outer BC for Karman Vortex Street");

            // Noslip BC
            for(int i = 0; i < grid.imaxb(); i++){
                //bottom
                u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
                v_velocity.at(i).at(0) = 0;
                //top
                u_velocity.at(i).at(jmax) = -u_velocity.at(i).at(jmax-1);
                v_velocity.at(i).at(jmax) = 0;
            }

            u_inflow = 1.0;
            v_inflow = 0.0;

            // Inflow BC (from left)
            for(int j = 1; j < grid.jmaxb()-1; j++) {
                v_velocity.at(0).at(j) = v_inflow;
                u_velocity.at(0).at(j) = u_inflow;
            }

            // Outflow BC



            break;

        case 3: printf("BC for flow over step");

            // Noslip BC
            for(int i = 0; i < grid.imaxb(); i++){
                //bottom
                u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
                v_velocity.at(i).at(0) = 0;
                //top
                u_velocity.at(i).at(jmax+1) = -u_velocity.at(i).at(jmax);
                v_velocity.at(i).at(jmax+1) = 0;
            }

            u_inflow = 1.0;
            v_inflow = 0.0;
            break;


    }



    //set boundary values for heat transfer
    TD = val_TD;
    kappa = val_kappa;
    heat_flux = val_heat_flux;
};

 */