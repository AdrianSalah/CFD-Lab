#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

void boundaryvalues(int imax, int jmax, Grid& grid, matrix<double>& F,
    matrix<double>& G) {
    
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

    /* ------ implementation of boundary values ws1 ------ */

 

    static matrix<double> pres;
    grid.pressure(pres);


    static matrix<double> temp;
    grid.temperature(temp);
    // -----Boundary conditions initializaion----- //
    
    for (int i = 1; i < grid.imaxb()-1; i++) {
        for (int j = 1; j < grid.jmaxb()-1; j++) {
            // NO slip boundary confitions
            if (grid.cell(i, j)._cellType == NOSLIP) {
                //B_NE
                if (grid.cell(i, j)._nbNorth->_cellType > 1  && grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i+1).at(j - 1);
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
                    pres.at(i).at(j) = (pres.at(i).at(j + 1) + pres.at(i -1).at(j)) / 2;
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
                else if (grid.cell(i, j)._nbNorth->_cellType > 1 ) {
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1);
                    u_velocity.at(i).at(j) = - u_velocity.at(i).at(j + 1);
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


    grid.set_velocity(u_velocity, velocity_type::U);
    grid.set_velocity(v_velocity, velocity_type::V);

    grid.set_pressure(pres);


    grid.set_temperature(temp);
    
}
