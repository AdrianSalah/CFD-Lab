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
