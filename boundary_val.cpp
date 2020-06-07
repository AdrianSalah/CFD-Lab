#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

void boundaryvalues(int imax,
                    int jmax,
                    Grid& grid,
                    double& v_inflow,
                    double& u_inflow,
                    matrix<double>& F,
                    matrix<double>& G,
                    double& T_h,
                    double& T_c,
                    double& dx,
                    double& dy,
                    double &kappa,
                    double &heat_flux,
                    double &beta,
                    double &delta_t,
                    double &GX,
                    double &GY,
                    int scenarioSpec) {
    
    // VELOCITY - Declaration and Initialisation

    static matrix<double> u_velocity;
    static matrix<double> v_velocity;

    grid.velocity(u_velocity, velocity_type::U);
    grid.velocity(v_velocity, velocity_type::V);

    // PRESSURE - Declaration and Initialisation
    static matrix<double> pres;
    grid.pressure(pres);

    // Temperature - Declaration and Initialisation
    static matrix<double> temp;
    grid.temperature(temp);

    // ----- Boundary conditions NO SLIP inner cells ----- //

    for (int i = 1; i < grid.imaxb() - 1; i++) {
        for (int j = 1; j < grid.jmaxb() - 1; j++) {
            // NO slip boundary confitions
            if (grid.cell(i, j)._cellType == NOSLIP) {
                //B_NE
                if (grid.cell(i, j)._nbNorth->_cellType > 1  && grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1); // Was commented by Issa, but is actually reasonable
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i+1).at(j - 1); // Was commented by Issa, but is actually reasonable
                    pres.at(i).at(j) = (pres.at(i).at(j + 1) + pres.at(i + 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j + 1) + temp.at(i + 1).at(j)) / 2;
                    //F[i][j] = u_velocity.at(i).at(j);
                    //G[i][j] = v_velocity.at(i).at(j);
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
                //B_NW
                else if (grid.cell(i, j)._nbNorth->_cellType > 1 && grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0; // Was commented by Issa, but is actually reasonable
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j + 1);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i - 1).at(j - 1); // Was commented by Issa, but is actually reasonable
                    pres.at(i).at(j) = (pres.at(i).at(j + 1) + pres.at(i - 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j + 1) + temp.at(i - 1).at(j)) / 2;
                    F[i - 1][j] = u_velocity.at(i - 1).at(j) - beta * delta_t / 2 * (temp.at(i - 1).at(j) + temp.at(i).at(j)) * GX;
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
                //B_SW
                else if (grid.cell(i, j)._nbSouth->_cellType > 1 && grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0; // Was commented by Issa, but is actually reasonable
                    v_velocity.at(i).at(j - 1) = 0; // Was commented by Issa, but is actually reasonable
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j - 1); // Here was a BIG ISSUE potentially
                    v_velocity.at(i).at(j) = -v_velocity.at(i - 1).at(j);
                    pres.at(i).at(j) = (pres.at(i).at(j - 1) + pres.at(i - 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j - 1) + temp.at(i - 1).at(j)) / 2;
                    F[i - 1][j] = u_velocity.at(i - 1).at(j) - beta * delta_t / 2 * (temp.at(i - 1).at(j) + temp.at(i).at(j)) * GX;
                    G[i][j - 1] = v_velocity.at(i).at(j - 1) - beta * delta_t / 2 * (temp.at(i).at(j - 1) + temp.at(i).at(j)) * GY;
                }
                //B_SE
                else if (grid.cell(i, j)._nbSouth->_cellType > 1 && grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j - 1) = 0; // Was commented by Issa, but is actually reasonable
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j - 1); // Was commented by Issa, but is actually reasonable
                    v_velocity.at(i).at(j) = -v_velocity.at(i + 1).at(j);
                    pres.at(i).at(j) = (pres.at(i).at(j - 1) + pres.at(i + 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j - 1) + temp.at(i + 1).at(j)) / 2;
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                    G[i][j - 1] = v_velocity.at(i).at(j - 1) - beta * delta_t / 2 * (temp.at(i).at(j - 1) + temp.at(i).at(j)) * GY;
                }
                //B_N
                else if (grid.cell(i, j)._nbNorth->_cellType > 1 ) {
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1); // Was commented by Issa, but is actually reasonable
                    u_velocity.at(i).at(j) = - u_velocity.at(i).at(j + 1);
                    pres.at(i).at(j) = pres.at(i).at(j + 1);
                    temp.at(i).at(j) = temp.at(i).at(j + 1);
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
                //B_E
                else if (grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i + 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i + 1).at(j - 1); // Was commented by Issa, but is actually reasonable
                    pres.at(i).at(j) = pres.at(i + 1).at(j);
                    temp.at(i).at(j) = temp.at(i + 1).at(j);
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                }
                //B_S
                else if (grid.cell(i, j)._nbSouth->_cellType > 1) {
                    v_velocity.at(i).at(j - 1) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j - 1); // Was commented by Issa, but is actually reasonable
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j - 1);
                    pres.at(i).at(j) = pres.at(i).at(j - 1);
                    temp.at(i).at(j) = temp.at(i).at(j - 1);
                    G[i][j - 1] = v_velocity.at(i).at(j - 1) - beta * delta_t / 2 * (temp.at(i).at(j - 1) + temp.at(i).at(j)) * GY;
                }
                //B_W
                else if (grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i - 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i - 1).at(j - 1); // Was commented by Issa, but is actually reasonable
                    pres.at(i).at(j) = pres.at(i - 1).at(j);
                    temp.at(i).at(j) = temp.at(i - 1).at(j);
                    F[i - 1][j] = u_velocity.at(i - 1).at(j) - beta * delta_t / 2 * (temp.at(i - 1).at(j) + temp.at(i).at(j)) * GX;
                }
                else if(grid.cell(i, j)._nbEast->_cellType == NOSLIP && grid.cell(i, j)._nbWest->_cellType == NOSLIP&&
                    grid.cell(i, j)._nbNorth->_cellType == NOSLIP&& grid.cell(i, j)._nbSouth->_cellType == NOSLIP) {
                    u_velocity.at(i).at(j) = 0; // EAST
                    u_velocity.at(i - 1).at(j) = 0; // WEST
                    v_velocity.at(i).at(j) = 0; // NORTH
                    v_velocity.at(i).at(j - 1) = 0; // SOUTH

                    pres.at(i).at(j)=0;
                    temp.at(i).at(j)=0;
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
            }
        }
    }

    /* ---- BC outer cells ---- */

    // bottom and top
    for(int i = 0; i < grid.imaxb(); i++){
        //noslip bottom: checked and added pres and temp [Adrian]
        if(grid.cell(i,0)._cellType == NOSLIP and grid.cell(i,0)._nbNorth->_cellType == FLUID) {
            u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
            v_velocity.at(i).at(0) = 0;

            //G.at(i).at(0) = v_velocity.at(i).at(0);
            G.at(i).at(0) = v_velocity.at(i).at(0) - beta * delta_t / 2 * (temp.at(i).at(0) + temp.at(i).at(1)) * GY;
            pres.at(i).at(0) = pres.at(i).at(1);
            temp.at(i).at(0) = temp.at(i).at(1);

        }
        //noslip top: checked and added pres and temp [Adrian]
        if(grid.cell(i, grid.jmaxb() - 1)._cellType == NOSLIP and grid.cell(i, grid.jmaxb() - 1)._nbSouth->_cellType == FLUID) {
            u_velocity.at(i).at(grid.jmaxb() - 1) = -u_velocity.at(i).at(grid.jmaxb() - 2);
            v_velocity.at(i).at(grid.jmaxb() - 1) = 0; // Solid cell
            v_velocity.at(i).at(grid.jmaxb() - 2) = 0; // Fluid cell below

            G.at(i).at(grid.jmaxb() - 1) = v_velocity.at(i).at(grid.jmaxb() - 1); // Solid cell
            

            //G.at(i).at(grid.jmaxb() - 2) = v_velocity.at(i).at(grid.jmaxb() - 2);

            G.at(i).at(grid.jmaxb() - 2) = v_velocity.at(i).at(grid.jmaxb() - 2)
                - beta * delta_t / 2 * (temp.at(i).at(grid.jmaxb() - 2) + temp.at(i).at(grid.jmaxb() - 1)) * GY; // Fluid cell below
            
            pres.at(i).at(grid.jmaxb() - 1) = pres.at(i).at(grid.jmaxb() - 2);
            temp.at(i).at(grid.jmaxb() - 1) = temp.at(i).at(grid.jmaxb() - 2);
        }
    }

    // left and right
    for(int j = 0; j < grid.jmaxb(); j++){
        //noslip left: checked and added pres and temp [Adrian]
        if(grid.cell(0,j)._cellType == NOSLIP and grid.cell(0,j)._nbEast->_cellType == FLUID) {
            u_velocity.at(0).at(j) = 0;
            v_velocity.at(0).at(j) = -v_velocity.at(1).at(j);

            //F.at(0).at(j) = u_velocity.at(0).at(j);

            F.at(0).at(j) = u_velocity.at(0).at(j) - beta * delta_t / 2 * (temp.at(0).at(j) + temp.at(1).at(j)) * GX;
            pres.at(0).at(j) = pres.at(1).at(j);
            temp.at(0).at(j) = temp.at(1).at(j);
        }

        //noslip right: checked and added pres and temp [Adrian]
        if(grid.cell(grid.imaxb() - 1,j)._cellType == NOSLIP and grid.cell(grid.imaxb() - 1, j)._nbWest->_cellType == FLUID) {
            u_velocity.at(grid.imaxb() - 1).at(j) = 0; // Solid cell
            u_velocity.at(grid.imaxb() - 2).at(j) = 0; // Fluid cell
            v_velocity.at(grid.imaxb() - 1).at(j) = -v_velocity.at(grid.imaxb() - 2).at(j);

            F.at(grid.imaxb() - 1).at(j) = u_velocity.at(grid.imaxb() - 1).at(j); // Solid cell (no need to calculate it)

            //F.at(grid.imaxb() - 2).at(j) = u_velocity.at(grid.imaxb() - 2).at(j);

            F.at(grid.imaxb() - 2).at(j) = u_velocity.at(grid.imaxb() - 2).at(j)
                - beta * delta_t / 2 * (temp.at(grid.imaxb() - 2).at(j) + temp.at(grid.imaxb() - 1).at(j)) * GX; // Fluid cell (still no need to calculate it)
            
            pres.at(grid.imaxb() - 1).at(j) = pres.at(grid.imaxb() - 2).at(j);
            temp.at(grid.imaxb() - 1).at(j) = temp.at(grid.imaxb() - 2).at(j);
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
        if(grid.cell(grid.imaxb() - 1, j)._cellType == OUTFLOW and grid.cell(grid.imaxb() - 1, j)._nbWest->_cellType == FLUID){
            v_velocity.at(grid.imaxb() - 2).at(j) = v_velocity.at(grid.imaxb() - 1).at(j); //neumann BC also for v_velocity!
            u_velocity.at(grid.imaxb() - 2).at(j) = u_velocity.at(grid.imaxb() - 1).at(j);
        }
    }

    if(scenarioSpec == 5 or scenarioSpec == 6) {

        /* ----  Dirichlet BC Temperature ---- */
        // Natural Convection and Fluid Trap
        for (int j = 1; j < grid.jmaxb() - 1; j++) { //check indexing! // MODIFIED THE CODE
            //T_h left wall
            temp.at(0).at(j) = 2 * T_h - temp.at(1).at(j);
            //T_c right wall
            temp.at(grid.imaxb() - 1).at(j) = 2 * T_c - temp.at(grid.imaxb() - 2).at(j);
        }

        /* ---- Neumann BC Temperature ---- */
        // Natural Convection and Fluid Trap
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            //bottom wall
            temp.at(i).at(0) = temp.at(i).at(1) + dy * heat_flux / kappa;
            //top wall
            temp.at(i).at(grid.jmaxb() - 1) = temp.at(i).at(grid.jmaxb() - 2) + dy * heat_flux / kappa;
        }
    }
        // Rayleigh-Benard Convection
    else if (scenarioSpec == 7) {

        for (int i = 1; i < grid.imaxb() - 1; i++) {
            //T_h bottom wall
            temp.at(i).at(0) = 2 * T_h - temp.at(i).at(1);
            //T_c top wall
            temp.at(i).at(grid.jmaxb() - 1) = 2 * T_c - temp.at(i).at(grid.jmaxb() - 2);
        }

        // Rayleigh-Benard Convection
        for (int j = 1; j < grid.jmaxb() - 1; j++) {
            //left wall
            temp.at(0).at(j) = temp.at(1).at(j) + dx * heat_flux / kappa;
            //right wall
            temp.at(grid.imaxb() - 1).at(j) = temp.at(grid.imaxb() - 2).at(j) + dx * heat_flux / kappa;
        }
    }

        /* ---- Freeslip BC for lid driven cavity ---- */
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            if (grid.cell(i, grid.jmaxb() - 1)._cellType == FREESLIP) {
                v_velocity.at(i).at(grid.jmaxb() - 1) = 0;
                u_velocity.at(i).at(grid.jmaxb() - 1) = 2 - u_velocity.at(i).at(grid.jmaxb() - 2);
            }

        }



    grid.set_velocity(u_velocity, velocity_type::U);
    grid.set_velocity(v_velocity, velocity_type::V);

    grid.set_pressure(pres);

    grid.set_temperature(temp);
}


void spec_boundary_val(Grid& grid, 
        double &u_inflow,
        double &v_inflow,
        double& T_c,
        double& T_h,
        double &kappa,
        double &heat_flux,
        matrix<double>& U,
        matrix<double>& V,
        matrix<double>& P,
        matrix<double>& T,   
        matrix<double>& F,
        matrix<double>& G,
        int il,
        int ir,
        int jb,
        int jt)
{
        grid.velocity(U, velocity_type::U, il, ir, jb, jt);
        grid.velocity(V, velocity_type::V, il, ir, jb, jt);
        grid.pressure(P, il, ir, jb, jt);
        grid.temperature(T, il, ir, jb, jt);
        // left boundary
        if (il == 0) {
            for (int j = jb+1; j < jt; j++) {
                U[2][j+1] = 0;
                V[1][j+2] = -V[2][j+2];
                F[2][j+1] = U[2][j + 1];
                G[1][j+2] = V[1][j + 2];
                P[1][j+1] = P[2][j+1];
                T[1][j+1] = T[2][j+1];
            }
        }
        //right boundary
        if (ir == grid.imaxb()-1) {
            for (int j = jb+1; j < jt; j++) {
                U[ir+2][j+1] = 0;
                V[ir+1][j+2] = -V[ir][j+2];
                F[ir+2][j+1] = U[ir+2][j+1];
                G[ir+1][j+2] = V[ir+1][j+2];
                P[ir+1][j+1] = P[ir][j+1];
                T[ir+1][j+1] = T[ir][j+1];
            }
        }
        //bottom boundary
        if (jb == 0) {
            for (int i = il+1; i < ir; i++) {
                V[i+1][2] = 0;
                U[i+2][1] = -U[i+2][2];
                F[i+2][1] = U[i+2][1];
                G[i+1][2] = V[i+1][2];
                P[i+1][1] = P[i+1][2];
                T[i+1][1] = T[i+1][2];
            }
        }
        //top boundary free slip
        if (jt == grid.jmaxb()-1) {
            for (int i = il+1; i < ir; i++) {
                V[i+1][jt+1] = 0;
                U[i+2][jt+1] = 2-U[i+2][jt];
                F[i+2][jt+1] = U[i+2][jt+1];
                G[i+1][jt+2] = V[i+1][jt+2];
                P[i+1][jt+1] = P[i+1][jt];
                T[i+1][jt+1] = T[i+1][jt];
            }
        }
        grid.set_velocity(U, velocity_type::U, il, ir, jb, jt);
        grid.set_velocity(V, velocity_type::V, il, ir, jb, jt);
        grid.set_pressure(P, il, ir, jb, jt);
        grid.set_temperature(T, il, ir, jb, jt);
}

void assign_ptr_nbcells(Grid &grid){

    //store pointers to neighbours for inner cells
    for (int i = 1; i < grid.imaxb() - 1; i++) {
        for(int j = 1; j < grid.jmaxb() - 1; j++) {
            grid.cell(i, j)._nbNorth = &grid.cell(i, j + 1);
            grid.cell(i, j)._nbEast = &grid.cell(i + 1, j);
            grid.cell(i, j)._nbWest = &grid.cell(i - 1, j);
            grid.cell(i, j)._nbSouth = &grid.cell(i, j - 1);

            // Incrementing the quantity of FLUID cells
            if (grid.cell(i, j)._cellType == FLUID)
                grid.increment_fluid_cells();
        }
    }
    //neighbour edges
    //bottom left
    grid.cell(0,0)._nbNorth = &grid.cell(0,1);
    grid.cell(0,0)._nbEast = &grid.cell(1,0);
    //top right
    grid.cell(grid.imaxb() - 1, grid.jmaxb() - 1)._nbSouth = &grid.cell(grid.imaxb() - 1,grid.jmaxb() - 2);
    grid.cell(grid.imaxb() - 1, grid.jmaxb() - 1)._nbWest = &grid.cell(grid.imaxb() - 2,grid.jmaxb() - 1);
    //top left
    grid.cell(0,grid.jmaxb() - 1)._nbEast = &grid.cell(1, grid.jmaxb() - 1);
    grid.cell(0, grid.jmaxb() - 1)._nbSouth = &grid.cell(0, grid.jmaxb() - 2);
    //bottom right
    grid.cell(grid.imaxb() - 1, 0)._nbNorth = &grid.cell(grid.imaxb() - 1, 1);
    grid.cell(grid.imaxb() - 1, 0)._nbWest = &grid.cell(grid.imaxb() - 2, 0);

    for(int i = 1; i < grid.imaxb() - 1; i++){
        //bottom
        grid.cell(i, 0)._nbNorth = &grid.cell(i, 1);
        grid.cell(i, 0)._nbEast = &grid.cell(i + 1, 0);
        grid.cell(i, 0)._nbWest = &grid.cell(i - 1, 0);
        //top
        grid.cell(i, grid.jmaxb() - 1)._nbSouth = &grid.cell(i, grid.jmaxb() - 2);
        grid.cell(i, grid.jmaxb() - 1)._nbEast = &grid.cell(i + 1, grid.jmaxb() - 1);
        grid.cell(i, grid.jmaxb() - 1)._nbWest = &grid.cell(i - 1, grid.jmaxb() - 1);
    }
    for(int j = 1; j <  grid.jmaxb() - 1; j++){
        //left
        grid.cell(0, j)._nbNorth = &grid.cell(0, j + 1);
        grid.cell(0, j)._nbSouth = &grid.cell(0, j - 1);
        grid.cell(0, j)._nbEast = &grid.cell(1, j);
        //right
        grid.cell(grid.imaxb() - 1, j)._nbNorth = &grid.cell(grid.imaxb() - 1, j + 1);
        grid.cell(grid.imaxb() - 1, j)._nbSouth = &grid.cell(grid.imaxb() - 1, j - 1);
        grid.cell(grid.imaxb() - 1, j)._nbWest = &grid.cell(grid.imaxb() - 2, j);
    }
}
void boundary_val_sor(Grid& grid) {
    static matrix<double> pres;
    grid.pressure(pres);

    // Iterate over internal cells
    for (int i = 1; i < grid.imax() - 1; ++i) {
        for (int j = 1; j < grid.jmax() - 1; ++j) {
            // Setting pressure for all boundary (NOSLIP) cells
            if (grid.cell(i, j)._cellType == NOSLIP) {

                // Checking whether the neighbor cells belongs to FLUID, and taking its pressure

                // NORTH-EAST-ERN-cell
                if (((grid.cell(i, j)._nbNorth)->_cellType == FLUID) &&
                    ((grid.cell(i, j)._nbEast)->_cellType == FLUID))
                {
                    pres[i][j] = ((grid.cell(i, j)._nbNorth)->pressure() + (grid.cell(i, j)._nbEast)->pressure())/2;
                }

                // SOUTH-EAST-ERN-cell
                else if (((grid.cell(i, j)._nbSouth)->_cellType == FLUID) &&
                    ((grid.cell(i, j)._nbEast)->_cellType == FLUID))
                {
                    pres[i][j] = ((grid.cell(i, j)._nbSouth)->pressure() + (grid.cell(i, j)._nbEast)->pressure())/2;
                }

                // NORTH-WEST-ERN-cell
                else if (((grid.cell(i, j)._nbNorth)->_cellType == FLUID) &&
                    ((grid.cell(i, j)._nbWest)->_cellType == FLUID))
                {
                    pres[i][j] = ((grid.cell(i, j)._nbNorth)->pressure() + (grid.cell(i, j)._nbWest)->pressure()) / 2;
                }

                // SOUTH-WEST-ERN-cell
                else if (((grid.cell(i, j)._nbSouth)->_cellType == FLUID) &&
                    ((grid.cell(i, j)._nbWest)->_cellType == FLUID))
                {
                    pres[i][j] = ((grid.cell(i, j)._nbSouth)->pressure() + (grid.cell(i, j)._nbWest)->pressure()) / 2;
                }


                // NORTHERN-cell
                else if ((grid.cell(i, j)._nbNorth)->_cellType == FLUID)
                    pres[i][j] = (grid.cell(i, j)._nbNorth)->pressure() ;

                // EASTERN-cell
                else if ((grid.cell(i, j)._nbEast)->_cellType == FLUID)
                    pres[i][j] = (grid.cell(i, j)._nbEast)->pressure();

                // SOUTHERN-cell
                else if ((grid.cell(i, j)._nbSouth)->_cellType == FLUID)
                    pres[i][j] = (grid.cell(i, j)._nbSouth)->pressure();

                // WESTERN-cell (cowboy)
                else if ((grid.cell(i, j)._nbWest)->_cellType == FLUID)
                    pres[i][j] = (grid.cell(i, j)._nbWest)->pressure();

            }
        }
    }
    grid.set_pressure(pres);
}