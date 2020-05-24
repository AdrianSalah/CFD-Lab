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
                    double &heat_flux) {
    
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

    // ----- Boundary conditions NO SLIP inner cells ----- //

    for (int i = 1; i < imax-1; i++) {
        for (int j = 1; j < jmax - 1; j++) {
            // NO slip boundary confitions
            if (grid.cell(i, j)._cellType == NOSLIP) {
                //B_NE
                if (grid.cell(i, j)._nbNorth->_cellType > 1  && grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i+1).at(j - 1);
                    pres.at(i).at(j) = (pres.at(i).at(j + 1) + pres.at(i + 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j + 1) + temp.at(i + 1).at(j)) / 2;
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
                    temp.at(i).at(j) = (temp.at(i).at(j + 1) + temp.at(i - 1).at(j)) / 2;
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
                    temp.at(i).at(j) = (temp.at(i).at(j - 1) + temp.at(i - 1).at(j)) / 2;
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
                    temp.at(i).at(j) = (temp.at(i).at(j - 1) + temp.at(i + 1).at(j)) / 2;
                    F[i][j] = u_velocity.at(i).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                //B_N
                else if (grid.cell(i, j)._nbNorth->_cellType > 1 ) {
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1);
                    u_velocity.at(i).at(j) = - u_velocity.at(i).at(j + 1);
                    pres.at(i).at(j) = pres.at(i).at(j + 1);
                    temp.at(i).at(j) = temp.at(i).at(j + 1);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                //B_E
                else if (grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i + 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i + 1).at(j - 1);
                    pres.at(i).at(j) = pres.at(i + 1).at(j);
                    temp.at(i).at(j) = temp.at(i + 1).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                //B_S
                else if (grid.cell(i, j)._nbSouth->_cellType > 1) {
                    v_velocity.at(i).at(j - 1) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j - 1);
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j - 1);
                    pres.at(i).at(j) = pres.at(i).at(j - 1);
                    temp.at(i).at(j) = temp.at(i).at(j - 1);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                //B_W
                else if (grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i - 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i - 1).at(j - 1);
                    pres.at(i).at(j) = pres.at(i - 1).at(j);
                    temp.at(i).at(j) = temp.at(i - 1).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
                else if(grid.cell(i, j)._nbEast->_cellType == NOSLIP && grid.cell(i, j)._nbWest->_cellType == NOSLIP&&
                    grid.cell(i, j)._nbNorth->_cellType == NOSLIP&& grid.cell(i, j)._nbSouth->_cellType == NOSLIP) {
                    pres.at(i).at(j)=0;
                    temp.at(i).at(j)=0;
                    F[i][j] = u_velocity.at(i).at(j);
                    G[i][j] = v_velocity.at(i).at(j);
                }
            }
        }
    }

    /* ---- BC outer cells ---- */

    // bottom and top
    for(int i = 0; i < imax; i++){
        //noslip bottom
        if(grid.cell(i,0)._cellType == NOSLIP and grid.cell(i,0)._nbNorth->_cellType == FLUID) {
            u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
            v_velocity.at(i).at(0) = 0;

            G.at(i).at(0) = v_velocity.at(i).at(0);
        }
        //noslip top
        if(grid.cell(i,jmax-1)._cellType == NOSLIP and grid.cell(i, jmax-1)._nbSouth->_cellType == FLUID) {
            u_velocity.at(i).at(jmax-1) = -u_velocity.at(i).at(jmax - 2);
            v_velocity.at(i).at(jmax-1) = 0;

            G.at(i).at(jmax-1) = v_velocity.at(i).at(jmax-1);
        }
    }

    // left and right
    for(int j = 0; j < jmax; jmax++){
        //noslip left
        if(grid.cell(0,j)._cellType == NOSLIP and grid.cell(0,j)._nbEast->_cellType == FLUID) {
            u_velocity.at(0).at(j) = 0;
            v_velocity.at(0).at(j) = -v_velocity.at(1).at(j);

            F.at(0).at(j) = u_velocity.at(0).at(j);
        }

        //noslip right
        if(grid.cell(imax-1,j)._cellType == NOSLIP and grid.cell(imax-1, j)._nbWest->_cellType == FLUID) {
            u_velocity.at(imax - 1).at(j) = 0;
            v_velocity.at(imax - 1).at(j) = -v_velocity.at(imax-2).at(j);

            F.at(imax-1).at(j) = u_velocity.at(imax-1).at(j);
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
        if(grid.cell(imax-1, j)._cellType == OUTFLOW and grid.cell(imax-1, j)._nbWest->_cellType == FLUID){
            v_velocity.at(imax-2).at(j) = v_velocity.at(imax-1).at(j); //neumann BC also for v_velocity!
            u_velocity.at(imax-2).at(j) = u_velocity.at(imax-1).at(j);
        }
    }

    /* ----  Dirichlet BC Temperature ---- */
    // Natural Convection and Fluid Trap

    for(int j = 1; j < jmax; j++){ //check indexing!
        //T_h left wall
        temp.at(0).at(j) = 2*T_h - temp.at(1).at(j);
        //T_c right wall
        temp.at(imax-1).at(j) = 2*T_c - temp.at(imax-2).at(j);
    }

    // Rayleigh-Benard Convection
    /*
    for(int i = 1; i < imax; i++){
        //T_h bottom wall
        temp.at(i).at(0) = 2*T_h - temp.at(i).at(1);
        //T_c top wall
        temp.at(i).at(jmax-1) = 2*T_c - temp.at(i).at(jmax-2);
    }
     */

    /* ---- Neumann Boundary ---- */

    // Natural Convection and Fluid Trap
    for(int i = 1; i < imax; i++){
        //bottom wall
        temp.at(i).at(0) = temp.at(i).at(1) + dy * heat_flux / kappa;
        //top wall
        temp.at(i).at(jmax-1) = temp.at(i).at(jmax-2) + dy * heat_flux / kappa;
    }


    /*
    // Rayleigh-Benard Convection
    for(int j = 1; j < jmax; j++){
        //left wall
        temp.at(0).at(j) = temp.at(1).at(j) + dx * heat_flux / kappa;
        //right wall
        temp.at(imax-1).at(j) = temp.at(imax-2).at(j) + dx * heat_flux /kappa;
    }
     */



    grid.set_velocity(u_velocity, velocity_type::U);
    grid.set_velocity(v_velocity, velocity_type::V);

    grid.set_pressure(pres);


    grid.set_temperature(temp);
}


void spec_boundary_val(double &u_inflow,
        double &v_inflow,
        double& T_c,
        double& T_h,
        double &kappa,
        double &heat_flux,
        double val_u_inflow,
        double val_v_inflow,
        double val_T_c,
        double val_T_h,
        double val_kappa,
        double val_heat_flux) {

    u_inflow = val_u_inflow;
    v_inflow = val_v_inflow;

    T_c = val_T_c;
    T_h = val_T_h;
    kappa = val_kappa;
    heat_flux = val_heat_flux;

}

void assign_ptr_nbcells(int *imax, int *jmax, Grid &grid){

    //store pointers to neighbours for inner cells
    for(int j = 1; j < *jmax-1; j++) {
        for (int i = 1; i < *imax - 1; i++) {
            grid.cell(i, j)._nbNorth = &grid.cell(i, j + 1);
            grid.cell(i, j)._nbEast = &grid.cell(i + 1, j);
            grid.cell(i, j)._nbWest = &grid.cell(i - 1, j);
            grid.cell(i, j)._nbSouth = &grid.cell(i, j - 1);
        }
    }
    //neighbour edges
    //bottom left
    grid.cell(0,0)._nbNorth = &grid.cell(0,1);
    grid.cell(0,0)._nbEast = &grid.cell(1,0);
    //top right
    grid.cell(*imax-1, *jmax-1)._nbSouth = &grid.cell(*imax-1,*jmax-2);
    grid.cell(*imax-1, *jmax-1)._nbWest = &grid.cell(*imax-2,*jmax-1);
    //top left
    grid.cell(0,*jmax-1)._nbEast = &grid.cell(1, *jmax-1);
    grid.cell(0, *jmax-1)._nbSouth = &grid.cell(0, *jmax-2);
    //bottom right
    grid.cell(*imax-1, 0)._nbNorth = &grid.cell(*imax-1, 1);
    grid.cell(*imax-1, 0)._nbWest = &grid.cell(*imax-2, 0);

    for(int i = 1; i < *imax-1; i++){
        //bottom
        grid.cell(i, 0)._nbNorth = &grid.cell(i, 1);
        grid.cell(i, 0)._nbEast = &grid.cell(i + 1, 0);
        grid.cell(i, 0)._nbWest = &grid.cell(i-1, 0);
        //top
        grid.cell(i, *jmax-1)._nbSouth  = &grid.cell(i,*jmax-2);
        grid.cell(i, *jmax-1)._nbEast = &grid.cell(i+1, *jmax-1);
        grid.cell(i, *jmax-1)._nbWest = &grid.cell(i-1, *jmax-1);
    }
    for(int j = 1; j <  *jmax-1; j++){
        //left
        grid.cell(0, j)._nbNorth = &grid.cell(0, j+1);
        grid.cell(0, j)._nbSouth = &grid.cell(0, j-1);
        grid.cell(0, j)._nbEast = &grid.cell(1, j);
        //right
        grid.cell(*imax-1, j)._nbNorth = &grid.cell(*imax-1, j+1);
        grid.cell(*imax-1, j)._nbSouth = &grid.cell(*imax-1, j-1);
        grid.cell(*imax-1, j)._nbWest = &grid.cell(*imax-2, j);
    }
}