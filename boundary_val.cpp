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
    double& dx,
    double& dy,
    double& beta,
    double& delta_t,
    double& GX,
    double& GY) {
    
    // VELOCITY - Declaration and Initialisation

    static matrix<double> u_velocity;
    static matrix<double> v_velocity;

    grid.velocity(u_velocity, velocity_type::U);
    grid.velocity(v_velocity, velocity_type::V);

    // PRESSURE - Declaration and Initialization
    static matrix<double> pres;
    grid.pressure(pres);

    // TEMPERATURE - Declaration and Initialization
    static matrix<double> temp;
    grid.temperature(temp);

    // CONCENTRATION - Declaration and Initialization
    static matrix<double> conc_A;
    grid.concentration(conc_A, ID::A);

    static matrix<double> conc_B, conc_C, conc_D;
    grid.concentration(conc_B, ID::B);
    grid.concentration(conc_C, ID::C);
    grid.concentration(conc_D, ID::D);

    // ----- Boundary conditions of NOSLIP inner cells ----- //

    for (int i = 1; i < grid.imaxb() - 1; i++) {
        for (int j = 1; j < grid.jmaxb() - 1; j++) {
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
                    conc_A.at(i).at(j) = (conc_A.at(i).at(j + 1) + conc_A.at(i + 1).at(j)) / 2;
                    conc_B.at(i).at(j) = (conc_B.at(i).at(j + 1) + conc_B.at(i + 1).at(j)) / 2;
                    conc_C.at(i).at(j) = (conc_C.at(i).at(j + 1) + conc_C.at(i + 1).at(j)) / 2;
                    conc_D.at(i).at(j) = (conc_D.at(i).at(j + 1) + conc_D.at(i + 1).at(j)) / 2;
                    //F[i][j] = u_velocity.at(i).at(j);
                    //G[i][j] = v_velocity.at(i).at(j);
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
                //B_NW
                else if (grid.cell(i, j)._nbNorth->_cellType > 1 && grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0; 
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j + 1);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i - 1).at(j - 1); 
                    pres.at(i).at(j) = (pres.at(i).at(j + 1) + pres.at(i - 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j + 1) + temp.at(i - 1).at(j)) / 2;
                    conc_A.at(i).at(j) = (conc_A.at(i).at(j + 1) + conc_A.at(i - 1).at(j)) / 2;
                    conc_B.at(i).at(j) = (conc_B.at(i).at(j + 1) + conc_B.at(i - 1).at(j)) / 2;
                    conc_C.at(i).at(j) = (conc_C.at(i).at(j + 1) + conc_C.at(i - 1).at(j)) / 2;
                    conc_D.at(i).at(j) = (conc_D.at(i).at(j + 1) + conc_D.at(i - 1).at(j)) / 2;
                    F[i - 1][j] = u_velocity.at(i - 1).at(j) - beta * delta_t / 2 * (temp.at(i - 1).at(j) + temp.at(i).at(j)) * GX;
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
                //B_SW
                else if (grid.cell(i, j)._nbSouth->_cellType > 1 && grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0; 
                    v_velocity.at(i).at(j - 1) = 0; 
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j - 1); 
                    v_velocity.at(i).at(j) = -v_velocity.at(i - 1).at(j);
                    pres.at(i).at(j) = (pres.at(i).at(j - 1) + pres.at(i - 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j - 1) + temp.at(i - 1).at(j)) / 2;
                    conc_A.at(i).at(j) = (conc_A.at(i).at(j - 1) + conc_A.at(i - 1).at(j)) / 2;
                    conc_B.at(i).at(j) = (conc_B.at(i).at(j - 1) + conc_B.at(i - 1).at(j)) / 2;
                    conc_C.at(i).at(j) = (conc_C.at(i).at(j - 1) + conc_C.at(i - 1).at(j)) / 2;
                    conc_D.at(i).at(j) = (conc_D.at(i).at(j - 1) + conc_D.at(i - 1).at(j)) / 2;
                    F[i - 1][j] = u_velocity.at(i - 1).at(j) - beta * delta_t / 2 * (temp.at(i - 1).at(j) + temp.at(i).at(j)) * GX;
                    G[i][j - 1] = v_velocity.at(i).at(j - 1) - beta * delta_t / 2 * (temp.at(i).at(j - 1) + temp.at(i).at(j)) * GY;
                }
                //B_SE
                else if (grid.cell(i, j)._nbSouth->_cellType > 1 && grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j - 1) = 0; 
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j - 1);
                    v_velocity.at(i).at(j) = -v_velocity.at(i + 1).at(j);
                    pres.at(i).at(j) = (pres.at(i).at(j - 1) + pres.at(i + 1).at(j)) / 2;
                    temp.at(i).at(j) = (temp.at(i).at(j - 1) + temp.at(i + 1).at(j)) / 2;
                    conc_A.at(i).at(j) = (conc_A.at(i).at(j - 1) + conc_A.at(i + 1).at(j)) / 2;
                    conc_B.at(i).at(j) = (conc_B.at(i).at(j - 1) + conc_B.at(i + 1).at(j)) / 2;
                    conc_C.at(i).at(j) = (conc_C.at(i).at(j - 1) + conc_C.at(i + 1).at(j)) / 2;
					conc_D.at(i).at(j) = (conc_D.at(i).at(j - 1) + conc_D.at(i + 1).at(j)) / 2;
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                    G[i][j - 1] = v_velocity.at(i).at(j - 1) - beta * delta_t / 2 * (temp.at(i).at(j - 1) + temp.at(i).at(j)) * GY;
                }
                //B_N
                else if (grid.cell(i, j)._nbNorth->_cellType > 1 ) {
                    v_velocity.at(i).at(j) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j + 1); 
                    u_velocity.at(i).at(j) = - u_velocity.at(i).at(j + 1);
                    pres.at(i).at(j) = pres.at(i).at(j + 1);
                    temp.at(i).at(j) = temp.at(i).at(j + 1);
                    conc_A.at(i).at(j) = conc_A.at(i).at(j + 1);
                    conc_B.at(i).at(j) = conc_B.at(i).at(j + 1);
                    conc_C.at(i).at(j) = conc_C.at(i).at(j + 1);
                    conc_D.at(i).at(j) = conc_D.at(i).at(j + 1);
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
                //B_E
                else if (grid.cell(i, j)._nbEast->_cellType > 1) {
                    u_velocity.at(i).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i + 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i + 1).at(j - 1); 
                    pres.at(i).at(j) = pres.at(i + 1).at(j);
                    temp.at(i).at(j) = temp.at(i + 1).at(j);
                    conc_A.at(i).at(j) = conc_A.at(i + 1).at(j);
                    conc_B.at(i).at(j) = conc_B.at(i + 1).at(j);
                    conc_C.at(i).at(j) = conc_C.at(i + 1).at(j);
                    conc_D.at(i).at(j) = conc_D.at(i + 1).at(j);
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                }
                //B_S
                else if (grid.cell(i, j)._nbSouth->_cellType > 1) {
                    v_velocity.at(i).at(j - 1) = 0;
                    u_velocity.at(i - 1).at(j) = -u_velocity.at(i - 1).at(j - 1); 
                    u_velocity.at(i).at(j) = -u_velocity.at(i).at(j - 1);
                    pres.at(i).at(j) = pres.at(i).at(j - 1);
                    temp.at(i).at(j) = temp.at(i).at(j - 1);
                    conc_A.at(i).at(j) = conc_A.at(i).at(j - 1);
                    conc_B.at(i).at(j) = conc_B.at(i).at(j - 1);
                    conc_C.at(i).at(j) = conc_C.at(i).at(j - 1);
                    conc_D.at(i).at(j) = conc_D.at(i).at(j - 1);
                    G[i][j - 1] = v_velocity.at(i).at(j - 1) - beta * delta_t / 2 * (temp.at(i).at(j - 1) + temp.at(i).at(j)) * GY;
                }
                //B_W
                else if (grid.cell(i, j)._nbWest->_cellType > 1) {
                    u_velocity.at(i - 1).at(j) = 0;
                    v_velocity.at(i).at(j) = -v_velocity.at(i - 1).at(j);
                    v_velocity.at(i).at(j - 1) = -v_velocity.at(i - 1).at(j - 1); 
                    pres.at(i).at(j) = pres.at(i - 1).at(j);
                    temp.at(i).at(j) = temp.at(i - 1).at(j);
                    conc_A.at(i).at(j) = conc_A.at(i - 1).at(j);
                    conc_B.at(i).at(j) = conc_B.at(i - 1).at(j);
                    conc_C.at(i).at(j) = conc_C.at(i - 1).at(j);
                    conc_D.at(i).at(j) = conc_D.at(i - 1).at(j);
                    F[i - 1][j] = u_velocity.at(i - 1).at(j) - beta * delta_t / 2 * (temp.at(i - 1).at(j) + temp.at(i).at(j)) * GX;
                }
                else if(
                    grid.cell(i, j)._nbEast->_cellType == NOSLIP && grid.cell(i, j)._nbWest->_cellType == NOSLIP&&
                    grid.cell(i, j)._nbNorth->_cellType == NOSLIP && grid.cell(i, j)._nbSouth->_cellType == NOSLIP)
                {
                    u_velocity.at(i).at(j) = 0;     // EAST
                    u_velocity.at(i - 1).at(j) = 0; // WEST
                    v_velocity.at(i).at(j) = 0;     // NORTH
                    v_velocity.at(i).at(j - 1) = 0; // SOUTH

                    pres.at(i).at(j)=0;
                    temp.at(i).at(j)=0;
                    conc_A.at(i).at(j) = 0;
                    conc_B.at(i).at(j) = 0;
                    conc_C.at(i).at(j) = 0;
                    conc_D.at(i).at(j) = 0;
                    F[i][j] = u_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i + 1).at(j)) * GX;
                    G[i][j] = v_velocity.at(i).at(j) - beta * delta_t / 2 * (temp.at(i).at(j) + temp.at(i).at(j + 1)) * GY;
                }
            }
        }
    }

    /* ---- BC outer cells ---- */

    // bottom and top
    for(int i = 0; i < grid.imaxb(); i++)
    {
        //noslip bottom: checked and added pres and temp [Adrian]
        if(grid.cell(i,0)._cellType == NOSLIP and grid.cell(i,0)._nbNorth->_cellType == FLUID)
        {
            u_velocity.at(i).at(0) = -u_velocity.at(i).at(1);
            v_velocity.at(i).at(0) = 0;

            //G.at(i).at(0) = v_velocity.at(i).at(0);
            G.at(i).at(0) = v_velocity.at(i).at(0) - beta * delta_t / 2 * (temp.at(i).at(0) + temp.at(i).at(1)) * GY;
            pres.at(i).at(0) = pres.at(i).at(1);
            temp.at(i).at(0) = temp.at(i).at(1);
            conc_A.at(i).at(0) = conc_A.at(i).at(1);
            conc_B.at(i).at(0) = conc_B.at(i).at(1);
            conc_C.at(i).at(0) = conc_C.at(i).at(1);
            conc_D.at(i).at(0) = conc_D.at(i).at(1);
        }
        //noslip top: checked and added pres and temp [Adrian]
        if(grid.cell(i, grid.jmaxb() - 1)._cellType == NOSLIP and grid.cell(i, grid.jmaxb() - 1)._nbSouth->_cellType == FLUID)
        {
            u_velocity.at(i).at(grid.jmaxb() - 1) = -u_velocity.at(i).at(grid.jmaxb() - 2);
            v_velocity.at(i).at(grid.jmaxb() - 1) = 0; // Solid cell
            v_velocity.at(i).at(grid.jmaxb() - 2) = 0; // Fluid cell below

            G.at(i).at(grid.jmaxb() - 1) = v_velocity.at(i).at(grid.jmaxb() - 1); // Solid cell
            

            //G.at(i).at(grid.jmaxb() - 2) = v_velocity.at(i).at(grid.jmaxb() - 2);

            G.at(i).at(grid.jmaxb() - 2) = v_velocity.at(i).at(grid.jmaxb() - 2)
                - beta * delta_t / 2 * (temp.at(i).at(grid.jmaxb() - 2) + temp.at(i).at(grid.jmaxb() - 1)) * GY; // Fluid cell below
            
            pres.at(i).at(grid.jmaxb() - 1) = pres.at(i).at(grid.jmaxb() - 2);
            temp.at(i).at(grid.jmaxb() - 1) = temp.at(i).at(grid.jmaxb() - 2);
            conc_A.at(i).at(grid.jmaxb() - 1) = conc_A.at(i).at(grid.jmaxb() - 2);
            conc_B.at(i).at(grid.jmaxb() - 1) = conc_B.at(i).at(grid.jmaxb() - 2);
            conc_C.at(i).at(grid.jmaxb() - 1) = conc_C.at(i).at(grid.jmaxb() - 2);
            conc_D.at(i).at(grid.jmaxb() - 1) = conc_D.at(i).at(grid.jmaxb() - 2);
        }
    }

    // left and right
    for(int j = 0; j < grid.jmaxb(); j++)
    {
        // noslip left: checked and added pres and temp [Adrian]
        if(grid.cell(0,j)._cellType == NOSLIP and grid.cell(0,j)._nbEast->_cellType == FLUID)
        {
            u_velocity.at(0).at(j) = 0;
            v_velocity.at(0).at(j) = -v_velocity.at(1).at(j);

            // F.at(0).at(j) = u_velocity.at(0).at(j);

            F.at(0).at(j) = u_velocity.at(0).at(j) - beta * delta_t / 2 * (temp.at(0).at(j) + temp.at(1).at(j)) * GX;
            pres.at(0).at(j) = pres.at(1).at(j);
            temp.at(0).at(j) = temp.at(1).at(j);
            conc_A.at(0).at(j) = conc_A.at(1).at(j);
            conc_B.at(0).at(j) = conc_B.at(1).at(j);
            conc_C.at(0).at(j) = conc_C.at(1).at(j);
            conc_D.at(0).at(j) = conc_D.at(1).at(j);
        }

        // noslip right: checked and added pres and temp [Adrian]
        if(grid.cell(grid.imaxb() - 1,j)._cellType == NOSLIP and grid.cell(grid.imaxb() - 1, j)._nbWest->_cellType == FLUID)
        {
            u_velocity.at(grid.imaxb() - 1).at(j) = 0; // Solid cell
            u_velocity.at(grid.imaxb() - 2).at(j) = 0; // Fluid cell
            v_velocity.at(grid.imaxb() - 1).at(j) = -v_velocity.at(grid.imaxb() - 2).at(j);

            F.at(grid.imaxb() - 1).at(j) = u_velocity.at(grid.imaxb() - 1).at(j); // Solid cell (no need to calculate it)

            // F.at(grid.imaxb() - 2).at(j) = u_velocity.at(grid.imaxb() - 2).at(j);

            F.at(grid.imaxb() - 2).at(j) = u_velocity.at(grid.imaxb() - 2).at(j)
                - beta * delta_t / 2 * (temp.at(grid.imaxb() - 2).at(j) + temp.at(grid.imaxb() - 1).at(j)) * GX; // Fluid cell (still no need to calculate it)
            
            pres.at(grid.imaxb() - 1).at(j) = pres.at(grid.imaxb() - 2).at(j);
            temp.at(grid.imaxb() - 1).at(j) = temp.at(grid.imaxb() - 2).at(j);
            conc_A.at(grid.imaxb() - 1).at(j) = conc_A.at(grid.imaxb() - 2).at(j);
            conc_B.at(grid.imaxb() - 1).at(j) = conc_B.at(grid.imaxb() - 2).at(j);
            conc_C.at(grid.imaxb() - 1).at(j) = conc_C.at(grid.imaxb() - 2).at(j);
            conc_D.at(grid.imaxb() - 1).at(j) = conc_D.at(grid.imaxb() - 2).at(j);
        }

        // inflow left
        if(grid.cell(0,j)._cellType == INFLOW and grid.cell(0,j)._nbEast->_cellType == FLUID)
        {
                v_velocity.at(0).at(j) = v_inflow;
                u_velocity.at(0).at(j) = u_inflow;
        }

        //inflow from right
        //...

        //outflow to left
        //...

        //outflow to right
        if(grid.cell(grid.imaxb() - 1, j)._cellType == OUTFLOW and grid.cell(grid.imaxb() - 1, j)._nbWest->_cellType == FLUID)
        {
            v_velocity.at(grid.imaxb() - 2).at(j) = v_velocity.at(grid.imaxb() - 1).at(j); //neumann BC also for v_velocity!
            u_velocity.at(grid.imaxb() - 2).at(j) = u_velocity.at(grid.imaxb() - 1).at(j);
        }
    }

    // ---- Freeslip BC for lid driven cavity ---- //
    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        if (grid.cell(i, grid.jmaxb() - 1)._cellType == FREESLIP)
        {
            v_velocity.at(i).at(grid.jmaxb() - 1) = 0;
            u_velocity.at(i).at(grid.jmaxb() - 1) = 2 - u_velocity.at(i).at(grid.jmaxb() - 2);
        }
    }

    grid.set_velocity(u_velocity, velocity_type::U);
    grid.set_velocity(v_velocity, velocity_type::V);

    grid.set_pressure(pres);
    grid.set_temperature(temp);
    grid.set_concentration(conc_A, ID::A);
    grid.set_concentration(conc_B, ID::B);
    grid.set_concentration(conc_C, ID::C);
    grid.set_concentration(conc_D, ID::D);
}


void spec_boundary_val(
    int imax,
    int jmax,
    Grid& grid,
    double& v_inflow,
    double& u_inflow,
    double& T_h,
    double& T_c,
    double& C_inject_A,
    double& C_inject_B,
    double& C_inject_C,
    double& C_inject_D,
    double& dx,
    double& dy,
    double& kappa,
    double& heat_flux,
    double& beta,
    double& delta_t,
    double& GX,
    double& GY,
    int scenarioSpec,
    double& time,
    double& t_end)

{
    // TEMPERATURE - Declaration and Initialization
    static matrix<double> temp;
    grid.temperature(temp);

    // CONCENTRATION - Declaration and Initialization
    static matrix<double> conc_A;
    static matrix<double> conc_B, conc_C, conc_D;
    grid.concentration(conc_A, ID::A);
    grid.concentration(conc_B, ID::B);
    grid.concentration(conc_C, ID::C);
    grid.concentration(conc_D, ID::D);


    // Plane shear, Step over flow, Karmann Vortex (with chemical injection)
    if (scenarioSpec == 2 || scenarioSpec == 3 || scenarioSpec == 4)
    {
        // ----  Dirichlet BC Concentration ---- //
        // We set boundary values for all boundary cells equal zero
        // (no chemical components present in the feed flow)
        // LEFT Wall
        for (int j = 1; j < grid.jmaxb() - 1; j++)
            conc_A.at(0).at(j) = 2 * 0 - conc_A.at(1).at(j);

        // ---- Neumann BC Concentration ---- //
        // RIGHT Wall
        for (int j = 1; j < grid.jmaxb() - 1; j++)
            conc_A.at(grid.imaxb() - 1).at(j) = conc_A.at(grid.imaxb() - 2).at(j);

        for (int i = 1; i < grid.imaxb() - 1; i++)
        {
            // BOTTOM Wall
            conc_A.at(i).at(0) = conc_A.at(i).at(1);
            
            // TOP Wall
            conc_A.at(i).at(grid.jmaxb() - 1) = conc_A.at(i).at(grid.jmaxb() - 2);  
        }

        // Only these cells are non-zero: C_injection point at left wall (INFLOW)
        // and if current time is <50% of t_end
        if (time < t_end * 0.5)
        {
            for (int j = 12; j < 13; j++)
                conc_A.at(0).at(j) = 2 * C_inject_A - conc_A.at(1).at(j);
        }
    }


    // Natural Convection and Fluid Trap
    if(scenarioSpec == 5 or scenarioSpec == 6)
    {
        // ----  Dirichlet BC Temperature ---- //
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        { 
            // T_h LEFT Wall
            temp.at(0).at(j) = 2 * T_h - temp.at(1).at(j);

            // T_c RIGHT Wall
            temp.at(grid.imaxb() - 1).at(j) = 2 * T_c - temp.at(grid.imaxb() - 2).at(j);    
        }

        // ---- Neumann BC Temperature ---- //
        for (int i = 1; i < grid.imaxb() - 1; i++)
        {
            // BOTTOM Wall
            temp.at(i).at(0) = temp.at(i).at(1) + dy * heat_flux / kappa;

            // TOP Wall
            temp.at(i).at(grid.jmaxb() - 1) = temp.at(i).at(grid.jmaxb() - 2) + dy * heat_flux / kappa;
        }
    }
    
    // Rayleigh-Benard Convection
    else if (scenarioSpec == 7)
    {
        for (int i = 1; i < grid.imaxb() - 1; i++)
        {
            // T_h BOTTOM Wall
            temp.at(i).at(0) = 2 * T_h - temp.at(i).at(1);

            // T_c TOP Wall
            temp.at(i).at(grid.jmaxb() - 1) = 2 * T_c - temp.at(i).at(grid.jmaxb() - 2);
        }

        // Rayleigh-Benard Convection
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        {
            // LEFT Wall
            temp.at(0).at(j) = temp.at(1).at(j) + dx * heat_flux / kappa;

            // RIGHT Wall
            temp.at(grid.imaxb() - 1).at(j) = temp.at(grid.imaxb() - 2).at(j) + dx * heat_flux / kappa;
        }
    }

    grid.set_temperature(temp);
    grid.set_concentration(conc_A, ID::A);
    grid.set_concentration(conc_B, ID::B);
    grid.set_concentration(conc_C, ID::C);
    grid.set_concentration(conc_D, ID::D);
}


void assign_ptr_nbcells(Grid &grid){

    // store pointers to neighbours for inner cells
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
    // BOTTOM-LEFT
    grid.cell(0,0)._nbNorth = &grid.cell(0,1);
    grid.cell(0,0)._nbEast = &grid.cell(1,0);
    // TOP-RIGHT
    grid.cell(grid.imaxb() - 1, grid.jmaxb() - 1)._nbSouth = &grid.cell(grid.imaxb() - 1,grid.jmaxb() - 2);
    grid.cell(grid.imaxb() - 1, grid.jmaxb() - 1)._nbWest = &grid.cell(grid.imaxb() - 2,grid.jmaxb() - 1);
    // TOP-LEFT
    grid.cell(0,grid.jmaxb() - 1)._nbEast = &grid.cell(1, grid.jmaxb() - 1);
    grid.cell(0, grid.jmaxb() - 1)._nbSouth = &grid.cell(0, grid.jmaxb() - 2);
    // BOTTOM-RIGHT
    grid.cell(grid.imaxb() - 1, 0)._nbNorth = &grid.cell(grid.imaxb() - 1, 1);
    grid.cell(grid.imaxb() - 1, 0)._nbWest = &grid.cell(grid.imaxb() - 2, 0);

    for(int i = 1; i < grid.imaxb() - 1; i++){
        // BOTTOM
        grid.cell(i, 0)._nbNorth = &grid.cell(i, 1);
        grid.cell(i, 0)._nbEast = &grid.cell(i + 1, 0);
        grid.cell(i, 0)._nbWest = &grid.cell(i - 1, 0);
        // TOP
        grid.cell(i, grid.jmaxb() - 1)._nbSouth = &grid.cell(i, grid.jmaxb() - 2);
        grid.cell(i, grid.jmaxb() - 1)._nbEast = &grid.cell(i + 1, grid.jmaxb() - 1);
        grid.cell(i, grid.jmaxb() - 1)._nbWest = &grid.cell(i - 1, grid.jmaxb() - 1);
    }
    for(int j = 1; j <  grid.jmaxb() - 1; j++){
        // LEFT
        grid.cell(0, j)._nbNorth = &grid.cell(0, j + 1);
        grid.cell(0, j)._nbSouth = &grid.cell(0, j - 1);
        grid.cell(0, j)._nbEast = &grid.cell(1, j);
        // RIGHT
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