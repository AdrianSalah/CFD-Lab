#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"
#include "grid.hpp"


//The boundary values of the problem are set. 
void boundaryvalues(
    int imax,
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
    double& GY);


// The values of boundary conditions are specified
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
    double& t_end);

// Store pointers to neighbour cells for grid
void assign_ptr_nbcells(Grid &grid);

void boundary_val_sor(Grid& grid);

#endif
