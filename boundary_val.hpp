#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

/**
 * The boundary values of the problem are set.
 */
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
                    double &heat_flux);

/*
 * The values of boundary conditions are specified
 */
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
                       double val_heat_flux);

/*
 * Store pointers to neighbour cells for grid
 */

void assign_ptr_nbcells(Grid &grid);

#endif
