#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax, int jmax, Grid& grid, double& v_inflow, double& u_inflow, matrix<double>& F,
                    matrix<double>& G, double& TD, double& kappa, double& heat_flux);

/*
 * The values of boundary conditions are specified
 */
void spec_boundary_val(double &u_inflow, double &v_inflow, double& TD, double &kappa, double &heat_flux, double val_u_inflow, double val_v_inflow, double val_TD, double val_kappa, double val_heat_flux);




#endif
