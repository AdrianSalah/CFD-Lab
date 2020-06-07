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
                    double &heat_flux,
                    double &beta,
                    double &delta_t,
                    double& GX,
                    double& GY,
                    int scenarioSpec);

/*
 * The values of boundary conditions are specified
 */
void spec_boundary_val(Grid& grid,
    double& u_inflow,
    double& v_inflow,
    double& T_c,
    double& T_h,
    double& kappa,
    double& heat_flux,
    matrix<double>& U,
    matrix<double>& V,
    matrix<double>& P,
    matrix<double>& T,
    matrix<double>& F,
    matrix<double>& G,
    int il,
    int ir,
    int jb,
    int jt);

/*
 * Store pointers to neighbour cells for grid
 */

void assign_ptr_nbcells(Grid &grid);

void boundary_val_sor(Grid& grid);

#endif
