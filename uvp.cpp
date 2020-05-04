#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>

//my includes
#include <algorithm>
#include <iostream>

// Determines the value of F and G
void calculate_fg(
        double Re,
        double GX,
        double GY,
        double alpha,
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        matrix<double> &F,
        matrix<double> &G)
{
}

// This operation computes the right hand side of the pressure poisson equation.
void calculate_rs(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &F,
        matrix<double> &G,
        matrix<double> &RS)
{
}

//used in max_abs_velocity function to get index of absolute maximum value
static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}


// Changed the function signature: removed imax and jmax - they are already contained in Grid
// and can be acessed with imax() and jmax() methods. And jmax value is not used here at all.
// Used *pointer to access the abs_max_value directly
// Avoided creating additional vectors in the loop to improve performance and make code more readable
// [Oleg]

double max_abs_velocity(int imax, int jmax, Grid& grid, velocity_type type) {
    matrix<double> current_velocity; //matrix of current velocity U or V on grid
    grid.velocity(current_velocity, type); //assigns velocity U or V to current_velocity
    
    // Vector of maximum velocity values in every row (including boundaries, i.e. imaxb):
    std::vector<double> max_abs_value_per_row(grid.imaxb(), 0);

    for (int i = 0; i < grid.imaxb(); ++i)
        max_abs_value_per_row.at(i) = *std::max_element(current_velocity.at(i).begin(), current_velocity.at(i).end(), abs_compare);

    return *std::max_element(max_abs_value_per_row.begin(), max_abs_value_per_row.end());
}


// Determines the maximal time step size
void calculate_dt(double Re, double tau, double *dt, double dx, double dy, int imax, int jmax, Grid &grid) {

    //maximum absolute values for U, V on grid for current time step
    double max_abs_U = max_abs_velocity(imax, jmax, grid, velocity_type::U);
    double max_abs_V = max_abs_velocity(imax, jmax, grid, velocity_type::V);

    //first stability conditon

    double condition1 = 0.5 * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));  


    //check if CFL stability conditions are too small, then just use first stability condition
    //I am not sure if that's the best way to do it
    
    // Decreased the tolerance and replaced logical OR operator to AND ("&&"), which should be here
    // Because one value may be in "normal" range, even if another is close to zero.
    // Improved code readability and unnecessary initialization of variables (CFL1, CFL2, ec).
    // [Oleg]
    
    if (max_abs_V < 1e-06 && max_abs_U < 1e-06)
        *dt = tau * condition1;
    else
        *dt = tau * std::min(condition1,
                             std::min((dx / max_abs_U), (dy / max_abs_V)));
}





void calculate_uv(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        matrix<double> &F,
        matrix<double> &G)
{
}
