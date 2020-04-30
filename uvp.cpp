#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>

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

// Determines the maximal time step size
void calculate_dt(double Re, double tau, double *dt, double dx, double dy, int imax, int jmax, Grid &grid)
{
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
