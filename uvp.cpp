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
    matrix<double> u;
    matrix<double> v;

    grid.velocity(u, velocity_type::U);
    grid.velocity(v, velocity_type::V);

    // -----Boundary values for F and G----- //
    // F[0, j] = u[0, j]        j = 1...jmax    LEFT
    // F[imax, j] = u[imax, j]  j = 1...jmax    RIGHT
    //
    // G[i, 0] = v[i, 0]        i = 1...imax    BOTTOM
    // G[i, jmax] = v[i, jmax]  i = 1...imax    TOP


    for (int j = 1; j <= jmax; j++) {
        F[0][j] = u[0][j];
        F[imax][j] = u[imax][j];
    }

    for (int i = 1; i <= imax; i++) {
        G[i][0] = v[i][0];
        G[i][jmax] = v[i][jmax];
    }

    // -----F function initialisation----- //

    double d2_u_dx2;
    double d2_u_dy2;
    double d_u2_dx;
    double d_uv_dy;


    for (int i = 1; i < imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            d2_u_dx2 = 1 / (dx * dx) * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]);

            d2_u_dy2 = 1 / (dy * dy) * (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]);

            d_u2_dx = 1 / (4 * dx) * (
                (u[i][j] + u[i + 1][j]) * (u[i][j] + u[i + 1][j]) -
                (u[i - 1][j] + u[i][j]) * (u[i - 1][j] + u[i][j])
                ) +

                alpha / (4 * dx) * (
                    abs(u[i][j] + u[i + 1][j]) * (u[i][j] - u[i + 1][j]) -
                    abs(u[i - 1][j] + u[i][j]) * (u[i - 1][j] - u[i][j])
                    );

            d_uv_dy = 1 / (4 * dy) * (
                (v[i][j] + v[i + 1][j]) * (u[i][j] + u[i][j + 1]) -
                (v[i][j - 1] + v[i + 1][j - 1]) * (u[i][j - 1] + u[i][j])
                ) +
                
                alpha / (4 * dy) * (
                    abs(v[i][j] + v[i + 1][j]) * (u[i][j] - u[i][j + 1]) -
                    abs(v[i][j - 1] + v[i + 1][j - 1]) * (u[i][j - 1] - u[i][j])
                    );

            // To check whether GX should be divided by density
            F.at(i).at(j) = u[i][j] + dt * (1 / Re * (d2_u_dx2 + d2_u_dy2) - d_u2_dx - d_uv_dy + GX);
        }
    }

    // -----G function initialisation----- //

    double d2_v_dx2;
    double d2_v_dy2;
    double d_v2_dy;
    double d_uv_dx;


    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j < jmax; j++)
        {
            d2_v_dx2 = 1 / (dx * dx) * (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]);

            d2_v_dy2 = 1 / (dy * dy) * (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]);

            d_v2_dy = 1 / (4 * dy) * (
                (v[i][j] + v[i][j + 1]) * (v[i][j] + v[i][j + 1]) -
                (v[i][j - 1] + v[i][j]) * (v[i][j - 1] + v[i][j])
                ) +

                alpha / (4 * dy) * (
                    abs(v[i][j] + v[i][j + 1]) * (v[i][j] - v[i][j + 1]) -
                    abs(v[i][j - 1] + v[i][j]) * (v[i][j - 1] - v[i][j])
                    );

            d_uv_dx = 1 / (4 * dx) * (
                (u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) -
                (u[i - 1][j] + u[i - 1][j + 1]) * (v[i - 1][j] + v[i][j])
                ) +

                alpha / (4 * dx) * (
                    abs(u[i][j] + u[i][j + 1]) * (v[i][j] - v[i + 1][j]) -
                    abs(u[i - 1][j] + u[i - 1][j + 1]) * (v[i - 1][j] - v[i][j])
                    );

            // To check whether GY should be divided by density
            G.at(i).at(j) = v[i][j] + dt * (1 / Re * (d2_v_dx2 + d2_v_dy2) - d_uv_dx - d_v2_dy + GY);
        }
    }
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
    for (int i = 1; i < imax; i++)
    {
        for (int j = 1; j < jmax; j++)
        {
            RS.at(i).at(j) = 1 / dt * (
                (F.at(i).at(j) - F.at(i - 1).at(j)) / dx +
                (G.at(i).at(j) - G.at(i).at(j - 1)) / dy
                );
        }
    }
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
