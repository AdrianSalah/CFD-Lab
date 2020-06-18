#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <algorithm>


// Determines the values of F and G
void calculate_fg(
        double Re,
        double beta,
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
    static matrix<double> u;
    static matrix<double> v;
    static matrix<double> T;


    grid.velocity(u, velocity_type::U);
    grid.velocity(v, velocity_type::V);
    grid.temperature(T);

    // ----- Boundary values for F and G ----- //
    // F[0, j] = u[0, j]        j = 1...jmax    LEFT
    // F[imax, j] = u[imax, j]  j = 1...jmax    RIGHT
    //
    // G[i, 0] = v[i, 0]        i = 1...imax    BOTTOM
    // G[i, jmax] = v[i, jmax]  i = 1...imax    TOP


    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        F.at(0).at(j) = u.at(0).at(j);
        F.at(grid.imaxb() - 2).at(j) = u.at(grid.imaxb() - 2).at(j);
    }

    for (int i = 1; i < grid.imaxb() - 1; i++) {
        G.at(i).at(0) = v.at(i).at(0);
        G.at(i).at(grid.jmaxb() - 2) = v.at(i).at(grid.jmaxb() - 2);
    }

    // ----- F function initialisation ----- //

    static double d2_u_dx2;
    static double d2_u_dy2;
    static double d_u2_dx;
    static double d_uv_dy;

    // ------ Discretisation of differential operators of F ----- //
    for (int i = 1; i < grid.imaxb() - 1; i++) {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
            if (grid.cell(i, j)._cellType == FLUID && grid.cell(i, j)._nbEast->_cellType == FLUID)
            {
                //second derivative with respect to x
                d2_u_dx2 = 1 / (dx * dx) * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]);
                //second derivative with respect to y
                d2_u_dy2 = 1 / (dy * dy) * (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]);
                //first derivative of (u*u) with respect to x
                d_u2_dx = 1 / (4 * dx) * (
                        (u[i][j] + u[i + 1][j]) * (u[i][j] + u[i + 1][j]) -
                        (u[i - 1][j] + u[i][j]) * (u[i - 1][j] + u[i][j])
                ) +

                          alpha / (4 * dx) * (
                                  abs(u[i][j] + u[i + 1][j]) * (u[i][j] - u[i + 1][j]) -
                                  abs(u[i - 1][j] + u[i][j]) * (u[i - 1][j] - u[i][j])
                          );
                //first derivative of (u*v) with respect to y
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

                F.at(i).at(j) -= 0.5* beta * dt * (T.at(i).at(j) + T.at(i+1).at(j)) * GX;
            }
    }

    // ----- G function initialisation ----- //

    static double d2_v_dx2;
    static double d2_v_dy2;
    static double d_v2_dy;
    static double d_uv_dx;

    // ------ Discretisation of differential operators of G ----- //

    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
            if (grid.cell(i, j)._cellType == FLUID && grid.cell(i, j)._nbNorth->_cellType == FLUID)
            {
            //second derivative of v with respect to x
            d2_v_dx2 = 1 / (dx * dx) * (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]);
            //second derivative of v with respect to y
            d2_v_dy2 = 1 / (dy * dy) * (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]);
            //first derivative of (v*v) with respect to y
            d_v2_dy = 1 / (4 * dy) * (
                    (v[i][j] + v[i][j + 1]) * (v[i][j] + v[i][j + 1]) -
                    (v[i][j - 1] + v[i][j]) * (v[i][j - 1] + v[i][j])
            ) +

                      alpha / (4 * dy) * (
                              abs(v[i][j] + v[i][j + 1]) * (v[i][j] - v[i][j + 1]) -
                              abs(v[i][j - 1] + v[i][j]) * (v[i][j - 1] - v[i][j])
                      );
            //first derivative of (u*v) with respect to x
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

            G.at(i).at(j) -= 0.5 * beta * dt  * (T.at(i).at(j) + T.at(i).at(j+1)) * GY;

            }
    }
}


// Calculates temperature
void calculate_temp(
    double Re,
    double Pr,
    double alpha,
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    Grid& grid)
{
    static matrix<double> u;
    static matrix<double> v;
    static matrix<double> T_old;
    static matrix<double> T_new;

    grid.velocity(u, velocity_type::U);
    grid.velocity(v, velocity_type::V);
    grid.temperature(T_old);
    grid.temperature(T_new);

    static double duT_dx;
    static double dvT_dy;
    static double d2_T_dx2;
    static double d2_T_dy2;

    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        {
            if (grid.cell(i, j)._cellType > 1) {
                duT_dx = (u[i][j] * (T_old[i][j] + T_old[i + 1][j]) - u[i - 1][j] * (T_old[i - 1][j] + T_old[i][j]))
                    + (alpha * std::abs(u[i][j]) * (T_old[i][j] - T_old[i + 1][j])) - (alpha * std::abs(u[i - 1][j]) * (T_old[i - 1][j] - T_old[i][j]));
                duT_dx *= 0.5 / dx;

                dvT_dy = (v[i][j] * (T_old[i][j] + T_old[i][j + 1]) - v[i][j - 1] * (T_old[i][j - 1] + T_old[i][j]))
                    + (alpha * std::abs(v[i][j]) * (T_old[i][j] - T_old[i][j + 1])) - (alpha * std::abs(v[i][j - 1]) * (T_old[i][j - 1] - T_old[i][j]));
                dvT_dy *= 0.5 / dy;

                d2_T_dx2 = (T_old[i + 1][j] - 2 * T_old[i][j] + T_old[i - 1][j]) / (dx * dx);

                d2_T_dy2 = (T_old[i][j + 1] - 2 * T_old[i][j] + T_old[i][j - 1]) / (dy * dy);

                // Explicit euler to solve the temperature equation
                // TODO: add source/sink component Q_temp

                T_new[i][j] = T_old[i][j] + dt * ((1 / (Pr * Re)) * (d2_T_dx2 + d2_T_dy2) - duT_dx - dvT_dy);
            }
        }
    }
    grid.set_temperature(T_new);
}


// Calculates concentration
void calculate_concentration(
    double Re,
    double Pr_diffusion,
    double alpha,
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    Grid& grid)
{
    static matrix<double> u;
    static matrix<double> v;
    static matrix<double> C_old;
    static matrix<double> C_new;

    grid.velocity(u, velocity_type::U);
    grid.velocity(v, velocity_type::V);
    grid.concentration(C_old, ID::A);
    grid.concentration(C_new, ID::A);

    static double duC_dx;
    static double dvC_dy;
    static double d2_C_dx2;
    static double d2_C_dy2;

    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        {
            if (grid.cell(i, j)._cellType > 1) {
                duC_dx = (u[i][j] * (C_old[i][j] + C_old[i + 1][j]) - u[i - 1][j] * (C_old[i - 1][j] + C_old[i][j]))
                    + (alpha * std::abs(u[i][j]) * (C_old[i][j] - C_old[i + 1][j])) - (alpha * std::abs(u[i - 1][j]) * (C_old[i - 1][j] - C_old[i][j]));
                duC_dx *= 0.5 / dx;

                dvC_dy = (v[i][j] * (C_old[i][j] + C_old[i][j + 1]) - v[i][j - 1] * (C_old[i][j - 1] + C_old[i][j]))
                    + (alpha * std::abs(v[i][j]) * (C_old[i][j] - C_old[i][j + 1])) - (alpha * std::abs(v[i][j - 1]) * (C_old[i][j - 1] - C_old[i][j]));
                dvC_dy *= 0.5 / dy;

                d2_C_dx2 = (C_old[i + 1][j] - 2 * C_old[i][j] + C_old[i - 1][j]) / (dx * dx);

                d2_C_dy2 = (C_old[i][j + 1] - 2 * C_old[i][j] + C_old[i][j - 1]) / (dy * dy);

                // Explicit euler to solve the concentration equation
                // TODO: add source/sink component Q_concentration

                C_new[i][j] = C_old[i][j] + dt * ((1 / (Pr_diffusion * Re)) * (d2_C_dx2 + d2_C_dy2) - duC_dx - dvC_dy);
            }
        }
    }
    grid.set_concentration(C_new, ID::A);
}


// Calculatesright hand side of the pressure Poisson equation.
void calculate_rs(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &F,
        matrix<double> &G,
        matrix<double> &RS,
        Grid& grid)
{
    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        {
            RS.at(i).at(j) = 1 / dt * (
                    (F.at(i).at(j) - F.at(i - 1).at(j)) / dx +
                    (G.at(i).at(j) - G.at(i).at(j - 1)) / dy
            );
        }
    }
}


// Function helps to compare the magnitude of two values
static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}


// Calculates maximum absolute velocity across entire domain
double max_abs_velocity(
    int imax,
    int jmax,
    Grid& grid,
    velocity_type type)
{
    static matrix<double> current_velocity; //matrix of current velocity U or V on grid
    grid.velocity(current_velocity, type); //assigns velocity U or V to current_velocity

    // Vector of maximum velocity values in every row (including boundaries, i.e. imaxb):
    static std::vector<double> max_abs_value_per_row(grid.imaxb(), 0);

    // Resetting the values to zeros
    std::fill(max_abs_value_per_row.begin(), max_abs_value_per_row.end(), 0);


    for (int i = 0; i < grid.imaxb(); ++i) {
        max_abs_value_per_row.at(i) = *std::max_element(current_velocity.at(i).begin(), current_velocity.at(i).end(),
                                                        abs_compare);
    }
    //maximum velocity value in grid
    return abs(*std::max_element(max_abs_value_per_row.begin(), max_abs_value_per_row.end(), abs_compare));
}


// Calculates the value of timestep dt, considering stability conditions
void calculate_dt(
    double Re,
    double Pr,
    double Pr_diffusion,
    double tau,
    double* dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    Grid& grid)
{
    // Maximum absolute values for U, V on grid for current time step
    static double max_abs_U;
    max_abs_U = max_abs_velocity(grid.imaxb(), grid.jmaxb(), grid, velocity_type::U);

    static double max_abs_V;
    max_abs_V = max_abs_velocity(grid.imaxb(), grid.jmaxb(), grid, velocity_type::V);

    // Explicit time-stepping stability conditions:
    static double condition_temperature;        // Temperature equation
    static double condition_concentration;      // Concentration equation
    static double condition_viscousity;         // 
    static double condition_CFL;                //
    static double condition_common;             //

    condition_viscousity = 0.5 * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));

    // Pr=nu/alpha so Re*Pr= 1/alpha
    condition_temperature = 0.5 * Pr * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));

    // Pr_diffusion=nu/diffusion_coefficient, so Re*Pr= 1/diffusion_coefficient
    condition_concentration = 0.5 * Pr_diffusion * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));

    // Set CFL stability parameter to high number, if the value of u and v velocities are velow the tolerance
    if (max_abs_V < 1e-06 && max_abs_U < 1e-06) // error tolerance used 1e-06
        condition_CFL = INFINITY;
    else
        condition_CFL = std::min((dx / max_abs_U), (dy / max_abs_V));

    condition_common = std::min(std::min(condition_viscousity, condition_CFL),
        std::min(condition_concentration, condition_temperature));

    *dt = tau * condition_common;
}


// Calculates velocities u and v
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
    static matrix<double> u_velocity;
    static matrix<double> v_velocity;
    static matrix<double> pressure;

    grid.velocity(u_velocity, velocity_type::U);
    grid.velocity(v_velocity, velocity_type::V);
    grid.pressure(pressure);


    for(int i = 1; i < grid.imaxb() - 1; i++)
    {
        for(int j = 1; j < grid.jmaxb() - 1; j++)
            if (grid.cell(i, j)._cellType == FLUID && grid.cell(i, j)._nbEast->_cellType == FLUID)
            {
                u_velocity.at(i).at(j) = F.at(i).at(j) - dt/ dx * (pressure.at(i+1).at(j) - pressure.at(i).at(j));
            }
    }

    for (int i = 1; i < grid.imaxb() - 1; i++){
        for(int j = 1; j < grid.jmaxb() - 1; j++)
            if (grid.cell(i, j)._cellType == FLUID && grid.cell(i, j)._nbNorth->_cellType == FLUID)
            {
                v_velocity.at(i).at(j) = G.at(i).at(j) - dt/ dy * (pressure.at(i).at(j+1) - pressure.at(i).at(j));
            }
    }

    grid.set_velocity(u_velocity, velocity_type::U);
    grid.set_velocity(v_velocity, velocity_type::V);
}


// Initializes F, G and RS matrices
void init_fgrs(int imax,
            int jmax,
            matrix<double> &F,
            matrix<double> &G,
            matrix<double> &RS,
            double FI,
            double GI,
            double RSI,
            Grid& grid)

{
    F.resize(grid.imaxb());
    G.resize(grid.imaxb());
    RS.resize(grid.imaxb());

    for (int i = 0; i < grid.imaxb(); i++)
    {
        F.at(i).resize(grid.jmaxb(), FI);
        G.at(i).resize(grid.jmaxb(), GI);
        RS.at(i).resize(grid.jmaxb(), RSI);
    }
}


// Initializes u, v, p, T and c if the cells belongs to FLUID-cells
void init_uvptc(
    int imax,
    int jmax,
    matrix<double> U,
    matrix<double> V,
    matrix<double> P,
    matrix<double> T,
    matrix<double> C,
    double UI,
    double VI,
    double PI,
    double TI,
    double CI,
    Grid& grid)
{
    grid.velocity(U, velocity_type::U);
    grid.velocity(V, velocity_type::V);
    grid.pressure(P);
    grid.temperature(T);
    grid.concentration(C, ID::A);

    for (int i = 0; i < grid.imaxb(); i++) {
        for (int j = 0; j < grid.jmaxb(); j++) {
            if (grid.cell(i, j)._cellType == FLUID) {
                U.at(i).at(j) = UI;
                V.at(i).at(j) = VI;
                P.at(i).at(j) = PI;
                T.at(i).at(j) = TI;
                C.at(i).at(j) = CI;
            }
            else if (grid.cell(i, j)._cellType == NOSLIP) {
                U.at(i).at(j) = 0;
                V.at(i).at(j) = 0;
                P.at(i).at(j) = 0;
                // T.at(i).at(j) = 0;
                // C.at(i).at(j) = 0;
            }
        }
    }

    grid.set_velocity(U, velocity_type::U);
    grid.set_velocity(V, velocity_type::V);
    grid.set_pressure(P);
    grid.set_temperature(T);
    grid.set_concentration(C, ID::A);
}