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
        matrix<double> &G,
        matrix<double> &u,
        matrix<double> &v)
{

    // ----- Boundary values for F and G ----- //
    // F[0, j] = u[0, j]        j = 1...jmax    LEFT
    // F[imax, j] = u[imax, j]  j = 1...jmax    RIGHT
    //
    // G[i, 0] = v[i, 0]        i = 1...imax    BOTTOM
    // G[i, jmax] = v[i, jmax]  i = 1...imax    TOP

    for (int j = 2; j <= jmax+1; j++) {
        F[2][j] = u[2][j];
        F[imax+1][j] = u[imax+1][j];
    }

    for (int i = 2; i <= imax+1; i++) {
        G[i][2] = v[i][2];
        G[i][jmax+2] = v[i][jmax+1];
    }

    // ----- F function initialisation ----- //

    static double d2_u_dx2;
    static double d2_u_dy2;
    static double d_u2_dx;
    static double d_uv_dy;

    // ------ Discretisation of differential operators of G ----- //

    for (int i = 3; i < imax+2; i++) //shouldn't be the index from for (int i = 1; i < imax-1; i++)
    {
        for (int j = 2; j <= jmax+1; j++)
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
                    (v[i-1][j+1] + v[i][j+1]) * (u[i][j] + u[i][j + 1]) -
                    (v[i-1][j] + v[i][j]) * (u[i][j - 1] + u[i][j])
            ) +

                      alpha / (4 * dy) * (
                              abs(v[i-1][j+1] + v[i][j+1]) * (u[i][j] - u[i][j + 1]) -
                              abs(v[i-1][j] + v[i][j]) * (u[i][j - 1] - u[i][j])
                      );

            // To check whether GX should be divided by density
            F.at(i).at(j) = u[i][j] + dt * (1 / Re * (d2_u_dx2 + d2_u_dy2) - d_u2_dx - d_uv_dy + GX);
        }
    }

    // ----- G function initialisation ----- //

    static double d2_v_dx2;
    static double d2_v_dy2;
    static double d_v2_dy;
    static double d_uv_dx;

    // ------ Discretisation of differential operators of F ----- //

    for (int i = 2; i <= imax+1; i++)
    {
        for (int j = 3; j < jmax+2; j++)
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
                    (u[i+1][j-1] + u[i+1][j]) * (v[i][j] + v[i + 1][j]) -
                    (u[i][j-1] + u[i][j]) * (v[i - 1][j] + v[i][j])
            ) +

                      alpha / (4 * dx) * (
                              abs(u[i+1][j-1] + u[i+1][j]) * (v[i][j] - v[i + 1][j]) -
                              abs(u[i ][j-1] + u[i][j]) * (v[i - 1][j] - v[i][j])
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
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            RS.at(i).at(j) = 1 / dt * (
                    (F.at(i+2).at(j+1) - F.at(i+1).at(j+1)) / dx +
                    (G.at(i+1).at(j+2) - G.at(i+1).at(j+1)) / dy
            );
        }
    }
}



static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}



double max_abs_velocity(int imax, int jmax, Grid& grid, velocity_type type) {
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



void calculate_dt(double Re, double tau, double* dt, double dx, double dy, int imax, int jmax, Grid& grid, matrix<double> U, matrix<double> V ) {

    //maximum absolute values for U, V on grid for current time step
    static double max_abs_U=-1.0;
    for (int i = 0; i < imax + 5; i++) {
        for (int j = 0; j < jmax + 4; j++) {
            max_abs_U = std::max(max_abs_U, U.at(i).at(j));
        }
    }

    static double max_abs_V = -1.0;
    for (int i = 0; i < imax + 4; i++) {
        for (int j = 0; j < jmax + 5; j++) {
            max_abs_V = std::max(max_abs_V, V.at(i).at(j));
        }
    }

    //first stability conditon
    static double condition1;
    condition1 = 0.5 * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));


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
        matrix<double> &G,
        matrix<double>& u_velocity,
        matrix<double>& v_velocity,
        matrix<double>& pressure)
{


    for(int i = 2; i <= imax ; i++){
        for(int j = 2; j <= jmax+1; j++){
            u_velocity.at(i+1).at(j) = F.at(i+1).at(j) - dt/ dx * (pressure.at(i+1).at(j) - pressure.at(i).at(j));
        }
    }

    for (int i = 2; i <= imax+1; i++){
        for(int j = 2; j <= jmax ; j++){
            v_velocity.at(i).at(j+1) = G.at(i).at(j+1) - dt/ dy * (pressure.at(i).at(j+1) - pressure.at(i).at(j));
        }
    }

}

void init_fgrs(int imax,
            int jmax,
            matrix<double> &F,
            matrix<double> &G,
            matrix<double> &RS,
            double FI,
            double GI,
            double RSI
  ){
    F.resize(imax + 5);
    G.resize(imax + 4);
    RS.resize(imax + 2);

    for (int i = 0; i < imax + 5; i++) {

        F.at(i).resize(jmax + 4, FI);
    }
    for (int i = 0; i < imax + 4; i++) {

        G.at(i).resize(jmax + 5, GI);
    }
    for (int i = 0; i < imax + 2; i++) {

        RS.at(i).resize(jmax + 2, RSI);
    }
}
void init_uvpt(int imax,
    int jmax,
    matrix<double>& U,
    matrix<double>& V,
    matrix<double>& P,
    double UI,
    double VI,
    double PI
) {
    U.resize(imax + 5);
    V.resize(imax + 4);
    P.resize(imax + 4);

    for (int i = 0; i < imax + 5; i++) {

        U.at(i).resize(jmax + 4, UI);
    }
    for (int i = 0; i < imax + 4; i++) {

        V.at(i).resize(jmax + 5, VI);
    }
    for (int i = 0; i < imax + 4; i++) {

        P.at(i).resize(jmax + 4, PI);
    }
}
void init_uvpd(int imax,
    int jmax,
    matrix<double>& Ud,
    matrix<double>& Vd,
    matrix<double>& Pd,
    double UI,
    double VI,
    double PI
) {
    Ud.resize(imax + 2);
    Vd.resize(imax + 2);
    Pd.resize(imax + 2);

    for (int i = 0; i < imax + 2; i++) {

        Ud.at(i).resize(jmax + 2, UI);
    }
    for (int i = 0; i < imax + 2; i++) {

        Vd.at(i).resize(jmax + 2, VI);
    }
    for (int i = 0; i < imax + 2; i++) {

        Pd.at(i).resize(jmax + 2, PI);
    }
}