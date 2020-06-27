#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <algorithm>


static double R_gas_const = 8.314;    // Gas constant, J/(mol*K)


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
    double* Pr_diffusion,
    double alpha,
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    Grid& grid,
    ID  id)
{
    static matrix<double> u;
    static matrix<double> v;
    static matrix<double> C_old;
    static matrix<double> C_new;

    grid.velocity(u, velocity_type::U);
    grid.velocity(v, velocity_type::V);
    grid.concentration(C_old, id);
    grid.concentration(C_new, id);

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

                C_new[i][j] = C_old[i][j] + dt * ((1 / (Pr_diffusion[id] * Re)) * (d2_C_dx2 + d2_C_dy2) - duC_dx - dvC_dy);
            }
        }
    }
    grid.set_concentration(C_new, id);
}

void calculate_chem_kinetics(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        const bool* is_product,
        const double* stoichiometric_coeff,
        const double* homogeneous_reaction_coeff,
        const double* adsorption_coeff,
        const double* heat_capacity,
        const double& reaction_rate_constant_factor,
        const double& activation_energy_forward,
        const double& activation_energy_reverse,
        const double& activation_energy_catalyst,
        const double& surface_development_coeff,
        const double& vacant_centers_defficiency_coeff,
        const double& reaction_heat_effect_Q)
{
    static matrix<double> C_A;
    static matrix<double> C_B;
    static matrix<double> C_C;
    static matrix<double> C_D;

    static matrix<double> T;
    static matrix<double> P;

    static double Q;        // Total heat effect

    // Defines the sign of the component in the chemical reaction
    // +1 - for product
    // -1 - for initial components
    int signs[4];

    for(auto i = 0; i < 4; i++)
        signs[i] = is_product[i] ? 1 : (-1);

    // Get concentration of each component from the grid
    grid.concentration(C_A, ID::A);
    grid.concentration(C_B, ID::B);
    grid.concentration(C_C, ID::C);
    grid.concentration(C_D, ID::D);

    // Get temperature and pressure from the grid
    grid.temperature(T);
    grid.pressure(P);

    static double moles_total;

    double fwd_homog_acting_conc[4];
    double rev_homog_acting_conc[4];
    double fwd_homog_react_const;
    double rev_homog_react_const;
    double fwd_homog_react_rate;
    double rev_homog_react_rate;
    double homog_intencity(0);

    double fwd_heter_acting_surf[4];
    double rev_heter_acting_surf[4];
    double fwd_heter_react_const;
    double rev_heter_react_const;
    double fwd_heter_react_rate;
    double rev_heter_react_rate;
    double heter_intencity(0);

    double total_react_intencity(0);
    double fixed_point_acting_conc_new;
    double fixed_point_delta_conc;
    double initial_conc;

    // Defines the characteristic dimension of the catalyst surface
    // For our 2D case it is linear:
    // if FLUID cell is located to the NORTH or SOUTH of the CATALYST block, then dS = dx
    // if FLUID cell is located to the WEST of EAST of the CATALYST block, then dS = dy
    static double dS_heter;
    static double dS_homog;

    // Surface fraction defines the ratio of the catalyst surface with adsorbed molecules
    // of the corresponding component, to the total surface of the catalyst
    double heter_surface_fraction[4];
    static double surface_fraction_denominator;

    // Reduced heat capacity of the gas mixture
    static double reduced_heat_capacity;

    // Internal time steps to simulate chemical reaction
    int chem_time_steps = 100;
    double chem_dt = dt / chem_time_steps;

    // Define the accuracy and max iterations for chemical reaction fixed point solver
    double epsilon_chem = INFINITY;
    double max_iterations = 100;
    double internal_iterations = 0;

    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        {
            // Reaction is computed only for FLUID cell
            if (grid.cell(i, j)._cellType > 1)
            {
                for (int internal_time_step = 0; internal_time_step < chem_time_steps; ++internal_time_step)
                {
                    // Compute molar fractions of each component
                    moles_total = C_A[i][j] + C_B[i][j] + C_C[i][j];

                    // Molar heat capacity at constant volume (Cv = R*i/2), where i - degrees of freedom of gas component
                    // One of the model approximations. In general case, a polytropic process should be considered
                    // Compute reduced heat capacity (of the gas mix)
                    reduced_heat_capacity =
                        (heat_capacity[ID::A] * C_A[i][j]
                        + heat_capacity[ID::B] * C_B[i][j]
                        + heat_capacity[ID::C] * C_C[i][j]) / moles_total;

                    // --------------------------------------------- //
                    // ----- HOMOGENEOUS NON-CATALYST REACTION ----- //
                    // --------------------------------------------- //

                    // (reaction occurs across the VOLUME of a cell)

                    // Define the characteristic dimension (total volume of a cell)
                    dS_homog = dx * dy;

                    // ----- IMPLICIT EQUATION FIXED POINT SOLVER ----- //
                    initial_conc = C_A[i][j];

                    fwd_homog_acting_conc[ID::A] = C_A[i][j];
                    fwd_homog_acting_conc[ID::B] = C_B[i][j];

                    // Constant of FORWARD reaction
                    fwd_homog_react_const = reaction_rate_constant_factor *
                        exp(-activation_energy_forward / (R_gas_const * T[i][j]));

                    epsilon_chem = INFINITY;
                    internal_iterations = 0;

                    // Loop while relative error is larger than 1e-6 AND max iterations not exceeded
                    while ((epsilon_chem >= 1e-6) && (internal_iterations < max_iterations))
                    {
                        fixed_point_delta_conc =
                            chem_dt * stoichiometric_coeff[ID::A] * fwd_homog_react_const
                            * std::pow(fwd_homog_acting_conc[ID::A], homogeneous_reaction_coeff[ID::A])
                            * std::pow(fwd_homog_acting_conc[ID::B], homogeneous_reaction_coeff[ID::B]);

                        fixed_point_acting_conc_new = C_A[i][j] - fixed_point_delta_conc;

                        // Compute relative error
                        epsilon_chem =
                            std::abs((fixed_point_acting_conc_new - fwd_homog_acting_conc[ID::A]) / fwd_homog_acting_conc[ID::A]);

                        // Update acting concentration of component A
                        fwd_homog_acting_conc[ID::A] = fixed_point_acting_conc_new;

                        // Update acting concentration of component B:
                        fwd_homog_acting_conc[ID::B] =
                            (C_B[i][j] - stoichiometric_coeff[ID::B] / stoichiometric_coeff[ID::A] * fixed_point_delta_conc);

                        // Increment internal iterations
                        ++internal_iterations;
                    }

                    fwd_homog_acting_conc[ID::C] =
                        (C_C[i][j] + stoichiometric_coeff[ID::C] / stoichiometric_coeff[ID::A] * fixed_point_delta_conc);

                    // Update concentrations of components for FORWARD reaction
                    // Set to zero, if the concentration tends to negative values because of the stability problem
                    C_A[i][j] = (fwd_homog_acting_conc[ID::A] >= 0) ? fwd_homog_acting_conc[ID::A] : 0;
                    C_B[i][j] = (fwd_homog_acting_conc[ID::B] >= 0) ? fwd_homog_acting_conc[ID::B] : 0;
                    C_C[i][j] = (fwd_homog_acting_conc[ID::C] >= 0) ? fwd_homog_acting_conc[ID::C] : 0;

                    // Update temperature change due to heat release/absorption
                    T[i][j] +=
                        dS_homog * (C_A[i][j] - initial_conc) / stoichiometric_coeff[ID::A]
                        * (-reaction_heat_effect_Q) / reduced_heat_capacity;


                    // ----- IMPLICIT EQUATION FIXED POINT SOLVER ----- //

                    // Updates ACTING MOLAR CONCENTRATIONS for REVERSE reaction
                    initial_conc = C_C[i][j];

                    rev_homog_acting_conc[ID::C] = C_C[i][j];

                    // Constant of REVERSE reaction
                    rev_homog_react_const = reaction_rate_constant_factor *
                        exp(-activation_energy_reverse / (R_gas_const * T[i][j]));

                    epsilon_chem = INFINITY;
                    internal_iterations = 0;

                    // Loop while relative error is larger than 1e-6 AND max iterations not exceeded
                    while ((epsilon_chem >= 1e-6) && (internal_iterations < max_iterations))
                    {
                        fixed_point_delta_conc =
                            chem_dt * stoichiometric_coeff[ID::C] * rev_homog_react_const
                            * std::pow(rev_homog_acting_conc[ID::C], homogeneous_reaction_coeff[ID::C]);

                        fixed_point_acting_conc_new = C_C[i][j] - fixed_point_delta_conc;

                        // Compute relative error
                        epsilon_chem =
                            std::abs((fixed_point_acting_conc_new - rev_homog_acting_conc[ID::C]) / rev_homog_acting_conc[ID::C]);

                        // Update acting concentration of component C
                        rev_homog_acting_conc[ID::C] = fixed_point_acting_conc_new;

                        // Increment internal iterations
                        ++internal_iterations;
                    }

                    rev_homog_acting_conc[ID::A] =
                        (C_A[i][j] + stoichiometric_coeff[ID::A] / stoichiometric_coeff[ID::C] * fixed_point_delta_conc);
                    rev_homog_acting_conc[ID::B] =
                        (C_B[i][j] + stoichiometric_coeff[ID::B] / stoichiometric_coeff[ID::C] * fixed_point_delta_conc);

                    // Update concentrations of components for REVERSE reaction
                    // Set to zero, if the concentration tends to negative values because of the stability problem
                    C_A[i][j] = (rev_homog_acting_conc[ID::A] >= 0) ? rev_homog_acting_conc[ID::A] : 0;
                    C_B[i][j] = (rev_homog_acting_conc[ID::B] >= 0) ? rev_homog_acting_conc[ID::B] : 0;
                    C_C[i][j] = (rev_homog_acting_conc[ID::C] >= 0) ? rev_homog_acting_conc[ID::C] : 0;

                    // Update temperature change due to heat release/absorption (mind the minus sign)
                    T[i][j] -=
                        dS_homog * (C_C[i][j] - initial_conc) / stoichiometric_coeff[ID::C]
                        * (-reaction_heat_effect_Q) / reduced_heat_capacity;


                    // NOTE the difference between stoichiometric and homogenous reaction coefficients
                    // Normally, they are different, because the former is the analytical value,
                    // the latter - taken exlusively from experiment
                    rev_homog_react_rate =
                        rev_homog_react_const * std::pow(rev_homog_acting_conc[ID::C], homogeneous_reaction_coeff[ID::C]);

                    // Total effect from FORWARD AND REVERSE reactions, mols/sec
                    homog_intencity = (fwd_homog_react_rate - rev_homog_react_const) * dS_homog;

                    // ------------------------------------------- //
                    // ----- HETEROGENEOUS CATALYST REACTION ----- //
                    // ------------------------------------------- //

                    // Reaction occurs along the SURFACE of the catalyst block
                    // If has at least one adjacent CATALYST block
                    if ((grid.cell(i, j)._nbNorth->_cellType == CellType::CATALYST) ||
                        (grid.cell(i, j)._nbSouth->_cellType == CellType::CATALYST) ||
                        (grid.cell(i, j)._nbWest->_cellType == CellType::CATALYST) ||
                        (grid.cell(i, j)._nbEast->_cellType == CellType::CATALYST))

                    {
                        // Define the characteristic dimension (total area of the reaction surface)
                        dS_heter =
                                (grid.cell(i, j)._nbNorth->_cellType == CellType::CATALYST) * dx
                                + (grid.cell(i, j)._nbSouth->_cellType == CellType::CATALYST) * dx
                                + (grid.cell(i, j)._nbWest->_cellType == CellType::CATALYST) * dy
                                + (grid.cell(i, j)._nbEast->_cellType == CellType::CATALYST) * dy;

                        surface_fraction_denominator = 1
                                                       + adsorption_coeff[ID::A] * C_A[i][j]
                                                       + adsorption_coeff[ID::B] * C_B[i][j]
                                                       + adsorption_coeff[ID::C] * C_C[i][j];

                        // Compute adsorbed surface fraction of each component
                        heter_surface_fraction[ID::A] = adsorption_coeff[ID::A] * C_A[i][j] / surface_fraction_denominator;
                        heter_surface_fraction[ID::B] = adsorption_coeff[ID::B] * C_B[i][j] / surface_fraction_denominator;
                        heter_surface_fraction[ID::C] = adsorption_coeff[ID::C] * C_C[i][j] / surface_fraction_denominator;

                        // Updates surface fraction coefficients for FORWARD reaction
                        // taking into consideration the type of component:
                        // only INITIAL COMPONENTS will be taken into account for computation
                        fwd_heter_acting_surf[ID::A] = heter_surface_fraction[ID::A];
                        fwd_heter_acting_surf[ID::B] = heter_surface_fraction[ID::B];

                        // Constant of the rate of FORWARD CATALYST reaction
                        fwd_heter_react_const = reaction_rate_constant_factor *
                                                exp(-activation_energy_catalyst / (R_gas_const * T[i][j]));

                        // Vacant centers defficiency coefficient can be also considered here... (look up theory if needed)
                        fwd_heter_react_rate = fwd_heter_react_const
                                               * std::pow(fwd_heter_acting_surf[ID::A], stoichiometric_coeff[ID::A])
                                               * std::pow(fwd_heter_acting_surf[ID::B], stoichiometric_coeff[ID::B]);

                        // Updates surface fraction coefficients for REVERSE reaction
                        // taking into consideration the type of component:
                        // Constant of the rate of REVERSE CATALYST reaction
                        rev_heter_acting_surf[ID::C] = heter_surface_fraction[ID::C];

                        rev_heter_react_const = reaction_rate_constant_factor *
                                                exp(-(activation_energy_catalyst - reaction_heat_effect_Q) / (R_gas_const * T[i][j]));

                        rev_heter_react_rate = rev_heter_react_const
                                               * std::pow(rev_heter_acting_surf[ID::C], stoichiometric_coeff[ID::C]);

                        // Total effect from FORWARD AND REVERSE CATALYST reactions, mols/sec
                        heter_intencity = (fwd_heter_react_rate - rev_heter_react_const) * surface_development_coeff * dS_heter * chem_dt;
                    }

                    // Update concentration of the components
                    C_A[i][j] += -stoichiometric_coeff[ID::A] * heter_intencity;
                    C_B[i][j] += -stoichiometric_coeff[ID::B] * heter_intencity;
                    C_C[i][j] += stoichiometric_coeff[ID::C] * heter_intencity;
                    
                    // Check if the values are not negative
                    C_A[i][j] = (C_A[i][j] >= 0) ? C_A[i][j] : 0;
                    C_B[i][j] = (C_B[i][j] >= 0) ? C_B[i][j] : 0;
                    C_C[i][j] = (C_C[i][j] >= 0) ? C_C[i][j] : 0;

                    // Update temperature change due to heat release/absorption
                    T[i][j] -= dS_homog * heter_intencity / stoichiometric_coeff[ID::C]
                        * (-reaction_heat_effect_Q) / reduced_heat_capacity;
                }
            }
        }
    }

    grid.set_concentration(C_A, ID::A);
    grid.set_concentration(C_B, ID::B);
    grid.set_concentration(C_C, ID::C);
    grid.set_concentration(C_D, ID::D);

    grid.set_temperature(T);
}



/* --- Original Implementation of calculate_chem_kinetics kept to prevent bugs in cleaned up version --- */
// Calculates chemical kinetics
// INTERNAL loop may be required for equation integration with very small time steps
// depending on the rate of the reaction. Might be needed for fast catalytic reactions.
/*
void calculate_chem_kinetics(
    double dt,
    double dx,
    double dy,
    int imax,
    int jmax,
    Grid& grid,
    const bool& component_A_is_the_product,
    const bool& component_B_is_the_product,
    const bool& component_C_is_the_product,
    const bool& component_D_is_the_product,
    const double& stoichiometric_coeff_a,
    const double& stoichiometric_coeff_b,
    const double& stoichiometric_coeff_c,
    const double& stoichiometric_coeff_d,
    const double& homogeneous_reaction_coeff_a,
    const double& homogeneous_reaction_coeff_b,
    const double& homogeneous_reaction_coeff_c,
    const double& homogeneous_reaction_coeff_d,
    const double& adsorption_coeff_A,
    const double& adsorption_coeff_B,
    const double& adsorption_coeff_C,
    const double& adsorption_coeff_D,
    const double& heat_capacity_A,
    const double& heat_capacity_B,
    const double& heat_capacity_C,
    const double& heat_capacity_D,
    const double& reaction_rate_constant_factor,
    const double& activation_energy_forward,
    const double& activation_energy_reverse,
    const double& activation_energy_catalyst,
    const double& surface_development_coeff,
    const double& vacant_centers_defficiency_coeff,
    const double& reaction_heat_effect_Q)
{
    static matrix<double> C_A;
    static matrix<double> C_B;
    static matrix<double> C_C;
    static matrix<double> C_D;

    static matrix<double> T;
    static matrix<double> P;

    static double Q;        // Total heat effect

    // Defines the sign of the component in the chemical reaction
    // +1 - for product
    // -1 - for initial components
    int sign_A = (component_A_is_the_product) ? 1 : (-1);
    int sign_B = (component_B_is_the_product) ? 1 : (-1);
    int sign_C = (component_C_is_the_product) ? 1 : (-1);
    int sign_D = (component_D_is_the_product) ? 1 : (-1);

    // Get concentration of each component from the grid
    grid.concentration(C_A, ID::A);
    grid.concentration(C_B, ID::B);
    grid.concentration(C_C, ID::C);
    grid.concentration(C_D, ID::D);

    // Get temperature and pressure from the grid
    grid.temperature(T);
    grid.pressure(P);

    static double moles_total;

    double fwd_homog_acting_conc[4];
    double rev_homog_acting_conc[4];
    double fwd_homog_react_const;
    double rev_homog_react_const;
    double fwd_homog_react_rate;
    double rev_homog_react_rate;
    double homog_intencity(0);

    double fwd_heter_acting_surf[4];
    double rev_heter_acting_surf[4];
    double fwd_heter_react_const;
    double rev_heter_react_const;
    double fwd_heter_react_rate;
    double rev_heter_react_rate;
    double heter_intencity(0);

    double total_react_intencity(0);

    // Defines the characteristic dimension of the catalyst surface
    // For our 2D case it is linear:
    // if FLUID cell is located to the NORTH or SOUTH of the CATALYST block, then dS = dx
    // if FLUID cell is located to the WEST of EAST of the CATALYST block, then dS = dy
    static double dS_heter;
    static double dS_homog;

    // Surface fraction defines the ratio of the catalyst surface with adsorbed molecules
    // of the corresponding component, to the total surface of the catalyst
    double heter_surface_fraction[4];
    static double surface_fraction_denominator;

    // Reduced heat capacity of the gas mixture
    static double reduced_heat_capacity;

    int chem_time_steps = 50;
    double chem_dt = dt / chem_time_steps;

    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        {
            for (int k = 0; k < chem_time_steps; ++k)
            {
                // Reaction is computed only for FLUID cell
                if (grid.cell(i, j)._cellType > 1)
                {
                    //// Compute molar fractions of each component
                    moles_total = C_A[i][j] + C_B[i][j] + C_C[i][j] + C_D[i][j];

                    // Molar heat capacity at constant volume (Cv = R*i/2), where i - degrees of freedom of gas component
                    // One of the model approximations. In general case, a polytropic process should be considered
                    // Compute reduced heat capacity (of the gas mix)
                    reduced_heat_capacity = (heat_capacity_A * C_A[i][j] + heat_capacity_B * C_B[i][j]
                        + heat_capacity_C * C_C[i][j] + heat_capacity_D * C_D[i][j]) / moles_total;

                    // HOMOGENEOUS NON-CATALYST REACTION
                    // Reaction occurs across the VOLUME of a cell

                    // Define the characteristic dimension (total volume of a cell)
                    dS_homog = dx * dy;

                    // Updates ACTING MOLAR CONCENTRATIONS for FORWARD reaction
                    // taking into consideration the type of component:
                    // only INITIAL COMPONENTS will be taken into account for computation
                    // i.e. the products will NOT be included into equation (multiplied by x1)
                    fwd_homog_acting_conc[0] = component_A_is_the_product ? 1 : C_A[i][j];
                    fwd_homog_acting_conc[1] = component_B_is_the_product ? 1 : C_B[i][j];
                    fwd_homog_acting_conc[2] = component_C_is_the_product ? 1 : C_C[i][j];
                    fwd_homog_acting_conc[3] = component_D_is_the_product ? 1 : C_D[i][j];

                    // Constant of FORWARD reaction
                    fwd_homog_react_const = reaction_rate_constant_factor *
                        exp(-activation_energy_forward / (R_gas_const * T[i][j]));

                    fwd_homog_react_rate = fwd_homog_react_const
                        * std::pow(fwd_homog_acting_conc[0], homogeneous_reaction_coeff_a)
                        * std::pow(fwd_homog_acting_conc[1], homogeneous_reaction_coeff_b)
                        * std::pow(fwd_homog_acting_conc[2], homogeneous_reaction_coeff_c)
                        * std::pow(fwd_homog_acting_conc[3], homogeneous_reaction_coeff_d);

                    // Updates ACTING MOLAR CONCENTRATIONS for REVERSE reaction 
                    // taking into consideration the type of component:
                    // only PRODUCTS will be taken into account for computation
                    // i.e. the initial components will NOT be included into equation (multiplied by x1)
                    rev_homog_acting_conc[0] = component_A_is_the_product ? C_A[i][j] : 1;
                    rev_homog_acting_conc[1] = component_B_is_the_product ? C_B[i][j] : 1;
                    rev_homog_acting_conc[2] = component_C_is_the_product ? C_C[i][j] : 1;
                    rev_homog_acting_conc[3] = component_D_is_the_product ? C_D[i][j] : 1;

                    // Constant of the REVERSE reaction
                    rev_homog_react_const = reaction_rate_constant_factor *
                        exp(-activation_energy_reverse / (R_gas_const * T[i][j]));

                    // NOTE the difference between stoichiometric and homogenous reaction coefficients
                    // Normally, they are different, because the former is the analytical value,
                    // the latter - taken exlusively from experiment
                    rev_homog_react_rate = rev_homog_react_const
                        * std::pow(rev_homog_acting_conc[0], homogeneous_reaction_coeff_a)
                        * std::pow(rev_homog_acting_conc[1], homogeneous_reaction_coeff_b)
                        * std::pow(rev_homog_acting_conc[2], homogeneous_reaction_coeff_c)
                        * std::pow(rev_homog_acting_conc[3], homogeneous_reaction_coeff_d);

                    // Total effect from FORWARD AND REVERSE reactions, mols/sec
                    homog_intencity = (fwd_homog_react_rate - rev_homog_react_const) * dS_homog;

                    // HETEROGENEOUS CATALYST REACTION
                    // Reaction occurs along the SURFACE of the catalyst block
                    // If has at least one adjacent CATALYST block
                    if ((grid.cell(i, j)._nbNorth->_cellType == CellType::CATALYST) ||
                        (grid.cell(i, j)._nbSouth->_cellType == CellType::CATALYST) ||
                        (grid.cell(i, j)._nbWest->_cellType == CellType::CATALYST) ||
                        (grid.cell(i, j)._nbEast->_cellType == CellType::CATALYST))

                    {
                        // Define the characteristic dimension (total area of the reaction surface)
                        dS_heter =
                            (grid.cell(i, j)._nbNorth->_cellType == CellType::CATALYST) * dx
                            + (grid.cell(i, j)._nbSouth->_cellType == CellType::CATALYST) * dx
                            + (grid.cell(i, j)._nbWest->_cellType == CellType::CATALYST) * dy
                            + (grid.cell(i, j)._nbEast->_cellType == CellType::CATALYST) * dy;

                        surface_fraction_denominator = 1
                            + adsorption_coeff_A * C_A[i][j]
                            + adsorption_coeff_B * C_B[i][j]
                            + adsorption_coeff_C * C_C[i][j]
                            + adsorption_coeff_D * C_D[i][j];

                        // Compute adsorbed surface fraction of each component
                        heter_surface_fraction[0] = adsorption_coeff_A * C_A[i][j] / surface_fraction_denominator;
                        heter_surface_fraction[1] = adsorption_coeff_B * C_B[i][j] / surface_fraction_denominator;
                        heter_surface_fraction[2] = adsorption_coeff_C * C_C[i][j] / surface_fraction_denominator;
                        heter_surface_fraction[3] = adsorption_coeff_D * C_D[i][j] / surface_fraction_denominator;

                        // Updates surface fraction coefficients for FORWARD reaction 
                        // taking into consideration the type of component:
                        // only INITIAL COMPONENTS will be taken into account for computation
                        // i.e. the products will NOT be included into equation (multiplied by x1)
                        fwd_heter_acting_surf[0] = component_A_is_the_product ? 1 : heter_surface_fraction[0];
                        fwd_heter_acting_surf[1] = component_B_is_the_product ? 1 : heter_surface_fraction[1];
                        fwd_heter_acting_surf[2] = component_C_is_the_product ? 1 : heter_surface_fraction[2];
                        fwd_heter_acting_surf[3] = component_D_is_the_product ? 1 : heter_surface_fraction[3];

                        // Constant of the rate of FORWARD CATALYST reaction
                        fwd_heter_react_const = reaction_rate_constant_factor *
                            exp(-activation_energy_catalyst / (R_gas_const * T[i][j]));

                        // Vacant centers defficiency coefficient can be also considered here... (look up theory if needed)
                        fwd_heter_react_rate = fwd_heter_react_const
                            * std::pow(fwd_heter_acting_surf[0], stoichiometric_coeff_a)
                            * std::pow(fwd_heter_acting_surf[1], stoichiometric_coeff_b)
                            * std::pow(fwd_heter_acting_surf[2], stoichiometric_coeff_c)
                            * std::pow(fwd_heter_acting_surf[3], stoichiometric_coeff_d);

                        // Updates surface fraction coefficients for REVERSE reaction
                        // taking into consideration the type of component:
                        // only PRODUCTS will be taken into account for computation
                        // i.e. the initial components will NOT be included into equation (multiplied by x1)
                        rev_heter_acting_surf[0] = component_A_is_the_product ? heter_surface_fraction[0] : 1;
                        rev_heter_acting_surf[1] = component_B_is_the_product ? heter_surface_fraction[1] : 1;
                        rev_heter_acting_surf[2] = component_C_is_the_product ? heter_surface_fraction[2] : 1;
                        rev_heter_acting_surf[3] = component_D_is_the_product ? heter_surface_fraction[3] : 1;

                        // Constant of the rate of REVERSE CATALYST reaction
                        rev_heter_react_const = reaction_rate_constant_factor *
                            exp(-(activation_energy_catalyst - reaction_heat_effect_Q) / (R_gas_const * T[i][j]));

                        rev_heter_react_rate = rev_heter_react_const
                            * std::pow(rev_heter_acting_surf[0], stoichiometric_coeff_a)
                            * std::pow(rev_heter_acting_surf[1], stoichiometric_coeff_b)
                            * std::pow(rev_heter_acting_surf[2], stoichiometric_coeff_c)
                            * std::pow(rev_heter_acting_surf[3], stoichiometric_coeff_d);

                        // Total effect from FORWARD AND REVERSE CATALYST reactions, mols/sec
                        heter_intencity = (fwd_heter_react_rate - rev_heter_react_const) * surface_development_coeff * dS_heter;
                    }

                    // Total effect from HOMOGENEOUS reaction AND HETEROGENEOUS CATALYST reaction, mols/sec
                    total_react_intencity = (homog_intencity + heter_intencity) * chem_dt;

                    // Update concentration of the components
                    C_A[i][j] += sign_A * stoichiometric_coeff_a * total_react_intencity;
                    C_B[i][j] += sign_B * stoichiometric_coeff_b * total_react_intencity;
                    C_C[i][j] += sign_C * stoichiometric_coeff_c * total_react_intencity;
                    C_D[i][j] += sign_D * stoichiometric_coeff_d * total_react_intencity;

                    // Update temperature change due to heat release/absorption
                    T[i][j] += total_react_intencity * (-reaction_heat_effect_Q) / reduced_heat_capacity;
                }
            }
        }
    }

    grid.set_concentration(C_A, ID::A);
    grid.set_concentration(C_B, ID::B);
    grid.set_concentration(C_C, ID::C);
    grid.set_concentration(C_D, ID::D);

    grid.set_temperature(T);
}
 */


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
    double* Pr_diffusion,
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
    static double Pr_diffusion_max=INFINITY;             //

    for (int i = 0; i < LAST; i++) {
        Pr_diffusion_max = std::min(Pr_diffusion_max, Pr_diffusion[i]);
    }

    condition_viscousity = 0.5 * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));

    // Pr=nu/alpha so Re*Pr= 1/alpha
    condition_temperature = 0.5 * Pr * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));

    // Pr_diffusion=nu/diffusion_coefficient, so Re*Pr= 1/diffusion_coefficient
    condition_concentration = 0.5 * Pr_diffusion_max * Re * (dx * dx) * (dy * dy) / ((dx * dx) + (dy * dy));

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
    double UI,
    double VI,
    double PI,
    double TI,
    double *CI,
    Grid& grid)
{
    grid.velocity(U, velocity_type::U);
    grid.velocity(V, velocity_type::V);
    grid.pressure(P);
    grid.temperature(T);

    matrix<double> C_A, C_B, C_C, C_D;
    grid.concentration(C_A, ID::A);
    grid.concentration(C_B, ID::B);
    grid.concentration(C_C, ID::C);
    grid.concentration(C_D, ID::D);


    for (int i = 0; i < grid.imaxb(); i++) {
        for (int j = 0; j < grid.jmaxb(); j++) {
            if (grid.cell(i, j)._cellType == FLUID) {
                U.at(i).at(j) = UI;
                V.at(i).at(j) = VI;
                P.at(i).at(j) = PI;
                T.at(i).at(j) = TI;
                C_A.at(i).at(j) = CI[ID::A];
                C_B.at(i).at(j) = CI[ID::B];
                C_C.at(i).at(j) = CI[ID::C];
                C_D.at(i).at(j) = CI[ID::D];

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
    grid.set_concentration(C_A, ID::A);
    grid.set_concentration(C_B, ID::B);
    grid.set_concentration(C_C, ID::C);
    grid.set_concentration(C_D, ID::D);
}