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


// The function is used to smooth the gradients of temperature across full domain.
// If relative delta for a particular cell is higher than the threshold value soecified upfront, then the function
// will "smooth" the value at this cell with respect to other adjacent cells.
void smooth_temp(
    int imax,
    int jmax,
    Grid& grid,
    double time)
{

    static matrix<double> T;

    grid.temperature(T);
    double rel = 0.04;
    for (int i = 1; i < grid.imaxb() - 1; i++)
    {
        for (int j = 1; j < grid.jmaxb() - 1; j++)
        {

            if (grid.cell(i, j)._cellType == 4) {
                double rel_diff = std::abs(T[i][j] - T[i + 1][j]) / std::max(T[i][j], T[i + 1][j])
                                + std::abs(T[i][j] - T[i][j + 1]) / std::max(T[i][j], T[i][j + 1])
                                + std::abs(T[i][j] - T[i - 1][j]) / std::max(T[i][j], T[i - 1][j])
                                + std::abs(T[i][j] - T[i][j - 1]) / std::max(T[i][j], T[i][j - 1]);
                if (rel_diff > 4 * rel) {
                    T[i][j] = (T[i][j] + T[i][j + 1] + T[i + 1][j] + T[i - 1][j] + T[i][j - 1]) / 5;
                }
            }
        }
    }
    grid.set_temperature(T);
}
// Compute kinetics of HOMOGENEOUS NON-CATALYST REACTION
// (reaction occurs across the entire VOLUME of a cell) 
void homogeneous_noncatalyst_reaction(
    double& C_A_cell,
    double& C_B_cell,
    double& C_C_cell,
    double& T_cell,
    const double* stoichiometric_coeff,
    const double* homogeneous_reaction_coeff,
    const double& activation_energy_forward,
    const double& activation_energy_reverse,
    const double& reaction_rate_constant_factor,
    const double& reaction_heat_effect_Q,
    const double& reduced_heat_capacity,
    const double& chem_dt,
    const double& dS_homog,
    const int& max_fixed_point_iterations)

{
    static double fwd_homog_acting_conc_deficit;
    static double fwd_homog_acting_conc_excess[4];
    static double fwd_homog_acting_conc[4];
    static double rev_homog_acting_conc[4];
    static double fwd_homog_react_const;
    static double rev_homog_react_const;

    static double fixed_point_acting_conc_new;
    static double fixed_point_delta_conc;
    static double initial_conc;

    static double* fwd_homog_acting_conc_deficit_alias;      // Alias for concentration
    static double fwd_homog_stoich_coeff_deficit;
    static double homog_react_coeff_deficit;

    static double* fwd_homog_acting_conc_excess_alias[4];    // Alias for concentration
    static double fwd_homog_stoich_coeff_excess[4];
    static double homog_react_coeff_excess[4];

    static double reaction_final_conc;

    // It makes sense to compute the reaction if there are initial components available
    // Otherwise skip the forward reaction computation
    if (C_A_cell > 0 && C_B_cell > 0)
    {
        // Determine which component is taken in deficiency, defining aliases
        // The program can be extended for extra components
        if (C_A_cell / stoichiometric_coeff[ID::A] <= C_B_cell / stoichiometric_coeff[ID::B])
        {
            initial_conc = C_A_cell;

            fwd_homog_acting_conc_deficit = C_A_cell;
            fwd_homog_acting_conc_deficit_alias = &C_A_cell;
            fwd_homog_stoich_coeff_deficit = stoichiometric_coeff[ID::A];
            homog_react_coeff_deficit = homogeneous_reaction_coeff[ID::A];

            fwd_homog_acting_conc_excess[0] = C_B_cell;
            fwd_homog_acting_conc_excess_alias[0] = &C_B_cell;
            fwd_homog_stoich_coeff_excess[0] = stoichiometric_coeff[ID::B];
            homog_react_coeff_excess[0] = homogeneous_reaction_coeff[ID::B];
        }
        else
        {
            initial_conc = C_B_cell;

            fwd_homog_acting_conc_deficit = C_B_cell;
            fwd_homog_acting_conc_deficit_alias = &C_B_cell;
            fwd_homog_stoich_coeff_deficit = stoichiometric_coeff[ID::B];
            homog_react_coeff_deficit = homogeneous_reaction_coeff[ID::B];

            fwd_homog_acting_conc_excess[0] = C_A_cell;
            fwd_homog_acting_conc_excess_alias[0] = &C_A_cell;
            fwd_homog_stoich_coeff_excess[0] = stoichiometric_coeff[ID::A];
            homog_react_coeff_excess[0] = homogeneous_reaction_coeff[ID::A];
        }

        // Constant of FORWARD reaction
        fwd_homog_react_const = reaction_rate_constant_factor *
            exp(-activation_energy_forward / (R_gas_const * T_cell));

        static double epsilon_chem;
        static int internal_iterations;
        epsilon_chem = INFINITY;
        internal_iterations = 0;

        // Loop while relative error is larger than 1e-6 AND max iterations not exceeded
        while ((epsilon_chem >= 1e-6) && (internal_iterations < max_fixed_point_iterations))
        {
            fixed_point_delta_conc =
                chem_dt * fwd_homog_stoich_coeff_deficit * fwd_homog_react_const
                * std::pow(fwd_homog_acting_conc_deficit, homog_react_coeff_deficit)
                * std::pow(fwd_homog_acting_conc_excess[0], homog_react_coeff_excess[0]);

            // If the reaction rate is fast enough to completely deplete the component
            // which may result in having negative concentration of the components,
            // then it can be considered the reaction has been completed
            reaction_final_conc = *fwd_homog_acting_conc_deficit_alias - fixed_point_delta_conc;
            if (reaction_final_conc < 0)
            {
                fixed_point_delta_conc = *fwd_homog_acting_conc_deficit_alias;
                fixed_point_acting_conc_new = 0;
            }
            else
                fixed_point_acting_conc_new = reaction_final_conc;

            // Compute relative error
            epsilon_chem =
                std::abs((fixed_point_acting_conc_new - fwd_homog_acting_conc_deficit)
                    / fwd_homog_acting_conc_deficit);

            // Update acting concentration of the component in deficit
            fwd_homog_acting_conc_deficit = fixed_point_acting_conc_new;

            // Update acting concentration of the components in excess:
            fwd_homog_acting_conc_excess[0] =
                (*fwd_homog_acting_conc_excess_alias[0]
                    - fwd_homog_stoich_coeff_excess[0] / fwd_homog_stoich_coeff_deficit * fixed_point_delta_conc);

            // Increment internal iterations
            ++internal_iterations;
        }

        fwd_homog_acting_conc[ID::C] =
            (C_C_cell + stoichiometric_coeff[ID::C] / fwd_homog_stoich_coeff_deficit * fixed_point_delta_conc);

        // Update concentrations of components for FORWARD reaction
        *fwd_homog_acting_conc_deficit_alias = fwd_homog_acting_conc_deficit;
        *fwd_homog_acting_conc_excess_alias[0] = fwd_homog_acting_conc_excess[0];
        C_C_cell = (fwd_homog_acting_conc[ID::C] > 1e-10) ? fwd_homog_acting_conc[ID::C] : 0;

        
        // Update temperature change due to heat release/absorption
        //T_cell +=
        //    dS_homog * (*fwd_homog_acting_conc_deficit_alias - initial_conc) / fwd_homog_stoich_coeff_deficit
        //    * (-reaction_heat_effect_Q) / reduced_heat_capacity;
        
    }


    // ++++++++ TURNED OFF TO ISOLATE FORWARD HOMOGENEOUS REACTION FOR DEBUGGING ++++++++ //
    /*

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
    while ((epsilon_chem >= 1e-6) && (internal_iterations < max_fixed_point_iterations))
    {
        fixed_point_delta_conc =
            chem_dt * stoichiometric_coeff[ID::C] * rev_homog_react_const
            * std::pow(rev_homog_acting_conc[ID::C], homogeneous_reaction_coeff[ID::C]);

        // If the reaction rate is fast enough to completely deplete the component
        // which may result in having negative concentration of the components,
        // then it can be considered the reaction has been completed
        reaction_final_conc = C_C[i][j] - fixed_point_delta_conc;
        if (reaction_final_conc < 0)
        {
            fixed_point_delta_conc = C_C[i][j];
            fixed_point_acting_conc_new = 0;
        }
        else
            fixed_point_acting_conc_new = reaction_final_conc;

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
    C_A[i][j] = rev_homog_acting_conc[ID::A];
    C_B[i][j] = rev_homog_acting_conc[ID::B];
    C_C[i][j] = rev_homog_acting_conc[ID::C];


    // Update temperature change due to heat release/absorption (mind the minus sign)
    T_cell -=
        dS_homog * (C_C_cell - initial_conc) / stoichiometric_coeff[ID::C]
        * (-reaction_heat_effect_Q) / reduced_heat_capacity;
    */
};



// Calculates chemical kinetics
// INTERNAL loop may be required for equation integration with very small time steps
// depending on the rate of the reaction. Might be needed for fast catalytic reactions.
void calculate_chem_kinetics(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        const int* is_product,
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

    static double Q;        // Total heat effect

    // Get concentration of each component from the grid
    grid.concentration(C_A, ID::A);
    grid.concentration(C_B, ID::B);
    grid.concentration(C_C, ID::C);
    grid.concentration(C_D, ID::D);

    // Get temperature and pressure from the grid
    grid.temperature(T);

    static double moles_total;

    double fwd_heter_acting_surf[4];
    double rev_heter_acting_surf[4];
    double fwd_heter_react_const;
    double rev_heter_react_const;
    double fwd_heter_react_rate;
    double rev_heter_react_rate;
    double heter_intencity(0);

    // Define the characteristic dimension of the catalyst surface
    // For our 2D case it is linear:
    // if FLUID cell is located to the NORTH or SOUTH of the CATALYST block, then dS = dx
    // if FLUID cell is located to the WEST of EAST of the CATALYST block, then dS = dy
    static double dS_heter;

    // Define the characteristic dimension (total volume of a cell)
    static double dS_homog = dx * dy;

    // Surface fraction defines the ratio of the catalyst surface with adsorbed molecules
    // of the corresponding component, to the total surface of the catalyst
    double heter_surface_fraction[4];
    static double surface_fraction_denominator;

    // Reduced heat capacity of the gas mixture
    static double reduced_heat_capacity;

    // Internal time steps to simulate chemical reaction
    int chem_time_steps = 100;
    double chem_dt = dt / chem_time_steps;

    // Define max iterations for chemical reaction fixed point solver
    int max_fixed_point_iterations = 100;

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

                    
                    homogeneous_noncatalyst_reaction(
                        C_A[i][j], C_B[i][j], C_C[i][j], T[i][j], stoichiometric_coeff, homogeneous_reaction_coeff,
                        activation_energy_forward, activation_energy_reverse,  reaction_rate_constant_factor,
                        reaction_heat_effect_Q, reduced_heat_capacity, chem_dt, dS_homog, max_fixed_point_iterations);

                    // Check if the values are not negative
                    C_A[i][j] = (C_A[i][j] >= 1e-10) ? C_A[i][j] : 0;
                    C_B[i][j] = (C_B[i][j] >= 1e-10) ? C_B[i][j] : 0;
                    C_C[i][j] = (C_C[i][j] >= 1e-10) ? C_C[i][j] : 0;
                    

                    
                    // ++++++++ TURNED OFF TO ISOLATE HOMOGENEOUS REACTION FOR DEBUGGING ++++++++ //
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
                        heter_intencity = (fwd_heter_react_rate - rev_heter_react_rate) * surface_development_coeff * dS_heter * chem_dt;
                    }

                    // Update concentration of the components
                    C_A[i][j] += -stoichiometric_coeff[ID::A] * heter_intencity;
                    C_B[i][j] += -stoichiometric_coeff[ID::B] * heter_intencity;
                    C_C[i][j] += stoichiometric_coeff[ID::C] * heter_intencity;
                   
                    
                    // Check if the values are not negative
                    C_A[i][j] = (C_A[i][j] >= 1e-10) ? C_A[i][j] : 0;
                    C_B[i][j] = (C_B[i][j] >= 1e-10) ? C_B[i][j] : 0;
                    C_C[i][j] = (C_C[i][j] >= 1e-10) ? C_C[i][j] : 0;



                    
                    // Update temperature change due to heat release/absorption
                    //T[i][j] -= dS_homog * heter_intencity / stoichiometric_coeff[ID::C]
                    //    * (-reaction_heat_effect_Q) / reduced_heat_capacity;
                                    
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