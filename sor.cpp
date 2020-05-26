#include "sor.hpp"
#include <cmath>

void sor(
        double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        Grid& grid,
        matrix<double> &RS,
        double *res
) {
    static int i,j;
    static double rloc;
    static double coeff;
    coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    static matrix<double> P;

    // MAY NOT BE NEEDED. TO REMOVE LATER.
    // Setting pressure for boundary cells inside the spatial domain before SOR algorithm
    // grid.set_pressure_for_internal_boundaries();

    // Getting pressure over the whole domain
    grid.pressure(P);


    // Setting pressure for boundary cells inside the spatial domain after the SOR algorithm
    grid.set_pressure_for_internal_boundaries();


    // Set boundary values for the outmost cells of the domain
    for (i = 1; i < grid.imaxb() - 1; i++) {
        // BOTTOM
        if (grid.cell(i, 0)._cellType == NOSLIP)
            P.at(i).at(0) = P.at(i).at(1);

        // TOP
        if (grid.cell(i, grid.jmaxb() - 1)._cellType == NOSLIP)
            P.at(i).at(grid.jmaxb() - 1) = P.at(i).at(grid.jmaxb() - 2);
    }

    for (j = 1; j < grid.jmaxb() - 1; j++) {
        // LEFT
        if (grid.cell(0, j)._cellType == NOSLIP)
            P.at(0).at(j) = P.at(1).at(j);

        // RIGHT
        if (grid.cell(grid.imaxb() - 1, j)._cellType == NOSLIP)
            P.at(grid.imaxb() - 1).at(j) = P.at(grid.imaxb() - 2).at(j);
    }
    
    
    /* SOR iteration for FLUID-cells only*/
    for(i = 1; i < grid.imaxb() - 1; i++) {
        for(j = 1; j < grid.jmaxb() - 1; j++) {
            if (grid.cell(i, j)._cellType == FLUID)
            {
                P.at(i).at(j) = (1.0 - omg) * P.at(i).at(j)
                    + coeff * (
                        (P.at(i + 1).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                        + (P.at(i).at(j + 1) + P.at(i).at(j - 1)) / (dy * dy)
                        - RS.at(i).at(j)
                    );
            }
        }
    }

    /* compute the residual for FLUID-cells only*/
    rloc = 0;
    for(i = 1; i < grid.imaxb() - 1; i++) {
        for(j = 1; j < grid.jmaxb() - 1; j++) {
            if (grid.cell(i, j)._cellType == FLUID)
            {
                rloc += ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                    + (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) / (dy * dy) - RS.at(i).at(j)) *
                    ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                        + (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) / (dy * dy) - RS.at(i).at(j));
            }
        }
    }

    // Compute residual of SOR iteration for FLUID-cells only
    // (now dividing by the quantity of fluid cells instead of by (imax*jmax)
    rloc = rloc / (grid.get_fluid_cells_quantity());
    rloc = sqrt(rloc);
    /* set residual */
    *res = rloc;


    // Setting pressure for the entire domain after the SOR algorithm
    grid.set_pressure(P);

}

