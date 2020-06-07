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
        double *res,
        int il,
        int ir,
        int jb,
        int jt,
        int my_rank
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
    grid.pressure(P, il, ir, jb, jt);


    // Setting pressure for boundary cells inside the spatial domain after the SOR algorithm
    //boundary_val_sor(grid);


    // Set boundary values for the outmost cells of the domain
    if (il == 0) {
        for (int j = 2; j < jt-jb+1; j++) {
            P[1][j] = P[2][j];
        }
    }
    if (ir == grid.imaxb() - 1) {
        for (int j = 2; j < jt - jb + 1; j++) {
            P[ir-il+1][j] = P[ir - il ][j];
        }
    }
    if (jb == 0) {
        for (int i = 2; i < ir-il+1; i++) {
            P[i][1] = P[i][2];
        }
    }
    if (jt == grid.jmaxb() - 1) {
        for (int i = 2; i < ir - il + 1; i++) {
            P[i][jt-jb+1] = P[i][jt-jb];
        }
    }
    
    
    /* SOR iteration for FLUID-cells only*/
    for(i = 1+(il==0); i < ir-il+2 -(ir==(grid.imaxb() - 1)); i++) {
        for(j = 1+(jb==0); j <jt-jb+2 -(jt==(grid.jmaxb() - 1)); j++) {
            //if (grid.cell(i, j)._cellType == FLUID)
            {
                P.at(i).at(j) = (1.0 - omg) * P.at(i).at(j)
                    + coeff * (
                        (P.at(i + 1).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                        + (P.at(i).at(j + 1) + P.at(i).at(j - 1)) / (dy * dy)
                        - RS.at(i-1).at(j-1)
                    );
            }
        }
    }

    /* compute the residual for FLUID-cells only*/
        rloc = 0;
        for (i = 1 + (il == 0); i < ir - il + 2 - (ir == (grid.imaxb() - 1)); i++) {
            for (j = 1 + (jb == 0); j < jt - jb + 2 - (jt == (grid.jmaxb() - 1)); j++) {
                //if (grid.cell(i, j)._cellType == FLUID)
                {
                    rloc += ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                        + (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) / (dy * dy) - RS.at(i-1).at(j-1)) *
                        ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                            + (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) / (dy * dy) - RS.at(i-1).at(j-1));
                }
            }
        }

        // Compute residual of SOR iteration for FLUID-cells only
        // (now dividing by the quantity of fluid cells instead of by (imax*jmax)
        rloc = rloc / ((imax) * (jmax));
        rloc = sqrt(rloc);
        /* set residual */
        *res = rloc;

    // Setting pressure for the entire domain after the SOR algorithm
    grid.set_pressure(P, il, ir, jb, jt);

}

