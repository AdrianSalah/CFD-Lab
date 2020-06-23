#include "sor.hpp"
#include <cmath>
#include <iostream>
#include <mpi.h>

void sor(
        double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        Grid& grid,
        matrix<double> &RS,
        matrix<double> &P,
        double *res,
        double *res_temp,
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

    // MAY NOT BE NEEDED. TO REMOVE LATER.
    // Setting pressure for boundary cells inside the spatial domain before SOR algorithm
    // grid.set_pressure_for_internal_boundaries();

    // Getting pressure over the whole domain
    //grid.pressure(P, il, ir, jb, jt);

    // Setting pressure for boundary cells inside the spatial domain after the SOR algorithm
    //boundary_val_sor(grid);


    // Set boundary values for the outmost cells of the domain
    // LEFT Wall
    if (il == 0) {
        for (int j = 1+ (jb == 0); j < (jt - jb + 2)- (jt == grid.jmaxb() - 1); j++) {
            P[1][j] = P[2][j];
        }
    }
    // RIGHT Wall
    if (ir == (grid.imaxb() - 1)) {
        for (int j = 1 + (jb == 0); j < (jt - jb + 2) - (jt == grid.jmaxb() - 1); j++) {
            P[ir - il + 1][j] = P[ir - il][j];
        }
    }
    // BOTTOM Wall
    if (jb == 0) {
        for (int i = 1 + (il==0); i < (ir - il + 2) - (ir == grid.imaxb() - 1); i++) {
            P[i][1] = P[i][2];
        }
    }
    // TOP Wall
    if (jt == (grid.jmaxb() - 1)) {
        for (int i = 1 + (il == 0); i < (ir - il + 2) - (ir == grid.imaxb() - 1); i++) {
            P[i][jt - jb + 1] = P[i][jt - jb];
        }
    }
    
    /* SOR iteration for FLUID-cells only*/
    for (i = 1 + (il == 0); i < (ir - il + 2) - (ir == (grid.imaxb() - 1)); i++) {
        for (j = 1 + (jb == 0); j < (jt - jb + 2) - (jt == (grid.jmaxb() - 1)); j++) {
            //if (grid.cell(i, j)._cellType == FLUID)
            {
                P.at(i).at(j) = (1.0 - omg) * P.at(i).at(j)
                    + coeff * (
                        (P.at(i + 1).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                        + (P.at(i).at(j + 1) + P.at(i).at(j - 1)) / (dy * dy)
                        - RS.at(i - 1).at(j - 1)
                    );
            }
        }
    }

    /* compute the residual for FLUID-cells only*/
        rloc = 0;
        for (i = 1 + (il == 0); i < (ir - il + 2) - (ir == (grid.imaxb() - 1)); i++) {
            for (j = 1 + (jb == 0); j < (jt - jb + 2) - (jt == (grid.jmaxb() - 1)); j++) {
                //if (grid.cell(i, j)._cellType == FLUID)
                {
                    rloc += ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                        + (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) / (dy * dy) - RS.at(i - 1).at(j - 1)) *
                        ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) / (dx * dx)
                            + (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) / (dy * dy) - RS.at(i - 1).at(j - 1));
                }
            }
        }

        // NOT VALID ANYMORE (Compute residual of SOR iteration for FLUID-cells only)
        // NOT VALID ANYMORE (now dividing by the quantity of fluid cells instead of by (imax*jmax)
        
        
        // Residual is divided by the number of cells inside each chunk
        rloc = rloc / ((imax) * (jmax));
        rloc = sqrt(rloc);
        /* set residual */
        *res_temp = rloc;
        MPI_Allreduce(res_temp, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (il == 0) {
            for (int j = 1 + (jb == 0); j < (jt - jb + 2) - (jt == grid.jmaxb() - 1); j++) {
                P[1][j] = P[2][j];
            }
        }
        // RIGHT Wall
        if (ir == (grid.imaxb() - 1)) {
            for (int j = 1 + (jb == 0); j < (jt - jb + 2) - (jt == grid.jmaxb() - 1); j++) {
                P[ir - il + 1][j] = P[ir - il][j];
            }
        }
        // BOTTOM Wall
        if (jb == 0) {
            for (int i = 1 + (il == 0); i < (ir - il + 2) - (ir == grid.imaxb() - 1); i++) {
                P[i][1] = P[i][2];
            }
        }
        // TOP Wall
        if (jt == (grid.jmaxb() - 1)) {
            for (int i = 1 + (il == 0); i < (ir - il + 2) - (ir == grid.imaxb() - 1); i++) {
                P[i][jt - jb + 1] = P[i][jt - jb];
            }
        }
    // Setting pressure for the entire domain after the SOR algorithm
    grid.set_pressure(P, il, ir, jb, jt);

}

