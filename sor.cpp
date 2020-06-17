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
    matrix<double>& P
) {
    static int i,j;
    static double rloc;
    static double coeff;
    coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

    /* SOR iteration */
    for(i = 2; i <= imax+1; i++) {
        for(j = 2; j<=jmax+1; j++) {
            P.at(i).at(j) = (1.0-omg)*P.at(i).at(j)
                            + coeff*(( P.at(i+1).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)+P.at(i).at(j-1))/(dy*dy) - RS.at(i-1).at(j-1));
        }
    }

    /* compute the residual */
    rloc = 0;
    for(i = 2; i <= imax+1; i++) {
        for(j = 2; j <= jmax+1; j++) {
            rloc += ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i-1).at(j-1))*
                    ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i-1).at(j-1));
        }
    }
    rloc = rloc/(imax*jmax);
    rloc = sqrt(rloc);
    /* set residual */
    *res = rloc;


    /* set boundary values */
    for(i = 2; i <= imax+1; i++) {
        P.at(i)[1] = P.at(i)[2];
        P.at(i).at(jmax+2) = P.at(i).at(jmax+1);
    }
    for(j = 2; j <= jmax+1; j++) {
        P[1].at(j) = P[2].at(j);
        P.at(imax+2).at(j) = P.at(imax+1).at(j);
    }

}

