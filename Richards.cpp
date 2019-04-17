#include "Richards.h"
#include "global.h"



//========================= Richards equation solver  ===============================
arma::vec Richards( double tau_0,  double lambda, double dx, double dy, double qIn, int p){

    double tau[Ny][Nx],
            KSouth[Ny][Nx],
            KNorth[Ny][Nx],
            KEast[Ny][Nx],
            KWest[Ny][Nx];


    arma::sp_mat Mat(Nx * Ny, Nx * Ny);

    arma::vec F(Ny * Nx),
            delta_theta(Ny*Nx);


    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

            tau[j][i] = tau_0 * pow(Swf[p][j][i], -lambda);


            if (j > 0) {
                KSouth[j][i] = pow(K[p][j - 1][i] * K[p][j][i], 0.5);
//                KSouth[j][i] = (K[j - 1][i] + K[j][i])  * 0.5;

            } else {
                KSouth[j][i] = K[p][j][i];
            }

            if (j < Ny - 1) {
                KNorth[j][i] = pow(K[p][j + 1][i] * K[p][j][i], 0.5);
//                KNorth[j][i] = (K[j + 1][i] + K[j][i]) * 0.5;

            } else {
                KNorth[j][i] = K[p][j][i];
            }

            if (i > 0) {
                KWest[j][i] = pow(K[p][j][i - 1] * K[p][j][i], 0.5);
//                KWest[j][i] = (K[j][i - 1] + K[j][i]) * 0.5;

            } else {
                KWest[j][i] = 0.;
            }

            if (i < Nx - 1) {
                KEast[j][i] = pow(K[p][j][i + 1] * K[p][j][i], 0.5);
//                KEast[j][i] = (K[j][i + 1] + K[j][i]) * 0.5;

            } else {
                KEast[j][i] = 0.;
            }

        }
    }

    arma::uword k = 0;

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx;++i) {

            if (j == 0) { //======== BOTTOM BOUNDARY ==================================

                Mat(k, k) = 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                              + KNorth[j][i] / pow(dy, 2.));  // (i,j)

                if (i > 0 && i < Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[p][j][i] - psi[p][j + 1][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else if (i == 0) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[p][j][i] - psi[p][j + 1][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else if (i == Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][0] - psi[p][j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[p][j][i] - psi[p][j + 1][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                }

            } else if (j == Ny - 1) { //======== TOP BOUNDARY ==================================

                Mat(k, k) = 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                              + KSouth[j][i] / pow(dy, 2.));  // (i,j)

                if (i > 0 && i < Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy, 2.)
                            + (qIn - KSouth[j][i]) / dy);
                } else if (i == 0) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy, 2.)
                            + (qIn - KSouth[j][i]) / dy);
                } else if (i == Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][0] - psi[p][j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy, 2.)
                            + (qIn - KSouth[j][i]) / dy);
                }


            } else {

                Mat(k, k) = 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                              + (KSouth[j][i] + KNorth[j][i]) / pow(dy, 2.));  // (i,j)

                if (i == 0) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[p][j][i] - psi[p][j + 1][i]) / pow(dy, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else if (i == Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][0] - psi[p][j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[p][j][i] - psi[p][j + 1][i]) / pow(dy, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else {
                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[p][j][i] - psi[p][j + 1][i]) / pow(dy, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                }
            }

            if (i > 0) {
                Mat(k, k - 1) = -tau[j][i - 1] * KWest[j][i] / pow(dx, 2.);  // (i-1,j)
            } else {
                Mat(k, (j + 1) * Nx - 1) = -tau[j][Nx - 1] * KWest[j][i] / pow(dx, 2.); // (i-1,j)
            }

            if (i < Nx - 1) {
                Mat(k, k + 1) = -tau[j][i + 1] * KEast[j][i] / pow(dx, 2.);  // (i + 1, j);
            } else {
                Mat(k, j * Nx) = -tau[j][0] * KEast[j][i] / pow(dx, 2.);  // (i + 1, j)
            }

            if (j > 0) {
                Mat(k, k - Nx) = -tau[j - 1][i] * KSouth[j][i] / pow(dy, 2.); // (i,j-1)
            }

            if (j < Ny - 1) {
                Mat(k, k + Nx) = -tau[j + 1][i] * KNorth[j][i] / pow(dy, 2.); // (i,j+1)
            }


            k += 1;
        }
    }

    delta_theta = arma::spsolve(Mat, F, "superlu");


    return delta_theta;

}


