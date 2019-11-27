#include "Richards.h"
#include "global.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

typedef Eigen::Triplet<double> Tr;
using namespace Eigen;



//========================= Richards equation solver  ===============================
Eigen::VectorXd Richards(const double& dx, const double& dy, const int& p, const int& M, double dy_top[])
{

    double  KSouth[M][Nx],
            KNorth[M][Nx],
            KEast[M][Nx],
            KWest[M][Nx];


    std::vector<Tr> coeffs;

    Eigen::VectorXd F(Ny * Nx);

    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < Nx; ++i) {


            if (j > 0) {
                KSouth[j][i] = sqrt(K[p][j - 1][i] * K[p][j][i]); // Geometric mean
//                KSouth[j][i] = (K[p][j - 1][i] + K[p][j][i])  * 0.5; // Arithmetic mean

            } else {
                KSouth[j][i] = K[p][j][i];
            }

            if (j < M - 1) {
                KNorth[j][i] = sqrt(K[p][j + 1][i] * K[p][j][i]); // Geometric mean
//                KNorth[j][i] = (K[p][j + 1][i] + K[p][j][i]) * 0.5; // Arithmetic mean

            } else {
                KNorth[j][i] = K[p][j][i];
            }

            if (i > 0) {
                KWest[j][i] = sqrt(K[p][j][i - 1] * K[p][j][i]);  // Geometric mean
//                KWest[j][i] = (K[p][j][i - 1] + K[p][j][i]) * 0.5;  // Arithmetic mean

            } else {
                KWest[j][i] = sqrt(K[p][j][Nx - 1] * K[p][j][i]); // Periodic boundary condition, geometric mean
//                KWest[j][i] = (K[p][j][Nx - 1] + K[p][j][i]) * 0.5; // Periodic boundary condition, arithmetic mean

            }

            if (i < Nx - 1) {
                KEast[j][i] = sqrt(K[p][j][i + 1] * K[p][j][i]); // Geometric mean
//                KEast[j][i] = (K[p][j][i + 1] + K[p][j][i]) * 0.5; // Arithmetic mean



            } else {
                KEast[j][i] = sqrt(K[p][j][0] * K[p][j][i]); // Periodic boundary condition, geometric mean
//                KEast[j][i] = (K[p][j][0] + K[p][j][i]) * 0.5; // Periodic boundary condition, arithmetic mean

            }

        }
    }

    int k = 0;

    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < Nx;++i) {

            if (j == 0) { //======== BOTTOM BOUNDARY ==================================

                coeffs.emplace_back(Tr(k,k,1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                                         + KNorth[j][i] / pow(dy, 2.))));
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

                    coeffs.emplace_back(Tr(k, k, 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                             + KSouth[j][i] / pow(dy_top[i], 2.))));

                if (i > 0 && i < Nx - 1) {

                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy_top[i], 2.)
                            + (qIn[i] - KSouth[j][i]) / dy);

                } else if (i == 0) {

                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy_top[i], 2.)
                            + (qIn[i] - KSouth[j][i]) / dy);

                } else if (i == Nx - 1) {

                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][0] - psi[p][j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy_top[i], 2.)
                            + (qIn[i] - KSouth[j][i]) / dy);
                }

            } else {  //======== MIDDLE CELLS (j > 0 and j < Ny-1) ==================================


                coeffs.emplace_back(Tr(k,k, 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                              + (KSouth[j][i] + KNorth[j][i]) / pow(dy, 2.))));
                if (i == 0) {

                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[p][j][i + 1] - psi[p][j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[p][j][i] - psi[p][j + 1][i]) / pow(dy, 2.)
                            + KSouth[j][i] * (psi[p][j - 1][i] - psi[p][j][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);

                } else if (i == Nx - 1) {

                    F(k) = (-KWest[j][i] * (psi[p][j][i] - psi[p][j][i - 1]) / pow(dx, 2.)
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
                coeffs.emplace_back(Tr(k,k-1, -tau[j][i - 1] * KWest[j][i] / pow(dx, 2.)));
            } else {
                coeffs.emplace_back(Tr(k,(j + 1) * Nx - 1, -tau[j][Nx - 1] * KWest[j][i] / pow(dx, 2.)));
            }

            if (i < Nx - 1) {
                coeffs.emplace_back(Tr(k,k+1, -tau[j][i + 1] * KEast[j][i] / pow(dx, 2.)));
            } else {
                coeffs.emplace_back(Tr(k, j*Nx, -tau[j][0] * KEast[j][i] / pow(dx, 2.)));
            }

            if (j > 0) {
                coeffs.emplace_back(Tr(k, k-Nx, -tau[j - 1][i] * KSouth[j][i] / pow(dy, 2.)));
            }

            if (j < Ny - 1) {
                coeffs.emplace_back(Tr(k, k+Nx, -tau[j + 1][i] * KNorth[j][i] / pow(dy, 2.)));
            }


            k += 1;
        }
    }

    Eigen::SparseMatrix<double> Mat(Nx * M, Nx * M);
    Mat.setFromTriplets(coeffs.begin(), coeffs.end());


   Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
   cg.compute(Mat);


   Eigen::VectorXd delta_theta = cg.solve(F);
   return delta_theta;


}


