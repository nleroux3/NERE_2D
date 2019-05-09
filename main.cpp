#include <iostream>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <random>
#include <tuple>
#include <omp.h>
#include "Richards.h"
#include "hysteresis.h"
#include "KFun.h"
#include "melt.h"
#include "global.h"
#include <Eigen/Sparse>
#include <fstream>

double  alpha[Ny][Nx],
        n[Ny][Nx],
        porosity[Ny][Nx],
        Ks[Ny][Nx],
        psi[2][Ny][Nx],
        K[2][Ny][Nx],
        alphaW[Ny][Nx],
        Swf[2][Ny][Nx],
        tau[Ny][Nx],
        qIn[Nx],
        Ts[Nx],
        T[Ny][Nx],
        Keff[Ny][Nx];


int main() {

    clock_t tStart = clock();


    double  Lx = 0.5, // [m]
            Ly_ini = 0.5,
            nu = 1.7e-6,
            tf = 3600. ,
            thetaR_dry = 0.02,
            qIn_ini = 0. ,  // Input flux [m/s]
            epsilon = 1e-6,
            tolerance = 0.05, // relative truncation error tolerances
            r_min = 0.1,
            Ts_ini = 0., // Initial snow surface temperature assumed constant [C]
            T_ini = -5., // Initial snow internal tempeature assumed constant [C]
            EPS = 1e-10,
            s = 0.6,
            dt_new,
            relative_error,
            dt_ini = 0.1,
            theta_ini = epsilon,
            dt,
            Hn = 100., // Input flux for melt [W/m2]
            rhoI = 917., // density of ice [kg/m3]
            error = 0,
            dryRho_ini = 300,
            Dgrain_ini = 1e-3;


    const int tPlot = 60; // Time interval for plotting [s]


    int     wrc[Ny][Nx], //  0 : initial condition is T, 1: initial condition is theta
            M = Ny; // Number of vertical numerical cells that change through time due to melting

    double  x[Ny][Nx],
            y[Ny][Nx],
            dryRho[Ny][Nx], // Dry snow density [kg m-3]
            Dgrain[Ny][Nx],
            theta[Ny][Nx],
            thetaR[Ny][Nx],
            theta_s[Ny][Nx],
            theta_r[Ny][Nx],
            theta_p[Ny][Nx],
            theta_new[Ny][Nx], // predicted theta^(n+1) from Euler
            gamma[Ny][Nx],
            melt_layer[Nx],
            Ly[Nx];   // Height of the snowpack that changes through time due to melting




    std::random_device rd;
//    std::mt19937 gen(5489u);
//    std::mt19937 gen{1};
    std::mt19937 gen{rd()};

    std::ofstream tau_output;

    tau_output.open("tau.dat");

    std::ofstream theta_output;

    theta_output.open("theta.dat");

    std::ofstream grain_output;

    grain_output.open("grain.dat");

    std::ofstream density_output;

    density_output.open("density.dat");

    std::ofstream outflow_output;

    outflow_output.open("outflow.dat");

    std::ofstream height_output;

    height_output.open("height_snowpack.dat");

    double dx =  Lx / double(Nx);
    double dy =  Ly_ini / double(Ny);


    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

            if (i == 0) {
                x[j][i] = (Lx / double(Nx)) * 0.5;
            } else {
                x[j][i] = (Lx / double(Nx)) + x[j][i - 1];
            }
            if (j == 0) {
                y[j][i] = (Ly_ini / double(Ny)) * 0.5;
            } else {
                y[j][i] = (Ly_ini / double(Ny)) + y[j - 1][i];
            }

            std::normal_distribution<double> r{dryRho_ini, dryRho_ini * 0.05};
            std::normal_distribution<double> d{Dgrain_ini, Dgrain_ini * 0.05};

            dryRho[j][i] = r(gen)  ;
            Dgrain[j][i] = d(gen)  ;

//            dryRho[j][i] = dryRho_ini  ;
//            Dgrain[j][i] = Dgrain_ini  ;

//            if (y[j][i] > 0.3) {
//                std::normal_distribution<double> o{0.5e-3, 0.5e-3 * 0.02};
//                Dgrain[j][i] = o(gen)  ;
//            }

            double permeability = 3. * pow(Dgrain[j][i] * 0.5, 2.) * exp(-0.013 * dryRho[j][i]);

            gamma[j][i] = 2.;

            tau[j][i] = 0.92 - 0.079 * gamma[j][i] - 1.9e-3 * qIn_ini*3600*1e3 + 3e8 * permeability + 4.3e-4 * gamma[j][i] * qIn_ini*3600*1e3 - 1.5e8 * gamma[j][i] * permeability  ;

            tau_output << tau[j][i] << " " ;

            porosity[j][i] = 1. - dryRho[j][i] / rhoI;

            theta[j][i] = theta_ini;

            alpha[j][i] = 4.4e6 * pow(dryRho[j][i] / Dgrain[j][i], -0.98);

            alphaW[j][i] = gamma[j][i] * alpha[j][i];

            n[j][i] = 1. + 2.7e-3 * pow(dryRho[j][i] / Dgrain[j][i], 0.61);

            Ks[j][i] = 9.81 / nu * 3. * pow(Dgrain[j][i] * 0.5, 2.) * exp(-0.013 * dryRho[j][i]);

            thetaR[j][i] = 0.;

            wrc[j][i] = 1;

            theta_s[j][i] = 0.9 * porosity[j][i];

            theta_r[j][i] = thetaR[j][i];

            auto outputs = hysteresis(theta[j][i], theta[j][i], theta_s[j][i], theta_r[j][i], 0.9 * porosity[j][i],
                                       thetaR[j][i], wrc[j][i], thetaR_dry, i, j);

            std::tie(thetaR[j][i], theta_s[j][i], theta_r[j][i], wrc[j][i]) = outputs;


            psi[0][j][i] = psi[1][j][i];
            Swf[0][j][i] = Swf[1][j][i];

            K[0][j][i] = KFun(Ks[j][i], 0, i, j);

            qIn[i] = qIn_ini;
            Ts[i] = Ts_ini;
            T[j][i] = T_ini;
            Ly[i] = Ly_ini;

            Keff[j][i] = 1e-9;

            melt_layer[i] = dy; // depth of the melting layer (used in Richards, changes if Hn > 0)
        }
    }

    tau_output.close();

    // =================================================================================================
    // ============================ Main loop through time =============================================
    // =================================================================================================

    dt = dt_ini;
    double time = 0.;

    while (time <  tf) {


        //=================================================================================
        //============================== Full step hydrology ==============================
        //=================================================================================



        if (Hn > 0.) {

            int melt_flag = 0;

            for (int i = 0; i < Nx; i++) {


                qIn[i] = melt(Hn, y[M - 1][i], Ts[i], T[M - 1][i], Ly[i], Keff[M - 1][i], dryRho[M - 1][i],
                              theta[M - 1][i], dt);

                melt_layer[i] = Ly[i] - (y[M - 2][i] + dy / 2.); // depth of the melting layer. Used to compute flux if melt_flag = 1

                if (Ly[i] <= (y[M - 2][i] + dy / 2. + 0.002e-3) and
                    M > 2) { // the top layer disappears if it gets less than 0.01 mm
                    melt_flag = 1;
                }

            }

            if (melt_flag == 1) { // Merging of layers

                for (int i = 0; i < Nx; i++) {
                    qIn[i] += melt_layer[i] * (dryRho[M - 1][i] / 1000. + theta[M - 1][i]);
                }

                M = M - 1;

            }
        }


        bool NotConverged = true; // Condition for convergence of the first and second order solutions of Heun's method

        while (NotConverged) {

            //==================== Compute theta_p (intermediate value) (forward Euler) =========================
            Eigen::VectorXd delta_theta_p = Richards(dx, dy, 0, M, melt_layer);


            int k = Nx*M - 1;

            for (int j = M-1; j >= 0; --j) {
                for (int i = Nx-1; i >= 0; --i) {

                    theta_p[j][i] = theta[j][i] + dt * delta_theta_p(k); // new theta predicted using Euler

                    if (theta_p[j][i] > (0.9 * porosity[j][i] - epsilon) and j > 0) { // If it goes above saturation, move water to cell below
                        delta_theta_p(k-Nx) += (theta_p[j][i] - (0.9 * porosity[j][i] - epsilon))/dt;
                        theta_p[j][i] = 0.9 * porosity[j][i] - epsilon;
                    }



                    auto outputs = hysteresis(theta[j][i], theta_p[j][i], theta_s[j][i], theta_r[j][i],
                                              0.9 * porosity[j][i],
                                              thetaR[j][i], wrc[j][i], thetaR_dry, i, j);


                    K[1][j][i] = KFun(Ks[j][i], 1, i, j);

                    k -= 1;

                }
            }

            //==================== Compute theta_new (forward Heun) =========================


            Eigen::VectorXd delta_theta = Richards(dx, dy, 1, M, melt_layer);


            k = Nx*M - 1;

            for (int j = M-1; j >= 0; --j) {
                for (int i = Nx-1; i >= 0; --i) {

                    theta_new[j][i] =
                            theta[j][i] + dt * 0.5 * (delta_theta(k) + delta_theta_p(k)); // new theta using Heun

                    if (theta_new[j][i] > (0.9 * porosity[j][i] - epsilon)) { // If it goes above saturation, move water to cell below
                        delta_theta(k-Nx) +=  (theta_new[j][i] - (0.9 * porosity[j][i] - epsilon))/ (0.5 * dt) ;
                        theta_new[j][i] = 0.9 * porosity[j][i] - epsilon;
                    }

                    //==================== Calculate error  =========================
                    error = abs(theta_p[j][i] - theta_new[j][i]) / theta_new[j][i];

                    if (k == 0) {
                        relative_error = error;
                    } else {
                        relative_error = std::max(relative_error, error);
                    }

                    k -= 1;
                }
            }

            //==================== Compare relative error and tolerance =========================
            if (relative_error <= tolerance) { // time step accepted

                dt_new = std::min(dt * s * sqrt(tolerance / std::max(relative_error, EPS)), dt_ini);

                NotConverged = false;

            } else { // condition rejected, still in while loop

                dt = dt * std::max(s * sqrt(tolerance / std::max(relative_error, EPS)), r_min);

            }

        } // end while convergence loop

        //==================== Calculate psi_new =========================

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < Nx; ++i) {

                auto outputs = hysteresis(theta[j][i], theta_new[j][i], theta_s[j][i], theta_r[j][i], 0.9 * porosity[j][i],
                                          thetaR[j][i], wrc[j][i], thetaR_dry, i, j);


                std::tie(thetaR[j][i], theta_s[j][i], theta_r[j][i], wrc[j][i]) = outputs;

                psi[0][j][i] = psi[1][j][i];
                theta[j][i] = theta_new[j][i];
                Swf[0][j][i] = Swf[1][j][i];
                K[0][j][i]  =  KFun(Ks[j][i], 0, i, j);

            }
        }

        //=================================================================================
        //================================ Update time  ===================================
        //=================================================================================

        time = time + dt;

        //=================================================================================
        //==============================  write output files  =============================
        //=================================================================================

        if (time > 0 && int(time) % tPlot == 0) {
            if (int(time-dt) != int(time) ) {

                std::cout << int(time / tPlot) << std::endl;
                height_output << M << " ";
                for (int i = 0; i < Nx; ++i) {
                    for (int j = 0; j < M; j++) {

                        theta_output << theta[j][i] << " ";
                        grain_output << Dgrain[j][i] << " ";
                        density_output << dryRho[j][i] << " " ;
                    }

                    outflow_output << K[0][0][i] << " " ;
                    height_output << Ly[i] << " ";

                }


                theta_output << std::endl;
                grain_output << std::endl;
                density_output << std::endl;
                outflow_output << std::endl;
                height_output << std::endl;

            }

        }

        dt = dt_new;

    } // end of time loop


    theta_output.close();
    grain_output.close();
    density_output.close();

    printf("Time taken: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

}






