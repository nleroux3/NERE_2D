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
        Swf[2][Ny][Nx];


int main() {

    clock_t tStart = clock();


    double  Lx = 0.5, // [m]
            Ly_ini = 0.5,
            nu = 1.7e-6,
            lambda = 0,
            tau_0 = 0.1,
            gamma = 2.,
            tf = 3600. ,
            thetaR_dry = 0.02,
            qIn = 10e-3/3600.  ,  // Input flux [m/s]
            epsilon = 1e-6,
            time = 0.,
            tolerance = 0.01, // relative truncation error tolerances
            r_min = 0.1,
            EPS = 1e-10,
            s = 0.6,
            dt_new,
            relative_error,
            dt_ini = 0.1,
            theta_ini = epsilon,
            dt,
            rhoI = 917., // density of ice [kg/m3]
            error = 0,
            dryRho_ini = 300.,
            Dgrain_ini = 1e-3;


    const int tPlot = 60; // Time interval for plotting [s]


    int     wrc[Ny][Nx],
            wrc_new[Ny][Nx]; //  0 : initial condition is T, 1: initial condition is theta

    double  x[Ny][Nx],
            y[Ny][Nx],
            dryRho[Ny][Nx], // Dry snow density [kg m-3]
            Dgrain[Ny][Nx],
            theta[Ny][Nx],
            thetaR[Ny][Nx],
            theta_s[Ny][Nx],
            theta_r[Ny][Nx],
            theta_s_new[Ny][Nx],
            theta_r_new[Ny][Nx],
            theta_p[Ny][Nx],
            theta_new[Ny][Nx], // predicted theta^(n+1) from Euler
            thetaR_new[Ny][Nx];


    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(95., 105.);


    std::ofstream theta_output;

    theta_output.open("theta.dat");


    std::ofstream grain_output;

    grain_output.open("grain.dat");

    std::ofstream density_output;

    density_output.open("density.dat");

    std::ofstream outflow_output;

    outflow_output.open("outflow.dat");

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


            dryRho[j][i] = dryRho_ini  * dist(mt)/100.;
            Dgrain[j][i] = Dgrain_ini  * dist(mt)/100.;

            if (y[j][i] > 0.3){
//                dryRho[j][i] = 417.  * dist(mt)/100.;
                Dgrain[j][i] = 0.7e-3  * dist(mt)/100.;
            }

            porosity[j][i] = 1. - dryRho[j][i] / rhoI;

            theta[j][i] = theta_ini;

            alpha[j][i] = 4.4e6 * pow(dryRho[j][i] / Dgrain[j][i], -0.98);

            alphaW[j][i] = gamma * alpha[j][i];

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



        }
    }


    dt = dt_ini;

    // =================================================================================================
    // ============================ Main loop through time =============================================
    // =================================================================================================



    while (time < tf) {



        //=================================================================================
        //============================== Full step hydrology ==============================
        //=================================================================================


        //==================== Compute theta_p (intermediate value) (forward Euler) =========================


        Eigen::VectorXd delta_theta_p = Richards(tau_0, lambda, dx, dy, qIn, 0);

        again_hydrology:

        int k = 0;

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {


                theta_p[j][i] = theta[j][i] + dt * delta_theta_p(k); // new theta predicted using Euler

                if (theta_p[j][i] > (0.9 * porosity[j][i]) ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta)
                    dt = dt * 0.8;
//                    std::cout<< 111 << std::endl;
                    goto again_hydrology;
                }

                if (theta_p[j][i] <= thetaR[j][i]){
                    if (abs(theta_p[j][i] - thetaR[j][i]) > 1e-10  ) {
//                    std::cout <<  theta_p[j][i] - thetaR[j][i] << " " << 222 << std::endl;

                    dt = dt * 0.8;
                    goto again_hydrology;
                    } else {
                        thetaR[j][i] = theta_p[j][i] - 1e-10;
                    }


                }

                compute_psi:

                auto outputs = hysteresis(theta[j][i], theta_p[j][i], theta_s[j][i], theta_r[j][i], 0.9 * porosity[j][i],
                                          thetaR[j][i], wrc[j][i], thetaR_dry, i, j);

                double thetaR1;

                std::tie(thetaR1, std::ignore, std::ignore, std::ignore) = outputs;


                if (theta_p[j][i] <= thetaR1 ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta
                    dt = dt * 0.8;
//                    std::cout<< 333 << std::endl;

                    goto again_hydrology;

                }

                K[1][j][i]  =  KFun(Ks[j][i], 1, i, j);

                k += 1;

            }
        }

        //==================== Compute theta_new (forward Heun) =========================

        Eigen::VectorXd delta_theta = Richards(tau_0, lambda,dx, dy, qIn, 1);


        k = 0;

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {

                theta_new[j][i] = theta[j][i] + dt * 0.5 * (delta_theta(k) + delta_theta_p(k)) ; // new theta using Heun

                if (theta_new[j][i] > (0.9 * porosity[j][i])  ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta)
                    dt = dt * 0.8;
//                    std::cout<< 555 << std::endl;

                    goto again_hydrology;
                } else if (theta_new[j][i] <= thetaR[j][i]) {
                    if (abs(theta_new[j][i] - thetaR[j][i]) > 1e-10  ) {
//                    std::cout << 222 << std::endl;

                    dt = dt * 0.8;
                    goto again_hydrology;
                    } else {
                        thetaR[j][i] = theta_new[j][i] - 1e-10;
                    }

                }


                //==================== Calculate error  =========================
                error = abs(theta_p[j][i] - theta_new[j][i]) / theta_new[j][i];


                if (k == 0) {
                    relative_error = error;
                } else {
                    relative_error = std::max(relative_error, error);
                }

                k += 1;
            }
        }

        //==================== Compare relative error and tolerance =========================


        if (relative_error <= tolerance) { // time step accepted

            dt_new = std::min(dt * s * sqrt(tolerance / std::max(relative_error,EPS)), dt_ini);

        } else { // condition rejected

            dt = dt * std::max(s * sqrt(tolerance / std::max(relative_error,EPS)), r_min);
//            std::cout<< 666 << std::endl;

            goto again_hydrology;

        }

        //==================== Calculate psi_new =========================

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {


                auto outputs = hysteresis(theta[j][i], theta_new[j][i], theta_s[j][i], theta_r[j][i], 0.9 * porosity[j][i],
                                          thetaR[j][i], wrc[j][i], thetaR_dry, i, j);


                std::tie(thetaR_new[j][i], theta_s_new[j][i], theta_r_new[j][i], wrc_new[j][i]) = outputs;


                if (theta_new[j][i] <= thetaR_new[j][i] ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta)
                    if (abs(theta_new[j][i] - thetaR_new[j][i]) > 1e-10  ) {
//                        std::cout << 222 << std::endl;

                        dt = dt * 0.8;
                        goto again_hydrology;
                    } else {
                        thetaR_new[j][i] = theta_new[j][i] - 1e-10;
                    }

                }

            }
        }



        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {

                theta[j][i] = theta_new[j][i];

                thetaR[j][i] = thetaR_new[j][i];

                psi[0][j][i] = psi[1][j][i];

                Swf[0][j][i] = Swf[1][j][i];

                wrc[j][i] = wrc_new[j][i];

                theta_s[j][i] = theta_s_new[j][i];

                theta_r[j][i] = theta_r_new[j][i];

                K[0][j][i]  =  KFun(Ks[j][i], 0, i, j);
            }
        }


        //=================================================================================
        //================================ Update time  ===================================
        //=================================================================================

        time = time + dt;

        //=================================================================================
        //================ Compute mass balance and write output files  ===================
        //=================================================================================


        if (time > 0 && int(time) % tPlot == 0) {
            if (int(time-dt) != int(time) ) {

                std::cout << int(time / tPlot) << std::endl;

                for (int i = 0; i < Nx; ++i) {
                    for (int j = 0; j < Ny; j++) {

                        theta_output << theta[j][i] << " ";
                        grain_output << Dgrain[j][i] << " ";
                        density_output << dryRho[j][i] << " " ;
                    }

                    outflow_output << K[0][0][i] << " " ;
                }



                density_output << std::endl;

                theta_output << std::endl;
                grain_output << std::endl;
                density_output << std::endl;
                outflow_output << std::endl;
            }

        }



        dt = dt_new;


    } // end of time loop


    theta_output.close();
    grain_output.close();
    density_output.close();

    printf("Time taken: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

}






