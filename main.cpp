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
        Swf[2][Ny][Nx],
        tau[Ny][Nx];


int main() {

    clock_t tStart = clock();


    double  Lx = 0.5, // [m]
            Ly_ini = 0.5,
            nu = 1.7e-6,
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
            dryRho_ini = 300,
            Dgrain_ini = 1e-3;


    const int tPlot = 60; // Time interval for plotting [s]


    int     wrc[Ny][Nx]; //  0 : initial condition is T, 1: initial condition is theta

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
            gamma[Ny][Nx];



    std::random_device rd;
//    std::mt19937 gen(5489u);
//    std::mt19937 gen{1};
    std::mt19937 gen{rd()};


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

    double a = -7.3e-4 * qIn * 1e3 * 3600.  + 0.63, // parameters for the parametric equation relating tau and gamma
           b = 3.5e-3 * qIn * 1e3 * 3600.   - 3.08,
           c = -4.6e-3 * qIn * 1e3 * 3600.  + 4.1;

    double p = 535.3, // parameters for the parametric equation relating gamma and grain size
           q = 1.43;

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

            std::normal_distribution<double> r{dryRho_ini, dryRho_ini * 0.1};
            std::normal_distribution<double> d{Dgrain_ini, Dgrain_ini * 0.1};


            dryRho[j][i] = r(gen)  ;
            Dgrain[j][i] = d(gen)  ;

//            if (y[j][i] > 0.1){
//                std::normal_distribution<double> o{1.463e-3, 1.463e-3 * 0.1};
//                Dgrain[j][i] = o(gen)  ;
//                std::normal_distribution<double> q{472., 472. * 0.1};
//                dryRho[j][i] = q(gen)  ;
//
//            }

            gamma[j][i] = p * Dgrain[j][i] + q;

            tau[j][i] = a * pow(gamma[j][i], 2.) + b * gamma[j][i] + c;

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



        Eigen::VectorXd delta_theta_p = Richards(dx, dy, qIn, 0);


        bool NotConverged = true; // Condition for convergence of the first and second order solutions of Heun's method

        while (NotConverged) {

            int k = 0;


            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {


                    theta_p[j][i] = theta[j][i] + dt * delta_theta_p(k); // new theta predicted using Euler

                    if (theta_p[j][i] > (0.9 * porosity[j][i]) || isnan(theta_p[j][i])) { // If it goes above saturation, stop model

                        std::cout << "NaN" << std::endl;
                        return 0;
                    }


                    auto outputs = hysteresis(theta[j][i], theta_p[j][i], theta_s[j][i], theta_r[j][i],
                                              0.9 * porosity[j][i],
                                              thetaR[j][i], wrc[j][i], thetaR_dry, i, j);


                    K[1][j][i] = KFun(Ks[j][i], 1, i, j);

                    k += 1;

                }
            }

            //==================== Compute theta_new (forward Heun) =========================

            Eigen::VectorXd delta_theta = Richards(dx, dy, qIn, 1);


            k = 0;

            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {

                    theta_new[j][i] =
                            theta[j][i] + dt * 0.5 * (delta_theta(k) + delta_theta_p(k)); // new theta using Heun

                    if (theta_new[j][i] > (0.9 * porosity[j][i]) || isnan(theta_new[j][i])) { // If it goes above saturation, stop model

                        std::cout << "NaN " << theta_new[j][i] << std::endl;
                        return 0;

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

                dt_new = std::min(dt * s * sqrt(tolerance / std::max(relative_error, EPS)), dt_ini);

                NotConverged = false;

            } else { // condition rejected, still in while loop

                dt = dt * std::max(s * sqrt(tolerance / std::max(relative_error, EPS)), r_min);

            }

        } // end while loop

        //==================== Calculate psi_new =========================

        for (int j = 0; j < Ny; ++j) {
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






