#include <iostream>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <math.h>
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
#include <vector>

int Nx, Ny;


std::vector<std::vector<std::vector<double>>> psi,
                                              K,
                                              Swf;

std::vector<std::vector<double>> alpha,
                                 n,
                                 porosity,
                                 Ks,
                                 alphaW,
                                 tau,
                                 T,
                                 Keff;

std::vector<double> qIn,
                    Ts;

int main() {

    clock_t tStart = clock();


    double Lx, // [m]
            Ly_ini,
            nu = 1.7e-6,
            tf ,
            thetaR_dry = 0.02,
            qIn_ini ,  // Input flux [m/s]
            epsilon = 1e-6,
            tolerance = 0.01, // relative truncation error tolerances
            r_min = 0.1,
            Ts_ini = 0., // Initial snow surface temperature assumed constant [C]
            T_ini = 0., // Initial snow internal tempeature assumed constant [C]
            EPS = 1e-10,
            s = 0.6,
            dt_new,
            relative_error,
            dt_ini ,
            theta_ini,
            dt,
            Hn_ini = 0., // Input flux for melt [W/m2]
            rhoI = 917., // density of ice [kg/m3]
            error = 0,
            dryRho_ini,
            Dgrain_ini ,
            Qflux,
            Qbottom = 0., // Heat flux at the bottom of the snowpack [w/m2]
            Hn,
            sigma;

    int tPlot; // Time interval for plotting [s]



    // read inputs from input file

    std::ifstream infile;

    infile.open("../inputs.txt");

    if (!infile) {
        std::cerr << "Unable to open file inputs.txt";
        exit(1);   // call system to stop
    }

    std::string str;

    int line = 0;
    while (std::getline (infile,str)){
        line += 1;
        if (line == 1){
             std::stringstream geek(str);
             geek >> Lx;
        } else if (line == 2){
            std::stringstream geek(str);
            geek >> Ly_ini;
        } else if (line == 3) {
            std::stringstream geek(str);
            geek >> Nx;
        } else if (line == 4) {
            std::stringstream geek(str);
            geek >> Ny;
        } else if (line == 5) {
            std::stringstream geek(str);
            geek >> tPlot;
        } else if (line == 6) {
            std::stringstream geek(str);
            geek >> tf;
        } else if (line == 7) {
            std::stringstream geek(str);
            geek >> dt_ini;
        } else if (line == 8) {
            std::stringstream geek(str);
            geek >> sigma;
        } else if (line == 9) {
            std::stringstream geek(str);
            geek >> qIn_ini;
            qIn_ini = qIn_ini / 3600. ; // [m/hr] to [m/s]
        } else if (line == 10) {
            std::stringstream geek(str);
            geek >> Hn_ini;
        } else if (line == 12) {
            std::stringstream geek(str);
            geek >> dryRho_ini >> theta_ini >> Dgrain_ini;
            theta_ini = std::max(theta_ini,epsilon);
        } else if (line == 13) {
            std::stringstream geek(str);
            geek >> T_ini;
        }
    }

    infile.close();

    alpha.resize(Ny, std::vector<double>(Nx));
    n.resize(Ny, std::vector<double>(Nx));
    porosity.resize(Ny, std::vector<double>(Nx));
    Ks.resize(Ny, std::vector<double>(Nx));
    alphaW.resize(Ny, std::vector<double>(Nx));
    tau.resize(Ny, std::vector<double>(Nx));
    T.resize(Ny, std::vector<double>(Nx));
    Keff.resize(Ny, std::vector<double>(Nx));
    Ts.resize(Ny);
    qIn.resize(Ny);

    psi.resize(2);
    Swf.resize(2);
    K.resize(2);

    for(int i=0; i<2; i++)
    {
        psi[i].resize(Ny);
        K[i].resize(Ny);
        Swf[i].resize(Ny);

        for(int j=0; j < Ny;j++)
        {
            psi[i][j].resize(Nx);
            Swf[i][j].resize(Nx);
            K[i][j].resize(Nx);
        }
    }


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
            Ly[Nx],   // Height of the snowpack that changes through time due to melting
            T_new[Ny][Nx];




    std::random_device rd;
//    std::mt19937 gen(5489u);
//    std::mt19937 gen{1};
    std::mt19937 gen{rd()};

    std::ofstream tau_output;


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

    std::ofstream T_output;

    T_output.open("T.dat");

    double dx =  Lx / double(Nx);
    double dy =  Ly_ini / double(Ny);

    // Initialize the model inputs and parameters

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

            std::normal_distribution<double> r{dryRho_ini, dryRho_ini * sigma};
            std::normal_distribution<double> d{Dgrain_ini, Dgrain_ini * sigma};

            dryRho[j][i] = r(gen)  ;
            Dgrain[j][i] = d(gen)  ;


//            if (y[j][i] > 0.3) {
//                std::normal_distribution<double> o{0.5e-3, 0.5e-3 * 0.02};
//                Dgrain[j][i] = o(gen)  ;
//            }

            gamma[j][i] = 2.;

            tau[j][i] = std::max(0.748 + 0.346 * gamma[j][i] - 1.32e-3 * qIn_ini*3600*1e3 + 92.4 * Dgrain[j][i] -5.02e-4 * gamma[j][i] * qIn_ini*3600*1e3 - 308. * gamma[j][i] * Dgrain[j][i] + 0.976 * qIn_ini*3600*1e3 * Dgrain[j][i], 0.);

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


            qIn[i] = qIn_ini ;

            Ts[i] = Ts_ini;
            T[j][i] = T_ini;
            Ly[i] = Ly_ini;

            double Kw = 0.569;

            Keff[j][i] = (2.5e-6 * pow(dryRho[j][i], 2.) - 1.23e-4 * dryRho[j][i] + 0.024) * (1. - theta[j][i]) +
                    theta[j][i] * Kw;


            melt_layer[i] = dy; // depth of the melting layer (used in Richards, changes if Hn > 0)
        }
    }


    // =================================================================================================
    // ============================ Main loop through time =============================================
    // =================================================================================================

    dt = dt_ini;
    double time = 0.;

    while (time <  tf) {


        //=================================================================================
        //======================== Full step thermodynamics ================================
        //=================================================================================

        Hn = Hn_ini;

        for (int j = 0; j < M; ++j) {
            for (int i = 0 ; i < Nx ; ++i) {

                if (i > 0 and i < Nx-1) {

                    if (j > 0 && j < M - 1) {
                        Qflux = 1. / pow(dy, 2) * (
                                (Keff[j + 1][i] + Keff[j][i])* 0.5 * (T[j + 1][i] - T[j][i]) +
                                (Keff[j - 1][i] + Keff[j][i])* 0.5 * (T[j - 1][i] - T[j][i])) +
                                1. / pow(dx, 2) * (
                                        (Keff[j][i + 1] + Keff[j][i]) * 0.5 * (T[j][i + 1] - T[j][i]) +
                                        (Keff[j][i - 1] + Keff[j][i]) * 0.5 * (T[j][i - 1] - T[j][i]))  ;

                    } else if (j == 0) {

                        Qflux = 1. / pow(dy, 2) * (
                                (Keff[j + 1][i] + Keff[j][i]) * 0.5 * (T[j + 1][i] - T[j][i]) +
                                Qbottom) +
                                1. / pow(dx, 2) * (
                                        (Keff[j][i + 1] + Keff[j][i]) * 0.5 * (T[j][i+  1] - T[j][i]) +
                                        (Keff[j][i - 1] + Keff[j][i]) * 0.5 * (T[j][i - 1] - T[j][i]))  ;


                    } else if (j == M - 1) { // top layer
                        Qflux = 1. / melt_layer[i] * (
                                Hn +
                                (Keff[j - 1][i] + Keff[j][i]) * 0.5 * (T[j - 1][i] - T[j][i]) / dy) +
                                1. / pow(dx, 2) * (
                                        (Keff[j][i + 1] + Keff[j][i]) * 0.5 * (T[j][i + 1] - T[j][i]) +
                                        (Keff[j][i - 1] + Keff[j][i]) * 0.5 * (T[j][i - 1] - T[j][i]));
                    }

                } else if (i == 0) { // West boundary, periodic boundary conditions

                    if (j > 0 && j < M - 1) {
                        Qflux = 1. / pow(dy, 2) * (
                                (Keff[j + 1][i] + Keff[j][i]) * 0.5 * (T[j + 1][i] - T[j][i]) +
                                (Keff[j - 1][i] + Keff[j][i]) * 0.5 * (T[j - 1][i] - T[j][i])) +
                                1. / pow(dx, 2) * (
                                        (Keff[j][i+1] + Keff[j][i]) * 0.5 * (T[j][i + 1] - T[j][i]) +
                                        (Keff[j][Nx - 1] + Keff[j][i]) * 0.5 * (T[j][Nx - 1] - T[j][i]))  ;

                    } else if (j == 0) {

                        Qflux = 1. / pow(dy, 2) * (
                                (Keff[j + 1][i] + Keff[j][i]) * 0.5 * (T[j + 1][i] - T[j][i]) +
                                Qbottom) +
                                1. / pow(dx, 2) * (
                                        (Keff[j][i+1] + Keff[j][i]) * 0.5 * (T[j][i + 1] - T[j][i]) +
                                        (Keff[j][Nx - 1] + Keff[j][i]) * 0.5 * (T[j][Nx - 1] - T[j][i]))  ;


                    } else if (j == M - 1) { // top layer
                        Qflux = 1. / melt_layer[i] * (
                                Hn +
                                (Keff[j - 1][i] + Keff[j][i]) * 0.5 * (T[j - 1][i] - T[j][i]) / dy) +
                                1. / pow(dx, 2) *
                                ((Keff[j][i + 1] + Keff[j][i]) * 0.5 * (T[j][i + 1] - T[j][i]) +
                                (Keff[j][Nx - 1] + Keff[j][i]) * 0.5 * (T[j][Nx - 1] - T[j][i]));
                    }

                } else if (i == Nx - 1 ){ // East boundary, periodic boundary conditions

                    if (j > 0 && j < M - 1) {
                        Qflux = 1. / pow(dy, 2) * (
                                (Keff[j + 1][i] + Keff[j][i]) * 0.5 * (T[j + 1][i] - T[j][i]) +
                                (Keff[j - 1][i] + Keff[j][i]) * 0.5 * (T[j - 1][i] - T[j][i])) +
                                1. / pow(dx, 2) * (
                                        (Keff[j][0] + Keff[j][i]) * 0.5 * (T[j][0] - T[j][i]) +
                                        (Keff[j][i - 1] + Keff[j][i]) * 0.5 * (T[j][i - 1] - T[j][i]))  ;

                    } else if (j == 0) {

                        Qflux = 1. / pow(dy, 2) * (
                                (Keff[j + 1][i] + Keff[j][i]) * 0.5 * (T[j + 1][i] - T[j][i]) +
                                Qbottom) +
                                1. / pow(dx, 2) * (
                                        (Keff[j][0] + Keff[j][i]) * 0.5 * (T[j][0] - T[j][i]) +
                                        (Keff[j][i - 1] + Keff[j][i]) * 0.5 * (T[j][i - 1] - T[j][i]))  ;


                    } else if (j == M - 1) {  // Top layer

                        Qflux = 1. / melt_layer[i] * (
                                Hn +
                                (Keff[j - 1][i] + Keff[j][i]) * 0.5 * (T[j - 1][i] - T[j][i]) / dy) +
                                1. / pow(dx, 2) * ((Keff[j][0] + Keff[j][i]) * 0.5 * (T[j][0] - T[j][i]) +
                                (Keff[j][i - 1] + Keff[j][i]) * 0.5 * (T[j][i - 1] - T[j][i]));

                    }
                }



                double rhoCp = rhoI * Cpi * (1. - porosity[j][i]) + rhoW * Cpw * theta[j][i];

                T_new[j][i] = T[j][i] + dt * Qflux / rhoCp;

            }

        }

        //=================================================================================
        //======================== Melt and upper BC for heat =============================
        //=================================================================================


        if (Hn > 0.) {

            int melt_flag = 0;

            for (int i = 0; i < Nx; i++) {

                    if (T_new[M-1][i] >= 0. and T[M-1][i] < 0.){
                        Hn = Hn_ini  + dryRho[M-1][i] * Cpi * T[M-1][i] * melt_layer[i] / dt;
                        T_new[M-1][i] = 0.;
                    } else if (T_new[M-1][i] >= 0. and T[M-1][i] >= 0.) {
                        T_new[M-1][i] = 0.;
                    } else if (T_new[M-1][i] == 0.) {
                        Hn = 0.;
                    }

                if (T[M - 1][i] == 0.) { // Upper layer is at 0oC, melt occurs

                    qIn[i] = melt(Hn, dryRho[M - 1][i], theta[M - 1][i], Ly[i], y[M-1][i], Ly_ini, dt);

                    melt_layer[i] = Ly[i] - (y[M - 2][i] + dy / 2.); // depth of the melting layer. Used to compute flux if melt_flag = 1

                    if (Ly[i] <= (y[M - 2][i] + dy / 2. + 0.002e-3) and
                        M > 2) { // the top layer disappears if it gets less than 0.01 mm
                        melt_flag = 1;
                    }

                }
            }

            if (melt_flag == 1) { // Merging of layers

                for (int i = 0; i < Nx; i++) {
                    qIn[i] += melt_layer[i] * (dryRho[M - 1][i] / 1000. + theta[M - 1][i]);
                }

                M = M - 1;

            }
        }




        //=================================================================================
        //============================== Full step hydrology ==============================
        //=================================================================================



        bool NotConverged = true; // Condition for convergence of the first and second order solutions of Heun's method

        while (NotConverged) {

            //==================== Compute theta_p (intermediate value) (forward Euler) =========================
            Eigen::VectorXd delta_theta_p = Richards(dx, dy, 0, M, melt_layer);


            int k = Nx*M - 1;

            for (int j = M-1; j >= 0; --j) {
                for (int i = Nx-1; i >= 0; --i) {

                    theta_p[j][i] = theta[j][i] + dt * delta_theta_p(k); // new theta predicted using Euler

                    // If it goes above saturation, move to boundary drying curve. Avoid model divergence
                    if (theta_p[j][i] > (0.9 * porosity[j][i] - epsilon) and j > 0) {

                        delta_theta_p(k-Nx) += (theta_p[j][i] - (0.9 * porosity[j][i] - epsilon))/dt;
                        theta_p[j][i] = 0.9 * porosity[j][i] - epsilon;
                        wrc[j][i] = 5; // Move directly to drying curve


                    } else if (theta_p[j][i] > (0.9 * porosity[j][i] - epsilon) and j == 0 ) {
                        theta_p[j][i] = 0.9 * porosity[j][i] - epsilon;
                        wrc[j][i] = 5; // Move directly to drying curve
                    }


//                     Compute capillary pressure and hydraulic conductivity for the intermediate water content (theta_p)
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


                    // If it goes above saturation, move to boundary drying curve. Avoid model divergence
                    if (theta_new[j][i] > (0.9 * porosity[j][i] - epsilon) and j > 0) { // If it goes above saturation, stop model
                        delta_theta(k-Nx) +=  (theta_new[j][i] - (0.9 * porosity[j][i] - epsilon))/ (0.5 * dt) ;
                        theta_new[j][i] = 0.9 * porosity[j][i] - epsilon;
                        wrc[j][i] = 5; // Move directly to drying curve
                    }  else if (theta_new[j][i] > (0.9 * porosity[j][i] - epsilon) and j == 0 ) {
                        theta_new[j][i] = 0.9 * porosity[j][i] - epsilon;
                        wrc[j][i] = 5; // Move directly to drying curve
                    }

                    //==================== Calculate error between intermediate value and final theta from Heun's method  =========================
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

        //=================================================================================
        //======================== Refreezing and updating variables  =====================
        //=================================================================================

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < Nx; ++i) {

                if (T_new[j][i] < 0.) {

                    double rhoCp = rhoI * Cpi * (1. - porosity[j][i]) + rhoW * Cpw * theta[j][i];

                    double theta_f = (-rhoCp * T_new[j][i]) / (Lf * rhoW); // total amount of liquid water that can refreeze

                    if (theta_f + epsilon < theta_new[j][i]) { //   not all liquid water content refreezes

                        T_new[j][i] = 0.;
                        theta_new[j][i] -= theta_f;
                        double thetaI = (1. - porosity[j][i]) + theta_f * rhoW / rhoI;
                        porosity[j][i] = 1. - thetaI;
                        dryRho[j][i] = thetaI * rhoI;

                    } else  {  // all liquid water content refreezes,  i.e. freeze = water_content

                        T_new[j][i] += (theta_new[j][i] - epsilon) * Lf * rhoW / rhoCp;
                        double thetaI = (1. - porosity[j][i]) + (theta_new[j][i] - epsilon) * rhoW / rhoI;
                        theta_new[j][i] = epsilon;
                        porosity[j][i] = 1. - thetaI;
                        dryRho[j][i] = thetaI * rhoI;

                    }

                }
            }
        }



        //==================== Update variables at the end of the time step =========================

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < Nx; ++i) {

                alpha[j][i] = 4.4e6 * pow(dryRho[j][i] / Dgrain[j][i], -0.98);
                alphaW[j][i] = gamma[j][i] * alpha[j][ i];
                n[j][i] = 1. + 2.7e-3 * pow(dryRho[j][i] / Dgrain[j][i], 0.61);
                Ks[j][i] = 9.81 / nu * 3. * pow(Dgrain[j][i] * 0.5, 2.) * exp(-0.013 * dryRho[j][i]);

                double Kw = 0.569;

                Keff[j][i] = (2.5e-6 * pow(dryRho[j][i], 2.) - 1.23e-4 * dryRho[j][i] + 0.024) * (1. - theta_new[j][i]) +
                             theta_new[j][i] * Kw;

                auto outputs = hysteresis(theta[j][i], theta_new[j][i], theta_s[j][i], theta_r[j][i], 0.9 * porosity[j][i],
                                          thetaR[j][i], wrc[j][i], thetaR_dry, i, j);

                std::tie(thetaR[j][i], theta_s[j][i], theta_r[j][i], wrc[j][i]) = outputs;

                psi[0][j][i] = psi[1][j][i];
                theta[j][i] = theta_new[j][i];
                Swf[0][j][i] = Swf[1][j][i];
                K[0][j][i]  =  KFun(Ks[j][i], 0, i, j);

                T[j][i] = T_new[j][i];

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

                std::cout << "Time = " << int(time / tPlot) << " min" << std::endl;
                height_output << M << " ";

                for (int i = 0; i < Nx; ++i) {
                    for (int j = 0; j < M; j++) {

                        theta_output << theta[j][i] << " ";
                        grain_output << Dgrain[j][i] << " ";
                        density_output << dryRho[j][i] << " " ;
                        T_output << T[j][i] << " " ;
                    }

                    outflow_output << K[0][0][i] << " " ;
                    height_output << Ly[i] << " ";

                }


                theta_output << std::endl;
                grain_output << std::endl;
                density_output << std::endl;
                outflow_output << std::endl;
                height_output << std::endl;
                T_output << std::endl;

            }

        }

        dt = dt_new;

    } // end of time loop


    theta_output.close();
    grain_output.close();
    density_output.close();
    height_output.close();
    T_output.close();

    printf("Time taken: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

}






