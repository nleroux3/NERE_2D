#include <iostream>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <tuple>
#include <omp.h>
#include <armadillo>
//#include <eigen3/Dense>

const int Nx = 50,  // Number of cells in lateral direction
          Ny = 50;  // Number of cells in vertical direction

double KFun(double , double , double);

std::tuple<double, double, double, double, int, double>
hyst_simple(double , double , double , double ,
            double , double , double , double , double ,
            int , double , double );


arma::vec Richards(double , double[Ny][Nx], double, double[Ny][Nx],
                double , double, double[Ny][Nx], double, int);


int main() {

    clock_t tStart = clock();



    double  Lx = 0.2, // [m]
            Ly_ini = 0.2,
            nu = 1.7e-6,
            lambda = 0,
            tau_0 = 0.1,
            gamma = 2.,
            tf = 7000. ,
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
            Dgrain[Ny][Nx], // [m]
            alpha[Ny][Nx],
            n[Ny][Nx],
            porosity[Ny][Nx],
            Ks[Ny][Nx],
            theta[Ny][Nx],
            psi[Ny][Nx],
            K[Ny][Nx],
            alphaW[Ny][Nx],
            thetaR[Ny][Nx],
            theta_s[Ny][Nx],
            theta_r[Ny][Nx],
            theta_s_new[Ny][Nx],
            theta_r_new[Ny][Nx],
            Swf[Ny][Nx],
            theta_p[Ny][Nx],
            theta_new[Ny][Nx], // predicted theta^(n+1) from Euler
            psi_new[Ny][Nx],
            Swf_new[Ny][Nx],
            K_new[Ny][Nx],
            thetaR_new[Ny][Nx];


    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(95., 105.);


    std::ofstream theta_output;

    theta_output.open("theta.dat");

    std::ofstream mass_output;

    mass_output.open("mass_conservation.dat");

    std::ofstream grain_output;

    grain_output.open("grain.dat");

    std::ofstream density_output;

    density_output.open("density.dat");

    std::ofstream outflow_output;

    outflow_output.open("outflow.dat");

    double dx =  Lx / double(Nx);
    double dy =  Ly_ini / double(Ny);


    arma::uword k = 0;

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

//            if (y[j][i] > 0.1){
//                dryRho[j][i] = 417.  * dist(mt)/100.;
//                Dgrain[j][i] = 0.406e-3  * dist(mt)/100.;
//            }

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

            Swf[j][i] = theta[j][i] / (0.9 * porosity[j][i]);


            auto outputs = hyst_simple(theta[j][i], theta[j][i], theta_s[j][i], theta_r[j][i], n[j][i], alpha[j][i], alphaW[j][i], 0.9 * porosity[j][i],
                                       thetaR[j][i], wrc[j][i], psi[j][i], thetaR_dry);

            std::tie(psi[j][i], thetaR[j][i], theta_s[j][i], theta_r[j][i], wrc[j][i], Swf[j][i]) = outputs;

            K[j][i] = KFun(Swf[j][i], Ks[j][i], n[j][i]);


            k += 1;

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


        arma::vec delta_theta_p = Richards(tau_0, Swf, lambda, K, dx, dy, psi, qIn, Ny);

        again_hydrology:

        k = 0;

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {

                theta_p[j][i] = theta[j][i] + dt * delta_theta_p(k); // new theta predicted using Euler

                if (theta_p[j][i] > (0.9 * porosity[j][i]) ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta)
                    dt = dt * 0.8;

                    goto again_hydrology;
                }

                if (theta_p[j][i] <= thetaR[j][i] and abs(theta_p[j][i] - thetaR[j][i])  >= 1e-6 ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta
                    dt = dt * 0.8;

                    goto again_hydrology;

                } else if (theta_p[j][i] <= thetaR[j][i] and abs(theta_p[j][i] - thetaR[j][i])  < 1e-6 ){
                    theta_p[j][i] = thetaR[j][i] + 1e-6;

                }

                compute_psi:

                auto outputs = hyst_simple(theta[j][i], theta_p[j][i], theta_s[j][i], theta_r[j][i],
                                           n[j][i], alpha[j][i], alphaW[j][i],
                                           0.9 * porosity[j][i], thetaR[j][i],
                                           wrc[j][i], psi[j][i], thetaR_dry);

                double thetaR1;


                std::tie(psi_new[j][i], thetaR1, std::ignore, std::ignore, std::ignore, Swf_new[j][i]) = outputs;

                if (theta_p[j][i] <= thetaR1 and abs(theta_p[j][i] - thetaR1)  >= 1e-6 ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta
                    dt = dt * 0.8;

                    goto again_hydrology;

                } else if (theta_p[j][i] <= thetaR1 and abs(theta_p[j][i] - thetaR1)  < 1e-6 ){
                    theta_p[j][i] = thetaR1 + 1e-6;

                    goto compute_psi;

                }

                K_new[j][i] = KFun(Swf_new[j][i], Ks[j][i], n[j][i]);

                k += 1;

            }
        }

        //==================== Compute theta_new (forward Heun) =========================

        arma::vec delta_theta = Richards(tau_0, Swf_new, lambda, K_new, dx, dy, psi_new, qIn, Ny);


        k = 0;

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {

                theta_new[j][i] = theta[j][i] + dt * 0.5 * (delta_theta(k) + delta_theta_p(k)) ; // new theta using Heun

                if (theta_new[j][i] > (0.9 * porosity[j][i])  ) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta)
                    dt = dt * 0.8;

                    goto again_hydrology;
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

            dt_new = dt * std::max(s * sqrt(tolerance / std::max(relative_error,EPS)), r_min);

            goto again_hydrology;

        }

        //==================== Calculate psi_new =========================

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {

                compute_psi1:

                auto outputs = hyst_simple(theta[j][i], theta_new[j][i], theta_s[j][i], theta_r[j][i], n[j][i],
                                           alpha[j][i], alphaW[j][i],
                                           0.9 * porosity[j][i], thetaR[j][i],
                                           wrc[j][i], psi[j][i], thetaR_dry);

                std::tie(psi_new[j][i], thetaR_new[j][i], theta_s_new[j][i], theta_r_new[j][i], wrc_new[j][i],
                         Swf_new[j][i]) = outputs;


                if (theta_new[j][i] <= thetaR_new[j][i] and abs(theta_new[j][i] - thetaR_new[j][i]) >= 1e-6) { // If it goes above saturation, excess of water goes to the cell below (included in delta_theta)
                    dt = dt * 0.8;//                    std::cout << 666 << std::endl;

                    goto again_hydrology;

                } else if (theta_new[j][i] <= thetaR_new[j][i] and abs(theta_new[j][i] - thetaR_new[j][i]) < 1e-6) {
                    theta_new[j][i] = thetaR_new[j][i] + 1e-6;

                    goto compute_psi1;
                }

            }
        }



        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {

                theta[j][i] = theta_new[j][i];

                thetaR[j][i] = thetaR_new[j][i];

                psi[j][i] = psi_new[j][i];

                Swf[j][i] = Swf_new[j][i];

                wrc[j][i] = wrc_new[j][i];

                theta_s[j][i] = theta_s_new[j][i];
                theta_r[j][i] = theta_r_new[j][i];

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

                    outflow_output << K[0][i] << " " ;
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


inline double KFun(double S, double Ks, double n) {
    double K,
            m;

    m = 1. - 1. / n;

    K = Ks * pow(S, 0.5) * pow(1. - pow(1. - pow(S, 1. / m), m), 2);

    return K;
}

//========================= SIMPLE HYSTERESIS ===============================
std::tuple<double, double, double, double, int, double>
hyst_simple(double theta_1, double theta_2, double theta_s, double theta_r, double n, double alpha, double alpha_wetting, double thetaS, double thetaR,
            int wrc_ini, double psi, double thetaR_dry) {

    double  m,
            Swf,
            Sd_Pid,
            Pw,
            Pd,
            Si_Pdi,
            Pdi,
            Pid;

    int wrc;


    m = 1. - 1. / n;

    if (wrc_ini == 1) {  // On main wetting curve

        if (theta_2 >= theta_1) {   // Stays on main wetting curve

            wrc = 1;
            thetaR = 0.;
            theta_s = thetaS;
            theta_r = thetaR;

            Swf = theta_2 / thetaS;

            psi = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n) / alpha_wetting;

         } else { // Moves to scanning drying curve

            wrc = 3;
            Pid = psi;

            thetaR = theta_1 / (1. + (thetaS / thetaR_dry - 1.) * theta_1 / thetaS);
            theta_r = thetaR;

            Swf = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                   (thetaS / thetaR_dry - 1.) * (theta_2 / thetaS - thetaR / thetaS)));

            Sd_Pid = pow(1. + pow(-alpha * Pid, n), -m);

            theta_s = (theta_1 - thetaR * (1.-Sd_Pid)) / Sd_Pid;

             psi = -pow(pow((theta_2 - thetaR) / (theta_s - thetaR), -1. / m) - 1., 1. / n) /
                  alpha;


        }

    }


    if (wrc_ini == 3) { // On drying scanning curve

        wrc = 3;

        psi = -pow(pow((theta_2 - thetaR) / (theta_s - thetaR), -1. / m) - 1., 1 / n) /
              alpha;

        Swf = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                               (thetaS /
                                                                                                                thetaR_dry -
                                                                                                                1.) *
                                                                                                               (theta_2 /
                                                                                                                thetaS -
                                                                                                                thetaR /
                                                                                                                thetaS)));

        Pw = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n) / alpha_wetting;

        if (abs(psi) < abs(Pw) || theta_2 > theta_s) { // Moves to main wetting curve
            wrc = 1;
            thetaR = 0.;
            psi = Pw;
        }
    }


    if (wrc_ini == 4) { //  On wetting scanning scanning

        wrc = 4;

        psi = -pow(pow((theta_2 - theta_r) / (thetaS - theta_r), -1. / m) - 1., 1. / n) /
              alpha_wetting;

        Swf = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                               (thetaS /
                                                                                                                thetaR_dry -
                                                                                                                1.) *
                                                                                                               (theta_2 /
                                                                                                                thetaS -
                                                                                                                thetaR /
                                                                                                                thetaS)));

        Pw = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n) / alpha_wetting;
        Pd = -pow(pow((theta_2 - thetaR_dry) / (thetaS - thetaR_dry), -1. / m) - 1., 1. / n) / alpha;


        if (abs(psi) < abs(Pw)) { //  Moves to main wetting curve
            wrc = 1;
            thetaR = 0;
            psi = Pw;

        } else if (abs(psi) > abs(Pd) or theta_2 < theta_r) { //  Moves to main wetting curve
            wrc = 2;
            thetaR = thetaR_dry;
            psi = Pd;
        }


    }


    if (wrc_ini == 2) { // On main drying curve

        wrc = 2;

        thetaR = thetaR_dry;
        theta_s = thetaS;
        theta_r = thetaR;


        if (theta_2 > theta_1) { // WRC moves to wetting scanning (di: drainage to imbibition)

            wrc = 4;

            Pdi = psi;

            Swf = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                                   (thetaS /
                                                                                                                    thetaR_dry -
                                                                                                                    1.) *
                                                                                                                   (theta_2 /
                                                                                                                    thetaS -
                                                                                                                    thetaR /
                                                                                                                    thetaS)));

            Si_Pdi = pow(1. + pow(-alpha_wetting * Pdi, n), -m);

            theta_r = (theta_1 - thetaS * Si_Pdi)/(1. - Si_Pdi);

            psi = -pow(pow((theta_2 - theta_r) / (thetaS - theta_r), -1. / m) - 1., 1. / n) /
                  alpha_wetting;

            Pw = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n) / alpha_wetting;

            if (abs(psi) < abs(Pw)) { // Moves to main wetting curve
                wrc = 1;
                psi = Pw;
            }


        } else { // Stays on main drying curve

            Swf = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                                   (thetaS /
                                                                                                                    thetaR_dry -
                                                                                                                    1.) *
                                                                                                                   (theta_2 /
                                                                                                                    thetaS -
                                                                                                                    thetaR /
                                                                                                                    thetaS)));

            psi = -pow(pow((theta_2 - thetaR) / (thetaS - thetaR), -1. / m) - 1., 1. / n) / alpha;

        }
    }

    return std::make_tuple(psi, thetaR, theta_s, theta_r, wrc, Swf);
}

//========================= Richards equation solver  ===============================
arma::vec Richards( double tau_0, double Swf[Ny][Nx], double lambda, double K[Ny][Nx],
                double dx, double dy, double psi[Ny][Nx], double qIn, const int M){

    double tau[M][Nx],
           KSouth[M][Nx],
           KNorth[M][Nx],
           KEast[M][Nx],
           KWest[M][Nx];


    arma::sp_mat Mat(Nx * M, Nx * M);

    arma::vec F(M * Nx),
              delta_theta(M*Nx);



    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < Nx; ++i) {

            tau[j][i] = tau_0 * pow(Swf[j][i], -lambda);


            if (j > 0) {
                KSouth[j][i] = pow(K[j - 1][i] * K[j][i], 0.5);
//                KSouth[j][i] = (K[j - 1][i] + K[j][i])  * 0.5;

            } else {
                KSouth[j][i] = K[j][i];
            }

            if (j < M - 1) {
                KNorth[j][i] = pow(K[j + 1][i] * K[j][i], 0.5);
//                KNorth[j][i] = (K[j + 1][i] + K[j][i]) * 0.5;

            } else {
                KNorth[j][i] = K[j][i];
            }

            if (i > 0) {
                KWest[j][i] = pow(K[j][i - 1] * K[j][i], 0.5);
//                KWest[j][i] = (K[j][i - 1] + K[j][i]) * 0.5;

            } else {
                KWest[j][i] = 0.;
            }

            if (i < Nx - 1) {
                KEast[j][i] = pow(K[j][i + 1] * K[j][i], 0.5);
//                KEast[j][i] = (K[j][i + 1] + K[j][i]) * 0.5;

            } else {
                KEast[j][i] = 0.;
            }

        }
    }

    arma::uword k = 0;

    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < Nx;++i) {

            if (j == 0) { //======== BOTTOM BOUNDARY ==================================

                Mat(k, k) = 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                              + KNorth[j][i] / pow(dy, 2.));  // (i,j)

                if (i > 0 && i < Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][i + 1] - psi[j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[j][i] - psi[j + 1][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else if (i == 0) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][i + 1] - psi[j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[j][i] - psi[j + 1][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else if (i == Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][0] - psi[j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[j][i] - psi[j + 1][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                }

            } else if (j == M - 1) { //======== TOP BOUNDARY ==================================

                Mat(k, k) = 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                              + KSouth[j][i] / pow(dy, 2.));  // (i,j)

                if (i > 0 && i < Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][i + 1] - psi[j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[j - 1][i] - psi[j][i]) / pow(dy, 2.)
                            + (qIn - KSouth[j][i]) / dy);
                } else if (i == 0) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][i + 1] - psi[j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[j - 1][i] - psi[j][i]) / pow(dy, 2.)
                            + (qIn - KSouth[j][i]) / dy);
                } else if (i == Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][0] - psi[j][i]) / pow(dx, 2.)
                            + KSouth[j][i] * (psi[j - 1][i] - psi[j][i]) / pow(dy, 2.)
                            + (qIn - KSouth[j][i]) / dy);
                }


            } else {

                Mat(k, k) = 1. + tau[j][i] * ((KWest[j][i] + KEast[j][i]) / pow(dx, 2.)
                                              + (KSouth[j][i] + KNorth[j][i]) / pow(dy, 2.));  // (i,j)

                if (i == 0) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][i + 1] - psi[j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[j][i] - psi[j + 1][i]) / pow(dy, 2.)
                            + KSouth[j][i] * (psi[j - 1][i] - psi[j][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else if (i == Nx - 1) {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][Nx - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][0] - psi[j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[j][i] - psi[j + 1][i]) / pow(dy, 2.)
                            + KSouth[j][i] * (psi[j - 1][i] - psi[j][i]) / pow(dy, 2.)
                            + (KNorth[j][i] - KSouth[j][i]) / dy);
                } else {
                    F(k) = (-KWest[j][i] * (psi[j][i] - psi[j][i - 1]) / pow(dx, 2.)
                            + KEast[j][i] * (psi[j][i + 1] - psi[j][i]) / pow(dx, 2.)
                            - KNorth[j][i] * (psi[j][i] - psi[j + 1][i]) / pow(dy, 2.)
                            + KSouth[j][i] * (psi[j - 1][i] - psi[j][i]) / pow(dy, 2.)
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

            if (j < M - 1) {
                Mat(k, k + Nx) = -tau[j + 1][i] * KNorth[j][i] / pow(dy, 2.); // (i,j+1)
            }


            k += 1;
        }
    }

        delta_theta = arma::spsolve(Mat, F, "superlu");


        return delta_theta;

}


