//
// Created by Nicolas Leroux on 2019-04-16.
//
#include <iostream>
#include <tuple>
#include <math.h>

#include "global.h"

//========================= SIMPLE HYSTERESIS ===============================
std::tuple<double, double, double, int>
hysteresis(const double& theta_1, const double& theta_2,  double theta_s,  double theta_r, const double& thetaS,  double thetaR,
           const int& wrc_ini, const double& thetaR_dry, const int& i, const int& j) {

    double  m,
            Sd_Pid,
            Pw,
            Pd,
            Si_Pdi,
            Pdi,
            Pid;

    int wrc;


    m = 1. - 1. / n[j][i];

    if (wrc_ini == 1) {  // On main wetting curve

        if (theta_2 >= theta_1) {   // Stays on main wetting curve

            wrc = 1;
            thetaR = 0.;
            theta_s = thetaS;
            theta_r = thetaR;

            Swf[1][j][i] = theta_2 / thetaS;

            psi[1][j][i] = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n[j][i]) / alphaW[j][i];

        } else { // Moves to scanning drying curve

            wrc = 3;
            Pid = psi[0][j][i];

            thetaR = std::min(theta_1 / (1. + (thetaS / thetaR_dry - 1.) * theta_1 / thetaS) , theta_2 - 1e-10);
            theta_r = thetaR;

            Swf[1][j][i] = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                                   (thetaS / thetaR_dry - 1.) * (theta_2 / thetaS - thetaR / thetaS)));

            Sd_Pid = pow(1. + pow(-alpha[j][i]* Pid, n[j][i]), -m);

            theta_s = (theta_1 - thetaR * (1.-Sd_Pid)) / Sd_Pid;

            psi[1][j][i] = -pow(pow((theta_2 - thetaR) / (theta_s - thetaR), -1. / m) - 1., 1. / n[j][i]) /
                  alpha[j][i];


        }

    }


    if (wrc_ini == 3) { // On drying scanning curve

        wrc = 3;

        psi[1][j][i] = -pow(pow((theta_2 - thetaR) / (theta_s - thetaR), -1. / m) - 1., 1 / n[j][i]) /
              alpha[j][i];

        Swf[1][j][i] = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                               (thetaS /
                                                                                                                thetaR_dry -
                                                                                                                1.) *
                                                                                                               (theta_2 /
                                                                                                                thetaS -
                                                                                                                thetaR /
                                                                                                                thetaS)));

        Pw = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n[j][i]) / alphaW[j][i];

        if (abs(psi[1][j][i]) < abs(Pw) || theta_2 > theta_s) { // Moves to main wetting curve
            wrc = 1;
            thetaR = 0.;
            psi[1][j][i] = Pw;
        }
    }


    if (wrc_ini == 4) { //  On wetting scanning scanning

        wrc = 4;

        psi[1][j][i] = -pow(pow((theta_2 - theta_r) / (thetaS - theta_r), -1. / m) - 1., 1. / n[j][i]) /
              alphaW[j][i];

        Swf[1][j][i] = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                               (thetaS /
                                                                                                                thetaR_dry -
                                                                                                                1.) *
                                                                                                               (theta_2 /
                                                                                                                thetaS -
                                                                                                                thetaR /
                                                                                                                thetaS)));

        Pw = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n[j][i]) / alphaW[j][i];
        Pd = -pow(pow((theta_2 - thetaR_dry) / (thetaS - thetaR_dry), -1. / m) - 1., 1. / n[j][i]) / alpha[j][i];


        if (abs(psi[1][j][i]) < abs(Pw)) { //  Moves to main wetting curve
            wrc = 1;
            thetaR = 0;
            psi[1][j][i] = Pw;

        } else if (abs(psi[1][j][i]) > abs(Pd) or theta_2 < theta_r) { //  Moves to main wetting curve
            wrc = 2;
            thetaR = thetaR_dry;
            psi[1][j][i] = Pd;
        }


    }


    if (wrc_ini == 2) { // On main drying curve

        wrc = 2;

        thetaR = thetaR_dry;
        theta_s = thetaS;
        theta_r = thetaR;


        if (theta_2 > theta_1) { // WRC moves to wetting scanning (di: drainage to imbibition)

            wrc = 4;

            Pdi = psi[0][j][i];

            Swf[1][j][i] = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                                   (thetaS /
                                                                                                                    thetaR_dry -
                                                                                                                    1.) *
                                                                                                                   (theta_2 /
                                                                                                                    thetaS -
                                                                                                                    thetaR /
                                                                                                                    thetaS)));

            Si_Pdi = pow(1. + pow(-alphaW[j][i] * Pdi, n[j][i]), -m);

            theta_r = std::min((theta_1 - thetaS * Si_Pdi)/(1. - Si_Pdi), theta_2 - 1e-10);

            psi[1][j][i] = -pow(pow((theta_2 - theta_r) / (thetaS - theta_r), -1. / m) - 1., 1. / n[j][i]) /
                  alphaW[j][i];

            Pw = -pow(pow(theta_2 / thetaS, -1. / m) - 1., 1. / n[j][i]) / alphaW[j][i];

            if (abs(psi[1][j][i]) < abs(Pw)) { // Moves to main wetting curve
                wrc = 1;
                psi[1][j][i] = Pw;
            }


        } else { // Stays on main drying curve

            Swf[1][j][i] = 0.5 * ((theta_2 / thetaS - thetaR / thetaS) + sqrt(pow(theta_2 / thetaS - thetaR / thetaS, 2.) + 4. /
                                                                                                                   (thetaS /
                                                                                                                    thetaR_dry -
                                                                                                                    1.) *
                                                                                                                   (theta_2 /
                                                                                                                    thetaS -
                                                                                                                    thetaR /
                                                                                                                    thetaS)));

            psi[1][j][i] = -pow(pow((theta_2 - thetaR) / (thetaS - thetaR), -1. / m) - 1., 1. / n[j][i]) / alpha[j][i];

        }
    }

    return std::make_tuple(thetaR, theta_s, theta_r, wrc);
}