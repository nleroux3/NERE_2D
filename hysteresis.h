//
// Created by Nicolas Leroux on 2019-04-16.
//

#pragma once
#include <tuple>


std::tuple<double, double, double, int>
hysteresis(const double& theta_1, const double& theta_2,  double theta_s,  double theta_r, const double& thetaS,  double thetaR,
           const int& wrc_ini, const double& thetaR_dry, const int& i, const int& j);