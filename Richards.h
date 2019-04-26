//
// Created by Nicolas Leroux on 2019-04-16.
//

#pragma once
#include <iostream>
#include <Eigen/Sparse>

Eigen::VectorXd Richards( const double& tau_0,  const double& lambda, const double& dx, const double& dy, const double& qIn, const int& p);