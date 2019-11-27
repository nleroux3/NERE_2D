//
// Created by Nicolas Leroux on 2019-04-16.
//

#pragma once
#include <iostream>
#include "global.h"
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>

Eigen::VectorXd Richards(const double& dx, const double& dy, const int& p, const int&  M, double dy_top[]);