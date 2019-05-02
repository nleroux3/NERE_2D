//
// Created by Nicolas Leroux on 2019-04-16.
//
#include "KFun.h"
#include "global.h"
#include <math.h>

double KFun(const double& Ks, const int& p, const int& i, const int& j) {

    double K,
           m;

    m = 1. - 1. / n[j][i];

    K = Ks * sqrt(Swf[p][j][i]) * pow(1. - pow(1. - pow(Swf[p][j][i], 1. / m), m), 2);

    return K;
}