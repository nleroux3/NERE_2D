//
// Created by Nicolas Leroux on 2019-04-16.
//
#include "KFun.h"
#include "global.h"
#include <math.h>

double KFun( double Ks, int p, int i, int j) {

    double K,
           m;

    m = 1. - 1. / n[j][i];

    K = Ks * pow(Swf[p][j][i], 0.5) * pow(1. - pow(1. - pow(Swf[p][j][i], 1. / m), m), 2);

    return K;
}