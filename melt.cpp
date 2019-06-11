#include "melt.h"
#include "global.h"
#include <iostream>


double melt(double& Hn,  double dryRho, double theta, double& Ly, double& y, double Ly_ini, double dt) {

    double Vn;


//    double  Ts_ini = T, // Temperature at the beginning of the time step
//            Ly_ini = Ly;
//
//    Ts = (Ly-y) / Keff * Hn + T;
//
//    if (Ts >= 0. and Ts_ini < 0.) {
//        T = 0.;
//        Hn = Hn - (0. - T) * Keff / (Ly-y) ;
//        Vn = Hn / (Lf * (dryRho + theta * rhoW)) ;
//    }
//
//    if (Ts >= 0. and Ts_ini >= 0.) {
//        Ts = 0.;
//        Vn = Hn / (Lf * (dryRho + theta * rhoW)) ;
//    }
//
//    if (Ts < 0.) {
//        Vn = 0. ;
//    }


    Vn = Hn / (Lf * (dryRho + theta * rhoW)) ;

    // Infiltration rate
     double Qin = Vn * (dryRho / rhoW + theta);

    // Change of snowdepth
    Ly = Ly_ini - Vn * dt;

    // Moving mesh
    y = y - (Ly_ini - Ly) * 0.5 ;


    return Qin;


}