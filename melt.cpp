#include "melt.h"
#include "global.h"
#include <iostream>


double melt(double& Hn,  double dryRho, double theta, double& Ly, double& y, double Ly_ini, double dt) {

    double Vn;



    Vn = Hn / (Lf * (dryRho + theta * rhoW)) ;

    // Infiltration rate
     double Qin = Vn * (dryRho / rhoW + theta);

    // Change of snowdepth
    Ly = Ly_ini - Vn * dt;

    // Moving mesh
    y = y - (Ly_ini - Ly) * 0.5 ;


    return Qin;


}