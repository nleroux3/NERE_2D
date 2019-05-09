#include "melt.h"
#include "global.h"


double melt(double& Hn, double& y, double& Ts, double T, double& Ly, double Keff, double dryRho,
            double theta, double dt) {

    double Vn;


    double  Ts_ini = Ts, // Temperature at the beginning of the time step
            Ly_ini = Ly;

    Ts = (Ly-y) / Keff * Hn + T;

    if (Ts >= 0. and Ts_ini < 0.) {
        Ts = 0.;
        Hn = Hn - (0. - T) * Keff / (Ly-y) ;
        Vn = Hn / (Lf * (dryRho + theta * rhoW)) ;
    }

    if (Ts >= 0. and Ts_ini >= 0.) {
        Ts = 0.;
        Vn = Hn / (Lf * (dryRho + theta * rhoW)) ;
    }

    if (Ts < 0.) {
        Vn = 0. ;
    }

    // Infiltration rate
     double Qin = Vn * (dryRho / rhoW + theta);

    // Change of snowdepth
    Ly = Ly_ini - Vn * dt;

    // Moving mesh
    y = y - (Ly_ini - Ly) * 0.5 ;


    return Qin;


}