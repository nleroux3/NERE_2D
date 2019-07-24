//
// Created by Nicolas Leroux on 2019-04-16.
//
#include <iostream>
#include <vector>
#pragma once

extern int Nx, Ny;

const double  rhoI = 917., // density of ice [kg/m3]
              Cpi = 2110., // Heat capacity of ice [J/kg/K]
              Lf = 338000., // Latent heat of fusion [J/kg]
              rhoW = 1000., // density of water [kg/m3]
              Cpw = 4181.3;  // Heat capacity of water [J/kg/K]


extern std::vector<std::vector<std::vector<double>>> psi,
                                                     K,
                                                     Swf;


extern std::vector<std::vector<double>> alpha,
                                        n,
                                        porosity,
                                        Ks,
                                        alphaW,
                                        tau,
                                        T,
                                        Keff;

extern std::vector<double> qIn,
                           Ts;

