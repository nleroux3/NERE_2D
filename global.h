//
// Created by Nicolas Leroux on 2019-04-16.
//
#pragma once
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
using namespace Eigen;

const int Nx = 1,  // Number of cells in lateral direction
          Ny = 2000;  // Number of cells in vertical direction


extern double  alpha[Ny][Nx],
               n[Ny][Nx],
               porosity[Ny][Nx],
               Ks[Ny][Nx],
               psi[2][Ny][Nx],
               K[2][Ny][Nx],
               alphaW[Ny][Nx],
               Swf[2][Ny][Nx],
               tau[Ny][Nx],
               qIn[Nx],
               Ts[Nx],
               T[Ny][Nx],
               Keff[Ny][Nx]; // thermal conductivity [W/(m K)];


const double  rhoI = 917., // density of ice [kg/m3]
              Cpi = 2110., // Heat capacity of ice [J/kg/K]
              Lf = 338000., // Latent heat of fusion [J/kg]
              rhoW = 1000., // density of water [kg/m3]
              Cpw = 4181.3;  // Heat capacity of water [J/kg/K]