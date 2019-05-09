//
// Created by Nicolas Leroux on 2019-04-16.
//
#pragma once
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
using namespace Eigen;

const int Nx = 100,  // Number of cells in lateral direction
          Ny = 100;  // Number of cells in vertical direction


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
