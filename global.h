//
// Created by Nicolas Leroux on 2019-04-16.
//
#pragma once

const int Nx = 50,  // Number of cells in lateral direction
          Ny = 50;  // Number of cells in vertical direction


extern double  alpha[Ny][Nx],
               n[Ny][Nx],
               porosity[Ny][Nx],
               Ks[Ny][Nx],
               psi[2][Ny][Nx],
               K[2][Ny][Nx],
               alphaW[Ny][Nx],
               Swf[2][Ny][Nx];