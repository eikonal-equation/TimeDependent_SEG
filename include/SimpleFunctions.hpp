/*
* ==============================================================================
*
*  Copyright (C) 2019  Marc Aurèle Gilles
*  Copyright (C) 2019  Elliot Cartee, Qianli Song, Lexiao Lai
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
* ------------------------------------------------------------------------------
*
* The primary purpose in distributing this source code is to enable readers to
* reproduce the numerical results reported in the manuscript
* "Time-Dependent Surveillance-Evasion Games" by E. Cartee, L. Lai, Q. Song, 
* and A. Vladimirsky. See <https://arxiv.org/abs/1903.01332>.
*
* The code can be found at 
*       <https://github.com/eikonal-equation/TimeDependent_SEG>.
* Please see README.md for instructions on configuring/running this program.
*
* ------------------------------------------------------------------------------
*
* File: SimpleFunctions.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This file contains inline definitions of many small/simple 
* functions used for speed or pointwise observability.
*
* ==============================================================================
*/

#ifndef SIMPLEFUNCTIONS_HPP
#define SIMPLEFUNCTIONS_HPP

#include "GlobalConfiguration.hpp"

/** ----- Constant functions -------------------------------------------------*/
inline double constant(double aX) {
  return 1;
}

inline double constant2(double aX, double aY) {
  return 1;
}

inline double constant3(double aX, double aY, double aT) {
  return 1;
}

inline double constantcost3(int aI, int aJ, int aK) {
  return 1;
}

/** ----- Non-constant speed function ----------------------------------------*/
inline double wavy(double aX, double aY) {
  return 1 + 0.5*cos(aX*10*PI)*cos(aY*10*PI);
}

/** ------ Pointwise Observability -------------------------------------------*/
inline double exponential(double aPosX, double aPosY, double aX, double aY) {
  return exp(-20*((aPosX - aX)*(aPosX - aX) + (aPosY - aY)*(aPosY - aY)));
}

inline double inverseSquared(double aPosX, double aPosY, double aX, double aY) {
  return 1/(((aPosX - aX)*(aPosX - aX) + (aPosY - aY)*(aPosY - aY)) + 0.1);
}

inline double smallInverseSquared(double aPosX, double aPosY, 
                                  double aX, double aY) {
  return 1/(30*((aPosX - aX)*(aPosX - aX) + (aPosY - aY)*(aPosY - aY)) + 0.1);
}

#endif
