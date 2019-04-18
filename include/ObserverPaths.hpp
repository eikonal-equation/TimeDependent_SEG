/*
* ==============================================================================
*
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
* File: ObserverPaths.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*
* Description: This file contains all of the position and velocity functions for
* all of the observer trajectories.
*
* ==============================================================================
*/

#ifndef OBSERVERPATHS_HPP
#define OBSERVERPATHS_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"

/*==============================================================================
      Simple Paths
==============================================================================*/

/** ------ Simple path 1 -----------------------------------------------------*/
inline double xsimple1(double t){
	return 0.3 + 0.12*cos(2.0*PI*t);
}
inline double ysimple1(double t){
	return 0.8 + 0.12*sin(2.0*PI*t);
}
inline double xsimpleder1(double t){
  return -0.12*2.0*PI*sin(2.0*PI*t);
}
inline double ysimpleder1(double t){
  return 0.12*2.0*PI*cos(2.0*PI*t);
}

/** ------ Simple path 2 -----------------------------------------------------*/
inline double xsimple2(double t){
	return 0.9 - 0.08*cos(2.0*PI*t);
}
inline double ysimple2(double t){
  return 0.1 - 0.08*sin(2.0*PI*t);
}
inline double xsimpleder2(double t){
  return 0.08*2.0*PI*sin(2.0*PI*t);
}
inline double ysimpleder2(double t){
  return -0.08*2.0*PI*cos(2.0*PI*t);
}

/*==============================================================================
      Tilted Figure Eight Shaped Paths
==============================================================================*/

/** ------ Figure eight constants --------------------------------------------*/
constexpr double r1 = (1.0/3-0.05)/(1+1.0/SQRT2);
constexpr double r2 = SQRT2/6;

/** ------ Tilted Figure Eight path 1 ----------------------------------------*/
inline double xtilted1(double t) {
  if (t<1) {
    return 0.05 +r1 + r1*cos(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return 0.5+r2*cos(0.75*PI-PI*(t-1));
  } else {
    return 0.05+r1+r1*cos(-PI/4+PI*(t-3));
  }
}
inline double ytilted1(double t){
  if (t<1) {
    return 0.95-r1+r1*sin(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return 0.5+r2*sin(0.75*PI-PI*(t-1));
  } else {
    return 0.95-r1+r1*sin(-PI/4+PI*(t-3));
  }
}
inline double xtiltedder1(double t){
  if (t<1) {
    return -r1*PI*sin(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return r2*PI*sin(0.75*PI-PI*(t-1));
  } else {
    return -r1*PI*sin(-PI/4+PI*(t-3));
  }
}
inline double ytiltedder1(double t){
  if (t<1) {
    return PI*r1*cos(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return -PI*r2*cos(0.75*PI-PI*(t-1));
  } else {
    return PI*r1*cos(-PI/4+PI*(t-3));
  }
}

/** ------ Tilted Figure Eight path 2 ----------------------------------------*/
inline double xtilted2(double t){
  if (t<1) {
    return 0.95-r1+r1*sin(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return 0.5+r2*sin(0.75*PI-PI*(t-1));
  } else {
    return 0.95-r1+r1*sin(-PI/4+PI*(t-3));
  }
}
inline double ytilted2(double t){
  if (t<1) {
    return 0.05+r1 + r1*cos(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return 0.5+r2*cos(0.75*PI-PI*(t-1));
  } else {
    return 0.05+r1+r1*cos(-PI/4+PI*(t-3));
  }
}
inline double xtiltedder2(double t){
  if (t<1) {
    return PI*r1*cos(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return -PI*r2*cos(0.75*PI-PI*(t-1));
  } else {
    return PI*r1*cos(-PI/4+PI*(t-3));
  }
}
inline double ytiltedder2(double t){
  if (t<1) {
    return -r1*PI*sin(0.75*PI+PI*t);
  } else if (t>=1 && t<3) {
    return r2*PI*sin(0.75*PI-PI*(t-1));
  } else {
    return -r1*PI*sin(-PI/4+PI*(t-3));
  }
}

/**=============================================================================
      Circular paths for maze example
==============================================================================*/

/** ------ Maze circle constants ---------------------------------------------*/
constexpr double xc = 0.4;
constexpr double yc = 0.525;
constexpr double rad = 0.25;

/** ------ Maze circle 1 -----------------------------------------------------*/
inline double xmazecircle1(double t) {
  return xc + rad*cos(t);
}
inline double ymazecircle1(double t) {
  return yc + rad*sin(t);
}
inline double xmazecircleder1(double t) {
  return -rad*sin(t);
}
inline double ymazecircleder1(double t) {
  return rad*cos(t);
}

/** ------ Maze circle 2 -----------------------------------------------------*/
inline double xmazecircle2(double t) {
  return xc + rad*cos(t + PI/2);
}
inline double ymazecircle2(double t) {
  return yc + rad*sin(t + PI/2);
}
inline double xmazecircleder2(double t) {
  return -rad*sin(t + PI/2);
}
inline double ymazecircleder2(double t) {
  return rad*cos(t + PI/2);
}

/** ------ Maze circle 3 -----------------------------------------------------*/
inline double xmazecircle3(double t) {
  return xc + rad*cos(t + PI);
}
inline double ymazecircle3(double t) {
  return yc + rad*sin(t + PI);
}
inline double xmazecircleder3(double t) {
  return -rad*sin(t + PI);
}
inline double ymazecircleder3(double t) {
  return rad*cos(t + PI);
}

/** ----- Maze circle 4 ------------------------------------------------------*/
inline double xmazecircle4(double t) {
  return xc + rad*cos(t + 3*PI/2);
}
inline double ymazecircle4(double t) {
  return yc + rad*sin(t + 3*PI/2);
}
inline double xmazecircleder4(double t) {
  return -rad*sin(t + 3*PI/2);
}
inline double ymazecircleder4(double t) {
  return rad*cos(t + 3*PI/2);
}

#endif