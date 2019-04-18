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
* File: CStationaryObstacle.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This class implements the CObstacle interface for a stationary 
* rectangular obstacle.
* (See also CStationaryObstacle.cpp)
*
* ==============================================================================
*/
#include "CStationaryObstacle.hpp"

/*==============================================================================
  Constructor
==============================================================================*/
CStationaryObstacle::CStationaryObstacle(const double Left, const double Right, 
                                         const double Bottom, const double Top) {
  fLeft = Left;
  fRight = Right;
  fBottom = Bottom;
  fTop = Top;
}

/*==============================================================================
  Destructor
==============================================================================*/
CStationaryObstacle::~CStationaryObstacle(){}

/*==============================================================================
  Check if (aX,aY,aT) is in Obstacle
==============================================================================*/
bool CStationaryObstacle::isIn(const double aX, const double aY, 
                               const double aT) const {
  return (aX >= fLeft && aX <= fRight && aY >= fBottom && aY <= fTop);
}

/*==============================================================================
  Check if (aX,aY,aT) is on boundary of Obstacle
==============================================================================*/
bool CStationaryObstacle::isOnBoundary(const double aX, const double aY, 
                                       const double aT) const {
  if ((nearly_equals(aX,fLeft) || nearly_equals(aX,fRight)) && (aY >= fBottom && aY <= fTop)) {
    return 1;
  } else if ((aX >= fLeft && aX <= fRight) && (nearly_equals(aY,fBottom) || nearly_equals(aY,fTop))) {
    return 1;
  } else {
  	return 0;
  }
}

/*==============================================================================
  Helper function for computing whether aX is within distance epsilon of aY
==============================================================================*/
constexpr double epsilon = 1e-8;
bool CStationaryObstacle::nearly_equals(const double aX, const double aY) const {
  return ((aX + epsilon >= aY) && (aX - epsilon <= aY));
}