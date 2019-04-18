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
* File: CStationaryObstacle.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This class implements the CObstacle interface for a stationary 
* rectangular obstacle.
*
* ==============================================================================
*/

#ifndef CSTATIONARYOBSTACLE_HPP
#define CSTATIONARYOBSTACLE_HPP

#include "GlobalConfiguration.hpp"
#include "CObstacle.hpp"

class CStationaryObstacle : public CObstacle {
  public:
    /* Returns true if (aX,aY,aT) is on the boundary, false otherwise */
    virtual bool isOnBoundary(const double aX, const double aY, 
                              const double aT = 0) const override;

    /** Returns true if (aX,aY,aT) is inside the obstacle, false otherwise */
    virtual bool isIn(const double aX, const double aY, 
                      const double aT = 0) const override;

    /** Constructors */
    CStationaryObstacle() = default;
    CStationaryObstacle(const double aLeft, const double aRight, 
                        const double aBottom, const double aTop);

    /** Destructor */
    virtual ~CStationaryObstacle();

  private:
    /** The position of the obstacle : [fLeft,fRight]x[fTop,fBottom] */
    double fLeft, fRight, fTop, fBottom; 

    /** Returns whether aX is within distance 10^(-8) of aY */
    bool nearly_equals(const double aX, const double aY) const;
};

#endif
