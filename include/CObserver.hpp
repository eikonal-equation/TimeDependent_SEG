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
* File: CObserver.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This is an abstract class interface for the observers.
* Sub-classes will need to implement functions for keeping track of the 
* observer's position, and for calculating the pointwise observability function.
*
* ==============================================================================
*/

#ifndef COBSERVER_HPP
#define COBSERVER_HPP

/** ------ Project-specific header files -------------------------------------*/
#include "CTerrain.hpp"

class CObserver
{
  public:
    /** ------ Getters -------------------------------------------------------*/
    /** Pointwise observability for gridpoint (aI,aJ,aK) */
    virtual double getObservabilityGrid(const int aI, const int aJ, 
                                        const int aK) const = 0;

    /** Pointwise observability at physical location (aX,aY,aT) */
    virtual double getObservabilityPhysical(const double aX, const double aY, 
                                            const double aT) const = 0;

    /** Path followed by the observer */
    virtual double getPosX(const double aT) const = 0;
    virtual double getPosY(const double aT) const = 0;

    /** ------ Destructor ----------------------------------------------------*/
    virtual ~CObserver() {};

  protected:
    /** ------Fields ---------------------------------------------------------*/
    /** Path (x(t),y(t)) followed by observer */
  	std::function<double(double)> fXpath, fYpath;
    /** The pointwise observability function */
    std::function<double(double,double,double,double)> fObservability;
    /** Shared pointer to terrain function */
    std::shared_ptr<CTerrain> fTerrain;
};

#endif
