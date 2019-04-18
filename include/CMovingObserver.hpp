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
* File: CMovingObserver.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This is the class for moving Observers. It implements the 
* CObserver interface, and calculates the pointwise observability in the domain.
*
* ==============================================================================
*/

#ifndef CMOVINGOBSERVER_HPP
#define CMOVINGOBSERVER_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <limits>

/** ----- Project-specific header files --------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "MemoryAllocations.hpp"
#include "CTimeGrid.hpp"
#include "CTerrain.hpp"
#include "SimpleFunctions.hpp"
#include "CObserver.hpp"

class CMovingObserver : public CObserver
{
  public:
    /* Default constructor */
    CMovingObserver() = default;

    /**
     * This is the main constructor for moving observers. 
     * It initializes the properties and computes the visibility information.
     * @param aXpath function x(t), physical location of the observer at time t
     * @param aYpath function y(t), physical location of the observer at time t
     * @param aTerrain CTerrain object containing the obstacles.
     * @param aObservability  the pointwise observability function. 
     *        Defaults to inverse Squared.
     * @param aN size of grid on which visibility is computed. Default is 51.
     * @param aNt number of timesteps. Default is 101.
     * @param aPhysMin lower bound of physical domain. Defaults to 0.
     * @param aPhysMax upper bound of physical domain. Defaults to 1.
     * @param aTime time interval of problem. Defaults to 4.
     The domain is [aPhysMin, aPhysMax]x[aPhysMin,aPhysMax]x[0,aTime]. 
     * @param aXder function x'(t), the observer's horizontal velocity. 
     *        Defaults to 1. Only relevant if alpha < 2*pi.
     * @param aYder function y'(t), the observer's vertical velocity. 
     *        Defaults to 1. Only relevant if alpha < 2*pi.
     * @param aAlpha the angle of the visible sector. Defaults to 2pi.
     */
    CMovingObserver(const std::function<double(double)> aXpath, 
                    const std::function<double(double)> aYpath, 
                    const std::shared_ptr<CTerrain> aTerrain, 
                    const std::function<double(double,double,double,double)> aObservability = &inverseSquared,
                    const int aN = 51, const int aNt = 101, 
                    const double aPhysMin = 0, const double aPhysMax = 1, 
                    const double aTime = 4, 
                    const std::function<double(double)> aXder = &constant, 
                    const std::function<double(double)> aYder = &constant, 
                    const double aAlpha = 2*PI);

    /** ------ Getters -------------------------------------------------------*/
    /* Pointwise observability for gridpoint (aI,aJ,aK) */
    virtual double getObservabilityGrid(const int aI, const int aJ, 
                                        const int aK) const override;

    /* Pointwise observability at physical location (aX,aY,aT) */
    virtual double getObservabilityPhysical(const double aX, const double aY, 
                                            const double aT) const override;

    /* Path followed by the observer */
    double getPosX(double aT) const override;
    double getPosY(double aT) const override;

  private:
    /** ------ Fields --------------------------------------------------------*/
    /* Number of gridpoints */
    int fN;
    int fNt;
    /* Physical size of domain */
    double fPhysMin;
    double fPhysMax;
    double fTime;

    /* The derivative (x'(t),y'(t)) of the trajectory */
    std::function<double(double)> fXder, fYder;
    /* Angle of visible sector for anisotropic observers */
    double fAlpha;

    /* TimeGrid object containing the visibility of each gridpoint */
    std::shared_ptr<CTimeGrid> fVisibility;
    /* TimeGrid which contains the signed distance to nearest obstacle */
    std::shared_ptr<CTimeGrid> fPhi;

    /** ------ Functions for computing visibility ----------------------------*/
    /* This is the main function for computing visibility */
    virtual void compute_visibility();
    /* Helper function that computes signed distance from the obstacles */
    virtual void compute_signed_distance();
    /* Given signed distance function, computes line of sight information */
    virtual void compute_line_of_sight();
};

/** ------ Inline definitions of functions -----------------------------------*/
/* Pointwise observability for gridpoint (aI,aJ,aK) */
inline double CMovingObserver::getObservabilityGrid(const int aI, const int aJ, 
                                                    const int aK) const {
  const double x = fVisibility->xGridToPhysical(aI);
  const double y = fVisibility->yGridToPhysical(aJ);
  const double t = fVisibility->tGridToPhysical(aK);
  if (fTerrain->isNotInObstacles(x,y,t)) {
    const double vis = fVisibility->getValue(aI,aJ,aK);
    const double obs = fObservability(getPosX(t), getPosY(t), x, y);
    return vis*obs + OBSERVABILITY_IN_SHADOWS;
  } else {
    return std::numeric_limits<double>::infinity();
  }
}

/* Path followed by the observer */
inline double CMovingObserver::getPosX(const double aT) const {
  return fXpath(aT);
}
inline double CMovingObserver::getPosY(const double aT) const {
  return fYpath(aT);
}

#endif
