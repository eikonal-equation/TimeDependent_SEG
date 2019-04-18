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
* File: CTimeDependentTracer.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This is the trajectory tracer class for time-dependent problems.
* Given the solution to the time-dependent HJB equation, this class computes the
* optimal trajectory. 
* It assumes a 2D regular grid with equal spacing in both directions.
* It uses trilinear interpolation of the HJB solution.
* The target and source are assumed to be single points.
* This implementation computes a grid search over a uniform grid of directions.
* It depends on the CTimeGrid class and is called by CAdversarialPlan.
*
* ==============================================================================
*/

#ifndef TIMEDEPENDENTTRACER_HPP
#define TIMEDEPENDENTTRACER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <vector>
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CTimeGrid.hpp"
#include "CTerrain.hpp"

/** ------ Data structures and typedefs --------------------------------------*/
/* Struct for a single Gridpoint */
struct Point2d {
  double x, y;
  Point2d(double a, double b) {
    x = a;
    y = b;
  }
};
/* Typedef for path objects */
typedef std::vector<Point2d> Path;

/** ------ Main Class definition ---------------------------------------------*/
class CTimeDependentTracer {
  public:
    /* ------ Constructors ---------------------------------------------------*/
    CTimeDependentTracer() = default;

    /**
     * The CTimeDependentTracer Constructor
     * This is the main constructor.
     * @param aGrid CTimeGrid object containing the value function of the PDE.
     * @param aSource physical coordinates of the source.
     * @param aTarget physical coordinates of the target.
     * @param aTerrain a shared pointer to the terrain object.
     * @param aTimeFactor the number of tracer steps per PDE timestep.
     */
    CTimeDependentTracer(std::shared_ptr<CTimeGrid> aGrid, 
                         const std::vector<double> aSource, 
                         const std::vector<double> aTarget, 
                         const std::shared_ptr<CTerrain> aTerrain, 
                         const int aTimeFactor);

    /*-- Main ----------------------------------------------------------------*/
    /**
    * A path printing function.
    * Computes the optimal path and prints it to file specified by aFilename
    * @param aFilename a string containing the name of the file.
    */
    void printOptimalPathToFile(const std::string aFilename = "optpath") const;

    /**
    * The helper function for actually computing the optimal path
    */
    Path getOptimalPath() const;

  private:
    /* The number of tracer steps per PDE timestep.
       tau = dt/fTimeFactor */
    int fTimeFactor;

    /* Source and Target X and Y coordinates */
    std::vector<double> fSource;
    std::vector<double> fTarget;
    
    /* Timegrid containing value function of the PDE */
    std::shared_ptr<CTimeGrid> fGrid;
    /* Shared pointer to terrain object */
    std::shared_ptr<CTerrain> fTerrain;

    /** ------ Helpers -------------------------------------------------------*/
    /* Computes physical distance from (aX,aY) to the target */
    double distance_to_target(const double aX, const double aY) const;
};

/** ------ Inline definition of function -------------------------------------*/
inline double CTimeDependentTracer::distance_to_target(const double aX, 
                                                       const double aY) const {
    return std::sqrt(std::pow(fTarget[0] - aX,2) + std::pow(fTarget[1] - aY,2));
}

#endif
