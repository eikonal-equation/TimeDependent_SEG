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
* File: CTimeDependentHjbSolver.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*
* Description: This class is for solving the time-dependent HJB equation.
* It assumes that fGrid contains the speed and cost functions.
* It assumes a 2D regular grid with equal spacing in both directions.
* The target is assumed to be single points.
* Updates from both 5-point and 9-point stencils are possible, 
*     depends on setting in GlobalConfiguration.hpp
* It depends on the CTimeGrid class and is called by CAdversarialPlan.
*
* ==============================================================================
*/

#ifndef CTIMEDEPENDENTHJBSOLVER_HPP
#define CTIMEDEPENDENTHJBSOLVER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>

/** ------ Project-specific header files -------------------------------------*/
#include "CTimeGrid.hpp"
#include "GlobalConfiguration.hpp"
#include "CTerrain.hpp"

class CTimeDependentHjbSolver
{
  public:
    /** Constructor
    * @param aGrid  shared pointer to a 3D TimeGrid representing value function
    * @param aTarget  vector containing physical coordinates of the target
    * @param aTerrain  shared pointer to the terrain object
    */
    CTimeDependentHjbSolver(std::shared_ptr<CTimeGrid> aGrid, 
                            std::vector<double> aTarget, 
                            std::shared_ptr<CTerrain> aTerrain);

    /** This is the main function for solving the HJB equation */
    void solveHJB(); 

  private:
    /* Shared pointer to a TimeGrid object which contains 
     * - The approximation u[i,j] of the value function 
     * - The speed f(x,t) of travel in the domain
     * - The poinwise expected observability K_lambda 
     */
    std::shared_ptr<CTimeGrid> fGrid;

    /* Physical coordinates of the target */
    std::vector<double> fTarget;

    /* Shared pointer to the terrain object */
    std::shared_ptr<CTerrain> fTerrain;

    /* Helper function that computes HJB update using 5-point stencil */
    double update_5pt(const int aI, const int aJ, const int aK) const;

    /* Helper function that computes HJB update using 9-point stencil */
    double update_9pt(const int aI, const int aJ, const int aK) const;

    /* Helper function for evaluating gradient on 9-point stencil */
    void eval_gradient(const int i_diff, const int i_direct, const double v_diag, 
                       const double v_dir, double (&u)[2]) const;

    /* Computes physical distance from (aX, aY) to the target */
    double distance_to_target(const double aX, const double aY) const;
};

/* Computes physical distance from (aX, aY) to the target */
inline double CTimeDependentHjbSolver::distance_to_target(const double aX, 
                                                          const double aY) const {
  return std::sqrt(std::pow(aX-fTarget[0],2) + std::pow(aY-fTarget[1],2));
}

#endif