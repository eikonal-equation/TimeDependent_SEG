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
* File: CTimeDependentTracer.cpp
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
* (see also CTimeDependentTracer.hpp)
*
* ==============================================================================
*/
#include "CTimeDependentTracer.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath> 
#include <fstream>
#include <iostream>
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "CTimeGrid.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;

/*==============================================================================
  Constructor
==============================================================================*/
CTimeDependentTracer::CTimeDependentTracer(shared_ptr<CTimeGrid> aGrid, 
                                           const vector<double> aSource, 
                                           const vector<double> aTarget, 
                                           const shared_ptr<CTerrain> aTerrain, 
                                           const int aTimeFactor) {
  /* Copy arguments */
  fGrid = aGrid;
  fSource = aSource;
  fTarget = aTarget;
  fTerrain = aTerrain;
  fTimeFactor = aTimeFactor;
}

/*==============================================================================
  Compute optimal path
==============================================================================*/
Path CTimeDependentTracer::getOptimalPath() const {
  /* Pseudo-timestep for path computation */
  const double tau = fGrid->getDt() / fTimeFactor;

  /* Add source point to path */
  Path path;
  path.push_back(Point2d(fSource[0], fSource[1]));
  double dist_to_target = distance_to_target(fSource[0], fSource[1]);
  double speed = fGrid->getSpeedPhysical(path[0].x, path[0].y);
  int n_steps = 1;

  /** ======== Main Loop =====================================================*/
  while (dist_to_target > tau * speed) {
    const double curr_x = path[path.size() - 1].x;
    const double curr_y = path[path.size() - 1].y;
    double min_val = INF;
    double best_x = -1;
    double best_y = -1;
    double best_theta = -1;

    /* Grid search over directions */
    for (int i = 0; i < N_THETA; ++i) {
      const double theta = 2 * PI * (i) / N_THETA;
      const double new_x = curr_x + cos(theta) * tau * speed;
      const double new_y = curr_y + sin(theta) * tau * speed;
      const double new_t = n_steps * tau;

      double temp_val;
      if (fTerrain->isNotInObstacles(new_x,new_y,new_t) && fGrid->getMinX() <= new_x && new_x <= fGrid->getMaxX() && fGrid->getMinY() <= new_y && new_y <= fGrid->getMaxY()) {
        temp_val = fGrid->getValuePhysical(new_x, new_y, n_steps*tau);
      } else {
        temp_val = INF;
      }

      if(temp_val < min_val) {
        best_theta = theta;
        min_val = temp_val;
      }
    }
    best_x = curr_x + cos(best_theta) * tau * speed;
    best_y = curr_y + sin(best_theta) * tau * speed;

    /* Check if it is better to stay at current point */
    const double value_at_curr = fGrid->getValuePhysical(curr_x, curr_y, n_steps*tau);
    if (min_val >= value_at_curr) {
      best_x = curr_x;
      best_y = curr_y;
    }

    /* Add new point to path */
    path.push_back(Point2d(best_x, best_y));
    dist_to_target = distance_to_target(best_x, best_y);
    speed = fGrid->getSpeedPhysical(best_x, best_y);
    n_steps++;
  }

  /* Add target as final point of path. */
  path.push_back(Point2d(fTarget[0], fTarget[1]));
  return path;
}

/*==============================================================================
  Compute optimal path and write to file
==============================================================================*/
void CTimeDependentTracer::printOptimalPathToFile(const string aFilename) const {
  ofstream out("output/" + aFilename, ios::binary);
  Path path = getOptimalPath();
  for (int i = 0; i < path.size(); ++i) {
    out.write((char*) &path[i].x, sizeof(double));
    out.write((char*) &path[i].y, sizeof(double));
  }
}