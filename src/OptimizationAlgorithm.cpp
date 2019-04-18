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
* File: OptimizationAlgorithm.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This file contains the implementation of the projected 
* subgradient descent method. 
* 
* This method is used to compute min_{lambda} -G(lambda),
* which is an equivalent problem to max_{lambda} G(lambda).
* Internally it solves the convex problem with subgradients, but it outputs
* supergradients to the command line to match the physical problem.
*
* ==============================================================================
*/
#include "OptimizationAlgorithm.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <iostream>
#include <algorithm>

/** ------ Project-specific header files -------------------------------------*/
#include "WriteToFile.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace Eigen;
using namespace std;

/*==============================================================================
  Projected Subgradient Method
==============================================================================*/
Eigen::VectorXd projectedSubgradientMethod(function<void(VectorXd, double&, VectorXd&)> aValueSubgradientFun, 
                                           function<void(VectorXd&)> aProjectionFun, 
                                           const VectorXd aX0, 
                                           const int aMaxIter, 
                                           const string aFilename) {
  const int dimension = aX0.size();
  const int max_iter = aMaxIter;

  cout << endl << "Beginning projected supergradient ascent method..." << endl;

  /* Project initial guess. */
  VectorXd current_x = aX0;
  aProjectionFun(current_x);
  cout << "Initial Lambda: ";
  io::printEigenVec(current_x);
  cout << endl;

  /* Calculate initial value */
  double current_value;
  VectorXd current_subgradient = VectorXd::Zero(dimension);
  aValueSubgradientFun(current_x, current_value, current_subgradient);

  /* Initialize best values */
  VectorXd best_x = current_x;
  double best_value = current_value;
  double stag_value = current_value;

  /* Set scaled initial step size */
  const double base_step_size = 1.0/current_subgradient.norm();
  double step_size = base_step_size;

  /* Number of steps without improvement in objective */
  int stagnated_steps = 0;
  /* Used to store the values for all iterations. */
  vector<double> allvalues;
  /* Used to store the subgradients for all iterations. */
  vector<double> allsubgradients;
  /* Used to store the subgradients for all iterations. */
  vector<double> alliterates;
  /*============= Main Loop ==================================================*/
  for (int iter = 0; iter < aMaxIter; ++iter) {
    /* The subgradient step */
    VectorXd temp_x = current_x - step_size*current_subgradient;

    /* Projection */
    aProjectionFun(temp_x);

    /* Get new subgradient and value. */
    cout << "Iteration: " << iter << " | ";
    aValueSubgradientFun(temp_x, current_value, current_subgradient);
    if (USE_SUM_SUBGRADIENTS) {
      current_value = 0;
      for (int i = 0; i < dimension; ++i) {
        current_value += temp_x[i]*current_subgradient[i];
      }
    }
    current_x = temp_x;
    cout << "                 Supergradient: ";
    io::printEigenVec(-current_subgradient);
    cout << " | ";

    /** Check for stagnation. 
     *  Don't start before STEPS_BEFORE_STAGNATION iterations 
     * (until we have seen enough decrease in the stepsize) */
    if ((current_value + TOLERANCE > stag_value) && (iter > STEPS_BEFORE_STAGNATION)) {
      stagnated_steps++;
    } else {
      stagnated_steps = 0;
      stag_value = current_value;
    }
    cout << "Stagnated Steps: " << stagnated_steps << " | Stepsize: "<< step_size << endl;

    /* Update best value seen so far */
    if (current_value <= best_value) {
      best_value = current_value;
      best_x = current_x;
    }

    /* If stagnated for too long, end iterations */
    if (stagnated_steps >= MAX_STAGNATED_STEPS) {
      cout << "Iteration has stagnated for " << stagnated_steps 
           << " steps, with tolerance of " << TOLERANCE
           << ". Ending iterations." << endl;
      break;
    }

    /* Print current iteration */
    if (aFilename.compare("") != 0) {
      allvalues.push_back(current_value);
      for (int i = 0; i < dimension; ++i) {
        allsubgradients.push_back(current_subgradient[i]);
        alliterates.push_back(current_x[i]);
      }
    }

    /* Update step size */
    if (SQRT_STEP_SIZE) {
      step_size = base_step_size/(sqrt(iter) + 1);
    } else {
      step_size = base_step_size/(iter + 1);
    }
  }

  /**--------------------------------------------------------------------------/
  * Save all iterations to file
  * --------------------------------------------------------------------------*/
  if (aFilename.compare("") != 0) {
    io::writeVectorToFile(aFilename+"IteratesValues", allvalues);
    io::writeVectorToFile(aFilename+"IteratesSubgradients", allsubgradients);
    io::writeVectorToFile(aFilename+"Iterates", alliterates);
  }
  return best_x;
}