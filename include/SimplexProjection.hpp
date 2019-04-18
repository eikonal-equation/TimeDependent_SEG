/*
* ==============================================================================
*
*  Copyright (C) 2019  Marc Aurèle Gilles
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
* File: SimplexProjection.hpp
*
* Authors: Marc Aurèle Gilles
*
* Description: Function that performs the simplex projection.
*
* ==============================================================================
*/

#ifndef SIMPLEXPROJECTION_HPP
#define SIMPLEXPROJECTION_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <vector>
#include <Eigen>

/** ----- Project-specific header files --------------------------------------*/
#include "GlobalConfiguration.hpp"

/**
 * A function which performs the orthogonal projection onto the probability 
 * simplex S = { x \in \mathcal{R}^n | x_i \geq 0, \sum_i x_i = 1 }
 * This function performs the operation:
 * x = argmin_{ sum(x_i) = 1, x_i >= 0  } \| x - \lambda \|_2
 * As described by Wang & Carreria-Perpinan, 2013.
 * @param aLambda Eigen::VectorXd of which is projected onto the probability 
 * simplex.
*/
void simplexProjection(Eigen::VectorXd& aLambda) {
  int nonzeros = 0;
  double cumsum=0;
  double threshold = 0;

  /** Cast to STL for sorting. */
  std::vector<double> SortedLambda(aLambda.data(), aLambda.data() + aLambda.size());
  std::sort(SortedLambda.begin(), SortedLambda.end(), std::greater<double>());

  for (int i = 0; i <  SortedLambda.size(); ++i) {
    cumsum += SortedLambda[i];
    if (SortedLambda[i] + 1.0/(i + 1)*(1 - cumsum) <= 0) {
      cumsum = cumsum - SortedLambda[i];
      break;
    }
    nonzeros++;
  }

  const double kappa = (1.0/(nonzeros))*(1 - cumsum);
  for (int i = 0; i < SortedLambda.size(); ++i) {
      aLambda[i] = std::max((double) aLambda[i] + kappa, (double) 0.0);
  }
}

#endif
