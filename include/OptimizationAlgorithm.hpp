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
* Description: This file contains the declaration of the projected subgradient
* descent method.
*
* ==============================================================================
*/

#ifndef OPTIMIZATIONALGORITHM_HPP
#define OPTIMIZATIONALGORITHM_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <vector>
#include <string>
#include <Eigen>

/** ----- Project-specific header files --------------------------------------*/
#include "GlobalConfiguration.hpp"

/**
 * A projected subgradient descent method.
 * It solves the optimization problem : min_{x \in S} f(x) by iterating:
 *
 * for k=0:maxIter
 *  x_{k+1} = Pi(x_k - s* \partial f x_{k} )
 * end
 *
 * where Pi is an orthogonal projection onto the set S 
 * and \partial f x is a subgradient of f. 
 * The stepsize s is by default 1/(norm(\partial f(aX0))*(k+1)).
 * This implementation is however specialized to our specific problem, 
 *      and it is reflected in a couple of ways:
 * Convergence is computed by a stagnation coefficient, 
 *      but only after a fixed number of iterations has been attempted. 
 * This is because we have a step rule where for the first few iterations, 
 *      we may be taking steps which are way too large.
 * @param aValueSubgradientFun 
 *    A function which returns the value and subgradient. 
 *    aValueSubgradientFun should take three arguments: 
 *        First argument is the input x, 
 *        Second argument is the output f(x) 
 *        Third argument is the output \partial f (x).
 * @param aProjectionFun 
 *    A function which performs the orthogonal projection onto the set S.
 *    aProjectionFun should have one argument that is both input and output.
 * @param aX0 a Eigen::VectorXd containing the initial guess of the method.
 * @param aMaxIter integer, the maximum number of iterations
 * @param aFilename string, prints the iterates to this field.
 *
 * This method is used to compute min_{lambda} -G(lambda),
 * which is an equivalent problem to max_{lambda} G(lambda).
 * Internally it solves the convex minimization problem with subgradients, 
 * but it actually prints supergradients to the command line.
*/
Eigen::VectorXd projectedSubgradientMethod(
    std::function<void(Eigen::VectorXd,double&,Eigen::VectorXd&)> aValueSubgradientFun,
    std::function<void(Eigen::VectorXd&)> aProjectionFun, Eigen::VectorXd aX0,
    const int aMaxiter = 100, const std::string aFilename = "");

#endif
