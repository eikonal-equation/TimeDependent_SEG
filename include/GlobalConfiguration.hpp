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
* File: GlobalConfiguration.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This file contains the definitions of various numerical constants
* and types. It also contains all global parameters for all of the numerics.
*
* ==============================================================================
*/

#ifndef GLOBAL_CONFIG_HPP
#define GLOBAL_CONFIG_HPP

/*-- STL ---------------------------------------------------------------------*/
#include <cassert>
#include <functional>
#include <limits>

/*-----------------------------------------------------------------------------/
/-- Numerical Constants
/-----------------------------------------------------------------------------*/
constexpr double INF = std::numeric_limits<double>::max();
constexpr double PI = 3.141592653589793;
constexpr double SQRT2 = 1.41421356237309504880168872420969807;

/*-----------------------------------------------------------------------------/
/-- Gridpoint status for FMM
/-----------------------------------------------------------------------------*/
enum status_t {FAR, CONSIDERED, ACCEPTED};

/* ------ Shorter version of function<double(double)> ------------------------*/
typedef double (*Function1D)(double);

/*-----------------------------------------------------------------------------/
/-- Global parameters & flags
/-----------------------------------------------------------------------------*/
/** ------ Constants for time-dependent HJB solver ---------------------------*/
/* Flag determining whether or not to use Lagrangian updates near the target */
constexpr bool LAGRANGIAN_NEAR_TARGET = true;
/* Flag determining whether to use the nine point stencil or not */
constexpr bool NINE_POINT_STENCIL = true;
/* (physical) distance in which to compare with Lagrangian updates */
constexpr double LAGRANGIAN_THRESHOLD = 0.02;
/* Terminal conditions at all non-target points */
constexpr double TERMINAL_VALUE = 10;
/* Placeholder for value inside the obstacles */
constexpr double VALUE_INSIDE_OBSTACLES = 1e8;

/** ------ Constants for observability computation ---------------------------*/
/* Pointwise observability in shadow zones (referred to as sigma in paper) */
constexpr double OBSERVABILITY_IN_SHADOWS = 0.1;

/** ------ Constants for subgradient computation -----------------------------*/
/* Flag determining whether to determine pointwise observability 
* at non-grid locations using rounding or interpolation.
* (Note: only affects subgradient computations, not HJB solver or tracer) */
constexpr bool INTERPOLATE_OBSERVABILITY = true;

/** ------ Constants for outer loop ------------------------------------------*/
/* Flag for scaling outer loop step sizes as 1/sqrt(n_iter) instead */
constexpr bool SQRT_STEP_SIZE = false;
/* Flag for using sum of subgradients instead for the outer loop */
constexpr bool USE_SUM_SUBGRADIENTS = false;
/* Steps of outer loop before tracking stagnation */
constexpr int STEPS_BEFORE_STAGNATION = 30;
/* Maximum number of stagnated steps before terminating outer loop */
constexpr int MAX_STAGNATED_STEPS = 100;
/* Tolerance in outer loop iterations */
constexpr double TOLERANCE = 1e-6;

/** ------ Parameters for multiple path algorithm. Fixed by experience. ------*/
/* The initial size of the perturbation. */
constexpr double INIT_STEP_SIZE = 1e-4; 
/* Probabilities below this threshold are set to 0. */
constexpr double ZERO_THRESHOLD = 5e-3; 
/* How close the characterization based on the residual is to 0. */
constexpr double RESIDUAL_THRESHOLD = 3e-3; 
/* Parameter which gives a little observability to observer positions which 
 * should have none at all (because their probability was set to 0). 
 * This is essentially for tiebreaking. 
 * Can be set to zero to turn this option off. */
constexpr double ALMOST_ZERO_OBSERVABILITY = 1e-6; 

/** ------ Parameters for trajectory tracing ---------------------------------*/
/* N_THETA is the number of directions to search in.
 * Should be divisible by 4 so that the 4 grid direction are searched. */
  const int N_THETA = 360;

#endif
