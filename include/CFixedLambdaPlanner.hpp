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
* File: CFixedLambdaPlanner.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This is the abstract evader planner class, meant as an interface
* for CAdversarialPlan.hpp. The main function that children of this class should
* implement is getValueAndSubgradient, which is called in the
* subgradient/perturbation methods of CAdversarialPlan.
*
* ==============================================================================
*/

#ifndef CFIXEDLAMBDAPLANNER_HPP
#define CFIXEDLAMBDAPLANNER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <Eigen>
#include <memory>
#include "boost/multi_array.hpp"

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "MemoryAllocations.hpp"

class CFixedLambdaPlanner{
  public:
    /** Destructor */
    virtual ~CFixedLambdaPlanner() {};

    /* Computes the value and subgradient of the objective function.
    * This is the main function to be implemented by children.
    * @param aLambda (input), the probability of each observer location.
    * @param aValue (output), the value of the objective function at lambda
    * @param aSubgradient (output), the subgradient of the objective function.
    */
    virtual void getValueAndSubgradient(const Eigen::VectorXd aLambda, 
                                        double& aValue, 
                                        Eigen::VectorXd& aSubgradient) = 0;
    
    /* Returns the number of possible observer patrol trajectories. */
    virtual int getNumberOfObservers() const = 0;

    /* This saves the current strategy to file. 
     * Useful since this object is continually reused in CAdversarialPlan. */
    virtual void saveStrategyToFile(const std::string filename) const = 0;

    /* Saves the state of all variables to file */
    virtual void writeAllToFile(const std::string filename) const = 0;
};

#endif
