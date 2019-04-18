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
* File: CAdversarialPlan.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: Computes the Nash equilibrium of a surveillance-evasion game.
* First the probability distribution is computed using a projected subgradient
* descent method, where the subgradient is the vector of partial costs
* (integrated observability from each position) at the source and the value of
* the function is the scalarized cost (expected integrated observability) at the
* source. Then the optimal probability distribution is perturbed to obtain a set
* of different paths which form a Nash equilibrium.
* It uses the class CFixedLambdaPlanner.
*
* ==============================================================================
*/

#ifndef CADVERSARIALPLAN_HPP
#define CADVERSARIALPLAN_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <Eigen>

/** ----- Project-specific header files --------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CFixedLambdaPlanner.hpp"

class CAdversarialPlan
{
  public:
    CAdversarialPlan() = default;

    /** The CAdversarialPlan Constructor. 
    * @param aFixedLambdaPlan, the fixed lambda planner object */
    CAdversarialPlan(std::shared_ptr<CFixedLambdaPlanner> aFixedLambdaPlan);

    /**
     * This computes the Nash Equilibrium for the whole problem.
     * It finds the optimal probability distribution for the observer
     *    via the project subgradient descent method.
     * Then it uses path tracing to find the optimal strategy for the evader.
     * @param aFilename 
     *    String containing the name of the file for which to save all the data.
     *    Defaults to the empty string.
     * @param aMaxIter 
     *    Maximum number of subgradient method iterations. Defaults to 300.
     */
    void findNash(const std::string aFilename = "", const int aMaxIter = 300);

    /**
     * This function computes the optimal strategy for the observer.
     * This function perturbs the given probability distribution aLambda, 
     *    (which should be the optimal strategy for the observer) 
     *    in order to find the optimal path(s) associated with aLambda.
     * @param aFilename 
     *    String containing the name of the file for which to save all the data.
     *    Defaults to the empty string.
     */
    void findMultiplePaths(const Eigen::VectorXd aLambda, 
                           const std::string aFilename = "");

    /** Computes the convex part of the pareto front using the scalarization.
    * NOTE: Only works for 2 observers.
    * @param aFilename 
    *    String containing the name of the file for which to save all the data.
    *    Defaults to the empty string.
    * @param aGridSize integer, contains the number of lambda values used.
    */
    void getParetoFront(const std::string aFilename = "", 
                        const int aGridSize = 100);


  private:
    /** The CFixedLambdaPlanner object that is used */
    std::shared_ptr<CFixedLambdaPlanner> fFixedLambdaPlanner;

    /**
     * A helper function for the findMultiplePaths function.
     * Given a long Eigen::VectorXd, and a std::vector of boolean 
     *    indicating which indices of the long vector should be kept, 
     *    it returns the corresponding shortened Eigen::VectorXd aShort.
     * @param aShort output, 
     *    Eigen:vector of length M containing the entries in aLong 
     *    which have corresponding "1" entry in aKeep, where M = sum(aKeep).
     * @param aLong input, Eigen:vector of length N.
     * @param aKeep input, std::vector<bool> of length N with M true/1 entries.
     */
    static void get_short_vector(Eigen::VectorXd& aShort, Eigen::VectorXd aLong, 
                               std::vector<bool> aKeep);

    /**
     * Inverse operation of GetShortVector.
     */
    static void get_long_vector(Eigen::VectorXd aShort, Eigen::VectorXd& aLong, 
                              std::vector<bool> aKeep);

    /**
     * Uses the Library QuadProg++ to solve a quadratic program.
     * @param aA, Matrix of Partial costs: Used to perform the least square.
     * @param aRhs, Right hand side: constant entry equal to Value.
     */
    Eigen::VectorXd call_quad_prog(Eigen::MatrixXd aA, Eigen::VectorXd aRhs);
};

#endif
