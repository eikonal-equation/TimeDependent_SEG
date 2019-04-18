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
* File: CTimeDependentPlanner.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*
* Description: This is an implementation of the CFixedLambdaPlanner class for
* time-dependent surveillance-evasion games. It should handle the time-dependent
* HJB solver and time-dependent path tracer to compute values and subgradients
* of the objective function.
*
* ==============================================================================
*/

#ifndef CTIMEDEPENDENTPLANNER_HPP
#define CTIMEDEPENDENTPLANNER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <vector>
#include <string>
#include <Eigen>

/** ------ Project-specific header files -------------------------------------*/
#include "CFixedLambdaPlanner.hpp"
#include "GlobalConfiguration.hpp"
#include "MemoryAllocations.hpp"
#include "CObserver.hpp"
#include "CTimeGrid.hpp"
#include "CTimeDependentHjbSolver.hpp"
#include "CTimeDependentTracer.hpp"
#include "CTerrain.hpp"

class CTimeDependentPlanner : public CFixedLambdaPlanner
{
  public:
    /** Default constructor and destructor */
    virtual ~CTimeDependentPlanner();
    CTimeDependentPlanner() = default;

    /**
     * The CTimeDependentPlanner Constructor.
     * This is the main constructor.
     * @param aObserverPointers a vector of shared pointers to observer objects
     * @param aSource a vector containing the location of the source.
     * @param aTarget a vector containing the location of the target.
     * @param aTerrain a pointer to a CTerrain object specifying the obstacles.
     * @param aObservability function double(double,double,double,double) 
     *      The (unobstructed) pointwise observability function K(x,y).
     *      The signature should have four inputs: (x,y,observer x, observer y).
     * @param aSpeed function double(double,double) 
     *      Defines the speed of the evader at all points in the domain.
     *      Signature should be v(x,y).
     * @param aN  Number of grid points in each space dimension. Defaults to 51.
     * @param aNt  Number of grid points in the time dimension. Defaults to 101.
     * @param aNobs  Number of observers
     * @param aTimeFactor Ratio of path tracer timesteps to HJB timesteps.
     * @param aMinPhys 
     *      Used to define the domain [aMinPhys,aMaxPhys]x[aMinPhys,aMaxPhys]. 
     *      Defaults to 0.
     * @param aMaxPhys 
     *      Used to define the domain [aMinPhys,aMaxPhys]x[aMinPhys,aMaxPhys]. 
     *      Defaults to 1.
     * @param aTime Used to define the time interval [0,aTime]. Defaults to 4.
     */
    CTimeDependentPlanner(const std::vector<std::shared_ptr<CObserver>> aObserverPointers, 
                          const std::vector<double> aSource,
                          const std::vector<double> aTarget,
                          const std::shared_ptr<CTerrain> aTerrain,
                          const std::function<double(double,double,double,double)> aObservability = &inverseSquared,
                          const std::function<double(double,double)> aSpeed = &constant2,
                          const int aN = 51, const int aNt = 101, 
                          const int aNobs = 2, const int aTimeFactor = 1,
                          const double aMinPhys = 0, const double aMaxPhys = 1, 
                          const double aTime = 4);

    /**
     * Computes the expected pointwise observability
     *    for a given probability distribution lambda of observer locations
     *    for a specific point on the grid.
     * @param aLambda A probability distribution over the observer positions.
     * @param aI, aJ, aK are the grid coordinates of the point in the domain.
     * NOTE: Uses grid coordinates, not physical coordinates!
     */
    virtual double getExpectedObservability(const Eigen::VectorXd aLambda, 
                                            const int aI, 
                                            const int aJ, 
                                            const int aK) const;

    /**
     * A helper function for the subgradient method.
     * For a given probability distribution over the observer locations aLambda,
     * this function computes the optimal expected cost 
     * (i.e. the solution of the scalarized HJB at the source) 
     * which is the value of the objective function, 
     * and the vector of partial costs at the source,
     * which is the subgradient of the objective function. 
     * It uses the CTimeDependentHjbSolver and CTimeDependentTracer objects. 
     * NOTE: It returns the negative of the quantities expected, 
     *    This is because subgradient descent is meant for a convex function, 
     *    while our objective function is concave.
     * @param aLambda A probability distribution over the observer positions.
     * @param aValue output, the value of the objective function.
     * @param aSubgradient output, the subgradient of the objective function.
     */
    virtual void getValueAndSubgradient(const Eigen::VectorXd aLambda, 
                                        double& aValue, 
                                        Eigen::VectorXd& aSubgradient) override;

    /**
     * I/O function which writes the current state of all variables to file.
     * @param aFilename 
     *    String containing the name of the file for which to save all the data.
     *    Defaults to the empty string.
     */
    virtual void writeAllToFile(const std::string aFilename = "") const override;

    /* Writes current path to file.
     * @param aFilename 
     *    String containing the name of the file for which to save all the data.
     *    Defaults to the empty string.
     */
    virtual void saveStrategyToFile(const std::string filename) const override;

    /** Returns the number of observers. */
    virtual int getNumberOfObservers() const override;

  private:
    /* Size of the grid. */
    int fN;
    /* Number of timesteps. */
    int fNt;
    /* Number of sub time slices for path tracing. */
    int fTimeFactor;
    /* Number of observer trajectories. */
    int fNumObservers;
    /* Total physical time. */
    double fTime;
    /* Vector of pointers to CObserver objects. */
    std::vector<std::shared_ptr<CObserver>> fObserverPointers;
    /* Pointer to CTimeGrid with the properties of the HJB equation. */
    std::shared_ptr<CTimeGrid> fExpectedObsGrid;
    /* Pointer to CHjbSolver Object, responsible for solving the HJB PDE. */
    std::shared_ptr<CTimeDependentHjbSolver> fHjbSolver;
    /* Vector of length 2 containing the physical coordinates of the target */
    std::vector<double> fTarget;
    /* Vector of length 2 containing the physical coordinates of the source */
    std::vector<double> fSource;
    /* Pointer to fTracer object which computes the optimal trajectory 
     * through gradient descent on the value function of the scalarized HJB 
     * (contained in fExpectedObsGrid) */
    std::shared_ptr<CTimeDependentTracer> fTracer; 
    /* Pointer to Terrain object which owns the obstacles */
    std::shared_ptr<CTerrain> fTerrain;
};

#endif
