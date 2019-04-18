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
* File: CTimeDependentPlanner.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*
* Description: This is an implementation of the CFixedLambdaPlanner class for
* time-dependent surveillance-evasion games. It should handle the time-dependent
* HJB solver and time-dependent path tracer to compute values and subgradients
* of the objective function.
* (See also CTimeDependentPlanner.cpp)
*
* Internally it computes -G(lambda), but will actually print G(lambda) to the 
* command line to match the physical problem.
*
* ==============================================================================
*/
#include "CTimeDependentPlanner.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <iostream>
#include <cmath>
#include <functional>
#include <Eigen>

/** ------ Project-specific header files -------------------------------------*/
#include "SimpleFunctions.hpp"
#include "WriteToFile.hpp"
#include "OptimizationAlgorithm.hpp"
#include "CObserver.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace Eigen;
using namespace memory;
using namespace std;
using namespace std::placeholders;

/*==============================================================================
  Destructor
==============================================================================*/
CTimeDependentPlanner::~CTimeDependentPlanner() {}

/*==============================================================================
  Constructor
==============================================================================*/
CTimeDependentPlanner::CTimeDependentPlanner(const vector<shared_ptr<CObserver>> aObserverPointers, 
                                             const vector<double> aSource, 
                                             const std::vector<double> aTarget, 
                                             const shared_ptr<CTerrain> aTerrain, 
                                             const function<double(double,double,double,double)> aObservability, 
                                             const std::function<double(double, double)> aSpeed,  
                                             const int aN, const int aNt, 
                                             const int aNobs, const int aTimeFactor, 
                                             const double aMinPhys, 
                                             const double aMaxPhys, 
                                             const double aTime) {
  /* Copy arguments */
  fN = aN;
  fNt = aNt;
  fTime = aTime;

  /* Initialize source,target and observer locations. */
  fSource = aSource;
  fTarget = aTarget;
  fNumObservers = aObserverPointers.size();
  fTimeFactor = aTimeFactor;
  fObserverPointers = aObserverPointers;
  fTerrain = aTerrain;

  /* Initialize the Expected value grid. */
  const std::function<double(double, double)> notInsideObstacles = 
      std::bind(&CTerrain::isNotInObstacles, fTerrain, _1, _2, 0);
  fExpectedObsGrid = make_shared<CTimeGrid>(&constantcost3, notInsideObstacles, fN, fNt, aMinPhys, aMaxPhys, aTime);

  /* Initialize HJB and fTracer Objects. */
  fHjbSolver = make_shared<CTimeDependentHjbSolver>(fExpectedObsGrid, fTarget, fTerrain);
  fTracer = make_shared<CTimeDependentTracer>(fExpectedObsGrid, fSource, fTarget, fTerrain, fTimeFactor);
}

/*==============================================================================
  Pointwise expected observability function
==============================================================================*/
double CTimeDependentPlanner::getExpectedObservability(const VectorXd aLambda, 
                                                       const int aI, 
                                                       const int aJ, 
                                                       const int aK) const {
  double cost = 0;
  for (int k = 0; k < fNumObservers; ++k) {
    cost += fObserverPointers[k]->getObservabilityGrid(aI,aJ,aK)*aLambda[k];
  }
  return cost;
}

/*=======================================================================================================================
  Write everything to file
=======================================================================================================================*/
void CTimeDependentPlanner::writeAllToFile(const string aFilename) const {
  /* Write Expected visibility grid to file. */
  fExpectedObsGrid->writeGridToFile(aFilename);

  /* Write Path to File */
  fTracer->printOptimalPathToFile(aFilename + "OptPath");

  /* This file will contain source and target locations */
  vector<double> locations = {fSource[0], fSource[1], fTarget[0], fTarget[1]};
  io::writeVectorToFile<double>(aFilename + "locations", locations);

  /* This file will contain the number of observers, nx, nt and time_factor */
  vector<int> configs = {fNumObservers, fNt, fN, fTimeFactor};
  io::writeVectorToFile<int>(aFilename + "Configs", configs);

  /* Write observer paths to file */
  for(int k = 0; k < fNumObservers; ++k){
    memory::array2D_t<double> observer_locations = memory::allocateArray2D<double>(fNt,2);
    for(int i=0; i < fNt; ++i) {
      const double t = fExpectedObsGrid->tGridToPhysical(i);
      observer_locations[i][0] = fObserverPointers[k]->getPosX(t);
      observer_locations[i][1] = fObserverPointers[k]->getPosY(t);
    }
    io::writeToFile2D<double>(aFilename + "Observer" + std::to_string(k+1) + "locations", observer_locations);
  }
}

/*==============================================================================
  Save current path to file
==============================================================================*/
void CTimeDependentPlanner::saveStrategyToFile(const std::string aFilename) const {
  fTracer->printOptimalPathToFile(aFilename);
}

/** Returns the number of observers. */
int CTimeDependentPlanner::getNumberOfObservers() const {
  return fNumObservers;
}


/*==============================================================================
  Helper function for projected subgradient method. 
  Computes value and subgradient of G(lambda)
==============================================================================*/
void CTimeDependentPlanner::getValueAndSubgradient(const VectorXd aLambda, 
                                                   double& aValue, 
                                                   VectorXd& aSubgradient) {
  /* Get source location */
  const int xSource = fExpectedObsGrid->xPhysicalToGrid(fSource[0]);
  const int ySource = fExpectedObsGrid->yPhysicalToGrid(fSource[1]);

  /* Update Cost function. */
  if (fNumObservers > 0) {
    std::function<double(int, int, int)> CombinedCost = 
        std::bind(&CTimeDependentPlanner::getExpectedObservability, this, aLambda, _1, _2, _3);
    fExpectedObsGrid->setCost(CombinedCost);
  } else {
    fExpectedObsGrid->setCost(&constantcost3);
  }

  /* Run Time-dependent HJB Solver. */
  fHjbSolver->solveHJB();

  /* Value is the solution of time-dependent HJB evaluated at the source. */
  aValue = -fExpectedObsGrid->getValue(xSource, ySource, 0);

  /* Trace optimal path */
  Path path = fTracer->getOptimalPath();
  aSubgradient = VectorXd::Constant(fNumObservers,0);

  

  /* Compute subgradient along path */
  const double tau = fExpectedObsGrid->getDt() / fTimeFactor; 
  for (int i = 0; i < path.size() - 1; i++) {
    double timestep = tau;
    if (i == path.size() - 2) {
      /* Special case for end of path */
      const double speed = fExpectedObsGrid->getSpeedPhysical(fTarget[0],fTarget[1]);
      const double dist = sqrt(pow(path[i].x - path[i+1].x,2) + pow(path[i].y - path[i+1].y,2));
      timestep = dist/speed;
    }
    for (int k = 0; k < fNumObservers; k++) {
      const double observability = fObserverPointers[k]->getObservabilityPhysical(path[i].x, path[i].y, (i+1)*tau);
      const double cost_tmp = timestep * observability;
      aSubgradient[k] -= cost_tmp;
    }
  }
  path.clear();

  /* Print value and weighted sum of subgradients */
  double sum_costs = 0;
  for (int k = 0; k < fNumObservers; ++k) {
    sum_costs += aLambda[k] * aSubgradient[k];
  }
  cout << "Lambda: ";
  io::printEigenVec(aLambda);
  cout << " | G(Lambda): " << -aValue << endl; /* Use -aValue to get G(lambda) */
}