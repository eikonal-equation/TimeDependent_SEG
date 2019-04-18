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
*  along with this program. If not, see <http://www.gnu.org/licenses/>.
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
* File: CAdversarialPlan.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: Computes the Nash equilibrium of a surveillance-evasion game.
* First the optimal probability distribution for the observer is computed using 
* a projected subgradient descent method, where the subgradient is the vector of
* partial costs (integrated observability from each position) at the source and 
* the value of the function is the scalarized cost (expected integrated 
* observability) at the source. Then the optimal probability distribution is 
* perturbed to obtain a set of different paths which form a Nash equilibrium.
* It uses the class CFixedLambdaPlanner.
* (also: see CAdversarialPlan.hpp)
*
* ==============================================================================
*/
#include "CAdversarialPlan.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <functional>
#include <iostream>
#include <Eigen>
#include "QuadProg++.hh"

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "SimpleFunctions.hpp"
#include "WriteToFile.hpp"
#include "SimplexProjection.hpp"
#include "OptimizationAlgorithm.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace Eigen;
using namespace memory;
using namespace std;
using namespace std::placeholders;

/** ------ Constructor -------------------------------------------------------*/
CAdversarialPlan::CAdversarialPlan(shared_ptr<CFixedLambdaPlanner> aFixedLambdaPlan) {
  fFixedLambdaPlanner = aFixedLambdaPlan;
}

/** -- Find Nash equilibrium of SE game using projected subgradient descent --*/
void CAdversarialPlan::findNash(const string aFilename, const int aMaxIter) {
  /* Define Objective function & subgradient Projected Subgradient Descent */
  std::function<void(VectorXd, double&, VectorXd&)> val_and_subgrad = 
      std::bind(&CFixedLambdaPlanner::getValueAndSubgradient, fFixedLambdaPlanner, _1, _2, _3);

  /* Define Projection for Projected Subgradient Descent */
  std::function<void(VectorXd&)> projection = simplexProjection;

  /* Initial value of Lambda is uniform probability */
  const int num_observers = fFixedLambdaPlanner->getNumberOfObservers();
  VectorXd X0 = VectorXd::Constant(num_observers, 1.0 / num_observers);

  /* Find best strategy for observer */
  VectorXd Lambda = projectedSubgradientMethod(val_and_subgrad, projection, X0, aMaxIter, aFilename);

  /* Find best strategy for evader. */
  if (num_observers > 0) {
    findMultiplePaths(Lambda, aFilename);
  }

  /* Write results to file */
  fFixedLambdaPlanner->writeAllToFile(aFilename);
}

/** Find optimal trajectory (trajectories) corresponding to observer 
 *  probability distribution aLambda, and print the results to aFilename */
void CAdversarialPlan::findMultiplePaths(const VectorXd aLambda, const string aFilename) {
  cout << endl << "Finding optimal strategy for Evader..." << endl;
  /* Number of possible observer trajectories */
  const int num_observers = aLambda.size();
  
  /* Variables for holding value and subgradient */
  VectorXd lambda = aLambda;
  double value;

  /* Set small probabilities to zero. */
  int non_zeros = 0; 
  std::vector<bool> keep(num_observers);
  for (int k = 0; k < num_observers; ++ k) {
    if (lambda[k] < ZERO_THRESHOLD) {
      lambda[k] = 0;
      keep[k] = false;
    } else {
      keep[k] =  true;
      non_zeros += 1;
    }
  }
  VectorXd reduced_lambda(non_zeros);
  VectorXd direction(non_zeros);

  /* Maximum numer of optimal Evader trajectories */
  const int max_number_paths = non_zeros + 1;

  /* Project the vector with zero-ed components onto probability simplex. */
  get_short_vector(reduced_lambda, lambda, keep);
  simplexProjection(reduced_lambda);
  get_long_vector(reduced_lambda, lambda, keep);

  /* Get New Value and subgradient from tresholded probability.
   * Recall that this function really returns the negative of the value */
  VectorXd subgradient(num_observers);
  fFixedLambdaPlanner->getValueAndSubgradient(lambda, value, subgradient);
  cout << "Total cost of best trajectory found: " << -value << endl;

  /** Matrix of Partial costs: Used to perform the least square. */
  MatrixXd A = MatrixXd::Zero(non_zeros,max_number_paths);

  /** Right hand side: constant entry equal to Value. */
  MatrixXd rhs = MatrixXd::Constant(non_zeros, 1, lambda.transpose()*subgradient);

  /* Store probabilities and costs. */
  array2D_t<double> costs = memory::allocateArray2D<double>(max_number_paths, num_observers + 1);
  array2D_t<double> probabilities = memory::allocateArray2D<double>(max_number_paths, num_observers);

  /* Partial Costs */
  for (int i = 0; i < num_observers; ++i) {
    costs[0][i] = subgradient[i];
  }
  cout << "Partial Costs of best trajectory found: ";
  io::printEigenVec(-subgradient);
  cout << endl;

  /* In last entry of Costs matrix, store expected costs. */
  costs[0][num_observers] = value;

  /* Store Probability distribution */
  for (int i = 0; i < num_observers; ++ i) {
    probabilities[0][i] = lambda[i];
  }

  /* Store first strategy. */
  fFixedLambdaPlanner->saveStrategyToFile(aFilename + "path0");

  /*----------- LOOP STARTS HERE ---------------------------------------------*/
  int k;
  for (k = 1; k < max_number_paths; ++k) {
    /* Extract appropriate component of subgradients to form A matrix */
    int i = 0; 
    int j = 0;
    while(j < non_zeros) {
      if (keep[i]== true) {
        A(j,(k-1)) = subgradient(i);
        ++j;
      }
      ++i;
    }
    cout << "--------------- Partial Cost Matrix ---------------" << endl;
    cout << -A.block(0,0,non_zeros,k) << endl;
    cout << "--------------- Partial Cost Matrix ---------------" << endl;

    /* Solve least square problem */
    VectorXd theta = call_quad_prog(A.block(0,0,non_zeros,k),rhs);

    /* The perturbation direction is the direction of the subproblem. */
    direction = A.block(0,0,non_zeros,k)*theta - rhs;
    cout << "Norm of residual at step " << k << " is " 
         << direction.norm()/(fabs(value)*non_zeros) << endl;
    /* If residual is small, then we are done. */
    if (direction.norm()/(fabs(value)*non_zeros) < RESIDUAL_THRESHOLD) {
      cout << "Small residual. Terminating..." << endl;
      k++;
      break;
    }

    if (theta[k-1] < 1e-8) {
      cout << "WARNING: Last trajectory found was not used. Terminating... " << endl;
      break;
    }

    cout << "Partial cost matrix: " << endl << -A.block(0,0,non_zeros,k) << endl;
    cout << "Theta: ";
    io::printEigenVec(theta);
    cout << endl;
    cout << "Lambda: ";
    io::printEigenVec(lambda);
    cout << endl;

    VectorXd temp_lambda = lambda;

    double stepsize = INIT_STEP_SIZE;
    bool new_path_found = false;
    while ((new_path_found == false)  && (stepsize < 1)) {
      /* Perturb and project back onto probabilty simplex */
      get_short_vector(reduced_lambda, temp_lambda, keep);
      reduced_lambda = reduced_lambda - direction*stepsize;
      simplexProjection(reduced_lambda);
      get_long_vector(reduced_lambda, temp_lambda, keep);

      /* Add tiny values for zeros */
      for (int i = 0 ; i < subgradient.size(); ++i) {
        if (temp_lambda[i] == 0){
          temp_lambda[i] = ALMOST_ZERO_OBSERVABILITY;
        }
      }

      /* Compute new subgradient & value */
      fFixedLambdaPlanner->getValueAndSubgradient(temp_lambda, value, subgradient);
      cout << "Partial costs found: ";
      io::printEigenVec(-subgradient);
      cout << endl;

      /* Check if similar to previous subgradients */
      int num_paths = 0;
      for (int j = 0; j < k; ++j) {
          double diff = 0;
          for (int i = 0; i < subgradient.size(); ++i) {
            /* Compute L^2 distance in the vector of individual costs. */
            diff += (subgradient[i] - costs[j][i])*(subgradient[i] - costs[j][i]);
          }
          if (sqrt(diff) > 0.5*subgradient.norm()) {
            /* If the difference is greater than this threshold, 
             * we have found a new path. */
            num_paths++;
          }
      }
      /* If the current subgradient is different from all others, 
       * then we have a new path. */
      if (num_paths == k) {
        new_path_found = true;
      } else {
        cout << "New trajectory not found. Increasing perturbation size at step " << k << endl;
        stepsize = stepsize*2;
      }
    }

    cout << "G(lambda): " << -value << endl;

    /** Store Partial Costs */
    for (int i = 0; i < num_observers; ++i) {
      costs[k][i] = subgradient[i];
    }

    /** Store Expected cost */
    costs[k][num_observers] = value;

    /** Store Probabilities */
    for(int i = 0; i < num_observers; ++i) {
      probabilities[k][i] = temp_lambda[i];
    }

    /** Save current strategy to file */
    fFixedLambdaPlanner->saveStrategyToFile(aFilename + "path" + to_string(k));
  }

  /** Save to file */
  io::writeToFile2D<double>(aFilename + "Costs", costs);
  io::writeToFile2D<double>(aFilename + "Probabilities", probabilities);

  /** Recompute theta */
  VectorXd theta = call_quad_prog(A.block(0,0,non_zeros,k-1),rhs);
  cout << "Trajectory tracing complete." << endl;
  cout << "Theta : ";
  io::printEigenVec(theta);
  cout << endl;
  vector<double> theta_vec(theta.data(), theta.data() + theta.size());
  io::writeVectorToFile(aFilename + "Theta", theta_vec);

  /** Recompute once more with input aLambda.
   *  This way, writing to file uses the best possible lambda */
  fFixedLambdaPlanner->getValueAndSubgradient(aLambda, value, subgradient);
}

void CAdversarialPlan::get_short_vector(VectorXd& aShort, VectorXd aLong, 
                                        vector<bool> aKeep) {
  int j = 0;
  for (int k = 0; k < aLong.size(); ++k) {
    if (aKeep[k] == true) {
      aShort[j] = aLong[k];
      ++j;
    }
  }
}

void CAdversarialPlan::get_long_vector(VectorXd aShort, VectorXd& aLong, 
                                       std::vector<bool> aKeep) {
  int j = 0;
  for (int k = 0; k < aLong.size(); ++k) {
    if (aKeep[k] == true) {
      aLong[k] = aShort[j];
      ++j;
    }
  }
}

/** Compute the Pareto front. Only for the 2 observer case */
void CAdversarialPlan::getParetoFront(const std::string aFilename, 
                                      const int aMaxIter) {
  /* Check that we really only have two observers */
  assert(fFixedLambdaPlanner->getNumberOfObservers() == 2);

  cout << endl << "Computing Pareto Front..." << endl;

  /* Calculate Pareto front */
  vector<double> PF;
  for (int i = 0; i <= aMaxIter; ++i) {
    VectorXd temp_x(2);
    temp_x(0) = float(i)/aMaxIter;
    temp_x(1) = 1 - float(i)/aMaxIter;
    double temp_value;
    VectorXd subgradient(2);
    cout << "Pareto Front iteration " << i + 1;
    cout << " of " << aMaxIter + 1 << ":" << endl << "    ";

    fFixedLambdaPlanner->getValueAndSubgradient(temp_x, temp_value, subgradient);
    PF.push_back(subgradient[0]);
    PF.push_back(subgradient[1]);
  }
  /* Write Pareto front to file */
  io::writeVectorToFile(aFilename + "PF", PF);
}

/* Solve the Quadratic Program.
*  The problem is in the form:
*  min 0.5 * x G x + g0 x
*  s.t.
*      CE^T x + ce0 = 0
*      CI^T x + ci0 >= 0
*/
VectorXd CAdversarialPlan::call_quad_prog(MatrixXd aA, VectorXd aRhs) {
  const int dim = aA.cols();
  const int num_evaders = 1;
  MatrixXd G = aA.transpose()*aA;
  VectorXd g0 =  -1*aA.transpose()*aRhs;
  MatrixXd CE(dim,num_evaders);
  for (int i = 0; i < dim/num_evaders; i++) {
    CE.block(i*num_evaders,0,num_evaders,num_evaders) = MatrixXd::Identity(num_evaders,num_evaders);
  }
  VectorXd ce0 = VectorXd::Constant(num_evaders,-1);
  MatrixXd CI(dim,2*dim);
  CI.block(0,0,dim,dim) = MatrixXd::Identity(dim,dim);
  CI.block(0,dim,dim,dim) = -1*MatrixXd::Identity(dim,dim);
  VectorXd ci0(2*dim);
  ci0.head(dim) = VectorXd::Constant(dim,0);
  ci0.tail(dim) = VectorXd::Constant(dim,1);
  VectorXd theta(dim);
  double d = QuadProgPP::solve_quadprog(G,g0,CE,ce0,CI,ci0,theta);
  return theta;
}