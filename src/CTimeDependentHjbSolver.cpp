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
* File: CTimeDependentHjbSolver.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (Initial implementation of 8pt solver written by Tristan Reynoso)
*
* Description: This class is for solving the time-dependent HJB equation.
* It assumes that fGrid contains the speed and cost functions.
* It assumes a 2D regular grid with equal spacing in both directions.
* The target is assumed to be single points.
* Updates from both 5-point and 9-point stencils are possible, 
*     depends on setting in GlobalConfiguration.hpp
* It depends on the CTimeGrid class and is called by CAdversarialPlan.
* (See also CTimeDependentHjbSolver.hpp)
*
* ==============================================================================
*/
#include "CTimeDependentHjbSolver.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath> 
#include <string>
#include <iostream>
#include <algorithm>

/** ------ Project-specific header files -------------------------------------*/
#include "CTimeGrid.hpp"
#include "GlobalConfiguration.hpp"
#include "CTerrain.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace memory;
using namespace std;

/*==============================================================================
  Constructor
==============================================================================*/
CTimeDependentHjbSolver::CTimeDependentHjbSolver(std::shared_ptr<CTimeGrid> aValue, 
                                                 std::vector<double> aTarget, 
                                                 std::shared_ptr<CTerrain> aTerrain) {
  fGrid = aValue;
  fTarget = aTarget;
  fTerrain = aTerrain;
}

/*==============================================================================
  Main function for solving time-dependent HJB equation
==============================================================================*/
void CTimeDependentHjbSolver::solveHJB() {
  /* Get grid dimensions */
  const int nx = fGrid->getGridSizeX();
  const int ny = fGrid->getGridSizeY();
  const int nt = fGrid->getGridSizeT();
  const double h  = fGrid->getH();
  const double dt = fGrid->getDt();

  /* Enforce terminal conditions on value function */
  for (int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      fGrid->setValue(i,j,nt-1,TERMINAL_VALUE);
    }
  }

  /* Enforce value of zero at target */
  const int i_target = fGrid->xPhysicalToGrid(fTarget[0]);
  const int j_target = fGrid->yPhysicalToGrid(fTarget[1]);
  for (int k = 0; k < nt; ++k) {
    fGrid->setValue(i_target,j_target,k,0);
  }

  /* Evolve value function backward in time */
  for(int k = nt - 1; k > 0; --k) {
    for(int i = 0; i < nx; ++i) {
      for(int j = 0; j < ny; ++j) {
        const double t = fGrid->tGridToPhysical(k);
        const double x = fGrid->xGridToPhysical(i);
        const double y = fGrid->yGridToPhysical(j);
        double u_new;
        if(i == i_target && j == j_target) {
          u_new = 0;
        } else if(fTerrain->isNotInObstacles(x,y)) {
          if (NINE_POINT_STENCIL) {
            u_new = update_9pt(i,j,k);
          } else {
            u_new = update_5pt(i,j,k);
          }

          /* Compute Lagrangian update and compare to PDE value near target */
          const double dist_to_target = distance_to_target(x,y);
          if (LAGRANGIAN_NEAR_TARGET && (dist_to_target <= LAGRANGIAN_THRESHOLD)) {
            const double K = fGrid->getCost(i,j,k);
            const double f = fGrid->getSpeed(i,j);
            const double u_lagrange = K*dist_to_target/f;
            if (u_lagrange < u_new) {
              u_new = u_lagrange;
            } 
          }
        } else {
          u_new = VALUE_INSIDE_OBSTACLES;
        }
        fGrid->setValue(i,j,k-1,u_new);
      }
    }
  }
  return;
}

/*==============================================================================
  Nine-point stencil HJB update 
==============================================================================*/
double CTimeDependentHjbSolver::update_9pt(const int aI, const int aJ, 
                                           const int aK) const {
  /* Grid sizes */
  const int nx = fGrid->getGridSizeX();
  const int ny = fGrid->getGridSizeY();
  const int nt = fGrid->getGridSizeT();

  /* Get local variables */
  const double x = fGrid->xGridToPhysical(aI);
  const double y = fGrid->yGridToPhysical(aJ);
  const double t = fGrid->tGridToPhysical(aK);
  const double K = fGrid->getCost(aI,aJ,aK);
  const double f = fGrid->getSpeed(aI,aJ);
  const double u = fGrid->getValue(aI,aJ,aK);
  const double h  = fGrid->getH();
  const double dt = fGrid->getDt();

  /* Array encoding nine point stencil */
  const int stencil[8][2] = {{1,0}, {1,1}, {0,1}, {-1,1}, {-1,0}, {-1,-1}, {0,-1}, {1,-1}};
  
  double u_x, u_y;
  /** Calculate upwind discretization of gradient (8 neighbors) 
   *  (See Diagram below) */
  /*
     d[3]    d[2]   d[1]
          o----o----o
          |\   |   /|
          | \  |  / |
          |  \ | /  |
     d[4] o----u----o d[0]
          |  / | \  |
          | /  |  \ |
          |/   |   \|
          o----o----o
     d[5]   d[6]   d[7]
                             */

  /* Calculate the difference in value with all 8 neighbors */
  double dneighbors[8];
  for (int it = 0; it < 8; ++it) {
    const int i1 = stencil[it][0];
    const int j1 = stencil[it][1];
    if (fTerrain->isNotInObstacles(x + i1*h, y + j1*h) && aI +i1 >= 0 && aI + i1 <= nx-1 && aJ + j1 >= 0 && aJ + j1 <= ny-1) {
      dneighbors[it]  = (u - fGrid->getValue(aI + stencil[it][0], aJ + stencil[it][1], aK))/h;
    } else {
      dneighbors[it] = 0;
    }
  }

  /** i_direct, i_diagonal = index of the smallest direct/diagonal neighbor
   *   i_direct ∈ {0, 2, 4, 6}; i_diagonal ∈ {1, 3, 5, 7} */
  int i_direct = 0, i_diagonal = 1;
  double min_direct = dneighbors[0], min_diagonal = dneighbors[1];
  for (int i = 1; i < 4; ++i) {
    if (dneighbors[i*2] > min_direct) {
      i_direct = i*2;
      min_direct = dneighbors[i*2];
    }
    if (dneighbors[i*2 + 1] > min_diagonal) {
      i_diagonal = i*2 + 1;
      min_diagonal = dneighbors[i*2 + 1];
    }
  }

  /* i_diff indicates if the direct and diagonal neighbors chosen are adjacent */
  /* i_diff also helps encode (+)(-) relations between u_x and u_y in different simplices */
  int i_diff = i_diagonal - i_direct;
  if (i_diff == 7) {
    i_diff = -1;
  }
 
  /* If i_diff equals 1 or -1, then the direct and diagonal neighbors chosen are adjacent */
  if (i_diff * i_diff == 1) {
    double grad_u[2];
    eval_gradient(i_diff, i_direct, min_diagonal, min_direct, grad_u);
    u_x = grad_u[0];
    u_y = grad_u[1];
  } else {
    /* The case where the direct and diagonal neighbors chosen are not adjacent */
    double nom_dir[2], nom_diag[2];
    int diff_dir, diff_diag;

    /* Diagonal neighbor nominates a simplex */
    if (dneighbors[(i_diagonal+1)%8]>=dneighbors[i_diagonal-1]) {
      diff_diag = -1;
      eval_gradient(-1, (i_diagonal+1)%8, min_diagonal, dneighbors[(i_diagonal+1)%8], nom_diag);
    } else {
      diff_diag = 1;
      eval_gradient(1, i_diagonal-1, min_diagonal, dneighbors[i_diagonal-1], nom_diag);
    }

    /* Direct neighbor nominates a simplex */
    if (dneighbors[i_direct+1]>=dneighbors[(i_direct-1)%8]) {
      diff_dir = 1;
      eval_gradient(1, i_direct, dneighbors[i_direct+1], min_direct, nom_dir);
    } else {
      diff_dir = -1;
      eval_gradient(-1, i_direct, dneighbors[(i_direct-1)%8], min_direct, nom_dir);
    }

    /* Compare the two candidate simplices nominated by direct and diagonal neighbors */
    if (nom_dir[0] * nom_dir[0] + nom_dir[1] * nom_dir[1] >= nom_diag[0] * nom_diag[0] + nom_diag[1] * nom_diag[1]) {
      u_x = nom_dir[0];
      u_y = nom_dir[1];
      i_diagonal = i_direct + diff_dir;
    } else {
      u_x = nom_diag[0];
      u_y = nom_diag[1];
      i_direct = i_diagonal - diff_diag;
    }
  }

  double u_new;
  /* Check upwinding condition */
  if (((i_direct == 0 || i_direct == 4) && abs(u_y) > abs(u_x)) || ((i_direct == 2 || i_direct == 6) && abs(u_y) < abs(u_x))) {
    /* Upwinding condition is not satisfied, use smaller of one-sided updates */
    const double u_new1 = K*dt + dt*f/h*(u - dneighbors[i_direct]*h) + (1 - dt*f/h)*u;
    const double u_new2 = K*dt + dt*f/(h*sqrt(2))*(u - dneighbors[i_diagonal]*h) + (1 - dt*f/(h*sqrt(2)))*u;
    if (u_new1 >= u_new2){
      u_new = u_new2;
    } else {
      u_new = u_new1;
    }
  } else {
    /* Upwinding condition is satisfied, use Euler step */
    u_new = u + K*dt - f*sqrt(u_x*u_x + u_y*u_y)*dt;
  }
  return u_new;
}

/*==============================================================================
  Helper function for computing gradient using nine-point stencil
================================================================================         

               d[3]    d[2]   d[1]
                    o----o----o
                    |\ 3 | 2 /|
                    | \  |  / |
                    | 4\ | /1 |
               d[4] o----u----o d[0]
                    | 5/ | \8 |
                    | /  |  \ |
                    |/ 6 | 7 \|
                    o----o----o
               d[5]   d[6]   d[7]
        
        Compute the gradients in 8 simplices as follows:

        Simplex 1:  u_x = - d[0];
                    u_y = - d[1] - u_x;

        Simplex 2:  u_y = - d[2];
                    u_x = - d[1] - u_y;

        Simplex 3:  u_y = - d[2];
                    u_x = d[3] + u_y;

        Simplex 4:  u_x = d[4];
                    u_y = - d[3] + u_x;

        Simplex 5:  u_x = d[4];
                    u_y = d[5] - u_x;

        Simplex 6:  u_y = d[6];
                    u_x = d[5] - u_y;

        Simplex 7:  u_y = d[6];
                    u_x = - d[7] + u_y;

        Simplex 8:  u_x = - d[0];
                    u_y = d[7] + u_x;
                                            */
void CTimeDependentHjbSolver::eval_gradient(const int i_diff, const int i_direct, 
                                            const double v_diag, 
                                            const double v_dir, 
                                            double (&u)[2]) const {
  const int div = (i_direct+2) / 4;
  if (v_dir < 0) {
    u[0] = 0;
    u[1] = 0;
  }
  else if (v_diag <= 0) {
    if (i_direct == 2 || i_direct == 6) {
      u[1] = pow(-1, div) * v_dir;
      u[0] = 0;
    } else {
      u[0] = pow(-1, div+1) * v_dir;
      u[1] = 0;
    }
  } else {
    const double u_xy = v_diag;
    if (i_direct == 2 || i_direct == 6) {
      u[1] = pow(-1, div) * v_dir;
      if (abs(v_diag) < abs(v_dir)) {
        u[0] = 0;
      } else {
        u[0] = i_diff * pow(-1, div + 1) * (u_xy + pow(-1, div + 1)*u[1]);
      }
    } else {
      u[0] = pow(-1, div+1) * v_dir;
      if (abs(v_diag) < abs(v_dir)) {
        u[1] = 0;
      } else {
        u[1] = i_diff * pow(-1, div + 1) * (u_xy + pow(-1, div)*u[0]);
      }
    }
  }
}

/*==============================================================================
  Five-point stencil HJB update
==============================================================================*/
double CTimeDependentHjbSolver::update_5pt(const int aI, const int aJ, 
                                           const int aK) const {
  /* Grid sizes */
  const int nx = fGrid->getGridSizeX();
  const int ny = fGrid->getGridSizeY();

  /* Get local variables */
  const double x = fGrid->xGridToPhysical(aI);
  const double y = fGrid->yGridToPhysical(aJ);
  const double u = fGrid->getValue(aI,aJ,aK);
  const double h = fGrid->getH();
          
  /* Calculate directional derivatives of value function */
  double dplusx, dminusx, dplusy, dminusy;
  if (fTerrain->isNotInObstacles(x+h,y) && aI < nx-1) {
    dplusx  = (u-fGrid->getValue(aI+1,aJ,aK))/h;
  } else {
    dplusx = 0;
  }
  if (fTerrain->isNotInObstacles(x-h,y) && aI > 0) {
    dminusx = (u-fGrid->getValue(aI-1,aJ,aK))/h;
  } else {
    dminusx = 0;
  }
  if (fTerrain->isNotInObstacles(x,y+h) && aJ < ny-1) {
    dplusy  = (u-fGrid->getValue(aI,aJ+1,aK))/h;
  } else {
    dplusy = 0;
  }
  if (fTerrain->isNotInObstacles(x,y-h) && aJ > 0) {
    dminusy = (u-fGrid->getValue(aI,aJ-1,aK))/h;
  } else {
    dminusy = 0;
  }

  /* Calculate partial derivatives using upwind discretization */
  double u_x, u_y;
  if (dplusx >= dminusx && dplusx > 0) {
    u_x = -dplusx;
  } else if (dminusx >= dplusx && dminusx > 0) {
    u_x = dminusx;
  } else {
    u_x = 0;
  }
  if (dplusy >= dminusy && dplusy > 0) {
    u_y = -dplusy;
  } else if (dminusy >= dplusy && dminusy > 0) {
    u_y = dminusy;
  } else {
    u_y = 0;
  }

  /* Euler step */
  const double K = fGrid->getCost(aI,aJ,aK);
  const double f = fGrid->getSpeed(aI,aJ);
  const double dt = fGrid->getDt();
  return u + K*dt - f*sqrt(u_x*u_x + u_y*u_y)*dt;
}