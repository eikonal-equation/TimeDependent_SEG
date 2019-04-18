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
* File: CMovingObserver.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This is the class for moving Observers. It implements the 
* CObserver interface, and calculates the pointwise observability in the domain.
* 
* We first calculate a "signed distance function" Phi that is 
*     -positive outside the obstacles
*     -negative inside the obstacles
* This is calculated using FMM to solve an Eikonal equation, with the boundaries
* of the obstacles as the boundary conditions for the Eikonal.
* Then the values inside the obstacles are flipped as a post-processing step.
* Then we calculate a visibility function Psi for each gridpoint via the 
* quasi-variational technique in (Tsai, Cheng, Osher, Burchard, Sapiro 2004).
* Then we loop over all of the gridpoints, labeling them as visible or not based
* on the values of Psi and including any angular restrictions on the observer.
*
* (See also CMovingObserver.hpp)
*
* ==============================================================================
*/
#include "CMovingObserver.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CFMM.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace memory;
using namespace std;
using namespace std::placeholders;

/*==============================================================================
  Constructor
==============================================================================*/
CMovingObserver::CMovingObserver(const std::function<double(double)> aXpath, 
                                 const std::function<double(double)> aYpath, 
                                 const std::shared_ptr<CTerrain> aTerrain, 
                                 const std::function<double(double, double,double,double)> aObservability, 
                                 const int aN, const int aNt, 
                                 const double aPhysMin, const double aPhysMax, 
                                 const double aTime, 
                                 const std::function<double(double)> aXder, 
                                 const std::function<double(double)> aYder, 
                                 const double aAlpha) {
  /* Copy arguments */
  fXpath = aXpath;
  fYpath = aYpath;
  fTerrain = aTerrain;
  fObservability = aObservability;
  fN = aN;
  fNt = aNt;
  fPhysMin = aPhysMin;
  fPhysMax = aPhysMax;
  fTime = aTime;
  fXder = aXder;
  fYder = aYder;
  fAlpha = aAlpha;
  
  /* Compute shadow zones */
  compute_visibility();
}

/*
================================================================================
COMPUTING THE SHADOW ZONE.
================================================================================

DESCRIPTION:

The shadow zone is computed in three steps:

1) Signed distance function
For every gridpoint, we first compute the distance to the nearest obstacle boundary. 
This is computed by solving an Eikonal equation using FMM
  where all gridpoints on the boundaries initialized to 0. 
We then flip the sign of the gridpoints inside obstacles. 
We end up with a signed distance function Phi, whose value is:
  positive outside the obstacles
  negative inside the obstacles 
  zero if the gridpoint is on the boundary

2) Line of sight 
For each gridpoint, we compute the minimum of the signed distance function Phi 
  along the line connecting the observer and the gridpoint.
We end up with a function Psi, whose value is:
  negative if there is no direct line of sight from the observer to the gridpoint
  non-negative if there is a direct line of sight from the observer to the gridpoint

3) Angular restrictions
Gridpoints are then set to be visible if:
  a) The angle between the observer's heading 
        and the line connecting the observer and gridpoint is less than fAlpha
  b) There is a direct line of sight from the observer to the gridpoint 
        (i.e. Psi >= 0)

In implementation terms:
compute_visibility() is the main function that computes all of these steps
compute_signed_distance() is a helper function for step (1)
compute_line_of_sight() is a helper function for step (2)

(Note: The implementation assumes a square grid and stationary obstacles)

================================================================================
*/

/*==============================================================================
  Signed distance function Phi
    computed via Fast Marching Method
==============================================================================*/
void CMovingObserver::compute_signed_distance() {
  /* Get grid sizes */
  const int nx = fN;
  const int ny = fN;
  const int nt = fNt;

  /* Phi will store the signed distance function */
  fPhi = make_shared<CTimeGrid>(&constantcost3, &constant2, fN, 1, fPhysMin, fPhysMax, fTime);
  
  /** boundary[i][j]==true if (i,j) is on the boundary.
  * This stores all points that will be initialized as ACCEPTED for CFMM. */
  array2D_t<bool> boundary = allocateArray2D<bool>(nx, ny);

  /** in_obstacles[i][j]==1 if (i,j) is inside obstacles.
  * This stores all points whose signs of Phi values will be flipped. */
  array2D_t<bool> in_obstacles = allocateArray2D<bool>(nx, ny);

  /** Go through all gridpoints to check if they are on the boundary or inside obstacles */
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      const double x = fPhi->xGridToPhysical(i);
      const double y = fPhi->yGridToPhysical(j);
      if(fTerrain->isNotOnBoundary(x,y,0)) {
        boundary[i][j] = false;
      } else {
        boundary[i][j] = true; 
      }
      if(fTerrain->isNotInObstacles(x,y,0)) {
        in_obstacles[i][j] = false;
      } else {
        in_obstacles[i][j] = true;
      }
    }
  }

  /** Run FMM to solve the Eikonal 
   * with all points on the boundary initialized as 0 and ACCEPTED */
  CFMM FMMobject = CFMM(fPhi);
  FMMobject.march(boundary);

  /** Flip the Phi value of gridpoints inside obstacles */
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (in_obstacles[i][j]) {
        const double tmp = fPhi->getValue(i, j, 0);
        fPhi->setValue(i, j, 0, -tmp);
      }
    }
  }
}


/*==============================================================================
  Line of sight function Psi
    computed via taking minimum of signed distance function Phi 
    along line connecting observer to gridpoint
    (see Tsai, Cheng, Osher, Burchard, Sapiro 2004)
==============================================================================*/
void CMovingObserver::compute_line_of_sight() {
  /* Stencils used for dynamically handling the four quadrants */
  constexpr int stencil1[4][2] = {{1,1},{-1, -1}, {1, -1}, {-1, 1}};
  constexpr int stencil2[4][2] = {{1,0},{-1, 0}, {0, 1}, {0, -1}};

  /* Get grid sizes */
  const int nx = fVisibility->getGridSizeX();
  const int ny = fVisibility->getGridSizeY();
  const int nt = fVisibility->getGridSizeT();

  /* ------ Main Loop starts here --------------------------------------------*/
  for(int step = 0; step < nt; ++step){
    /** Initialization */
    for(int i = 0; i < nx; ++i) {
      for(int j = 0; j < ny; ++j) {
        fVisibility->setValue(i,j,step, fPhi->getValue(i,j,0));
      }
    }

    /** Obtain integer grid position of observer */
    const double t = fVisibility->tGridToPhysical(step);
    const int i_observer = fVisibility->xPhysicalToGrid(getPosX(t));
    const int j_observer = fVisibility->yPhysicalToGrid(getPosY(t));

    /** First update Psi along the diagonals */
    for (int k = 0; k < 4; ++k) {
      int i = i_observer;
      int j = j_observer;

      
      while ((i + stencil1[k][0] >= 0) && (i + stencil1[k][0] < nx) && (j + stencil1[k][1] >= 0) && (j + stencil1[k][1] < ny)) {
        i = i + stencil1[k][0];
        j = j + stencil1[k][1];
        const double tmp = min(fPhi->getValue(i, j, 0), fVisibility->getValue(i-stencil1[k][0], j-stencil1[k][1], step));
        fVisibility -> setValue(i, j, step, tmp);
      }
    }

    /** Update Psi values on horizontal and vertical lines */
    for (int k = 0; k < 4; ++k) {
      int i = i_observer;
      int j = j_observer;

      while (( i + stencil2[k][0] >= 0 ) && (i + stencil2[k][0] < nx) && (j + stencil2[k][1] >= 0 ) && (j + stencil2[k][1] < ny)) {
        i = i + stencil2[k][0];
        j = j + stencil2[k][1];
        const double tmp = min(fPhi->getValue(i, j, 0), fVisibility->getValue(i-stencil2[k][0], j-stencil2[k][1], step));
        fVisibility -> setValue(i, j, step, tmp);
      }
    }

    /** Update Psi values everywhere else */
    for (int k = 0; k < 4; ++k) {
      for (int i = i_observer + stencil1[k][0]; (i >= 0) && (i < nx); i = i + stencil1[k][0]) {
        for (int j = j_observer + stencil1[k][1]; (j >= 0) && (j < ny); j = j + stencil1[k][1]) {
          const double slope = (j-j_observer)*1.0/(i-i_observer);
          const double intercept = j*1.0 - slope * (i*1.0);
          if (fabs(i-i_observer) > fabs(j-j_observer)) {
            const double r = fabs(slope * (i-stencil1[k][0]) + intercept - j);
            const double tmp1 = (1-r)*fVisibility->getValue(i-stencil1[k][0], j, step) + r*fVisibility->getValue(i-stencil1[k][0], j-stencil1[k][1], step);
            const double tmp2 = fPhi->getValue(i, j, 0);
            fVisibility -> setValue(i, j, step, min(tmp1, tmp2));
          } else if (fabs(i-i_observer) < fabs(j-j_observer)) {
            const double r = fabs((j*1.0-stencil1[k][1]-intercept)/slope - i);
            const double tmp1 = (1-r)*fVisibility->getValue(i, j-stencil1[k][1], step) + r*fVisibility->getValue(i-stencil1[k][0], j-stencil1[k][1], step);
            const double tmp2 = fPhi->getValue(i, j, 0);
            fVisibility -> setValue(i, j, step, min(tmp1, tmp2));
          }
        }
      }
    }
  }
}

/*==============================================================================
  Compute Visibility
    Main function for computing visibility
==============================================================================*/
void CMovingObserver::compute_visibility() {
  /* Allocate arrays and compute line of sight */
  fVisibility = make_shared<CTimeGrid>(&constantcost3, &constant2, fN, fNt, fPhysMin, fPhysMax, fTime);
  compute_signed_distance();
  compute_line_of_sight();

  fPhi->writeGridToFile("temp_phi");
  fVisibility->writeGridToFile("temp_psi");

  /* Loop through gridpoints to determine visibility */
  const int nx = fVisibility->getGridSizeX();
  const int ny = fVisibility->getGridSizeY();
  const int nt = fVisibility->getGridSizeT();
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nt; ++k) {
        const double x = fVisibility->xGridToPhysical(i);
        const double y = fVisibility->yGridToPhysical(j);
        const double t = fVisibility->tGridToPhysical(k);

        /* Compute angle between line from observer to gridpoint, and observer's heading */
        const double x_dist = x-getPosX(t);
        const double y_dist = y-getPosY(t);
        const double numerator = x_dist * fXder(t) + y_dist * fYder(t);
        const double norm1 = sqrt(x_dist * x_dist + y_dist * y_dist);
        const double norm2 = sqrt(fXder(t) * fXder(t) + fYder(t) * fYder(t));
        const double cos_theta = numerator/(norm1 * norm2);

        /* Mark invisible if no line of sight, or angle restriction is violated */
        if (fVisibility->getValue(i,j,k) < 0 || cos_theta < cos(fAlpha/2)) {
          fVisibility->setValue(i,j,k,0);
        } else {
          fVisibility->setValue(i,j,k,1);
        }
      }
    }
  }
}

/*==============================================================================
  Interpolated observability function
    Returns observability at the physical location (aX,aY,aT).
    Only used for subgradient computation/path tracer. 
    The HJB solver uses getObservabilityGrid() instead.
==============================================================================*/
double CMovingObserver::getObservabilityPhysical(const double aX, 
                                                 const double aY, 
                                                 const double aT) const {
  if (INTERPOLATE_OBSERVABILITY) {
    /** Trilinear interpolation of getObservabilityGrid() */
    const double h  = fVisibility->getH();
    const double dt = fVisibility->getDt();
    const int nx = fVisibility->getGridSizeX();
    const int ny = fVisibility->getGridSizeY();
    const int nt = fVisibility->getGridSizeT();

    /* Find logical coordinates of "lower-left corner" */
    const int i = floor((aX - fVisibility->getMinX()) / h);
    const int j = floor((aY - fVisibility->getMinY()) / h);
    const int k = floor((aT - fVisibility->getMinT()) / dt);

    /* Interpolation coefficients */
    const double r_x = (aX - i*h) / h;
    const double r_y = (aY - j*h) / h;
    const double r_t = (aT - k*dt) / dt;

    /* Find gridpoints for neighbors (used to handle boundary issues) */
    int i1 = i + 1;
    int j1 = j + 1;
    int k1 = k + 1;
    if (i == nx-1) { /* Check if point is on right side of grid */
      i1 = i;
    }
    if (j == ny-1) { /* Check if point is on top side of grid */
      j1 = j;
    }
    if (k == nt-1) { /* Check if point is at final timestep */
      k1 = k;
    }

    /** Trilinear interpolation */
    /* Interpolate along t-axis */
    const double K1 = (1-r_t)*getObservabilityGrid(i,j,k)   + (r_t)*getObservabilityGrid(i,j,k1);
    const double K2 = (1-r_t)*getObservabilityGrid(i1,j,k)  + (r_t)*getObservabilityGrid(i1,j,k1);
    const double K3 = (1-r_t)*getObservabilityGrid(i,j1,k)  + (r_t)*getObservabilityGrid(i,j1,k1);
    const double K4 = (1-r_t)*getObservabilityGrid(i1,j1,k) + (r_t)*getObservabilityGrid(i1,j1,k1);

    /* Interpolate along y-axis */
    const double K5 = (1-r_y)*K1 + (r_y)*K3;
    const double K6 = (1-r_y)*K2 + (r_y)*K4;

    /* Interpolate along x-axis */
    return (1-r_x)*K5 + (r_x)*K6;
  } else {
    /* Round to nearest grid coordinate */
    const int i = fVisibility->xPhysicalToGrid(aX);
    const int j = fVisibility->yPhysicalToGrid(aY);
    const int k = fVisibility->tPhysicalToGrid(aT);
    return getObservabilityGrid(i,j,k);
  }
}