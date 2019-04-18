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
* File: CTimeGrid.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This class is used to handle both the logical and physical
* representation of a 3D regular grid with fixed spacing in all 3 dimensions.
* It is used by CFMM, CTimeDependentHjbSolver, CTimeDependentTracer, and 
* CMovingObserver. This class is an interface for between the underlying
* data structures and the rest of the code. 
* Only the grid values are stored as an array. The speed and cost
* functions are stored only as function pointers, which are treated as inputs.
*
* (See also CTimeGrid.hpp)
*
* ==============================================================================
*/
#include "CTimeGrid.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <string>
#include "boost/multi_array.hpp"

/** ------ Project-specific header files -------------------------------------*/
#include "SimpleFunctions.hpp"
#include "WriteToFile.hpp"
#include "MemoryAllocations.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace memory;
using namespace std;

/*==============================================================================
  Constructor
================================================================================*/
CTimeGrid::CTimeGrid(function<double(int, int, int)> aCost, 
                     function<double(double, double)> aSpeed, 
                     const int aN, const int aNt, 
                     const double aPhysMin, const double aPhysMax, 
                     const double aTime) {
  assert(aPhysMax >= aPhysMin);
  assert(aTime > 0.0);
  fMinX = aPhysMin;
  fMinY = aPhysMin;
  fMinT = 0.0;
  fMaxX = aPhysMax;
  fMaxY = aPhysMax;
  fMaxT = aTime;
  fValues = make_shared<array3D_t<double>>(allocateArray3D<double>(aN,aN,aNt));
  fH = (aPhysMax - aPhysMin)/(aN-1);
  fDt = aTime / (aNt - 1);
  fCost = aCost;
  fSpeed = aSpeed;
}

/*==============================================================================
  Write grid to file
==============================================================================*/
void CTimeGrid::writeGridToFile(const string aFilename) const {
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  const int nt = getGridSizeT();
  array3D_t<double> CostArray = allocateArray3D<double>(nx, ny, nt);
  for (int i=0; i < nx; ++i) {
    for (int j=0; j < ny; ++j) {
      for (int k=0; k < nt; ++k) {
        CostArray[i][j][k] = getCost(i,j,k);
      }
    }
  }
  io::writeToFile3D<double>(aFilename + "Value", (*fValues));
  io::writeToFile3D<double>(aFilename + "Cost",  CostArray);
}

/*==============================================================================
  Setters
==============================================================================*/
void CTimeGrid::setSpeed(function<double(double, double)> aSpeed) {
  fSpeed = aSpeed;
}
void CTimeGrid::setCost(function<double(int, int, int)> aCost) {
  fCost = aCost;
}

/*==============================================================================
  Returns value function at the physical location (aX,aY,aT) 
      using trilinear interpolation
==============================================================================*/
double CTimeGrid::getValuePhysical(const double aX, const double aY, 
                                   const double aT) const {
  /* Grid sizes  */
  const int nx = fValues->shape()[0];
  const int ny = fValues->shape()[1];
  const int nt = fValues->shape()[2];

  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);
  const int k = floor((aT - fMinT) / fDt);

  /* Find interpolation coefficients */
  const double r_x = (aX - i*fH) / fH;
  const double r_y = (aY - j*fH) / fH;
  const double r_t = (aT - k*fDt) / fDt;

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
  const double u1 = (1-r_t)*getValue(i,j,k)   + (r_t)*getValue(i,j,k1);
  const double u2 = (1-r_t)*getValue(i1,j,k)  + (r_t)*getValue(i1,j,k1);
  const double u3 = (1-r_t)*getValue(i,j1,k)  + (r_t)*getValue(i,j1,k1);
  const double u4 = (1-r_t)*getValue(i1,j1,k) + (r_t)*getValue(i1,j1,k1);

  /* Interpolate along y-axis */
  const double u5 = (1-r_y)*u1 + (r_y)*u3;
  const double u6 = (1-r_y)*u2 + (r_y)*u4;

  /* Interpolate along x-axis */
  return (1-r_x)*u5 + (r_x)*u6;
}