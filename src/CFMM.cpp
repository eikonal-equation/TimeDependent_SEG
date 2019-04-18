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
* File: CFMM.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This is a class for solving stationary Eikonal equations using
* the Fast Marching Method. 
* It works on a 2D regular grid and a 5 point stencil. 
* It assumes equal spacing in both horizontal and vertical directions.
* The target set is determined by aBoundary. 
* The solution to the Eikonal is stored in the fGrid field.
* The FMM implementation uses the boost::heap and CTimeGrid classes.
*
* ==============================================================================
*/
#include "CFMM.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CTimeGrid.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;
using namespace memory;

CFMM::CFMM(shared_ptr<CTimeGrid> aGrid) {
  /** Initialize properties and allocate arrays. */
  fGrid = aGrid;
  const int nx = fGrid->getGridSizeX();
  const int ny = fGrid->getGridSizeY();
  fStatus = make_unique<array2D_t<status_t>>(allocateArray2D<status_t>(nx,ny));
  fHeapPointers = make_unique<array2D_t<handle_t>>(allocateArray2D<handle_t>(nx,ny));
}

/** The CFMM::march function has been adapted to use more general boundary 
 * conditions in order to compute the signed distance function used in 
 * visibility computations */
void CFMM::march(array2D_t<bool> aBoundary) {
  /* Initialize heap */
  CFMMHeap_t CFMMHeap;
  initialize_cfmm(CFMMHeap, aBoundary);

  /* perform Dijkstra's algorithm to update all grid points */
  while (!CFMMHeap.empty()) {
    CHeapGP current_GP = CFMMHeap.top();
    CFMMHeap.pop();
    (*fStatus)[current_GP.fI][current_GP.fJ] = ACCEPTED;
    update_neighbors(CFMMHeap, current_GP.fI, current_GP.fJ);
  }
}

/** The CFMM::initialize_cfmm function has been adapted to use more general 
 * boundary conditions in order to compute the signed distance function used in 
 * visibility computations */
void CFMM::initialize_cfmm(CFMMHeap_t& aCFMMHeap, array2D_t<bool> aBoundary) {
  for (int i = 0; i < fGrid->getGridSizeX(); ++i) {
    for (int j = 0; j < fGrid->getGridSizeY(); ++j) {
      if (aBoundary[i][j]) {
        /* All points on the border are set to ACCEPTED. */
        fGrid->setValue(i, j, 0, 0);
        (*fStatus)[i][j] = ACCEPTED;
        aCFMMHeap.push(CHeapGP(i, j, 0));
      } else {
        /* All other grid points are set to FAR. */
        fGrid->setValue(i, j, 0, INF);
        (*fStatus)[i][j] = FAR;
      }
    }
  }
}

/*----------------------------------------------------------------------------//
/-- Helper functions
/-----------------------------------------------------------------------------*/
void CFMM::update_neighbors(CFMMHeap_t& aCFMMHeap, const int aCurrent_i, const int aCurrent_j) {
  /* Five-point Stencil */
  constexpr int stencil[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};

  /* Iterate over neighbors */
  for (int k = 0; k < 4; ++ k) {
    const int nbr_i = aCurrent_i + stencil[k][0];
    const int nbr_j = aCurrent_j + stencil[k][1];

    if (in_domain(nbr_i,nbr_j)) {
      /* If accepted nothing to do, otherwise, try to update */
      if ((*fStatus)[nbr_i][nbr_j] != ACCEPTED) {
        const bool gp_was_updated = update_gp(nbr_i, nbr_j)  ;
        /* If it was already considered and was updated then update heap */
        if ((*fStatus)[nbr_i][nbr_j] == CONSIDERED) {
          if (gp_was_updated) {
            aCFMMHeap.update((*fHeapPointers)[nbr_i][nbr_j],CHeapGP(nbr_i, nbr_j, fGrid->getValue(nbr_i,nbr_j,0)));
          }
        } else {
          /* Else add to heap */
          (*fStatus)[nbr_i][nbr_j] = CONSIDERED;
          (*fHeapPointers)[nbr_i][nbr_j] = aCFMMHeap.push(CHeapGP(nbr_i, nbr_j, fGrid->getValue(nbr_i,nbr_j,0)));
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------/
/-- Compute functions
/-----------------------------------------------------------------------------*/
double CFMM::compute_from_two_neighbors(const double aEffectiveSpeed, 
                                        const double aVerticalNeighbor, 
                                        const double aHorizontalNeighbor) {
  double newVal;

  /* If speed is zero, then value function should be set to infinity */
  if(aEffectiveSpeed == 0) {
    newVal = INF;
  }

  const double h = fGrid->getH();
  const double diff = aVerticalNeighbor - aHorizontalNeighbor;
  const double rat1 = h/aEffectiveSpeed;

  if (fabs(diff) >= rat1) {
    /* Upwinding condition is not satisfied. Use one-sided update instead. */
    if (aVerticalNeighbor < aHorizontalNeighbor) {
      newVal = aVerticalNeighbor + rat1;
    } else {
      newVal = aHorizontalNeighbor + rat1;
    }
  } else {
    /* Upwinding condition is satisfied. Use two-sided update */
    const double rat = SQRT2*rat1;
    const double discriminant = (rat - diff)*(rat + diff);
    assert(discriminant >= 0);

    newVal = .5*((aVerticalNeighbor + aHorizontalNeighbor) + sqrt(discriminant));
  }

  return newVal;
}

/*-----------------------------------------------------------------------------/
/-- Update gridpoint
/-----------------------------------------------------------------------------*/
bool CFMM::update_gp(const int aI, const int aJ) {
  /* Compute smaller vertical and horizontal neighbors. */
  const double val_h = smaller_h_neighbor(aI,aJ);
  const double val_v = smaller_v_neighbor(aI,aJ);

  /* Compute effective speed f = speed/cost */
  const double f = fGrid->getSpeed(aI,aJ)/fGrid->getCost(aI,aJ,0); 
  assert(!isnan(f));

  /* If speed is 0, then value has to be infinity there. */
  if (f == 0) {
    fGrid->setValue(aI,aJ,0,INF);
    return false;
  }

  /* Compute update based on smaller horizontal and vertical neighbors */
  const double new_value = compute_from_two_neighbors(f, val_h, val_v);
  if (new_value < fGrid->getValue(aI,aJ,0)) {
    fGrid->setValue(aI, aJ, 0, new_value);
    return true;
  } else {
    return false;
  }
}
