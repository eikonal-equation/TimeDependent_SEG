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
* File: CFMM.hpp
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

#ifndef CFMM_HPP
#define CFMM_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <boost/heap/binomial_heap.hpp>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CTimeGrid.hpp"

class CFMM
{
  public:
    CFMM() = default;

    /**
     * The main CFMM Constructor
     * @param aGrid shared pointer to CTimeGrid object 
     *        on which the Eikonal equation is solved.
     */
    CFMM(std::shared_ptr<CTimeGrid> aGrid);

    /**
     * This function actually computes the solution to the Eikonal equation
     *      using the Fast Marching Method.
     * The solution is stored in fGrid (accessible by getValue/setValue).
     * @param aBoundary 2D array defining the boundary conditions.
     *        if aBorder[i][j]==1, then u(x[i],y[j])=0
     */
    void march(memory::array2D_t<bool> aBoundary);

  private:
    /* A GridPoint Class to be used by the heap. */
    class CHeapGP{
      public:
        int fI;
        int fJ;
        double fValue;
        CHeapGP(int aI, int aJ, double aValue): fI(aI), fJ(aJ), fValue(aValue) {};
    };

    /* Struct comparison which is required by boost::heap */
    struct compare_CHeapGP
    {
      bool operator()(const CHeapGP& aPoint1, const CHeapGP& aPoint2) const
      {
        return aPoint1.fValue > aPoint2.fValue;
      }
    };

    /* Typedef Heap types to make it easier to read. */
    typedef boost::heap::binomial_heap<CHeapGP, boost::heap::compare<compare_CHeapGP> > CFMMHeap_t;
    typedef typename boost::heap::binomial_heap<CHeapGP, boost::heap::compare<compare_CHeapGP> >::handle_type handle_t;

    /* A shared pointer to a CTimeGrid object
    *  which contains all of the information (cost, speed, physical grid)
    *  about the Eikonal PDE */
    std::shared_ptr<CTimeGrid> fGrid;

    /* The status (far/considered/accepted) of each grid point. */
    std::unique_ptr<memory::array2D_t<status_t>> fStatus;
    /* The grid of backpointers for the heap */
    std::unique_ptr<memory::array2D_t<handle_t>> fHeapPointers;


    /**
     * FMM Initialization.
     * This function sets the status of all nodes 
     *    and adds the boundary of the domain to the heap.
     * @param aCFMMHeap_t& aCFMMHeap Boost::heap passed by value 
     *    it is initialized by this function.
     * @param aBoundary 2D array defining the boundary conditions.
     *        if aBoundary[i][j]==true, then u(x[i],y[j])=0
     */
    void initialize_cfmm(CFMMHeap_t& aCFMMHeap, 
                        const memory::array2D_t<bool> aBoundary);

    /** A helper function which updates neighbors and adds them to the heap.
    * @param aCFMMHeap boost heap object
    * @param aCurrent_i int x logical coordinate of grid point
    * @param aCurrent_j int y logical coordinate of grid point
    */
    void update_neighbors(CFMMHeap_t& aCFMMheap, 
                         const int aCurrent_i, const int aCurrent_j);

    /**
     * A helper function to figure out if grid point (aI,aJ) is in the domain.
     * @param aI int x logical coordinate
     * @param aJ int y logical coordinate
     */
    bool in_domain(const int aI, const int aJ) const;

    /** Returns the value of the smaller of the two horizontal neighbors.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    double smaller_h_neighbor(const int aI, const int aJ) const;

    /** Returns the value of the smaller of the two vertical neighbors.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    double smaller_v_neighbor(const int aI, const int aJ) const;

    /** A compute function which returns the Eikonal update.
    * @param aEffectiveSpeed (= speed / cost) of fGrid.
    * @param aSmallerVerticalNeighbor value of the vertical neighbor.
    * @param aSmallerHorizontalNeighbor value of the horizontal neighbor.
    */
    double compute_from_two_neighbors(const double aEffectiveSpeed, 
                                      const double aVerticalNeighbor, 
                                      const double aHorizontalNeighbor);

    /** Update a gridpoint. Compute the update, assign, and return if updated.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    bool update_gp(const int aI, const int aJ);
};

/*-----------------------------------------------------------------------------/
/-- Inline function definitions
/-----------------------------------------------------------------------------*/

/* Compute if gridpoints [i,j] is inside the domain */
inline
bool CFMM::in_domain(const int aI, const int aJ) const {
  return (aI >= 0 && aI < fGrid->getGridSizeX() && aJ >=0 && aJ < fGrid->getGridSizeY());
}

/* Returns the value of the smaller of the two horizontal neighbors */
inline
double CFMM::smaller_h_neighbor(const int aI, const int aJ) const {
  int bestI;
  if (aI == 0) {
    bestI = aI + 1;
  } else if (aI == fGrid->getGridSizeX() - 1) {
    bestI = aI - 1;
  } else if (fGrid->getValue(aI - 1, aJ, 0) < fGrid->getValue(aI + 1, aJ, 0)) {
    bestI = aI - 1;
  } else {
    bestI = aI + 1;
  }
  return fGrid->getValue(bestI, aJ, 0);
}

/* Returns the value of the smaller of the two vertical neighbors */
inline
double CFMM::smaller_v_neighbor(const int aI, const int aJ) const {
  int bestJ;
  if (aJ == 0) {
    bestJ = aJ + 1;
  } else if (aJ == fGrid->getGridSizeY() - 1) {
    bestJ = aJ - 1;
  } else if (fGrid->getValue(aI, aJ - 1, 0) < fGrid->getValue(aI, aJ + 1, 0)) {
    bestJ = aJ - 1;
  } else{
    bestJ = aJ + 1;
  }
  return fGrid->getValue(aI, bestJ, 0);
}

#endif
