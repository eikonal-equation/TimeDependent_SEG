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
* File: CTimeGrid.hpp
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
* ==============================================================================
*/

#ifndef CTIMEGRID_HPP
#define CTIMEGRID_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <cmath>
#include <string>
#include "boost/multi_array.hpp"

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"
#include "SimpleFunctions.hpp"

class CTimeGrid
{
  private:
    /** Properties of grid */
    /* Grid spacing: assuming to be the same in the two spatial variables */
    double fH; 
    /* Length of timestep */
    double fDt; 
    /* Grid of values. Stored as a pointer to a 3D boost::multi_array */
    std::shared_ptr<memory::array3D_t<double>> fValues; 
    /* Physical bounds on the grid */
    double fMinX;
    double fMinY;
    double fMinT;
    double fMaxX;
    double fMaxY;
    double fMaxT;

    /** Speed function, stored as a function pointer. 
        This has to be given as an input.  
        Computes speed f at physical location (x,y).
        Signature is fSpeed(x,y).  */
    std::function<double(double, double)> fSpeed; 

    /** Cost function, stored as a function pointer. 
        This has to be given as an input.  
        Computes cost K at gridpoint (i,j,k).
        Signature is fCost(i,j,k).  */
    std::function<double(int, int, int)> fCost;

  public:
    /** ========================================================================
    *    Constructors
    * ========================================================================*/
    /** Default constructor */
    CTimeGrid() = default;

    /**
     * This is the main CTimeGrid constructor. 
     * It initializes the grid and computes the grid spacing
     * for a square grid with aN points in both spatial dimensions,
     * and aNt points in the time dimension.
     * @param aCost cost function pointer. Signature is aCost(i,j,k).
     * @param aSpeed speed function pointer. Defaults to a constant speed.
     * @param aN an integer. The number of spatial grid points. Defaults to 51.
     * @param aNt an integer. The number of time steps. Defaults to 100.
     * @param aPhysMin a real number. The lower bound on the physical grid.
     * @param aPhysMax a real number. The upper bound on the physical grid.
     * @param aTime a real number. The length of the time Interval [0, aTime]
     * The domain is [aPhysMin,aPhysMax] x [aPhysMin,aPhysMax] x [0,aTime].
     */
    CTimeGrid(std::function<double(int, int, int)> aCost, 
              std::function<double(double, double)> aSpeed = &constant2, 
              const int aN = 51, const int aNt = 301, 
              const double aPhysMin = 0, const double aPhysMax = 1, 
              const double aTime = 1);

    /** ========================================================================
    *    Setters 
    *=========================================================================*/
    void setValue(const int aI, const int aJ, const int aK, 
                  const double aValue);
    void setSpeed(std::function<double(double, double)> aSpeed);
    void setCost(std::function<double(int, int, int)> aCost);

    /** ========================================================================
    *    Getters
    * ========================================================================*/
    /** Logical-coordinate functions */
    double getValue(const int aI, const int aJ, const int aK) const;
    double getSpeed(const int aI, const int aJ) const;
    double getCost(const int aI, const int aJ, const int aK) const;

    /** Physical-coordinate functions */
    double getValuePhysical(const double aX, const double aY, 
                            const double aT) const;
    double getSpeedPhysical(const double aX, const double aY) const;

    /** Grid parameters */
    int getGridSizeX() const;
    int getGridSizeY() const;
    int getGridSizeT() const;
    double getH() const;
    double getDt() const;
    double getMinX() const;
    double getMinY() const;
    double getMinT() const;
    double getMaxX() const;
    double getMaxY() const;
    double getMaxT() const;

    /** ========================================================================
    *    Other
    * ========================================================================*/
    /** Mapping back and forth between grid and physical */
    double xGridToPhysical(const int aI) const;
    double yGridToPhysical(const int aJ) const;
    double tGridToPhysical(const int aK) const;
    int xPhysicalToGrid(const double aX) const;
    int yPhysicalToGrid(const double aY) const;
    int tPhysicalToGrid(const double aT) const;

    /**
     * A I/O member which prints the value and cost grids to file.
     * @param aFilename a string which contains the prefix to the names 
     *      of the files to which the grids will be printed. 
     * The (cost/value) grids will be printed to
     *      files called "aFilename"+(Cost/Value)
     */
    void writeGridToFile(const std::string aFilename) const;
};

/** ============================================================================
*    Inline function definitions
* ============================================================================*/
/** Inline definitions must be in the header file. 
 * These functions are used frequently. */

/** ------ Inline definition of setter function ------------------------------*/
inline void CTimeGrid::setValue(const int aI, const int aJ, const int aK, 
                                const double aValue) {
  (*fValues)[aI][aJ][aK] = aValue;
}

/** ------ Inline definition of getters --------------------------------------*/
inline double CTimeGrid::getValue(const int aI, const int aJ, 
                                  const int aK) const {
  return (*fValues)[aI][aJ][aK];
}

inline double CTimeGrid::getSpeed(const int aI, const int aJ) const {
  return fSpeed(xGridToPhysical(aI), yGridToPhysical(aJ));
}

inline double CTimeGrid::getCost(const int aI, const int aJ, 
                                 const int aK) const {
    return fCost(aI,aJ,aK);
}

inline double CTimeGrid::getSpeedPhysical(const double aX, 
                                          const double aY) const {
  return fSpeed(aX, aY);
}

inline int CTimeGrid::getGridSizeX() const{
  return fValues->shape()[0];
}
inline int CTimeGrid::getGridSizeY() const{
  return fValues->shape()[1];
}
inline int CTimeGrid::getGridSizeT() const{
  return fValues->shape()[2];
}
inline double CTimeGrid::getH() const{
  return fH;
}
inline double CTimeGrid::getDt() const{
  return fDt;
}
inline double CTimeGrid::getMinX() const{
  return fMinX;
}
inline double CTimeGrid::getMinY() const{
  return fMinY;
}
inline double CTimeGrid::getMinT() const{
  return fMinT;
}
inline double CTimeGrid::getMaxX() const{
  return fMaxX;
}
inline double CTimeGrid::getMaxY() const{
  return fMaxY;
}
inline double CTimeGrid::getMaxT() const{
  return fMaxT;
}

/** ------ Inline definition of coordinate functions -------------------------*/
inline double CTimeGrid::xGridToPhysical(const int aI) const {
  return fMinX + (double)aI * fH;
}

inline double CTimeGrid::yGridToPhysical(const int aJ) const {
  return fMinY + (double)aJ * fH;
}

inline double CTimeGrid::tGridToPhysical(const int aK) const {
  return fMinT + (double)aK * fDt;
}

inline int CTimeGrid::xPhysicalToGrid(const double aX) const {
  return std::round((aX - fMinX)/fH);
}

inline int CTimeGrid::yPhysicalToGrid(const double aY) const {
  return std::round((aY - fMinY)/fH);
}

inline int CTimeGrid::tPhysicalToGrid(const double aT) const {
  return std::round((aT - fMinT)/fDt);
}

#endif
