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
* File: WriteToFile.hpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This file contains helper functions for writing multi-dimensional
* Boost arrays and vectors to file, and for printing nicely-formatted vectors to
* the command line.
*
* ==============================================================================
*/

#ifndef WRITE_TO_FILE_HPP
#define WRITE_TO_FILE_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen>

/** ----- Project-specific header files --------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "MemoryAllocations.hpp"

namespace io {

/** 
* This function writes the 2D Boost array aArray to a file with name aFilename
*/
template <class T>
void writeToFile2D(std::string aFilename, memory::array2D_t<T> aArray) {
  const int aDim0 = aArray.shape()[0];
  const int aDim1 = aArray.shape()[1];
  std::ofstream out("output/" + aFilename, std::ios::binary);

  for (int i = 0; i < aDim0; ++i) {
    for (int j = 0; j < aDim1; ++j) {
      out.write((char*) &aArray[i][j], sizeof(T));
    }
  }
}

/** 
* This function writes the 3D Boost array aArray to a file with name aFilename
*/
template <class T>
void writeToFile3D(std::string aFilename, memory::array3D_t<T> aArray) {
  const int aDim0 = aArray.shape()[0];
  const int aDim1 = aArray.shape()[1];
  const int aDim2 = aArray.shape()[2];
  std::ofstream out("output/" + aFilename, std::ios::binary);

  for (int i = 0; i < aDim0; ++i) {
    for (int j = 0; j < aDim1; ++j) {
      for (int k = 0; k < aDim2; ++k) {
        out.write((char*) &aArray[i][j][k], sizeof(T));
      }
    }
  }
}

/** 
* This function writes the 1D std::vector aVec to a file with name aFilename
*/
template <class T>
void writeVectorToFile(std::string aFilename, std::vector<T> aVec) {
  std::ofstream out("output/" + aFilename);
  for (int i = 0; i < aVec.size(); ++i) {
    out.write((char*) &aVec[i], sizeof(T));
  }
}

/** 
* This function writes the 2D std::vector aVec to a file with name aFilename
*/
template <class T>
void writeVectorToFile2D(std::string aFilename, 
                              std::vector<std::vector<T>> aVec) {
  std::ofstream out("output/" + aFilename);
  for (int i = 0 ; i < aVec.size(); ++i) {
    for (int j = 0; j < aVec[i].size(); ++j) {
      out.write((char*) &aVec[i][j], sizeof(T));
    }
  }
}

/**
 * Prints Eigen::VectorXd aEigenVec to the command line
 */
inline void printEigenVec(Eigen::VectorXd aEigenVec) {
  const int dimension = aEigenVec.size();
  std::cout << "(";
  for (int k = 0; k < dimension; ++k) {
    std::cout << aEigenVec[k]; 
    if (k < dimension - 1) {
      std::cout << ", ";
    }
  }
  std::cout << ")";
}

} /* namespace io */

#endif
