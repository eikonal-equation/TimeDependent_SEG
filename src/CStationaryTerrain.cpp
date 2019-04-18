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
* File: CStationaryTerrain.cpp
*
* Authors: Elliot Cartee, Qianli Song, Lexiao Lai
*   (based on code by Marc Aurèle Gilles)
*
* Description: This is the class for stationary terrain. It implements the 
* CTerrain interface, and is responsible for figuring out whether a given point
* in the domain is contained in any of the obstacles.
* (See also CStationaryTerrain.hpp)
*
* ==============================================================================
*/
#include "CStationaryTerrain.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <vector>

/** ------ Project-specific header files -------------------------------------*/
#include "CStationaryObstacle.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;

/*==============================================================================
  Destructor
==============================================================================*/
CStationaryTerrain::~CStationaryTerrain() {}

/*==============================================================================
  Default (i.e. no arguments) Constructor
==============================================================================*/
CStationaryTerrain::CStationaryTerrain() {
  fNumObstacles = 0;
}

/*==============================================================================
  Constructor for 2D vector of doubles
==============================================================================*/
CStationaryTerrain::CStationaryTerrain(vector<vector<double>> aSquareObstacles) {
  fObstacles.resize(aSquareObstacles.size());
  fNumObstacles = aSquareObstacles.size();
  for (int i = 0; i < fNumObstacles; ++i) {
    fObstacles[i] = make_shared<CStationaryObstacle>(aSquareObstacles[i][0],
        aSquareObstacles[i][1],aSquareObstacles[i][2],aSquareObstacles[i][3]);
  }
}

/*==============================================================================
  Constructor for vector of pointers to obstacles
==============================================================================*/
CStationaryTerrain::CStationaryTerrain(vector<shared_ptr<CObstacle>> aObstacles) {
  fObstacles = aObstacles;
  fNumObstacles = aObstacles.size();
}