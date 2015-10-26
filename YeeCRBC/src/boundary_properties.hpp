/*
  Copyright 2015 John LaGrone

  This file is part of the Yee CRBC Library from RBCPACK.

  The Yee CRBC Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or (at your 
  option) any later version.

  The Yee CRBC Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the Yee CRBC Library.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
  Contains a class that holds enums used in many of the other CRBC classes
*/


#ifndef CRBC_BOUNDARY_PROPERTIES_H_
#define CRBC_BOUNDARY_PROPERTIES_H_

namespace crbc
{

class BoundaryProperties
{

  public:

    /// Enumeration for boundary types
    enum Boundary {
        NONE = -1, /**< no boundary type                               */
        DIR = 0,   /**< homogeneous Dirichlet                          */
        NEUM,      /**< homogeneous Neumann (in the normal direction)  */
        CRBC       /**< Complete radition boundary condition           */
    };

    /// Enumeration of the valid coordinate directions
    enum Coordinate {
        X1 = 0,     /**< x1-coordinate direction  */
        X2,         /**< x2-coordinate direction  */
        X3          /**< x3-coordinate direction  */
    };

    /// Enumeration of the valid boundary sides
    enum Side {
        X1Left = 0,     /**< Left side in the x1-coordinate direction  */
        X1Right,        /**< Right side in the x1-coordinate direction */
        X2Left,         /**< Left side in the x2-coordinate direction  */
        X2Right,        /**< Right side in the x2-coordinate direction */
        X3Left,         /**< Left side in the x3-coordinate direction  */
        X3Right         /**< Right side in the x3-coordinate direction */
    };

  }; // end class

} // end namespace

#endif //#define CRBC_BOUNDARY_PROPERTIES_H_
