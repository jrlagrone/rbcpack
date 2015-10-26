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
  This file allows CMake to generate a configuration header to save the definitions
  that provide optional features 
*/

#ifndef CRBC_CONFIG_H_
#define CRBC_CONFIG_H_

// this defines the field data type --- 0 = float, 1 = double
/* #undef USE_CRBC_SCALAR_FLOAT */

// this defines the coefficient data type --- 0 = float, 1 = double
/* #undef USE_CRBC_COEF_FLOAT */

// The following two are mutually excluside and allow the indexing type to be
// changed from the default of int. Set the first to 1 to use long and the second
// to 1 to use long long. If both are set to 1, the type will be long
/* #undef USE_CRBC_INDEX_LONG */
/* #undef USE_CRBC_INDEX_LONG_LONG */

// Enable HDF5 functionality (requires HDF5)
/* #undef USE_HDF5_RESTARTS */

// Enable OpenMP functionality (requires OpenMP)
/* #undef USE_OPENMP */

#endif // CRBC_CONFIG_H_
