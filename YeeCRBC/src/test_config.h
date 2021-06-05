/* This is the default configuration file. If building with Cmake, this will
   be overwritten with the choices make in Cmake.
*/

#ifndef CRBC_CONFIG_H_
#define CRBC_CONFIG_H_

// this defines the field data type --- 0 = float, 1 = double
#define USE_CRBC_SCALAR_FLOAT 0

// this defines the coefficient data type --- 0 = float, 1 = double
#define USE_CRBC_COEF_FLOAT 0

// The following two are mutually excluside and allow the indexing type to be
// changed from the default of int. Set the first to 1 to use long and the second
// to 1 to use long long. If both are set to 1, the type will be long
#define USE_CRBC_INDEX_LONG 0
#define USE_CRBC_INDEX_LONG_LONG 0

// Enable HDF5 functionality (requires HDF5)
#define USE_HDF5_RESTARTS 0

// Enable OpenMP functionality (requires OpenMP)
#define USE_OPENMP 0

#endif // CRBC_CONFIG_H_
