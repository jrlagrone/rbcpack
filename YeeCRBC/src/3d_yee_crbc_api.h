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

/* C/C++ interface to the DAB/CRBC updates for the 3D Yee scheme

*/

/// \example yee.c

#ifndef THREE_D_YEE_CRBC_API_H
#define THREE_D_YEE_CRBC_API_H

// load the configuration file
#include "config.h"

#if USE_CRBC_SCALAR_FLOAT
  typedef float ScalarType;
#else
  typedef double ScalarType;
#endif

#if USE_CRBC_COEF_FLOAT
  typedef float CoefType;
#else
  typedef double CoefType;
#endif

#if USE_CRBC_INDEX_LONG
  typedef long int IndexType;
#elif USE_CRBC_INDEX_LONG_LONG
  typedef long long IndexType;
#else
  typedef int IndexType;
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct YeeCrbcUpdater;
typedef struct YeeCrbcUpdater YeeCrbcUpdater;

/// define enumerations for the valid boundary types
typedef enum CRBC_Boundaries
{
  CRBC_PEC,  /**< Perfect Electric Conductor */
  // PMC,  /**< Perfect Magnetic Conductor */
  CRBC_CRBC /**< Complete Radiation Boundary Condition (Implemented as a DAB) */
} CRBC_Boundaries_t;

/// define the valid Field Component Types
typedef enum CRBC_Fields
{
  CRBC_Ex, /**< x-component of the Electric Field */
  CRBC_Ey, /**< y-component of the Electric Field */
  CRBC_Ez /*,*/ /**< z-component of the Electric Field */
//  Hx, /**< x-component of the Magnetic Field */
//  Hy, /**< y-component of the Magnetic Field */
//  Hz /**< z-component of the Magnetic Field */
} CRBC_Fields_t;

/// define valid sides
typedef enum CRBC_Side {                                                            
  CRBC_XLeft,  /**< Left side in the x-coordinate direction */
  CRBC_XRight, /**< Right side in the x-coordinate direction */
  CRBC_YLeft,  /**< Left side in the y-coordinate direction */
  CRBC_YRight, /**< Right side in the y-coordinate direction */
  CRBC_ZLeft,  /**< Left side in the z-coordinate direction */
  CRBC_ZRight  /**< Right side in the z-coordinate direction */
} CRBC_Side_t;                                                                 

/// \brief create a new updater object with default settings
YeeCrbcUpdater* CRBC_new_updater (const CoefType T,
                             const CoefType h[3],  
                             const CoefType dt,                  
                             const CoefType c,
                             const CRBC_Boundaries_t boundaries[6]);

/// \brief create a new updater object and specify the number of recursions
YeeCrbcUpdater* CRBC_new_updater_p (const CoefType T,
                               const CoefType h[3],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC_Boundaries_t boundaries[6],
                               const int P);

/// \brief create a new updater object and specify the a tolerance for the 
///        recursions
YeeCrbcUpdater* CRBC_new_updater_tol (const CoefType T,
                               const CoefType h[3],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC_Boundaries_t boundaries[6],
                               const int Pmax,
                               const double tol);

/// \brief delete an updater object
void CRBC_delete_updater (YeeCrbcUpdater *U);

/// \brief initialize a component on a face
int CRBC_init_face (YeeCrbcUpdater *U,
               const CRBC_Side_t side,
               const CRBC_Fields_t comp, 
               const IndexType low_index[3],
               const IndexType high_index[3],
               const double delta);

/// \brief input new data into the boundary updater from the solver
void CRBC_load_face_data (YeeCrbcUpdater *U,
                            const CRBC_Side_t side,
                            const CRBC_Fields_t comp,
                            const IndexType P[3],
                            const ScalarType *val);

/// \brief return the updated boundary values to the solver
ScalarType CRBC_get_new_face_vals (YeeCrbcUpdater *U,
                                     const CRBC_Side_t side,
                                     const CRBC_Fields_t comp,
                                     const IndexType P[3]);

/// \brief compute the new boundary values
int CRBC_compute_updates (YeeCrbcUpdater *U);

/// \brief get the reflection coefficient for the given side
double CRBC_get_reflection_coef (YeeCrbcUpdater *U, const CRBC_Side_t side);

/// \brief get the maximum reflection coefficient
double CRBC_get_max_reflection_coef (YeeCrbcUpdater *U);

/// \brief get the the number of recursions on a side
int CRBC_get_num_recursions (YeeCrbcUpdater *U, const CRBC_Side_t side);


/// \brief get the components we expect to be inputted into the updater
void CRBC_get_component_inputs (YeeCrbcUpdater *U, 
                           const CRBC_Side_t side, 
                           CRBC_Fields_t *comp, 
                           int *n);

/// \brief get the index bounds we expect for input
void CRBC_get_input_extents (YeeCrbcUpdater *U, 
                       const CRBC_Side_t side, 
                       const CRBC_Fields_t comp, 
                       IndexType low[3], 
                       IndexType high[3]);

/// \brief get the index bounds we can output
void CRBC_get_output_extents (YeeCrbcUpdater *U,
                         const CRBC_Side_t side, 
                         const CRBC_Fields_t comp, 
                         IndexType low[3], 
                         IndexType high[3]);

/// \brief return the wave speed
CoefType CRBC_get_c (YeeCrbcUpdater *U);

/// \brief return the final simulation time
CoefType CRBC_get_T (YeeCrbcUpdater *U);

/// \brief return the time step size
CoefType CRBC_get_dt (YeeCrbcUpdater *U); 

/// \brief return the number of faces for a specific component
int CRBC_get_num_faces (YeeCrbcUpdater *U, const CRBC_Fields_t comp);

/// \brief return the number of edges for a specific component
int CRBC_get_num_edges (YeeCrbcUpdater *U, const CRBC_Fields_t comp);

/// \brief return the number of corners for a specific component
int CRBC_get_num_corners (YeeCrbcUpdater *U, const CRBC_Fields_t comp);


#if USE_HDF5_RESTARTS
/// \brief load data from a previous save point
int CRBC_restart (YeeCrbcUpdater **U, const char *fname);

/// \brief save the current state so the simulation may be restarted
int CRBC_save_state (YeeCrbcUpdater *U, const char *fname);
#endif

#ifdef __cplusplus
}
#endif

#endif // 3D_YEE_CRBC_API_H
