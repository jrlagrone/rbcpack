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

/* C/C++ interface to the DAB/CRBC updates for the 2D wave or TE/TM simulations

*/

/// \example yee_TM.c

#ifndef TWO_D_CRBC_API_H
#define TWO_D_CRBC_API_H

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

struct CrbcUpdater2d;
typedef struct CrbcUpdater2d CrbcUpdater2d;

/// define enumerations for the valid boundary types
typedef enum CRBC2d_Boundaries
{
  CRBC2d_DIR,    /**< Homogeneous Dirichlet */
  CRBC2d_NEUM,   /**< Homogeneous Neumann */
  CRBC2d_CRBC,   /**< Complete Radiation Boundary Condition (Implemented as a DAB) */
  CRBC2d_NONE,   /**< No boundary condition */
} CRBC2d_Boundaries_t;


/// define valid sides
typedef enum CRBC2d_Side {                                                            
  CRBC2d_XLeft = 0,  /**< Left side in the x-coordinate direction */
  CRBC2d_XRight, /**< Right side in the x-coordinate direction */
  CRBC2d_YLeft,  /**< Left side in the y-coordinate direction */
  CRBC2d_YRight, /**< Right side in the y-coordinate direction */
} CRBC2d_Side_t;                                                                 

/// \brief create a new updater object with default settings
CrbcUpdater2d* CRBC2d_new_updater (const CoefType T,
                             const CoefType h[2],  
                             const CoefType dt,                  
                             const CoefType c,
                             const CRBC2d_Boundaries_t boundaries[4]);

/// \brief create a new updater object and specify the number of recursions
CrbcUpdater2d* CRBC2d_new_updater_p (const CoefType T,
                               const CoefType h[2],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC2d_Boundaries_t boundaries[4],
                               const int P);

/// \brief create a new updater object and specify the a tolerance for the 
///        recursions
CrbcUpdater2d* CRBC2d_new_updater_tol (const CoefType T,
                               const CoefType h[2],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC2d_Boundaries_t boundaries[4],
                               const int Pmax,
                               const double tol);

/// \brief delete an updater object
void CRBC2d_delete_updater (CrbcUpdater2d *U);

/// \brief initialize a boundary side
int CRBC2d_init_face (CrbcUpdater2d *U,
               const CRBC2d_Side_t side,
               const IndexType low_index[2],
               const IndexType high_index[2],
               const double delta);

/// \brief input new data into the boundary updater from the solver
void CRBC2d_load_face_data (CrbcUpdater2d *U,
                            const CRBC2d_Side_t side,
                            const IndexType P[2],
                            const ScalarType *val);

/// \brief return the updated boundary values to the solver
ScalarType CRBC2d_get_new_face_vals (CrbcUpdater2d *U,
                                     const CRBC2d_Side_t side,
                                     const IndexType P[2]);

/// \brief compute the new boundary values
int CRBC2d_compute_updates (CrbcUpdater2d *U);

/// \brief get the reflection coefficient for the given side
double CRBC2d_get_reflection_coef (CrbcUpdater2d *U, const CRBC2d_Side_t side);

/// \brief get the maximum reflection coefficient
double CRBC2d_get_max_reflection_coef (CrbcUpdater2d *U);

/// \brief get the the number of recursions on a side
int CRBC2d_get_num_recursions (CrbcUpdater2d *U, const CRBC2d_Side_t side);


/// \brief get the components we expect to be inputted into the updater
void CRBC2d_get_component_inputs (CrbcUpdater2d *U, 
                           const CRBC2d_Side_t side, 
                           int *n);

/// \brief get the index bounds we expect for input
void CRBC2d_get_input_extents (CrbcUpdater2d *U, 
                       const CRBC2d_Side_t side, 
                       IndexType low[2], 
                       IndexType high[2]);

/// \brief get the index bounds we can output
void CRBC2d_get_output_extents (CrbcUpdater2d *U,
                         const CRBC2d_Side_t side, 
                         IndexType low[2], 
                         IndexType high[2]);

/// \brief return the wave speed
CoefType CRBC2d_get_c (CrbcUpdater2d *U);

/// \brief return the final simulation time
CoefType CRBC2d_get_T (CrbcUpdater2d *U);

/// \brief return the time step size
CoefType CRBC2d_get_dt (CrbcUpdater2d *U); 

/// \brief return the number of faces for a specific component
int CRBC2d_get_num_sides (CrbcUpdater2d *U);

/// \brief return the number of corners for a specific component
int CRBC2d_get_num_corners (CrbcUpdater2d *U);


#if USE_HDF5_RESTARTS
/// \brief load data from a previous save point
int CRBC2d_restart (CrbcUpdater2d **U, const char *fname);

/// \brief save the current state so the simulation may be restarted
int CRBC2d_save_state (CrbcUpdater2d *U, const char *fname);
#endif

#ifdef __cplusplus
}
#endif

#endif // TWO_D_CRBC_API_H
