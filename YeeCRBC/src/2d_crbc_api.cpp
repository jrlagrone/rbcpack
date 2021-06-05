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

/* definitions for the c interface class for the CRBC updates in 2d */

#include "2d_crbc_api.h"
#include "crbc_updates.hpp"
#include "boundary_properties.hpp"
#include <iostream>

// rename crbc::CrbcUpdates<2, ScalarType, IndexType> to something shorter
typedef crbc::CrbcUpdates<2, ScalarType, IndexType, CoefType> cxx_updater;

// rename crbc::BoundaryProperties to something shorter
typedef crbc::BoundaryProperties cxx_bounds;


// simple function to convert between the input boundary types and those
// that the underlying interface expects
cxx_bounds::Boundary convert_bt(const CRBC2d_Boundaries_t &in)
{

  cxx_bounds::Boundary out = cxx_bounds::DIR;

  switch (in)
  {
    case (CRBC2d_DIR):
      out = cxx_bounds::DIR;
      break;
    case (CRBC2d_NEUM):
      out = cxx_bounds::NEUM;
      break;
    case (CRBC2d_CRBC):
      out = cxx_bounds::CRBC;
      break;
    case (CRBC2d_NONE):
      out = cxx_bounds::NONE;
      break;
  }

  return out;

}


/// initialize only the parameters that are common to all of
/// the faces in all circumstances.
///
/// \param[in] T              Simulation length (time)
/// \param[in] h              Grid spacings, e.g. hx, hy, hz
/// \param[in] dt             Time step size
/// \param[in] c              Wave Speed
/// \param[in] boundaries     Boundary conditions
///
/// \return pointer to a CRBC boundary updater object
CrbcUpdater2d* CRBC2d_new_updater (const CoefType T,
                             const CoefType h[2],  
                             const CoefType dt,                  
                             const CoefType c,
                             const CRBC2d_Boundaries_t boundaries[4]) 
{
  crbc::BoundaryProperties::Boundary tmp_bounds[4];
  for (int i =0; i<4; i++)
    tmp_bounds[i] = convert_bt(boundaries[i]);

  return reinterpret_cast<CrbcUpdater2d*> (new cxx_updater (T, h, dt, c, tmp_bounds));
}

/// initialize only the parameters that are common to all of
/// the faces and set the number of recursions to be P in all directions
///
/// \param[in] T              Simulation length (time)
/// \param[in] h              Grid spacings, e.g. hx, hy
/// \param[in] dt             Time step size
/// \param[in] c              Wave Speed
/// \param[in] boundaries     Boundary conditions
/// \param[in] P              Number of recursions to use
///
/// \return pointer to a CRBC boundary updater object
CrbcUpdater2d* CRBC2d_new_updater_p (const CoefType T,
                               const CoefType h[2],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC2d_Boundaries_t boundaries[4],
                               const int P) 
{
  crbc::BoundaryProperties::Boundary tmp_bounds[4];
  for (int i =0; i<4; i++)
    tmp_bounds[i] = convert_bt(boundaries[i]);

  return reinterpret_cast<CrbcUpdater2d*> (new cxx_updater (T, h, dt, c, tmp_bounds, P));
}

/// initialize only the parameters that are common to all of
/// the faces and set a maximum number or recursions and a tolerance
///
/// \param[in] T              Simulation length (time)
/// \param[in] h              Grid spacings, e.g. hx, hy, hz
/// \param[in] dt             Time step size
/// \param[in] c              Wave Speed
/// \param[in] boundaries     Boundary conditions
/// \param[in] Pmax           Maximum number of recursions to use
/// \param[in] tol            Tolarence for the recursions. This is used to
///                           bound the reflection coefficient.
///
/// \return pointer to a CRBC boundary updater object
CrbcUpdater2d* CRBC2d_new_updater_tol (const CoefType T,
                               const CoefType h[2],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC2d_Boundaries_t boundaries[4],
                               const int Pmax,
                               const double tol) 
{
  crbc::BoundaryProperties::Boundary tmp_bounds[4];
  for (int i =0; i<4; i++)
    tmp_bounds[i] = convert_bt(boundaries[i]);

  return reinterpret_cast<CrbcUpdater2d*> (new cxx_updater (T, h, dt, c, tmp_bounds, Pmax, tol));
}

/// Deallocate a CRBC boundary updater object
/// \param[inout] U CRBC boundary updater object
void CRBC2d_delete_updater (CrbcUpdater2d *U)
{
  delete reinterpret_cast<cxx_updater*>(U);
  U = NULL;
}


/// Routine to initialize parameters on each face that might allow a solver
/// to access the new data using it's native (physical) indexing. We assume
/// the solution is initially compactly supported away from the boundaries
/// so we initialize the data to 0.
///
///
/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side       The side boundary side this face is on
/// \param[in] low_index  The lower bounds for the "physical" indexing in
///                       the solver. We expect, for example, to be given 
///                       (i,j) in 2D where i is the x index, j the y.
///                       We do not need to know the logical indexing. This 
///                       should include the last layer of points that the solver
///                       is able to update on it's own.
/// \param[in] high_index The upper bounds for the "physical" indexing.
/// \param[in] delta      The minimum superation from this boundary face and
///                       any sources, scatterers, or other inhomogeneities
///
/// \return Error flag, 0 is success all others are failures.
int CRBC2d_init_face (CrbcUpdater2d *U,
               const CRBC2d_Side_t side,
               const IndexType low_index[2],
               const IndexType high_index[2],
               const double delta)
{

  try {
    reinterpret_cast<cxx_updater*>(U)->init_face(side, low_index, high_index, delta);
  } catch (...) {
    return -1;
  }

  return 0;
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side  The boundary side to load the data in
/// \param[in] P     The physical (solver) coordinate to load the value into
/// \param[in] val   The new value from the solver at point P
void CRBC2d_load_face_data (CrbcUpdater2d *U,
                            const CRBC2d_Side_t side,
                            const IndexType P[2],
                            const ScalarType *val)
{
  reinterpret_cast<cxx_updater*>(U)->load_face_data(side, P, *val);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side  The boundary side to get data from
/// \param[in] P     The physical (solver) coordinate to get the data at
///
/// \return   The new value at the coordinate P
ScalarType CRBC2d_get_new_face_vals (CrbcUpdater2d *U,
                                     const CRBC2d_Side_t side,
                                     const IndexType P[2])
{
  return reinterpret_cast<cxx_updater*>(U)->get_new_face_vals (side, P);
}

/// function to run the updates. This is the function that should be called
/// after all of the new values have been loaded from the solver.
///
/// \param[inout] U       CRBC boundary updater object
///
/// \return Error flag, 0 is success all others are failures.
int CRBC2d_compute_updates (CrbcUpdater2d *U)
{
  try {
    reinterpret_cast<cxx_updater*>(U)->compute_updates();
  } catch (...) {
    return -1;
  }

  return 0;
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side  The boundary side to return a reflection coefficient 
///                  for.
/// \return  The reflection coefficient for the given side. If the side is
///          not of CRBC type or the cosine parameters have not yet been
///          computed, this returns -1
double CRBC2d_get_reflection_coef (CrbcUpdater2d *U, const CRBC2d_Side_t side)
{
  return reinterpret_cast<cxx_updater*>(U)->get_reflection_coef(side);
}

/// \brief routine to get the maximum reflection coefficient
///
/// \param[inout] U       CRBC boundary updater object
///
/// \return  The maximum reflection coefficient over all sides. If no side
///          is of CRBC type or the cosine parameters have not yet been
///          computed, this returns -1
double CRBC2d_get_max_reflection_coef (CrbcUpdater2d *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_max_reflection_coef();
}

/// \param[in] side The boundary side to return the number of recursions for
/// 
/// \param[inout] U       CRBC boundary updater object
///
/// \return  The number of recursions used on the provided side. If the side
///          is not of CRBC type, the default number of recursions (5) or
///          the number of recursions provided as a constructor argument
///          will be returned even though recursions will not be used ...
int CRBC2d_get_num_recursions (CrbcUpdater2d *U, const CRBC2d_Side_t side)
{
  return reinterpret_cast<cxx_updater*>(U)->get_num_recursions (side);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side The boundary side to get indexing extents for
/// \param[out] low array containing the lower index bounds in each direction
/// \param[out] high array containing the upper index bounds in each direction
void CRBC2d_get_input_extents (CrbcUpdater2d *U, 
                               const CRBC2d_Side_t side, 
                               IndexType low[2], 
                               IndexType high[2])
{
  reinterpret_cast<cxx_updater*>(U)->get_input_extents (side, low, high);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side The boundary side to get indexing extents for
/// \param[out] low array containing the lower index bounds in each direction
/// \param[out] high array containing the upper index bounds in each direction
void CRBC2d_get_output_extents (CrbcUpdater2d *U,
                         const CRBC2d_Side_t side, 
                         IndexType low[2], 
                         IndexType high[2])
{
  reinterpret_cast<cxx_updater*>(U)->get_output_extents (side, low, high);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \return the wave speed
CoefType CRBC2d_get_c (CrbcUpdater2d *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_c();
}

/// \param[inout] U       CRBC boundary updater object
///
/// \return the final time
CoefType CRBC2d_get_T (CrbcUpdater2d *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_T();
}

/// \param[inout] U       CRBC boundary updater object
///
/// \return Time step size
CoefType CRBC2d_get_dt (CrbcUpdater2d *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_dt();
}

/// \param[inout] U       CRBC boundary updater object
///
/// \return The number of faces for the given component
int CRBC2d_get_num_sides (CrbcUpdater2d *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_num_faces();
}

/// \param[inout] U       CRBC boundary updater object
/// 
/// \return The number of corners for the given component
int CRBC2d_get_num_corners (CrbcUpdater2d *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_num_edges();
}


#if USE_HDF5_RESTARTS
/// \param[inout
int CRBC2d_restart (CrbcUpdater2d **U, const char *fname)
{
  std::string name(fname);

  // first create a new object
  try {
    *U = reinterpret_cast<CrbcUpdater2d*> (new cxx_updater(fname));
  } catch (...) {
    return -1;
  }

  return 0;
}
int CRBC2d_save_state (CrbcUpdater2d *U, const char *fname)
{
  std::string name(fname);

  try {
    reinterpret_cast<cxx_updater*>(U)->save_state(name);
  } catch(...) {
    return -1;
  }

  return 0;
}
#endif


