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

/* definitions for the c yee interface class for the CRBC updates in 3d */

#include "3d_yee_api.hpp"
#include "3d_yee_crbc_api.h"
#include <iostream>

// rename yee_crbc::YeeCrbcUpdates3d to something shorter
typedef yee_crbc::YeeCrbcUpdates3d cxx_updater;


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
YeeCrbcUpdater* CRBC_new_updater (const CoefType T,
                             const CoefType h[3],  
                             const CoefType dt,                  
                             const CoefType c,
                             const CRBC_Boundaries_t boundaries[6]) 
{
  yee_crbc::Boundaries tmp_bounds[6];
  for (int i =0; i<6; i++)
    tmp_bounds[i] = static_cast<yee_crbc::Boundaries>(boundaries[i]);

  return reinterpret_cast<YeeCrbcUpdater*> (new yee_crbc::YeeCrbcUpdates3d (T, h, dt, c, tmp_bounds));
}

/// initialize only the parameters that are common to all of
/// the faces and set the number of recursions to be P in all directions
///
/// \param[in] T              Simulation length (time)
/// \param[in] h              Grid spacings, e.g. hx, hy, hz
/// \param[in] dt             Time step size
/// \param[in] c              Wave Speed
/// \param[in] boundaries     Boundary conditions
/// \param[in] P              Number of recursions to use
///
/// \return pointer to a CRBC boundary updater object
YeeCrbcUpdater* CRBC_new_updater_p (const CoefType T,
                               const CoefType h[3],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC_Boundaries_t boundaries[6],
                               const int P) 
{
  yee_crbc::Boundaries tmp_bounds[6];
  for (int i =0; i<6; i++)
    tmp_bounds[i] = static_cast<yee_crbc::Boundaries>(boundaries[i]);

  return reinterpret_cast<YeeCrbcUpdater*> (new cxx_updater (T, h, dt, c, tmp_bounds, P));
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
YeeCrbcUpdater* CRBC_new_updater_tol (const CoefType T,
                               const CoefType h[3],  
                               const CoefType dt,                  
                               const CoefType c,
                               const CRBC_Boundaries_t boundaries[6],
                               const int Pmax,
                               const double tol) 
{
  yee_crbc::Boundaries tmp_bounds[6];
  for (int i =0; i<6; i++)
    tmp_bounds[i] = static_cast<yee_crbc::Boundaries>(boundaries[i]);

  return reinterpret_cast<YeeCrbcUpdater*> (new cxx_updater (T, h, dt, c, tmp_bounds, Pmax, tol));
}

/// Deallocate a CRBC boundary updater object
/// \param[inout] U CRBC boundary updater object
void CRBC_delete_updater (YeeCrbcUpdater *U)
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
/// \param[in] comp       The field component being initialized on this face
/// \param[in] low_index  The lower bounds for the "physical" indexing in
///                       the solver. We expect, for example, to be given 
///                       (i,j,k) in 3D where i is the x index, j the y, and
///                       k the z index. We do not need to know the logical 
///                       indexing. This should include the last layer of 
///                       points that the solver is able to update on it's 
///                       own.
/// \param[in] high_index The upper bounds for the "physical" indexing.
/// \param[in] delta      The minimum superation from this boundary face and
///                       any sources, scatterers, or other inhomogeneities
///
/// \return Error flag, 0 is success all others are failures.
int CRBC_init_face (YeeCrbcUpdater *U,
               const CRBC_Side_t side,
               const CRBC_Fields_t comp, 
               const IndexType low_index[3],
               const IndexType high_index[3],
               const double delta)
{
  return reinterpret_cast<cxx_updater*>(U)->init_face(static_cast<yee_crbc::Side>(side), static_cast<yee_crbc::Fields>(comp), low_index, high_index, delta);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side  The boundary side to load the data in
/// \param[in] comp  The field component being initialized on this face
/// \param[in] P     The physical (solver) coordinate to load the value into
/// \param[in] val   The new value from the solver at point P
void CRBC_load_face_data (YeeCrbcUpdater *U,
                            const CRBC_Side_t side,
                            const CRBC_Fields_t comp,
                            const IndexType P[3],
                            const ScalarType *val)
{
  reinterpret_cast<cxx_updater*>(U)->load_face_data(static_cast<yee_crbc::Side>(side), static_cast<yee_crbc::Fields>(comp), P, *val);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side  The boundary side to get data from
/// \param[in] comp  The field component being initialized on this face
/// \param[in] P     The physical (solver) coordinate to get the data at
///
/// \return   The new value at the coordinate P
ScalarType CRBC_get_new_face_vals (YeeCrbcUpdater *U,
                                     const CRBC_Side_t side,
                                     const CRBC_Fields_t comp,
                                     const IndexType P[3])
{
  return reinterpret_cast<cxx_updater*>(U)->get_new_face_vals (static_cast<yee_crbc::Side>(side), static_cast<yee_crbc::Fields>(comp), P);
}

/// function to run the updates. This is the function that should be called
/// after all of the new values have been loaded from the solver.
///
/// \param[inout] U       CRBC boundary updater object
///
/// \return Error flag, 0 is success all others are failures.
int CRBC_compute_updates (YeeCrbcUpdater *U)
{
  return reinterpret_cast<cxx_updater*>(U)->compute_updates();
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side  The boundary side to return a reflection coefficient 
///                  for.
/// \return  The reflection coefficient for the given side. If the side is
///          not of CRBC type or the cosine parameters have not yet been
///          computed, this returns -1
double CRBC_get_reflection_coef (YeeCrbcUpdater *U, const CRBC_Side_t side)
{
  return reinterpret_cast<cxx_updater*>(U)->get_reflection_coef(static_cast<yee_crbc::Side>(side));
}

/// \brief routine to get the maximum reflection coefficient
///
/// \param[inout] U       CRBC boundary updater object
///
/// \return  The maximum reflection coefficient over all sides. If no side
///          is of CRBC type or the cosine parameters have not yet been
///          computed, this returns -1
double CRBC_get_max_reflection_coef (YeeCrbcUpdater *U)
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
int CRBC_get_num_recursions (YeeCrbcUpdater *U, const CRBC_Side_t side)
{
  return reinterpret_cast<cxx_updater*>(U)->get_num_recursions (static_cast<yee_crbc::Side>(side));
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in]  side   The boundary side
/// \param[out] comp   The field component(s) we expect to be inputs. This should
///                    be allocated to hold at least 3 ints.
/// \param[out] n      The number of field components we expect as input
void CRBC_get_component_inputs (YeeCrbcUpdater *U, 
                           const CRBC_Side_t side, 
                           CRBC_Fields_t *comp, 
                           int *n)
{
  yee_crbc::Fields tmp_fields[3];

  reinterpret_cast<cxx_updater*>(U)->get_component_inputs (static_cast<yee_crbc::Side>(side), tmp_fields, *n);

  for (int i=0; i<*n; i++)
    comp[i] = static_cast<CRBC_Fields>(tmp_fields[i]);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side The boundary side to get indexing extents for
/// \param[in] comp  The field component being initialized on this face
/// \param[out] low array containing the lower index bounds in each direction
/// \param[out] high array containing the upper index bounds in each direction
void CRBC_get_input_extents (YeeCrbcUpdater *U, 
                       const CRBC_Side_t side, 
                       const CRBC_Fields_t comp, 
                       IndexType low[3], 
                       IndexType high[3])
{
  reinterpret_cast<cxx_updater*>(U)->get_input_extents (static_cast<yee_crbc::Side>(side), static_cast<yee_crbc::Fields>(comp), low, high);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] side The boundary side to get indexing extents for
/// \param[in] comp  The field component being initialized on this face
/// \param[out] low array containing the lower index bounds in each direction
/// \param[out] high array containing the upper index bounds in each direction
void CRBC_get_output_extents (YeeCrbcUpdater *U,
                         const CRBC_Side_t side, 
                         const CRBC_Fields_t comp, 
                         IndexType low[3], 
                         IndexType high[3])
{
  reinterpret_cast<cxx_updater*>(U)->get_output_extents (static_cast<yee_crbc::Side>(side), static_cast<yee_crbc::Fields>(comp), low, high);
}

/// \param[inout] U       CRBC boundary updater object
///
/// \return the wave speed
CoefType CRBC_get_c (YeeCrbcUpdater *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_c();
}

/// \param[inout] U       CRBC boundary updater object
///
/// \return the final time
CoefType CRBC_get_T (YeeCrbcUpdater *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_T();
}

/// \param[inout] U       CRBC boundary updater object
///
/// \return Time step size
CoefType CRBC_get_dt (YeeCrbcUpdater *U)
{
  return reinterpret_cast<cxx_updater*>(U)->get_dt();
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] comp  The field component to get number of faces for
///
/// \return The number of faces for the given component
int CRBC_get_num_faces (YeeCrbcUpdater *U, const CRBC_Fields_t comp)
{
  return reinterpret_cast<cxx_updater*>(U)->get_num_faces(static_cast<yee_crbc::Fields>(comp));
}

/// \param[inout] U       CRBC boundary updater object
///
/// \param[in] comp  The field component to get number of edges for
///
/// \return The number of edges for the given component
int CRBC_get_num_edges (YeeCrbcUpdater *U, const CRBC_Fields_t comp)
{
  return reinterpret_cast<cxx_updater*>(U)->get_num_edges(static_cast<yee_crbc::Fields>(comp));
}

/// \param[inout] U       CRBC boundary updater object
/// 
/// \param[in] comp  The field component to get number of corners for
///
/// \return The number of corners for the given component
int CRBC_get_num_corners (YeeCrbcUpdater *U, const CRBC_Fields_t comp)
{
  return reinterpret_cast<cxx_updater*>(U)->get_num_corners(static_cast<yee_crbc::Fields>(comp));
}


#if USE_HDF5_RESTARTS
/// \param[inout
int CRBC_restart (YeeCrbcUpdater **U, const char *fname)
{
  std::string name(fname);

  // first create a new object or overwrite the old one by assignment
  *U = reinterpret_cast<YeeCrbcUpdater*> (new yee_crbc::YeeCrbcUpdates3d ());

  return reinterpret_cast<cxx_updater*>(*U)->restart(name);
}
int CRBC_save_state (YeeCrbcUpdater *U, const char *fname)
{
  std::string name(fname);
  return reinterpret_cast<cxx_updater*>(U)->save_state(name);
}
#endif


