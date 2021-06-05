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

/* definitions for the c++ yee interface class for the CRBC updates in 3d */

#include "3d_yee_api.hpp"
#include <iostream>

#if USE_HDF5_RESTARTS
  #include "util/hdf5_interface.hpp"
#endif

namespace yee_crbc
{

/// Constructor --- initialize only the parameters that are common to all of
/// the faces in all circumstances.
/// \param[in] T              Simulation length (time)
/// \param[in] h              Grid spacings, e.g. hx, hy, hz
/// \param[in] dt             Time step size
/// \param[in] c              Wave Speed
/// \param[in] boundaries     Boundary conditions
YeeCrbcUpdates3d::YeeCrbcUpdates3d(const CoefType &T,
                                   const CoefType h[3],  
                                   const CoefType &dt,                  
                                   const CoefType &c,
                                   const Boundaries boundaries[6])
{

  typedef crbc::CrbcUpdates<3, ScalarType, IndexType, CoefType> CRBC;

  int i;

  useTol = false;
  Pmax = 5;
  tol = -1;

  // copy in the boundaries
  for (i=0; i<6; i++)
    input_boundaries[i] = boundaries[i];

  // figure out the boundaries for each field component
  compute_internal_bounds();

  // initialize the boundary updater for each of the require components
  for (i=0; i<3; i++) {
    if (require_component[i])
      updater[i] = CRBC (T, h, dt, c, internal_bounds[i]);
  }
  initialized = true;
}

/// Constructor --- initialize only the parameters that are common to all of
/// the faces and set the number of recursions to be P in all directions
/// \param[in] T              Simulation length (time)
/// \param[in] h              Grid spacings, e.g. hx, hy, hz
/// \param[in] dt             Time step size
/// \param[in] c              Wave Speed
/// \param[in] boundaries     Boundary conditions
/// \param[in] P              Number of recursions to use
YeeCrbcUpdates3d::YeeCrbcUpdates3d(const CoefType &T,
                                   const CoefType h[3],  
                                   const CoefType &dt,                  
                                   const CoefType &c,
                                   const Boundaries boundaries[6],
                                   const int &P)
{

  typedef crbc::CrbcUpdates<3, ScalarType, IndexType, CoefType> CRBC;

  int i;

  useTol = false;
  Pmax = P;
  tol = -1;

  // copy in the boundaries
  for (i=0; i<6; i++)
    input_boundaries[i] = boundaries[i];

  // figure out the boundaries for each field component
  compute_internal_bounds();

  // initialize the boundary updater for each of the require components
  for (i=0; i<3; i++) {
    if (require_component[i])
      updater[i] = CRBC (T, h, dt, c, internal_bounds[i], P);
  }
  initialized = true;
}

/// Constructor --- initialize only the parameters that are common to all of
/// the faces and set the number of recursions to be P in all directions
/// \param[in] T              Simulation length (time)
/// \param[in] h              Grid spacings, e.g. hx, hy, hz
/// \param[in] dt             Time step size
/// \param[in] c              Wave Speed
/// \param[in] boundaries     Boundary conditions
/// \param[in] Pmax           Maximum number of recursions to use
/// \param[in] tol            Tolerance for the recursions. This is used as a 
///                           bound for the reflection coefficient.
YeeCrbcUpdates3d::YeeCrbcUpdates3d(const CoefType &T,
                                   const CoefType h[3],  
                                   const CoefType &dt,                  
                                   const CoefType &c,
                                   const Boundaries boundaries[6],
                                   const int &Pmax,
                                   const double &tol)
{

  typedef crbc::CrbcUpdates<3, ScalarType, IndexType, CoefType> CRBC;

  int i;

  useTol = true;
  this->tol = tol;
  this->Pmax = Pmax;

  // copy in the boundaries
  for (i=0; i<6; i++)
    input_boundaries[i] = boundaries[i];

  // figure out the boundaries for each field component
  compute_internal_bounds();

  // initialize the boundary updater for each of the require components
  for (i=0; i<3; i++) {
    if (require_component[i])
      updater[i] = CRBC (T, h, dt, c, internal_bounds[i]);
  }
  initialized = true;
}

/// \brief Initialize face
/// 
/// Routine to initialize parameters on each face that might allow a solver
/// to access the new data using it's native (physical) indexing. We assume
/// the solution is initially compactly supported away from the boundaries
/// so we initialize the data to 0.
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
int YeeCrbcUpdates3d::init_face(const Side &side,
                  const Fields &comp, 
                  const IndexType low_index[3],
                  const IndexType high_index[3],
                  const double &delta)
{
  int flag = 0;

  // we don't currently throw any exceptions (\todo we should) however, 
  // memory is being allocated here so it is possible that an exception could
  // be thrown. We want to make sure we notice this ...
  try {
    if (useTol) {
      updater[comp].init_face(side, low_index, high_index, delta, Pmax, tol);
    } else {
      updater[comp].init_face(side, low_index, high_index, delta);
    }
  } catch(...) {
    std::cerr << "Calling init_face failed" << std::endl;
    flag = -1;
  }
  return flag;
}

/// \brief routine to load data into faces
///
/// \param[in] side  The boundary side to load the data in
/// \param[in] comp  The field component being initialized on this face
/// \param[in] P     The physical (solver) coordinate to load the value into
/// \param[in] val   The new value from the solver at point P
void YeeCrbcUpdates3d::load_face_data(const Side &side,
                               const Fields &comp,
                               const IndexType P[3],
                               const ScalarType &val)
{
  updater[comp].load_face_data(side, P, val);
}

/// \brief routine to get updated values from the specified face
///
/// \param[in] side  The boundary side to get data from
/// \param[in] comp  The field component being initialized on this face
/// \param[in] P     The physical (solver) coordinate to get the data at
///
/// \return   The new value at the coordinate P
ScalarType YeeCrbcUpdates3d::get_new_face_vals(const Side &side,
                                        const Fields &comp,
                                        const IndexType P[3]) const
{
  return updater[comp].get_new_face_vals(side, P);
}

/// function to run the updates. This is the function that should be called
/// after all of the new values have been loaded from the solver.
int YeeCrbcUpdates3d::compute_updates()
{
  int i, flag = 0;

  // again, we're not throwing any errors ....
  try {
    for (i=0; i<3; i++) {
      if (require_component[i])
        updater[i].compute_updates();
    }
  } catch (...) {
    std::cerr << "Error calling compute_updates" << std::endl;
    flag = -1;
  }
  return flag;
}

/// \brief routine to get the reflection coefficient for the provided side
///
/// \param[in] side  The boundary side to return a reflection coefficient 
///                  for.
/// \return  The reflection coefficient for the given side. If the side is
///          not of CRBC type or the cosine parameters have not yet been
///          computed, this returns -1
double YeeCrbcUpdates3d::get_reflection_coef(const Side &side) const
{
  return updater[update_index].get_reflection_coef(side);
}

/// \brief routine to get the maximum reflection coefficient
///
/// \return  The maximum reflection coefficient over all sides. If no side
///          is of CRBC type or the cosine parameters have not yet been
///          computed, this returns -1
double YeeCrbcUpdates3d::get_max_reflection_coef() const 
{
  return updater[update_index].get_max_reflection_coef();
}

/// routine to get the number of recursions used on the given side
/// \param[in] side The boundary side to return the number of recursions for
/// 
/// \return  The number of recursions used on the provided side. If the side
///          is not of CRBC type, the default number of recursions (5) or
///          the number of recursions provided as a constructor argument
///          will be returned even though recursions will not be used ...
int YeeCrbcUpdates3d::get_num_recursions(const Side &side) const
{
  return updater[update_index].get_num_recursions(side);
}

/// \brief routine to return the field components we expect to be inputs on the
///        given side
///
/// \param[in]  side   The boundary side
/// \param[out] comp   The field component(s) we expect to be inputs. This should
///                    be allocated to hold at least 3 ints.
/// \param[out] n      The number of field components we expect as input
void YeeCrbcUpdates3d::get_component_inputs(const Side &side, Fields *comp, int &n) const
{
  int i; 
  n = 0;

  for (i = 0; i<3; i++) {
    if ((require_component[i]) && (input_boundaries[side] == CRBC)) {
      comp[n] = static_cast<Fields>(i);
      ++n;
    }
  }
}

/// \brief get the indexing bounds we expect for data input on the given side
///
/// \param[in] side The boundary side to get indexing extents for
/// \param[in] comp  The field component being initialized on this face
/// \param[out] low array containing the lower index bounds in each direction
/// \param[out] high array containing the upper index bounds in each direction
void YeeCrbcUpdates3d::get_input_extents(const Side &side, const Fields &comp, 
                                     IndexType low[3], IndexType high[3]) const
{
  updater[comp].get_input_extents(side, low, high);
}

/// \brief  get the indexing bounds we can generate output data for on the given side
///
/// \param[in] side The boundary side to get indexing extents for
/// \param[in] comp  The field component being initialized on this face
/// \param[out] low array containing the lower index bounds in each direction
/// \param[out] high array containing the upper index bounds in each direction
void YeeCrbcUpdates3d::get_output_extents(const Side &side, const Fields &comp, 
                                     IndexType low[3], IndexType high[3]) const
{
  updater[comp].get_output_extents(side, low, high);
}

/// return the wave speed
/// \return the wave speed
CoefType YeeCrbcUpdates3d::get_c() const
{
  return updater[update_index].get_c();
}

/// return the simulation time
/// \return the final time
CoefType YeeCrbcUpdates3d::get_T() const
{
  return updater[update_index].get_T();
}

/// get the time step size
/// \return Time step size
CoefType YeeCrbcUpdates3d::get_dt() const
{
  return updater[update_index].get_dt();
}

/// get the number of faces
///
/// \param[in] comp  The field component to get number of faces for
///
/// \return The number of faces for the given component
int YeeCrbcUpdates3d::get_num_faces(const Fields &comp) const
{
  return updater[comp].get_num_faces();
}

/// get the number of edges
///
/// \param[in] comp  The field component to get number of edges for
///
/// \return The number of edges for the given component
int YeeCrbcUpdates3d::get_num_edges(const Fields &comp) const
{
  return updater[comp].get_num_edges();
}

/// get the number of corners
///
/// \param[in] comp  The field component to get number of corners for
///
/// \return The number of corners for the given component
int YeeCrbcUpdates3d::get_num_corners(const Fields &comp) const
{
  return updater[comp].get_num_corners();
}


/**
  \brief Function to compute the actual boundary conditions for each component
         base on the boundary conditions for the whole problem

  To do this, we are making the assumption that the yee cells are set up such
  that the E-field components are on the edges of the cell. It is possible
  instead to define the cell so the H-field components are on the edge of the 
  cell and the E-field components are at the centers of the faces.

  At this time, we only have support for the configuration where the E-fields
  are located on the edges and the H-field components are at the centers of the
  faces. Furthermore, we require that the computational domain be comprised of 
  full Yee cells (again, it is possible in some circumstances to use half a cell
  at the boundaries)

  \todo Add support for the Yee cells with the H-field components on the edges.
        In general this should work now by just calling Hx --> Ex, but this 
        becomes confusing ...

*/
void YeeCrbcUpdates3d::compute_internal_bounds()
{

  int i, j, direction;

  // set default requirements
  for (i=0; i<3; i++)
    require_component[i] = false;
  update_index = -1;

  // default the boundaries to be of type none
  for (i=0; i<6; i++) {
    for (j=0; j<3; j++) {
      internal_bounds[j][i] = crbc::BoundaryProperties::NONE;
    }
  }

  // loop over the input boundaries
  for (i=0; i<6; i++) {

    // get the normal direction
    direction = i / 2;

    if (input_boundaries[i] == CRBC) {

      // set all the boundaries to be CRBC even if we won't need one of the
      // components
    
      internal_bounds[0][i] = crbc::BoundaryProperties::CRBC;
      internal_bounds[1][i] = crbc::BoundaryProperties::CRBC;
      internal_bounds[2][i] = crbc::BoundaryProperties::CRBC;

      // note that we require updates for the tangential components
      for (j=0; j<3; j++) {
        if (j != direction)
          require_component[j] = true;
      }

    } else if (input_boundaries[i] == PEC) {

      // The normal component is half a grid point removed from the face so
      // the correct boundary type is 0 Neumann
      internal_bounds[direction][i] = crbc::BoundaryProperties::NEUM;

      // the tangential faces are 0 Dirichlet
      for (j=0; j<3; j++) {
        if (j != direction)
          internal_bounds[j][i] = crbc::BoundaryProperties::DIR;
      }
    }

  }

  // find a component that will be used and store it's index so we can use it 
  // as a reference later
  for (j=0; j<3; j++) {
    if (require_component[j]) {
      update_index = j;
      break;
    }
  }
}

#if USE_HDF5_RESTARTS
/**
  \brief Restart the boundaries from a previously saved state.

  This will look for a file called <fname>.h5 which should contain all the
  information needed to reinitialize this class data. Additionally, it will
  load the state information for each of the require updater components by
  looking for the appropriate files (at least one of <fname>_ex.h5, 
  <fname>_ey.h5, or <fname>_ez.h5)

*/
int YeeCrbcUpdates3d::restart (const std::string &fname)
{
  typedef crbc::CrbcUpdates<3, ScalarType, IndexType, CoefType> CRBC;

  int tmp_req_comp[3], tmp_useTol, tmp_in_bounds[6], tmp_bounds[18];
  std::string dsname, aname, gname;

  // create a data set to hold top level book keeping stuff
  gname = "/yee_3d_data";
  util_hdf5::ReadHDF5 hdf5_reader(fname);

  // read the attributes
  aname = "require_component";
  hdf5_reader.getAttributeFromGroup(gname, aname, tmp_req_comp);
  aname = "input_boundaries";
  hdf5_reader.getAttributeFromGroup(gname, aname, tmp_in_bounds);
  aname = "internal_boundaries";
  hdf5_reader.getAttributeFromGroup(gname, aname, tmp_bounds);
  aname = "update_index";
  hdf5_reader.getAttributeFromGroup(gname, aname, &update_index);
  aname = "useTol";
  hdf5_reader.getAttributeFromGroup(gname, aname, &tmp_useTol);
  aname = "tol";
  hdf5_reader.getAttributeFromGroup(gname, aname, &tol);
  aname = "Pmax";
  hdf5_reader.getAttributeFromGroup(gname, aname, &Pmax);

  // convert ints to boundary and bool types
  for (int i=0; i<3; i++) {
    require_component[i] = (tmp_req_comp[i] == 0) ? false : true;
    for (int j=0; j<6; j++)
      internal_bounds[i][j] = static_cast<crbc::BoundaryProperties::Boundary>(tmp_bounds[j + 6*i]);
  }
  for (int j=0; j<6; j++)
    input_boundaries[j] = static_cast<Boundaries>(tmp_in_bounds[j]);
  useTol = (tmp_useTol == 0) ? false : true;


  std::vector< std::string> comp_names;
  comp_names.push_back("_ex");
  comp_names.push_back("_ey");
  comp_names.push_back("_ez");

  // now load the data from each of the required updater objects
  for (int i=0; i<3; i++)
    if (require_component[i])
      updater[i] = CRBC(fname + comp_names.at(i));
  
  return 0;
}


/**
  \brief Save the current state.

  This will create a file called <fname>.h5 containing all the
  information needed to reinitialize this class data. Additionally, it will
  save the state information for each of the require updater components in
  the  files <fname>_ex.h5, <fname>_ey.h5, or <fname>_ez.h5 as required.

*/

int YeeCrbcUpdates3d::save_state (const std::string &fname) const
{

  int tmp_req_comp[3], tmp_useTol, tmp_in_bounds[6], tmp_bounds[18];
  std::string dsname, aname, gname;

  if (initialized) {

    // create an object to write HDF5 files
    util_hdf5::OutputHDF5 hdf5_writer(fname);

    // create a data set to hold top level book keeping stuff
    gname = "/yee_3d_data";
    hdf5_writer.addGroup(gname);

    // convert boundary and bool types to ints
    for (int i=0; i<3; i++) {
      tmp_req_comp[i] = (require_component[i]) ? 1 : 0;
      for (int j=0; j<6; j++)
        tmp_bounds[j + 6*i] = static_cast<int>(internal_bounds[i][j]);
    }
    for (int j=0; j<6; j++)
      tmp_in_bounds[j] = static_cast<int>(input_boundaries[j]);
    tmp_useTol = (useTol) ? 1 : 0;

    // save the attributes
    aname = "require_component";
    hdf5_writer.addAttributeToGroup(gname, aname, tmp_req_comp, 3);
    aname = "input_boundaries";
    hdf5_writer.addAttributeToGroup(gname, aname, tmp_in_bounds, 6);
    aname = "internal_boundaries";
    hdf5_writer.addAttributeToGroup(gname, aname, tmp_bounds, 18);
    aname = "update_index";
    hdf5_writer.addAttributeToGroup(gname, aname, &update_index, 1);
    aname = "useTol";
    hdf5_writer.addAttributeToGroup(gname, aname, &tmp_useTol, 1);
    aname = "tol";
    hdf5_writer.addAttributeToGroup(gname, aname, &tol, 1);
    aname = "Pmax";
    hdf5_writer.addAttributeToGroup(gname, aname, &Pmax, 1);


    std::vector< std::string> comp_names;
    comp_names.push_back("_ex");
    comp_names.push_back("_ey");
    comp_names.push_back("_ez");
    // now write the data from each of the required updater objects
    for (int i=0; i<3; i++)
      if (require_component[i])
        updater[i].save_state(fname + comp_names.at(i));

    return 0;

  } else {

    std::cerr << "Error in save_state(): there is nothing to save" << std::endl;
    return -1;

  }
  
}
#endif

} // end namespace



