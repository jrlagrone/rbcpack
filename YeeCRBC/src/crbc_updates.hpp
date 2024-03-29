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

/** This file contains the implementations of the main routines for the CRBC/DAB
    updates.

   The backend routines are largely dimension independent; however, this
   interface is set up for DIM <= 3. It seems to be clear on how to generalize
   everything to a higher dimension, but this is not done here primarily because 
   we believe the code would be too difficult to maintain. Additionally, we feel
   that in higher dimensions, one would likely be better served with a higher
   order implementation due to memory and CPU time requirements.
*/

#ifndef CRBC_UPDATES_H_
#define CRBC_UPDATES_H_
#include "crbc_data.hpp"
#include "grid/point.hpp"
#include "grid/box.hpp"
#include "grid/grid_array.hpp"
#include "boundary_properties.hpp"
#include "wave_equations.hpp"
#include "recursions.hpp"
#include "optimal_cosines/optimal_cosines.h"

#include<vector>

#if USE_HDF5_RESTARTS
  #include "util/hdf5_interface.hpp"
#endif

#include <string>
#include <sstream>

namespace crbc {

/// \class CrbcUpdates
/// \brief Class to handle the forward and backward recurions updates
///
/// This class ties together the recursion and wave equation updates to compute
/// the values for the auxilliary variables in the Double Absorbing Boundary
/// layer. 
///
//  the template parameters are: \n
/// \tparam DIM        The number of physical dimensions for the problem \n
/// \tparam DATA_TYPE  The type of data to use. This isn't checked, but it
///                     should be a float type. Note, we calculate the cosines
///                     using double precision.
/// \tparam IndexType  The data type to use for indexing. This isn't checked, but
///                     should be an integral type. This needs to be large
///                     enough to store the number of points each dimension.
/// \tparam CoefType   The data type used for the coefficients such as time step
///                    size, wave speed, etc.
///
/// \example wave_eq.cpp
/// The usage of the underlying C++ interface is demonstrated in this simple
/// example using the scalar wave equation.
template <int DIM = 3, 
          class DataType = double, 
          class IndexType = long int,
          class CoefType = DataType>
class CrbcUpdates : virtual public BoundaryProperties,
                    protected WaveEquationUpdates<DIM, DataType, IndexType, CoefType>,
                    protected CrbcRecursions<DIM, DataType, IndexType, CoefType>

{

  public:

    /// Constructor --- initialize only the bare minimum, to fully initialize
    ///                 assignment must be used
    CrbcUpdates ();

    /// Constructor --- initialize only the parameters that are common to all of
    /// the faces in all circumstances.
    /// \param[in] T              Simulation length (time)
    /// \param[in] h              Grid spacings, e.g. hx, hy, hz
    /// \param[in] dt             Time step size
    /// \param[in] c              Wave Speed
    /// \param[in] boundaries     Boundary conditions
    CrbcUpdates (const CoefType &T,
                 const CoefType h[DIM],  
                 const CoefType &dt,                  
                 const CoefType &c,
                 const Boundary boundaries[2*DIM]);

    /// Constructor --- initialize only the parameters that are common to all of
    /// the faces and set the number of recursions to be P in all directions
    /// \param[in] T              Simulation length (time)
    /// \param[in] h              Grid spacings, e.g. hx, hy, hz
    /// \param[in] dt             Time step size
    /// \param[in] c              Wave Speed
    /// \param[in] boundaries     Boundary conditions
    /// \param[in] P              Number of recursions to use
    CrbcUpdates (const CoefType &T,
                 const CoefType h[DIM],  
                 const CoefType &dt,                  
                 const CoefType &c,
                 const Boundary boundaries[2*DIM],
                 const int &P);

    /// Constructor --- initialize only the parameters that are common to all of
    /// the faces and set the number of recursions to be calculated based on the
    /// provided tolerance
    /// \param[in] T              Simulation length (time)
    /// \param[in] h              Grid spacings, e.g. hx, hy, hz
    /// \param[in] dt             Time step size
    /// \param[in] c              Wave Speed
    /// \param[in] boundaries     Boundary conditions
    /// \param[in] Pmax           Maximum number of recursions to use
    /// \param[in] tol            Tolerance (for the reflection coefficient)
    CrbcUpdates (const CoefType &T,
                 const CoefType h[DIM],  
                 const CoefType &dt,                  
                 const CoefType &c,
                 const Boundary boundaries[2*DIM],
                 const int &Pmax,
                 const double &tol);

    /// Constructor --- Load a state from an HDF5 file.
    ///
    /// @TODO This needs a lot more error checking than is present ...
    ///
    /// \param[in] fname input filename
    #if USE_HDF5_RESTARTS
    CrbcUpdates (const std::string &fname);
    #endif

    /// destructor
    virtual ~CrbcUpdates() {};

    /// \brief Initialize face with the default recursions.
    /// 
    /// Routine to initialize parameters on each face that might allow a solver
    /// to access the new data using it's native (physical) indexing. We assume
    /// the solution is initially compactly supported away from the boundaries
    /// so we initialize the data to 0. The default type of recursions is based
    /// upon which constructor was called.
    ///
    /// \param[in] side       The side boundary side this face is on
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
    void init_face(const int &side,
                   const IndexType low_index[DIM],
                   const IndexType high_index[DIM],
                   const double &delta);

    /// \brief Initialize face with specified number of recursions
    /// 
    /// Routine to initialize parameters on each face that might allow a solver
    /// to access the new data using it's native (physical) indexing. We assume
    /// the solution is initially compactly supported away from the boundaries
    /// so we initialize the data to 0.
    ///
    /// \param[in] side       The side boundary side this face is on
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
    /// \param[in] P          The number of recursions to use on this face.
    void init_face(const int &side,
                   const IndexType low_index[DIM],
                   const IndexType high_index[DIM],
                   const double &delta,
                   const int &P);

    /// \brief Initialize face using reflection tolerance
    /// 
    /// Routine to initialize parameters on each face that might allow a solver
    /// to access the new data using it's native (physical) indexing. We assume
    /// the solution is initially compactly supported away from the boundaries
    /// so we initialize the data to 0.
    ///
    /// \param[in] side       The side boundary side this face is on
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
    /// \param[in] Pmax       The maximum number of recursions to use on this 
    ///                       face.
    /// \param[in] tol        The desired error tolerance for this face
    void init_face(const int &side,
                   const IndexType low_index[DIM],
                   const IndexType high_index[DIM],
                   const double &delta,
                   const int &Pmax,
                   const double &tol);

    /// \brief Method to save the current state to an HDF5 file so the simulation
    /// can be restarted
    ///
    /// \param[in] file name (this will be overwritten if it exists)
    #if USE_HDF5_RESTARTS
    void save_state(const std::string &fname) const;
    #endif

    /// \brief routine to load data into faces
    ///
    /// \param[in] side  The boundary side to load the data in
    /// \param[in] P     The physical (solver) coordinate to load the value into
    /// \param[in] val   The new value from the solver at point P
    void load_face_data(const int &side,
                        const IndexType P[DIM],
                        const DataType &val);

    /// \brief routine to get updated values from the specified face
    ///
    /// \param[in] side  The boundary side to get data from
    /// \param[in] P     The physical (solver) coordinate to get the data at
    ///
    /// \return   The new value at the coordinate P
    inline DataType get_new_face_vals(const int &side,
                        const IndexType P[DIM]) const;

    /// \brief routine to get the array of point values from the specified face
    ///        at the requested recursion level of the auxiliary variables
    ///
    /// \param[in] side  The boundary side to get data from
    /// \param[in] p     The recursion level to return the data array for
    ///
    /// \return pointer to array containing all of the auxiliary variables
    ///         at the requested side and recursion level
    DataType* get_auxiliary_vars (const int &side, const IndexType &p);

    /// \brief routine to get the array of point values from the specified face
    ///        at the requested recursion level of the auxiliary variables
    ///
    /// \param[in] side  The boundary side to get data from
    /// \param[in] p     The recursion level to return the data array for
    ///
    /// \return const pointer to array containing all of the auxiliary variables
    ///         at the requested side and recursion level
    const DataType* get_auxiliary_vars (const int &side, const IndexType &p) const;

    /// \brief routine to get the auxilliary variables at the specified point 
    ///        and recursion level
    ///
    /// \param[in] side  The boundary side to get data from
    /// \param[in] P     The physical (solver) coordinate to get the data at
    /// \param[in] q     The recursion level to return the data array for
    ///
    /// \return aux. variable at the specified location
    inline DataType get_auxiliary_vars (const int &side, 
                                        const IndexType P[DIM], 
                                        const IndexType &q) const;

    /// \brief routine to set the auxilliary variables at the specified point 
    ///        and recursion level
    ///
    /// \param[in] side  The boundary side to get data from
    /// \param[in] P     The physical (solver) coordinate to get the data at
    /// \param[in] q     The recursion level to return the data array for
    /// \param[in] val   Value to set
    void set_auxiliary_vars (const int &side, 
                             const IndexType P[DIM], 
                             const IndexType &q,
                             const DataType &val);

    /// \brief routine to return the index of an edge by giving the index of
    ///        to intersecting faces
    ///
    /// \param[in] sidea  First side
    /// \param[in] sideb  Second side
    ///
    /// \return Index of the edge. Returns -1 if the edge is not used or does
    ///         not exist.
    int get_edge_index(const int &sidea, const int &sideb) const;
    
    /// \brief routine to get the auxilliary variables at the specified point 
    ///        and recursion level
    ///
    /// \param[in] edge  The boundart edge to get data from
    /// \param[in] P     The physical (solver) coordinate to get the data at
    /// \param[in] q     The recursion level to return the data array for
    ///
    /// \return Edge auxilliary variable at the specified location
    inline DataType get_edge_auxiliary_vars (const int &edge, 
                                        const IndexType P[DIM], 
                                        const IndexType q[2]) const;

    /// \brief routine to set the auxilliary variables at the specified point 
    ///        and recursion level
    ///
    /// \param[in] edge  The boundary side to get data from
    /// \param[in] P     The physical (solver) coordinate to get the data at
    /// \param[in] q     The recursion level to return the data array for
    /// \param[in] val   Value to set
    void set_edge_auxiliary_vars (const int &edge, 
                             const IndexType P[DIM], 
                             const IndexType q[2],
                             const DataType &val);

    /// \brief routine to return the order of the indexing for the auxiliary 
    ///        variables. 
    /// \return array of size DIM - 0 = x, 1 = y, 2 = z, etc.
    std::vector<int> get_auxiliary_order(const int &side) const;

    /// function to run the updates. This is the function that should be called
    /// after all of the new values have been loaded from the solver.
    void compute_updates();

    
    /// \brief routine to get the reflection coefficient for the provided side
    ///
    /// \param[in] side  The boundary side to return a reflection coefficient 
    ///                  for.
    /// \return  The reflection coefficient for the given side. If the side is
    ///          not of CRBC type or the cosine parameters have not yet been
    ///          computed, this returns -1
    double get_reflection_coef(const int &side) const {return ref_coefs[side];};

    /// \brief routine to get the maximum reflection coefficient
    ///
    /// \return  The maximum reflection coefficient over all sides. If no side
    ///          is of CRBC type or the cosine parameters have not yet been
    ///          computed, this returns -1
    double get_max_reflection_coef() const;

    /// routine to get the number of recursions used on the given side
    /// \param[in] side The boundary side to return the number of recursions for
    /// 
    /// \return  The number of recursions used on the provided side. If the side
    ///          is not of CRBC type, the default number of recursions (5) or
    ///          the number of recursions provided as a constructor argument
    ///          will be returned even though recursions will not be used ...
    int get_num_recursions(const int &side) const {return num_recursions[side];};

    /// \brief get the indexing bounds we expect for data input on the given side
    ///
    /// \param[in] side The boundary side to get indexing extents for
    /// \param[out] low array containing the lower index bounds in each direction
    /// \param[out] high array containing the upper index bounds in each direction
    void get_input_extents(const int &side, IndexType low[DIM], IndexType high[DIM]) const;

    /// \brief  get the indexing bounds we can generate output data for on the given side
    ///
    /// \param[in] side The boundary side to get indexing extents for
    /// \param[out] low array containing the lower index bounds in each direction
    /// \param[out] high array containing the upper index bounds in each direction
    void get_output_extents(const int &side, IndexType low[DIM], IndexType high[DIM]) const;

    /// \brief  get the indexing bounds we can generate output data for on the given edge
    ///
    /// \param[in] edge The boundary edge to get indexing extents for
    /// \param[out] low array containing the lower index bounds in each direction
    /// \param[out] high array containing the upper index bounds in each direction
    /// \param[out] plow array containing the lower index bounds for recursion levels
    /// \param[out] phigh array containing the upper index bounds for recursion levels
    void get_edge_extents(const int &edge, 
                          IndexType low[DIM], 
                          IndexType high[DIM],
                          IndexType plow[2],
                          IndexType phigh[2]) const;

    /// return the wave speed
    CoefType get_c() const {return c;};

    /// return the simulation time
    CoefType get_T() const {return T;};

    /// get the time step size
    CoefType get_dt() const {return dt;}; 

    /// get the boundary conditions
    Boundary* get_boundaries() const {return boundaries.data();};  

    /// get the number of faces
    int get_num_faces() const {return num_faces;};

    /// get the number of edges
    int get_num_edges() const {return num_edges;};
 
    /// get the number of corners
    int get_num_corners() const {return num_corners;};

    /// get approximate memory usage in MB
    double get_mem_usage() const;

  private:

    CRBC_Data<1, DIM, DataType, IndexType, CoefType> Faces[6];  //< storage objects for face data
    CRBC_Data<2, DIM, DataType, IndexType, CoefType> Edges[12]; //< storage objects for edge data
    CRBC_Data<3, DIM, DataType, IndexType, CoefType> Corners[8]; //< storage objects for corner data
    std::vector<CoefType> recurs_coefs_a[6]; //< storage for the 'a' recursion coefficients for each face
    std::vector<CoefType> recurs_coefs_ab[6]; //< storage for the '\f$\bar{a}\f$' recursion coefficients for each face
    std::vector<CoefType> recurs_coefs_sig[6]; //< storage for the '\f$\sigma\f$' recursion coefficients for each face
    std::vector<CoefType> recurs_coefs_sigb[6]; //< storage for the '\f$\var{\sigma}\f$' recursion coefficients for each face
    int num_faces; //< number of CRBC faces
    int num_edges; //< number of CRBC edges
    int num_corners; //< number of CRBC corners
    bool useTol; //< flag to indicate if the tolerance based constructor was used.
    bool have_face[6]; //< flag indicating whether we have set up the edge for each side
    bool have_new_input[6]; //< used to check if values have been updated (\todo remove? debug only?)
    bool have_cosines[6]; //< flag indicating wheter cosines were allocated or not
    bool initialized; //< flag to check if we are fully initialized
    int edge_pairs[12][2]; //< stores the indices for faces that intersect to for the edge
    int corner_triples[8][3]; //< stores the indices for edges that intersect to for the corner
    int num_recursions[6]; //< number of recursions on each face
    double ref_coefs[6]; //< maximum reflection coefficient on each face
    double delta[6]; //< face seperation from sources, scatterers, etc.
    grid::Point<2*DIM, Boundary> boundaries; //< boundary conditions
    grid::Point<DIM,CoefType> h; //< mesh widths
    CoefType dt; //< time step size
    CoefType c; //< wave speed
    CoefType T; //< Final simulation time
    int Pmax; //< maximum number of recursions when using tolerance based methods
    double tol; //< tolerance (for the reflection coefficient
    grid::Point<1,IndexType> zero; //< a zero vector --- used frequently for input/output

    /// Routine to calculate the wave equation updates. 
    /// 
    /// \tparam REC_DIM           Number of recursion directions, e.g. 1 for a 
    ///                           face, 2 for an edge, 3 for a corner
    /// \param[inout] crbc_data   A CRBC_Data array that we want to update the 
    ///                           data in using the wave equation
    template <int REC_DIM>
    void wave_equation_updates(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data);

    /// Routine to apply the recursion updates. 
    /// 
    /// \tparam REC_DIM           Number of recursion directions, e.g. 1 for a 
    ///                           face, 2 for an edge, 3 for a corner
    /// \param[inout] crbc_data   A CRBC_Data array that we want to update the 
    ///                           data in using the wave equation
    template <int REC_DIM>
    void recursion_updates(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data);

    /// routine to copy data from the face updates into the edge routines
    void copy_face_to_edge();
  
    /// routine to copy data from the edge updates into the corner routines
    void copy_edge_to_corner();

    /// routine to copy data from corner updates into the edge updates
    void copy_corner_to_edge();
   
    /// routine to copy data from the edge updates into the face updates
    void copy_edge_to_face();

    /// routine to find edges
    void find_edges();
 
    /// routine to find corners
    void find_corners();

    /// routine to step face
    void step_faces();

    /// routine to step edges
    void step_edges();

    /// routine to step corners
    void step_corners();

    /// routine to copy the current data to the old data and the new data to the
    /// current data arrays (via pointer swapping)
    void rotate_data();

    /// routine to set default values
    /// this sets all of the checks to false and the number of recursions
    /// to be 5
    void set_defaults();

    /// function to set up any edges and corners we have after the faces have 
    /// been initialized
    void initialize();

    /// function to check the cosine parameters
    void check_cos_params(double &eta, int &P);

    /// funtion to compute coefficients a, bar{a}, sigma, and bar{sigma} used 
    /// in the recursions for the given side
    void compute_crbc_coefficients(const int &side, const double *cosines);

    /// function to initialize face parameters
    void init_face_params(const int &side, 
                          const IndexType low_index[DIM], 
                          const IndexType high_index[DIM]);

    /// function to print out message is optimal cosines fails
    void print_opt_cos_info(const int &flag) const;

    /// function to get the CRBC data properties from hdf5 file
    #if USE_HDF5_RESTARTS
    template <int REC_DIM>
    void get_crbc_props(grid::Point<DIM, IndexType> &domain_low,
                     grid::Point<DIM, IndexType> &domain_high,
                     grid::Point<DIM, IndexType> &domain_order, 
                     grid::Point<REC_DIM, IndexType> &rec_low, 
                     grid::Point<REC_DIM, IndexType> &rec_high,
                     grid::Point<REC_DIM, IndexType> &rec_order,
                     grid::Point<REC_DIM, int> &normals, 
                     grid::Point<DIM*2, Boundary> &bounds,
                     util_hdf5::ReadHDF5 &hdf5_reader,
                     const std::string &gname) const;
    #endif
    

    /// function to convert integers to strings ... This is a in the C++11 
    /// standard, but we don't necessarily want to have to require c++11 just
    /// to get this function
    template <typename T> 
    std::string to_string( const T& n ) const
    {
        std::ostringstream s;
        s << n ;
        return s.str() ;
    }
    

}; // end CrbcUpdates class

//------------------------------------------------------------------------------
//                                 Definitions
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//                                 Contructors
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
CrbcUpdates<DIM, DataType, IndexType, CoefType>::CrbcUpdates()
{
  // set the defaults
  set_defaults();
}

template <int DIM, class DataType, class IndexType, class CoefType>
CrbcUpdates<DIM, DataType, IndexType, CoefType>::CrbcUpdates(const CoefType &T, 
                 const CoefType h[DIM], 
                 const CoefType &dt,      
                 const CoefType &c,  
                 const Boundary boundaries[2*DIM]) : 
                 WaveEquationUpdates<DIM, DataType, IndexType, CoefType>(h, c, dt)
{

  // set the defaults
  set_defaults();

  // save the inputs
  this->T = T;
  this->h = grid::Point<DIM,CoefType>(h);
  this->dt = dt;
  this->c = c;
  this->boundaries = grid::Point<2*DIM,Boundary>(boundaries);

  // use the boundaries to figure out where we need to compute updates
  for (int i=0; i<2*DIM; i++) {
    if (boundaries[i] == CRBC) {
      have_face[i] = true;
      ++num_faces;
    }
  }
  find_edges();
  find_corners();

}

template <int DIM, class DataType, class IndexType, class CoefType>
CrbcUpdates<DIM, DataType, IndexType, CoefType>::CrbcUpdates(const CoefType &T, 
                 const CoefType h[DIM], 
                 const CoefType &dt,      
                 const CoefType &c,  
                 const Boundary boundaries[2*DIM],
                 const int &P)  : 
                 WaveEquationUpdates<DIM, DataType, IndexType, CoefType>(h, c, dt)
{

  // set the defaults
  set_defaults();

  // save the inputs
  this->T = T;
  this->h = grid::Point<DIM,CoefType>(h);
  this->dt = dt;
  this->c = c;
  this->boundaries = grid::Point<2*DIM,Boundary>(boundaries);

  // save the number of recursions
  for (int i=0; i<6; i++)
    num_recursions[i] = P;

  // use the boundaries to figure out where we need to compute updates
  for (int i=0; i<6; i++) {
    if (boundaries[i] == CRBC) {
      have_face[i] = true;
      ++num_faces;
    }
  }
  find_edges();
  find_corners();

}

template <int DIM, class DataType, class IndexType, class CoefType>
CrbcUpdates<DIM, DataType, IndexType, CoefType>::CrbcUpdates (const CoefType &T,
                 const CoefType h[DIM],  
                 const CoefType &dt,                  
                 const CoefType &c,
                 const Boundary boundaries[2*DIM],
                 const int &Pmax,
                 const double &tol): 
                 WaveEquationUpdates<DIM, DataType, IndexType, CoefType>(h, c, dt)
{
  // set the defaults
  set_defaults();

  useTol = true;

  // save the inputs
  this->T = T;
  this->h = grid::Point<DIM,CoefType>(h);
  this->dt = dt;
  this->c = c;
  this->boundaries = grid::Point<2*DIM,Boundary>(boundaries);
  this->Pmax = Pmax;
  this->tol = tol;

  // use the boundaries to figure out where we need to compute updates
  for (int i=0; i<2*DIM; i++) {
    if (boundaries[i] == CRBC) {
      have_face[i] = true;
      ++num_faces;
    }
  }
  find_edges();
  find_corners();
}


//------------------------------------------------------------------------------
//                       Restart constructor ....
//
// TODO add checks to make sure this works!
//------------------------------------------------------------------------------
#if USE_HDF5_RESTARTS
template <int DIM, class DataType, class IndexType, class CoefType>
CrbcUpdates<DIM, DataType, IndexType, CoefType>::CrbcUpdates (const std::string &fname)
{
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData1d;
  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData2d;
  typedef CRBC_Data<3, DIM, DataType, IndexType, CoefType> CData3d;
  std::string dsname, aname, gname, gname2;
  int temp[6], temp2[2*DIM];

  // set the default values
  set_defaults();

  // create an object to read HDF files
  util_hdf5::ReadHDF5 hdf5_reader(fname);
  
  // first get the "book keeping" attributes
  gname = "/book_keeping";
  aname = "h";
  hdf5_reader.getAttributeFromGroup(gname, aname, h.data());
  aname = "dt";
  hdf5_reader.getAttributeFromGroup(gname, aname, &dt);
  aname = "c";
  hdf5_reader.getAttributeFromGroup(gname, aname, &c);
  aname = "T";
  hdf5_reader.getAttributeFromGroup(gname, aname, &T);
  aname = "num_faces";
  hdf5_reader.getAttributeFromGroup(gname, aname, &num_faces);
  aname = "num_edges";
  hdf5_reader.getAttributeFromGroup(gname, aname, &num_edges);
  aname = "num_corners";
  hdf5_reader.getAttributeFromGroup(gname, aname, &num_corners);
  aname = "num_recursions";
  hdf5_reader.getAttributeFromGroup(gname, aname, &(num_recursions[0]));
  aname = "delta";
  hdf5_reader.getAttributeFromGroup(gname, aname, &(delta[0]));
  aname = "boundaries";
  hdf5_reader.getAttributeFromGroup(gname, aname, &(temp2[0]));

  for (int i=0; i<2*DIM; i++)
    boundaries.set_value(i, static_cast<Boundary>(temp2[i]));

  aname = "have_face";
  hdf5_reader.getAttributeFromGroup(gname, aname, &(temp[0]));

  for (int i=0; i<6; i++) {

    if (temp[i] == 0) {
      have_face[i] = false;
    } else {
      have_face[i] = true;

      // calculate the recursion coefficients ... 
      double eta = static_cast<double>(delta[i] / (c*T));
      double emax;
      int flag;

      check_cos_params(eta, num_recursions[i]);

      // allocate memory for the cosines
      double *cosines = new double[2*num_recursions[i]];

      if (num_recursions[i] != 0) {
        flag = optimal_cosinesP(eta, num_recursions[i], cosines, &emax);
        ref_coefs[i] = emax;

        print_opt_cos_info(flag);
      }

      // calculate a, ab, sig, and sigb
      compute_crbc_coefficients(i, cosines);

      delete[] cosines;

      gname = "/Face_" + to_string(i);
      grid::Point<1, int> normals;
      grid::Point<DIM*2, Boundary> bounds;

      grid::Point<DIM, IndexType> domain_high, domain_low, domain_order;
      grid::Point<1, IndexType> rec_high, rec_low, rec_order;

      // get the indexing properties for this face ...
      get_crbc_props(domain_low, domain_high, domain_order, rec_low, rec_high, \
                     rec_order, normals, bounds, hdf5_reader, gname);

      grid::Box<DIM, IndexType> domain(domain_low, domain_high, domain_order);
      grid::Box<1, IndexType> recursions(rec_low, rec_high, rec_order);

      // initialize the face
      Faces[i] = CRBC_Data<1, DIM, DataType, IndexType, CoefType>(domain, recursions, normals, bounds); 

      // now load the cosines into the data structure
      Faces[i].calculate_coefficients(i, recurs_coefs_a[i].data(), 
                                     recurs_coefs_ab[i].data(),
                                     recurs_coefs_sig[i].data(), 
                                     recurs_coefs_sigb[i].data(),
                                     c*dt/h[i/2]);

      // create an iterator over the recursions
      grid::BoxIterator<1, IndexType> R = 
          grid::BoxIterator<1, IndexType> (Faces[i].get_rec_box());

      // read the data set for each recursion level
      for (R.begin(); R.notDone(); ++R) {
        dsname = gname + "/old/Recursion_Level_" +  to_string(R()[0]); 
        hdf5_reader.getDataSet(dsname, Faces[i](CData1d::OLD, R()).get_data());
       
        dsname = gname + "/current/Recursion_Level_" +  to_string(R()[0]); 
        hdf5_reader.getDataSet(dsname, Faces[i](CData1d::CUR, R()).get_data());

        dsname = gname + "/new/Recursion_Level_" +  to_string(R()[0]); 
        hdf5_reader.getDataSet(dsname, Faces[i](CData1d::NEW, R()).get_data());
      }
    } // end if
  } // end for over have_face ...

  // loop over the edges ...
  for (int i=0; i<num_edges; i++) {
    gname = "/Edge_" + to_string(i);
    grid::Point<2, int> normals;
    grid::Point<DIM*2, Boundary> bounds;
    
    grid::Point<DIM, IndexType> domain_high, domain_low, domain_order;
    grid::Point<2, IndexType> rec_high, rec_low, rec_order;

    // get the indexing properties for this edge ...
    get_crbc_props(domain_low, domain_high, domain_order, rec_low, rec_high, \
                     rec_order, normals, bounds, hdf5_reader, gname);

    grid::Box<DIM, IndexType> domain(domain_low, domain_high, domain_order);
    grid::Box<2, IndexType> recursions(rec_low, rec_high, rec_order);

    aname = "edge_pairs";
    hdf5_reader.getAttributeFromGroup(gname, aname, &(edge_pairs[i][0]));

    // initialize the edge
    Edges[i] = CRBC_Data<2, DIM, DataType, IndexType, CoefType>(domain, recursions, normals, bounds); 

    // now load the cosines into the data structure
    for (int j=0; j<2; j++)
      Edges[i].calculate_coefficients(normals[j], recurs_coefs_a[normals[j]].data(), 
                                     recurs_coefs_ab[normals[j]].data(),
                                     recurs_coefs_sig[normals[j]].data(), 
                                     recurs_coefs_sigb[normals[j]].data(),
                                     c*dt/h[normals[j]/2]);

    // create an iterator over the recursions
    grid::BoxIterator<2, IndexType> R = 
        grid::BoxIterator<2, IndexType> (Edges[i].get_rec_box());

    // read the data set for each recursion level
    for (R.begin(); R.notDone(); ++R) {
      dsname = gname + "/old/Recursion_Level_" + to_string(R()[0])
                     + "_" + to_string(R()[1]); 
      hdf5_reader.getDataSet(dsname, Edges[i](CData2d::OLD, R()).get_data());
     
      dsname = gname + "/current/Recursion_Level_" + to_string(R()[0])
                     + "_" + to_string(R()[1]);
      hdf5_reader.getDataSet(dsname, Edges[i](CData2d::CUR, R()).get_data());

      dsname = gname + "/new/Recursion_Level_" + to_string(R()[0])
                     + "_" + to_string(R()[1]);
      hdf5_reader.getDataSet(dsname, Edges[i](CData2d::NEW, R()).get_data());
    }
  } // loop over edges

  // loop over the corners ...
  for (int i=0; i<num_corners; i++) {
    gname = "/Corner_" + to_string(i);
    grid::Point<3, int> normals;
    grid::Point<DIM*2, Boundary> bounds;
   
    grid::Point<DIM, IndexType> domain_high, domain_low, domain_order;
    grid::Point<3, IndexType> rec_high, rec_low, rec_order;

    // get the indexing properties for this face ...
    get_crbc_props(domain_low, domain_high, domain_order, rec_low, rec_high, \
                     rec_order, normals, bounds, hdf5_reader, gname);

    grid::Box<DIM, IndexType> domain(domain_low, domain_high, domain_order);
    grid::Box<3, IndexType> recursions(rec_low, rec_high, rec_order);

    aname = "corner_triples";
    hdf5_reader.getAttributeFromGroup(gname, aname, &(corner_triples[i][0]));

    // initialize the edge
    Corners[i] = CRBC_Data<3, DIM, DataType, IndexType, CoefType>(domain, recursions, normals, bounds); 

    // now load the cosines into the data structure
    for (int j=0; j<3; j++)
      Corners[i].calculate_coefficients(normals[j], recurs_coefs_a[normals[j]].data(), 
                                     recurs_coefs_ab[normals[j]].data(),
                                     recurs_coefs_sig[normals[j]].data(), 
                                     recurs_coefs_sigb[normals[j]].data(),
                                     c*dt/h[normals[j]/2]);

    // create an iterator over the recursions
    grid::BoxIterator<3, IndexType> R = 
        grid::BoxIterator<3, IndexType> (Corners[i].get_rec_box());

    // read the data set for each recursion level
    for (R.begin(); R.notDone(); ++R) {
      dsname = gname + "/old/Recursion_Level_" + to_string(R()[0])
                     + "_" + to_string(R()[1]) + "_" + to_string(R()[2]); 
      hdf5_reader.getDataSet(dsname, Corners[i](CData3d::OLD, R()).get_data());
     
      dsname = gname + "/current/Recursion_Level_" + to_string(R()[0])
                     + "_" + to_string(R()[1]) + "_" + to_string(R()[2]); 
      hdf5_reader.getDataSet(dsname, Corners[i](CData3d::CUR, R()).get_data());

      dsname = gname + "/new/Recursion_Level_" + to_string(R()[0])
                     + "_" + to_string(R()[1]) + "_" + to_string(R()[2]);
      hdf5_reader.getDataSet(dsname, Corners[i](CData3d::NEW, R()).get_data());
    }
  } // loop over corners


  // initialize the wave equation updater stuff since we couldn't call the 
  // constructor directly like in the other constructors
  this->set_wave_speed(c);
  this->set_time_step(dt);
  this->set_grid_spacing(h);

  // finally, say we are intilized ...
  initialized = true;

}
#endif

//------------------------------------------------------------------------------
//                          Face initialization routines
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::init_face(const int &side,
                   const IndexType low_index[DIM],
                   const IndexType high_index[DIM],
                   const double &delta)
{
  if (useTol) {
    init_face(side, low_index, high_index, delta, Pmax, tol);
  } else {
    init_face(side, low_index, high_index, delta, num_recursions[side]);
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::init_face(const int &side,
                   const IndexType low_index[DIM],
                   const IndexType high_index[DIM],
                   const double &delta,
                   const int &P)
{
  int flag;
  double emax;

  // copy inputs
  num_recursions[side] = P;
  this->delta[side] = delta;

  // calculate eta
  double eta = static_cast<double> (delta / (c*T));

  check_cos_params(eta, num_recursions[side]);

  // allocate memory for the cosines
  double *cosines = new double[2*num_recursions[side]];

  // calculate cosines
  if (num_recursions[side] != 0) {
    flag = optimal_cosinesP(eta, num_recursions[side], cosines, &emax);
    ref_coefs[side] = emax;

    print_opt_cos_info(flag);
  }

  // calculate a, ab, sig, and sigb
  compute_crbc_coefficients(side, cosines);

  delete[] cosines;

  // now set up the crbc data structure ...
  init_face_params(side, low_index, high_index);

}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::init_face(const int &side,
                   const IndexType low_index[DIM],
                   const IndexType high_index[DIM],
                   const double &delta,
                   const int &Pmax,
                   const double &tol)
{
  int flag;
  double emax;
  int P, loc_Pmax;

  this->delta[side] = delta;

  // calculate eta
  double eta = static_cast<double> (delta / (c*T));

  loc_Pmax = Pmax;
  check_cos_params(eta, loc_Pmax);

  // allocate memory for the cosines
  double *cosines = new double[2*loc_Pmax];

    // calculate cosines
  if (loc_Pmax != 0) {
    flag = optimal_cosines(eta, Pmax, tol, cosines, &P, &emax);
    ref_coefs[side] = emax;

    print_opt_cos_info(flag); 
  }
  // calculate cosines

  // save the number of recursions
  num_recursions[side] = P;

  // calculate a, ab, sig, and sigb
  compute_crbc_coefficients(side, cosines);
 
  delete[] cosines;

  // now set up the crbc data structure ...
  init_face_params(side, low_index, high_index);
}

//------------------------------------------------------------------------------
//                         Routine to save the current state
//------------------------------------------------------------------------------
#if USE_HDF5_RESTARTS
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::save_state(const std::string &fname) const
{

  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData1d;
  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData2d;
  typedef CRBC_Data<3, DIM, DataType, IndexType, CoefType> CData3d;

  // there's no reason to dump any data if we haven't initialized it ...
  if (initialized) {

    std::string dsname, aname, gname, gname2;

    // create an object to write HDF5 files
    util_hdf5::OutputHDF5 hdf5_writer(fname);

    // create a data set to hold top level book keeping stuff
    gname = "/book_keeping";
    hdf5_writer.addGroup(gname);

    // save the input values and book keeping stuff as attributes 
    aname = "h";
    hdf5_writer.addAttributeToGroup(gname, aname, h.data(), 2*DIM);
    aname = "dt";
    hdf5_writer.addAttributeToGroup(gname, aname, &dt, 1);
    aname = "c";
    hdf5_writer.addAttributeToGroup(gname, aname, &c, 1);
    aname = "T";
    hdf5_writer.addAttributeToGroup(gname, aname, &T, 1);
    aname = "num_faces";
    hdf5_writer.addAttributeToGroup(gname, aname, &num_faces, 1);
    aname = "num_edges";
    hdf5_writer.addAttributeToGroup(gname, aname, &num_edges, 1);
    aname = "num_corners";
    hdf5_writer.addAttributeToGroup(gname, aname, &num_corners, 1);
    aname = "num_recursions";
    hdf5_writer.addAttributeToGroup(gname, aname, &(num_recursions[0]), 6);
    aname = "delta";
    hdf5_writer.addAttributeToGroup(gname, aname, &(delta[0]), 6);

    // convert boolean values to ints
    int hface[6];
    for (int i=0; i<6; i++)
      hface[i] = (have_face[i]) ? 1 : 0;

    aname = "have_face";
    hdf5_writer.addAttributeToGroup(gname, aname, &(hface[0]), 6);

    // change the boundary types to ints to save ...
    grid::Point<2*DIM, int> tmp_bounds;
    for (int i=0; i<2*DIM; i++)
       tmp_bounds.set_value(i, static_cast<int>(boundaries[i]));
    aname = "boundaries";
    hdf5_writer.addAttributeToGroup(gname, aname, tmp_bounds.data(), 2*DIM);

    // add a group for each of the faces we have data for and save the data
    // in them ...
    for (int i=0; i<6; i++) {
      if (have_face[i]) {
        gname = "/Face_" + to_string(i);
        hdf5_writer.addGroup(gname);

        // add the data and recursion extends as attributes
        aname = "data_extent_low";
        hdf5_writer.addAttributeToGroup(gname, aname, Faces[i].get_data_box().get_low().data(), DIM);
        aname = "data_extent_high";
        hdf5_writer.addAttributeToGroup(gname, aname, Faces[i].get_data_box().get_high().data(), DIM);
        aname = "data_ordering";
        hdf5_writer.addAttributeToGroup(gname, aname, Faces[i].get_data_box().get_ordering().data(), DIM);
        aname = "recursion_extent_low";
        hdf5_writer.addAttributeToGroup(gname, aname, Faces[i].get_rec_box().get_low().data(), 1);
        aname = "recursion_extent_high";
        hdf5_writer.addAttributeToGroup(gname, aname, Faces[i].get_rec_box().get_high().data(), 1);
        aname = "recursion_ordering";
        hdf5_writer.addAttributeToGroup(gname, aname, Faces[i].get_rec_box().get_ordering().data(), 1);
        
        // add the normal and boundary properties
        aname = "recursion_normals";
        hdf5_writer.addAttributeToGroup(gname, aname, Faces[i].get_sides().data(), 1);
        for (int j=0; j<2*DIM; j++)
          tmp_bounds.set_value(i, static_cast<int>(Faces[i].get_boundaries()[i]));
        aname = "boundaries";
        hdf5_writer.addAttributeToGroup(gname, aname, tmp_bounds.data(), 2*DIM);

        // create an iterator over the recursions
        grid::BoxIterator<1, IndexType> R = 
          grid::BoxIterator<1, IndexType> (Faces[i].get_rec_box());

        // create groups for the old, current, and new data
        gname2 = gname + "/new";
        hdf5_writer.addGroup(gname2);
        gname2 = gname + "/old";
        hdf5_writer.addGroup(gname2);
        gname2 = gname + "/current";
        hdf5_writer.addGroup(gname2);

        // Create a data set for each recursion level
        for (R.begin(); R.notDone(); ++R) {
          dsname = gname + "/old/Recursion_Level_" +  to_string(R()[0]); 
          hdf5_writer.addDataSet(dsname, \
                                Faces[i](CData1d::OLD, R()).get_data(), 
                                Faces[i](CData1d::OLD, R()).get_n());
       
          dsname = gname + "/current/Recursion_Level_" +  to_string(R()[0]); 
          hdf5_writer.addDataSet(dsname, 
                                Faces[i](CData1d::CUR, R()).get_data(), 
                                Faces[i](CData1d::CUR, R()).get_n());

          dsname = gname + "/new/Recursion_Level_" +  to_string(R()[0]); 
          hdf5_writer.addDataSet(dsname, 
                                Faces[i](CData1d::NEW, R()).get_data(), 
                                Faces[i](CData1d::NEW, R()).get_n());
        }
      }
    }
     
    // add a group for each of the edges we have data for and save the data
    // in them ...
    for (int i=0; i<num_edges; i++) {
      gname = "/Edge_" + to_string(i);
      hdf5_writer.addGroup(gname);

      // add the data and recursion extends as attributes
      aname = "data_extent_low";
      hdf5_writer.addAttributeToGroup(gname, aname, Edges[i].get_data_box().get_low().data(), DIM);
      aname = "data_extent_high";
      hdf5_writer.addAttributeToGroup(gname, aname, Edges[i].get_data_box().get_high().data(), DIM);
      aname = "data_ordering";
      hdf5_writer.addAttributeToGroup(gname, aname, Edges[i].get_data_box().get_ordering().data(), DIM);
      aname = "recursion_extent_low";
      hdf5_writer.addAttributeToGroup(gname, aname, Edges[i].get_rec_box().get_low().data(), 2);
      aname = "recursion_extent_high";
      hdf5_writer.addAttributeToGroup(gname, aname, Edges[i].get_rec_box().get_high().data(), 2);
      aname = "recursion_ordering";
      hdf5_writer.addAttributeToGroup(gname, aname, Edges[i].get_rec_box().get_ordering().data(), 2);

      // add the normal and boundary properties
      aname = "recursion_normals";
      hdf5_writer.addAttributeToGroup(gname, aname, Edges[i].get_sides().data(), 2);
      for (int j=0; j<2*DIM; j++)
        tmp_bounds.set_value(i, static_cast<int>(Edges[i].get_boundaries()[i]));
      aname = "boundaries";
      hdf5_writer.addAttributeToGroup(gname, aname, tmp_bounds.data(), 2*DIM);

      // add edge pairs
      aname = "edge_pairs";
      hdf5_writer.addAttributeToGroup(gname, aname, &(edge_pairs[i][0]), 2);

      // create an iterator over the recursions
      grid::BoxIterator<2, IndexType> R = 
          grid::BoxIterator<2, IndexType> (Edges[i].get_rec_box());

      // create groups for the old, current, and new data
      gname2 = gname + "/new";
      hdf5_writer.addGroup(gname2);
      gname2 = gname + "/old";
      hdf5_writer.addGroup(gname2);
      gname2 = gname + "/current";
      hdf5_writer.addGroup(gname2);

      // Create a data set for each recursion level
      for (R.begin(); R.notDone(); ++R) {
        dsname = gname + "/old/Recursion_Level_" + to_string(R()[0]) + "_" + to_string(R()[1]); 
        hdf5_writer.addDataSet(dsname, \
                                Edges[i](CData2d::OLD, R()).get_data(), 
                                Edges[i](CData2d::OLD, R()).get_n());
       
        dsname = gname + "/current/Recursion_Level_" + to_string(R()[0]) + "_" + to_string(R()[1]); 
        hdf5_writer.addDataSet(dsname, 
                                Edges[i](CData2d::CUR, R()).get_data(), 
                                Edges[i](CData2d::CUR, R()).get_n());

        dsname = gname + "/new/Recursion_Level_" + to_string(R()[0]) + "_" + to_string(R()[1]);
        hdf5_writer.addDataSet(dsname, 
                                Edges[i](CData2d::NEW, R()).get_data(), 
                                Edges[i](CData2d::NEW, R()).get_n());
      }
    }

    // add a group for each of the corners we have data for and save the data
    // in them ...
    for (int i=0; i<num_corners; i++) {
      gname = "/Corner_" + to_string(i);
      hdf5_writer.addGroup(gname);

      // add the data and recursion extends as attributes
      aname = "data_extent_low";
      hdf5_writer.addAttributeToGroup(gname, aname, Corners[i].get_data_box().get_low().data(), DIM);
      aname = "data_extent_high";
      hdf5_writer.addAttributeToGroup(gname, aname, Corners[i].get_data_box().get_high().data(), DIM);
      aname = "data_ordering";
      hdf5_writer.addAttributeToGroup(gname, aname, Corners[i].get_data_box().get_ordering().data(), DIM);
      aname = "recursion_extent_low";
      hdf5_writer.addAttributeToGroup(gname, aname, Corners[i].get_rec_box().get_low().data(), 3);
      aname = "recursion_extent_high";
      hdf5_writer.addAttributeToGroup(gname, aname, Corners[i].get_rec_box().get_high().data(), 3);
      aname = "recursion_ordering";
      hdf5_writer.addAttributeToGroup(gname, aname, Corners[i].get_rec_box().get_ordering().data(), 3);

      // add the normal and boundary properties
      aname = "recursion_normals";
      hdf5_writer.addAttributeToGroup(gname, aname, Corners[i].get_sides().data(), 3);
      for (int j=0; j<2*DIM; j++)
        tmp_bounds.set_value(i, static_cast<int>(Corners[i].get_boundaries()[i]));
      aname = "boundaries";
      hdf5_writer.addAttributeToGroup(gname, aname, tmp_bounds.data(), 2*DIM);

      // add corner triples
      aname = "corner_triples";
      hdf5_writer.addAttributeToGroup(gname, aname, &(corner_triples[i][0]), 3);

      // create an iterator over the recursions
      grid::BoxIterator<3, IndexType> R = 
          grid::BoxIterator<3, IndexType> (Corners[i].get_rec_box());

      // create groups for the old, current, and new data
      gname2 = gname + "/new";
      hdf5_writer.addGroup(gname2);
      gname2 = gname + "/old";
      hdf5_writer.addGroup(gname2);
      gname2 = gname + "/current";
      hdf5_writer.addGroup(gname2);

      // Create a data set for each recursion level
      for (R.begin(); R.notDone(); ++R) {
        dsname = gname + "/old/Recursion_Level_" + to_string(R()[0]) 
                       + "_" + to_string(R()[1]) + "_" + to_string(R()[2]); 
        hdf5_writer.addDataSet(dsname, \
                                Corners[i](CData3d::OLD, R()).get_data(), 
                                Corners[i](CData3d::OLD, R()).get_n());
       
        dsname = gname + "/current/Recursion_Level_" + to_string(R()[0]) 
                       + "_" + to_string(R()[1]) + "_" + to_string(R()[2]); 
        hdf5_writer.addDataSet(dsname, 
                                Corners[i](CData3d::CUR, R()).get_data(), 
                                Corners[i](CData3d::CUR, R()).get_n());

        dsname = gname + "/new/Recursion_Level_" + to_string(R()[0]) 
                       + "_" + to_string(R()[1]) + "_" + to_string(R()[2]); 
        hdf5_writer.addDataSet(dsname, 
                                Corners[i](CData3d::NEW, R()).get_data(), 
                                Corners[i](CData3d::NEW, R()).get_n());
      }
    }

  } else {
    std::cerr << "CrbcUpdates::save_state: Data has not been initialized." << std::endl;
  }

}
#endif


//------------------------------------------------------------------------------
//                           Routines to set/get data
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::load_face_data(const int &side,
                        const IndexType P[DIM],
                        const DataType &val)
{
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData;
  have_new_input[side] = true;
  Faces[side](CData::NEW, zero)(P) = val;
}

template <int DIM, class DataType, class IndexType, class CoefType>
DataType CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_new_face_vals(const int &side,
                        const IndexType P[DIM]) const
{
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData;

  return Faces[side](CData::CUR, zero)(P);

}

template <int DIM, class DataType, class IndexType, class CoefType>
DataType* CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_auxiliary_vars 
                                         (const int &side, const IndexType &p)
{
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData;
  DataType *temp = NULL;
  grid::Point<1, IndexType> P(&p);

  #ifdef CRBC_DEBUG
  if (have_face[side]) {
    if ((p >= 0) && (p<= num_recursions[side])) {
      temp = (Faces[side](CData::NEW, P)).get_data();
    } else {
      std::cerr << "error get_auxiliary_vars: recursion index p out of bounds " << std::endl;
    }
  } else {
    std::cerr << "error get_auxiliary_vars: no recursions on this face" << std::endl;
  }
  #else
  temp = (Faces[side](CData::NEW, P)).get_data();
  #endif

  return temp;
}

template <int DIM, class DataType, class IndexType, class CoefType>
const DataType* CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_auxiliary_vars 
                                    (const int &side, const IndexType &p) const
{
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData;
  DataType *temp = NULL;
  grid::Point<1, IndexType> P(&p);

  #ifdef CRBC_DEBUG
  if (have_face[side]) {
    if ((p >= 0) && (p<= num_recursions[side])) {
      temp = (Faces[side](CData::NEW, P)).get_data();
    } else {
      std::cerr << "error get_auxiliary_vars: recursion index p out of bounds " << std::endl;
    }
  } else {
    std::cerr << "error get_auxiliary_vars: no recursions on this face" << std::endl;
  }
  #else
  temp = (Faces[side](CData::NEW, P)).get_data();
  #endif

  return temp;
}

template <int DIM, class DataType, class IndexType, class CoefType>
std::vector<int> CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_auxiliary_order
                                                         (const int &side) const
{

  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData;
  std::vector<int> temp;

  grid::Point<DIM, IndexType> order = (Faces[side](CData::NEW, 0)).get_data_box().get_ordering();

  for (int i=0; i<DIM; ++i)
    temp.push_back(static_cast<int>(order[i]));

  return temp;
}

template <int DIM, class DataType, class IndexType, class CoefType>
DataType CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_auxiliary_vars 
                                       (const int &side, 
                                        const IndexType P[DIM], 
                                        const IndexType &q) const
{
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData;
  grid::Point<1, IndexType> Q(&q);

  return Faces[side](CData::CUR, Q)(P);
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::set_auxiliary_vars 
                            (const int &side, 
                             const IndexType P[DIM], 
                             const IndexType &q,
                             const DataType &val)
{
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData;
  grid::Point<1, IndexType> Q(&q);

  Faces[side](CData::CUR, Q)(P) = val;
}

template <int DIM, class DataType, class IndexType, class CoefType>
int CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_edge_index
                                                     (const int &sidea, 
                                                      const int &sideb) const
{
  int i, temp = -1;

  for (i=0; i<num_edges; ++i) {
    if (((edge_pairs[i][0] == sidea) && (edge_pairs[i][1] == sideb)) ||
        ((edge_pairs[i][1] == sidea) && (edge_pairs[i][0] == sideb))) {
      temp = i;
      break;
    }
  }
  return temp;
}

template <int DIM, class DataType, class IndexType, class CoefType>
inline DataType CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_edge_auxiliary_vars 
                                       (const int &edge, 
                                        const IndexType P[DIM], 
                                        const IndexType q[2]) const
{
  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData;
  grid::Point<2, IndexType> Q(q);

  return Edges[edge](CData::CUR, Q)(P);
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::set_edge_auxiliary_vars 
                            (const int &edge, 
                             const IndexType P[DIM], 
                             const IndexType q[2],
                             const DataType &val)
{
  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData;
  grid::Point<2, IndexType> Q(q);

  Edges[edge](CData::CUR, Q)(P) = val;
}


//------------------------------------------------------------------------------
//                     initialization for edges/corners
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::initialize()
{
  int i, j, ind;

  // set up the edges
  for (i=0; i<num_edges; i++) {

    // figure out the recursions indexing 
    grid::Point<2, int> normals;
    grid::Point<2, IndexType> high, low;

    for (j=0; j<2; j++) {
      normals.set_value(j, Faces[edge_pairs[i][j]].get_sides()[0]);
      low.set_value(j, Faces[edge_pairs[i][j]].get_rec_box().get_low()[0]);
      high.set_value(j, Faces[edge_pairs[i][j]].get_rec_box().get_high()[0]);
    }

    grid::Box<2, IndexType> rec_box(low, high); 

    // figure out the data indexing
    grid::Box<DIM, IndexType> data_box = Faces[edge_pairs[i][0]].get_data_box()
                           .intersection(Faces[edge_pairs[i][1]].get_data_box());

    // expand in the recursion directions by 1
    for (j=0; j<2; j++) {
      if (Faces[edge_pairs[i][j]].get_sides()[0] % 2 == 0) { // left side
        data_box = data_box.add_low(grid::Point<DIM, IndexType>::scaled_unit
                              (Faces[edge_pairs[i][j]].get_sides()[0] / 2, -1));
      } else { // right side
        data_box = data_box.add_high(grid::Point<DIM, IndexType>::unit
                                  (Faces[edge_pairs[i][j]].get_sides()[0] / 2));
      }
    }

    // figure out "new" boundaries --- set all the sides "opposite" the CRBC
    // boundaries on this edge to be of type none ...
    grid::Point<DIM*2, Boundary> bounds = boundaries;
    for (j=0; j<2; j++) {
      
      // get the normal
      ind = Faces[edge_pairs[i][j]].get_sides()[0];
      
      // calculate the normal in the opposite direction
      ind = 2*(ind / 2) + ((1 + (ind % 2)) % 2);      

      bounds.set_value(ind, NONE);
    }

    // initialize the edge
    Edges[i] = CRBC_Data<2, DIM, DataType, IndexType, CoefType>(data_box, rec_box, 
                                                      normals, bounds);

    // now load the cosines into the data structure
    for (j=0; j<2; j++)
      Edges[i].calculate_coefficients(normals[j], recurs_coefs_a[normals[j]].data(), 
                                     recurs_coefs_ab[normals[j]].data(),
                                     recurs_coefs_sig[normals[j]].data(), 
                                     recurs_coefs_sigb[normals[j]].data(),
                                     c*dt/h[normals[j]/2]);

  }

  // set up the corners
  for (i=0; i<num_corners; i++) {

    grid::Box<DIM, IndexType> data_box, xdata_box, ydata_box, zdata_box;

    // figure out the recursions indexing 
    grid::Point<3, int> normals;
    grid::Point<3, IndexType> high, low;

    // figure out the normal directions for the recursions. This is done with 
    // assumming that we are using 3 dimensions, so we know we'll have a normal
    // in the x, y, and z directions. Therefore, we just need to figure out the 
    // sign, but we will also put them in order.
    for (j=0; j<3; j++) {
      // find the x-normal and recursion ranges:
      if (edge_pairs[corner_triples[i][j]][0] / 2 == 0) {

        normals.set_value(0, edge_pairs[corner_triples[i][j]][0]);
        low.set_value(0, Faces[edge_pairs[corner_triples[i][j]][0]].get_rec_box().get_low()[0]);
        high.set_value(0, Faces[edge_pairs[corner_triples[i][j]][0]].get_rec_box().get_high()[0]);
        
      } else if (edge_pairs[corner_triples[i][j]][1] / 2 == 0) {

        normals.set_value(0, edge_pairs[corner_triples[i][j]][1]);
        low.set_value(0, Faces[edge_pairs[corner_triples[i][j]][1]].get_rec_box().get_low()[0]);
        high.set_value(0, Faces[edge_pairs[corner_triples[i][j]][1]].get_rec_box().get_high()[0]);
        
      } 
    }
    for (j=0; j<3; j++) {
      // find the y-normal and recursion ranges:
      if (edge_pairs[corner_triples[i][j]][0] / 2 == 1) {

        normals.set_value(1, edge_pairs[corner_triples[i][j]][0]);
        low.set_value(1, Faces[edge_pairs[corner_triples[i][j]][0]].get_rec_box().get_low()[0]);
        high.set_value(1, Faces[edge_pairs[corner_triples[i][j]][0]].get_rec_box().get_high()[0]);
        
      } else if (edge_pairs[corner_triples[i][j]][1] / 2 == 1) {

        normals.set_value(1, edge_pairs[corner_triples[i][j]][1]);
        low.set_value(1, Faces[edge_pairs[corner_triples[i][j]][1]].get_rec_box().get_low()[0]);
        high.set_value(1, Faces[edge_pairs[corner_triples[i][j]][1]].get_rec_box().get_high()[0]);
        
      } 
    }
    for (j=0; j<3; j++) {
      // find the z-normal and recursion ranges:
      if (edge_pairs[corner_triples[i][j]][0] / 2 == 2) {

        normals.set_value(2, edge_pairs[corner_triples[i][j]][0]);
        low.set_value(2, Faces[edge_pairs[corner_triples[i][j]][0]].get_rec_box().get_low()[0]);
        high.set_value(2, Faces[edge_pairs[corner_triples[i][j]][0]].get_rec_box().get_high()[0]);
        
      } else if (edge_pairs[corner_triples[i][j]][1] / 2 == 2) {

        normals.set_value(2, edge_pairs[corner_triples[i][j]][1]);
        low.set_value(2, Faces[edge_pairs[corner_triples[i][j]][1]].get_rec_box().get_low()[0]);
        high.set_value(2, Faces[edge_pairs[corner_triples[i][j]][1]].get_rec_box().get_high()[0]);
        
      }
    }

    // set up the recursion indexing box
    grid::Box<3, IndexType> rec_box(low, high); 

    // caluclate the data box
    for (j=0; j<3; j++) {

      ind = 3 - edge_pairs[corner_triples[i][j]][0]/2 - edge_pairs[corner_triples[i][j]][1]/2;
      if (ind == 0) { // x-direction
        xdata_box = Edges[corner_triples[i][j]].get_data_box();
      } else if (ind == 1) { // y-direction
        ydata_box = Edges[corner_triples[i][j]].get_data_box();
      } else if (ind == 2) { // z-direction
        zdata_box = Edges[corner_triples[i][j]].get_data_box();
      }
    }
    data_box = xdata_box.intersection(ydata_box);
    data_box = data_box.intersection(zdata_box);

    // expand in the recursion directions by 1
    for (j=0; j<3; j++) {
      if (normals[j] % 2 == 0) { // left side
        data_box = data_box.add_low(grid::Point<DIM, IndexType>::scaled_unit(j, -1));
      } else { // right side
        data_box = data_box.add_high(grid::Point<DIM, IndexType>::unit(j));
      }
    }

    // figure out "new" boundaries --- set all the sides "opposite" the CRBC
    // boundaries on this edge to be of type none ...
    grid::Point<DIM*2, Boundary> bounds;
    for (j=0; j<2; j++) {
      bounds.set_value(ind, NONE);
    }
    // initialize the corner
    Corners[i] = CRBC_Data<3, DIM, DataType, IndexType, CoefType>(data_box, rec_box, 
                                                      normals, bounds);
    
    // now load the cosines into the data structure
    for (j=0; j<3; j++)
      Corners[i].calculate_coefficients(normals[j], recurs_coefs_a[normals[j]].data(), 
                                     recurs_coefs_ab[normals[j]].data(),
                                     recurs_coefs_sig[normals[j]].data(), 
                                     recurs_coefs_sigb[normals[j]].data(),
                                     c*dt/h[normals[j]/2]);

  }

  initialized = true;

}

//------------------------------------------------------------------------------
//                          Routine to run the updates
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::compute_updates()
{

  int i;

  // make sure the edges and corners have been set up
  if (!initialized)
    initialize();

  // calculate all of the wave equation updates
  for (i=0; i<6; i++) {
    if (have_face[i]) {

      //TODO better warnings / failure checks
      if (!have_new_input[i])
        std::cerr << "New data has not been input, updates will likely produce unexpected results" << std::endl;
 
      wave_equation_updates<1>(Faces[i]);
    }
  }

  // apply the recursion updates to the faces
  step_faces();

  // apply the recursion updates the edges
  step_edges();

  // apply the recursion updates to the corners
  step_corners();

  // copy the updates from the corners into the edges
  copy_corner_to_edge();

  // copy the updates from the edges into the faces
  copy_edge_to_face();

  // copy the current values into the old values and the new values into the
  // current values
  rotate_data();

  // set the flags that check for new input to false
  for (i = 0; i<6; i++)
    have_new_input[i] = false;

}


//------------------------------------------------------------------------------
//                              return functions
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
double CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_max_reflection_coef() const
{

  double max;
  max = ref_coefs[0];
  for (int i=1; i<6; i++)
    max = (max < ref_coefs[i]) ? ref_coefs[i] : max;

  return max;

}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_input_extents
              (const int &side, IndexType low[DIM], IndexType high[DIM]) const
{
  int i;

  if (have_face[side]) {

    // get the interior side of the data extents
    grid::Box<DIM,IndexType> data_box = Faces[side].get_data_box().get_side(side/2, (side % 2 == 1));

    for (i=0; i<DIM; i++) {
      low[i] = data_box.get_low()[i];
      high[i] = data_box.get_high()[i];
    }

  } else { 

    // return indexing that won't go over any values
    for (i=0; i<DIM; i++) {
      low[i] = 0;
      high[i] = -1;
    }
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_output_extents
              (const int &side, IndexType low[DIM], IndexType high[DIM]) const
{
  int i;

  if (have_face[side]) {

    // get the data extents
    grid::Box<DIM,IndexType> data_box = Faces[side].get_data_box();

    // shrink the data box by one in the outward and inward normal directions
    // of the side (1 of these corresponds to input, the other to a "ghost"
    // recursion point that we don't update on the edges and corners
    data_box = data_box.add_high(grid::Point<DIM, IndexType>::scaled_unit(side/2, -1));
    data_box = data_box.add_low(grid::Point<DIM, IndexType>::unit(side/2));

    for (i=0; i<DIM; i++) {
      low[i] = data_box.get_low()[i];
      high[i] = data_box.get_high()[i];
    }

  } else { 

    // return indexing that won't go over any values
    for (i=0; i<DIM; i++) {
      low[i] = 0;
      high[i] = -1;
    }

  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_edge_extents
                         (const int &edge, 
                          IndexType low[DIM], 
                          IndexType high[DIM],
                          IndexType plow[2],
                          IndexType phigh[2]) const
{
  int i;

  if ((edge > -1) && (edge < num_edges)) {

    // get the data extents
    grid::Box<DIM,IndexType> data_box = Edges[edge].get_data_box();

    // shrink the data box by one in the outward and inward normal directions
    // of the side (1 of these corresponds to input, the other to a "ghost"
    // recursion point that we don't update on the edges and corners
    for (i=0; i<2; i++) {
      data_box = data_box.add_high(grid::Point<DIM, IndexType>::scaled_unit(edge_pairs[edge][i]/2, -1));
      data_box = data_box.add_low(grid::Point<DIM, IndexType>::unit(edge_pairs[edge][i]/2));
    }

    for (i=0; i<DIM; i++) {
      low[i] = data_box.get_low()[i];
      high[i] = data_box.get_high()[i];
    }

    // get the recursion extents
    grid::Box<2,IndexType> rec_box = Edges[edge].get_rec_box();

    for (i=0; i<2; i++) {
      plow[i] = rec_box.get_low()[i];
      phigh[i] = rec_box.get_high()[i];
    }

  } else { 

    // return indexing that won't go over any values
    for (i=0; i<DIM; i++) {
      low[i] = 0;
      high[i] = -1;
    }

    for (i=0; i<2; i++) {
      plow[i] = 0;
      phigh[i] = -1;
    }

  }
}

//------------------------------------------------------------------------------
//                                 Wave Equation
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
template <int REC_DIM>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::wave_equation_updates
               (CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data)
{

  // create an iterator over the recursion indices.
  grid::BoxIterator<REC_DIM, IndexType> i = 
      grid::BoxIterator<REC_DIM, IndexType> (crbc_data.get_rec_box());

  // loop over the recursion indices
  for (i.begin(); i.notDone(); ++i) {

    this->get_wave_eq_update(crbc_data.get_boundaries(),  \
                                    crbc_data(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType>::OLD, i()), \
                                    crbc_data(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType>::CUR, i()), \
                                    crbc_data(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType>::NEW, i()));

  }
}

//------------------------------------------------------------------------------
//                                 Recursions
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
template <int REC_DIM>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::recursion_updates
                (CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data)
{
  int i;

  for (i=0; i<REC_DIM; i++) {
    this -> template forward_recursions<REC_DIM> (crbc_data, \
                                                  crbc_data.get_sides()[i]/2);
    this -> template terminations<REC_DIM> (crbc_data, \
                                                  crbc_data.get_sides()[i]/2);
    this -> template backward_recursions<REC_DIM> (crbc_data, \
                                                  crbc_data.get_sides()[i]/2);
  }
  
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::step_faces() 
{
  int i;
  // loop over the possible faces
  for (i=0; i<6; i++)
    if (have_face[i])
      recursion_updates<1>(Faces[i]);
}

template <int DIM, class DataType, class IndexType, class CoefType>    
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::step_edges() 
{
  int i;

  grid::Point<2, IndexType> zero;
  grid::BoxIterator<DIM, IndexType> data(Edges[3].get_data_box());

  // copy in the data from the faces
  copy_face_to_edge();

  // wave equation updates
  for (i=0; i<num_edges; i++)
    wave_equation_updates<2>(Edges[i]);

  // loop over the edges
  for (i=0; i<num_edges; i++) 
    recursion_updates<2>(Edges[i]);

}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::step_corners() 
{
  int i;

  // copy in the data from the edges
  copy_edge_to_corner();

  // wave equation updates
  for (i=0; i<num_corners; i++)
    wave_equation_updates<3>(Corners[i]);

  // loop over the corners
  for (i=0; i<num_corners; i++) 
    recursion_updates<3>(Corners[i]);
}

//------------------------------------------------------------------------------
//                          Data Copying Routines
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::copy_face_to_edge()
{

  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData2d;
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData1d;
  int i, j, k;

  grid::Box<DIM, IndexType> data_box;
  grid::BoxIterator<DIM, IndexType> diter;
  grid::BoxIterator<1, IndexType> riter;

  // loop over all of the edges
  for (i=0; i<num_edges; i++) {

    // loop over the faces on the current edge
    for (j=0; j<2; j++) {

      // FIXME modify box indices so they don't overwrite precomputed WE updates
      // This isn't a huge issue, but it's a little bit nicer to do all of the
      // wave equation updates together and all of the recursion updates together
      // instead of alternating based on dimension
      // figure out the data ranges
      // data_box = (Edges[i].get_data_box()).intersection(Faces[edge_pairs[i][j]].get_data_box());

      // figure out the other side ...
      k = (j%2 + 1) % 2;
      data_box = Edges[i].get_data_box().get_side(edge_pairs[i][k]/2, (edge_pairs[i][k]%2 == 1));

      // create an iterator over the data we want to copy
      diter = grid::BoxIterator<DIM, IndexType>(data_box);

      // create an iterator over the face recursions
      riter = grid::BoxIterator<1, IndexType>(Faces[edge_pairs[i][j]].get_rec_box());

      // loop over the recurions
      for (riter.begin(); riter.notDone(); ++riter) {
        // loop over the data
        for (diter.begin(); diter.notDone(); ++diter) {
          
          // copy data from the face into the edge
          Edges[i](CData2d::NEW, grid::Point<2, IndexType>::scaled_unit(j, riter()[0]))(diter()) = \
            Faces[edge_pairs[i][j]](CData1d::NEW, riter())(diter());

        }
      }
    }
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::copy_edge_to_corner()
{

  typedef CRBC_Data<3, DIM, DataType, IndexType, CoefType> CData3d;
  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData2d;
  typedef grid::Point<3, IndexType> P;
  int i, j, k, ind1, ind2;

  grid::Box<DIM, IndexType> data_box;
  grid::BoxIterator<DIM, IndexType> diter;
  grid::BoxIterator<2, IndexType> riter;

  // loop over all of the corners
  for (i=0; i<num_corners; i++) {

    // loop over the edges on the current corner
    for (j=0; j<3; j++) {

      // FIXME modify box indices so they don't overwrite precomputed WE updates
      // This isn't a huge issue, but it's a little bit nicer to do all of the
      // wave equation updates together and all of the recursion updates together
      // instead of alternatiing based on dimension
      // figure out the data ranges
      data_box = Corners[i].get_data_box().intersection(Edges[corner_triples[i][j]].get_data_box());

      // shrink to region that was updated on the edges
      for (k=0; k<3; k++) {
        
        if ((Corners[i].get_sides()[k] != Edges[corner_triples[i][j]].get_sides()[0]) &&
            (Corners[i].get_sides()[k] != Edges[corner_triples[i][j]].get_sides()[1])) {

          if (Corners[i].get_sides()[k] % 2 == 0) { // left side
            data_box = data_box.add_low(grid::Point<DIM, IndexType>::unit(k));
          } else { 
            data_box = data_box.add_high(grid::Point<DIM, IndexType>::scaled_unit(k, -1));
          }

          break;
        }
      }

      // create an iterator over the data we want to copy
      diter = grid::BoxIterator<DIM, IndexType>(data_box);

      // create an iterator over the edge recursions
      riter = grid::BoxIterator<2, IndexType>(Edges[corner_triples[i][j]].get_rec_box());

      // figure out which components the edge recursions are in
      ind1 = Edges[corner_triples[i][j]].get_sides()[0]/2;
      ind2 = Edges[corner_triples[i][j]].get_sides()[1]/2;

      // loop over the recurions
      for (riter.begin(); riter.notDone(); ++riter) {
        // loop over the data
        for (diter.begin(); diter.notDone(); ++diter) {
          
          // copy data from the edge into the corner
          Corners[i](CData3d::NEW, P::scaled_unit(ind1, riter()[0]) + \
                      P::scaled_unit(ind2, riter()[1]))(diter()) = \
            Edges[corner_triples[i][j]](CData2d::NEW, riter())(diter());

        }
      }
    }
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::copy_edge_to_face()
{

  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData2d;
  typedef CRBC_Data<1, DIM, DataType, IndexType, CoefType> CData1d;
  int i, j, k;

  grid::Box<DIM, IndexType> data_box;
  grid::BoxIterator<DIM, IndexType> diter;
  grid::BoxIterator<1, IndexType> riter;

  // loop over all of the edges
  for (i=0; i<num_edges; i++) {

    // loop over the faces on the current edge
    for (j=0; j<2; j++) {

      // figure out the data ranges (only need to copy data in "interior")
      data_box = Edges[i].get_data_box().grow(-1).intersection(Faces[edge_pairs[i][j]].get_data_box());

      if (DIM == 3) {
        k = 3 - Edges[i].get_sides()[0]/2 - Edges[i].get_sides()[1]/2;
        data_box = data_box.add_high(grid::Point<DIM, IndexType>::unit(k));
        data_box = data_box.add_low(grid::Point<DIM, IndexType>::scaled_unit(k, -1));
      }

      // create an iterator over the data we want to copy
      diter = grid::BoxIterator<DIM, IndexType>(data_box);

      // create an iterator over the face recursions
      riter = grid::BoxIterator<1, IndexType>(Faces[edge_pairs[i][j]].get_rec_box());

      // loop over the recurions
      for (riter.begin(); riter.notDone(); ++riter) {
        // loop over the data
        for (diter.begin(); diter.notDone(); ++diter) {
          
          // copy data from the edge into the face
          Faces[edge_pairs[i][j]](CData1d::NEW, riter())(diter()) = \
            Edges[i](CData2d::NEW, grid::Point<2, IndexType>::scaled_unit(j, riter()[0]))(diter());

        }
      }
    }
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::copy_corner_to_edge()
{

  typedef CRBC_Data<3, DIM, DataType, IndexType, CoefType> CData3d;
  typedef CRBC_Data<2, DIM, DataType, IndexType, CoefType> CData2d;
  typedef grid::Point<3, IndexType> P;
  int i, j, ind1, ind2;

  grid::Box<DIM, IndexType> data_box;
  grid::BoxIterator<DIM, IndexType> diter;
  grid::BoxIterator<2, IndexType> riter;

  // loop over all of the corners
  for (i=0; i<num_corners; i++) {

    // loop over the edges on the current corner
    for (j=0; j<3; j++) {

      // figure out the data ranges (only need to copy data in "interior")
      data_box = Corners[i].get_data_box().grow(-1).intersection(Edges[corner_triples[i][j]].get_data_box());

      // create an iterator over the data we want to copy
      diter = grid::BoxIterator<DIM, IndexType>(data_box);

      // create an iterator over the edge recursions
      riter = grid::BoxIterator<2, IndexType>(Edges[corner_triples[i][j]].get_rec_box());

      // figure out which components the edge recursions are in
      ind1 = Edges[corner_triples[i][j]].get_sides()[0]/2;
      ind2 = Edges[corner_triples[i][j]].get_sides()[1]/2;

      // loop over the recurions
      for (riter.begin(); riter.notDone(); ++riter) {
        // loop over the data
        for (diter.begin(); diter.notDone(); ++diter) {

          // copy data from the corner into the edge
          Edges[corner_triples[i][j]](CData2d::NEW, riter())(diter()) = \
            Corners[i](CData3d::NEW, P::scaled_unit(ind1, riter()[0]) + \
                      P::scaled_unit(ind2, riter()[1]))(diter());

        }
      }
    }
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::rotate_data()
{
  int i;

  // rotate the face data
  for (i=0; i<6; i++)
    if (have_face[i])
      Faces[i].rotate_data();

  // rotate the edge data
  for (i=0; i<num_edges; i++)
    Edges[i].rotate_data();

  // rotate the corner data
  for (i=0; i<num_corners; i++)
    Corners[i].rotate_data();

}

//------------------------------------------------------------------------------
//                              Helper functions
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::set_defaults()
{

  int i;

  num_faces = 0;
  num_edges = 0;
  num_corners = 0;

  for (i = 0; i<6; i++) {
    have_face[i] = false;
    have_new_input[i] = false;
    num_recursions[i] = 5;
    ref_coefs[i] = -1.0;
    boundaries.set_value(i, NONE);
    have_cosines[i] = false;
    delta[i] = -1.0;
  }

  Pmax = 20;
  tol = 1e-2;

  useTol = false;
  initialized = false;

}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::find_edges()
{
  int i, j;

  // loop over sides
  for (i=0; i<2*DIM-1; i++) {
    // get a second side to check
    for (j=i+1; j<2*DIM; j++) {
      // make sure the sides are not parallel
      if (j/2 != i/2) {
        
        if ((boundaries[i] == CRBC) && (boundaries[j] == CRBC)) {
          edge_pairs[num_edges][0] = i;
          edge_pairs[num_edges][1] = j;
          num_edges++;
           
        }
      }
    }
  }
}
 
template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::find_corners()
{
  int i, j, k, x, y, z;
  int edge_dir[3];
  bool test;

  // just to avoid a "possibly unitialized warning" that shouldn't be able to 
  // occur ... this should force an out of bounds index if the warning is in 
  // fact true ...
  x = y = z = -1;

  // FIXME I believe this could be much more straightforward. In anycase, we 
  // only do this once and the expectation is that the updates will be orders
  // of magnitude more expensive, so there isn't really any reason to change it
  // if it isn't broken ...

  // loop over edges
  for (i=0; i<num_edges-2; i++) {
    // get a second edge to check
    for (j=i+1; j<num_edges-1; j++) {
      // get a third edge to check
      for (k=j+1; k<num_edges; k++) {
        // make sure all 3 edges are in different directions
        // calculate the directions along each edge
        edge_dir[0] = 3 - edge_pairs[i][0]/2 - edge_pairs[i][1]/2;
        edge_dir[1] = 3 - edge_pairs[j][0]/2 - edge_pairs[j][1]/2;
        edge_dir[2] = 3 - edge_pairs[k][0]/2 - edge_pairs[k][1]/2;

        if ((edge_dir[0] != edge_dir[1]) && (edge_dir[0] != edge_dir[2]) && \
                                            (edge_dir[1] != edge_dir[2])) {

          test = true;
          // figure out which index corresponds to each direction
          switch (edge_dir[0]) {
            case 0 : x = i; break;
            case 1 : y = i; break;
            case 2 : z = i; break;
          }
          switch (edge_dir[1]) {
            case 0 : x = j; break;
            case 1 : y = j; break;
            case 2 : z = j; break;
          }
          switch (edge_dir[2]) {
            case 0 : x = k; break;
            case 1 : y = k; break;
            case 2 : z = k; break;
          }

          // the edges in the x & y directions need to share a face in z
          if ((edge_pairs[x][0] != edge_pairs[y][0]) && \
              (edge_pairs[x][0] != edge_pairs[y][1]) && \
              (edge_pairs[x][1] != edge_pairs[y][0]) && \
              (edge_pairs[x][1] != edge_pairs[y][1])) {
            test = false;
          }

          // the edges in the x & z directions need to share a face in y
          if ((edge_pairs[x][0] != edge_pairs[z][0]) && \
              (edge_pairs[x][0] != edge_pairs[z][1]) && \
              (edge_pairs[x][1] != edge_pairs[z][0]) && \
              (edge_pairs[x][1] != edge_pairs[z][1])) {
            test = false;
          }

          // the edges in the y & z directions need to share a face in x
          if ((edge_pairs[y][0] != edge_pairs[z][0]) && \
              (edge_pairs[y][0] != edge_pairs[z][1]) && \
              (edge_pairs[y][1] != edge_pairs[z][0]) && \
              (edge_pairs[y][1] != edge_pairs[z][1])) {
            test = false;
          }

          if (test) {
            corner_triples[num_corners][0] = i;
            corner_triples[num_corners][1] = j;
            corner_triples[num_corners][2] = k;
            num_corners++;
          }
        }
      }
    }
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::check_cos_params(double &eta, int &P)
{

  if (eta > 0.1) {
    #ifdef CRBC_DEBUG
    std::cerr << "warning: eta > 0.1, using eta = 0.1 to compute optimal cosines "
              << "(this is a harder problem than requested)" << std::endl;
    #endif
    eta = 0.1;
  } else if (eta < 1e-7) {
    #ifdef CRBC_DEBUG
    std::cerr << "warning: eta < 1e-7, using eta = 1e-7 to compute optimal cosines "
              << "(this is an easier problem than requested)" << std::endl;
    #endif
    eta = 1e-7;
  }

  if (P > 40) {
    #ifdef CRBC_DEBUG
    std::cerr << "warning: P (or Pmax) > 40, trying P=40" << std::endl;
    #endif
    P = 40;
  } else if (P < 0) {
    #ifdef CRBC_DEBUG
    std::cerr << "warning: P (or Pmax) < 0, trying P=0 (Sommerfeld)" << std::endl;
    #endif
    P = 0;
  }

}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::compute_crbc_coefficients
                                       (const int &side, const double *cosines)
{

  int j;

  // allocate memory for coefficients
  recurs_coefs_a[side].reserve(num_recursions[side]);
  recurs_coefs_ab[side].reserve(num_recursions[side]);
  recurs_coefs_sig[side].reserve(num_recursions[side]);
  recurs_coefs_sigb[side].reserve(num_recursions[side]);

  for (j=0; j<num_recursions[side]; j++) {
    recurs_coefs_a[side].push_back (static_cast<CoefType> (cosines[2*j]));
    recurs_coefs_ab[side].push_back (static_cast<CoefType> (cosines[2*j+1]));
  }

  for (j=0; j<num_recursions[side]; j++) {
    recurs_coefs_sig[side].push_back (static_cast<CoefType> (\
         (0.5 * dt)*(1.0 - recurs_coefs_a[side][j]*recurs_coefs_a[side][j]) \
         / (T*recurs_coefs_a[side][j])));
    recurs_coefs_sigb[side].push_back (static_cast<CoefType> (\
         (0.5 * dt)*(1.0 - recurs_coefs_ab[side][j]*recurs_coefs_ab[side][j]) \
         / (T*recurs_coefs_ab[side][j])));
  }

  have_cosines[side] = true;
}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::init_face_params
                          (const int &side, 
                          const IndexType low_index[DIM], 
                          const IndexType high_index[DIM])
{

  grid::Point<DIM,IndexType> low_ind(low_index);
  grid::Point<DIM,IndexType> high_ind(high_index);
  // expand the domain box to add a ghost layer
  if (side % 2 == 0) { // left side
    low_ind.add_to_component(side/2, -1);   
  } else { // right side
    high_ind.add_to_component(side/2, 1);   
  }
  grid::Box<DIM, IndexType> domain = grid::Box<DIM, IndexType>(low_ind, high_ind);
  
  // set up a box over the recursions
  grid::Point<1, IndexType> low = grid::Point<1, IndexType>();
  grid::Point<1, IndexType> high = grid::Point<1, IndexType>::scaled_ones(num_recursions[side]);
  grid::Box<1, IndexType> recursions = grid::Box<1, IndexType>(low, high);
  
  grid::Point<1, int> normals = grid::Point<1, int>::scaled_ones(side);
  
  grid::Point<DIM*2, Boundary> bounds = boundaries;
  
  // set the opposite boundary to type NONE
  bounds.set_value(2*(side / 2) + ((1 + (side % 2)) % 2), NONE);
  
  Faces[side] = CRBC_Data<1, DIM, DataType, IndexType, CoefType>(domain, recursions, normals, bounds); 

  // now load the cosines into the data structure
  Faces[side].calculate_coefficients(side, recurs_coefs_a[side].data(), 
                                     recurs_coefs_ab[side].data(),
                                     recurs_coefs_sig[side].data(), 
                                     recurs_coefs_sigb[side].data(),
                                     c*dt/h[side/2]);

  have_new_input[side] = true;

}

template <int DIM, class DataType, class IndexType, class CoefType>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::print_opt_cos_info(const int &flag) const
{
  if (flag != 0) {
    std::cerr << "optimal_cosinesP failed with flag = " << flag;

    switch (flag) {
      case 1:
        std::cerr << " (eta out of range (eta < 1e-7 or eta > 0.1))" << std::endl;
        break;
      case 2:
        std::cerr << " (P out of range (P < 1 or P > 40))" << std::endl;
        break;
      case 4:
        std::cerr << " (Remez failed to converge)" << std::endl;
        break;
      case 5:
        std::cerr << " (Newton failed to converge)" << std::endl;
        break;
      case 6:
        std::cerr << " (Linear solve failed)" << std::endl;
        break;
    }
  }
}

 #if USE_HDF5_RESTARTS
template <int DIM, class DataType, class IndexType, class CoefType>
template <int REC_DIM>
void CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_crbc_props(grid::Point<DIM, IndexType> &domain_low,
                     grid::Point<DIM, IndexType> &domain_high,
                     grid::Point<DIM, IndexType> &domain_order, 
                     grid::Point<REC_DIM, IndexType> &rec_low, 
                     grid::Point<REC_DIM, IndexType> &rec_high,
                     grid::Point<REC_DIM, IndexType> &rec_order,
                     grid::Point<REC_DIM, int> &normals, 
                     grid::Point<DIM*2, Boundary> &bounds,
                     util_hdf5::ReadHDF5 &hdf5_reader,
                     const std::string &gname) const
{

  int l, temp2[2*DIM];
  std::string aname;

  aname = "data_extent_low";
  hdf5_reader.getAttributeFromGroup(gname, aname, domain_low.data());
  aname = "data_extent_high";
  hdf5_reader.getAttributeFromGroup(gname, aname, domain_high.data());
  aname = "data_ordering";
  hdf5_reader.getAttributeFromGroup(gname, aname, domain_order.data());
  aname = "recursion_extent_low";
  hdf5_reader.getAttributeFromGroup(gname, aname, rec_low.data());
  aname = "recursion_extent_high";
  hdf5_reader.getAttributeFromGroup(gname, aname, rec_high.data());
  aname = "recursion_ordering";
  hdf5_reader.getAttributeFromGroup(gname, aname, rec_order.data());
  aname = "recursion_normals";
  hdf5_reader.getAttributeFromGroup(gname, aname, normals.data());
  aname = "boundaries";
  hdf5_reader.getAttributeFromGroup(gname, aname, &(temp2[0]));
  for (l=0; l<2*DIM; l++)
    bounds.set_value(l, static_cast<Boundary>(temp2[l]));

}
#endif

template <int DIM, class DataType, class IndexType, class CoefType>
double CrbcUpdates<DIM, DataType, IndexType, CoefType>::get_mem_usage() const
{
  std::size_t i, mem_use = 0;

  for (i=0; i<6; ++i)
    if (have_face[i])
      mem_use += Faces[i].get_mem_usage();

  for (i=0; i<num_edges; ++i)
    mem_use += Edges[i].get_mem_usage();

  for (i=0; i<num_corners; ++i)
    mem_use += Corners[i].get_mem_usage();

  return (mem_use / ((double) 1024*1024));

}


} // end namespace
#endif // CRBC_UPDATES_H_
