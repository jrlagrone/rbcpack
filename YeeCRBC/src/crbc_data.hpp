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

/* This file defines a data class for holding the data needed to do the crbc/dab
   update operations
*/
#ifndef CRBC_DATA_H_
#define CRBC_DATA_H_

// include std::vector
#include <vector>

// include basic grid data objects
#include "grid/point.hpp"
#include "grid/box.hpp"
#include "grid/grid_array.hpp"
#include "boundary_properties.hpp"

#include <utility> // std::swap

namespace crbc {

// TODO: add generic components for non-maxwell wave equations
// and figure out if we can have a mix of E & H fields
/// Enumeration of the valid field component types.
enum Field {
  Fx1 = 0,     /**< x component of the field  */
  Fx2,         /**< y component of the field  */
  Fx3,         /**< z component of the field  */
};


/// Class to create a data structure for the CRBC/DAB auxiliary variables
/// based off of the Array data structure in the grid/grid_arry.hpp
/// 
/// the template parameters are: \n
/// \tparam REC_DIM   the number recursions directions needed. Typically, 1 
///                     corresponds to a face, 2 to an edge, 3 to a corner, etc.
/// \tparam DATA_DIM  the number of physical dimensions for the problem 
/// \tparam DataType  The type of data to use. This isn't checked, but it
///                     should be a float type. Note, we calculate the cosines
///                     using double precision.
/// \tparam IndexType The data type to use for indexing. This isn't checked, but
///                     should be an integral type. This needs to be large
///                     enough to store the number of points each dimension.
/// \tparam CoefType   The data type used for the coefficients such as time step
///                    size, wave speed, etc.
template <int REC_DIM,
          int DATA_DIM = 3, 
          class DataType = double, 
          class IndexType = long int,
          class CoefType = DataType>
class CRBC_Data : virtual public BoundaryProperties
{

  public:

    /// enumeration to make accessing the data at different times more readable
    enum Time {
      OLD = 0, //< indentifier for time t-1
      CUR,     //< indentifier for time t
      NEW      //< indentifier for time t+1
    };

    /// constructor --- default, mostly unitialized ....
    CRBC_Data();

    /// constructor --- creates data storage for auxiliary variables.
    /// The box domain provides the index constraints for the data 
    /// and the box rescursions provides the index constraints for
    /// the auxiliary variables. Additionally, we require the field component
    /// and outward normal direction(s) to be identified.
    ///
    /// \param[in] domain      A box defining the index bounds for the data
    /// \param[in] recursions  A box defining the bounds for the recursions
    ///                        indices
    /// \param[in] normals     The normal directions that the CRBC recursions 
    ///                        are defined on
    /// \param[in] bounds      The boundary conditions for the domain
    CRBC_Data(const grid::Box<DATA_DIM, IndexType> &domain, 
              const grid::Box<REC_DIM, IndexType> &recursions,
              const grid::Point<REC_DIM, int> &normals,
              const grid::Point<DATA_DIM*2, Boundary> &bounds);

    /// assignment operator
    CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType> & operator=
            (const CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType> &B);

    /// destructor
    virtual ~CRBC_Data() {};

    /// Routine to calculate the coefficients used in the recurion updates
    ///
    /// \param[in] normal  normal direction to calculate the parameters for
    /// \param[in] a       array containing a[j] = \f$ \frac{\cos \phi_j}{c} \f$
    /// \param[in] ab      array containing ab[j] = \f$ \frac{\cos \bar{\phi}_j}{c} \f$
    /// \param[in] sig     array containing sig[j] = 
    ///                    \f$ \frac{1}{cT} \frac{\sin^2 \phi_j}{\cos \phi_j} \f$
    /// \param[in] sigb    array containing sigb[j] = 
    ///                    \f$ \frac{1}{cT} \frac{\sin^2 \bar{\phi}_j}{\cos \bar{\phi}_j} \f$
    /// \param[in] ccfl    wave speed times the cfl condition --- c*dt/h, where
    ///                    h is the grid spacing in the given direction
    void calculate_coefficients(const int &normal,
                                const CoefType *a,
                                const CoefType *ab,
                                const CoefType *sig,
                                const CoefType *sigb,
                                const CoefType &ccfl);
    
    /// Routine to get the coefficients needed to calculate the forward 
    /// recursions
    ///
    /// \param[in] comp_ind  The index for the direction we want to get the 
    ///                       recursion coefficients in
    /// \param[in] p          The index of the auxilliary variable we want the 
    ///                       the coefficient for
    /// \param[in] i          The index of the term we want to coefficient for
    ///
    /// \return               The coefficient for the ith term of the forward
    ///                       recursion in the appropriate direction  
    inline CoefType get_fwd_coef(int &comp_ind, const IndexType &p, const int i) const;

    /// Routine to get the coefficients needed to calculate the backward 
    /// recursions
    ///
    /// \param[in] comp_ind   The index for the direction we want to get the 
    ///                       recursion coefficients in
    /// \param[in] p          The index of the auxilliary variable we want the 
    ///                       the coefficient for
    /// \param[in] i          The index of the term we want to coefficient for
    ///
    /// \return               The coefficient for the ith term of the forward
    ///                       recursion in the appropriate direction  
    inline CoefType get_bwd_coef(int &comp_ind, const IndexType &p, const int i) const;

    /// Routine to the wave speed times the cfl number in the direction requested
    ///
    /// \param[in] comp_ind   The index for the direction we want to get the 
    ///                       recursion coefficients in
    ///
    /// \return               c*dt/h[comp_ind]  
    inline CoefType get_ccfl(int &comp_ind) const;

    /// method to access the data array associated with the given recursion
    /// point index
    /// 
    /// \param[in] T         Time to get data reference at, either old, current,
    ///                      or new
    /// \param[in] rec_index The recursion index to get data at
    ///
    /// \return              A reference to a grid::Array object, which can be
    ///                      used to get or modify data
    inline grid::Array<DATA_DIM, DataType, IndexType> & operator()
                         (const Time &T,
                          const grid::Point<REC_DIM, IndexType> &rec_index);

    inline const grid::Array<DATA_DIM, DataType, IndexType> & operator()
                         (const Time &T,
                          const grid::Point<REC_DIM, IndexType> &rec_index) const;

    /// method to access the data array associated with the given recursion
    /// point index
    /// 
    /// \param[in] T         Time to get data reference at, either old, current,
    ///                      or new
    /// \param[in] rec_index The logical recursion index to get data at
    ///
    /// \return              A reference to a grid::Array object, which can be
    ///                      used to get or modify data
    inline grid::Array<DATA_DIM, DataType, IndexType> & operator()
                         (const Time &T,
                          const std::size_t &rec_index);

    inline const grid::Array<DATA_DIM, DataType, IndexType> & operator()
                         (const Time &T,
                          const std::size_t &rec_index) const;

    /// get normal direction(s)
    inline grid::Point<REC_DIM, int> get_sides() const {return sides;};

    /// get boundaries
    ///
    /// \return boundary conditions
    inline grid::Point<2*DATA_DIM, Boundary> get_boundaries() const {return bounds;};

    /// get the box for the recursions
    ///
    /// \return box defining indexing for the recursions
    inline grid::Box<REC_DIM, IndexType> get_rec_box() const {return rec_box;};

    /// get the box for the data
    ///
    /// \return box defining indexing for the underlying data
    inline grid::Box<DATA_DIM, IndexType> get_data_box() const {return data_box;};

    /// method to move current data -> old data and new data -> current data
    /// this is done by swapping pointers, the data is not actually moved
    void rotate_data();

    /// method to zero out the new data. We don't need to do this to actually
    /// complete the updates successfully, but it makes some debugging tasks
    /// easier (such as checking the the input data is correctly copied in)
    void zero_new_data();

    /// method to return the memory usage
    std::size_t get_mem_usage() const;
    

  private:

    bool dataAllocated;
    grid::Box<DATA_DIM, IndexType> data_box;
    grid::Box<REC_DIM, IndexType> rec_box;
    std::vector<grid::Array<DATA_DIM, DataType, IndexType> > data[3];
    grid::Point<REC_DIM, int> sides;
    grid::Point<DATA_DIM*2, Boundary> bounds;
    std::size_t n;
    DataType ccfl[REC_DIM];
    std::vector<CoefType> fwd_recursion_coefs[REC_DIM];
    std::vector<CoefType> bwd_recursion_coefs[REC_DIM];
    bool have_recursion_coefs[REC_DIM];

}; // end CRBC_Data class

//------------------------------------------------------------------------------
//                               definitions
//------------------------------------------------------------------------------
template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::CRBC_Data()
{

  for (int i=0; i<REC_DIM; i++)
    have_recursion_coefs[i] = false;
  dataAllocated = false;

}


template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::CRBC_Data
             (const grid::Box<DATA_DIM, IndexType> &domain, 
              const grid::Box<REC_DIM, IndexType> &recursions,
              const grid::Point<REC_DIM, int> &normal,
              const grid::Point<DATA_DIM*2, Boundary> &bounds)
{

  data_box = domain;
  rec_box = recursions;
  sides = normal;
  this->bounds = bounds;
  n = rec_box.get_npoints();
  data[OLD].assign(n, grid::Array<DATA_DIM, DataType, IndexType>());
  data[CUR].assign(n, grid::Array<DATA_DIM, DataType, IndexType>());
  data[NEW].assign(n, grid::Array<DATA_DIM, DataType, IndexType>());
  for (int i=0; i<REC_DIM; i++)
    have_recursion_coefs[i] = false;
  
  // figure out an ordering based on the normals. The idea here is that we are
  // calculating derivatives in the direction of the normals so we should 
  // make the strides shorter in these directions and hopefully gain some
  // runtime benefit
  grid::Point<DATA_DIM, IndexType> order = grid::Point<DATA_DIM, IndexType>();
  for (int i=0; i<DATA_DIM; i++)
    order.set_value(i,i);        // set the default order
  if (REC_DIM <= DATA_DIM)  
    for (int i=0; i<REC_DIM; i++) {  // reorder
      for (int j=0; j<DATA_DIM; j++)
        if (order[j] == sides[i]/2)
          order.set_value(j, order[i]);
      order.set_value(i, sides[i]/2);  
    }

  data_box.set_ordering(order);
  for (std::size_t i = 0; i<n; i++) {
    data[OLD].at(i) = grid::Array<DATA_DIM, DataType, IndexType>(data_box);
    data[CUR].at(i) = grid::Array<DATA_DIM, DataType, IndexType>(data_box);
    data[NEW].at(i) = grid::Array<DATA_DIM, DataType, IndexType>(data_box);
  }
  dataAllocated = true;

}
template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType> & 
            CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::operator=
            (const CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType> &B)
{

    dataAllocated = B.dataAllocated;
    data_box = B.data_box;
    rec_box = B.rec_box;
    data[0] = B.data[0];
    data[1] = B.data[1];
    data[2] = B.data[2];
    sides = B.sides;
    bounds = B.bounds;
    n = B.n;
    for (int i=0; i<REC_DIM; i++) {
      ccfl[i] = B.ccfl[i];
      have_recursion_coefs[i] = B.have_recursion_coefs[i];
      fwd_recursion_coefs[i] = B.fwd_recursion_coefs[i];
      bwd_recursion_coefs[i] = B.bwd_recursion_coefs[i];
    }    

  return *this;
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
void CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::calculate_coefficients
                               (const int &normal,
                                const CoefType *a,
                                const CoefType *ab,
                                const CoefType *sig,
                                const CoefType *sigb,
                                const CoefType &ccfl)
{

  int ind, i, n;

  // figure out the index for storage
  i = -1;
  for (i=0; i<REC_DIM; i++)
    if (sides[i] == normal)
      ind = i;

  // TODO add meaningful error checking / throw ...
  if (i != -1) {

    // save ccfl
    this->ccfl[ind] = ccfl;

    // allocate memory
    n = rec_box.get_high()[ind] - rec_box.get_low()[ind];
    fwd_recursion_coefs[ind].reserve(n*7);
    bwd_recursion_coefs[ind].reserve(n*7);
    have_recursion_coefs[ind] = true;

    // calculate the coefficients
    for (i=0; i<n; i++) {
      bwd_recursion_coefs[ind].push_back((a[i]-ccfl-sig[i])/(a[i] + ccfl + sig[i]));
      bwd_recursion_coefs[ind].push_back((a[i]+ccfl-sig[i])/(a[i] + ccfl + sig[i]));
      bwd_recursion_coefs[ind].push_back((-ab[i]+ccfl+sigb[i])/(a[i] + ccfl + sig[i]));
      bwd_recursion_coefs[ind].push_back((-ab[i]-ccfl+sigb[i])/(a[i] + ccfl + sig[i]));
      bwd_recursion_coefs[ind].push_back((-a[i]+ccfl-sig[i])/(a[i] + ccfl + sig[i]));
      bwd_recursion_coefs[ind].push_back((ab[i]+ccfl+sigb[i])/(a[i] + ccfl + sig[i]));
      bwd_recursion_coefs[ind].push_back((ab[i]-ccfl+sigb[i])/(a[i] + ccfl + sig[i]));

      fwd_recursion_coefs[ind].push_back((ab[i]-ccfl-sigb[i])/(ab[i] + ccfl + sigb[i]));
      fwd_recursion_coefs[ind].push_back((ab[i]+ccfl-sigb[i])/(ab[i] + ccfl + sigb[i]));
      fwd_recursion_coefs[ind].push_back((-a[i]+ccfl+sig[i])/(ab[i] + ccfl + sigb[i]));
      fwd_recursion_coefs[ind].push_back((-a[i]-ccfl+sig[i])/(ab[i] + ccfl + sigb[i]));
      fwd_recursion_coefs[ind].push_back((-ab[i]+ccfl-sigb[i])/(ab[i] + ccfl + sigb[i]));
      fwd_recursion_coefs[ind].push_back((a[i]+ccfl+sig[i])/(ab[i] + ccfl + sigb[i]));
      fwd_recursion_coefs[ind].push_back((a[i]-ccfl+sig[i])/(ab[i] + ccfl + sigb[i]));

      
    }
  }
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
CoefType CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::get_fwd_coef
                  (int &comp_ind, const IndexType &p, const int i) const
{
  // TODO add throw ... or some meaningful way to show that this failed ...
  if (have_recursion_coefs[comp_ind]) {
    return fwd_recursion_coefs[comp_ind][p*7+i];
  } else {
    std::cerr << "returning 0 for coeficient" << std::endl;
    return 0;
  }
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
CoefType CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::get_bwd_coef
                  (int &comp_ind, const IndexType &p, const int i) const
{
  // TODO add throw ... or some meaningful way to show that this failed ...
  if (have_recursion_coefs[comp_ind]) {
    return bwd_recursion_coefs[comp_ind][p*7+i];
  } else {
    return 0;
  }
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
CoefType CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::get_ccfl
                  (int &comp_ind) const
{
  // TODO add throw ... or some meaningful way to show that this failed ...
  if (have_recursion_coefs[comp_ind]) {
    return ccfl[comp_ind];
  } else {
    return 0;
  }
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
grid::Array<DATA_DIM, DataType, IndexType> & 
  CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::
  operator()(const Time &T, const grid::Point<REC_DIM, IndexType> &rec_index)
{
  return data[T][rec_box.get_lindex(rec_index)];
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
const grid::Array<DATA_DIM, DataType, IndexType> & 
  CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::
  operator()(const Time &T, const grid::Point<REC_DIM, IndexType> &rec_index) const
{
  return data[T][rec_box.get_lindex(rec_index)];
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
grid::Array<DATA_DIM, DataType, IndexType> & 
  CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::
  operator()(const Time &T, const std::size_t &rec_index)
{
  return data[T][rec_index];
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
const grid::Array<DATA_DIM, DataType, IndexType> & 
  CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::
  operator()(const Time &T, const std::size_t &rec_index) const
{
  return data[T][rec_index];
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
void CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::rotate_data()
{
  for (std::size_t i = 0; i<n; i++) {
    data[OLD][i].swap(data[CUR][i]);
    data[NEW][i].swap(data[CUR][i]);
  }
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
void CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::zero_new_data()
{
  for (std::size_t i = 0; i<n; i++) {
    data[NEW][i].zero();
  }
}

template <int REC_DIM, int DATA_DIM, class DataType, class IndexType, class CoefType>
std::size_t CRBC_Data<REC_DIM, DATA_DIM, DataType, IndexType, CoefType>::get_mem_usage() const
{
  std::size_t i, j, mem_use = 0;

  mem_use += sizeof(bool)*(REC_DIM+1);
  mem_use += data_box.get_mem_usage();
  mem_use += rec_box.get_mem_usage();
  for (i=0; i<3; ++i)
    for (j=0; j<data[i].size(); ++j)
      mem_use += data[i][j].get_mem_usage();
  mem_use += sides.get_mem_usage();
  mem_use += bounds.get_mem_usage();
  mem_use += sizeof(std::size_t);
  mem_use += REC_DIM*sizeof(DataType);
  for (i=0; i<REC_DIM; ++i) {
    mem_use += fwd_recursion_coefs[i].size()*sizeof(CoefType);
    mem_use += bwd_recursion_coefs[i].size()*sizeof(CoefType);
  }
  
  return mem_use;
}
    
} // end namespace

#endif  // CRBC_DATA_H_
