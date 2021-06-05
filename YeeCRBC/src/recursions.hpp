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

/* This file contains a class that implements CRBC recursions on a CRBC_Data
   data array 
*/

#ifndef RECURSIONS_H_
#define RECURSIONS_H
#include "grid/point.hpp"
#include "grid/box.hpp"
#include "grid/grid_array.hpp"
#include "boundary_properties.hpp"
#include "crbc_data.hpp"

#include<iostream>

#if USE_OPENMP
  #include <omp.h>
#endif

namespace crbc
{

/// Class to apply CRBC recursion updates to CRBC_Data objects
/// 
/// the template parameters are: \n
/// \tparam REC_DIM   the number recursions directions needed. Typically, 1 
///                     corresponds to a face, 2 to an edge, 3 to a corner, etc.
/// \tparam DATA_DIM  the number of physical dimensions for the problem 
/// \tparam DataType  The type of data to use. This isn't checked, but it
///                     should be a float type. Note, we calculate the cosines
///                     using double precision.
/// \tparam IndexType The data type to use for indexing. This isn't checked, but
///                     should be an integral type. This needs to be larger
///                     enough to store the number of points each dimension.
/// \tparam CoefType   The data type used for the coefficients such as time step
///                    size, wave speed, etc.
template <int DIM = 3, 
          class DataType = double, 
          class IndexType = long int,
          class CoefType = DataType>
class CrbcRecursions : virtual public BoundaryProperties
{

  public:

    /// Constructor --- default
    CrbcRecursions() {};

    /// destructor --- declared so we can make it virtual
    virtual ~CrbcRecursions() {};

    /// Routine to calculate the forward recursive updates in the provided 
    /// direction
    ///
    /// \tparam REC_DIM           Number of recursion directions, e.g. 1 for a 
    ///                           face, 2 for an edge, 3 for a corner
    /// \param[inout] crbc_data   A CRBC_Data array that we want to run 
    ///                           recursions over
    /// \param[in]    coord       The coordinate direction we want to run the 
    ///                           the recursions in. We use 0-based indexing, 
    ///                           so 0 corresponds to the first component, 1
    ///                           corresponds to the second component and so on.
    ///                           For example in 3D, 0 corresponds to x, 1 to y
    ///                           and 2 to z.
    template <int REC_DIM>
    void forward_recursions(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data,
                            const int &coord);

    /// Routine to calculate the backward recursive updates in the provided 
    /// direction
    ///
    /// \tparam REC_DIM           Number of recursion directions, e.g. 1 for a 
    ///                           face, 2 for an edge, 3 for a corner
    /// \param[inout] crbc_data   A CRBC_Data array that we want to run 
    ///                           recursions over
    /// \param[in]    coord       The coordinate direction we want to run the 
    ///                           the recursions in. We use 0-based indexing, 
    ///                           so 0 corresponds to the first component, 1
    ///                           corresponds to the second component and so on.
    ///                           For example in 3D, 0 corresponds to x, 1 to y
    ///                           and 2 to z.
    template <int REC_DIM>
    void backward_recursions(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data,
                             const int &coord);

    /// Routine to calculate the termination conditions in the provided
    /// direction
    ///
    /// \tparam REC_DIM           Number of recursion directions, e.g. 1 for a 
    ///                           face, 2 for an edge, 3 for a corner
    /// \param[inout] crbc_data   A CRBC_Data array that we want to run 
    ///                           recursions over
    /// \param[in]    coord       The coordinate direction we want to run the 
    ///                           the recursions in. We use 0-based indexing, 
    ///                           so 0 corresponds to the first component, 1
    ///                           corresponds to the second component and so on.
    ///                           For example in 3D, 0 corresponds to x, 1 to y
    ///                           and 2 to z.
    template <int REC_DIM>
    void terminations(CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data,
                              const int &coord);


}; // end class

//------------------------------------------------------------------------------
//                                 Definitions
//------------------------------------------------------------------------------
template <int DIM, class DataType, class IndexType, class CoefType>
template <int REC_DIM>
void CrbcRecursions<DIM, DataType, IndexType, CoefType>::forward_recursions(
                            CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data,
                            const int &coord)
{
  typedef CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> CRBC_Data;
  
  int comp_ind = 0;
  int side;
  int i;
  std::size_t j, l, ind, shift, rind, rshift;
  grid::Point<REC_DIM, IndexType> rec_nvec, r;
  grid::Point<DIM, IndexType> nvec, d;
  grid::Box<REC_DIM, IndexType> rec_box;
  grid::Box<DIM, IndexType> data_box;

  // figure out which component of the recursion indexing corresponds to the
  // direction we want to run the recursions in
  for (i = 0; i<REC_DIM; i++)
    comp_ind = (crbc_data.get_sides()[i] / 2 == coord) ? i : comp_ind;

  // save the boundary side -- this is the outward normal direction for the 
  // recursions
  side = crbc_data.get_sides()[comp_ind];

  // create a box that goes from p=1,...,P in the recursion direction and 
  // p = 0,...,P in all other directions
  rec_nvec = grid::Point<REC_DIM, IndexType>();
  rec_nvec.set_value(comp_ind, 1);
  rec_box = crbc_data.get_rec_box().add_low(rec_nvec);

  // figure out the data box
  data_box = crbc_data.get_data_box().get_side(coord, (side % 2 == 1));

  // shrink the box in the tangential direction(s) unless the boundary condition
  // is Neumann. We don't need to update Dirichlet boundaries because they will
  // already be correctly set to 0. CRBC type boundaries will be updated in an
  // intersecting layer (e.g. edges or corners). Finally, we never update the 
  // boundary type NONE.
  nvec = grid::Point<DIM, IndexType>();
  for (i = 0; i<DIM; i++) {
    if ((crbc_data.get_boundaries()[2*i] != NEUM) && (i != coord))
      nvec.set_value(i, 1);
  }
  data_box = data_box.add_low(nvec);

  nvec = grid::Point<DIM, IndexType>();
  for (i = 0; i<DIM; i++) {
    if ((crbc_data.get_boundaries()[2*i+1] != NEUM) && (i != coord))
      nvec.set_value(i, -1);
  }
  data_box = data_box.add_high(nvec);

  // make an outward normal unit vector for the spatial relatons
  nvec = grid::Point<DIM, IndexType>();
  if (side % 2 == 0) { // left side
    nvec.set_value(coord, -1);
  } else { // right side
    nvec.set_value(coord, 1);
  }

  // figure out what the normal vector is in terms of logical indexing
  rshift = crbc_data.get_rec_box().get_lstride(rec_nvec);

  // loop over the recursion points
  #if USE_OPENMP
    #pragma omp parallel default(shared) private(d, r, j, l, ind, rind)
    {
  #endif
  for (l = 0; l < rec_box.get_npoints(); ++l) {

    r = rec_box.get_pindex(l);

    rind = crbc_data.get_rec_box().get_lindex(r);

    // figure out what the normal vector is in terms of logical indexing
    shift = crbc_data(CRBC_Data::NEW, rind).get_box().get_lstride(nvec);
    
    // loop over data points
    #if USE_OPENMP
      #pragma omp for
    #endif
    for (j=0; j<data_box.get_npoints(); j++) {

      // get the physical index
      d = data_box.get_pindex(j);

      // now get the logical index in the original data box
      ind = crbc_data(CRBC_Data::NEW, rind).get_box().get_lindex(d);

      // update using the forward recursions. The update expression for the
      // x direction in 3D is as follows.
      //
      // u->set_new(i,j,k,p) = b[p-1][0] * u->get_cur(i,j,k,p) 
      //                     + b[p-1][1] * u->get_cur(i+1,j,k,p) 
      //                     + b[p-1][2] * u->get_cur(i+1,j,k,p-1)
      //                     + b[p-1][3] * u->get_cur(i,j,k,p-1)
      //                     + b[p-1][4] * u->get_new(i+1,j,k,p)
      //                     + b[p-1][5] * u->get_new(i+1,j,k,p-1) 
      //                     + b[p-1][6] * u->get_new(i,j,k,p-1),
      // where
      //   b[i][0] = (ab[i]-ccflx-sigb[i])/(ab[i] + ccflx + sigb[i]);
      //   b[i][1] = (ab[i]+ccflx-sigb[i])/(ab[i] + ccflx + sigb[i]);
      //   b[i][2] = (-a[i]+ccflx+sig[i])/(ab[i] + ccflx + sigb[i]);
      //   b[i][3] = (-a[i]-ccflx+sig[i])/(ab[i] + ccflx + sigb[i]);
      //   b[i][4] = (-ab[i]+ccflx-sigb[i])/(ab[i] + ccflx + sigb[i]);
      //   b[i][5] = (a[i]+ccflx+sig[i])/(ab[i] + ccflx + sigb[i]);
      //   b[i][6] = (a[i]-ccflx+sig[i])/(ab[i] + ccflx + sigb[i]);
      //
      // ccflx, is the c*dt/hx and the remaining terms come for the optimal
      // cosines.

      crbc_data(CRBC_Data::NEW, rind)[ind] = \
         crbc_data.get_fwd_coef(comp_ind, r[comp_ind]-1, 0) * crbc_data(CRBC_Data::CUR, rind)[ind] \
       + crbc_data.get_fwd_coef(comp_ind, r[comp_ind]-1, 1) * crbc_data(CRBC_Data::CUR, rind)[ind+shift] \
       + crbc_data.get_fwd_coef(comp_ind, r[comp_ind]-1, 2) * crbc_data(CRBC_Data::CUR, rind-rshift)[ind+shift] \
       + crbc_data.get_fwd_coef(comp_ind, r[comp_ind]-1, 3) * crbc_data(CRBC_Data::CUR, rind-rshift)[ind] \
       + crbc_data.get_fwd_coef(comp_ind, r[comp_ind]-1, 4) * crbc_data(CRBC_Data::NEW, rind)[ind+shift] \
       + crbc_data.get_fwd_coef(comp_ind, r[comp_ind]-1, 5) * crbc_data(CRBC_Data::NEW, rind-rshift)[ind+shift] \
       + crbc_data.get_fwd_coef(comp_ind, r[comp_ind]-1, 6) * crbc_data(CRBC_Data::NEW, rind-rshift)[ind];
    }
  }
  #if USE_OPENMP
    }
  #endif
}

template <int DIM, class DataType, class IndexType, class CoefType>
template <int REC_DIM>
void CrbcRecursions<DIM, DataType, IndexType, CoefType>::backward_recursions(
                            CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data,
                            const int &coord)
{

  typedef CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> CRBC_Data;

  int comp_ind = 0;
  int side;
  int i;
  std::size_t j, ind, shift, rind, rshift;
  grid::Point<REC_DIM, IndexType> rec_nvec;
  grid::Point<DIM, IndexType> nvec, d;
  grid::Box<REC_DIM, IndexType> rec_box;
  grid::Box<DIM, IndexType> data_box;
  grid::BoxIterator<REC_DIM, IndexType> rec;
  grid::BoxIterator<DIM, IndexType> data;
  

  // figure out which component of the recursion indexing corresponds to the
  // direction we want to run the recursions in
  for (i = 0; i<REC_DIM; i++)
    comp_ind = (crbc_data.get_sides()[i] / 2 == coord) ? i : comp_ind;

  // save the boundary side -- this is the outward normal direction for the 
  // recursions
  side = crbc_data.get_sides()[comp_ind];

  // create a box that goes from p=0,...,P-1 in the recursion direction and 
  // p = 0,...,P in all other directions
  rec_nvec = grid::Point<REC_DIM, IndexType>();
  rec_nvec.set_value(comp_ind, -1);
  rec_box = crbc_data.get_rec_box().add_high(rec_nvec);
  rec_nvec.set_value(comp_ind, 1);

  // create an iterator for the recursions 
  rec = grid::BoxIterator<REC_DIM, IndexType> (rec_box);

  // figure out the data box
  data_box = crbc_data.get_data_box().get_side(coord, (side % 2 == 0));

  // shrink the box in the tangential direction(s) unless the boundary condition
  // is Neumann. We don't need to update Dirichlet boundaries because they will
  // already be correctly set to 0. CRBC type boundaries will be updated in an
  // intersecting layer (e.g. edges or corners). Finally, we never update the 
  // boundary type NONE.
  nvec = grid::Point<DIM, IndexType>();
  for (i = 0; i<DIM; i++) {
    if ((crbc_data.get_boundaries()[2*i] != NEUM) && (i != coord))
      nvec.set_value(i, 1);
  }
  data_box = data_box.add_low(nvec);

  nvec = grid::Point<DIM, IndexType>();
  for (i = 0; i<DIM; i++) {
    if ((crbc_data.get_boundaries()[2*i+1] != NEUM) && (i != coord))
      nvec.set_value(i, -1);
  }
  data_box = data_box.add_high(nvec);

  // make an outward normal unit vector for the spatial relatons
  nvec = grid::Point<DIM, IndexType>();
  if (side % 2 == 0) { // left side
    nvec.set_value(coord, -1);
  } else { // right side
    nvec.set_value(coord, 1);
  }

  // figure out what the normal vector is in terms of logical indexing
  rshift = crbc_data.get_rec_box().get_lstride(rec_nvec);

  // loop over the recursion points
  for (rec.end(); rec.notDone(); --rec) {

    // figure out what the normal vector is in terms of logical indexing
    shift = crbc_data(CRBC_Data::NEW, rec()).get_box().get_lstride(nvec);

    rind = crbc_data.get_rec_box().get_lindex(rec());

    #if USE_OPENMP
      #pragma omp parallel for default(shared) private(d, j, ind)
    #endif
    for (j=0; j<data_box.get_npoints(); j++) {

      // get the physical index
      d = data_box.get_pindex(j);

      // now get the logical index in the original data box
      ind = crbc_data(CRBC_Data::NEW, rind).get_box().get_lindex(d);

      // update using the backward recursions. The update expression for the
      // x direction in 3D is as follows.
      //
      // u->set_new(i,j,k,p) = b[p][0] * u->get_cur(i,j,k,p) 
      //                     + b[p][1] * u->get_cur(i-1,j,k,p) 
      //                     + b[p][2] * u->get_cur(i-1,j,k,p+1)
      //                     + b[p][3] * u->get_cur(i,j,k,p+1)
      //                     + b[p][4] * u->get_new(i-1,j,k,p)
      //                     + b[p][5] * u->get_new(i-1,j,k,p+1) 
      //                     + b[p][6] * u->get_new(i,j,k,p+1),
      // where
      //   b[i][0] = (a[i]-ccflx-sig[i])/(a[i] + ccflx + sig[i]);
      //   b[i][1] = (a[i]+ccflx-sig[i])/(a[i] + ccflx + sig[i]);
      //   b[i][2] = (-ab[i]+ccflx+sigb[i])/(a[i] + ccflx + sig[i]);
      //   b[i][3] = (-ab[i]-ccflx+sigb[i])/(a[i] + ccflx + sig[i]);
      //   b[i][4] = (-a[i]+ccflx-sig[i])/(a[i] + ccflx + sig[i]);
      //   b[i][5] = (ab[i]+ccflx+sigb[i])/(a[i] + ccflx + sig[i]);
      //   b[i][6] = (ab[i]-ccflx+sigb[i])/(a[i] + ccflx + sig[i]);
      //
      // ccflx, is the c*dt/hx and the remaining terms come for the optimal
      // cosines.
      crbc_data(CRBC_Data::NEW, rind)[ind] = \
         crbc_data.get_bwd_coef(comp_ind, rec()[comp_ind], 0) * crbc_data(CRBC_Data::CUR, rind)[ind] \
       + crbc_data.get_bwd_coef(comp_ind, rec()[comp_ind], 1) * crbc_data(CRBC_Data::CUR, rind)[ind-shift] \
       + crbc_data.get_bwd_coef(comp_ind, rec()[comp_ind], 2) * crbc_data(CRBC_Data::CUR, rind+rshift)[ind-shift] \
       + crbc_data.get_bwd_coef(comp_ind, rec()[comp_ind], 3) * crbc_data(CRBC_Data::CUR, rind+rshift)[ind] \
       + crbc_data.get_bwd_coef(comp_ind, rec()[comp_ind], 4) * crbc_data(CRBC_Data::NEW, rind)[ind-shift] \
       + crbc_data.get_bwd_coef(comp_ind, rec()[comp_ind], 5) * crbc_data(CRBC_Data::NEW, rind+rshift)[ind-shift] \
       + crbc_data.get_bwd_coef(comp_ind, rec()[comp_ind], 6) * crbc_data(CRBC_Data::NEW, rind+rshift)[ind];
    }
  }
}

template <int DIM, class DataType, class IndexType, class CoefType>
template <int REC_DIM>
void CrbcRecursions<DIM, DataType, IndexType, CoefType>::terminations(
                            CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> &crbc_data,
                            const int &coord)
{

  typedef CRBC_Data<REC_DIM,DIM,DataType,IndexType,CoefType> CRBC_Data;

  int comp_ind = 0;
  int side;
  int i;
  std::size_t j, l, ind, shift, rind;
  grid::Point<REC_DIM, IndexType> rec_nvec, r;
  grid::Point<DIM, IndexType> nvec, d;
  grid::Box<REC_DIM, IndexType> rec_box;
  grid::Box<DIM, IndexType> data_box;
  grid::BoxIterator<DIM, IndexType> data;
  

  // figure out which component of the recursion indexing corresponds to the
  // direction we want to run the recursions in
  for (i = 0; i<REC_DIM; i++)
    comp_ind = (crbc_data.get_sides()[i] / 2 == coord) ? i : comp_ind;

  // save the boundary side -- this is the outward normal direction for the 
  // recursions
  side = crbc_data.get_sides()[comp_ind];

  // create a box that goes from p=P in the recursion direction and 
  // p = 0,...,P in all other directions
  rec_nvec = grid::Point<REC_DIM, IndexType>();
  rec_nvec.set_value(comp_ind, 1);
  rec_box = crbc_data.get_rec_box().get_side(comp_ind, false);

  // figure out the data box
  data_box = crbc_data.get_data_box().get_side(coord, (side % 2 == 0));

  // shrink the box in the tangential direction(s) unless the boundary condition
  // is Neumann. We don't need to update Dirichlet boundaries because they will
  // already be correctly set to 0. CRBC type boundaries will be updated in an
  // intersecting layer (e.g. edges or corners). Finally, we never update the 
  // boundary type NONE.
  nvec = grid::Point<DIM, IndexType>();
  for (i = 0; i<DIM; i++) {
    if ((crbc_data.get_boundaries()[2*i] != NEUM) && (i != coord))
      nvec.set_value(i, 1);
  }
  data_box = data_box.add_low(nvec);

  nvec = grid::Point<DIM, IndexType>();
  for (i = 0; i<DIM; i++) {
    if ((crbc_data.get_boundaries()[2*i+1] != NEUM) && (i != coord))
      nvec.set_value(i, -1);
  }
  data_box = data_box.add_high(nvec);

  // make an outward normal unit vector for the spatial relatons
  nvec = grid::Point<DIM, IndexType>();
  if (side % 2 == 0) { // left side
    nvec.set_value(coord, -1);
  } else { // right side
    nvec.set_value(coord, 1);
  }

  // loop over the recursion points
  #if USE_OPENMP
    #pragma omp parallel default(shared) private(d, r, j, l, ind, rind)
    {
  #endif
  for (l = 0; l < rec_box.get_npoints(); ++l) {

    r = rec_box.get_pindex(l);

    rind = crbc_data.get_rec_box().get_lindex(r);

    // figure out what the normal vector is in terms of logical indexing
    shift = crbc_data(CRBC_Data::NEW, rind).get_box().get_lstride(nvec);
    
    // loop over data points
    #if USE_OPENMP
      #pragma omp for
    #endif
    for (j=0; j<data_box.get_npoints(); j++) {

      // get the physical index
      d = data_box.get_pindex(j);

      // now get the logical index in the original data box
      ind = crbc_data(CRBC_Data::NEW, rind).get_box().get_lindex(d);

      // update using the termination conditions. The update expression for the
      // x direction in 3D is as follows. (this is Sommerfeld)
      //
      // u->set_new(i,j,k,p) = (u->get_cur(i,j,k,p) 
      //                       + u->get_cur(i-1,j,k,p) 
      //                       - u->get_new(i-1,j,k,p)
      //                       + ccflx*(u->get_new(i-1,j,k,p) 
      //                       - u->get_cur(i,j,k,p) 
      //                       + u->get_cur(i-1,j,k,p)))
      //                       / (1.0 + ccflx);
      //
      // ccflx, is the c*dt/hx 
      crbc_data(CRBC_Data::NEW, rind)[ind] = ( \
         crbc_data(CRBC_Data::CUR, rind)[ind] \
       + crbc_data(CRBC_Data::CUR, rind)[ind-shift] \
       - crbc_data(CRBC_Data::NEW, rind)[ind-shift] \
       + crbc_data.get_ccfl(comp_ind) * (crbc_data(CRBC_Data::NEW, rind)[ind-shift] \
       - crbc_data(CRBC_Data::CUR, rind)[ind] \
       + crbc_data(CRBC_Data::CUR, rind)[ind-shift] \
       )) / (1.0 + crbc_data.get_ccfl(comp_ind));
    }
  }
  #if USE_OPENMP
    }
  #endif
}

} // end namespace
#endif //define RECURSIONS_H
