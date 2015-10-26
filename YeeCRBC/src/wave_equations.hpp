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

/* This file contains a class that implements wave equation updates on a 
   grid/Array object
*/

#ifndef WAVE_EQUATIONS_H_
#define WAVE_EQUATIONS_H
#include "grid/point.hpp"
#include "grid/box.hpp"
#include "grid/grid_array.hpp"
#include "boundary_properties.hpp"

#if USE_OPENMP
  #include <omp.h>
#endif

namespace crbc {

/// Class to apply wave equation updates to grid::Array objects
/// 
/// the template parameters are: \n
/// \tparam DIM       the number recursions directions needed. Typically, 1 
///                     corresponds to a face, 2 to an edge, 3 to a corner, etc.
/// \tparam DataType  The type of data to use. This isn't checked, but it
///                     should be a float type. Note, we calculate the cosines
///                     using double precision.
/// \tparam IndexType The data type to use for indexing. This isn't checked, but
///                     should be an integral type. This needs to be larger
///                     enough to store the number of points each dimension.
template <int DIM = 3, 
          class DataType = double, 
          class IndexType = long int,
          class CoefType = DataType>
class WaveEquationUpdates : virtual public BoundaryProperties
{

  public:

    /// constructor --- empty
    WaveEquationUpdates() {init_params();};

    /// contructor --- load basic parameters, grid spacing, wave speed, and time
    /// step
    ///
    /// \param[in] grid_spacing the grid spacings in each dimension (hx, hy, ...)
    /// \param[in] c            the wave speed
    /// \param[in] dt           time step size
    WaveEquationUpdates(const grid::Point<DIM, CoefType> &grid_spacing,
                        const CoefType &c,
                        const CoefType &dt);

    /// method to set the wave speed
    /// \param[in] c wave speed
    void set_wave_speed(const CoefType &c) {this->c = c; calc_coefs();};

    /// method to set the time step
    /// \param[in] dt time step
    void set_time_step(const CoefType &dt) {this->dt = dt; calc_coefs();};
 
    /// method to set the grid spacings
    /// \param[in] h grid spacings
    void set_grid_spacing(const grid::Point<DIM, CoefType> h) {this->h = h; calc_coefs();};

    /// destrucotr, empty
    virtual ~WaveEquationUpdates() {};

    /// function to apply the wave equation update using the standard second
    /// order central differences in space and time. We assume that the
    /// boundaries are included in the grid
    /// 
    /// \param[in] boundaries Vector of boundary conditions
    /// \param[in] old_data Array of data at old time t_{n-1}
    /// \param[in] cur_data Array of data at current time t_{n}
    /// \param[out] new_data updated array of data at time t_{n+1}
    void get_wave_eq_update(const grid::Point<2*DIM, Boundary> &boundaries,
                       const grid::Array<DIM,DataType,IndexType> &old_data,
                       const grid::Array<DIM,DataType,IndexType> &cur_data,
                       grid::Array<DIM,DataType,IndexType> &new_data) const;

  protected:

    /// Routine to get the wave coefficients. This is protected so the coefs.
    /// can be used by a derived class
    ///
    /// \param[in] i index of coefficient to get
    /// \return coefs[i]
    inline CoefType get_coef(const int &i) const {return coefs[i];};

  private:
    
    CoefType dt, c;
    grid::Point<DIM, CoefType> h;
    grid::Point<DIM, CoefType> coefs;

    /// function to intilize all of the parameters to 1
    void init_params();

    /// function to calculate the coefficients used in the wave equation 
    /// updates (dt^2*c^2/h^2)
    void calc_coefs();

    /// function that actually does the wave equation updates. Again, we assume
    /// that the boundaries are included in the grid. This is virtual because
    /// we want to allow some flexibility to evolve other wave equations (e.g.
    /// the Klein-Gordon equation)
    /// 
    /// \param[in] boundaries Vector of boundary conditions
    /// \param[in] old_data Array of data at old time t_{n-1}
    /// \param[in] cur_data Array of data at current time t_{n}
    /// \param[out] new_data updated array of data at time t_{n+1}
    virtual void wave_eq_update(const grid::Point<2*DIM, Boundary> &boundaries,
                           const grid::Array<DIM,DataType,IndexType> &old_data,
                           const grid::Array<DIM,DataType,IndexType> &cur_data,
                           grid::Array<DIM,DataType,IndexType> &new_data) const;

    /// function to apply Neumann type boundary conditions.
    /// 
    /// \param[in] boundaries Vector of boundary conditions
    /// \param[in] old_data Array of data at old time t_{n-1}
    /// \param[in] cur_data Array of data at current time t_{n}
    /// \param[out] new_data updated array of data at time t_{n+1}
    virtual void apply_neumman_boundaries(
                             const grid::Point<2*DIM, Boundary> &boundaries,
                             const grid::Array<DIM,DataType,IndexType> &old_data,
                             const grid::Array<DIM,DataType,IndexType> &cur_data,
                             grid::Array<DIM,DataType,IndexType> &new_data) const;
  

}; // end class

//------------------------------------------------------------------------------
//                         constructor
//------------------------------------------------------------------------------
template<int DIM, class DataType, class IndexType, class CoefType>
WaveEquationUpdates<DIM, DataType, IndexType,CoefType>::WaveEquationUpdates
                       (const grid::Point<DIM, CoefType> &grid_spacing,
                        const CoefType &c,
                        const CoefType &dt)
{

  this->c = c;
  this->dt = dt;
  h = grid_spacing;

  // calculate the scalings used in the discretization
  calc_coefs();

}

//------------------------------------------------------------------------------
//                         Public Definitions
//------------------------------------------------------------------------------
template<int DIM, class DataType, class IndexType, class CoefType>
void WaveEquationUpdates<DIM, DataType, IndexType,CoefType>::get_wave_eq_update
                      (const grid::Point<2*DIM, Boundary> &boundaries,
                       const grid::Array<DIM,DataType,IndexType> &old_data,
                       const grid::Array<DIM,DataType,IndexType> &cur_data,
                       grid::Array<DIM,DataType,IndexType> &new_data) const
{
  wave_eq_update(boundaries, old_data, cur_data, new_data);
}

//------------------------------------------------------------------------------
//                         Private Definitions
//------------------------------------------------------------------------------
template<int DIM, class DataType, class IndexType, class CoefType>
void  WaveEquationUpdates<DIM, DataType, IndexType,CoefType>::calc_coefs()
{
  for (int i=0; i<DIM; i++)
    coefs.set_value(i, dt*dt*c*c/h[i]/h[i]);
}

template<int DIM, class DataType, class IndexType, class CoefType>
void  WaveEquationUpdates<DIM, DataType, IndexType,CoefType>::init_params()
{
  dt = c = static_cast<CoefType>(1.0);
  h = grid::Point<DIM, CoefType>::ones();
  
}

template<int DIM, class DataType, class IndexType, class CoefType>
void WaveEquationUpdates<DIM, DataType, IndexType,CoefType>::wave_eq_update(
                           const grid::Point<2*DIM, Boundary> &boundaries,
                           const grid::Array<DIM,DataType,IndexType> &old_data,
                           const grid::Array<DIM,DataType,IndexType> &cur_data,
                           grid::Array<DIM,DataType,IndexType> &new_data) const
{
  int j;
  std::size_t k, ind, shift[3];

  grid::Point<DIM, IndexType> e_j;

  // get the box indexing the data. TODO, we should really check that the new,
  // current, and old data arrays are defined over the same box
  grid::Box<DIM, IndexType> box = new_data.get_box();

  // shrink the box so we only loop over the interal points that we can apply
  // the full wave equation stencil to. Note that we have assumed that the 
  // boundaries have been included in the grid
  box = box.grow(-1);

  // grid::BoxIterator<DIM, IndexType> i;
  grid::Point<DIM, IndexType> i;

  // figure out what the normal vector is in terms of logical indexing
  for (j=0; j<DIM; j++) {
    // create a the standard unit vector e_j
    e_j = grid::Point<DIM, IndexType>::unit(j);
    shift[j] = new_data.get_box().get_lstride(e_j);
  }

  // loop over all of the data boxes
  #if USE_OPENMP
    #pragma omp parallel for default(shared) private(k, i, j, e_j, ind)
  #endif
  for (k=0; k<box.get_npoints(); ++k) {

    // get the physical point
    i = box.get_pindex(k);
   
    // now get the logical index in the original data box
    ind = new_data.get_box().get_lindex(i);

    // unew = 2*ucur - uold
    new_data[ind] = 2*cur_data[ind] - old_data[ind];  

    // loop over the dimensions
    for (j=0; j<DIM; j++) {

      // apply the 1d 3 point WE stencil in the current dimension
      // adds dt^2*c^2/h_{x_i}^2 *(u_{x_i-1} - 2u_{x_i} + u_{x_i+1})
      new_data[ind] += coefs[j]*(cur_data[ind - shift[j]] - 2*cur_data[ind] \
                                 + cur_data[ind + shift[j]] );
    }
  }

  // now apply boundary updates.
  // In our setting, Dirichlet boundaries should always be correct since we 
  // never update them. We will not update the NONE boundary type nor will we
  // update the CRBC type boundaries. Therefore, we only need to update the 
  // Neumann boundaries
  apply_neumman_boundaries(boundaries, old_data, cur_data, new_data);

}

template<int DIM, class DataType, class IndexType, class CoefType>
void WaveEquationUpdates<DIM, DataType, IndexType,CoefType>::apply_neumman_boundaries
                         (const grid::Point<2*DIM, Boundary> &boundaries,
                          const grid::Array<DIM,DataType,IndexType> &old_data,
                          const grid::Array<DIM,DataType,IndexType> &cur_data,
                          grid::Array<DIM,DataType,IndexType> &new_data) const
{
  int i,j,k,l,m, count, *neum_ind, *int_ind, int_count, den;
  int num_neum = 0;
  int intersects = 0;
  bool found_intersect = false, test;
  grid::Box<DIM, IndexType> w_box;
  grid::BoxIterator<DIM, IndexType> iter; 
  grid::Point<DIM, IndexType> e_j;

  // find how many boundary sides are Neumann
  for (i=0; i<2*DIM; i++)
    num_neum += (boundaries[i] == NEUM) ? 1 : 0;

  // calclate the indices of the neumann sides
  neum_ind = new int[num_neum];
  count = 0;
  for (i=0; i<2*DIM; i++) {
    if (boundaries[i] == NEUM) {
      neum_ind[count] = i;
      ++count;
    }
  }

  // figure out the maximum number of intersecting Neumann boundaries
  intersects = (num_neum > DIM) ? DIM : num_neum;
  int_ind = new int[intersects];

  // loop over the possible number of intersections
  for (i=0; i<intersects; i++) {
    
    // set flag so we can stop early if all of the intersections are empty
    if (i > 0)
      found_intersect = false;

    // loop over sides
    for (j=0; j<num_neum - i; j++) {
 
      // calulate the number of "intersections" we need to check. These are 
      // i-dimensional triangular numbers so we will use the recursive formula
      // t^{d}_{k} = \frac{k+d-1}{d} t^{d-1}_{k}
      // t^{0}_{k} = 1
      int_count = 1;
      den = 1;
      for (k=1; k<=i; k++) {
        int_count *= num_neum - j - k;
        den *= i - k + 1;
      }
      int_count /= den;

      // reset the index array
      int_ind[0] = j;
      for (k=1; k<=i; k++)
        int_ind[k] = j+k;
      if (i > 0)
        int_ind[i] = int_ind[i]-1;

      // now calculate the needed intersections
      for (k = 0; k<int_count; k++) {
        
        // the box at side j 
        w_box = new_data.get_box().get_side(neum_ind[j]/2, (neum_ind[j] % 2 == 0));
       
        // now find the indexes of the boxes to check for intersections ...
        for (l=i; l>0; l--) {
          
          if (int_ind[l] < num_neum - 1) {
            int_ind[l] = int_ind[l] + 1;
            for (m = l+1; m<= i; m++)
              int_ind[m] = int_ind[m-1] + 1;
            break;
          }
        }

        // calculate intersections
        for (l=1; l<=i; l++) {
          w_box = w_box.intersection(new_data.get_box().get_side(neum_ind[int_ind[l]]/2, (neum_ind[int_ind[l]] % 2 == 0)));
          if (w_box.isEmpty())
            break;
        }

        // if box is empty don't do anything
        if (w_box.isEmpty())
            continue;

        found_intersect = true;

        // shrink the box in the tangetial directions
        grid::Point<DIM, IndexType> diff = grid::Point<DIM, IndexType>::ones();
        for (l=0; l<=i; l++)
          diff.set_value(neum_ind[int_ind[l]]/2, 0);

        w_box = w_box.add_low(diff);
        diff.scale(-1);
        w_box = w_box.add_high(diff);

        // create an iterator
        iter = grid::BoxIterator<DIM, IndexType>(w_box);

        // loop over the points in the box and compute the differences from the time 
        // derivatives
        // unew = 2*ucur - uold
        for (iter.begin(); iter.notDone(); ++iter)
          new_data(iter()) = 2*cur_data(iter()) - old_data(iter()); 

        // loop over the dimensions
        for (l=0; l<DIM; l++) {
          // create a the standard unit vector e_j
          e_j = grid::Point<DIM, IndexType>::unit(l);

          test = false;
          for (m=0; m<=i; m++) {
            if (neum_ind[int_ind[m]]/2 == l) {
              test = true;
              if (neum_ind[int_ind[m]]%2 == 1)
                e_j.scale(-1);
            }
          }

          if (test) { 
            // there is a Neumann boundary in this dimension so we have
            // u(i + e_j) = u(i)

            for (iter.begin(); iter.notDone(); ++iter)
              new_data(iter()) += coefs[l]*(cur_data(iter() + e_j) \
                                 - cur_data(iter()));
           
          } else {
            // loop over the points in the box and apply the 1d 3 point WE stencil in 
            // the current dimension
            // adds dt^2*c^2/h_{x_i}^2 *(u_{x_i-1} - 2u_{x_i} + u_{x_i+1})
            for (iter.begin(); iter.notDone(); ++iter)
              new_data(iter()) += coefs[l]*(cur_data(iter() - e_j) \
                                 - 2*cur_data(iter()) \
                                 + cur_data(iter() + e_j) );
          }

        } // end l loop

      } // end k loop    

    } // end j loop    

    // break out of loop if we didn't find any non-empty intersections
    if (!found_intersect)
      break;
  } // end i loop

  delete[] neum_ind;
  delete[] int_ind;

}
} // end namespace

#endif //WAVE_EQUATIONS_H
