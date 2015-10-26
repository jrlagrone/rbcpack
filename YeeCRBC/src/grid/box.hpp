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

/* 
  Header file to define a box class and some associated helper classes
  for handing indexing bounds

*/
#include <cstddef>   // for std::size_t
#include <iostream>  // cout, etc.
#include "point.hpp"

#ifndef GRID_BOX_H_
#define GRID_BOX_H_

namespace grid {

/// Class to define a "box" that contains all of the indexing information for
/// a regular grid of points
///
/// \tparam DIM       Number of dimensions to make the box
/// \tparam DataType  Type of data to use. For the indended use of indexing, an
///                   integral type is appropriate. Other uses are untested and
///                   may cause problems.
template <int DIM, class DataType = long int>
class Box
{

  public:

    /// constructor --- defines an empty box
    Box() {empty = true;};

    /// constructor --- defines a rectangular box with the lower corner at
    /// low and upper corner at high using the default ordering (column major).
    ///
    /// \param low   The bottom left corner of the box
    /// \param high  The top right corner of the box
    Box(const Point<DIM, DataType> &low, const Point<DIM, DataType> &high);

    /// constructor --- defines a rectangular box with the lower corner at
    /// low and upper corner at high with the provided data ordering.
    /// 
    /// The point ordering gives the priority of the ordering, where
    /// 0 is first and (DIM-1) is last. For example, for DIM = 3, the point
    /// (0,1,2) corresponds to the default ordering (i + j*nx + k*nx*ny). The 
    /// point (2,0,1) would give the ordering (k + i*nz + j*nx*nz)
    ///
    /// \param low      The bottom left corner of the box
    /// \param high     The top right corner of the box
    /// \param ordering A point giving the ordering priority. This can be used
    ///                 to give column or row major ordering or something
    ///                 inbetween
    Box(const Point<DIM, DataType> &low, 
        const Point<DIM, DataType> &high, 
        const Point<DIM, DataType> &ordering);

    /// routine to return the lower corner
    inline Point<DIM, DataType> get_low() const {return low_corner;};

    /// routine to return the upper corner
    inline Point<DIM, DataType> get_high() const {return high_corner;};

    /// routine to return the ordering
    inline Point<DIM, DataType> get_ordering() const {return order;};

    /// routine to change the ordering
    void set_ordering(const Point<DIM, DataType> &ordering); 
    
    /// routine to get the number of points
    std::size_t get_npoints() const;

    /// routine to return the data strides
    Point<DIM-1, DataType> get_strides() const {return strides;};

    /// routine to return if the box is empty or not
    inline bool isEmpty() const {return empty;};

    /// function to increase/decrease the size of a box
    /// e.g. grow(-1) shrinks the box by adding one to each component of
    /// the lower point and subtracting one from the upper point
    /// This is useful to get indexing that only goes over the interior points.
    Box<DIM, DataType> grow(const DataType &i);

    /// function to increase/decrease the size of a box by changing the high 
    /// corner e.g. grow_high(-1) shrinks the box by subtracting one from the 
    /// upper point.
    Box<DIM, DataType> grow_high(const DataType &i);

    /// function to increase/decrease the size of a box by changing the high 
    /// corner e.g. grow_low(-1) shrinks the box by adding one to the lower
    /// point.
    Box<DIM, DataType> grow_low(const DataType &i);

    /// function to change the size of a box by adding the point P to the 
    /// high corner
    Box<DIM, DataType> add_high(const Point<DIM, DataType> &P);

    /// function to change the size of a box by adding the point P to the 
    /// low corner
    Box<DIM, DataType> add_low(const Point<DIM, DataType> &P);    

    /// function to partition the box into smaller disjoing boxes
    /// The components of n tell us how many divisions to make in each
    /// direction, the output is stored in part_boxes and should be sufficiently
    /// large to store <n,n> new boxes.
    void partition(const Point<DIM, DataType> &n, Box<DIM, DataType> *part_boxes);

    /// return the box along the given dimension
    /// \param i index corresponding to the requested dimension
    /// \param low true returns the left side in the requested dimension and
    ///            false returns the right side in the requested dimension
    Box<DIM, DataType> get_side(const int &i, const bool &low = true);

    /// return the box with the high and low value of the specified component
    /// set equal to the input index
    /// \param i   index corresponding to the requested dimension
    /// \param val index to in the dimension to take the slice at
    Box<DIM, DataType> get_slice(const int &i, const DataType &val);

    /// Calculate the box that is the intersection of this box and box B
    Box<DIM, DataType> intersection(const Box<DIM, DataType> &B);

    /// get the logical index of a physical point
    inline std::size_t get_lindex(const Point<DIM, DataType> &P) const;
    inline std::size_t get_lindex(const DataType P[DIM]) const;

    /// routine to return the logical data stride produced by adding the 
    /// provided vector
    ///
    /// \param[in] P vector to add to an index
    ///
    /// return shift in logical indexing
    inline std::size_t get_lstride(const Point<DIM, DataType> &P) const;
   
    /// get the physical index from a logical index
    inline Point<DIM, DataType> get_pindex(const std::size_t &i) const;

    /// get the memory usage in bytes
    std::size_t get_mem_usage() const;

  private:

    bool empty;
    Point<DIM, DataType> low_corner;
    Point<DIM, DataType> high_corner;
    Point<DIM, DataType> order;
    Point<DIM-1, DataType> strides;

};

template <int DIM, class DataType = long int>
class BoxIterator
{

  public:

    /// contructor --- defines an empty iterator
    BoxIterator();

    /// constructor --- associates the iterator to a specific box
    /// uses the default ordering (i + j*nx + k*nx*ny) where i, j, and k are
    /// the indices for the x, y, and z coordinates, respectively
    BoxIterator(const Box<DIM, DataType> &B);

    /// routine to return the data strides
    Point<DIM-1, DataType> get_strides() const {return strides;};

    /// start the iterations from the beginning
    inline void begin();

    /// start the iterations from the end
    inline void end();

    /// start the iterations from the provided point
    inline void begin(const Point<DIM, DataType> &P);

    /// function to check to see if the iterations have traversed the points
    /// in the box
    inline bool notDone() const {return !(done);};

    /// operator to increment the iteration forward
    inline void operator++();

    /// operator to increment the iteration backward
    inline void operator--();

    /// routine to return the current point
    inline Point<DIM, DataType> operator()() const {return pt;};

  private:

    bool done;
    Box<DIM, DataType> box;
    Point<DIM, DataType> pt;
    Point<DIM, DataType> ordering;
    Point<DIM-1, DataType> strides;

};

//
// definitions for Box class
//
template <int DIM, class DataType>
Box<DIM,DataType>::Box(const Point<DIM, DataType> &low, const Point<DIM, DataType> &high)
{
  bool sanity_check = true;
  
  // check the high point is above the low point
  for (int i=0; i<DIM; i++) {
    if (high[i] < low[i])
      sanity_check = false;
  }

  if (sanity_check) {

    low_corner = low;
    high_corner = high;
    empty = false;
    order = Point<DIM, DataType>();
    for (std::size_t i=0; i<DIM; i++)
      order.set_value(i,i);

    strides = Point<DIM-1, DataType>();
    for (int i=0; i<DIM-1; i++) {
      strides.set_value(i,high_corner[order[i]]-low_corner[order[i]]+1);
      if (i > 0)
        strides.set_value(i, strides[i]*strides[i-1]);
    }
  } else {
    empty = true;
  }
}

template <int DIM, class DataType>
Box<DIM,DataType>::Box(const Point<DIM, DataType> &low, 
                       const Point<DIM, DataType> &high,
                       const Point<DIM, DataType> &ordering)
{

  bool sanity_check = true;
  
  // check the high point is above the low point
  for (int i=0; i<DIM; i++) {
    if (high[i] < low[i])
      sanity_check = false;
  }

  if (sanity_check) {

    low_corner = low;
    high_corner = high;
    empty = false;
    order = ordering;

    strides = Point<DIM-1, DataType>();
    for (int i=0; i<DIM-1; i++) {
      strides.set_value(i,high_corner[order[i]]-low_corner[order[i]]+1);
      if (i > 0)
        strides.set_value(i, strides[i]*strides[i-1]);
    }
  } else {
    empty = true;
  }
}

template <int DIM, class DataType>
void Box<DIM,DataType>::set_ordering(const Point<DIM, DataType> &ordering)
{
  order = ordering;
  strides = Point<DIM-1, DataType>();
  for (int i=0; i<DIM-1; i++) {
    strides.set_value(i,high_corner[order[i]]-low_corner[order[i]]+1);
    if (i > 0)
      strides.set_value(i, strides[i]*strides[i-1]);
  }
}

template <int DIM, class DataType>
std::size_t Box<DIM,DataType>::get_npoints() const 
{
  std::size_t tmp;

  if (empty) {
    tmp = 0;
  } else {
    tmp = 1;
  }
  for (int i=0; i<DIM; i++) {
    tmp *= high_corner[i] - low_corner[i] + 1;
  }
  return tmp;
  
}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::grow(const DataType &i)
{
  if (empty)
    return Box<DIM, DataType>();

  Point<DIM, DataType> temp = Point<DIM, DataType>::scaled_ones(i);

  Point<DIM, DataType> low = low_corner;
  Point<DIM, DataType> high = high_corner;

  // subtract from lower corner
  low.subtract(temp);

  // add to high corner
  high.add(temp);

  // make sure the upper and lower bounds didn't cross
  for (int i=0; i<DIM; i++) {
    if (low[i] > high[i])
      return Box<DIM, DataType>(); // empty box;
  }

  return Box<DIM, DataType>(low, high, order);
}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::grow_high(const DataType &i)
{
  if (empty)
    return Box<DIM, DataType>();

  Point<DIM, DataType> temp = Point<DIM, DataType>::scaled_ones(i);

  Point<DIM, DataType> high = high_corner;

  // add to high corner
  high.add(temp);

  return Box<DIM, DataType>(low_corner, high, order);
}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::grow_low(const DataType &i)
{
  if (empty)
    return Box<DIM, DataType>();

  Point<DIM, DataType> temp = Point<DIM, DataType>::scaled_ones(i);

  Point<DIM, DataType> low = low_corner;

  // subtract from lower corner
  low.subtract(temp);

  return Box<DIM, DataType>(low, high_corner, order);
}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::add_high(const Point<DIM, DataType> &P)
{
  if (empty)
    return Box<DIM, DataType>();

  Point<DIM, DataType> temp = P;

  // add to high corner
  temp.add(high_corner);

  return Box<DIM, DataType>(low_corner, temp, order);
}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::add_low(const Point<DIM, DataType> &P)
{
  if (empty)
    return Box<DIM, DataType>();

  Point<DIM, DataType> temp = P;
  temp.add(low_corner);

  return Box<DIM, DataType>(temp, high_corner, order);
}

template <int DIM, class DataType>
void Box<DIM,DataType>::partition(const Point<DIM, DataType> &n, Box<DIM, DataType> *part_boxes)
{
  std::size_t end, i;
  int j, k;

  Point<DIM, DataType> low = low_corner;
  Point<DIM, DataType> high = low_corner;
  Point<DIM, DataType> diffs = Point<DIM, DataType>();
  Point<DIM, DataType> counts = Point<DIM, DataType>();

  for (i=0; i<DIM; i++) {
    diffs.set_value(i, (high_corner[i] - low_corner[i]) / n[i]);
  }

  end = 1;
  for (i=0; i<DIM; i++)
    end *= n[i];

  for (i=0; i < end; i++) {

    for (j=0; j<DIM; j++) {
      low.set_value(j, low_corner[j] + counts[j]*diffs[j]);
      if (counts[j] == n[j]-1) {
        high.set_value(j, high_corner[j]);
      } else {
        high.set_value(j, low_corner[j] + (counts[j]+1)*diffs[j]);
      }
    }

    for (j=0; j<DIM; j++) {
      if (counts[j] < n[j]-1) {
        counts.add_to_component(j,1);
        for (k=0; k<j; k++)
          counts.set_value(k, 0);
        break;
      }
    }

    part_boxes[i] = Box<DIM,DataType>(low, high, order);

  }

}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::get_side(const int &i, const bool &low)
{
  if (empty)
    return Box<DIM, DataType>();

  Point<DIM, DataType> high = high_corner;
  Point<DIM, DataType> lo = low_corner;

  if (low) {
    high.set_value(i, low_corner[i]);
  } else {
    lo.set_value(i, high_corner[i]);
  }

  return Box<DIM,DataType>(lo, high, order);
}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::get_slice(const int &i, const DataType &val)
{
  if (empty)
    return Box<DIM, DataType>();

  // make sure the slice location is valid
  if ((val > high_corner[i]) || (val < low_corner[i]))
    return Box<DIM, DataType>();

  Point<DIM, DataType> high = high_corner;
  Point<DIM, DataType> lo = low_corner;

  high.set_value(i, val);
  lo.set_value(i, val);

  return Box<DIM,DataType>(lo, high, order);
}

template <int DIM, class DataType>
Box<DIM, DataType> Box<DIM,DataType>::intersection(const Box<DIM, DataType> &B)
{
  Point<DIM, DataType> high = Point<DIM, DataType>();
  Point<DIM, DataType> low = Point<DIM, DataType>();

  DataType temp;

  Box<DIM, DataType> output;

  if ((B.isEmpty()) || (empty))
    return output = Box<DIM, DataType>();

  for (int i = 0; i<DIM; i++) {
  
    if ((B.get_high()[i] >= low_corner[i]) && (B.get_high()[i] <= high_corner[i])) {
      temp = B.get_high()[i];
    } else if ((high_corner[i] >= B.get_low()[i]) && (high_corner[i] <= B.get_high()[i])) {
      temp = high_corner[i];
    } else { // no intersection
      return output = Box<DIM, DataType>();
    }

    high.set_value(i, temp);

    if ((B.get_low()[i] >= low_corner[i]) && (B.get_low()[i] <= high_corner[i])) {
      temp = B.get_low()[i];
    } else if ((low_corner[i] >= B.get_low()[i]) && (low_corner[i] <= B.get_high()[i])) {
      temp = low_corner[i];
    } else { // no intersection
      return output = Box<DIM, DataType>();
    }

    low.set_value(i, temp);

  } 

  return output = Box<DIM, DataType>(low, high, order);

}

template <int DIM, class DataType>
std::size_t Box<DIM,DataType>::get_lindex(const Point<DIM, DataType> &P) const
{

  std::size_t index;
  index = 0;
  for (int i=DIM-1; i>0; i--)
    index += (P[order[i]] - low_corner[order[i]])*strides[i-1];
  index += P[order[0]] - low_corner[order[0]];

  return index;
}

template <int DIM, class DataType>
std::size_t Box<DIM,DataType>::get_lindex(const DataType P[DIM]) const
{

  std::size_t index;
  index = 0;
  for (int i=DIM-1; i>0; i--)
    index += (P[order[i]] - low_corner[order[i]])*strides[i-1];
  index += P[order[0]] - low_corner[order[0]];

  return index;
}

template <int DIM, class DataType>
std::size_t Box<DIM,DataType>::get_lstride(const Point<DIM, DataType> &P) const
{
  std::size_t index;
  index = 0;
  for (int i=DIM-1; i>0; i--)
    index += P[order[i]]*strides[i-1];
  index += P[order[0]];

  return index;
}

template <int DIM, class DataType>
Point<DIM, DataType> Box<DIM,DataType>::get_pindex(const std::size_t &i) const
{
  std::size_t temp = i;
  Point<DIM, DataType> index = Point<DIM, DataType>();
  for (int i=DIM-1; i>0; i--) {
    index.set_value(order[i], low_corner[order[i]] + temp / strides[i-1]);
    temp %= strides[i-1];
  }
  index.set_value(order[0], low_corner[order[0]] + temp);
  return index;
}

template <int DIM, class DataType>
std::size_t Box<DIM,DataType>::get_mem_usage() const
{
  return low_corner.get_mem_usage() + high_corner.get_mem_usage() + order.get_mem_usage() + \
         strides.get_mem_usage() + sizeof(bool);

}

//------------------------------------------------------------------------------
//
//                Definitions for the box iterator class
//
//------------------------------------------------------------------------------

template <int DIM, class DataType>
BoxIterator<DIM,DataType>::BoxIterator(const Box<DIM, DataType> &B)
{
  done = true;
  box = B;
  pt = Point<DIM, DataType>();

  ordering = Point<DIM, DataType>();
  
  for (int i=0; i<DIM; i++)
    ordering.set_value(i,i);
}

template <int DIM, class DataType>
BoxIterator<DIM,DataType>::BoxIterator()
{
  done = true;
  box = Box<DIM, DataType>();
  pt = Point<DIM, DataType>();

  ordering = Point<DIM, DataType>();
  
  for (int i=0; i<DIM; i++)
    ordering.set_value(i,i);
}

template <int DIM, class DataType>
void BoxIterator<DIM,DataType>::begin()
{
  done = (box.isEmpty()) ? true : false;
  pt = box.get_low();
  ordering = box.get_ordering();
  strides = box.get_strides();
}

template <int DIM, class DataType>
void BoxIterator<DIM,DataType>::end()
{
  done = (box.isEmpty()) ? true : false;
  pt = box.get_high();
  ordering = box.get_ordering();
  strides = box.get_strides();
}

template <int DIM, class DataType>
void BoxIterator<DIM,DataType>::begin(const Point<DIM, DataType> &P)
{
  done = (box.isEmpty()) ? true : false;
  pt = P;
  ordering = box.get_ordering();
  strides = box.get_strides();
}

template <int DIM, class DataType>
void BoxIterator<DIM,DataType>::operator++()
{

  int j, k;

  for (j=0; j<DIM; j++) {
    if (pt[ordering[j]] < box.get_high()[ordering[j]]) {
      pt.add_to_component(ordering[j],1);
      for (k=0; k<j; k++)
        pt.set_value(ordering[k], box.get_low()[ordering[k]]);
      break;
    } else if (j == DIM-1) {
      done = true;
    }
  }

}

template <int DIM, class DataType>
void BoxIterator<DIM,DataType>::operator--()
{

  int j, k;

  for (j=0; j<DIM; j++) {
    if (pt[ordering[j]] > box.get_low()[ordering[j]]) {
      pt.add_to_component(ordering[j],-1);
      for (k=0; k<j; k++)
        pt.set_value(ordering[k], box.get_high()[ordering[k]]);
      break;
    } else if (j == DIM-1) {
      done = true;
    }
  }

}

} // end namespace

#endif // GRID_BOX_H_
