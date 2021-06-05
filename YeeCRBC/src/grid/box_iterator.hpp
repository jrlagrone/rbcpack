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
#ifndef GRID_BOX_ITERATOR_H_
#define GRID_BOX_ITERATOR_H_

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <cstddef>   // for std::size_t
#include <iostream>  // cout, etc.
#include "point.hpp"
#include "box.hpp"

namespace grid {


template <int DIM, class DataType = long int>
class BoxIterator : public std::iterator<std::random_access_iterator_tag, DataType>
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

    std::size_t m_i;

    bool done;
    Box<DIM, DataType> box;
    Point<DIM, DataType> pt;
    Point<DIM, DataType> ordering;
    Point<DIM-1, DataType> strides;

};

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

#endif // GRID_BOX_ITERATOR_H_
