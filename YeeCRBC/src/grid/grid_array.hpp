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
  Header file to define an multidimensional data class

*/

#include "point.hpp"
#include "box.hpp"
#include <vector>
#include <iostream>

#ifndef GRID_ARRAY_H_
#define GRID_ARRAY_H_
namespace grid {

template <int DIM, class DataType = double, class IndexType = long int>
class Array
{

  public:

    /// default constructor --- does nothing, if used, the object will have to 
    /// initialized later using, e.g, the assignment operator.
    Array() {};

    /// constructor --- defines an array with to store values associated with 
    /// the points in the provided box. Initializes values to 0.
    Array(const Box<DIM, IndexType> &domain);

    /// constructor --- defines an array with to store values associated with 
    /// the points in the provided box and points the data array to an existing
    /// array
    Array(const Box<DIM, IndexType> &domain, DataType *data_array);

    /// copy constructor --- this should be used with care becuase the memory 
    /// footprint can be substantial
    Array(const Array<DIM, DataType, IndexType> &B);

    /// Assignment operator --- this should be used with care becuase the memory 
    /// footprint can be substantial
    Array<DIM, DataType, IndexType> & operator=(const Array<DIM, DataType, IndexType> &B);

    /// swap function --- this swaps two Arrays without copying the data
    void swap(Array<DIM, DataType, IndexType> &B);

    /// destructor
    virtual ~Array() {};

    /// get the indexing box
    Box<DIM, IndexType> get_box() {return box;};

    // get number of data points
    std::size_t get_n() const {return n;};

    /// define a method for accessing data using the "physical" point values in
    /// the defining box
    inline DataType & operator() (const Point<DIM, IndexType> &P);
    inline DataType operator() (const Point<DIM, IndexType> &P) const;
    inline DataType & operator() (const IndexType P[DIM]);
    inline DataType operator() (const IndexType P[DIM]) const;

    /// define a method for accessing data using the "logical" index.
    inline DataType & operator[] (const std::size_t &i);
    inline DataType operator[] (const std::size_t &i) const;

    /// get a pointer to the data
    DataType* get_data() {return data.data();};
    const DataType* get_data() const {return data.data();};

    /// get memory size in bytes
    std::size_t get_mem_usage() const;

    /// method to zero out the data
    void zero();

  private:

    Box<DIM, IndexType> box;
    std::vector<DataType> data;
    std::size_t n;
    Point<DIM, IndexType> order;

};

//
// definitions
//

template <int DIM, class DataType, class IndexType>
Array<DIM, DataType, IndexType>::Array(const Box<DIM, IndexType> &domain)
{

  order = Point<DIM, IndexType>();
  box = domain;
  n = box.get_npoints();
  for (int i=0; i<DIM; i++) {
    order.set_value(i,i);
  }

  data.assign(n, static_cast<DataType>(0));
  order = box.get_ordering();
}

template <int DIM, class DataType, class IndexType>
Array<DIM, DataType, IndexType>::Array(const Box<DIM, IndexType> &domain, DataType *data_array)
{

  order = Point<DIM, IndexType>();
  box = domain;
  n = box.get_npoints();
  for (int i=0; i<DIM; i++) {
    order.set_value(i,i);
  }

  data.reserve(n);
  for (size_t i =0; i<n; i++)
    data.push_back(data_array[i]);

  order = box.get_ordering();
}


template <int DIM, class DataType, class IndexType>
Array<DIM, DataType, IndexType>::Array(const Array<DIM, DataType, IndexType> &B)
{
  order = B.order;
  box = B.box;
  n = B.n;
  
  data = B.data;
}

template <int DIM, class DataType, class IndexType>
Array<DIM, DataType, IndexType> & Array<DIM, DataType, IndexType>::operator=(const Array<DIM, DataType, IndexType> &B)
{
  Box<DIM, IndexType> tmp_box = B.box;
  std::vector<DataType> tmp_data;
  std::size_t tmp_n = B.n;
  Point<DIM, IndexType> tmp_order = B.order;

  tmp_data = B.data;

  data = tmp_data;
  n = tmp_n;
  box = tmp_box;
  order = tmp_order;
  
  return *this;
}

template <int DIM, class DataType, class IndexType>
void Array<DIM, DataType, IndexType>::swap(Array<DIM, DataType, IndexType> &B)
{
  std::swap(order, B.order);
  std::swap(box, B.box);
  std::swap(n, B.n);
  std::swap(data, B.data);

}

template <int DIM, class DataType, class IndexType>
DataType & Array<DIM, DataType, IndexType>::operator()(const Point<DIM, IndexType> &P)
{
  return data[box.get_lindex(P)];
}

template <int DIM, class DataType, class IndexType>
DataType Array<DIM, DataType, IndexType>::operator()(const Point<DIM, IndexType> &P) const
{
  return data[box.get_lindex(P)];
}

template <int DIM, class DataType, class IndexType>
DataType & Array<DIM, DataType, IndexType>::operator()(const IndexType P[DIM])
{
  return data[box.get_lindex(P)];
}

template <int DIM, class DataType, class IndexType>
DataType Array<DIM, DataType, IndexType>::operator()(const IndexType P[DIM]) const
{
  return data[box.get_lindex(P)];
}

template <int DIM, class DataType, class IndexType>
DataType & Array<DIM, DataType, IndexType>::operator[](const std::size_t &i)
{
  return data[i];
}

template <int DIM, class DataType, class IndexType>
DataType  Array<DIM, DataType, IndexType>::operator[](const std::size_t &i) const
{
  return data[i];
}

template <int DIM, class DataType, class IndexType>
std::size_t Array<DIM, DataType, IndexType>::get_mem_usage() const
{
  std::size_t mem_use;
 
  mem_use = sizeof(DataType)*n + sizeof(std::size_t);
  mem_use += box.get_mem_usage();
  mem_use += order.get_mem_usage();

  return mem_use;
}

template <int DIM, class DataType, class IndexType>
void Array<DIM, DataType, IndexType>::zero() 
{
  data.assign(n, static_cast<DataType>(0));
}

} // end namespace

#endif // GRID_ARRAY_H_
