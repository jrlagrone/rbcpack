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
  Header file to define a point class. The primary use here will be to store 
  index values.

*/
#include <cstddef>  // for std::size_t

#ifndef GRID_POINT_H_
#define GRID_POINT_H_

namespace grid {

/// Class to define a dimension independent point
/// \tparam DIM       Number of dimensions
/// \tparam DataType  Data type to use, the intention is that this is a float or
///                  integer type, but this is not checked. Using another type
///                  such as a char will likely break the implementation.
template<int DIM, class DataType = long int>
class Point
{

  public:

    /// Constructor --- create a point initialized to zeros (or whatever 
    /// static_cast<DataType>(0) gives --- this may be undefined in which
    /// case this shouldn't be expected to work as expected)
    Point();

    /// Constructor --- create a point at the given coordinates
    Point(const DataType P[DIM]);

    /// assignment operator this = B
    inline Point<DIM, DataType> & operator = (const Point<DIM, DataType> &B);

    /// [] data access -- const protected
    inline DataType operator[](const int &i) const;

    /// get a pointer to the raw data
    DataType* data() {return &(pt[0]);};
    const DataType* data() const {return &(pt[0]);};

    /// == operator --- check to see if two points are the same
    inline bool operator==(const Point<DIM, DataType> &B) const;

    /// != operator --- check to see if two points are different
    inline bool operator!=(const Point<DIM, DataType> &B) const;

    /// Add points return new point = this + B
    inline Point<DIM, DataType> operator+(const Point<DIM, DataType> &B) const;

    /// subtract points return new point = this - B
    inline Point<DIM, DataType> operator-(const Point<DIM, DataType> &B) const;

    /// Add points this = this + B
    inline void add(const Point<DIM, DataType> &B);

    /// Subtract points this = this - B
    inline void subtract(const Point<DIM, DataType> &B);

    /// multiple a point by a scalar: this = b*this
    inline void scale(const DataType &b);

    /// add b to the jth component return = this(j) + b (zero based indexing)
    inline void add_to_component(const int &j, const DataType &b);  

    /// Set the value of a point
    inline void set_value(const DataType *P);

    /// Set the value of the jth component of a point to a
    inline void set_value(const int &j, const DataType &a);

    /// get the value of a point
    inline void get_value(DataType *P) const;

    /// return the jth standard unit vector e_j
    inline static Point<DIM, DataType> unit(const int &j);

    /// return the jth standard unit vector scaled by a: return = a*e_j
    inline static Point<DIM, DataType> scaled_unit(const int &j, const DataType &a);

    /// return a vector where every component has value 1
    inline static Point<DIM, DataType> ones();

    /// return a vector where every component has value a
    inline static Point<DIM, DataType> scaled_ones(const DataType &a);

    /// return the memory usage in bytes
    std::size_t get_mem_usage() const;

    
  private:

    DataType pt[DIM];
    
};

  // add the definitions since everything is templated
  template<int DIM, class DataType> 
  Point<DIM, DataType>::Point() {

    for (int i = 0; i< DIM; i++)
      pt[i] = static_cast<DataType>(0);

  }

  template<int DIM, class DataType> 
  Point<DIM, DataType>::Point(const DataType P[DIM]) {

    for (int i = 0; i< DIM; i++)
      pt[i] = P[i];

  }

  template<int DIM, class DataType> 
  Point<DIM, DataType> & 
  Point<DIM, DataType>::operator = (const Point<DIM, DataType> &B) {
      
    for (int i=0; i<DIM; i++)
      pt[i] = B.pt[i];

    return *this;

  }

  template<int DIM, class DataType>
  DataType Point<DIM, DataType>::operator[](const int &i) const
  {
    return pt[i];
  }

  template<int DIM, class DataType>
  inline bool Point<DIM, DataType>::operator==(const Point<DIM, DataType> &B) const
  {
    for (int i=0; i<DIM; i++) {
      if (pt[i] != B[i])
        return false;
    }
    return true;
  }

  template<int DIM, class DataType>
  inline bool Point<DIM, DataType>::operator!=(const Point<DIM, DataType> &B) const
  {
  
    int test = 0;
   
    for (int i=0; i<DIM; i++) {
      if (pt[i] != B[i])
        return true;
    }
    return false;
  }

  template<int DIM, class DataType>
  Point<DIM, DataType> Point<DIM, DataType>::operator+(const Point<DIM, DataType> &B) const
  {
    Point<DIM, DataType> temp = Point<DIM, DataType>(pt);
    for (int i = 0; i< DIM; i++)
      temp.add_to_component(i,B[i]);

    return temp;
  }

  template<int DIM, class DataType>
  Point<DIM, DataType> Point<DIM, DataType>::operator-(const Point<DIM, DataType> &B) const
  {
    Point<DIM, DataType> temp = Point<DIM, DataType>(pt);
    for (int i = 0; i< DIM; i++)
      temp.add_to_component(i,-B[i]);

    return temp;
  }

  template<int DIM, class DataType>
  void Point<DIM, DataType>::add(const Point<DIM, DataType> &B)
  {
    for (int i = 0; i< DIM; i++)
      pt[i] += B[i];
  }

  template<int DIM, class DataType>
  void Point<DIM, DataType>::subtract(const Point<DIM, DataType> &B)
  {
    for (int i = 0; i< DIM; i++)
      pt[i] -= B[i];
  }

  template<int DIM, class DataType>
  void Point<DIM, DataType>::scale(const DataType &b)
  {
    for (int i = 0; i< DIM; i++)
      pt[i] *= b;
  }

  template<int DIM, class DataType>
  void Point<DIM, DataType>::add_to_component(const int &j, const DataType &b)
  {
    pt[j] += b;
  }  

  template<int DIM, class DataType> 
  void Point<DIM, DataType>::set_value(const DataType *P) {

    for (int i = 0; i< DIM; i++)
      pt[i] = P[i];

  }

  template<int DIM, class DataType> 
  void Point<DIM, DataType>::set_value(const int &j, const DataType &a) {
      pt[j] = a;
  }

  template<int DIM, class DataType> 
  void Point<DIM, DataType>::get_value(DataType *P) const {

    for (int i = 0; i< DIM; i++)
      P[i] = pt[i];

  }

  template<int DIM, class DataType>
  Point<DIM, DataType> 
  Point<DIM, DataType>::unit(const int &j)
  {
    DataType temp[DIM];
    for (int i = 0; i< DIM; i++)
      temp[i] = 0;
    temp[j] = 1;
    return Point<DIM, DataType>(temp);
  }

  template<int DIM, class DataType>
  Point<DIM, DataType> 
  Point<DIM, DataType>::scaled_unit(const int &j, const DataType &a)
  {
    DataType temp[DIM];
    for (int i = 0; i< DIM; i++)
      temp[i] = 0;
    temp[j] = a;
    return Point<DIM, DataType>(temp);
  }

  template<int DIM, class DataType>
  Point<DIM, DataType> 
  Point<DIM, DataType>::ones()
  {
    DataType temp[DIM];
    for (int i = 0; i< DIM; i++)
      temp[i] = 1;
    return Point<DIM, DataType>(temp);
  }

  template<int DIM, class DataType>
  Point<DIM, DataType> 
  Point<DIM, DataType>::scaled_ones(const DataType &a)
  {
    DataType temp[DIM];
    for (int i = 0; i< DIM; i++)
      temp[i] = a;
    return Point<DIM, DataType>(temp);
  }

  template<int DIM, class DataType>
  std::size_t Point<DIM, DataType>::get_mem_usage() const
  {
    return sizeof(DataType)*DIM;
  }


} // end namespace

#endif // GRID_POINT_H_
