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

/* C++ interface to the DAB/CRBC updates for the 3D Yee scheme

*/

#ifndef THREE_D_YEE_API_HPP_
#define THREE_D_YEE_API_HPP_

#if USE_CRBC_SCALAR_FLOAT
  typedef float ScalarType;
#else
  typedef double ScalarType;
#endif

#if USE_CRBC_COEF_FLOAT
  typedef float CoefType;
#else
  typedef double CoefType;
#endif

#if USE_CRBC_INDEX_LONG
  typedef long int IndexType;
#elif USE_CRBC_INDEX_LONG_LONG
  typedef long long IndexType;
#else
  typedef int IndexType;
#endif


// c++ interface
#include "crbc_updates.hpp"
#include "boundary_properties.hpp"

namespace yee_crbc
{

/// define enumerations for the valid boundary types
enum Boundaries
{
  PEC,  /**< Perfect Electric Conductor */
  // PMC,  /**< Perfect Magnetic Conductor */
  CRBC /**< Complete Radiation Boundary Condition (Implemented as a DAB) */
};

/// define the valid Field Component Types
enum Fields
{
  Ex, /**< x-component of the Electric Field */
  Ey, /**< y-component of the Electric Field */
  Ez /*,*/ /**< z-component of the Electric Field */
//  Hx, /**< x-component of the Magnetic Field */
//  Hy, /**< y-component of the Magnetic Field */
//  Hz /**< z-component of the Magnetic Field */
};

/// define valid sides
enum Side {                                                            
  XLeft,  //*< Left side in the x-coordinate direction */
  XRight, //*< Right side in the x-coordinate direction */
  YLeft,  //*< Left side in the y-coordinate direction */
  YRight, //*< Right side in the y-coordinate direction */
  ZLeft,  //*< Left side in the z-coordinate direction */
  ZRight  //*< Right side in the z-coordinate direction */
};                                                                 


/// Class to handle the DAB/CRBC updates for the Yee scheme
class YeeCrbcUpdates3d 
{

  public:

    YeeCrbcUpdates3d() {initialized = false;};

    YeeCrbcUpdates3d (const CoefType &T,
                      const CoefType h[3],  
                      const CoefType &dt,                  
                      const CoefType &c,
                      const Boundaries boundaries[6]);

    YeeCrbcUpdates3d (const CoefType &T,
                      const CoefType h[3],  
                      const CoefType &dt,                  
                      const CoefType &c,
                      const Boundaries boundaries[6],
                      const int &P);

    YeeCrbcUpdates3d (const CoefType &T,
                      const CoefType h[3],  
                      const CoefType &dt,                  
                      const CoefType &c,
                      const Boundaries boundaries[6],
                      const int &Pmax,
                      const double &tol);

    bool isEmpty() const {return !initialized;};

    virtual ~YeeCrbcUpdates3d() {};

    int init_face(const Side &side,
                  const Fields &comp, 
                  const IndexType low_index[3],
                  const IndexType high_index[3],
                  const double &delta);

    void load_face_data(const Side &side,
                               const Fields &comp,
                               const IndexType P[3],
                               const ScalarType &val);

    ScalarType get_new_face_vals(const Side &side,
                                        const Fields &comp,
                                        const IndexType P[3]) const;

    int compute_updates();

    double get_reflection_coef(const Side &side) const;

    double get_max_reflection_coef() const;

    int get_num_recursions(const Side &side) const;

    void get_component_inputs(const Side &side, Fields *comp, int &n) const;

    void get_input_extents(const Side &side, const Fields &comp, IndexType low[3], IndexType high[3]) const;

    void get_output_extents(const Side &side, const Fields &comp, IndexType low[3], IndexType high[3]) const;

    CoefType get_c() const;

    CoefType get_T() const;

    CoefType get_dt() const; 

    int get_num_faces(const Fields &comp) const;

    int get_num_edges(const Fields &comp) const;
 
    int get_num_corners(const Fields &comp) const;


    #if USE_HDF5_RESTARTS
      int restart (const std::string &fname);
      int save_state (const std::string &fname) const;
    #endif

  private:
   
    crbc::CrbcUpdates<3, ScalarType, IndexType> updater[3];  
    bool initialized;
    bool require_component[3];
    Boundaries input_boundaries[6];
    crbc::BoundaryProperties::Boundary internal_bounds[3][6];
    int update_index;
    bool useTol;
    double tol;
    int Pmax;

    void compute_internal_bounds();
    

}; // end class

} // end namespace

#endif // 3D_YEE_API_H
