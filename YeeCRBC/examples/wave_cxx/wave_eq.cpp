/**
    3-D FDTD code with CRBC/DAB boundary conditions

    This C++ code implements the finite difference time-domain solution of the 
    wave equation using standard second order, centered diffences. The grid is 
    terminated using Double Absorbing Boundaries (DAB).

    This is meant to serve as an example of how one might use the underlying
    C++ interface to the CRBC/DAB library. This interface is templated and 
    supports up to three dimensions. It is also templated on the internal
    indexing and data types.
    
    This code is freely distributed with the goal of illustrating the CRBC/DAB
    library usage and research purposes. We provide no warranties of its 
    performance or suitability for any purpose.

*/ 

// define some macros for indexing
#define i2d(i,j) (i + (j)*imax[0])
#define i3d(i,j,k) (i + imax[0]*((j) + (k)*imax[1]))
#define PI std::atan(1.0)*4

// Includes

// This file contains the most generic C++ interface to the CRBC/DAB boundary
// conditions.
#include <crbc_updates.hpp>

// C++ strings
#include <string>

// So we can append ints to strings
#include <sstream>

// for file output
#include <fstream>

// Needed for C++ output
#include <iostream>

// We use this to get std::swap so we can swap pointers for arrays and vectors
// without having to copy data. Note that as of C++11, std::swap is suppose to 
// move to <utility>. We are currently trying to avoid using C++11 for
// portability reasons
#include <algorithm>

// This has the declarations of the 'sqrt' and 'fabs' functions
#include <cmath>

#if USE_OPENMP
  #include <omp.h>
#endif


// We begin by creating a class to handle the wave equation updates along with
// all of the data. We use templates so that the same code can be used to handle
// 1, 2, or 3 dimensions.
template <int DIM>
class wave_equation
{

  
  private:

    // -------------------------------------------------------------------------
    // first we'll define all the data storage that we will need
    // -------------------------------------------------------------------------

    // arrays to hold old, current, and new wave equation values used in the
    // time stepping
    std::vector<double> old_data, cur_data, new_data;

    double c;       // wave speed
    double dt;      // time step
    double h[DIM];  // grid spacings

    int ntsteps;   // number of time steps
    std::size_t n; // number of grid points
    int imax[DIM]; // indexing limits in each direction
    int out_freq;  // output frequency, in time steps
    int src[DIM];  // source location

    bool save_output; // flag indicating whether we should save the field vals.

    std::string base_name; // base file name for output

    // Boundary conditions. We will use the boundaries enumerations provided
    // by the BoundaryProperties class in the CRBC library.
    crbc::BoundaryProperties::Boundary boundaries[2*DIM];
    
    // a CRBC updater object to handle the boundaries
    // We declare this using the template parameters to indicate we our
    // calculations with be done in DIM dimensions, we will use doubles to store
    // the field values and use ints for indexing.
    crbc::CrbcUpdates<DIM, double, int> boundary_updater;

    // -------------------------------------------------------------------------
    // now define the private member functions
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // function to update the internal grid points
    // -------------------------------------------------------------------------
    void time_step() {

      int i, j, k;
      
      double b[DIM];

      // precompute coefficients
      for (i=0; i<DIM; ++i)
        b[i] = dt*dt*c*c/(h[i]*h[i]);

      // we'll explicitly differentiate each of the possible dimensions here for
      // simplicity, but not that it is possible to do this in a truly dimension
      // independent manner.
      //
      // loop over the internal grid points and apply the wave equation
      // update in the appropriate number of dimensions
      switch (DIM) {
      
        case 1:

          //   unew = 2ucur - uold + dt^2*c^2/h_{x}^2 *(u_{x_i-1} - 2u_{x_i} + u_{x_i+1})
          for (i=1; i<imax[0]-1; ++i) {
            new_data[i] = 2.0*cur_data[i] - old_data[i] \
              + b[0]*(cur_data[i-1] - 2.0*cur_data[i] + cur_data[i+1]);
          }
          break;

        case 2:

          //   unew = 2ucur - uold + dt^2*c^2/h_{x}^2 *(u_{x_i-1} - 2u_{x_i} + u_{x_i+1})
          //                       + dt^2*c^2/h_{y}^2 *(u_{y_i-1} - 2u_{y_i} + u_{y_i+1})
          for (j=1; j<imax[1]-1; ++j) {
            for (i=1; i<imax[0]-1; ++i) {

              new_data[i2d(i,j)] = 2.0*cur_data[i2d(i,j)] - old_data[i2d(i,j)] \
                + b[0]*(cur_data[i2d(i-1,j)] - 2.0*cur_data[i2d(i,j)]+ cur_data[i2d(i+1,j)]) \
                + b[1]*(cur_data[i2d(i,j-1)] - 2.0*cur_data[i2d(i,j)] + cur_data[i2d(i,j+1)]);
            }
          }
          break;

        case 3:

          //   unew = 2ucur - uold + dt^2*c^2/h_{x}^2 *(u_{x_i-1} - 2u_{x_i} + u_{x_i+1})
          //                       + dt^2*c^2/h_{y}^2 *(u_{y_i-1} - 2u_{y_i} + u_{y_i+1})
          //                       + dt^2*c^2/h_{z}^2 *(u_{z_i-1} - 2u_{z_i} + u_{z_i+1})
          for (k=1; k<imax[2]-1; ++k) {
            for (j=1; j<imax[1]-1; ++j) {
              for (i=1; i<imax[0]-1; ++i) {
                new_data[i3d(i,j,k)] = 2.0*cur_data[i3d(i,j,k)] - old_data[i3d(i,j,k)] \
                  + b[0]*(cur_data[i3d(i-1,j,k)] - 2.0*cur_data[i3d(i,j,k)] + cur_data[i3d(i+1,j,k)]) \
                  + b[1]*(cur_data[i3d(i,j-1,k)] - 2.0*cur_data[i3d(i,j,k)] + cur_data[i3d(i,j+1,k)]) \
                  + b[2]*(cur_data[i3d(i,j,k-1)] - 2.0*cur_data[i3d(i,j,k)] + cur_data[i3d(i,j,k+1)]);
              }
            }
          }
          break;

        default:
          // In practice one should probably specify the exceptions for better error
          // handling, but for a simple example this should be fine.
          std::cerr << " Unsupported dimension in the time_step() function" << std::endl;
          throw;

      }

    };  // end time_step()

    // -------------------------------------------------------------------------
    // function to apply the source
    // -------------------------------------------------------------------------
    void apply_driving(const int &tstep)
    {
    
      switch (DIM)
      {
 
        case 1: 
          new_data[src[0]] = sin(c*tstep*PI*dt);
          break;
        case 2:
          new_data[i2d(src[0], src[1])] += 5.*sin(c*tstep*PI*dt);
          break;
        case 3:
          new_data[i3d(src[0], src[1], src[2])] += 100.*sin(c*tstep*PI*dt);
          break;
        default:
          // simple error handling
          std::cerr << " Unsupported dimension in the apply_driving() function" << std::endl;
          throw;
      }

    }; // end apply_driving

    // -------------------------------------------------------------------------
    // function to apply boundary updates. For this, we will demonstrate 
    // using the CRBC library and homogeneous Dirichlet boundaries. Note the 
    // CRBC library also supports homogenous Neumann boundaries.
    // -------------------------------------------------------------------------
    void step_boundaries()
    {

      int i, j, k, l;
      int low_ind[DIM], high_ind[DIM], ind[DIM];

      // since we initializaed the fields to be 0, Dirichlet boundaries are
      // automatically enforced. So we need only apply the CRBC type boundaries.
      // If supporting Neumann boundaries, these should be applied where 
      // possible before applying the CRBC library updates
      
      // First we need to copy the new values into the CRBC updater object
      // loop over the boundary faces:
      for (l=0; l<2*DIM; ++l) {

        // check to see if this is a face that the CRBC updater is handling
        if (boundaries[l] == crbc::BoundaryProperties::CRBC) {

          // get the indices the updater object expects as input from this face.
          // Note that these values are inclusive
          boundary_updater.get_input_extents(l, low_ind, high_ind);

          // Copy the updated values into CRBC library object. Again, we note that
          // we explicitly handle each of the dimensions for simplicity.
          switch (DIM) {

            case 1:

              for (i=low_ind[0]; i<=high_ind[0]; ++i)
                boundary_updater.load_face_data(l, &i, new_data[i]);

              break;

            case 2:

              for (j=low_ind[1]; j<=high_ind[1]; ++j) {
                ind[1] = j;
                for (i=low_ind[0]; i<=high_ind[0]; ++i) {
                  ind[0] = i;
                  boundary_updater.load_face_data(l, ind, new_data[i2d(i,j)]);
               }
              }

              break;

            case 3:

              for (k=low_ind[2]; k<=high_ind[2]; ++k) {
                ind[2] = k;
                for (j=low_ind[1]; j<=high_ind[1]; ++j) {
                  ind[1] = j;
                  for (i=low_ind[0]; i<=high_ind[0]; ++i) {
                    ind[0] = i;
                    boundary_updater.load_face_data(l, ind, new_data[i3d(i,j,k)]);
                  }
                }
              }

              break;

            default:
              // simple error handling
              std::cerr << " Unsupported dimension in the step_boundaries() function" 
                        << std::endl;
              throw;


          }
        } //end if
      } // end loop over faces

      // Now we can tell the CRBC library boundary updater to compute the updates
      boundary_updater.compute_updates();

      // Finally, we need to copy the updated boundary values from the CRBC library
      // We again loop over all of the faces, but this time we request the output
      // data extents from the updater object and then copy the new values into 
      // the appropriate locations
      for (l=0; l<2*DIM; ++l) {

        // check to see if this is a face that the CRBC updater is handling
        if (boundaries[l] == crbc::BoundaryProperties::CRBC) {

          // get the indices the updater object expects as input from this face.
          // Note that these values are inclusive
          boundary_updater.get_output_extents(l, low_ind, high_ind);

          // Copy the updated values from the crbc updater object.
          switch (DIM) {

            case 1:

              for (i=low_ind[0]; i<=high_ind[0]; ++i)
                new_data[i] = boundary_updater.get_new_face_vals(l, &i);

              break;

            case 2:

              for (j=low_ind[1]; j<=high_ind[1]; ++j) {
                ind[1] = j;
                for (i=low_ind[0]; i<=high_ind[0]; ++i) {
                  ind[0] = i;
                  new_data[i2d(i,j)] = boundary_updater.get_new_face_vals(l, ind);
               }
              }

              break;

            case 3:

              for (k=low_ind[2]; k<=high_ind[2]; ++k) {
                ind[2] = k;
                for (j=low_ind[1]; j<=high_ind[1]; ++j) {
                  ind[1] = j;
                  for (i=low_ind[0]; i<=high_ind[0]; ++i) {
                    ind[0] = i;
                    new_data[i3d(i,j,k)] = boundary_updater.get_new_face_vals(l, ind);
                  }
                }
              }

              break;

            default:
              // simple error handling
              std::cerr << " Unsupported dimension in the step_boundaries() function" 
                        << std::endl;
              throw;


          }
        } //end if
      } // end loop over faces

    }; // end step_boundaries()

    // -------------------------------------------------------------------------
    // Function to write output
    // -------------------------------------------------------------------------
    void write_output(std::string &fname) {

      int i;
      std::size_t j, cells;

      // calculate the number of cells
      cells = 1;
      for (i=0; i<DIM; i++)
        cells *= (imax[i]-1);

      // open output file
      std::ofstream outfile;
      outfile.open(fname.c_str());

      // write out the basic VTK header info
      outfile << "# vtk DataFile Version 3.0" << std::endl;
      outfile << "vtk output" << std::endl;
      outfile << "ASCII" << std::endl;
      outfile << "DATASET RECTILINEAR_GRID" << std::endl;

      // set the dimensions. note that this needs to be in 3D regardless of the
      // actual dimension.
      outfile << "DIMENSIONS " << imax[0] << " ";
      if (DIM > 1) {
        outfile << imax[1] << " ";
        if (DIM > 2) {
          outfile << imax[2] << std::endl;
        } else {
          outfile << 1 << std::endl;
        }
      } else {
        outfile << "1 1" << std::endl;
      }

      // save the coordinates
      outfile << "X_COORDINATES " << imax[0] << " float" << std::endl;
      for (i=0; i<imax[0]; i++)
        outfile << i*h[0] << std::endl;

      outfile << "Y_COORDINATES ";
      if (DIM > 1) {
        outfile << imax[1] << " float" << std::endl;
        for (i=0; i<imax[1]; i++)
          outfile << i*h[1] << std::endl;
      } else {
        outfile << 1 << " float" << std::endl;
        outfile << 0 << std::endl;
      }

      outfile << "Z_COORDINATES ";
      if (DIM > 2) {
        outfile << imax[2] << " float" << std::endl;
        for (i=0; i<imax[2]; i++)
          outfile << i*h[2] << std::endl;
      } else {
        outfile << 1 << " float" << std::endl;
        outfile << 0 << std::endl;
      }

      // set up a cell and field
      outfile << "CELL_DATA " << cells << std::endl;
      outfile << "POINT_DATA " << n << std::endl;
      outfile << "FIELD FieldData 1" << std::endl;
      outfile << "wave 1 " << n << " double" << std::endl;

      // now actually write the data, just a new line after every entry because
      // some systems have a character per line limit that we don't want to hit
      for (j=0; j<n; j++)
        outfile << new_data[j] << std::endl;

      // close file
      outfile.close();

    }; // end write output


  public:

    // -------------------------------------------------------------------------
    // constructor:
    //
    // Here, we'll specify the number of grid points in each direction and the
    // spacings as well as the wave speed, boundary conditions and the number of 
    // time steps. Finally the output frequency, file name base, and the
    // grid point to place the source are inputs.
    // -------------------------------------------------------------------------
    wave_equation (const int grid_points[DIM],
                   const double grid_spacing[DIM],
                   const double &c,
                   const crbc::BoundaryProperties::Boundary bounds[2*DIM],
                   const int ntsteps,
                   const std::string &base_file_name,
                   const int &out_freq,
                   const int src_location[DIM])
    {

      int i, l, low[DIM], high[DIM];
      double T, delta, tol;

      // first we'll save the inputs
      this->c = c;
      for (i=0; i<DIM; ++i) {
        imax[i] = grid_points[i];
        h[i] = grid_spacing[i];
        src[i] = src_location[i];
      }
      for (i=0; i<2*DIM; ++i)
        boundaries[i] = bounds[i];
      this->ntsteps = ntsteps;
      base_name = base_file_name;
      this->out_freq = out_freq;
      save_output = false;

      // calculate the time step size (0.99 * cfl)
      dt = 0.0;
      for (i=0; i<DIM; ++i)
        dt += 1.0/(h[i]*h[i]);
      dt = 0.99 / (c * sqrt(dt));

      // calculate the total simulation time
      T = dt * ntsteps;

      // calculate the total number of grid points
      n = 1;
      for (i=0; i<DIM; ++i)
        n *= imax[i];

      // initialize the data storage to be zero
      new_data.assign(n, 0.0);
      cur_data.assign(n, 0.0);
      old_data.assign(n, 0.0);

      // Now initialize the CRBC updater object (by assignment)
      boundary_updater = crbc::CrbcUpdates<DIM, double, int> (T, h, dt, c, boundaries);
      // Note we can change the default number of recursions from 5 to, say, 7 by
      // instead using the following constructor:
      // boundary_updater = crbc::CrbcUpdates<DIM, double, int> (T, h, dt, c, boundaries, 7);

      // now we need to set up the parameters for each of the faces. This can be 
      // done in 3 different ways. The main difference is how we choose the number
      // of recursions: we can use the defaults, we can specify the number of 
      // recursions, or we can specify a tolerance and the number of recursions
      // can be determined based on this tolerance.
      // We begin by looping over all of the possible sides:
      for (l=0; l<2*DIM; ++l) {

        // check to see if this is a face that the CRBC updater is handling.
        // NOTE that the sides are assumed to be in the following order
        //   left side in x  := 0
        //   right side in x := 1
        //   left side in y  := 2
        //   right side in y := 3
        //   left side in z  := 4
        //   right side in z := 5
        if (boundaries[l] == crbc::BoundaryProperties::CRBC) {

          // we need to calculate the minimum seperation between the boudary 
          // face and the source. This is just the distance from the source to
          // the boundary in the direction normal to the boundary.
          //
          // note if the boundary is on the left this distance is just the
          // appropriate component of the source location time the grid spacing
          // hence the (l%2)
          delta = std::abs(src[l/2] - (l%2)*imax[l/2]) * h[l/2];

          // handle the different dimensions explicitly
          switch (DIM) {

            case 1:

              // For the 1D case we will illustrate the initilizer where we
              // specify the number of recursions because using 0 recursions
              // corresponds to using the Sommerfeld radiation boundary conditions
              // which is exact in this case. The CRBC boundary updater needs to
              // know the index of the point on the boundary as well as the index
              // of the point immediately interior to the boundary. So if this is
              // the left boundary, the extents are [0,1]. For the right boundary,
              // the extents are [imax[0]-2, imax[0]-1] (-1 from 0 based indexing).
              // The updater expects these to be inclusive;
              if (l == 0 ) { // left side
                low[0] = 0;
                high[0] = 1;
              } else { // right side
                low[0] = imax[0]-2;
                high[0] = imax[0]-1;
              }
          
              boundary_updater.init_face(l, low, high, delta, 0);

              break;
 
            case 2:

              // For the 2D case we will illustrate the default initializer.
              // In 2D the updater needs to know the indexing extents for the
              // line of points on the boundary as well as the parallel line
              // of points immediately interior to the boundary.
              
              if (l == 0) { 
                // left boundary in x, need [0,1] in x, all in y
                low[0] = 0;
                low[1] = 0;
                high[0] = 1;
                high[1] = imax[1] - 1;
              } else if (l == 1) { 
                // right boundary in x, need [imax[0]-2, imax[0]-1] in x, all y
                low[0] = imax[0]-2;
                low[1] = 0;
                high[0] = imax[0] - 1;
                high[1] = imax[1] - 1;    
              } else if (l == 2) {               
                // left boundary in y, need [0,1] in y, all in x
                low[0] = 0;
                low[1] = 0;
                high[0] = imax[0] - 1;
                high[1] = 1;
              } else {
                // right boundary in y, need [imax[1]-2, imax[1]-1] in y, all x
                low[0] = 0;
                low[1] = imax[1]-2;
                high[0] = imax[0] - 1;
                high[1] = imax[1] - 1;   
              }

              // call initializer
              boundary_updater.init_face(l, low, high, delta);
              break;

            case 3:

              // For the 3D case, we will use the tolerance based initializer.
              // This tolerance controls the reflection coefficient of the 
              // boundary. In general, this usually provides a reasonable estimate
              // of the relative error due to the boundary. Some though should 
              // be given to it's choice, using too tight of a tolerance results
              // in more work being done for little or no accuracy benefit and
              // too loose of a tolerance results in the boundary error dominating.
              // A good choice is typicallyon the order of the expected
              // discretization error. Here we'll just choose 1e-2.
              tol = 1e-2;

              // In 3D the updater needs to know the indexing extents for the
              // plane of points on the boundary as well as the parallel plane
              // of points immediately interior to the boundary.
              
              if (l == 0) { 
                // left boundary in x, need [0,1] in x, all in y, z
                low[0] = 0;
                low[1] = 0;
                low[2] = 0;
                high[0] = 1;
                high[1] = imax[1] - 1;
                high[2] = imax[2] - 1;
              } else if (l == 1) { 
                // right boundary in x, need [imax[0]-2, imax[0]-1] in x, all y, z
                low[0] = imax[0]-2;
                low[1] = 0;
                low[2] = 0;
                high[0] = imax[0] - 1;
                high[1] = imax[1] - 1; 
                high[2] = imax[2] - 1;   
              } else if (l == 2) {               
                // left boundary in y, need [0,1] in y, all in x, z
                low[0] = 0;
                low[1] = 0;
                low[2] = 0;
                high[0] = imax[0] - 1;
                high[1] = 1;
                high[2] = imax[2] - 1;
              } else if (l == 3) {
                // right boundary in y, need [imax[1]-2, imax[1]-1] in y, all x, z
                low[0] = 0;
                low[1] = imax[1]-2;
                low[2] = 0;
                high[0] = imax[0] - 1;
                high[1] = imax[1] - 1; 
                high[2] = imax[2] - 1;  
              } else if (l == 4) {               
                // left boundary in z, need [0,1] in z, all in x, y
                low[0] = 0;
                low[1] = 0;
                low[2] = 0;
                high[0] = imax[0] - 1;
                high[1] = imax[1] - 1;
                high[2] = 1;
              } else {
                // right boundary in z, need [imax[2]-2, imax[2]-1] in z, all x, y
                low[0] = 0;
                low[1] = 0;
                low[2] = imax[2]-2;
                high[0] = imax[0] - 1;
                high[1] = imax[1] - 1; 
                high[2] = imax[2] - 1;  
              }

              // call initializer and limit the number of recursions to at most 20
              boundary_updater.init_face(l, low, high, delta, 20, tol);
              break;

            default:
              // simple error handling
              std::cerr << " Unsupported dimension in the constructor" << std::endl;
              throw;
          }
        }
      } // end loop over sides

      // Now we'll print out some information from the boundary updater.
      std::cout << "Recursion properties by face (reflection coef =-1 means that no updates are performed):" << std::endl;
      std::cout << "  Left side in x:" << std::endl;
      std::cout << "    recursions      = " << boundary_updater.get_num_recursions(0) << std::endl;
      std::cout << "    reflection coef = " << boundary_updater.get_reflection_coef(0) << std::endl;
      std::cout << "  Right side in x:" << std::endl;
      std::cout << "    recursions      = " << boundary_updater.get_num_recursions(1) << std::endl;
      std::cout << "    reflection coef = " << boundary_updater.get_reflection_coef(1) << std::endl;
      if (DIM > 1) {
        std::cout << "  Left side in y:" << std::endl;
        std::cout << "    recursions      = " << boundary_updater.get_num_recursions(2) << std::endl;
        std::cout << "    reflection coef = " << boundary_updater.get_reflection_coef(2) << std::endl;
        std::cout << "  Right side in y:" << std::endl;
        std::cout << "    recursions      = " << boundary_updater.get_num_recursions(3) << std::endl;
        std::cout << "    reflection coef = " << boundary_updater.get_reflection_coef(3) << std::endl;
      }
      if (DIM > 2) {
        std::cout << "  Left side in z:" << std::endl;
        std::cout << "    recursions      = " << boundary_updater.get_num_recursions(4) << std::endl;
        std::cout << "    reflection coef = " << boundary_updater.get_reflection_coef(4) << std::endl;
        std::cout << "  Right side in z:" << std::endl;
        std::cout << "    recursions      = " << boundary_updater.get_num_recursions(5) << std::endl;
        std::cout << "    reflection coef = " << boundary_updater.get_reflection_coef(5) << std::endl;
      }
      std::cout << "The maximum reflection coeficient is " << boundary_updater.get_max_reflection_coef()
                << std::endl;

      
      std::cout << "The CRBC library is updating" << std::endl;
      std::cout << "  " << boundary_updater.get_num_faces() << " faces" << std::endl;
      std::cout << "  " << boundary_updater.get_num_edges() << " edges" << std::endl;
      std::cout << "  " << boundary_updater.get_num_corners() << " corners" << std::endl;

      

    } // end constructor

    // -------------------------------------------------------------------------
    // function to turn the output file generation on/off
    // -------------------------------------------------------------------------
    void set_write_output(const bool &write_out) {save_output = write_out;};

    // -------------------------------------------------------------------------
    // function to run the simulation
    // -------------------------------------------------------------------------
    void run() {

      int t;
      std::ostringstream s;
      std::string fname;

      // time step:
      for (t=0; t<ntsteps; ++t) {

        // time step the interior
        time_step();
 
        // apply the driving term
        apply_driving(t);
   
        // update the boundaries
        step_boundaries();

        // generate output if needed
        if ((t % out_freq == 0) && (save_output)) {
 
          // strncpy (fname, base_name, sizeof(fname));
          s.str("");
          s.clear();
          s << base_name << "_" << t/out_freq << ".vtk";
          fname = s.str();
          write_output(fname);

        }

        // swap the storage vectors
        std::swap(cur_data, old_data);
        std::swap(cur_data, new_data);

      }

    }; //end run
    

}; // end wave_equation class

//------------------------------------------------------------------------------
//                           MAIN ROUTINE
//------------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  bool write_output_files = false;

  // read in input to see if we should write output files. If output files are
  // enabled this program writes out the 1D results 500 times which takes
  // approximately 7 MB of space. The 2D simulation generates 50 files, taking
  // up approximately 15 MB, and the 3D simulation generates 35 files taking up 
  // about 120 MB. 
  // By default, the file output is turned off.
  // There is only one option, so for simplicity we'll just assume that if we
  // recieve any command line option, then we should enable output instead of
  // actually parsing and identifying a specific option.
  if (argc > 1) {
    std::cout << "This program will generate output files." << std::endl;
    write_output_files = true;
  }

  std::cout << "1D simulation ..." << std::endl << std::endl;
  // first run the 1D wave equations
  int npoints_1d[] = {1000};
  double grid_spacing_1d[] = {0.01};
  double c = 1.0;
  int ntsteps = 5000;
  std::string bname = "1d_wave";
  int out_freq = 10;
  int src_location_1d[] = {350};
  crbc::BoundaryProperties::Boundary bounds_1d[2];

  bounds_1d[0] = crbc::BoundaryProperties::CRBC;
  bounds_1d[1] = crbc::BoundaryProperties::CRBC;

  // create a 1D simulation
  wave_equation<1> wave1d(npoints_1d, grid_spacing_1d, c, bounds_1d, ntsteps, bname, out_freq, src_location_1d);
  wave1d.set_write_output(write_output_files);

  // run
  try {
    wave1d.run();
  } catch (...) {
    std::cerr << "something with wrong ..." << std::endl;
  }

  std::cout << std::endl << std::endl  << "2D simulation ..." << std::endl << std::endl;
  // now do a 2d simulation
  int npoints_2d[] = {200, 200};
  double grid_spacing_2d[] = {0.01, 0.01};
  ntsteps = 500;
  bname = "2d_wave";
  int src_location_2d[] = {75, 120};
  crbc::BoundaryProperties::Boundary bounds_2d[4];

  bounds_2d[0] = crbc::BoundaryProperties::CRBC;
  bounds_2d[1] = crbc::BoundaryProperties::CRBC;
  bounds_2d[2] = crbc::BoundaryProperties::CRBC;
  bounds_2d[3] = crbc::BoundaryProperties::CRBC;

  // create a 2D simulation
  wave_equation<2> wave2d(npoints_2d, grid_spacing_2d, c, bounds_2d, ntsteps, bname, out_freq, src_location_2d);
  wave2d.set_write_output(write_output_files);

  // run
  try {
    wave2d.run();
  } catch (...) {
    std::cerr << "something with wrong ..." << std::endl;
  }

  std::cout << std::endl << std::endl << "3D simulation ..." << std::endl << std::endl;
  // now do a 3d simulation
  int npoints_3d[] = {75, 75, 75};
  double grid_spacing_3d[] = {0.01, 0.01, 0.01};
  ntsteps = 350;
  bname = "3d_wave";
  int src_location_3d[] = {30, 60, 40};
  crbc::BoundaryProperties::Boundary bounds_3d[6];

  bounds_3d[0] = crbc::BoundaryProperties::CRBC;
  bounds_3d[1] = crbc::BoundaryProperties::CRBC;
  bounds_3d[2] = crbc::BoundaryProperties::CRBC;
  bounds_3d[3] = crbc::BoundaryProperties::CRBC;
  bounds_3d[4] = crbc::BoundaryProperties::CRBC;
  bounds_3d[5] = crbc::BoundaryProperties::CRBC;

  // create a 3D simulation
  wave_equation<3> wave3d(npoints_3d, grid_spacing_3d, c, bounds_3d, ntsteps, bname, out_freq, src_location_3d);
  wave3d.set_write_output(write_output_files);

  // run
  try {
    wave3d.run();
  } catch (...) {
    std::cerr << "something with wrong ..." << std::endl;
  }

  return 0;
}

