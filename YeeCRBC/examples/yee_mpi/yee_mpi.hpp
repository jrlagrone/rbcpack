/*
  Copyright 2016 John LaGrone

  This file is part of RBCPACK.

  The Yee CRBC Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or (at your 
  option) any later version.

  The RBCPACK Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the RBCPACK Library.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
   3-D FDTD code with CRBC/DAB boundary conditions

   This C++ code implements the finite difference time-domain solution of 
   Maxwell's equation using standard second order, centered diffences. The grid 
   is terminated using Double Absorbing Boundaries (DAB).

   This is meant to serve as an example of how one might use the underlying
   C++ interface to the CRBC/DAB library with MPI. Note that we do not
   officially support MPI in the library, but it is possible to use it in an MPI
   environment in some situtations. This implementation is meant to demonstrate
   such situation.

   Note that the library should ideally be compiled with MPI compilers.

*/ 

#ifndef CRBC_YEE_MPI_H_
#define CRBC_YEE_MPI_H_

// Includes

// This file contains the most generic C++ interface to the CRBC/DAB boundary
// conditions.
#include "crbc_updates.hpp"

// Header file for exact solution routines
#include "solutions.hpp"

// for std::vector
#include <vector>

// include MPI for parallel implementation
#include <mpi.h>

/// \class yee_updater
///
/// This class implements a Yee scheme solver with DAB boundaries with MPI. 
/// We note that the class assumes that the domain is a perfect cube and will
/// only utilize a perfect cube number of processes.
class yee_updater
{

public:

  /// constructor
  /// \param[in] comm      MPI communicator
  /// \param[in] nprocs    the number of processes to use in each direction
  /// \param[in] w         the approximate domain width (may be changed slighty due to discretization)
  /// \param[in] T         the total simulation time
  /// \param[in] CRBC_T    CRBC time parameter (usually CRBC_T = T)
  /// \param[in] CRBC_tol  the tolerance for the DAB boundaries
  /// \param[in] io_t      approximately how often to generate output
  /// \param[in] skip      the stride to use when sampling errors
  /// \param[in] eps       permittivity
  /// \param[in] mu        permeability
  /// \param[in] gamma     roughty 1/variance of the Gaussian pulse used for ICs
  /// \param[in] tau       how far in the past the source pulse was turned on (> 0)
  /// \param[in] dab_wt    weight factor for DAB in load balancing 
  yee_updater(MPI_Comm comm,          
              const int &nprocs,     
	      const double &w,        
	      const double &h,        
	      const double &T,        
	      const double &CRBC_T,   
	      const double &CRBC_tol, 
              const double &io_t = 0.05,
              const int &skip = 1,    
	      const double &eps = 1.0,      
	      const double &mu  = 1.0,
              const double &gamma = 160.0,
              const double &tau = 0.35,
              const int &dab_wt = 3);
 
  /// run the simulation
  void run();

  /// free the communicators created internally
  void free_comms();

  /// function to display the approximate memory usage in MB
  void print_mem_use() const;

  /// function to display timing information
  void print_timing_data() const;

private:

  // storage for field values
  std::vector<double> E[3], H[3];

  double Hcoef, Ecoef;

  // storage for mpi messages
  std::vector<double> E_sbuf[6];
  std::vector<double> E_rbuf[6];
  std::vector<double> DAB_sbuf[6];
  std::vector<double> DAB_rbuf[6];
  std::vector<double> DAB_corner_sbuf[12];
  std::vector<double> DAB_corner_rbuf[12];

  double eps, mu, gamma, tau, io_t, c;
  double T, dt, h, Etime, Htime;
  double tol, CRBC_T;
  double coord[3], domain_width;
  int CRBC_P[6];
  int nprocs, nprocs_cubed;
  int nx, ny, nz;
  int maxn;
  int ntsteps;
  int skip;
  int dab_wt;
  bool isBoundaryProc;
  crbc::BoundaryProperties::Boundary procBounds[6];
  int MPI_DIR[6];
  int my_id, cart_rank[3];
  std::vector<int> send_dirs, send_mpi_dirs, send_sides[4], corner_mpi_dirs[12];
  std::array<int, 3> send_corners[12];
  crbc::CrbcUpdates<3, double, int> bound_upd[3];
  MPI_Comm grid_comm, glob_comm;
  std::vector<MPI_Request> send_req, recv_req;
  std::vector<MPI_Request> DAB_send_req, DAB_recv_req;
  std::vector<int> DAB_props, rDAB_props;
  std::vector<double> DAB_refl, rDAB_refl;
  double scatter_radius;
  double create_comm_t,  alloc_mem_t, init_dab_t, step_E_t, step_inner_H_t, \
         step_outer_H_t, step_DAB_t, send_DAB_t, recv_DAB_t, send_E_t, \
         recv_E_t, sol_t, load_init_conds_t, calc_norm_t, calc_err_t, calc_params_t;

  double dab_mem_use;

  maxwell_solutions::MW_FreeSpace sol_obj;

  /// function to calculate parameters and do some basic load balancing
  void calc_params();

  /// function to create the internal mpi comm. It also labels the boundaries.
  void create_mpi_comm(MPI_Comm comm);

  /// function to set up the solution routines
  void init_solutions();

  /// function to allocate memory
  void allocate_memory();

  /// function that sets up the DAB boundary updaters
  void init_DAB();

  /// evolve the E fields
  void step_E();

  /// evolve the interior H fields
  void step_inner_H();

  /// evolve the H fields on the process boundaries
  void step_outer_H();

  /// update the DAB layers
  void step_DAB();

  /// copy the DAB updates back into the interior
  void copy_DAB();

  /// often used loop that copies data from the DAB to the interior
  void get_dab_vals_loop(std::vector<double> &buffer,
                         crbc::CrbcUpdates<3, double, int> &updater,
                         const int &side,
                         const int *low,
                         const int *high,
                         const int *plow,
                         const int *phigh,
                         const bool &isedge = false);

  /// often used loop that copies data from the interior to the DAB
  void set_dab_vals_loop(std::vector<double> &buffer,
                         crbc::CrbcUpdates<3, double, int> &updater,
                         int &count,
                         const int &side,
                         const int *low,
                         const int *high,
                         const int *plow,
                         const int *phigh,
                         const bool &isedge = false);

  /// function to indentify the sides and edges that need to be sent to update
  /// the DAB layer
  void calc_DAB_send_params();

  /// send DAB values between processes
  void send_DAB();
 
  /// recieve DAB values from neighboring processes
  void recv_DAB();
 
  /// send E field values to neighboring processes
  void send_E();

  /// recieve E field values form neighboring processes
  void recv_E();

  /// load the initial conditions
  void load_initial_conds(); 

  /// calculate the norm at the current time 
  double calc_norm();
  
  /// calculate the error at the current time
  double calc_error();

  void writeExField(int id);
                
};

#endif // CRBC_YEE_MPI_H_
