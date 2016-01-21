/**
   3-D FDTD code with CRBC/DAB boundary conditions

   This C++ code implements the finite difference time-domain solution of 
   Maxwell's equation using standard second order, centered diffences. The grid 
   is terminated using Double Absorbing Boundaries (DAB).

   This is specialized to simulate scattering from a sphere. Note the we only
   simulate the scattered solution by enforcing dirichlet data on an internal
   boundary using an exact solution.

   This is meant to serve as an example of how one might use the underlying
   C++ interface to the CRBC/DAB library with MPI. This interface is templated 
   and supports up to three dimensions. It is also templated on the internal
   indexing and data types.
    
   This code is freely distributed with the goal of illustrating the CRBC/DAB
   library usage and research purposes. We provide no warranties of its 
   performance or suitability for any purpose.

*/ 


// Includes

// Header file containing the class declarations.
#include "yee_mpi.hpp"

// Needed for C++ output
#include <iostream>

// We use this to get std::swap so we can swap pointers for arrays and vectors
// without having to copy data. Note that as of C++11, std::swap is suppose to 
// move to <utility>.
#include <algorithm>

// This has the declarations of the 'sqrt' and 'fabs' functions
#include <cmath>

/*******************************************************************************
                            constructor
*******************************************************************************/
yee_updater::yee_updater(MPI_Comm comm,          
                         const int &nprocs,     
	                 const double &w,        
	                 const double &h,        
	                 const double &T,        
	                 const double &CRBC_T,   
	                 const double &CRBC_tol, 
                         const double &io_t, 
                         const int &skip,    
	                 const double &eps,      
	                 const double &mu,
                         const double &gamma,
                         const double &tau,
                         const int &dab_wt)
{

  // save inputs
  this->T = T;
  this->CRBC_T = CRBC_T;
  tol = CRBC_tol;
  domain_width = w;
  this->nprocs = nprocs;
  this->h = h;
  this->io_t = io_t;
  this->skip = skip;
  this->eps = eps;
  this->mu = mu;
  this->gamma = gamma;
  this->tau = tau;
  this->dab_wt = dab_wt;
  

  // make a copy of the communicator
  MPI_Comm_dup(comm, &glob_comm);

  // initialize timers
  calc_params_t = create_comm_t = alloc_mem_t = init_dab_t = step_E_t = \
                  step_inner_H_t = step_outer_H_t = step_DAB_t =  send_DAB_t = \
                  recv_DAB_t = send_E_t = recv_E_t = sol_t = load_init_conds_t = \
                  calc_norm_t = calc_err_t = 0.0;

  // create cartesian communicator
  create_mpi_comm(glob_comm);

  // calculate grid and time step size
  if (grid_comm != MPI_COMM_NULL) {
    try {
      calc_params();
    } catch (const std::exception& e) {
      std::cout << "id = " << my_id 
        << " failed in calc_params() --- a standard exception was caught, with message '"
        << e.what() << "'" << std::endl;
      MPI_Abort(glob_comm, -2);
    } 

    // wait for all processes to finish calculating parameters
    MPI_Barrier(grid_comm);
  }

  // allocate memory
  if (grid_comm != MPI_COMM_NULL) {
    try {
      allocate_memory();
    } catch (const std::exception& e) {
      std::cout << "id = " << my_id 
        << " failed in allocate_memory() --- a standard exception was caught, with message '"
        << e.what() << "'" << std::endl;
      MPI_Abort(glob_comm, -2);
    } 
    
    // wait for all processes to finish allocating memory
    MPI_Barrier(grid_comm);
  }

  // initialize DAB updaters
  if (grid_comm != MPI_COMM_NULL) {
    try {
      init_DAB();
    } catch (const std::exception& e) {
      std::cout << "id = " << my_id 
        << " failed in init_DAB() --- a standard exception was caught, with message '"
        << e.what() << "'" << std::endl;
      MPI_Abort(glob_comm, -2);
    } 
     
    // wait for all processes to set up DAB updaters (if needed)
    MPI_Barrier(grid_comm);
  }

  // set up the solution routine
  if (grid_comm != MPI_COMM_NULL) {
    try {
      init_solutions();
    } catch (const std::exception& e) {
      std::cout << "id = " << my_id 
        << " failed in init_solutions() --- a standard exception was caught, with message '"
        << e.what() << "'" << std::endl;
      MPI_Abort(glob_comm, -2);
    } 
     
    // wait for all processes to set up solution routines
    MPI_Barrier(grid_comm);
  }

} // end constructor


/*******************************************************************************
                     function to run the simulation
*******************************************************************************/
void yee_updater::run()
{

  double loc_norm, glob_norm, loc_err, glob_err;
  int tskip; 
  int tstep;
  std::vector<double> err;
  std::vector<double> time;

  // figure out how many time steps should be taken between outputs
  tskip = io_t/dt;

  // reserve memory for errors 
  if (my_id==0) {
    err.reserve((ntsteps / tskip) + 1);
    time.reserve((ntsteps / tskip) + 1);
  }
   
  // load initial conditions
  if (grid_comm != MPI_COMM_NULL) {
    try {
      load_initial_conds();
    } catch (const std::exception& e) {
      std::cout << "id = " 
        << my_id << " failed in load_initial_conds() --- a standard exception was caught, with message '"
        << e.what() << "'" << std::endl;
      MPI_Abort(glob_comm, -2);
    } 
    MPI_Barrier(grid_comm);
  }

  // calculate norm of intitial conditions
  if (grid_comm != MPI_COMM_NULL) {
    try {
      loc_norm = calc_norm();
    } catch (const std::exception& e) {
      std::cout << "id = " << my_id 
        << " failed in calc_norm() --- a standard exception was caught, with message '"
        << e.what() << "'" << std::endl;
        MPI_Abort(glob_comm, -2);
    } 
    MPI_Barrier(grid_comm);
  }

  // use mpi_reduce to calculate global norm
  glob_norm = 0.0;
  if (grid_comm != MPI_COMM_NULL) {
    if (MPI_Reduce(&loc_norm, &glob_norm, 1, MPI_DOUBLE, MPI_SUM, 0, grid_comm) != MPI_SUCCESS)
      std::cerr << "MPI_Reduce failed, norm calculation" << std::endl;
  }
    
  if (my_id == 0) {
    glob_norm = std::sqrt(glob_norm);
  }

  // time step
  if (grid_comm != MPI_COMM_NULL) {

    for (tstep = 0; tstep < ntsteps; ++tstep) {

      // generate output
      if (tstep % tskip == 0) {

        // calculate error
	loc_err = calc_err();

	glob_err = 0.0;

	if (MPI_Reduce(&loc_err, &glob_err, 1, MPI_DOUBLE, MPI_SUM, 0, grid_comm) != MPI_SUCCESS)
	  std::cerr << "MPI_Reduce failed, err calculation" << std::endl;

	if (my_id==0) {
	  std::cout << "tstep = " << tstep << "	T (E) = " << Etime 
	      << "	err = " << std::sqrt(glob_err)
              << "	rel err = " << std::sqrt(glob_err)/glob_norm << std::endl;
	  err.push_back(std::sqrt(glob_err));
	  time.push_back(Etime);
	}
      } // output

      // update E fields
      step_E();

      // Send the E fields
      send_E();
     
      // update the DAB
      step_DAB();

      // update the current E time
      Etime += dt;

      // Send the DAB values
      send_DAB();

      // update the H fields
      step_inner_H();

      // wait for the E field sends to complete
      recv_E();

      // update the boundary H fields
      step_outer_H();

      // wait for the DAB sends to complete
      recv_DAB();

      // increment H time
      Htime += dt;

    } // end time stepping

    // calculate final error
    loc_err = calc_err();
  
    glob_err = 0.0;

    if (MPI_Reduce(&loc_err, &glob_err, 1, MPI_DOUBLE, MPI_SUM, 0, grid_comm) != MPI_SUCCESS)
      std::cerr << "MPI_Reduce failed, err calculation" << std::endl;

    if (my_id == 0) {
      err.push_back(std::sqrt(glob_err));
      time.push_back(Etime);

      std::cout << "tstep = " << tstep << "	T (E) = " << Etime 
	<< "	err = " << std::sqrt(glob_err)
        << "	rel err = " << std::sqrt(glob_err)/glob_norm << std::endl;
 
      std::cout << std::endl << std::endl;
	
      // print out all the errors again that are easier to import
      std::cout << "time, error, relative error," << std::endl;
      for (unsigned i=0; i<err.size(); ++i)
        std::cout << time[i] << ", " << err[i] << ", "  << err[i]/glob_norm 
          << std::endl;
    }

  } // end grid_comm != MPI_COMM_NULL
} // end run()


/*******************************************************************************
                   function to free the communicators
*******************************************************************************/
void yee_updater::free_comms()
{
  MPI_Comm_free(&glob_comm);
  MPI_Comm_free(&grid_comm);
}

/*******************************************************************************
               function to print approximate memory usage
*******************************************************************************/
void yee_updater::print_mem_use() const 
{
  if (grid_comm != MPI_COMM_NULL) {

    double tot_mem_use, dab_buff, ebuff, fields;

    // calculate the size of the field vectors in MB
    fields = sizeof(double)*(Ex.capacity() + Ey.capacity() + Ez.capacity() + 
          Hx.capacity() + Hy.capacity() + Hz.capacity()) / ((double) 1024*1024);

    // calculate the size of the E send/recieve buffers
    ebuff =  sizeof(double)*(N_sbuf.capacity() + S_sbuf.capacity() + 
             E_sbuf.capacity() + W_sbuf.capacity() + U_sbuf.capacity() + 
             D_sbuf.capacity() + N_rbuf.capacity() + S_rbuf.capacity() + 
             E_rbuf.capacity() + W_rbuf.capacity() + U_rbuf.capacity() + 
             D_rbuf.capacity()) / ((double) 1024*1024);

    // calculate the size of the DAB buffers
    dab_buff = 0;
    for (int i=0; i<4; ++i)
      dab_buff += sizeof(double)*(DAB_N_sbuf[i].capacity() + DAB_S_sbuf[i].capacity()
                + DAB_E_sbuf[i].capacity() + DAB_W_sbuf[i].capacity() 
                + DAB_U_sbuf[i].capacity() + DAB_D_sbuf[i].capacity() + 
                  DAB_N_rbuf[i].capacity() + DAB_S_rbuf[i].capacity() + 
                  DAB_E_rbuf[i].capacity() + DAB_W_rbuf[i].capacity() + 
                  DAB_U_rbuf[i].capacity() + DAB_D_rbuf[i].capacity()) 
                / ((double) 1024*1024);

    // calculate total mem use (dab_mem_use is calculate in init_DAB())
    tot_mem_use = dab_buff + ebuff + fields + dab_mem_use;
  
    // colect everything on 1 process
    std::vector<double> send, recv;
    send.push_back(dab_buff);
    send.push_back(ebuff);
    send.push_back(fields);
    send.push_back(dab_mem_use);
    send.push_back(tot_mem_use);
    recv.assign(nprocs_cubed*5, 0);

    if (MPI_Gather(send.data(), 5, MPI_DOUBLE, recv.data(), 5, MPI_DOUBLE, 0, grid_comm) != MPI_SUCCESS)
      std::cerr << "MPI_Gather failed " << std::endl;

    if (my_id == 0) {

      std::cout << " , DAB Buffers, E Buffers, Fields, DAB, Total" << std::endl;
      for (int l=0; l<nprocs_cubed; ++l) {

        std::cout << "proc " << l << ", ";
        for (int i=0; i<5; ++i)
           std::cout << recv[5*l + i] << ", ";

        std::cout << std::endl;   
      }

      std::cout << "Total, ";
      for (int i=0; i<5; ++i) {
        double tot = 0;
        for (int l=0; l<nprocs_cubed; ++l)
          tot += recv[5*l + i];
        std::cout << tot << ", ";
      }
      std::cout << std::endl;
           
    } // end my_id == 0
  } // end if grid_comm != MPI_COMM_NULL
} // end print_mem_use()
