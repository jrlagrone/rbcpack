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

// We use this for std::accumulate
#include <numeric>

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

  // compute the wave speed;
  c = 1.0 / std::sqrt(eps*mu); 

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
      std::cerr << "id = " << my_id 
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
      std::cerr << "id = " << my_id 
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
      std::cerr << "id = " << my_id 
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
      std::cerr << "id = " << my_id 
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
      std::cerr << "id = " 
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
      std::cerr << "id = " << my_id 
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

/*******************************************************************************
               function to print timing data
*******************************************************************************/
void yee_updater::print_timing_data() const 
{

  if (grid_comm != MPI_COMM_NULL) {

    std::vector<double> timers_send, timers_recv;

    timers_send.reserve(15);
    timers_recv.assign(15*nprocs_cubed, 0.0);

    // save all the local timers to a vector
    timers_send.push_back(calc_params_t);
    timers_send.push_back(create_comm_t);
    timers_send.push_back(alloc_mem_t);
    timers_send.push_back(init_dab_t);
    timers_send.push_back(step_E_t);
    timers_send.push_back(step_inner_H_t);
    timers_send.push_back(step_outer_H_t);
    timers_send.push_back(step_DAB_t);
    timers_send.push_back(send_DAB_t);
    timers_send.push_back(recv_DAB_t);
    timers_send.push_back(send_E_t);
    timers_send.push_back(recv_E_t);
    timers_send.push_back(load_init_conds_t);
    timers_send.push_back(calc_norm_t);
    timers_send.push_back(calc_err_t);

    // gather all timer data
    MPI_Gather(timers_send.data(), timers_send.size(), MPI_DOUBLE, timers_recv.data(), timers_send.size(), MPI_DOUBLE, 0, grid_comm);

    // print out timer data
    if (my_id == 0) {
      std::cout << std::endl << std::endl;
      std::cout << " ,"
		<< " Calculating Parameters,"
		<< " Creating Communicators,"
		<< " Allocating Memory,"
		<< " Initializing DABs,"
		<< " Stepping E,"
		<< " Stepping Inner H,"
		<< " Stepping Outer H,"
		<< " Stepping DAB,"
		<< " Sending DAB,"
		<< " Receiving DAB,"
		<< " Sending E,"
		<< " Receiving E,"
		<< " Loading Initial Condition,"
		<< " Calculating Norm,"
		<< " Calculating Error,"
		<< " Total "
		<< std::endl;

      timers_send.assign(16, 0.0);
      for (int i = 0; i<nprocs_cubed; ++i) {
	double sum = std::accumulate(timers_recv.begin() + 15*i, timers_recv.begin() + 15*(i+1), 0.0);

        // print timer data for each process
	std::cout << " process " << i << ",";
        for (int j=0; j<15; ++j)
          std::cout << timers_recv[15*i+j] << ",";
        std::cout << sum << std::endl;

        // update total times
        for (int j=0; j<15; ++j)
	  timers_send[j] += timers_recv[15*i + j];
	timers_send[15] += sum;
      }

      std::cout << " Total " << ",";
      for (int j=0; j<15; ++j)
	std::cout << timers_send[j] << ",";
      std::cout << timers_send[15] << std::endl << std::endl;

      // compute percentages 
      timers_send.assign(16, 0.0);
      for (int i = 0; i<nprocs_cubed; ++i) {
	std::cout << " process " << i << ",";
        for (int j=0; j<15; ++j)
          std::cout << 100*timers_recv[15*i + j] / std::accumulate(timers_recv.begin() + 15*i, timers_recv.begin() + 15*i + 15, 0.0) << ",";
        std::cout << "100" << std::endl;

        // total percentages
        for (int j=0; j<15; ++j)
	  timers_send[j] += timers_recv[15*i+j] / std::accumulate(timers_recv.begin() + 15*i, timers_recv.begin() + 15*i + 15, 0.0);
	timers_send[15] += 1;
      }

      std::cout << " Average " << ",";
      for (int j=0; j<15; ++j)
        std::cout << 100*timers_send[j]/nprocs_cubed << ",";
      std::cout << std::endl << std::endl;
    } // end if id == 0
  } // end comm check
} // end print_timing_data


/*******************************************************************************
               function to calculate parameters
*******************************************************************************/
void yee_updater::calc_params() 
{

  double eta, emax, cosines[40];
  int x, i, j, k, rem, n[3], P;
  double t1, t2;

  t1 = MPI_Wtime();  // start timer

  // \delta/(cT) where delta is the seperation of the boundary from the source
  // we assume the source is centered in the domain
  eta = (domain_width / 2.0 ) / (c*CRBC_T);

  // figure out the expected number of recursoins
  optimal_cosines(eta, 20, tol, cosines, &P, &emax);

  // compute the left-most point in each coordinate direction
  coord[0] = -domain_width/2.0;
  coord[1] = -domain_width/2.0;
  coord[2] = -domain_width/2.0;

  // Do some basic load balancing. The idea that using P auxilliary variables
  // is roughly equal (in terms of FLOPS) to doing ~3*P Yee cell updates, so
  // we will make the processes with DAB layers have fewer Yee updates to
  // compute but otherwise distribute the points as evenly as possible.
  // Note that the difference between the DAB and Yee memory access patterns
  // probably plays a role here, but we're ignoring it.
  //
  // Also note that we are overlaping the processor domains by a factor of h.
  // This is certainly less memory efficient, but it makes the message 
  // passing and updates a bit more straightforward
  if (nprocs == 1) {
    
    // if there's only 1 MPI process, we just compute the number of grid points
    // like normal
    x = domain_width/h + 1;
    for (i=0; i<3; ++i)
      n[i] = x;

  } else {

    // we start by calculating the total number of "grid points" assuming that 
    // each additional auxilliary variable in the DAB counts for dab_wt grid 
    // points. We additionally assume that each of the boundaries is in fact
    // a DAB, that is we assume the free space problem. We also assume the 
    // grid overlaps by h so processes share 2 grid points with their neighbors.
    //
    // So we want the domain to look something like
    // ......       ...........       ...........       ......
    //     ...........       ...........       ...........
    //
    // So if we had a single process we would have 
    //   (domain_width/h + 1) points
    // plus there are 2 DABS, which give
    //   (2*dab_wt*P) points
    // plus we have the extra points due to overlapping the grid, which gives
    //   ((nprocs-1)*2) points
    // putting this all together, we get that we need a total of 
    //   (domain_width/h + 2*dab_wt*P + 2*nprocs -1)
    // We try to divide this up evenly and calculate how many points are left
    //
    // IMPORTANT:
    // This can fail by making processes have too few points on the boundary or
    // or even assigning a negative number of points on the boundary processes.
    // This isn't really accounted for, but we attempt to abort if we see it ...
    x = ((int) (domain_width/h + 2*dab_wt*P + 2*(nprocs) - 1)) / nprocs;
    rem = ((int) (domain_width/h + 2*dab_wt*P + 2*(nprocs) - 1)) % nprocs;

    // Next we allocate points to processes in each direction bases on whether
    // they are on the boundary or not
    for (i=0; i<3; ++i) {
      if ((cart_rank[i] == 0) || (cart_rank[i] == nprocs-1)) {

        // if the process is on the boundary, we subtract of dab_wt*P points
        // to account for the DAB layer. We additionally calculate the left
        // most coordinate in this direction. Note that on the left side, we 
        // have already correctly set this value so we only do it if it is the
        // right-most process
	n[i] = x - dab_wt*P;
	if (cart_rank[i] == nprocs-1) {
	  coord[i] += (x - dab_wt*P - 2)*h + (cart_rank[i]-1)*(x-2)*h;
	}
      } else {
 
        // otherwise, we just assign the number of points as is and calculate
        // the left-most point of the process in the current direction
	n[i] = x;
	coord[i] += (x - dab_wt*P - 2)*h + (cart_rank[i]-1)*(x-2)*h;
      }
    }

    // now account for any left over points
    if (nprocs == 2) {
      // if there are only 2 processes per direction, just add the extra point(s)
      // to the left process and shift the right process' coordinates accordingly.
      for (i=0; i<3; ++i)
	if (cart_rank[i] == 0) {
	  n[i] += rem;
	} else {
	  coord[i] += rem*h;
	}
    } else {
   
      // otherwise we only add extra points to the interior processes. We do 
      // this by looping over the interior processes from left to right and
      // add one to the current process and the coordinates by h for all of the
      // processes to the right and repeat until we have no remaining points.
      int r[3];
      r[0] = rem;
      r[1] = rem;
      r[2] = rem;

      // loop over the number of remaining points just to make sure we iterate
      // enough times
      for (k=0; k<rem; ++k) {
        // loop over the interior processes
	for (j=1; j<nprocs-1; ++j) {
          // loop over the directions
	  for (i=0; i<3; ++i) {
            // if we have points left add one to the current process and shift
            // the process the coordinates for the processes to the right
	    if (r[i] > 0) {
	      if (cart_rank[i] == j) {
		n[i]++;
		r[i]--;
	      }
	      if (cart_rank[i] > j)
		coord[i] += h;
	    }
	  }
	}
      } // end for k
    }
  }

  // do a basic check to make sure that the grid partitioning is somewhat 
  // reasonable
  if ((n[0] < 3) || (n[1] < 3) || (n[2] < 3)) {
    std::cerr << "Grid partitioning failed. Try decreasing h, dab_wt, and/or nprocs" << std::endl;
      MPI_Abort(glob_comm, -3);
  }

  // save the number of grid points in each direction
  nx = n[0];
  ny = n[1];
  nz = n[2];

  // calculate the time step size and number of time steps
  dt = 0.99 / sqrt(3.0/(h*h));
  ntsteps = T / dt;
  Etime = 0;
  Htime = dt/2.0;

  // update timer
  t2 = MPI_Wtime();
  calc_params_t += t2-t1;
}


/*******************************************************************************
               function to create MPI communicator
*******************************************************************************/
void yee_updater::create_mpi_comm(MPI_Comm comm) {
  int periods[3];
  periods[0] = 0; // not periodic
  periods[1] = 0;
  periods[2] = 0;
  int reorder = 1;
  double t1, t2;

  t1 = MPI_Wtime(); // start timer

  int dims[3];
  dims[0] = nprocs;
  dims[1] = nprocs;
  dims[2] = nprocs;

  nprocs_cubed = nprocs*nprocs*nprocs;
  
  // create a cartesian communicator with nprocs in each direction
  if (MPI_Cart_create(comm, 3, dims, periods, reorder, &grid_comm) != MPI_SUCCESS) {
    std::cerr << "MPI_Cart_create failed" << std::endl;
  }

  // figure out neighboring processes
  NORTH = MPI_PROC_NULL;
  EAST = MPI_PROC_NULL;
  SOUTH = MPI_PROC_NULL;
  WEST = MPI_PROC_NULL;
  UP = MPI_PROC_NULL;
  DOWN = MPI_PROC_NULL;

  cart_rank[0] = -1;
  cart_rank[1] = -1;
  cart_rank[2] = -1;

  if (grid_comm != MPI_COMM_NULL) {
    if (MPI_Comm_rank(grid_comm, &my_id) != MPI_SUCCESS)
      std::cerr << "MPI_Comm_rank failed" << std::endl;
    if (MPI_Cart_coords(grid_comm, my_id, 3, cart_rank) != MPI_SUCCESS)
      std::cerr << "MPI_Cart_coords failed" << std::endl;

    // figure the ids of the processes we might need to send data to
    if (MPI_Cart_shift(grid_comm, 2, -1, &UP, &DOWN) != MPI_SUCCESS)
      std::cerr << "MPI_Cart_shift failed" << std::endl;

    if (MPI_Cart_shift(grid_comm, 1, -1, &NORTH, &SOUTH) != MPI_SUCCESS)
      std::cerr << "MPI_Cart_shift failed" << std::endl;

    if (MPI_Cart_shift(grid_comm, 0, -1, &EAST, &WEST) != MPI_SUCCESS)
      std::cerr << "MPI_Cart_shift failed" << std::endl;

  }

  // figure out which processes are on the boundary
  isBoundaryProc = false;

  // if a process doesn't have a neighbor on at least 1 side, its on the boundary
  if ((grid_comm != MPI_COMM_NULL) && ((NORTH == MPI_PROC_NULL) || 
       (EAST == MPI_PROC_NULL) || (SOUTH == MPI_PROC_NULL) || 
       (WEST == MPI_PROC_NULL) || (UP == MPI_PROC_NULL) || 
       (DOWN == MPI_PROC_NULL)))
    isBoundaryProc = true;

  // label the boundaries. Use type NONE for interior sides and CRBC for the
  // exterior boundaries. To do a wave guide, e.g., one might change the type
  // to PEC on the appropriate sides
  for (int i=0; i<6; ++i)
    procBounds[i] = crbc::BoundaryProperties::NONE;

  if ((grid_comm != MPI_COMM_NULL) && (UP == MPI_PROC_NULL))
    procBounds[5] = crbc::BoundaryProperties::CRBC;
   
  if ((grid_comm != MPI_COMM_NULL) && (DOWN == MPI_PROC_NULL))
    procBounds[4] = crbc::BoundaryProperties::CRBC;

  if ((grid_comm != MPI_COMM_NULL) && (NORTH == MPI_PROC_NULL))
    procBounds[3] = crbc::BoundaryProperties::CRBC;

  if ((grid_comm != MPI_COMM_NULL) && (SOUTH == MPI_PROC_NULL))
    procBounds[2] = crbc::BoundaryProperties::CRBC;

  if ((grid_comm != MPI_COMM_NULL) && (EAST == MPI_PROC_NULL))
    procBounds[1] = crbc::BoundaryProperties::CRBC;

  if ((grid_comm != MPI_COMM_NULL) && (WEST == MPI_PROC_NULL))
    procBounds[0] = crbc::BoundaryProperties::CRBC;

  // stop timer
  t2 = MPI_Wtime();
  create_comm_t += t2-t1;
}


/*******************************************************************************
               function to initialize solutions
*******************************************************************************/
void yee_updater::init_solutions() 
{
  
  double src_loc[3];
  double hloc[3];

  // place the source at (0,0,0) plus a small perturbation to decrease the 
  // chances of coinciding with a grid point
  src_loc[0] = 1.43e-7;
  src_loc[1] = -5.3423e-7;
  src_loc[2] = 9.012e-8;

  // set the grid spacing to be the same in all directions
  hloc[0] = h;
  hloc[1] = h;
  hloc[2] = h;
  

  sol_obj = maxwell_solutions::MW_FreeSpace(gamma, tau, eps, mu, src_loc);
  sol_obj.set_grid_spacing(hloc);

}

/*******************************************************************************
               function to allocate memory
*******************************************************************************/
void yee_updater::allocate_memory()
{

  int m;
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  // figure out the largest number number of grid points possible
  m = (nx > ny) ? nx : ny;
  m = (m > nz) ? m : nz;
  if (MPI_Allreduce(&m, &maxn, 1, MPI_INT, MPI_MAX, grid_comm) != MPI_SUCCESS)
    std::cerr << "MPI_Allreduce failed";

  // allocate Fields and initialize to 0
  Ex.assign((nx-1)*ny*nz, 0.0);
  Ey.assign(nx*(ny-1)*nz, 0.0);
  Ez.assign(nx*ny*(nz-1), 0.0);
  Hx.assign(nx*(ny-1)*(nz-1), 0.0);
  Hy.assign((nx-1)*ny*(nz-1), 0.0);
  Hz.assign((nx-1)*(ny-1)*nz, 0.0);

  // compute update coefficients
  Hcoef = (dt/mu)/h;
  Ecoef = (dt/eps)/h;

  // allocate buffer space for the worst possible case. We could use slightly
  // less memory here by figuring out the exact sizes we need.
  if (procBounds[0] == crbc::BoundaryProperties::NONE) {
    W_sbuf.assign(2*maxn*maxn, 0.0);
    W_rbuf.assign(2*maxn*maxn, 0.0);
  }
  if (procBounds[1] == crbc::BoundaryProperties::NONE) {
    E_sbuf.assign(2*maxn*maxn, 0.0);
    E_rbuf.assign(2*maxn*maxn, 0.0);
  }
  if (procBounds[2] == crbc::BoundaryProperties::NONE) {
    S_sbuf.assign(2*maxn*maxn, 0.0);
    S_rbuf.assign(2*maxn*maxn, 0.0);
  }
  if (procBounds[3] == crbc::BoundaryProperties::NONE) {
    N_sbuf.assign(2*maxn*maxn, 0.0);
    N_rbuf.assign(2*maxn*maxn, 0.0);
  }
  if (procBounds[4] == crbc::BoundaryProperties::NONE) {
    D_sbuf.assign(2*maxn*maxn, 0.0);
    D_rbuf.assign(2*maxn*maxn, 0.0);
  }
  if (procBounds[5] == crbc::BoundaryProperties::NONE) {
    U_sbuf.assign(2*maxn*maxn, 0.0);
    U_rbuf.assign(2*maxn*maxn, 0.0);
  }   

  // stop timer
  t2 = MPI_Wtime();
  alloc_mem_t += t2-t1;
}

/*******************************************************************************
                Initialize the DAB updaters
*******************************************************************************/
void yee_updater::init_DAB()
{

  double delta;
  int l,m;
  int low[3], high[3];
  double htmp[3];
  htmp[0] = htmp[1] = htmp[2] = h;
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  // storage for DAB boundary properties that we might want to print out or need
  // to use elsewhere
  if (my_id == 0) {
    rDAB_props.assign(15*nprocs_cubed, 0.0);
    rDAB_refl.assign(10*nprocs_cubed, 0.0);
  }
  DAB_props.assign(15, 0);
  DAB_refl.assign(10, 0.0);
  DAB_refl[6] = coord[0];
  DAB_refl[7] = coord[1];
  DAB_refl[8] = coord[2];


  if (isBoundaryProc) {

    // we're assuming this is a free space problem, so we initialize DAB updaters
    // for all 3 E components. In a wave guide, e.g., where there are no DAB
    // edges or corners, we only need updaters for the tangential components.
    // There's no harm in having an updater for the normal components when it's
    // not needed, but it represents unnecessary work
    //
    // We initialize the updaters in 3D with double field values and ints for 
    // indexing (and by default doubles for coeficients) and provide the run
    // time, grid spacing, time step size, wave speed, and boundary configuration
    // on each boundary process
    bound_upd_Ex = crbc::CrbcUpdates<3, double, int> (CRBC_T, htmp, dt, c, procBounds);
    bound_upd_Ey = crbc::CrbcUpdates<3, double, int> (CRBC_T, htmp, dt, c, procBounds);
    bound_upd_Ez = crbc::CrbcUpdates<3, double, int> (CRBC_T, htmp, dt, c, procBounds);

    // loop over and set up each of the possible boundary sides
    //
    // We are dealing with the MPI by overlapping the DAB layer for tangential
    // components, the grids are overlapped by a factor of h so neighboring
    // processes share 2 physical points. The same thing is done for the normal
    // component, but it looks slightly different because of the staggered grid.
    // While they both overlap by a factor of h to share 2 grid points, in the 
    // Yee scheme it is 1 cell for tangential components and 2 cells for normal
    // components.
    //
    // Tangential components            Normal Components
    //   ---------                      ---x-----x-----x--- 
    //   x   x   x                      |     |     |     |
    //   ---------                      ---x-----x-----x---
    //       ---------                        ---x-----x-----x---
    //       x   x   x                        |     |     |     |
    //       ---------                        ---x-----x-----x---
    //
    // Note that the use of normal and tangential components here is somewhat
    // confusing because it is in reference to the boundaries with neighboring
    // processes, NOT the phyiscal boundary. We consider the direction in which
    // the message passing needs to take place as the normal direction. This 
    // results in the following, e.g.
    //
    // With our implementation, each process has components with the following 
    // indexing bounds:
    // Ex(0:nx-2, 0:ny-1, 0:nz-1) located at ((i+1/2)*h, j*h, k*h)
    // Ey(0:nx-1, 0:ny-2, 0:nz-1) located at (i*h, (j+1/2)*h, k*h)
    // Ez(0:nx-1, 0:ny-1, 0:nz-2) located at (i*h, j*h, (k+1/2)*h)
    //
    // So, for example, we'll consider the right boundary face in the x
    // direction. Then, we potentially need to pass information North/South 
    // (in the y-direction) or Up/Down (in the z-direction). For the case of 
    // needed to pass information in the North/South direction, the Ey 
    // component is normal to the interface between the two processes and Ex and
    // Ez are tangential. The tangential components are already overlap the way
    // we want because we overlapped the grids for the Yee scheme, therefore, we
    // tell the DAB updater the actual data extents for the points:
    //
    // For Ex the proper extents are [nx-3, nx-2] x [0, ny-1] x [0, nz-1]
    // because we include all the points in the y and z directions and the point
    // in x on the physical boundary and it's neighbor to the left.
    //
    // Similary, for Ez the extents are [nx-2, nx-1] x [0, ny-1] x [0, nz-2]
    //
    // For Ey, if we do the same thing, we would get the extents
    // [nx-2, nx-1] x [0, ny-2] x [0, nz-1], but this does not overlap the grids
    // by 2 points. To correct this, we tell the DAB layer that the extents are
    // greater than the actual data range by subtracting 1 from the lower y
    // extent if their is a neighboring process in the DOWN direction to get
    // [nx-2, nx-1] x [-1, ny-2] x [0, nz-1]
    // 
    // NOTE: the DAB updater considers the extents to be inclusive
    for (l=0; l<6; ++l) {

      if (procBounds[l] == crbc::BoundaryProperties::CRBC) {
        
	// figure out seperation from source
	delta = domain_width / 2.0;

	// loop over field components
	for (m=0; m<3; ++m) {

          // generic extents that are close to correct (some of the indices are
          // off by a factor of 1, which depends on the component)
	  low[0] = 0;
	  low[1] = 0;
	  low[2] = 0;
	  high[0] = nx - 1;
	  high[1] = ny - 1;
	  high[2] = nz - 1; 

	  if (l == 0) { 
	    // left boundary in x, need [0,1] in x, all in y, z
	    high[0] = 1;

	    // adjust based on field component
	    if (m == 1) { // Ey
	      high[1]--;
	      if (SOUTH != MPI_PROC_NULL)
		low[1]--;
	    } else if (m == 2) { // Ez
              high[2]--;
	      if (DOWN != MPI_PROC_NULL)
		low[2]--;
	    }
              
	  } else if (l == 1) { 
	    // right boundary in x, need [nx-2, nx-1] in x, all y, z for 
            // tangential and [nx-3, nx-2] in x, all y, z for normal
	    low[0] = nx-2;
    
	    // adjust based on field component
	    if (m == 0) {
	      high[0]--;
	      low[0]--;
	    } else if (m == 1) { // Ey
	      high[1]--;
	      if (SOUTH != MPI_PROC_NULL)
		low[1]--;
	    } else if (m == 2) { // Ez
	      high[2]--;
	      if (DOWN != MPI_PROC_NULL)
		low[2]--;
	    } 

	  } else if (l == 2) {               
	    // left boundary in y, need [0,1] in y, all in x, z
	    high[1] = 1;

	    // adjust based on field component
	    if (m == 0) { // Ex
	      high[0]--;
	      if (WEST != MPI_PROC_NULL)
		low[0]--;
	    } else if (m == 2) { // Ez
	      high[2]--;
	      if (DOWN != MPI_PROC_NULL)
		low[2]--;
	    }
	  } else if (l == 3) {
	    // right boundary in y, need [ny-2, ny-1] in y, all x, z for 
            // tangential and [ny-3, ny-2] in y, all x, z
	    low[1] = ny-2;
    
	    // adjust based on field component
	    if (m == 1) {
	      high[1]--;
	      low[1]--;
	    } else if (m == 0) { // Ex
	      high[0]--;
	      if (WEST != MPI_PROC_NULL)
		low[0]--;
	    } else if (m == 2) { // Ez
	      high[2]--;
	      if (DOWN != MPI_PROC_NULL)
		low[2]--;
	    } 
	  } else if (l == 4) {               
	    // left boundary in z, need [0,1] in z, all in x, y
	    high[2] = 1;

	    // adjust based on field component
	    if (m == 0) { // Ex
	      high[0]--;
	      if (WEST != MPI_PROC_NULL)
		low[0]--;
	    } else if (m == 1) { // Ey
	      high[1]--;
	      if (SOUTH != MPI_PROC_NULL)
		low[1]--;
	    }
	  } else if (l == 5) {
	    // right boundary in z, need [nz-2, nz-1] in z, all x, y for 
            // tangential and [nz-3, nz-2] in z, all x, y for normal
	    low[2] = nz-2;
    
	    // adjust based on field component
	    if (m == 2) {
	      high[2]--;
	      low[2]--;
	    } else if (m == 0) { // Ex
	      high[0]--;
	      if (WEST != MPI_PROC_NULL)
		low[0]--;
	    } else if (m == 1) { // Ey
	      high[1]--;
	      if (SOUTH != MPI_PROC_NULL)
		low[1]--;
	    } 
	  }

	  // call initializer and limit the number of recursions to at most 20
	  if (m == 0) {
	    bound_upd_Ex.init_face(l, low, high, delta, 20, tol);
	  } else if (m == 1) {
	    bound_upd_Ey.init_face(l, low, high, delta, 20, tol);
	  } else {
	    bound_upd_Ez.init_face(l, low, high, delta, 20, tol);
	  }

	} // end loop over components
      }
    } // end loop over sides

    // now get some properties from the updaters that we may be interested in
    // at a later time
    dab_mem_use = bound_upd_Ex.get_mem_usage() + bound_upd_Ey.get_mem_usage() \
                + bound_upd_Ez.get_mem_usage();

    // get number of recursions and reflection coefficients
    for (l = 0; l<6; ++l) {
      if (procBounds[l] == crbc::BoundaryProperties::CRBC) {
	DAB_props.at(l) = bound_upd_Ex.get_num_recursions(l);
	DAB_refl.at(l) = bound_upd_Ex.get_reflection_coef(l);
      }
    }
    DAB_refl[9] = dab_mem_use;

    // get info about the domain configuration 
    DAB_props.at(6) = bound_upd_Ex.get_num_faces();
    DAB_props.at(7) = bound_upd_Ex.get_num_edges();
    DAB_props.at(8) = bound_upd_Ex.get_num_corners();
    DAB_props.at(9) = bound_upd_Ey.get_num_faces();
    DAB_props.at(10) = bound_upd_Ey.get_num_edges();
    DAB_props.at(11) = bound_upd_Ey.get_num_corners();
    DAB_props.at(12) = bound_upd_Ez.get_num_faces();
    DAB_props.at(13) = bound_upd_Ez.get_num_edges();
    DAB_props.at(14) = bound_upd_Ez.get_num_corners();

  }

 /*
  // print out DAB info
  if (MPI_Gather(DAB_props.data(), 15, MPI_INT, rDAB_props.data(), 15, MPI_INT, 0, grid_comm) != MPI_SUCCESS)
    std::cerr << "MPI_Gather failed" << std::endl;
  if (MPI_Gather(DAB_refl.data(), 10, MPI_DOUBLE, rDAB_refl.data(), 10, MPI_DOUBLE, 0, grid_comm) != MPI_SUCCESS)
    std::cerr << "MPI_Gather failed" << std::endl;
  if (my_id == 0) {

    std::cout << "Domain width = " << domain_width << std::endl << std::endl;
    double mem_total = 0;

    for (l=0; l<nprocs_cubed; ++l) {

      int gc[3];

      MPI_Cart_coords(grid_comm, l, 3, gc);

      std::cout << "Process " << l << " with logical coordinates (" << gc[0] << ", "
		<< gc[1] << ", " << gc[2] << ") for the Ex component" << std::endl;
      std::cout << "         with a corner at (" << rDAB_refl.at(10*l + 6) << ", "
		<< rDAB_refl.at(10*l + 7) << ", " << rDAB_refl.at(10*l + 8) << ")" << std::endl;

      std::cout << "  Left side in x:" << std::endl;
      std::cout << "    recursions      = " << rDAB_props.at(15*l) << std::endl;
      std::cout << "    reflection coef = " << rDAB_refl.at(10*l) << std::endl;
      std::cout << "  Right side in x:" << std::endl;
      std::cout << "    recursions      = " << rDAB_props.at(15*l + 1) << std::endl;
      std::cout << "    reflection coef = " << rDAB_refl.at(10*l + 1) << std::endl;
      std::cout << "  Left side in y:" << std::endl;
      std::cout << "    recursions      = " << rDAB_props.at(15*l + 2) << std::endl;
      std::cout << "    reflection coef = " << rDAB_refl.at(10*l + 2) << std::endl;
      std::cout << "  Right side in y:" << std::endl;
      std::cout << "    recursions      = " << rDAB_props.at(15*l + 3) << std::endl;
      std::cout << "    reflection coef = " << rDAB_refl.at(10*l + 3) << std::endl;
      std::cout << "  Left side in z:" << std::endl;
      std::cout << "    recursions      = " << rDAB_props.at(15*l + 4) << std::endl;
      std::cout << "    reflection coef = " << rDAB_refl.at(10*l + 4) << std::endl;
      std::cout << "  Right side in z:" << std::endl;
      std::cout << "    recursions      = " << rDAB_props.at(15*l + 5) << std::endl;
      std::cout << "    reflection coef = " << rDAB_refl.at(10*l + 5) << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 6) << " Ex faces" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 7) << " Ex edges" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 8) << " Ex corners" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 9) << " Ey faces" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 10) << " Ey edges" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 11) << " Ey corners" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 12) << " Ez faces" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 13) << " Ez edges" << std::endl;
      std::cout << "  " << rDAB_props.at(15*l + 14) << " Ez corners" << std::endl;
      std::cout << "  Approx DAB Memory Use = " << rDAB_refl.at(10*l + 9) << " MB" << std::endl;
      mem_total += rDAB_refl.at(10*l + 9);

    }
    std::cout << "  Approx Total DAB Memory Use = " << mem_total << " MB" << std::endl;
  }
*/

  // stop timer
  t2 = MPI_Wtime();
  init_dab_t += t2-t1;
}

