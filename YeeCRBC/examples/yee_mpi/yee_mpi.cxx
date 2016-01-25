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

// optional OpenMP
#if USE_OPENMP
  #include <omp.h>
#endif


// for output, namely fprintf and printf
#include <stdio.h>   

// store and modify output filenames
#include <string.h>  
 
// standard libraries, we need this for malloc, free, etc.
#include <stdlib.h>   


void yee_updater::writeExField(int id) {

    int i, j, k, npts, cells;

    char fileBaseName[100];   
    sprintf(fileBaseName, "Ex_Field_%05d.%05d.vtk", my_id, id);   

    // open the file and write the VTK header information
    FILE *f = fopen(fileBaseName, "w");
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "vtk output\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET RECTILINEAR_GRID\n");

    // set the dimensions
    fprintf(f, "DIMENSIONS %i %i %i\n", nx-1, ny, nz);

    // save the coordinates
    fprintf(f, "X_COORDINATES %i double \n", nx-1);
    for (i=0; i < nx-1; ++i)
      fprintf(f, "%f\n", coord[0] + h/2.0 + i*h);

    fprintf(f, "Y_COORDINATES %i double \n", ny);
    for (j=0; j < ny; ++j)
      fprintf(f, "%f\n", coord[1] + j*h);

    fprintf(f, "Z_COORDINATES %i double \n", nz);
    for (k=0; k < nz; ++k)
      fprintf(f, "%f\n", coord[2] + k*h);

    // set up a cell and field
    npts = (nx-1) * ny * nz;
    cells = (nx-2) * (ny-1) * (nz-1);
    fprintf(f, "CELL_DATA %i\n", cells);
    fprintf(f, "POINT_DATA %i\n", npts);
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "Ex 1 %i double\n", npts);

    // now write out the data
    for (k=0; k < nz; ++k) 
      for (j=0; j < ny; ++j)
	for (i=0; i < nx-1; ++i)
	  fprintf(f, "%f\n", Ex[i + (j + k*ny)*(nx-1)]);
  	

    // close the file
    fclose(f);
  }

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

    int count = 0;

    for (tstep = 0; tstep < ntsteps; ++tstep) {

      // generate output
      if (tstep % tskip == 0) {

        writeExField(count++);

        // calculate error
	loc_err = calc_error();

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
    loc_err = calc_error();
  
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
    ebuff = 0;    
    for (int i=0; i<6; ++i)
      ebuff +=  sizeof(double)*(E_sbuf[i].capacity() + E_rbuf[i].capacity()) \
            / ((double) 1024*1024);

    // calculate the size of the DAB buffers
    dab_buff = 0;
    for (int i=0; i<6; ++i)
      dab_buff += sizeof(double)*(DAB_sbuf[i].capacity() + DAB_rbuf[i].capacity()) \
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
  for (int i=0; i<6; ++i) {
    if (procBounds[i] == crbc::BoundaryProperties::NONE) {
      E_sbuf[i].assign(2*maxn*maxn, 0.0);
      E_rbuf[i].assign(2*maxn*maxn, 0.0);
    }
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

  // figure out the message passing configuration
  calc_DAB_send_params();

  // stop timer
  t2 = MPI_Wtime();
  init_dab_t += t2-t1;
}


/*******************************************************************************
                routine to update E fields
*******************************************************************************/
void yee_updater::step_E()
{
  int i,j,k;
  int nxm, nym, nzm;
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();
  
  #if USE_OPENMP
  #pragma omp parallel default(shared) private(i,j,k)
  {
  #endif

    // compute updates to Ex
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm; ++k) {
      for (j=1; j < nym; ++j) {
        for (i=0; i < nxm; ++i) {
	  Ex[i + (j + k*ny)*nxm] += Ecoef * ((Hz[i + (j + k*nym)*nxm] - Hz[i + (j-1 + k*nym)*nxm]) \
		- (Hy[i + (j + k*ny)*nxm] - Hy[i + (j + (k-1)*ny)*nxm]));
        }
      }
    }

    // compute updates to Ey
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm; ++k) {
      for (j=0; j < nym; ++j) {
        for (i=1; i < nxm; ++i) {
	  Ey[i + (j + k*nym)*nx] += Ecoef * ((Hx[i + (j + k*nym)*nx] - Hx[i + (j + (k-1)*nym)*nx]) \
		- (Hz[i + (j + k*nym)*nxm] - Hz[i-1 + (j + k*nym)*nxm]));
        }
      }
    }

    // compute updates to Ez
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k < nzm; ++k) {
      for (j=1; j < nym; ++j) {
        for (i=1; i < nxm; ++i) {
	  Ez[i + (j + k*ny)*nx] += Ecoef * ((Hy[i + (j + k*ny)*nxm] - Hy[i-1 + (j + k*ny)*nxm]) \
		- (Hx[i + (j + k*nym)*nx] - Hx[i + (j-1 + k*nym)*nx]));
        }
      }
    }

  #if USE_OPENMP
  } // end parallel region
  #endif

  // stop timer
  t2 = MPI_Wtime();
  step_E_t += t2-t1;
}

/*******************************************************************************
                routine to update interior H fields
*******************************************************************************/
void yee_updater::step_inner_H()
{
  int i,j,k;
  int nxm, nym, nzm;
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  #if USE_OPENMP
  #pragma omp parallel default(shared) private(i,j,k)
  {
  #endif

    // compute updates to Hx
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm-1; ++k) {
      for (j=1; j < nym-1; ++j) {
        for (i=1; i < nxm; ++i) {
	  Hx[i + (j + k*nym)*nx] += Hcoef * ((Ey[i + (j + (k+1)*nym)*nx] - Ey[i + (j + k*nym)*nx]) \
			- (Ez[i + (j+1 + k*ny)*nx] - Ez[i + (j + k*ny)*nx]));
        }
      }
    }

    // compute updates to Hy
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm-1; ++k) {
      for (j=1; j < nym; ++j) {
        for (i=1; i < nxm-1; ++i) {
	  Hy[i + (j + k*ny)*nxm] += Hcoef * ((Ez[i+1 + (j + k*ny)*nx] - Ez[i + (j + k*ny)*nx]) \
			- (Ex[i + (j + (k+1)*ny)*nxm] - Ex[i + (j + k*ny)*nxm]));
        }
      }
    }

    // compute updates to Hz
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm; ++k) {
      for (j=1; j < nym-1; ++j) {
        for (i=1; i < nxm-1; ++i) {
	  Hz[i + (j + k*nym)*nxm] += Hcoef * ((Ex[i + (j+1 + k*ny)*nxm] - Ex[i + (j + k*ny)*nxm]) \
			 - (Ey[i+1 + (j + k*nym)*nx] - Ey[i + (j + k*nym)*nx]));
        }
      }
    }

  #if USE_OPENMP
  }
  #endif

  // stop timer
  t2 = MPI_Wtime();
  step_inner_H_t += t2-t1;
}

/*******************************************************************************
                routine to update boundary H fields
*******************************************************************************/
void yee_updater::step_outer_H()
{
  int i,j,k;
  int nxm, nym, nzm;
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  #if USE_OPENMP
  #pragma omp parallel default(shared) private(i,j,k)
  {
  #endif

    // compute updates to Hx
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k < nzm; ++k) {
      for (j=0; j < nym; ++j) {
        for (i=0; i < nx; i+=nxm) {
	  Hx[i + (j + k*nym)*nx] += Hcoef * ((Ey[i + (j + (k+1)*nym)*nx] - Ey[i + (j + k*nym)*nx]) \
			 - (Ez[i + (j+1 + k*ny)*nx] - Ez[i + (j + k*ny)*nx]));
        }
      }
    }

    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k < nzm; ++k) {
      for (j=0; j < nym; j+=nym-1) {
        for (i=1; i < nxm; ++i) {
	  Hx[i + (j + k*nym)*nx] += Hcoef * ((Ey[i + (j + (k+1)*nym)*nx] - Ey[i + (j + k*nym)*nx]) \
			- (Ez[i + (j+1 + k*ny)*nx] - Ez[i + (j + k*ny)*nx]));
        }
      }
    }

    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k < nzm; k+=nzm-1) {
      for (j=1; j < nym-1; ++j) {
        for (i=1; i < nxm; ++i) {
	  Hx[i + (j + k*nym)*nx] += Hcoef * ((Ey[i + (j + (k+1)*nym)*nx] - Ey[i + (j + k*nym)*nx]) \
			- (Ez[i + (j+1 + k*ny)*nx] - Ez[i + (j + k*ny)*nx]));
        }
      }
    }

    // compute updates to Hy
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k < nzm; ++k) {
      for (j=0; j < ny; j+=nym) {
        for (i=0; i < nxm; ++i) {
	  Hy[i + (j + k*ny)*nxm] += Hcoef * ((Ez[i+1 + (j + k*ny)*nx] - Ez[i + (j + k*ny)*nx]) \
			 - (Ex[i + (j + (k+1)*ny)*nxm] - Ex[i + (j + k*ny)*nxm]));
        }
      }
    }
 
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k < nzm; k+=nzm-1) {
      for (j=1; j < nym; ++j) {
        for (i=0; i < nxm; ++i) {
	  Hy[i + (j + k*ny)*nxm] += Hcoef * ((Ez[i+1 + (j + k*ny)*nx] - Ez[i + (j + k*ny)*nx]) \
			- (Ex[i + (j + (k+1)*ny)*nxm] - Ex[i + (j + k*ny)*nxm]));
        }
      }
    }
    
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm-1; ++k) {
      for (j=1; j < nym; ++j) {
        for (i=0; i < nxm; i+=nxm-1) {
	  Hy[i + (j + k*ny)*nxm] += Hcoef * ((Ez[i+1 + (j + k*ny)*nx] - Ez[i + (j + k*ny)*nx]) \
			 - (Ex[i + (j + (k+1)*ny)*nxm] - Ex[i + (j + k*ny)*nxm]));
        }
      }
    }

    // compute updates to Hz
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k < nz; k+=nzm) {
      for (j=0; j < nym; ++j) {
        for (i=0; i < nxm; ++i) {
	  Hz[i + (j + k*nym)*nxm] += Hcoef * ((Ex[i + (j+1 + k*ny)*nxm] - Ex[i + (j + k*ny)*nxm]) \
			 - (Ey[i+1 + (j + k*nym)*nx] - Ey[i + (j + k*nym)*nx]));
        }
      }
    }

    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm; ++k) {
      for (j=0; j < nym; j+=nym-1) {
        for (i=0; i < nxm; ++i) {
	  Hz[i + (j + k*nym)*nxm] += Hcoef * ((Ex[i + (j+1 + k*ny)*nxm] - Ex[i + (j + k*ny)*nxm]) \
			  - (Ey[i+1 + (j + k*nym)*nx] - Ey[i + (j + k*nym)*nx]));
        }
      }
    }

    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=1; k < nzm; ++k) {
      for (j=1; j < nym-1; ++j) {
        for (i=0; i < nxm; i+=nxm-1) {
	  Hz[i + (j + k*nym)*nxm] += Hcoef * ((Ex[i + (j+1 + k*ny)*nxm] - Ex[i + (j + k*ny)*nxm]) \
			  - (Ey[i+1 + (j + k*nym)*nx] - Ey[i + (j + k*nym)*nx]));
        }
      }
    }
  
  #if USE_OPENMP
  }
  #endif

  // stop timer
  t2 = MPI_Wtime();
  step_outer_H_t += t2-t1;

}


/*******************************************************************************
                Function to update the DAB values
*******************************************************************************/
void yee_updater::step_DAB()
{

  int i,j,k,l;
  int ind[3];
  int nxm, nym;
  int low_ind[3], high_ind[3];
  nxm = nx-1;
  nym = ny-1;
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  if (isBoundaryProc) {

    // loop over the boundary faces
    for (l=0; l<6; ++l) {
      
      // check to see if the current face is of type CRBC
      if (procBounds[l] == crbc::BoundaryProperties::CRBC) {

        // Ex component
        // get the indices the updater object expects as input from this face.
        // Note that these values are inclusive
        bound_upd_Ex.get_input_extents(l, low_ind, high_ind);

        // Because we overlapped the grid the range may extend outside of the 
        // field arrays. To fix this, we simply change -1 -> 0 in the indexing
        // if it occurs.
        for (i=0;i<3; ++i)
          low_ind[i] = (low_ind[i] == -1) ? 0 : low_ind[i];

        // copy in the face values to the Ex faces
        for (k=low_ind[2]; k<=high_ind[2]; ++k) {
          for (j=low_ind[1]; j<=high_ind[1]; ++j) {
            for (i=low_ind[0]; i<=high_ind[0]; ++i) {
              ind[0] = i;
              ind[1] = j;
              ind[2] = k;
              bound_upd_Ex.load_face_data(l, ind, Ex[i + (j + k*ny)*nxm]);
            }
          }
        }
        
        // Ey component
        // get the indices the updater object expects as input from this face.
        // Note that these values are inclusive
        bound_upd_Ey.get_input_extents(l, low_ind, high_ind);

        // Because we overlapped the grid the range may extend outside of the 
        // field arrays. To fix this, we simply change -1 -> 0 in the indexing
        // if it occurs.
        for (i=0;i<3; ++i)
          low_ind[i] = (low_ind[i] == -1) ? 0 : low_ind[i];

        // copy in the face values to the Ex faces
        for (k=low_ind[2]; k<=high_ind[2]; ++k) {
          for (j=low_ind[1]; j<=high_ind[1]; ++j) {
            for (i=low_ind[0]; i<=high_ind[0]; ++i) {
              ind[0] = i;
              ind[1] = j;
              ind[2] = k;
              bound_upd_Ey.load_face_data(l, ind, Ey[i + (j + k*nym)*nx]);
            }
          }
        }

        // Ez component
        // get the indices the updater object expects as input from this face.
        // Note that these values are inclusive
        bound_upd_Ez.get_input_extents(l, low_ind, high_ind);

        // Because we overlapped the grid the range may extend outside of the 
        // field arrays. To fix this, we simply change -1 -> 0 in the indexing
        // if it occurs.
        for (i=0;i<3; ++i)
          low_ind[i] = (low_ind[i] == -1) ? 0 : low_ind[i];

        // copy in the face values to the Ex faces
        for (k=low_ind[2]; k<=high_ind[2]; ++k) {
          for (j=low_ind[1]; j<=high_ind[1]; ++j) {
            for (i=low_ind[0]; i<=high_ind[0]; ++i) {
              ind[0] = i;
              ind[1] = j;
              ind[2] = k;
              bound_upd_Ez.load_face_data(l, ind, Ez[i + (j + k*ny)*nx]);
            }
          }
        }
      } // end if crbc
    } // end for 

    // compute updates
    bound_upd_Ex.compute_updates();
    bound_upd_Ey.compute_updates();
    bound_upd_Ez.compute_updates();

    // now copy the updated values from the DAB back into the fields. We only
    // need to copy the tangential fields because the normal components should
    // already be correct from the Yee updates.

    // loop over the boundary faces
    for (l=0; l<6; ++l) {
      
      // check to see if the current face is of type CRBC
      if (procBounds[l] == crbc::BoundaryProperties::CRBC) {

        // Ex component
        if (l / 2 != 0) { // skip faces with x-normal

          // get the indices the updater object expects to output from this face.
          // Note that these values are inclusive
          bound_upd_Ex.get_output_extents(l, low_ind, high_ind);

          // Because we overlapped the grid the range may extend outside of the 
          // field arrays. To fix this, we simply change -1 -> 0 in the indexing
          // if it occurs.
          for (i=0;i<3; ++i)
            low_ind[i] = (low_ind[i] == -1) ? 0 : low_ind[i];

          // copy in the face values to the Ex faces
          for (k=low_ind[2]; k<=high_ind[2]; ++k) {
            for (j=low_ind[1]; j<=high_ind[1]; ++j) {
              for (i=low_ind[0]; i<=high_ind[0]; ++i) {
                ind[0] = i;
                ind[1] = j;
                ind[2] = k;
                Ex[i + (j + k*ny)*nxm] = bound_upd_Ex.get_new_face_vals(l, ind);
              }
            }
          }
        }
        
        // Ey component
        if (l / 2 != 1) { // skip faces with y-normal

          // get the indices the updater object expects as output from this face.
          // Note that these values are inclusive
          bound_upd_Ey.get_output_extents(l, low_ind, high_ind);

          // Because we overlapped the grid the range may extend outside of the 
          // field arrays. To fix this, we simply change -1 -> 0 in the indexing
          // if it occurs.
          for (i=0;i<3; ++i)
            low_ind[i] = (low_ind[i] == -1) ? 0 : low_ind[i];

          // copy in the face values to the Ex faces
          for (k=low_ind[2]; k<=high_ind[2]; ++k) {
            for (j=low_ind[1]; j<=high_ind[1]; ++j) {
              for (i=low_ind[0]; i<=high_ind[0]; ++i) {
                ind[0] = i;
                ind[1] = j;
                ind[2] = k;
                Ey[i + (j + k*nym)*nx] = bound_upd_Ey.get_new_face_vals(l, ind);
              }
            }
          }
        }

        // Ez component
        if (l / 2 != 2) { // skip faces with z-normal

          // get the indices the updater object expects as output from this face.
          // Note that these values are inclusive
          bound_upd_Ez.get_output_extents(l, low_ind, high_ind);

          // Because we overlapped the grid the range may extend outside of the 
          // field arrays. To fix this, we simply change -1 -> 0 in the indexing
          // if it occurs.
          for (i=0;i<3; ++i)
            low_ind[i] = (low_ind[i] == -1) ? 0 : low_ind[i];

          // copy in the face values to the Ex faces
          for (k=low_ind[2]; k<=high_ind[2]; ++k) {
            for (j=low_ind[1]; j<=high_ind[1]; ++j) {
              for (i=low_ind[0]; i<=high_ind[0]; ++i) {
                ind[0] = i;
                ind[1] = j;
                ind[2] = k;
                Ez[i + (j + k*ny)*nx] = bound_upd_Ez.get_new_face_vals(l, ind);
              }
            }
          }
        }
      } // end if crbc
    } // end for 
  } // end isBoundaryProc

  // stop timer
  t2 = MPI_Wtime();
  step_DAB_t += t2-t1;
}


/*******************************************************************************
                Function to get DAB values 
*******************************************************************************/
void yee_updater::get_dab_vals_loop(std::vector<double> &buffer,
                         crbc::CrbcUpdates<3, double, int> &updater,
                         const int &side,
                         const int *low,
                         const int *high,
                         const int *plow,
                         const int *phigh,
                         const bool &isedge) 
{

  int i,j,k,p,q, ind[3];

  if (isedge) {
    int pind[2];
    for (p=plow[0]; p<=phigh[0]; ++p) {
	pind[0] = p;
	for (q=plow[1]; q<=phigh[1]; ++q) {
	  pind[1] = q;
	  for (k=low[2]; k<=high[2]; ++k) {
	    ind[2] = k;
	    for (j=low[1]; j<=high[1]; ++j) {
	      ind[1] = j;
	      for (i=low[0]; i<=high[0]; ++i) {
		ind[0] = i;
		buffer.push_back(updater.get_edge_auxiliary_vars(side, ind, pind));
	      } // i
	    } // j
	  } // k
	} // q
      } // p
    } else {
      for (p=plow[0]; p<=phigh[0]; ++p) {
	for (k=low[2]; k<=high[2]; ++k) {
	  ind[2] = k;
	  for (j=low[1]; j<=high[1]; ++j) {
	    ind[1] = j;
	    for (i=low[0]; i<=high[0]; ++i) {
	      ind[0] = i;
	      buffer.push_back(updater.get_auxiliary_vars(side, ind, p));
	    } // i
	  } // j
	} // k
      } // p
    }
  }


/*******************************************************************************
                Function to set DAB values 
*******************************************************************************/
void yee_updater::set_dab_vals_loop(std::vector<double> &buffer,
                         crbc::CrbcUpdates<3, double, int> &updater,
                         int &count,
                         const int &side,
                         const int *low,
                         const int *high,
                         const int *plow,
                         const int *phigh,
                         const bool &isedge) 
{
    int i,j,k,p,q, ind[3];

    if (isedge) {
      int pind[2];
      for (p=plow[0]; p<=phigh[0]; ++p) {
	pind[0] = p;
	for (q=plow[1]; q<=phigh[1]; ++q) {
	  pind[1] = q;
	  for (k=low[2]; k<=high[2]; ++k) {
	    ind[2] = k;
	    for (j=low[1]; j<=high[1]; ++j) {
	      ind[1] = j;
	      for (i=low[0]; i<=high[0]; ++i) {
		ind[0] = i;
		updater.set_edge_auxiliary_vars(side, ind, pind, buffer[count++]);

                
	      } // i
	    } // j
	  } // k
	} // q
      } // p
    } else {
      for (p=plow[0]; p<=phigh[0]; ++p) {
	for (k=low[2]; k<=high[2]; ++k) {
	  ind[2] = k;
	  for (j=low[1]; j<=high[1]; ++j) {
	    ind[1] = j;
	    for (i=low[0]; i<=high[0]; ++i) {
	      ind[0] = i;
	      updater.set_auxiliary_vars(side, ind, p, buffer[count++]);
	    } // i
	  } // j
	} // k
      } // p
    }
  }


/*******************************************************************************
                Function identify what needs to be sent of the DAB
*******************************************************************************/
void yee_updater::calc_DAB_send_params()
{

  // note this only checks the cases possible in this implementation. A more 
  // generic implementation would need to handle more cases. In particular, we
  // can't have 2 parallel faces need to pass data from the same process since
  // we are assuming each direction has the same number of processes abd we do 
  // not need to pass any information if we only have 1 process.

  unsigned int l, m, tang_sides[4];

  // first identify the directions that we need to send
  for (l=0; l<6; ++l) { // loop over sides
    if (procBounds[l] == crbc::BoundaryProperties::CRBC) {
      for (m=0; m<6; ++m) { // loop over sides
        if (l/2 == m/2)
          continue; // skip parallel sides
        if (procBounds[m] == crbc::BoundaryProperties::NONE)
          send_dirs.push_back(m);
      }
    } 
  }

  // remove any duplicates in send_dirs
  std::sort(send_dirs.begin(), send_dirs.end());
  auto last = std::unique(send_dirs.begin(), send_dirs.end());
  send_dirs.erase(last, send_dirs.end());

  // loop over the send directions
  for (l=0; l<send_dirs.size(); ++l) {

    // figure out which sides need to send in the current direction
    // start by listing the tangential sides
    switch (send_dirs[l] / 2) {

      case 0: // outward normal is +/-x
        tang_sides[0] = 2; // left y
        tang_sides[1] = 3; // right y
        tang_sides[2] = 4; // left z
        tang_sides[3] = 5; // right z
        break;
      case 1: // outward normal is +/-y
        tang_sides[0] = 0; // left x
        tang_sides[1] = 1; // right x
        tang_sides[2] = 4; // left z
        tang_sides[3] = 5; // right z
        break;
      case 2: // outward normal is +/-z
        tang_sides[0] = 0; // left x
        tang_sides[1] = 1; // right x
        tang_sides[2] = 2; // left y
        tang_sides[3] = 3; // left z
        break;
      default: // shouldn't happen
        for (m=0; m<4; ++m)
          tang_sides[m] = -1;
        std::cerr << "invalid side" << std::endl;
        break;
    }

    // check to see if any of the tangetial sides are a DAB/CRBC boundary
    for (m=0; m<4; ++m) {
      if (procBounds[tang_sides[m]] == crbc::BoundaryProperties::CRBC) {
        send_sides[l].push_back(tang_sides[m]);
      }
    }
  }

  // finally change the directions from local side indices to the MPI ranks in 
  // the appropriate direction
  for (l=0; l<send_dirs.size(); ++l) {
    switch (send_dirs[l]) {
            
      case 0: // WEST
        send_mpi_dirs.push_back(WEST);
        break;
      case 1: // EAST
        send_mpi_dirs.push_back(EAST);
        break;
      case 2: // SOUTH
        send_mpi_dirs.push_back(SOUTH);
        break;
      case 3: // NORTH
        send_mpi_dirs.push_back(NORTH);
        break;
      case 4: // DOWN
        send_mpi_dirs.push_back(DOWN);
        break;
      case 5: // UP
        send_mpi_dirs.push_back(UP);
        break;
    }
  }
}

/*******************************************************************************
                Function to send DAB values 
*******************************************************************************/
void yee_updater::send_DAB()
{

  // notice that we could make this much shorter if we declared the DAB updaters
  // as an array instead of individually for each component. We hope this 
  // presentation is slightly clearer.

  unsigned int l,m, side, sidea;
  int edge;
  int low[3], high[3], plow[2], phigh[2];
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  if (isBoundaryProc) {

    // loop over the directions we need to send
    for (l=0; l<send_dirs.size(); ++l) {

      side = send_dirs[l];

      // clear the buffers for this direction
      DAB_sbuf[side].clear();
      DAB_rbuf[side].clear();
      
      // loop over the sides we need to send in the current direction
      for (m=0; m<send_sides[l].size(); ++m) {

        sidea = send_sides[l].at(m);

        // get data extents for Ex
        // This gets us the indices for all the points on the plane
        // parallel to the phyiscal boundary
        bound_upd_Ex.get_output_extents(sidea, low, high);

        // now we need to restrict to second to last line of points parallel to 
        // the boundary we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = ++low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = --high[side / 2];       
        }

        // copy data auxilliary data into the send buffer
        plow[0] = 0; // auxilliary index bounds
        phigh[0] = bound_upd_Ex.get_num_recursions(sidea);
        get_dab_vals_loop(DAB_sbuf[side], bound_upd_Ex, sidea, low, high, plow, phigh);

        // get data extents for Ey
        // This gets us the indices for all the points on the plane
        // parallel to the phyiscal boundary
        bound_upd_Ey.get_output_extents(sidea, low, high);

        // now we need to restrict to second to last line of points parallel to 
        // the boundary we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = ++low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = --high[side / 2];       
        }

        // copy data auxilliary data into the send buffer
        plow[0] = 0; // auxilliary index bounds
        phigh[0] = bound_upd_Ey.get_num_recursions(sidea);
        get_dab_vals_loop(DAB_sbuf[side], bound_upd_Ey, sidea, low, high, plow, phigh);

        // get data extents for Ez
        // This gets us the indices for all the points on the plane
        // parallel to the phyiscal boundary
        bound_upd_Ez.get_output_extents(sidea, low, high);

        // now we need to restrict to second to last line of points parallel to 
        // the boundary we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = ++low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = --high[side / 2];       
        }

        // copy data auxilliary data into the send buffer
        plow[0] = 0; // auxilliary index bounds
        phigh[0] = bound_upd_Ez.get_num_recursions(sidea);
        get_dab_vals_loop(DAB_sbuf[side], bound_upd_Ez, sidea, low, high, plow, phigh);

      } // end loop over sides in the current direction
      // now send any edge data
      if (send_sides[l].size() == 2) {

        // Ex
        // get edge index
        edge = bound_upd_Ex.get_edge_index(send_sides[l][0], send_sides[l][1]);

        // get edge data extents
        bound_upd_Ex.get_edge_extents(edge, low, high, plow, phigh);

        // now we need to restrict to second to last point parallel to 
        // the boundaries we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
        high[side / 2] = ++low[side / 2];
        } else { // right side in the appropriate direction
        low[side / 2] = --high[side / 2];       
        }

        // the true at the end tells the function plow, phigh are arrays of len 2
        get_dab_vals_loop(DAB_sbuf[side], bound_upd_Ex, edge, low, high, plow, phigh, true);

        // Ey
        // get edge index
        edge = bound_upd_Ey.get_edge_index(send_sides[l][0], send_sides[l][1]);

        // get edge data extents
        bound_upd_Ey.get_edge_extents(edge, low, high, plow, phigh);

        // now we need to restrict to second to last point parallel to 
        // the boundaries we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = ++low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = --high[side / 2];       
        }
        // the true at the end tells the function plow, phigh are arrays of len 2
        get_dab_vals_loop(DAB_sbuf[side], bound_upd_Ey, edge, low, high, plow, phigh, true);

        // Ez
        // get edge index
        edge = bound_upd_Ez.get_edge_index(send_sides[l][0], send_sides[l][1]);

        // get edge data extents
        bound_upd_Ez.get_edge_extents(edge, low, high, plow, phigh);

        // now we need to restrict to second to last point parallel to 
        // the boundaries we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = ++low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = --high[side / 2];       
        }
        // the true at the end tells the function plow, phigh are arrays of len 2
        get_dab_vals_loop(DAB_sbuf[side], bound_upd_Ez, edge, low, high, plow, phigh, true);

      } // end if 2 sides

      // now actually send the data. We use non-blocking sends in the hope that
      // the time it takes to compute the H-field updates will mask most or all
      // of the communication

      // create and save requests
      MPI_Request sreq, rreq;

      DAB_send_req.push_back(sreq);
      DAB_recv_req.push_back(rreq);

      // Send --- we use a tag of 1 for DAB send and 0 for E-field sends
      if (MPI_Isend(DAB_sbuf[side].data(), DAB_sbuf[side].size(), MPI_DOUBLE, \
          send_mpi_dirs[l], 1, grid_comm, &DAB_send_req.back()) != MPI_SUCCESS)
                std::cerr << "MPI_Isend failed" << std::endl;

      // make sure the recieve buffer is large enough
      if (DAB_rbuf[side].size() < DAB_sbuf[side].size())
        DAB_rbuf[side].assign(DAB_sbuf[side].size(), 0.0);

      // Recieve
      if (MPI_Irecv(DAB_rbuf[side].data(), DAB_rbuf[side].size(), MPI_DOUBLE, 
          send_mpi_dirs[l], 1, grid_comm, &DAB_recv_req.back()) != MPI_SUCCESS)
                std::cerr << "MPI_Isend failed" << std::endl;

    } // end loop over send directions

  } // end if boundary proc
  
  // stop timer
  t2 = MPI_Wtime();
  send_DAB_t += t2-t1;
} // end send_DAB

/*******************************************************************************
                Function to receive DAB values 
*******************************************************************************/
void yee_updater::recv_DAB()
{
  
  // note the recieve directions are the same as the send directions, so we do
  // almost exactly the same thing here as in send_DAB() except we copy from the
  // buffers to the DAB.
  unsigned int l,m, side, edge, sidea;
  int count;
  int low[3], high[3], plow[2], phigh[2];
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  // make sure all of the recieves are complete before we use the data
  if (MPI_Waitall(DAB_recv_req.size(), DAB_recv_req.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    std::cerr << "MPI_Waitall failed (DAB_recv_req)" << std::endl;

  DAB_recv_req.clear();

  // copy the values from the buffers
  if (isBoundaryProc) {

    // loop over the directions we need to send
    for (l=0; l<send_dirs.size(); ++l) {

      // used for indexing
      count = 0;

      side = send_dirs[l];

      // loop over the sides we need to send in the current direction
      for (m=0; m<send_sides[l].size(); ++m) {

        sidea = send_sides[l].at(m);

        // get data extents for Ex
        // This gets us the indices for all the points on the plane
        // parallel to the phyiscal boundary
        bound_upd_Ex.get_output_extents(sidea, low, high);

        // now we need to restrict to  last line of points parallel to 
        // the boundary we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = high[side / 2];       
        }

        // copy data auxilliary data into the send buffer
        plow[0] = 0; // auxilliary index bounds
        phigh[0] = bound_upd_Ex.get_num_recursions(sidea);
        set_dab_vals_loop(DAB_rbuf[side], bound_upd_Ex, count, sidea, low, high, plow, phigh);

        // get data extents for Ey
        // This gets us the indices for all the points on the plane
        // parallel to the phyiscal boundary
        bound_upd_Ey.get_output_extents(sidea, low, high);

        // now we need to restrict to last line of points parallel to 
        // the boundary we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = high[side / 2];       
        }

        // copy data auxilliary data into the send buffer
        plow[0] = 0; // auxilliary index bounds
        phigh[0] = bound_upd_Ey.get_num_recursions(sidea);
        set_dab_vals_loop(DAB_rbuf[side], bound_upd_Ey, count, sidea, low, high, plow, phigh);

        // get data extents for Ez
        // This gets us the indices for all the points on the plane
        // parallel to the phyiscal boundary
        bound_upd_Ez.get_output_extents(sidea, low, high);

        // now we need to restrict to last line of points parallel to 
        // the boundary we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = high[side / 2];       
        }

        // copy data auxilliary data into the send buffer
        plow[0] = 0; // auxilliary index bounds
        phigh[0] = bound_upd_Ez.get_num_recursions(sidea);
        set_dab_vals_loop(DAB_rbuf[side], bound_upd_Ez, count, sidea, low, high, plow, phigh);

      } // end loop over sides in the current direction

      // now receive any edge data
      if (send_sides[l].size() == 2) {

        // Ex
        // get edge index
        edge = bound_upd_Ex.get_edge_index(send_sides[l][0], send_sides[l][1]);
        // get edge data extents
        bound_upd_Ex.get_edge_extents(edge, low, high, plow, phigh);

        // now we need to restrict to last point parallel to 
        // the boundaries we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = high[side / 2];       
        }
        // the true at the end tells the function plow, phigh are arrays of len 2
        set_dab_vals_loop(DAB_rbuf[side], bound_upd_Ex, count, edge, low, high, plow, phigh, true);

        // Ey
        // get edge index
        edge = bound_upd_Ey.get_edge_index(send_sides[l][0], send_sides[l][1]);
        // get edge data extents
        bound_upd_Ey.get_edge_extents(edge, low, high, plow, phigh);

        // now we need to restrict to last point parallel to 
        // the boundaries we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = high[side / 2];       
        }
        // the true at the end tells the function plow, phigh are arrays of len 2
        set_dab_vals_loop(DAB_rbuf[side], bound_upd_Ey, count, edge, low, high, plow, phigh, true);

        // Ez
        // get edge index
        edge = bound_upd_Ez.get_edge_index(send_sides[l][0], send_sides[l][1]);
        // get edge data extents
        bound_upd_Ez.get_edge_extents(edge, low, high, plow, phigh);

        // now we need to restrict to last point parallel to 
        // the boundaries we are sending. side / 2 gives us the component of the
        // extents we need to modify and side % 2 tells us if we need to modify
        // the low or high extents
        if (side % 2 == 0) { // left side in the appropriate direction
          high[side / 2] = low[side / 2];
        } else { // right side in the appropriate direction
          low[side / 2] = high[side / 2];       
        }
        // the true at the end tells the function plow, phigh are arrays of len 2
        set_dab_vals_loop(DAB_rbuf[side], bound_upd_Ez, count, edge, low, high, plow, phigh, true);

      } // end if 2 sides

    } // end loop over send directions

  } // end if boundary proc

  // make sure all the sends have completed
  if (MPI_Waitall(DAB_send_req.size(), DAB_send_req.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS)
      std::cerr << "MPI_Waitall failed (DAB_send_req)" << std::endl;

  DAB_send_req.clear();
  
  // stop timer
  t2 = MPI_Wtime();
  recv_DAB_t += t2-t1;

} // end receive DAB


/*******************************************************************************
                Function to send E field values 
*******************************************************************************/
void yee_updater::send_E()
{
  // we just handle each case explicitly, but note that it wouldn't be terribly
  // difficult to greatly reduce the code here
 
  int i,j,k;
  int ns = maxn*maxn;
  int nxm, nym, nzm;
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;
  double t1, t2;

  // start timer
  t1 = MPI_Wtime();

  // send the last plane of points we were able to update for the tangential 
  // components
  if (procBounds[0] == crbc::BoundaryProperties::NONE) { // left side in x, send Ey, Ez

    MPI_Request sreq, rreq;

    send_req.push_back(sreq);
    recv_req.push_back(rreq);

    // copy Ey and Ez to buffer
    i = 1; // grid overlap
    for (k = 0; k<nz; ++k)
      for (j=0; j<nym; ++j)
        E_sbuf[0][j + k*nym] = Ey[i + (j + k*nym)*nx];

      for (k = 0; k<nzm; ++k)
        for (j=0; j<ny; ++j)
          E_sbuf[0][ns + j + k*ny] = Ez[i + (j + k*ny)*nx];

      if (MPI_Isend(E_sbuf[0].data(), 2*ns, MPI_DOUBLE, WEST, 0, grid_comm, &send_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

      if (MPI_Irecv(E_rbuf[0].data(), 2*ns, MPI_DOUBLE, WEST, 0, grid_comm, &recv_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

    }

    if (procBounds[1] == crbc::BoundaryProperties::NONE) { // right side in x, send Ey, Ez

      MPI_Request sreq, rreq;

      send_req.push_back(sreq);
      recv_req.push_back(rreq);

      // copy Ey and Ez to buffer
      i = nx-2; // overlap
      for (k = 0; k<nz; ++k)
        for (j=0; j<nym; ++j)
          E_sbuf[1][j + k*nym] = Ey[i + (j + k*nym)*nx];

      for (k = 0; k<nzm; ++k)
        for (j=0; j<ny; ++j)
          E_sbuf[1][ns + j + k*ny] = Ez[i + (j + k*ny)*nx];

      if (MPI_Isend(E_sbuf[1].data(), 2*ns, MPI_DOUBLE, EAST, 0, grid_comm, &send_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

      if (MPI_Irecv(E_rbuf[1].data(), 2*ns, MPI_DOUBLE, EAST, 0, grid_comm, &recv_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

    }

    if (procBounds[2] == crbc::BoundaryProperties::NONE) { // left side in y, send Ex, Ez

      MPI_Request sreq, rreq;

      send_req.push_back(sreq);
      recv_req.push_back(rreq);

      // copy Ex and Ez to buffer
      j = 1; // overlap
      for (k = 0; k<nz; ++k)
        for (i=0; i<nxm; ++i)
          E_sbuf[2][i + k*nxm] = Ex[i + (j + k*ny)*nxm];

      for (k = 0; k<nzm; ++k)
        for (i=0; i<nx; ++i)
          E_sbuf[2][ns + i + k*nx] = Ez[i + (j + k*ny)*nx];

      if (MPI_Isend(E_sbuf[2].data(), 2*ns, MPI_DOUBLE, SOUTH, 0, grid_comm, &send_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

      if (MPI_Irecv(E_rbuf[2].data(), 2*ns, MPI_DOUBLE, SOUTH, 0, grid_comm, &recv_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

    }

    if (procBounds[3] == crbc::BoundaryProperties::NONE) { // right side in y, send Ex, Ez

      MPI_Request sreq, rreq;

      send_req.push_back(sreq);
      recv_req.push_back(rreq);

      // copy Ex and Ez to buffer
      j = ny-2; // overlap
      for (k = 0; k<nz; ++k)
        for (i=0; i<nxm; ++i)
          E_sbuf[3][i + k*nxm] = Ex[i + (j + k*ny)*nxm];

      for (k = 0; k<nzm; ++k)
        for (i=0; i<nx; ++i)
          E_sbuf[3][ns + i + k*nx] = Ez[i + (j + k*ny)*nx];

      if (MPI_Isend(E_sbuf[3].data(), 2*ns, MPI_DOUBLE, NORTH, 0, grid_comm, &send_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

      if (MPI_Irecv(E_rbuf[3].data(), 2*ns, MPI_DOUBLE, NORTH, 0, grid_comm, &recv_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

    }

    if (procBounds[4] == crbc::BoundaryProperties::NONE) { // left side in z, send Ex, Ey

      MPI_Request sreq, rreq;

      send_req.push_back(sreq);
      recv_req.push_back(rreq);

      // copy Ex and Ey to buffer
      k = 1; // overlap
      for (j = 0; j<ny; ++j)
        for (i=0; i<nxm; ++i)
          E_sbuf[4][i + j*nxm] = Ex[i + (j + k*ny)*nxm];

      for (j = 0; j<nym; ++j)
        for (i=0; i<nx; ++i)
          E_sbuf[4][ns + i + j*nx] = Ey[i + (j + k*nym)*nx];

      if (MPI_Isend(E_sbuf[4].data(), 2*ns, MPI_DOUBLE, DOWN, 0, grid_comm, &send_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

      if (MPI_Irecv(E_rbuf[4].data(), 2*ns, MPI_DOUBLE, DOWN, 0, grid_comm, &recv_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

    }

    if (procBounds[5] == crbc::BoundaryProperties::NONE) { // right side in z, send Ex, Ey

      MPI_Request sreq, rreq;

      send_req.push_back(sreq);
      recv_req.push_back(rreq);

      // copy Ex and Ey to buffer
      k = nz-2; // overlap
      for (j = 0; j<ny; ++j)
        for (i=0; i<nxm; ++i)
          E_sbuf[5][i + j*nxm] = Ex[i + (j + k*ny)*nxm];

      for (j = 0; j<nym; ++j)
        for (i=0; i<nx; ++i)
          E_sbuf[5][ns + i + j*nx] = Ey[i + (j + k*nym)*nx];

      if (MPI_Isend(E_sbuf[5].data(), 2*ns, MPI_DOUBLE, UP, 0, grid_comm, &send_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

      if (MPI_Irecv(E_rbuf[5].data(), 2*ns, MPI_DOUBLE, UP, 0, grid_comm, &recv_req.back()) != MPI_SUCCESS)
        std::cerr << "MPI_Isend failed" << std::endl;

    }

  // stop timer
  t2 = MPI_Wtime();
  send_E_t += t2-t1;

}

/*******************************************************************************
                Function to receive E field values 
*******************************************************************************/
void yee_updater::recv_E()
{
  // again we just handle all the cases explicitly

  int i,j,k;
  int ns = maxn*maxn;
  int nxm, nym, nzm;
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;
  double t1, t2;
  t1 = MPI_Wtime();

  // make sure all of the recieves are complete
  if (MPI_Waitall(recv_req.size(), recv_req.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    std::cerr << "MPI_Waitall failed (recv_req)" << std::endl;

  // copy the values from the buffers
  recv_req.clear();

  if (procBounds[0] == crbc::BoundaryProperties::NONE) { // left side in x, recieve Ey, Ez

    // copy Ey and Ez from buffer
    i = 0;
    for (k = 0; k<nz; ++k)
      for (j=0; j<nym; ++j)
        Ey[i + (j + k*nym)*nx] = E_rbuf[0][j + k*nym];

    for (k = 0; k<nzm; ++k)
      for (j=0; j<ny; ++j)
        Ez[i + (j + k*ny)*nx] = E_rbuf[0][ns + j + k*ny];

  }

  if (procBounds[1] == crbc::BoundaryProperties::NONE) { // right side in x, recieve Ey, Ez

    // copy Ey and Ez from buffer
    i = nx-1;
    for (k = 0; k<nz; ++k)
      for (j=0; j<nym; ++j)
        Ey[i + (j + k*nym)*nx] = E_rbuf[1][j + k*nym];

    for (k = 0; k<nzm; ++k)
      for (j=0; j<ny; ++j)
        Ez[i + (j + k*ny)*nx] = E_rbuf[1][ns + j + k*ny];

  }

  if (procBounds[2] == crbc::BoundaryProperties::NONE) { // left side in y, recieve Ex, Ez

    // copy Ex and Ez from buffer
    j = 0; 
    for (k = 0; k<nz; ++k)
      for (i=0; i<nxm; ++i)
        Ex[i + (j + k*ny)*nxm] = E_rbuf[2][i + k*nxm];

    for (k = 0; k<nzm; ++k)
      for (i=0; i<nx; ++i)
        Ez[i + (j + k*ny)*nx] = E_rbuf[2][ns + i + k*nx];

  }

  if (procBounds[3] == crbc::BoundaryProperties::NONE) { // right side in y, recieve Ex, Ez

    // copy Ex and Ez from buffer
    j = ny-1;
    for (k = 0; k<nz; ++k)
      for (i=0; i<nxm; ++i)
        Ex[i + (j + k*ny)*nxm] = E_rbuf[3][i + k*nxm];

    for (k = 0; k<nzm; ++k)
      for (i=0; i<nx; ++i)
        Ez[i + (j + k*ny)*nx] = E_rbuf[3][ns + i + k*nx];

  }

  if (procBounds[4] == crbc::BoundaryProperties::NONE) { // left side in z, recieve Ex, Ey

    // copy Ex and Ey from buffer
    k = 0;
    for (j = 0; j<ny; ++j)
      for (i=0; i<nxm; ++i)
        Ex[i + (j + k*ny)*nxm] = E_rbuf[4][i + j*nxm];

    for (j = 0; j<nym; ++j)
      for (i=0; i<nx; ++i)
        Ey[i + (j + k*nym)*nx] = E_rbuf[4][ns + i + j*nx];

  }

  if (procBounds[5] == crbc::BoundaryProperties::NONE) { // right side in z, recieve Ex, Ey

    // copy Ex and Ey from buffer
    k = nz-1;
    for (j = 0; j<ny; ++j)
      for (i=0; i<nxm; ++i)
        Ex[i + (j + k*ny)*nxm] = E_rbuf[5][i + j*nxm];

      for (j = 0; j<nym; ++j)
        for (i=0; i<nx; ++i)
          Ey[i + (j + k*nym)*nx] = E_rbuf[5][ns + i + j*nx];

  }

  // make sure all the sends have completed
  if (MPI_Waitall(send_req.size(), send_req.data(), MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    std::cerr << "MPI_Waitall failed (send_req)" << std::endl;

  send_req.clear();

  t2 = MPI_Wtime();
  recv_E_t += t2-t1;
}


/*******************************************************************************
                Function to load intitial values
*******************************************************************************/
void yee_updater::load_initial_conds() {

  int i,j,k, nxm, nym, nzm;
  double x[3];
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;
  double t1, t2;
  t1 = MPI_Wtime();

  // set the solution routine to use the Yee schemes FD operator to compute the
  // curl in the solution. We do this to get a numerically div free initial 
  // condition
  sol_obj.set_derivative_method(maxwell_solutions::FD_YEE);

  #if USE_OPENMP
  #pragma omp parallel default(shared) private(i,j,k,x,sol_obj)
  {
  #endif

    // load Ex values
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k<nz; ++k) {
      for (j=0; j<ny; ++j) {
	for (i=0; i<nxm; ++i) {
          x[2] = coord[2] + h*k;
          x[1] = coord[1] + h*j;
	  x[0] = coord[0] + h*i + h/2.0;
	  Ex[i + (j + k*ny)*nxm] = sol_obj.get_Ex_solution(x, Etime);
	}
      }
    }

    // load Ey values
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k<nz; ++k) {
      for (j=0; j<nym; ++j) {
        for (i=0; i<nx; ++i) {
          x[2] = coord[2] + h*k;
          x[1] = coord[1] + h*j + h/2.0;
          x[0] = coord[0] + h*i;
          Ey[i + (j + k*nym)*nx] = sol_obj.get_Ey_solution(x, Etime);
        }
      }
    }

    // load Ez values
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k<nzm; ++k) {
      for (j=0; j<ny; ++j) {
        for (i=0; i<nx; ++i) {
          x[2] = coord[2] + h*k + h/2.0;
          x[1] = coord[1] + h*j;
          x[0] = coord[0] + h*i;
          Ez[i + (j + k*ny)*nx] = sol_obj.get_Ez_solution(x, Etime);
        }
      }
    }

    // load Hx values
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k<nzm; ++k) {
      for (j=0; j<nym; ++j) {
        for (i=0; i<nx; ++i) {
          x[2] = coord[2] + h*k+ h/2.0;
          x[1] = coord[1] + h*j+ h/2.0;
          x[0] = coord[0] + h*i;
          Hx[i + (j + k*nym)*nx] = sol_obj.get_Hx_solution(x, Htime);
        }
      }
    }

    // load Hy values
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k<nzm; ++k) {
      for (j=0; j<ny; ++j) {
        for (i=0; i<nxm; ++i) {
          x[2] = coord[2] + h*k+ h/2.0;
          x[1] = coord[1] + h*j;
          x[0] = coord[0] + h*i+ h/2.0;
          Hy[i + (j + k*ny)*nxm] = sol_obj.get_Hy_solution(x, Htime);
        }
      }
    }

    // load Hz values
    #if USE_OPENMP
    #pragma omp for collapse(3)
    #endif
    for (k=0; k<nz; ++k) {
      for (j=0; j<nym; ++j) {
        for (i=0; i<nxm; ++i) {
          x[2] = coord[2] + h*k;
          x[1] = coord[1] + h*j+ h/2.0;
          x[0] = coord[0] + h*i+ h/2.0;
	  Hz[i + (j + k*nym)*nxm] = sol_obj.get_Hz_solution(x, Htime);
	}
      }
    }

  #if USE_OPENMP
  }
  #endif

  t2 = MPI_Wtime();
  load_init_conds_t += t2-t1;
}


/*******************************************************************************
                Function to calculate norm (squared)
*******************************************************************************/
double yee_updater::calc_norm()
{
  int nxm, nym, nzm;
  double norm = 0.0, temp = 0.0;
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;

  double t1, t2;
  int i, j, k;
  t1 = MPI_Wtime();

  #if USE_OPENMP
  #pragma omp parallel default(shared) private(i,j,k)
  {
  #endif

    // load Ex values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nz; k+=skip) {
      for (j=0; j<ny; j+=skip) {
	for (i=0; i<nxm; i+=skip) {
          temp += Ex[i + (j + k*ny)*nxm] * Ex[i + (j + k*ny)*nxm];
	}
      }
    }

    // load Ey values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nz; k+=skip) {
      for (j=0; j<nym; j+=skip) {
        for (i=0; i<nx; i+=skip) {
          temp += Ey[i + (j + k*nym)*nx] * Ey[i + (j + k*nym)*nx];
        }
      }
    }

    // load Ez values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nzm; k+=skip) {
      for (j=0; j<ny; j+=skip) {
        for (i=0; i<nx; i+=skip) {
          temp += Ez[i + (j + k*ny)*nx] * Ez[i + (j + k*ny)*nx];
        }
      }
    }

    #if USE_OPENMP
    #pragma omp barrier
    #pragma omp single
    #endif
    norm = eps*temp;

    temp = 0.0;
    #if USE_OPENMP
    #pragma omp barrier
    #endif

    // load Hx values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nzm; k+=skip) {
      for (j=0; j<nym; j+=skip) {
        for (i=0; i<nx; i+=skip) {
          temp += Hx[i + (j + k*nym)*nx] * Hx[i + (j + k*nym)*nx];
        }
      }
    }

    // load Hy values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nzm; k+=skip) {
      for (j=0; j<ny; j+=skip) {
        for (i=0; i<nxm; i+=skip) {
          temp += Hy[i + (j + k*ny)*nxm] * Hy[i + (j + k*ny)*nxm];
        }
      }
    }

    // load Hz values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nz; k+=skip) {
      for (j=0; j<nym; j+=skip) {
        for (i=0; i<nxm; i+=skip) {
          temp += Hz[i + (j + k*nym)*nxm] * Hz[i + (j + k*nym)*nxm];
	}
      }
    }

  #if USE_OPENMP
  }
  #endif  
  
  norm += temp*mu;

  t2 = MPI_Wtime();
  calc_norm_t += t2-t1;

  return norm;
}

/*******************************************************************************
                Function to calculate error (squared)
*******************************************************************************/
double yee_updater::calc_error()
{
  int nxm, nym, nzm;
  double error = 0.0, temp = 0.0, sol, x[3];
  nxm = nx-1;
  nym = ny-1;
  nzm = nz-1;

  double t1, t2;
  int i, j, k;
  t1 = MPI_Wtime();

  // set the solution routine to use the exact expression for the curl
  sol_obj.set_derivative_method(maxwell_solutions::EXACT);

  #if USE_OPENMP
  #pragma omp parallel default(shared) private(i,j,k,x,sol,upd_obj)
  {
  #endif

    // load Ex values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nz; k+=skip) {
      for (j=0; j<ny; j+=skip) {
	for (i=0; i<nxm; i+=skip) {
          x[2] = coord[2] + h*k;
          x[1] = coord[1] + h*j;
	  x[0] = coord[0] + h*i + h/2.0;
          sol = sol_obj.get_Ex_solution(x, Etime) - Ex[i + (j + k*ny)*nxm];
          temp +=  sol * sol;
	}
      }
    }

    // load Ey values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nz; k+=skip) {
      for (j=0; j<nym; j+=skip) {
        for (i=0; i<nx; i+=skip) {
          x[2] = coord[2] + h*k;
          x[1] = coord[1] + h*j + h/2.0;
	  x[0] = coord[0] + h*i;
          sol = sol_obj.get_Ey_solution(x, Etime) - Ey[i + (j + k*nym)*nx];
          temp += sol * sol;
        }
      }
    }

    // load Ez values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nzm; k+=skip) {
      for (j=0; j<ny; j+=skip) {
        for (i=0; i<nx; i+=skip) {
          x[2] = coord[2] + h*k + h/2.0;
          x[1] = coord[1] + h*j;
	  x[0] = coord[0] + h*i;
          sol = sol_obj.get_Ez_solution(x, Etime) - Ez[i + (j + k*ny)*nx];
          temp += sol * sol;
        }
      }
    }

    #if USE_OPENMP
    #pragma omp barrier
    #pragma omp single
    #endif
    error = eps*temp;

    temp = 0.0;
    #if USE_OPENMP
    #pragma omp barrier
    #endif

    // load Hx values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nzm; k+=skip) {
      for (j=0; j<nym; j+=skip) {
        for (i=0; i<nx; i+=skip) {
          x[2] = coord[2] + h*k + h/2.0;
          x[1] = coord[1] + h*j + h/2.0;
	  x[0] = coord[0] + h*i;
          sol = sol_obj.get_Hx_solution(x, Htime) - Hx[i + (j + k*nym)*nx];
          temp += sol * sol;
        }
      }
    }

    // load Hy values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nzm; k+=skip) {
      for (j=0; j<ny; j+=skip) {
        for (i=0; i<nxm; i+=skip) {
          x[2] = coord[2] + h*k + h/2.0;
          x[1] = coord[1] + h*j;
	  x[0] = coord[0] + h*i + h/2.0;
          sol = sol_obj.get_Hy_solution(x, Htime) - Hy[i + (j + k*ny)*nxm];
          temp += sol * sol;
        }
      }
    }

    // load Hz values
    #if USE_OPENMP
    #pragma omp for reduction(+:temp) collapse(3)
    #endif
    for (k=0; k<nz; k+=skip) {
      for (j=0; j<nym; j+=skip) {
        for (i=0; i<nxm; i+=skip) {
          x[2] = coord[2] + h*k;
          x[1] = coord[1] + h*j + h/2.0;
	  x[0] = coord[0] + h*i + h/2.0;
          sol = sol_obj.get_Hz_solution(x, Htime) - Hz[i + (j + k*nym)*nxm];
          temp += sol * sol;
	}
      }
    }

  #if USE_OPENMP
  }
  #endif  
  
  error += temp*mu;

  t2 = MPI_Wtime();
  calc_err_t += t2-t1;

  return error;
}


