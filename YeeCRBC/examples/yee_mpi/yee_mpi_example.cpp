#include "yee_mpi.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

int main(int argc, char* argv[]) {

  int nprocs, id, procs_per_dim, n, p;
  double t, crbc_t, w, t1, t2;
  MPI_Comm comm = MPI_COMM_WORLD;

  // read commandline input
  // expect ./foo w n t crbc_t tol
  // where
  // w is the domain width
  // h is the grid spacing
  // t is the simulation time
  // crbc_t is the crbc time parameter
  // tol is the tolerance to determine the number of recurions
  if (argc < 6) {
    return 1;
  } else {

    w = std::atof(argv[1]);
    n = std::atoi(argv[2]);
    t = std::atof(argv[3]);
    crbc_t = std::atof(argv[4]);
    p = std::atoi(argv[5]);

  }

  // intialize MPI
  if ( MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Init" <<std::endl;
    return -1;
  }

  // start timer
  t1 = MPI_Wtime();

  if ( MPI_Comm_size(comm, &nprocs) != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Comm_size" <<std::endl;
    return -1;
  }
  if ( MPI_Comm_rank(comm, &id) != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Comm_rank" <<std::endl;
    return -1;
  }

  if (id == 0) 
    std::cout << "nprocs = " << nprocs << std::endl;

  // calculate the largest integer cube root
  procs_per_dim = floor(std::pow(nprocs + 1e-8, 1.0/3.0));

  if (id == 0) 
    std::cout << "using " << procs_per_dim << " per dimension" << std::endl;

  yee_updater solver(comm, procs_per_dim, w, n, t, crbc_t, p, 0.1, 1);

  try {
    solver.run();
    solver.print_timing_data();
    solver.free_comms();
  } catch (const std::exception& e) {
    std::cout << "id = " << id << " failed in solver.run() --- a standard exception was caught, with message '"
                  << e.what() << "'" << std::endl;
    MPI_Abort(comm, -2);
  } 

  MPI_Barrier(MPI_COMM_WORLD);
  t2 = MPI_Wtime();

  if (id == 0)
    std::cout << std::endl << "Wall time = " << t2-t1 << std::endl;

  MPI_Finalize();

  return 0;

} // end main
