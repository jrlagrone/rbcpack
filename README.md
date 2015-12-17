# RBCPack Repository #

The purpose of this library is to provide accurate and efficient radiation and absorbing boundary conditions to truncate the computational domain for the simulation of a variety of wave propagation problems in the time domain. The plan is for RBCPack to contain stand-alone components that implement complete radiation boundary conditions and/or double absorbing boundary layers for a variety of popular time-domain wave propagation solvers.

Notable features of RBCPack are

* *a priori* error estimates
* optimal selection of parameters given a desired tolerance

More details are available at [rbcpack.org](http://www.rbcpack.org)

### How do I get set up? ###

At this time, only the Yee/FDTD component and some basic utilities are available. After downloading, it is possible to build using CMake with most compilers. In some cases, it is necessary to tell the compiler to use c++11. Optional features require HDF5 and/or OpenMP. 