.. highlight:: c

***********************************************
2-D FDTD code with CRBC/DAB Boundary Conditions
***********************************************

This tutorial explains the code `yee_TM.c <https://bitbucket.org/rbcpack/rbcpack/src/default/YeeCRBC/examples/2d_yee/yee_TM.c>`_.

Introduction
============

This C code implements the finite difference time-domain solution of
Maxwell's equations (curl formulation) for a transverse-magnetic problem.
The grid is terminated using Double Absorbing Boundaries (DAB).

This is meant to serve as an example of how one might use the CRBC/DAB 
library. As such the Yee algorithm used is as simple as is practical to 
achieve this purpose. Furthermore, when it is an option, we have opted
for clarity rather than performance. We intend this to serve as a guide
to using the CRBC/DAB library.

What this program does
----------------------

For this example, we consider a set of transverse-magnetic Maxwell's equations in
a homogeneous, isotropic, dielectric material is given by

.. math::
  :nowrap:

  \begin{align}
    \frac{\partial \mathbf{H}}{\partial t} &= - \frac{1}{\mu}  \nabla \times \mathbf{E}, \\
    \frac{\partial \mathbf{E}}{\partial t}&= \frac{1}{\varepsilon} \nabla \times \mathbf{H},  
  \end{align}

where

.. math::
  :nowrap:

  \begin{align}
    E_x &= 0, \\
    E_y &= 0, \\
    H_z &= 0.
  \end{align}

.. _discretization:

To discretize these equations using Yee's scheme on a rectangular
domain :math:`[x_L, x_R] \times [y_L, y_R]` with mesh spacings
of :math:`\Delta x` and :math:`\Delta y`, in the 
:math:`x` and :math:`y` directions, respectively, we define

.. math::
  :nowrap:

  \begin{align}
    x_i &= x_L + i\Delta x, \\
    y_j &= y_L + j\Delta y
  \end{align}

We choose a time step size, :math:`\Delta t`, satisfying

.. math::
  :nowrap:

  \begin{align}
    \Delta t \leq \frac{1}{c \sqrt{\frac{1}{(\Delta x)^2} + \frac{1}{(\Delta y)^2}}},
  \end{align}

with the wave speed given by

.. math::
  :nowrap:

  \begin{align}
    c = \frac{1}{\sqrt{\varepsilon \mu}}.
  \end{align}

Letting 

.. math::
  :nowrap:

  \begin{align}
    t_n = n \Delta t,
  \end{align}

the fields are approximated on the staggered space time grids:

.. math::
  :nowrap:

  \begin{align}
    E_z^{i,j,n+\frac{1}{2}} & \sim E_z(x_i,y_j,t_{n+\frac{1}{2}}), \\
    H_x^{i,j+\frac{1}{2},n} & \sim H_x(x_i,y_{j+\frac{1}{2}},t_n), \\
    H_y^{i+\frac{1}{2},j,n} & \sim H_y(x_{i+\frac{1}{2}},y_j,t_n),
  \end{align}

where we require that the domain is terminated such that the :math:`E_z` 
component is on the boundaries. This corresponds to having an integer number of 
Yee cells that have :math:`E_z` points located on the corners.  

Finally, the fields are evolved with

.. math::
  :nowrap:

  \begin{align}
     H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n+1}  = H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n} 
     & - \frac{\Delta t}{\mu \Delta y} \left(E_z^{i,j+1,k+\frac{1}{2},n+\frac{1}{2}} 
     - E_z^{i,j,k+\frac{1}{2},n+\frac{1}{2}}    \right),  \\
     % % % %
     H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n+1}  = H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n}
     & + \frac{\Delta t}{\mu \Delta x} \left(E_z^{i+1,j,k+\frac{1}{2},n+\frac{1}{2}} 
     - E_z^{i,j,k+\frac{1}{2},n+\frac{1}{2}} \right), \\
     % % % %
     E_z^{i,j,k+\frac{1}{2},n+\frac{1}{2}}  = E_z^{i,j,k+\frac{1}{2},n-\frac{1}{2}} 
     & + \frac{\Delta t}{\varepsilon \Delta x}  \left(H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n} 
     - H_y^{i-\frac{1}{2},j,k+\frac{1}{2},n} \right)   \\
     & - \frac{\Delta t}{\varepsilon \Delta y}  \left(H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n} 
     - H_y^{i,j-\frac{1}{2},k-\frac{1}{2},n} \right). 
  \end{align}

We drive the simulation with a impulsive point source that takes the form of
a differentiated Gaussian. We place this source at the center of the computational
domain in the :math:`E_z` field and implement it as a soft source.

To make the problem more interesting, we place a rectangular and "L" shaped 
scatterer in the domain so that they form a channel. We define this scatterer to
be a perfect electric conductor. For a configuration of the layout see 
the :ref:`fig_scatterer`

.. _fig_scatterer:
.. figure:: scatterer.png
   :align: center
   :figwidth: 500 px
   :width: 450 px
   :alt: image of the scatterer in the wave guide.

   Domain configuration.
    
The commented program
=====================

Include files
-------------

We will require the following:

For output, namely ``printf()``, we need ::

  #include <stdio.h>

For the definitions of ``sqrt()``, ``fabs()``, etc. we use ::

  #include <math.h>   

We need the standard libraries for ``malloc()``, ``free()``, etc. ::

  #include <stdlib.h>   

To get the 2D interface to the CRBC/DAB library  
(see `2d_crbc_api.h <https://bitbucket.org/jlagrone/yee-crbc-testing/src/default/src/CRBC/2d_crbc_api.h?at=default/>`_) ::

  #include <2d_crbc_api.h>


Data structures
---------------

We define a struct to hold all of the data needed for the Yee scheme.
The only thing we do differently from what might be included in a typical Yee scheme
is to include an array of boundary types that the CRBC library requires. To do this we
utilize the type ``CRBC2d_Boundaries_t`` defined in the library interface. ::      

  typedef struct yeeData {

    // constants (MKS units)
    double pi;   
    double C;     // speed of light (m/s)
    double mu0;   // permeability (V*s/(A*m))
    double eps0;  // permittivity (F/m)

    // relative permeability
    double epsR; // 1.0 corresponds to a vacuum  

    // number of time steps 
    int ntsteps;    

    // grid size
    int imax;   
    int jmax;   

    // grid spacing in each direction
    double dx;   
    double dy;   

    // time step size, (we'll use 0.99*CFL)  
    double dt;  

    // boundary conditions
    CRBC2d_Boundaries_t boundaries[4];

    // source parameters 
    double tw;      // pulse width   
    double t0;      // delay   
    double amp;     // Amplitude 

    // specify how often to generate output (in time steps) 
    int save_freq;   

    // H & E Field components  
    double **Hx;   
    double **Hy;   
    double **Ez;   
  
    // permittivity, permeability 
    double epsilon;   
    double mu;   

    // a flag so we know if data has been allocated or not
    int flag;

    // update coefficients
    double CE, CH;

  } yeeData;
       

Function prototypes
-------------------

We declare the following function prototypes, which we will later define.      

A function to allocate memory, generate parameters, etc. ::

  void initialize(yeeData *d);  

A function to deallocate memory ::

  void finalize(yeeData *d);  


A function to setup CRBC parameters. This differs from a typical Yee implementation. ::

  void setup_crbc(yeeData *d, YeeCrbcUpdater **upd); 

A function to compute the E-field updates ::

  void computeE(yeeData *d, int *tstep);

A function to compute the H-field updates ::

  void computeH(yeeData *d);

A function to compute the boundary updates (E-field only). This also differs 
from a typical Yee implementation.  ::

  void computeBoundary(yeeData *d, YeeCrbcUpdater *upd);

Finally, a function to write Ez field output ::

  void writeEzField(yeeData *d, int id); 

       
Main routine
------------

First we declare the variables we need. 
We also perform a primitive command line input check to enable the Ez field to
be written out to files periodically. This feature is off by default because it
generates approximately 1.2 GB of files.
Note the that ``CrbcUpdater2d`` type is 
provided by the CRBC library and we need to declare it as a pointer because at 
this time its size is unknown. We will initialize it at a later point after 
we have set the simulation parameters. ::

  int main(int argc, char *argv[]) {   

    int tstep, i, j;
    double norm, enorm, hnorm;
    int write_output_files = 0;

    // read in input to see if we should write output files. If output files are
    // enabled this program writes out the Ez field 450 times. Each file is
    // file is approximately 2.6 MB, so this generates about 1.2 GB of data.
    // By default, the file output is turned off.
    // There is only one option, so for simplicity we'll just assume that if we
    // receive any command line option, then we should enable output instead of
    // actually parsing and identifying a specific option.
    if (argc > 1) {
      printf("The Ez field will be saved. \n");
      write_output_files = 1;
    }

    int tstep, i, j;
    double norm, enorm, hnorm;
  
    // declare a yee data structure
    yeeData data;

    // declare a boundary updater object.
    // This needs to be a pointer because the size is unknown at this time and will
    // be initialized later.
    CrbcUpdater2d *boundary_updater;

Next, we set the basic simulation parameters ::
  
  // relative permeability
  data.epsR = 1.0; // 1.0 corresponds to a vacuum  

  // number of time steps 
  data.ntsteps = 4500;   

  // grid size
  // NOTE: Field Value Locations
  //
  // Ez - at (i, j, k + 1/2)
  //  
  // Hx - at (i, j + 1/2, k + 1/2)
  // Hy - at (i + 1/2, j, k + 1/2)
  data.imax = 900;   
  data.jmax = 300;   

  // grid spacing in each direction
  data.dx = 5e-4;   
  data.dy = 5e-4;   

The boundary conditions are set using the ``CRBC2d_Boundaries_t`` 
enumeration type provided by the CRBC library.  At this time, 
the supported boundary types are 

* `CRBC2d_DIR`  --- Homogeneous Dirichlet
* `CRBC2d_NEUM` --- Homogeneous Neumann
* `CRBC2d_CRBC` --- Complete Radiation Boundary Condition (implemented as a DAB)
  

The library also provides enumerations that list the valid sides with the ``CRBC_Side_t``
type, they are

*  `CRBC_XLeft`
*  `CRBC_XRight`
*  `CRBC_YLeft`
*  `CRBC_YRight`
  

Here, we'll set the boundary conditions so that we a waveguide with parallel
PEC plates with normals in the y-direction. ::


  data.boundaries[CRBC2d_XLeft]  = CRBC2d_CRBC;
  data.boundaries[CRBC2d_XRight] = CRBC2d_CRBC;
  data.boundaries[CRBC2d_YLeft]  = CRBC2d_DIR;
  data.boundaries[CRBC2d_YRight] = CRBC2d_DIR;

Finally, we set the source parameters. For this example we will use an impulsive 
source that takes the form of a differentiated Gaussian. We will center the source in the 
Ex field and implement as a soft source. Additionally, we set the the frequency 
at which we generate output and call the ``initialize`` function to calculate the 
rest of the parameters such as the time step size. ::

  // source parameters 
  data.tw = 5E-11;           // pulse width   
  data.tO = 4.0 * data.tw;   // delay (this needs to be >= 4*tw to be smooth)
  data.amp = 1000;           // Amplitude 

  // specify how often to generate output (in time steps) 
  data.save_freq = 10;   

  // allocate memory and compute the remaining parameters
  initialize(&data);   

Now, we need to setup the boundary updater. The updater has 3 basic 
parameters: 

* delta --- The minimum separation of each boundary face from sources, scatterers, and other inhomogeneities.
* T     --- The total time being simulated.
* P     --- The number of recursions to use in approximating the boundary condition. (This can be chosen by providing a tolerance)

It additionally needs to know the wave speed, *c*, and the boundary conditions
on the rest of the domain. Finally, it needs to know some properties
of the discretization such as the time step size, mesh sizes, and the 
indexing for each face with a CRBC boundary.
Each of these is discussed in more detail in the ``setup_crbc`` function. 

Note that we need to pass the reference to the boundary_updater pointer because 
it is currently uninitialized. ::

  setup_crbc(&data, &boundary_updater); 

Finally, we can run the simulation by stepping the **H** fields followed by the 
**E** fields and then applying the boundary updates to the **E** fields from the
CRBC library as needed and repeating. We also generate some output at the 
requested intervals. ::

  // Start time stepping
  for (tstep=1; tstep <= data.ntsteps; ++tstep) {

    // compute the updates to the H-field
    computeH(&data);

    // compute the updates to the E-field
    computeE(&data, &tstep);

    // compute the updates to the boundaries of the E-field.
    computeBoundary(&data, boundary_updater);

    // print output if needed
    if (tstep % data.save_freq == 0) {

      // print out a norm to the screen.
      enorm = 0.0;
      hnorm = 0.0;

      // loop over Ez
      for (i=0; i<data.imax; ++i)
        for (j=0; j<data.jmax; ++j)
          enorm += data.Ez[i][j][k] * data.Ez[i][j][k];

      // loop over Hx
      for (i=0; i<data.imax; ++i)
        for (j=0; j<data.jmax-1; ++j)
          hnorm += data.Hx[i][j][k] * data.Hx[i][j][k];

      // loop over Hy
      for (i=0; i<data.imax-1; ++i)
        for (j=0; j<data.jmax; ++j)
          hnorm += data.Hy[i][j][k] * data.Hy[i][j][k];

      // compute norm of combined fields
      norm = sqrt(data.epsilon*enorm + data.mu*hnorm);
      
      // write to screen
      printf("tstep = %i \t |E| = %6.5e \t |H| = %6.5e \t l2 norm = %6.5e \n", \
             tstep, sqrt(enorm), sqrt(hnorm), norm);

      // write all of the Ez field data to a file so we can visualize it if
      // requested.
      if (write_output_files != 0) {
        writeEzField(&data, tstep / data.save_freq);
      }

    } // end output

After doing the updates, we demonstrate the ability to save and restart the
boundary updater. This capability requires the HDF5 library. When building the
CRBC library, HDF5 is optional, so this is only enabled when this option is turned
on.

This will write out the current state of the boundary library to 
several h5 files. In this case, it will save a file
called hdf5_save_test.h5 For demonstration
purposes, we delete the current boundary updater object to free all of the memory
associated to it and then create a new boundary updater object in its place using
the saved data half way through the simulation. :: 

    // test the restart capability (requires HDF5)
    #if USE_HDF5_RESTARTS
    if (tstep == data.ntsteps /2) {

      // This will write out the current state of the boundary library to 
      // several h5 files. In this case, it will save a file
      // called hdf5_save_test.h5 
      if (CRBC2d_save_state(boundary_updater, "hdf5_save_test") != 0)
        fprintf(stderr, "save_stat() failed \n");

      // delete the updater object to ensure that everything is working
      // correctly ...
      CRBC2d_delete_updater(boundary_updater);

      // Now restart. This is an alternative way to generate a boundary updater
      // object if a previously saved state exists.
      if (CRBC2d_restart(&boundary_updater, "hdf5_save_test") != 0)
        fprintf(stderr, "restart() failed \n");

    }  
    #endif


  } // finished time stepping

After the time stepping has been completed, we need to free the dynamically 
allocated memory before terminating the program. ::

    // free the memory associated with the solver
    finalize(&data);

    // free the memory associated with the boundary updater
    CRBC_delete_updater(boundary_updater);

    return 0;

  } // end main

Initialize function
-------------------

We use this function primarily to allocate the memory to store each of the field 
components. Note that we are using multidimensional arrays simply for clarity, 
but this should probably not be done in applications were performance is important.

Based on the :ref:`discretization <discretization>` presented, we have

.. math::
  :nowrap:

  \begin{align}  
    (imax) \cdot (jmax) \cdot (kmax-1) & & E_z \text{ field values} \\
    (imax) \cdot (jmax-1) \cdot (kmax-1) & & H_x \text{ field values} \\
    (imax-1) \cdot (jmax) \cdot (kmax-1) & & H_y \text{ field values} \\
  \end{align}

Additionally, we compute the constants that are used in the field update equations. ::

  void initialize(yeeData *d) {

    int i, j;
    int imax = d->imax;
    int jmax = d->jmax;

    d->flag = 0;

    d->pi = 3.14159265358979;   
    d->C = 2.99792458e8;   

    // time step size, (we'll use 0.99*CFL)  
    d->dt = 0.99/(d->C*sqrt(1.0/(d->dx*d->dx) + 1.0/(d->dy*d->dy))); 

    // calculate free space eps and mu
    d->mu0 = 4.0 * d->pi * 1.0E-7;   
    d->eps0 = 1.0 / (d->C * d->C * d->mu0);

    // calculate material epsilon and mu
    d->epsilon = d->epsR * d->eps0;
    d->mu = d->mu0;

    // calculate update coefficients
    d->CE = d->dt / d->epsilon;
    d->CH = d->dt / d->mu0;

    d->Ez = (double ***) malloc((imax) * sizeof(double));
    for (i=0; i < d->imax; i++) {
      d->Ez[i] = (double **) malloc((jmax) * sizeof(double));
      for (j=0; j < d->jmax; j++) {
        d->Ez[i][j] = 0.0;
      }
    }

    d->Hx = (double ***) malloc((imax) * sizeof(double));
    for (i=0; i < d->imax; i++) {
      d->Hx[i] = (double **) malloc((jmax-1) * sizeof(double));
      for (j=0; j < d->jmax-1; j++) {
        d->Hx[i][j] = 0.0;
      }
    }

    d->Hy = (double ***) malloc((imax-1) * sizeof(double));
    for (i=0; i < d->imax-1; i++) {
      d->Hy[i] = (double **) malloc((jmax) * sizeof(double));
      for (j=0; j < d->jmax; j++) {
        d->Hy[i][j] = 0.0;
      }
    }

    d->flag = 1;
  } 


Finalize function
-----------------

This function is used to free the memory allocated to store the field values. ::

  void finalize(yeeData *d) {

    int i, j;

    if (d->flag != 0) {

      for (i=0; i < d->imax; i++) {
        free(d->Ez[i]);
      }
      free(d->Ez);

      for (i=0; i < d->imax; i++) {
        free(d->Hx[i]);
      }
      free(d->Hx);

      for (i=0; i < d->imax-1; i++) {
        free(d->Hy[i]);
      }
      free(d->Hy);
    }

    d->flag = 0;
  } 

setup_crbc function
-------------------

This function is used to initialize the interface to the CRBC library. This 
happens in two parts, first we initialize the boundary updater object with the 
basic simulation parameters such as the discretization and wave speed. After we 
have an initialized boundary updater, we need to set up the parameters for each 
of the boundary faces that we want to use the CRBC library to update.

The first step is straightforward, and we simply provide the total simulation 
time, the grid spacings, the time step size, the wave speed and the boundary 
conditions in the format the CRBC library expects. Note that there are two 
alternative methods that we could use to create a new boundary updater. The 
first is to call ``CRBC_new_updater_p(...)``, which does the same thing but changes the
number of recursions used on the boundary from the default of 5 to the requested
number. The second is to call ``CRBC_new_updater_tol(...)``, in which the minimum 
number of recursions needed to meet the requested tolerance is chosen by the
library. ::

  void setup_crbc(yeeData *d, CrbcUpdater2d **upd) {

    int i, l;
    double T;
    double h[2];
    double delta;
    int low_index[2];
    int high_index[2];
    int n;
    CRBC2d_Side_t side;
  

    // First we will create a new boundary updater.
    h[0] = d->dx;
    h[1] = d->dy;
    T = d->dt * d->ntsteps; // total time is time step size * number of time steps
    *upd = CRBC2d_new_updater_tol(T, h, d->dt, d->C, d->boundaries);

    // alternatively one can call 
    // int P = 7;
    // upd = CRBC2d_new_updater_p(T, h, d->dt, d->C, d->boundaries, P);
    // This does exactly the same thing but changes the number of recursions from
    // the default of 5 to 7.
    //
    // or 
    // int Pmax = 15;
    // double tol = 1e-3;
    // upd = CRBC2d_new_updater_tol(T, h, d->dt, d->C, d->boundaries, Pmax, tol);
    // Will generate an updater that attempts to satsify the provided tolerance
    // and with fewer than Pmax recursions.

Before we precede, we'll print out the properties from the boundary updater to 
make sure they are correct. ::

  printf("The boundary updater was initialized with the following: \n");
  printf("  wave speed, c = %e \n", CRBC2d_get_c(*upd));
  printf("  total time, T = %e \n", CRBC2d_get_T(*upd));
  printf("  time step, dt = %e \n", CRBC2d_get_dt(*upd));
  printf("The boundary updater calculated the following: \n");
  printf("  %i edges \n", CRBC2d_get_num_sides(*upd));
  printf("  %i corners \n", CRBC2d_get_num_corners(*upd));

Now set up the faces. Start by looping over all of the possible faces. We
follow the order given in ``CRBC2d_Side_t``, so 

* l = 0 --> CRBC_XLeft
* l = 1 --> CRBC_XRight
* l = 2 --> CRBC_YLeft
* l = 3 --> CRBC_YRight

First we need to calculate the minimum distance between the boundary and
source, `delta`. Since we placed the source in the center of the domain, this is 
simply the distance from the boundary to the center of the domain in the
direction normal to the boundary. Additionally, we need to factor in the location
of the scatterer if it is present.
In general, if it is not possible to calculate
delta directly, using a lower bound for the separation is the safest thing to do,
but it may result in more work being done that is actually needed to achieve the
desired accuracy. ::


  for (l=0; l<4; ++l) {

    switch (l) {
   
      case 0: // left side in x
        delta = (d->imax / 2) * (d->dx);
        break;

      case 1: // right side in x
        // if the domain is big enough there is a scatterer ...
        if ((d->imax > 700) & (d->jmax > 250)) {
          delta = (d->imax - (d->imax + 250)) * (d->dx);
        } else {
          delta = (d->imax / 2) * (d->dx);
        }

      case 2: // left side in y

        // if the domain is big enough there is a scatterer
        if ((d->imax > 700) & (d->jmax > 250)) {
          delta = (d->jmax - 55) * (d->dy);
        } else {
          delta = (d->jmax / 2) * (d->dy);
        }
        break;

      case 3: // right side in y

        // if the domain is big enough there is a scatterer
        if ((d->imax > 700) & (d->jmax > 250)) {
          delta = (d->jmax - (d->jmax / 2 + 55)) * (d->dy);
        } else {
          delta = (d->jmax / 2) * (d->dy);
        }
        break;

    }

    // convert the index l to a CRBC2d_Side_t type
    switch (l) {
      case 0:
        side = CRBC2d_XLeft;
        break;
      case 1:
        side = CRBC2d_XRight;
        break;
      case 2:
        side = CRBC2d_YLeft;
        break;
      case 3:
        side = CRBC2d_YRight;
        break;
    }

The CRBC updater library attempts to communicate with the solvers "native"
indexing scheme, so we need to tell it the upper and lower indexing extents.
These extents need to include all of the requested variables
on the boundary face as well the adjacent parallel plane of points in the 
interior of the domain.

For the left boundary in x, this means we need to include the data 
points of the requested component at i=0 and i=1 and all possible 
indices for j. Therefore, the lower indexing extents are 0 for
all indices. The boundary updater considers these extents to be inclusive.

For the right boundary in x, this means we need to include the data
points of the requested component at i=imax-1 and i=imax-2 along with all of 
the possible indices for j.

For the left boundary in y, this means we need to include the data 
points of the requested component at j=0 and j=1 along with all of 
the possible indices for i.

Finally, for the right boundary in y, this means we need to include the data
points of the requested component at j=jmax-1 and j=jmax-2 along with all of 
the possible indices for i. ::

    if (d->boundaries[side] == CRBC2d_CRBC) {

      // Set the basic extents and then adjust them based on the side and field 
      // component. These should be inclusive extents;
      low_index[0] = 0;
      low_index[1] = 0;

      high_index[0] = d->imax-1; // minus 1 from 0-based indexing
      high_index[1] = d->jmax-1;

      switch (side) {

        case CRBC2d_XLeft:
          high_index[0] = 1;  
    
          break;

        case CRBC2d_XRight:
          low_index[0] = d->imax-2;

          break;

        case CRBC2d_YLeft:
          high_index[1] = 1;
    
          break;  

        case CRBC2d_YRight:
          low_index[1] = d->jmax-2;

          break;

      } // end switch side

Now, we can initialize the face: ::

      // now initialize the face.
      if (CRBC2d_init_face(*upd, side, low_index, high_index, delta) != 0)
      {
        fprintf(stderr, "Error: init_face(...) failed \n");
        exit(-1);
      }
    } // end if
  } // end for over possible boundary faces (l)

Finally, we print out some information about the recursions. Note that a reflection
coefficient of -1 indicates that we are not performing any updates on that face. ::

    // now we'll print out some information about the recursions
    printf("The faces were initialized with: \n");
    printf("  Left side in  x is using %i recursions "
           "with a reflection coefficient of %e \n"
           , CRBC2d_get_num_recursions(*upd, CRBC2d_XLeft), 
            CRBC2d_get_reflection_coef(*upd, CRBC2d_XLeft));
    printf("  Right side in x is using %i recursions "
           "with a reflection coefficient of %e \n"
           , CRBC2d_get_num_recursions(*upd, CRBC2d_XRight), 
           CRBC2d_get_reflection_coef(*upd, CRBC2d_XRight));
    printf("  Left side in  y is using %i recursions "
           "with a reflection coefficient of %e \n"
           , CRBC2d_get_num_recursions(*upd, CRBC2d_YLeft), 
           CRBC2d_get_reflection_coef(*upd, CRBC2d_YLeft));
    printf("  Right side in y is using %i recursions "
           "with a reflection coefficient of %e \n" 
           , CRBC2d_get_num_recursions(*upd, CRBC2d_YRight), 
           CRBC2d_get_reflection_coef(*upd, CRBC2d_YRight));
    printf("The maximum reflection coefficient is %e \n", 
           CRBC2d_get_max_reflection_coef(*upd));

  } // end setup_crbc     

computeE function
-----------------

This function computes the updates to the **E** field components using the Yee scheme. ::

  void computeE(yeeData *d, int *tstep) {

    int i, j, k;

    // compute updates to Ez
    for (i=1; i < d->imax-1; ++i) {
      for (j=1; j < d->jmax-1; ++j) {
         d->Ez[i][j] += d->CE * ((d->Hy[i][j] - d->Hy[i-1][j]) / d->dx \
                        - (d->Hx[i][j] - d->Hx[i][j-1]) / d->dy);
      }
    }

    // add the source term to the center point of the Ez field
    // This is a differentiated Gaussian applied as a "soft" source
    d->Ez[d->imax/2][d->jmax/2] += 2.0 * d->amp * d->CE \
                                * ((*tstep*d->dt - d->t0) / d->tw) \
                                * exp(-pow(((*tstep*d->dt - d->t0) / d->tw), 2.0));


    // add scatter if the domain is large enough
    // rectangular with L shapped channel
    if ((d->imax > 700) & (d->jmax > 250)) {

      for (i = d->imax/2 + 50; i <= d->imax/2 + 175; ++i) 
        for (j = d->jmax/2 + 45; j <= d->jmax/2 + 55; ++j)
          d->Ez[i][j] = 0.0;

      for (i = d->imax/2 + 50; i <= d->imax/2 + 250; ++i) 
        for (j = d->jmax/2 - 55; j <= d->jmax/2 - 45; ++j)
          d->Ez[i][j] = 0.0;

      for (i = d->imax/2 + 230; i <= d->imax/2 + 250; ++i) 
        for (j = d->jmax/2 - 45; j <= d->jmax/2 + 55; ++j)
          d->Ez[i][j] = 0.0;

    } 
  }

computeH function
-----------------

This function computes the updates to the **H** field components using the Yee scheme. ::

  void computeH(yeeData *d) {
 
    int i, j, k;

    // compute updates to Hx
    for (i=0; i < d->imax; ++i) {
      for (j=0; j < d->jmax-1; ++j) {
        d->Hx[i][j] += d->CH * (-(d->Ez[i][j+1] - d->Ez[i][j]) / d->dy);
      }
    }

    // compute updates to Hy
    for (i=0; i < d->imax-1; ++i) {
      for (j=0; j < d->jmax; ++j) {
        d->Hy[i][j] += d->CH * ((d->Ez[i+1][j] - d->Ez[i][j]) / d->dx);
      }
    }
  }

computeBoundary Function
------------------------

This function computes the boundary updates using the CRBC library. To do this, 
we first need to copy the values that have been updated by the Yee algorithm into
the CRBC updater. We start by looping over all of the possible sides. For each of 
these, the CRBC updater can return the index extents, which we use to copy 
the data into the CRBC library. ::

  void computeBoundary(yeeData *d, CrbcUpdater2d *upd) {

    int i, j, k, l, m;
    int low_index[2];
    int high_index[2];
    int index[2];
    int n;
    CRBC2d_Side_t side;

    // first we need to copy the values that have been updated by the Yee algorithm
    // into the crbc updater. First we will loop over all of the possible sides.
    for (m=0; m<4; ++m) {

      // convert index to CRBC2d_Side type
      switch (m) {
        case 0:
          side = CRBC2d_XLeft;
          break;
        case 1:
          side = CRBC2d_XRight;
          break;
        case 2:
          side = CRBC2d_YLeft;
          break;
        case 3:
          side = CRBC2d_YRight;
          break;
      }

      if (d->boundaries[side] == CRBC2d_CRBC) {

        // get the indices the updater expects as input.
        // These indices are inclusive.
        CRBC2d_get_input_extents(upd, side, low_index, high_index);

        // copy data into the updater
        for (i=low_index[0]; i<=high_index[0]; ++i) {
          index[0] = i;
          for (j=low_index[1]; j<=high_index[1]; ++j) {
            index[1] = j;
            CRBC2d_load_face_data(upd, side, index, &(d->Ez[i][j]));
          }
        }
      }
    } // end loop over sides

Now we can have the CRBC library compute the boundary updates. ::

  // now we can compute the updates to the boundaries
  if (CRBC2d_compute_updates(upd) != 0) {
    fprintf(stderr, "Error: compute_updates(...) failed \n");
    exit(-1);
  }

Finally, we need to copy the new values from the boundary updater into the arrays 
used by the Yee algorithm. To do this, we loop over all the possible sides and
get the output index extents from the boundary updater. We use these extent 
values to loop over the data arrays and copy the values. ::

    // Now copy the new boundary values into the arrays used by the Yee algorithm
    // loop over possible sides
    for (m=0; m<4; ++m) {

      // convert index to CRBC2d_Side type
      switch (m) {
        case 0:
          side = CRBC2d_XLeft;
          break;
        case 1:
          side = CRBC2d_XRight;
          break;
        case 2:
          side = CRBC2d_YLeft;
          break;
        case 3:
          side = CRBC2d_YRight;
          break;
      }
   
      if (d->boundaries[side] == CRBC2d_CRBC) {

        // get the indices the updater can output.
        // These indices are inclusive.
        CRBC2d_get_output_extents(upd, side, low_index, high_index);

        // copy data into the updater
        for (i=low_index[0]; i<=high_index[0]; ++i) {
          index[0] = i;
          for (j=low_index[1]; j<=high_index[1]; ++j) {
            index[1] = j;
            d->Ez[i][j] = CRBC2d_get_new_face_vals(upd, side, index);
          }
        }
      }
    } // end loop over sides
  

  } // end computeBoundary

writeEzField function
---------------------

We use this function to output the :math:`E_z` field component to an ASCII VTK 
file format that can be opened with visualization software such as ParaView. ::

  void writeEzField(yeeData *d, int id) {

    int i, j, k, n, cells;

    char step[10];   
    char fileBaseName[] = "Ez_Field_";   
    sprintf(step, "%d", id);   
    strcat(fileBaseName, step);   
    strcat(fileBaseName, ".vtk");   

    // open the file and write the VTK header information
    FILE *f = fopen(fileBaseName, "w");
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "vtk output\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET RECTILINEAR_GRID\n");

    // set the dimensions
    fprintf(f, "DIMENSIONS %i %i %i\n", d->imax, d->jmax, 1);

    // save the coordinates
    fprintf(f, "X_COORDINATES %i double \n", d->imax);
    for (i=0; i < d->imax; ++i)
      fprintf(f, "%f\n", i*d->dx);

    fprintf(f, "Y_COORDINATES %i double \n", d->jmax);
    for (j=0; j < d->jmax; ++j)
      fprintf(f, "%f\n", j*d->dy);

    fprintf(f, "Z_COORDINATES %i double \n", 1);
    fprintf(f, "0 \n");

    // set up a cell and field
    n = (d->imax) * d->jmax;
    cells = (d->imax-1) * (d->jmax-1);
    fprintf(f, "CELL_DATA %i\n", cells);
    fprintf(f, "POINT_DATA %i\n", n);
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "Ez 1 %i double\n", n);

    // now write out the data
    for (j=0; j < d->jmax; ++j)
      for (i=0; i < d->imax; ++i)
        fprintf(f, "%f\n", d->Ez[i][j]);
  	

    // close the file
    fclose(f);
  }

Output
======
A visualization of the 
`Ez field <https://drive.google.com/file/d/0B2tqyaSF82iOTnRRRE5iY3g4M0k/view?usp=sharing/>`_
shows on the top the error from a simulation using a tolerance of :math:`1e-3`,
in the middle the error using a tolerance of :math:`1e-5`, and on the :math:`E_z`
field on the bottom. The errors were generated by comparing to a simulation run
on a larger domain with the same parameters (a longer wave guide). Note that the
color scale for the :math:`E_z` field is chosen to be representative of the state
roughly midway through the simulation. Also note that the error scales are linear.
The files used to generate this movie can be generated by providing a command line
option at runtime, for instance ::

  ./yee_TM.x -output

The following is a truncated sample of the screen output. ::

  The boundary updater was initialized with the following: 
    wave speed, c = 2.997925e+08 
    total time, T = 5.253903e-09 
    time step, dt = 1.167534e-12 
  The boundary updater calculated the following: 
    2 edges 
    0 corners 
  The faces were initialized with: 
    Left side in  x is using 3 recursions with a reflection coefficient of 6.640539e-04 
    Right side in x is using 3 recursions with a reflection coefficient of 9.564357e-04 
    Left side in  y is using 5 recursions with a reflection coefficient of -1.000000e+00 
    Right side in y is using 5 recursions with a reflection coefficient of -1.000000e+00 
  The maximum reflection coefficient is 9.564357e-04 
  tstep = 10 	 |E| = 3.43084e-04 	 |H| = 1.71082e-06 	 l2 norm = 2.17261e-09 
  tstep = 20 	 |E| = 1.58550e-03 	 |H| = 8.76478e-06 	 l2 norm = 1.08993e-08 
  tstep = 30 	 |E| = 7.20991e-03 	 |H| = 4.08908e-05 	 l2 norm = 5.06106e-08 
  tstep = 40 	 |E| = 2.94279e-02 	 |H| = 1.70488e-04 	 l2 norm = 2.10223e-07 
  tstep = 50 	 |E| = 1.07180e-01 	 |H| = 6.34234e-04 	 l2 norm = 7.79229e-07 
  tstep = 60 	 |E| = 3.47840e-01 	 |H| = 2.10338e-03 	 l2 norm = 2.57505e-06 
  tstep = 70 	 |E| = 1.00479e+00 	 |H| = 6.21205e-03 	 l2 norm = 7.57841e-06 
  tstep = 80 	 |E| = 2.57999e+00 	 |H| = 1.63160e-02 	 l2 norm = 1.98361e-05 
  tstep = 90 	 |E| = 5.87854e+00 	 |H| = 3.80430e-02 	 l2 norm = 4.60941e-05 
  tstep = 100 	 |E| = 1.18606e+01 	 |H| = 7.85520e-02 	 l2 norm = 9.48659e-05 
  tstep = 110 	 |E| = 2.11351e+01 	 |H| = 1.43146e-01 	 l2 norm = 1.72350e-04 
  tstep = 120 	 |E| = 3.31680e+01 	 |H| = 2.29084e-01 	 l2 norm = 2.75115e-04 
  tstep = 130 	 |E| = 4.57490e+01 	 |H| = 3.19572e-01 	 l2 norm = 3.83233e-04 
  tstep = 140 	 |E| = 5.56463e+01 	 |H| = 3.84126e-01 	 l2 norm = 4.61343e-04 
  tstep = 150 	 |E| = 6.10928e+01 	 |H| = 3.90875e-01 	 l2 norm = 4.74384e-04 
  tstep = 160 	 |E| = 6.50698e+01 	 |H| = 3.31051e-01 	 l2 norm = 4.18581e-04 
  tstep = 170 	 |E| = 7.41876e+01 	 |H| = 2.54907e-01 	 l2 norm = 3.61089e-04 
  tstep = 180 	 |E| = 8.97395e+01 	 |H| = 2.76852e-01 	 l2 norm = 4.09416e-04 
  tstep = 190 	 |E| = 1.05730e+02 	 |H| = 3.77805e-01 	 l2 norm = 5.27587e-04 
  tstep = 200 	 |E| = 1.16836e+02 	 |H| = 4.54163e-01 	 l2 norm = 6.16494e-04 
  tstep = 210 	 |E| = 1.22165e+02 	 |H| = 4.73905e-01 	 l2 norm = 6.43711e-04 
  tstep = 220 	 |E| = 1.23991e+02 	 |H| = 4.51980e-01 	 l2 norm = 6.26765e-04 
  tstep = 230 	 |E| = 1.25525e+02 	 |H| = 4.16101e-01 	 l2 norm = 5.97566e-04 
  tstep = 240 	 |E| = 1.29196e+02 	 |H| = 3.83127e-01 	 l2 norm = 5.76410e-04 
  tstep = 250 	 |E| = 1.34828e+02 	 |H| = 3.56470e-01 	 l2 norm = 5.66249e-04 
  tstep = 260 	 |E| = 1.39206e+02 	 |H| = 3.39236e-01 	 l2 norm = 5.62311e-04 
  tstep = 270 	 |E| = 1.39147e+02 	 |H| = 3.37368e-01 	 l2 norm = 5.60769e-04 
  tstep = 280 	 |E| = 1.35141e+02 	 |H| = 3.48050e-01 	 l2 norm = 5.60297e-04 
  tstep = 290 	 |E| = 1.31040e+02 	 |H| = 3.59131e-01 	 l2 norm = 5.60460e-04 
  tstep = 300 	 |E| = 1.29613e+02 	 |H| = 3.63101e-01 	 l2 norm = 5.60736e-04 
  tstep = 310 	 |E| = 1.29709e+02 	 |H| = 3.62890e-01 	 l2 norm = 5.60760e-04 
  tstep = 320 	 |E| = 1.28959e+02 	 |H| = 3.64596e-01 	 l2 norm = 5.60618e-04 
  tstep = 330 	 |E| = 1.27304e+02 	 |H| = 3.68672e-01 	 l2 norm = 5.60620e-04 
  tstep = 340 	 |E| = 1.27330e+02 	 |H| = 3.69009e-01 	 l2 norm = 5.60951e-04 
  tstep = 350 	 |E| = 1.31880e+02 	 |H| = 3.58234e-01 	 l2 norm = 5.61482e-04 
  tstep = 360 	 |E| = 1.40134e+02 	 |H| = 3.35711e-01 	 l2 norm = 5.61694e-04 
  tstep = 370 	 |E| = 1.46015e+02 	 |H| = 3.16635e-01 	 l2 norm = 5.61037e-04 
  tstep = 380 	 |E| = 1.42905e+02 	 |H| = 3.24851e-01 	 l2 norm = 5.59847e-04 
  tstep = 390 	 |E| = 1.31251e+02 	 |H| = 3.57269e-01 	 l2 norm = 5.59401e-04 
  tstep = 400 	 |E| = 1.21396e+02 	 |H| = 3.82067e-01 	 l2 norm = 5.60288e-04 
  tstep = 410 	 |E| = 1.23271e+02 	 |H| = 3.79194e-01 	 l2 norm = 5.61457e-04 
  tstep = 420 	 |E| = 1.32254e+02 	 |H| = 3.57503e-01 	 l2 norm = 5.61674e-04 
  tstep = 430 	 |E| = 1.37921e+02 	 |H| = 3.41226e-01 	 l2 norm = 5.61020e-04 
  tstep = 440 	 |E| = 1.36918e+02 	 |H| = 3.43243e-01 	 l2 norm = 5.60391e-04 
  tstep = 450 	 |E| = 1.32606e+02 	 |H| = 3.54882e-01 	 l2 norm = 5.60320e-04 
  tstep = 460 	 |E| = 1.29885e+02 	 |H| = 3.62315e-01 	 l2 norm = 5.60653e-04 
  tstep = 470 	 |E| = 1.30267e+02 	 |H| = 3.61635e-01 	 l2 norm = 5.60887e-04 
  tstep = 480 	 |E| = 1.31233e+02 	 |H| = 3.59073e-01 	 l2 norm = 5.60812e-04 
  tstep = 490 	 |E| = 1.31569e+02 	 |H| = 3.58252e-01 	 l2 norm = 5.60849e-04 
  tstep = 500 	 |E| = 1.34372e+02 	 |H| = 3.51395e-01 	 l2 norm = 5.61282e-04 
  tstep = 510 	 |E| = 1.40548e+02 	 |H| = 3.34162e-01 	 l2 norm = 5.61448e-04 
  tstep = 520 	 |E| = 1.43762e+02 	 |H| = 3.23352e-01 	 l2 norm = 5.60700e-04 
  tstep = 530 	 |E| = 1.38384e+02 	 |H| = 3.38249e-01 	 l2 norm = 5.59763e-04 
  tstep = 540 	 |E| = 1.28248e+02 	 |H| = 3.65401e-01 	 l2 norm = 5.59833e-04 
  tstep = 550 	 |E| = 1.23479e+02 	 |H| = 3.77871e-01 	 l2 norm = 5.60742e-04 
  tstep = 560 	 |E| = 1.26905e+02 	 |H| = 3.70457e-01 	 l2 norm = 5.61296e-04 
  tstep = 570 	 |E| = 1.31474e+02 	 |H| = 3.58793e-01 	 l2 norm = 5.61087e-04 
  tstep = 580 	 |E| = 1.32941e+02 	 |H| = 3.54627e-01 	 l2 norm = 5.60818e-04 
  tstep = 590 	 |E| = 1.33929e+02 	 |H| = 3.52179e-01 	 l2 norm = 5.60962e-04 
  tstep = 600 	 |E| = 1.36585e+02 	 |H| = 3.45048e-01 	 l2 norm = 5.61063e-04 
  tstep = 610 	 |E| = 1.37870e+02 	 |H| = 3.40945e-01 	 l2 norm = 5.60694e-04 
  tstep = 620 	 |E| = 1.34936e+02 	 |H| = 3.48586e-01 	 l2 norm = 5.60276e-04 
  tstep = 630 	 |E| = 1.30379e+02 	 |H| = 3.60738e-01 	 l2 norm = 5.60391e-04 
  tstep = 640 	 |E| = 1.29171e+02 	 |H| = 3.64382e-01 	 l2 norm = 5.60877e-04 
  tstep = 650 	 |E| = 1.32037e+02 	 |H| = 3.57421e-01 	 l2 norm = 5.61158e-04 
  tstep = 660 	 |E| = 1.35418e+02 	 |H| = 3.48239e-01 	 l2 norm = 5.61036e-04 
  tstep = 670 	 |E| = 1.36695e+02 	 |H| = 3.44384e-01 	 l2 norm = 5.60787e-04 
  tstep = 680 	 |E| = 1.36283e+02 	 |H| = 3.45405e-01 	 l2 norm = 5.60688e-04 
  tstep = 690 	 |E| = 1.35641e+02 	 |H| = 3.47197e-01 	 l2 norm = 5.60701e-04 
  tstep = 700 	 |E| = 1.34757e+02 	 |H| = 3.49504e-01 	 l2 norm = 5.60614e-04 
  tstep = 710 	 |E| = 1.32396e+02 	 |H| = 3.55566e-01 	 l2 norm = 5.60424e-04 
  tstep = 720 	 |E| = 1.28737e+02 	 |H| = 3.64877e-01 	 l2 norm = 5.60397e-04 
  tstep = 730 	 |E| = 1.26412e+02 	 |H| = 3.70887e-01 	 l2 norm = 5.60668e-04 
  tstep = 740 	 |E| = 1.27580e+02 	 |H| = 3.68398e-01 	 l2 norm = 5.60949e-04 
  tstep = 750 	 |E| = 1.30980e+02 	 |H| = 3.59721e-01 	 l2 norm = 5.60811e-04 
  tstep = 760 	 |E| = 1.34000e+02 	 |H| = 3.50872e-01 	 l2 norm = 5.60083e-04 
  tstep = 770 	 |E| = 1.36110e+02 	 |H| = 3.43459e-01 	 l2 norm = 5.58810e-04 
  tstep = 780 	 |E| = 1.38204e+02 	 |H| = 3.35267e-01 	 l2 norm = 5.57107e-04 
  tstep = 790 	 |E| = 1.39483e+02 	 |H| = 3.29196e-01 	 l2 norm = 5.55379e-04 
  tstep = 800 	 |E| = 1.37480e+02 	 |H| = 3.33680e-01 	 l2 norm = 5.54316e-04 

