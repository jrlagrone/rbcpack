.. highlight:: c

***********************************************
3-D FDTD code with CRBC/DAB Boundary Conditions
***********************************************

This tutorial explains the code `yee.c <https://github.com/jrlagrone/rbcpack/blob/main/YeeCRBC/examples/3d_yee/yee.c>`_.

Introduction
============

This C code implements the finite difference time-domain solution of
Maxwell's equations (curl formulation) over a 3-D Cartesian lattice.
The grid is terminated using Double Absorbing Boundaries (DAB).

This is meant to serve as an example of how one might use the CRBC/DAB 
library. As such the Yee algorithm used is as simple as is practical to 
achieve this purpose. Furthermore, when it is an option, we have opted
for clarity rather than performance. We intend this to serve as a guide
to using the CRBC/DAB library.


What this program does
----------------------

For this example, we consider Maxwell's equations in a homogeneous, 
isotropic, dielectric material is given by

.. math::
  :nowrap:

  \begin{align}
    \frac{\partial \mathbf{H}}{\partial t} &= - \frac{1}{\mu}  \nabla \times \mathbf{E}, \\
    \frac{\partial \mathbf{E}}{\partial t}&= \frac{1}{\varepsilon} \nabla \times \mathbf{H},  
  \end{align}

subject to the constraints

.. math::
  :nowrap:

  \begin{align}
    \nabla \cdot \mathbf{E} &= 0, \\
    \nabla \cdot \mathbf{H} &= 0.
  \end{align}

.. _discretization:

To discretize these equations using Yee's scheme on a rectangular
domain :math:`[x_L, x_R] \times [y_L, y_R] \times [z_L, z_R]` with mesh spacings
of :math:`\Delta x`,  :math:`\Delta y`, and  :math:`\Delta z`, in the 
:math:`x`, :math:`y`, and :math:`z` directions, respectively, we define

.. math::
  :nowrap:

  \begin{align}
    x_i &= x_L + i\Delta x, \\
    y_j &= y_L + j\Delta y, \\
    z_k &= z_L + k\Delta z. 
  \end{align}

We choose a time step size, :math:`\Delta t`, satisfying

.. math::
  :nowrap:

  \begin{align}
    \Delta t \leq \frac{1}{c \sqrt{\frac{1}{(\Delta x)^2} + \frac{1}{(\Delta y)^2} + \frac{1}{(\Delta z)^2}}},
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
    E_x^{i+\frac{1}{2},j,k,n+\frac{1}{2}} & \sim E_x(x_{i+\frac{1}{2}},y_j,z_k,t_{n+\frac{1}{2}}), \\
    E_y^{i,j+\frac{1}{2},k,n+\frac{1}{2}} & \sim E_y(x_i,y_{j+\frac{1}{2}},z_k,t_{n+\frac{1}{2}}), \\
    E_z^{i,j,k+\frac{1}{2},n+\frac{1}{2}} & \sim E_z(x_i,y_j,z_{k+\frac{1}{2}},t_{n+\frac{1}{2}}), \\
    H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n} & \sim H_x(x_i,y_{j+\frac{1}{2}},z_{k+\frac{1}{2}},t_n), \\
    H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n} & \sim H_y(x_{i+\frac{1}{2}},y_j,z_{k+\frac{1}{2}},t_n), \\
    H_z^{i+\frac{1}{2},j+\frac{1}{2},k,n} & \sim H_z(x_{i+\frac{1}{2}},y_{j+\frac{1}{2}},z_k,t_n),
  \end{align}

where we require that the domain is terminated such that the tangential **E** 
components and the normal **H** component are located on the boundaries. This 
corresponds to having an integer number of Yee cells that match the 
illustrated :ref:`fig_yee_cell`

.. _fig_yee_cell:
.. figure:: cell.png
   :align: center
   :figwidth: 500 px
   :width: 450 px
   :alt: image of a Yee cell

   Spatial configuration of a Yee cell.

Finally, the fields are evolved with

.. math::
  :nowrap:

  \begin{align}
     H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n+1}  = H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n} 
     & +  \frac{\Delta t}{\mu \Delta z} \left(E_y^{i,j+\frac{1}{2},k+1,n+\frac{1}{2}} 
     - E_y^{i,j+\frac{1}{2},k,n+\frac{1}{2}} \right) \\
     & - \frac{\Delta t}{\mu \Delta y} \left(E_z^{i,j+1,k+\frac{1}{2},n+\frac{1}{2}} 
     - E_z^{i,j,k+\frac{1}{2},n+\frac{1}{2}}    \right),  \nonumber \\
     % % % %
     H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n+1}  = H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n}
     & + \frac{\Delta t}{\mu \Delta x} \left(E_z^{i+1,j,k+\frac{1}{2},n+\frac{1}{2}} 
     - E_z^{i,j,k+\frac{1}{2},n+\frac{1}{2}} \right) \\
     & -  \frac{\Delta t}{\mu \Delta z} \left(E_x^{i+\frac{1}{2},j,k+1,n+\frac{1}{2}} 
     - E_x^{i+\frac{1}{2},j,k,n+\frac{1}{2}} \right), \nonumber \\
     % % % %
     H_z^{i+\frac{1}{2},j+\frac{1}{2},k,n+1}  = H_z^{i+\frac{1}{2},j+\frac{1}{2},k,n} 
     & +  \frac{\Delta t}{\mu \Delta y} \left(E_x^{i+\frac{1}{2},j+1,k,n+\frac{1}{2}} 
     - E_x^{i+\frac{1}{2},j,k,n+\frac{1}{2}}\right)    \\
     & -  \frac{\Delta t}{\mu \Delta x} \left(E_y^{i+1,j+\frac{1}{2},k,n+\frac{1}{2}} 
     - E_y^{i,j+\frac{1}{2},k,n+\frac{1}{2}} \right), \nonumber \\
     % % % % %
     E_x^{i+\frac{1}{2},j,k,n+\frac{1}{2}}  = E_x^{i+\frac{1}{2},j,k,n-\frac{1}{2}}  
     & + \frac{\Delta t}{\varepsilon \Delta y} \left(H_z^{i+\frac{1}{2},j+\frac{1}{2},k,n}
     - H_z^{i+\frac{1}{2},j-\frac{1}{2},k,n} \right)  \\
     & - \frac{\Delta t}{\varepsilon \Delta z}  \left(H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n} 
     - H_y^{i+\frac{1}{2},j,k-\frac{1}{2},n}\right),  \nonumber \\
     % % % % %
     E_y^{i,j+\frac{1}{2},k,n+\frac{1}{2}}  = E_y^{i,j+\frac{1}{2},k,n-\frac{1}{2}} 
     & + \frac{\Delta t}{\varepsilon \Delta z}  \left(H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n} 
     - H_x^{i,j+\frac{1}{2},k-\frac{1}{2},n}\right) \\
     & - \frac{\Delta t}{\varepsilon \Delta x}  \left(H_z^{i+\frac{1}{2},j+\frac{1}{2},k,n} 
     - H_z^{i-\frac{1}{2},j+\frac{1}{2},k,n} \right),  \nonumber \\
     % % % % %
     E_z^{i,j,k+\frac{1}{2},n+\frac{1}{2}}  = E_z^{i,j,k+\frac{1}{2},n-\frac{1}{2}} 
     & + \frac{\Delta t}{\varepsilon \Delta x}  \left(H_y^{i+\frac{1}{2},j,k+\frac{1}{2},n} 
     - H_y^{i-\frac{1}{2},j,k+\frac{1}{2},n} \right)   \\
     & - \frac{\Delta t}{\varepsilon \Delta y}  \left(H_x^{i,j+\frac{1}{2},k+\frac{1}{2},n} 
     - H_y^{i,j-\frac{1}{2},k-\frac{1}{2},n} \right).  \nonumber
  \end{align}

We drive the simulation with a impulsive point source that takes the form of
a differentiated Gaussian. We place this source at the center of the computational
domain in the :math:`E_x` field and implement it as a soft source.


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

To get  the interface to the CRBC/DAB library specialized for the 3D Yee scheme 
(see `3d_yee_crbc_api.h <https://github.com/jrlagrone/rbcpack/blob/main/src/CRBC/3d_yee_crbc_api.h?at=default/>`_) ::

  #include <3d_yee_crbc_api.h>

Data structures
---------------

We define a struct to hold all of the data needed for the Yee scheme.
The only thing we do differently from what might be included in a typical Yee scheme
is to include an array of boundary types that the CRBC library requires. To do this we
utilize the type ``CRBC_Boundaries_t`` defined in the library interface. ::

  typedef struct yeeData {

    // constants (we'll use MKS units)
    double pi;   
    double C;     
    double mu0;   
    double eps0;  

    // relative permeability
    double epsR;

    // number of time steps 
    int ntsteps;    

    // grid sizes
    int imax;   
    int jmax;   
    int kmax;  

    // grid spacing in each direction
    double dx;   
    double dy;   
    double dz;  

    // time step size, (we'll use 0.99*CFL)  
    double dt;  

    // boundary conditions
    CRBC_Boundaries_t boundaries[6];

    // source parameters 
    double tw;        // pulse width   
    double t0;        // delay   
    double amp;       // Amplitude 

    // specify how often to generate output (in time steps) 
    int save_freq;   

    // H & E Field components  
    double ***Hx;   
    double ***Hy;   
    double ***Hz;   
    double ***Ex;   
    double ***Ey;   
    double ***Ez;   
       
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

Finally, a function to write Ex field output ::

  void writeExField(yeeData *d, int id); 

       
Main routine
------------

First we declare the variables we need. 
We also perform a primitive command line input check to enable the Ex field to
be written out to files periodically. This feature is off by default because it
generates approximately 350 MB of files.
Note the that ``YeeCrbcUpdater`` type is 
provided by the CRBC library and we need to declare it as a pointer because at 
this time its size is unknown. We will initialize it at a later point after 
we have set the simulation parameters. ::

  int main(int argc, char *argv[]) {   

    int tstep, i, j, k;
    double norm, enorm, hnorm;
    int write_output_files = 0;

    // read in input to see if we should write output files. If output files are
    // enabled this program writes out the Ex field 50 times. this generates about 
    // 350 MB of data.
    // By default, the file output is turned off.
    // There is only one option, so for simplicity we'll just assume that if we
    // receive any command line option, then we should enable output instead of
    // actually parsing and identifying a specific option.
    if (argc > 1) {
      printf("The Ex field will be saved. \n");
      write_output_files = 1;
    }
  
    // declare a yee data structure
    yeeData data;

    // declare a boundary updater object.
    // This needs to be a pointer because the size is unknown at this time and will
    // be initialized later.
    YeeCrbcUpdater *boundary_updater;

Next, we set the basic simulation parameters ::

  // relative permeability
  data.epsR = 1.0; // 1.0 corresponds to a vacuum  

  // number of time steps 
  data.ntsteps = 500;   

  // grid size 
  data.imax = 121;   
  data.jmax = 120;   
  data.kmax = 120;  

  // grid spacing in each direction
  data.dx = 1e-3;   
  data.dy = 1e-3;   
  data.dz = 1e-3;  

The boundary conditions are set using the ``CRBC_Boundaries_t`` 
enumeration type provided by the CRBC library.  At this time, 
the supported boundary types are 

 * `CRBC_PEC`  --- Perfect Electric Conductor
 * `CRBC_CRBC` --- Complete Radiation Boundary Condition (implemented as a DAB)

The library also provides enumerations that list the valid sides with the ``CRBC_Side_t``
type, they are

*  `CRBC_XLeft`
*  `CRBC_XRight`
*  `CRBC_YLeft`
*  `CRBC_YRight`
*  `CRBC_ZLeft`
*  `CRBC_ZRight`

Here, we'll set the boundary conditions so that we have parallel PEC plates with 
normals in the z-direction. ::

  // boundary conditions
  data.boundaries[CRBC_XLeft]  = CRBC_CRBC;
  data.boundaries[CRBC_XRight] = CRBC_CRBC;
  data.boundaries[CRBC_YLeft]  = CRBC_CRBC;
  data.boundaries[CRBC_YRight] = CRBC_CRBC;
  data.boundaries[CRBC_ZLeft]  = CRBC_PEC;
  data.boundaries[CRBC_ZRight] = CRBC_PEC;

Finally, we set the source parameters. For this example we will use an impulsive 
source that takes the form of a differentiated Gaussian. We will center in the 
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

      // loop over Ex
      for (i=0; i<data.imax-1; ++i)
        for (j=0; j<data.jmax; ++j)
          for (k=0; k<data.kmax; ++k)
            enorm += data.Ex[i][j][k] * data.Ex[i][j][k];

      // loop over Ey
      for (i=0; i<data.imax; ++i)
        for (j=0; j<data.jmax-1; ++j)
          for (k=0; k<data.kmax; ++k)
            enorm += data.Ey[i][j][k] * data.Ey[i][j][k];

      // loop over Ez
      for (i=0; i<data.imax; ++i)
        for (j=0; j<data.jmax; ++j)
          for (k=0; k<data.kmax-1; ++k)
            enorm += data.Ez[i][j][k] * data.Ez[i][j][k];

      // loop over Hx
      for (i=0; i<data.imax; ++i)
        for (j=0; j<data.jmax-1; ++j)
          for (k=0; k<data.kmax-1; ++k)
            hnorm += data.Hx[i][j][k] * data.Hx[i][j][k];

      // loop over Hy
      for (i=0; i<data.imax-1; ++i)
        for (j=0; j<data.jmax; ++j)
          for (k=0; k<data.kmax-1; ++k)
            hnorm += data.Hy[i][j][k] * data.Hy[i][j][k];

      // loop over Hz
      for (i=0; i<data.imax-1; ++i)
        for (j=0; j<data.jmax-1; ++j)
          for (k=0; k<data.kmax; ++k)
            hnorm += data.Hz[i][j][k] * data.Hz[i][j][k];

      // compute norm of combined fields
      norm = sqrt(data.epsilon*enorm + data.mu*hnorm);
      
      // write to screen
      printf("tstep = %i \t |E| = %6.5e \t |H| = %6.5e \t l2 norm = %16.5e \n", \
             tstep, sqrt(enorm), sqrt(hnorm), norm);

      // write all of the Ex field data to a file so we can visualize it if
      // requested.
      if (write_output_files != 0) {  
        writeExField(&data, tstep / data.save_freq);
      }

    } // end output

After doing the updates, we demonstrate the ability to save and restart the
boundary updater. This capability requires the HDF5 library. When building the
CRBC library, HDF5 is optional, so this is only enabled when this option is turned
on.

This writes the current state of the boundary library to several *h5* files. In
this case, it will save a basic parameter file called *hdf5_save_test.h5* and 
files for each of the three E-field components in the files `hdf5_save_test_ex.h5`,
`hdf5_save_test_ey.h5`, and `hdf5_save_test_ez.h5`. Depending on the boundary 
configuration, it may not save all of the boundary components. For demonstration
purposes, we delete the current boundary updater object to free all of the memory
associated to it and then create a new boundary updater object in its place using
the saved data half way through the simulation. :: 

    #if USE_HDF5_RESTARTS
    if (tstep == data.ntsteps /2) {

      // Save the current state
      if (CRBC_save_state(boundary_updater, "hdf5_save_test") != 0)
        fprintf(stderr, "save_stat() failed \n");

      // delete the updater object to ensure that everything is working
      // correctly ...
      CRBC_delete_updater(boundary_updater);

      // Now restart. This is an alternative way to generate a boundary updater
      // object if a previously saved state exists.
      if (CRBC_restart(&boundary_updater, "hdf5_save_test") != 0)
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
    (imax-1) \cdot (jmax) \cdot (kmax) & & E_x \text{ field values} \\
    (imax) \cdot (jmax-1) \cdot (kmax) & & E_y \text{ field values} \\
    (imax) \cdot (jmax) \cdot (kmax-1) & & E_z \text{ field values} \\
    (imax) \cdot (jmax-1) \cdot (kmax-1) & & H_x \text{ field values} \\
    (imax-1) \cdot (jmax) \cdot (kmax-1) & & H_y \text{ field values} \\
    (imax-1) \cdot (jmax-1) \cdot (kmax) & & H_z \text{ field values} 
  \end{align}

Additionally, we compute the constants that are used in the field update equations. ::

  void initialize(yeeData *d) {

    int i, j, k;
    int imax = d->imax;
    int jmax = d->jmax;
    int kmax = d->kmax;

    d->flag = 0;

    d->pi = 3.14159265358979;   
    d->C = 2.99792458e8;   

    // time step size, (we'll use 0.99*CFL)  
    d->dt = 0.99/(d->C*sqrt(1.0/(d->dx*d->dx) + 1.0/(d->dy*d->dy) + 1.0/(d->dz*d->dz))); 

    // calculate free space eps and mu
    d->mu0 = 4.0 * d->pi * 1.0E-7;   
    d->eps0 = 1.0 / (d->C * d->C * d->mu0);

    // calculate material epsilon and mu
    d->epsilon = d->epsR * d->eps0;
    d->mu = d->mu0;

    // calculate update coefficients
    d->CE = d->dt / d->epsilon;
    d->CH = d->dt / d->mu0;

    printf("CE = %e \n", d->CE);
    printf("CH = %e \n", d->CH);

    // allocate memory for the fields
    d->Ex = (double ***) malloc((imax-1) * sizeof(double));
    for (i=0; i < d->imax-1; i++) {
      d->Ex[i] = (double **) malloc((jmax) * sizeof(double));
      for (j=0; j < d->jmax; j++) {
        d->Ex[i][j] = (double *) malloc((kmax) * sizeof(double));
        for (k=0; k < d->kmax; k++) {
          d->Ex[i][j][k] = 0.0;
        }
      }
    }

    d->Ey = (double ***) malloc((imax) * sizeof(double));
    for (i=0; i < d->imax; i++) {
      d->Ey[i] = (double **) malloc((jmax-1) * sizeof(double));
      for (j=0; j < d->jmax-1; j++) {
        d->Ey[i][j] = (double *) malloc(kmax * sizeof(double));
        for (k=0; k < d->kmax; k++) {
          d->Ey[i][j][k] = 0.0;
        }
      }
    }

    d->Ez = (double ***) malloc((imax) * sizeof(double));
    for (i=0; i < d->imax; i++) {
      d->Ez[i] = (double **) malloc((jmax) * sizeof(double));
      for (j=0; j < d->jmax; j++) {
        d->Ez[i][j] = (double *) malloc((kmax-1) * sizeof(double));
        for (k=0; k < d->kmax-1; k++) {
          d->Ez[i][j][k] = 0.0;
        }
      }
    }

    d->Hx = (double ***) malloc((imax) * sizeof(double));
    for (i=0; i < d->imax; i++) {
      d->Hx[i] = (double **) malloc((jmax-1) * sizeof(double));
      for (j=0; j < d->jmax-1; j++) {
        d->Hx[i][j] = (double *) malloc((kmax-1) * sizeof(double));
        for (k=0; k < d->kmax-1; k++) {
          d->Hx[i][j][k] = 0.0;
        }
      }
    }

    d->Hy = (double ***) malloc((imax-1) * sizeof(double));
    for (i=0; i < d->imax-1; i++) {
      d->Hy[i] = (double **) malloc((jmax) * sizeof(double));
      for (j=0; j < d->jmax; j++) {
        d->Hy[i][j] = (double *) malloc((kmax-1) * sizeof(double));
        for (k=0; k < d->kmax-1; k++) {
          d->Hy[i][j][k] = 0.0;
        }
      }
    }

    d->Hz = (double ***) malloc((imax-1) * sizeof(double));
    for (i=0; i < d->imax-1; i++) {
      d->Hz[i] = (double **) malloc((jmax-1) * sizeof(double));
      for (j=0; j < d->jmax-1; j++) {
        d->Hz[i][j] = (double *) malloc((kmax) * sizeof(double));
        for (k=0; k < d->kmax; k++) {
          d->Hz[i][j][k] = 0.0;
        }
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

      // delete memory for the fields
      for (i=0; i < d->imax-1; i++) {
        for (j=0; j < d->jmax; j++) {
          free(d->Ex[i][j]);
        }
        free(d->Ex[i]);
      }
      free(d->Ex);

      for (i=0; i < d->imax; i++) {
        for (j=0; j < d->jmax-1; j++) {
          free(d->Ey[i][j]);
        }
        free(d->Ey[i]);
      }
      free(d->Ey);

      for (i=0; i < d->imax; i++) {
        for (j=0; j < d->jmax; j++) {
          free(d->Ez[i][j]);
        }
        free(d->Ez[i]);
      }
      free(d->Ez);

      for (i=0; i < d->imax; i++) {
        for (j=0; j < d->jmax-1; j++) {
          free(d->Hx[i][j]);
        }
        free(d->Hx[i]);
      }
      free(d->Hx);

      for (i=0; i < d->imax-1; i++) {
        for (j=0; j < d->jmax; j++) {
          free(d->Hy[i][j]);
        }
        free(d->Hy[i]);
      }
      free(d->Hy);

      for (i=0; i < d->imax-1; i++) {
        for (j=0; j < d->jmax-1; j++) {
          free(d->Hz[i][j]);
        }
        free(d->Hz[i]);
      }
      free(d->Hz);

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

  void setup_crbc(yeeData *d, YeeCrbcUpdater **upd) {

    int i, l;
    double T;
    double h[3];
    double delta;
    int low_index[3];
    int high_index[3];
    int n;
    CRBC_Fields_t components[3];
    CRBC_Side_t side;
  

    // First we will create a new boundary updater.
    h[0] = d->dx;
    h[1] = d->dy;
    h[2] = d->dz;
    T = d->dt * d->ntsteps; // total time is time step size * number of time steps
    *upd = CRBC_new_updater(T, h, d->dt, d->C, d->boundaries);

    // alternatively one can call 
    // int P = 7;
    // upd = CRBC_new_updater_p(T, h, d->dt, d->C, d->boundaries, P);
    // This does exactly the same thing but changes the number of recursions from
    // the default of 5 to 7.
    //
    // or 
    // int Pmax = 15;
    // double tol = 1e-3;
    // upd = CRBC_new_updater_tol(T, h, d->dt, d->C, d->boundaries, Pmax, tol);
    // Will generate an updater that attempts to satsify the provided tolerance
    // and with fewer than Pmax recursions.

Before we precede, we'll print out the properties from the boundary updater to 
make sure they are correct. ::

  printf("The boundary updater was initialized with the following: \n");
  printf("  wave speed, c = %e \n", CRBC_get_c(*upd));
  printf("  total time, T = %e \n", CRBC_get_T(*upd));
  printf("  time step, dt = %e \n", CRBC_get_dt(*upd));
  printf("The boundary updater calculated the following: \n");
  printf("  %i Ex faces \n", CRBC_get_num_faces(*upd, CRBC_Ex));
  printf("  %i Ey faces \n", CRBC_get_num_faces(*upd, CRBC_Ey));
  printf("  %i Ez faces \n", CRBC_get_num_faces(*upd, CRBC_Ez));
  printf("  %i Ex edges \n", CRBC_get_num_edges(*upd, CRBC_Ex));
  printf("  %i Ey edges \n", CRBC_get_num_edges(*upd, CRBC_Ey));
  printf("  %i Ez edges \n", CRBC_get_num_edges(*upd, CRBC_Ez));
  printf("  %i Ex corners \n", CRBC_get_num_corners(*upd, CRBC_Ex));
  printf("  %i Ey corners \n", CRBC_get_num_corners(*upd, CRBC_Ey));
  printf("  %i Ez corners \n", CRBC_get_num_corners(*upd, CRBC_Ez));

Now set up the faces. Start by looping over all of the possible faces. We
follow the order given in ``CRBC_Side_t``, so 

* l = 0 --> CRBC_XLeft
* l = 1 --> CRBC_XRight
* l = 2 --> CRBC_YLeft
* l = 3 --> CRBC_YRight
* l = 4 --> CRBC_ZLeft
* l = 5 --> CRBC_ZRight

First we need to calculate the minimum distance between the boundary and
source, `delta`. Since we placed the source in the center of the domain, this is 
simply the distance from the boundary to the center of the domain in the
direction normal to the boundary.  In general, if it is not possible to calculate
delta directly, using a lower bound for the separation is the safest thing to do,
but it may result in more work being done that is actually needed to achieve the
desired accuracy. ::

  for (l=0; l<6; ++l) {

    // calculate delta
    switch (l / 2) {
   
      case 0: // Faces with a normal in the x-direction
        delta = (d->imax / 2) * (d->dx);
        break;
    
      case 1: // faces with a normal in the y-direction
        delta = (d->jmax / 2) * (d->dy);
        break;

      case 2: // faces with a normal in the z-direction
        delta = (d->kmax / 2) * (d->dz);
        break;
    }

    // convert the index l to a CRBC_Side_t type
    switch (l) {
      case 0:
        side = CRBC_XLeft;
        break;
      case 1:
        side = CRBC_XRight;
        break;
      case 2:
        side = CRBC_YLeft;
        break;
      case 3:
        side = CRBC_YRight;
        break;
      case 4:
        side = CRBC_ZLeft;
        break;
      case 5:
        side = CRBC_ZRight;
        break;
    }


Next, we get the field components that the updater expects and set up each 
of these components on the current boundary face. The CRBC updater library 
attempts to communicate using the solvers "native" indexing scheme, so we need to
tell it the upper and lower indexing extents. These extents need to include all 
of the requested component variables on the boundary face (or just inside the 
boundary in the case of a normal component) as well the adjacent parallel plane 
of points in the interior of the domain. The boundary updater considers these
extents to be inclusive.

In particular, for the left boundary in the x-direction, this means we need to 
include the data points of the requested component at i=0 and i=1 and all possible 
indices for j and k. Therefore, the lower indexing extents are 0 for all indices.
The upper extents vary by component due to the staggered grid.

For the right boundary in the x-direction, this means we need to include the data
points of the requested component at i=imax-1 and i=imax-2 or 
i=imax-2 and i=imax-3 depending on the component along with all of 
the possible indices for j and k.

For the left boundary in the y-direction, this means we need to include the data 
points of the requested component at j=0 and j=1 along with all of 
the possible indices for i and k.

For the right boundary in the y-direction, this means we need to include the data
points of the requested component at j=jmax-1 and j=jmax-2 or 
j=jmax-2 and j=jmax-3 depending on the component along with all of 
the possible indices for i and k.

For the left boundary in the z-direction, this means we need to include the data 
points of the requested component at k=0 and k=1 along with all of 
the possible indices for j and k. 

For the right boundary in the z-direction, this means we need to include the data
points of the requested component at k=kmax-1 and k=kmax-2 or 
k=kmax-2 and k=kmax-3 depending on the component along with all of 
the possible indices for i and j. ::

    // get the field components that the updater expects. If this
    // is not a boundary that will be handled by the updater, this should 
    // return n = 0.
    CRBC_get_component_inputs(*upd, side, comps, &n);

    // Set up each component the updater has requested.
    for (i=0; i<n; ++i) {

      // Set the basic extents and then adjust them based on the side and field 
      // component. These should be inclusive extents;
      low_index[0] = 0;
      low_index[1] = 0;
      low_index[2] = 0;

      high_index[0] = d->imax-1; // minus 1 from 0-based indexing
      high_index[1] = d->jmax-1;
      high_index[2] = d->kmax-1; 

      switch (side) {

        case CRBC_XLeft:
          high_index[0] = 1;  
    
          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ey) {
            high_index[1]--;
          } else if (comps[i] == CRBC_Ez) {
            high_index[2]--;
          }

          break;

        case CRBC_XRight:
          low_index[0] = d->imax-2;

          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ex) {
            low_index[0]--;
            high_index[0]--;
          } else if (comps[i] == CRBC_Ey) {
            high_index[1]--;
          } else if (comps[i] == CRBC_Ez) {
            high_index[2]--;
          }

          break;

        case CRBC_YLeft:
          high_index[1] = 1;
    
          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ex) {
            high_index[0]--;
          } else if (comps[i] == CRBC_Ez) {
            high_index[2]--;
          }
          
          break;  

        case CRBC_YRight:
          low_index[1] = d->jmax-2;

          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ex) {
            high_index[0]--;
          } else if (comps[i] == CRBC_Ey) {
            low_index[1]--;
            high_index[1]--;
          } else if (comps[i] == CRBC_Ez) {
            high_index[2]--;
          }
     
          break;

        case CRBC_ZLeft:
          high_index[2] = 1;
    
          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ex) {
            high_index[0]--;
          } else if (comps[i] == CRBC_Ey) {
            high_index[1]--;
          }

          break;

        case CRBC_ZRight:
          low_index[2] = d->kmax-2;
    
          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ex) {
            high_index[0]--;
          } else if (comps[i] == CRBC_Ey) {
            high_index[1]--;
          } else if (comps[i] == CRBC_Ez) {
            low_index[2]--;
            high_index[2]--;
          } 

          break;

      } // end switch side

Now, we can initialize the face: ::

      // now initialize the face.
      if (CRBC_init_face(*upd, side, comps[i], low_index, high_index, delta) != 0)
      {
        fprintf(stderr, "Error: init_face(...) failed \n");
        exit(-1);
      }

    }  // end for over requested components (i)
  } // end for over possible boundary faces (l)

Finally, we print out some information about the recursions. Note that a reflection
coefficient of -1 indicates that we are not performing any updates on that face. ::

    // now we'll print out some information about the recursions
    printf("The faces were initialized with: \n");
    printf("  Left side in  x is using %i recursions "
           "with a reflection coefficient of %e \n"
           , CRBC_get_num_recursions(*upd, CRBC_XLeft), 
            CRBC_get_reflection_coef(*upd, CRBC_XLeft));
    printf("  Right side in x is using %i recursions "
           "with a reflection coefficient of %e \n"
           , CRBC_get_num_recursions(*upd, CRBC_XRight), 
           CRBC_get_reflection_coef(*upd, CRBC_XRight));
    printf("  Left side in  y is using %i recursions "
           "with a reflection coefficient of %e \n"
           , CRBC_get_num_recursions(*upd, CRBC_YLeft), 
           CRBC_get_reflection_coef(*upd, CRBC_YLeft));
    printf("  Right side in y is using %i recursions "
           "with a reflection coefficient of %e \n" 
           , CRBC_get_num_recursions(*upd, CRBC_YRight), 
           CRBC_get_reflection_coef(*upd, CRBC_YRight));
    printf("  Left side in  z is using %i recursions "
           "with a reflection coefficient of %e \n" 
           , CRBC_get_num_recursions(*upd, CRBC_ZLeft), 
           CRBC_get_reflection_coef(*upd, CRBC_ZLeft));
    printf("  Right side in z is using %i recursions "
           "with a reflection coefficient of %e \n" 
           , CRBC_get_num_recursions(*upd, CRBC_ZRight), 
           CRBC_get_reflection_coef(*upd, CRBC_ZRight));
    printf("The maximum reflection coefficient is %e \n", 
           CRBC_get_max_reflection_coef(*upd));

  } // end setup_crbc     


computeE function
-----------------

This function computes the updates to the **E** field components using the Yee scheme. ::

  void computeE(yeeData *d, int *tstep) {

    int i, j, k;

    // compute updates to Ex
    for (i=0; i < d->imax-1; ++i) {
      for (j=1; j < d->jmax-1; ++j) {
        for (k=1; k < d->kmax-1; ++k) {
           d->Ex[i][j][k] += d->CE * ((d->Hz[i][j][k] - d->Hz[i][j-1][k]) / d->dy \
                        - (d->Hy[i][j][k] - d->Hy[i][j][k-1]) / d->dz);
        }
      }
    }

    // compute updates to Ey
    for (i=1; i < d->imax-1; ++i) {
      for (j=0; j < d->jmax-1; ++j) {
        for (k=1; k < d->kmax-1; ++k) {
           d->Ey[i][j][k] += d->CE * ((d->Hx[i][j][k] - d->Hx[i][j][k-1]) / d->dz \
                          - (d->Hz[i][j][k] - d->Hz[i-1][j][k]) / d->dx);
        }
      }
    }

    // compute updates to Ez
    for (i=1; i < d->imax-1; ++i) {
      for (j=1; j < d->jmax-1; ++j) {
        for (k=0; k < d->kmax-1; ++k) {
           d->Ez[i][j][k] += d->CE * ((d->Hy[i][j][k] - d->Hy[i-1][j][k]) / d->dx \
                        - (d->Hx[i][j][k] - d->Hx[i][j-1][k]) / d->dy);
        }
      }
    }

    // add the source term to the center point of the Ex field
    // This is a differentiated Gaussian applied as a "soft" source
    d->Ex[d->imax/2][d->jmax/2][d->kmax/2] += 2.0 * d->amp * d->CE \
                                * ((*tstep*d->dt - d->tO) / d->tw) \
                                * exp(-pow(((*tstep*d->dt - d->tO) / d->tw), 2.0));

  }

computeH function
-----------------

This function computes the updates to the **H** field components using the Yee scheme. ::

  void computeH(yeeData *d) {

    int i, j, k;

    // compute updates to Hx
    for (i=0; i < d->imax; ++i) {
      for (j=0; j < d->jmax-1; ++j) {
        for (k=0; k < d->kmax-1; ++k) {
           d->Hx[i][j][k] += d->CH * ((d->Ey[i][j][k+1] - d->Ey[i][j][k]) / d->dz \
                        - (d->Ez[i][j+1][k] - d->Ez[i][j][k]) / d->dy);
        }
      }
    }

    // compute updates to Hy
    for (i=0; i < d->imax-1; ++i) {
      for (j=0; j < d->jmax; ++j) {
        for (k=0; k < d->kmax-1; ++k) {
           d->Hy[i][j][k] += d->CH * ((d->Ez[i+1][j][k] - d->Ez[i][j][k]) / d->dx \
                        - (d->Ex[i][j][k+1] - d->Ex[i][j][k]) / d->dz);
        }
      }
    }

    // compute updates to Hz
    for (i=0; i < d->imax-1; ++i) {
      for (j=0; j < d->jmax-1; ++j) {
        for (k=0; k < d->kmax; ++k) {
           d->Hz[i][j][k] += d->CH * ((d->Ex[i][j+1][k] - d->Ex[i][j][k]) / d->dy \
                        - (d->Ey[i+1][j][k] - d->Ey[i][j][k]) / d->dx);
        }
      }
    }

  }


computeBoundary Function
------------------------

This function computes the boundary updates using the CRBC library. To do this, 
we first need to copy the values that have been updated by the Yee algorithm into
the CRBC updater. We start by looping over all of the possible sides and getting
the list of components that the CRBC updater requests as input. For each of these
components, the CRBC updater can return the index extents, which we use to copy 
the data into the CRBC library. ::

  void computeBoundary(yeeData *d, YeeCrbcUpdater *upd) {

    int i, j, k, l, m;
    double ***field;
    int low_index[3];
    int high_index[3];
    int index[3];
    int n;
    CRBC_Fields_t components[3];
    CRBC_Side_t side;

    // first we need to copy the values that have been updated by the Yee algorithm
    // into the crbc updater. First we will loop over all of the possible sides.
    for (m=0; m<6; ++m) {

      // convert index to CRBC_Side type
      switch (m) {
        case 0:
          side = CRBC_XLeft;
          break;
        case 1:
          side = CRBC_XRight;
          break;
        case 2:
          side = CRBC_YLeft;
          break;
        case 3:
          side = CRBC_YRight;
          break;
        case 4:
          side = CRBC_ZLeft;
          break;
        case 5:
          side = CRBC_ZRight;
          break;
      }

      // Get the field components that the updater expects.
      CRBC_get_component_inputs(upd, side, components, &n);

      // loop over the components on this face
      for (l=0; l<n; l++) {

        /// get the array pointer for the correct field
        switch (components[l]) {
          case CRBC_Ex:
            field = d->Ex;
            break;
          case CRBC_Ey:
            field = d->Ey;
            break;
          case CRBC_Ez:
            field = d->Ez;
            break;
        }
 
        // get the indices the updater expects as input.
        // These indices are inclusive.
        CRBC_get_input_extents(upd, side, components[l], low_index, high_index);

        // copy data into the updater
        for (i=low_index[0]; i<=high_index[0]; ++i) {
          index[0] = i;
          for (j=low_index[1]; j<=high_index[1]; ++j) {
            index[1] = j;
            for (k=low_index[2]; k<=high_index[2]; ++k) {
              index[2] = k;
              CRBC_load_face_data(upd, side, components[l], index, &(field[i][j][k]));
            }
          }
        }

      } // end loop over sides
    } // end loop over sides

Now we can have the CRBC library compute the boundary updates. ::

  // now we can compute the updates to the boundaries
  if (CRBC_compute_updates(upd) != 0) {
    fprintf(stderr, "Error: compute_updates(...) failed \n");
    exit(-1);
  }

Finally, we need to copy the new values from the boundary updater into the arrays 
used by the Yee algorithm. To do this, we loop over all the possible sides and get
the components that are available from the boundary updater. We skip the normal 
components because the Yee algorithm is able to update these values on its own. 
Then, we get the output index extents from the boundary updater and use these 
values to loop over the data arrays and copy the values. ::

    // Now copy the new boundary values into the arrays used by the Yee algorithm
    // loop over possible sides
    for (m=0; m<6; ++m) {

      // convert index to CRBC_Side type
      switch (m) {
        case 0:
          side = CRBC_XLeft;
          break;
        case 1:
          side = CRBC_XRight;
          break;
        case 2:
          side = CRBC_YLeft;
          break;
        case 3:
          side = CRBC_YRight;
          break;
        case 4:
          side = CRBC_ZLeft;
          break;
        case 5:
          side = CRBC_ZRight;
          break;
      }

      // Get the field components that the updater expects.
      CRBC_get_component_inputs(upd, side, components, &n);

      // loop over the components on this face
      for (l=0; l<n; l++) {

        // If there is a normal component, the Yee solver should have already
        // updated all of the values correctly. Therefore, we don't have to
        // copy the new values, so we'll just skip them.
        if (side / 2 == components[l])
          continue;

        /// get the array pointer for the correct field
        switch (components[l]) {
          case CRBC_Ex:
            field = d->Ex;
            break;
          case CRBC_Ey:
            field = d->Ey;
            break;
          case CRBC_Ez:
            field = d->Ez;
            break;
        }

        // get the indices the updater can output.
        // These indices are inclusive.
        CRBC_get_output_extents(upd, side, components[l], low_index, high_index);

        // copy data into the updater
        for (i=low_index[0]; i<=high_index[0]; ++i) {
          index[0] = i;
          for (j=low_index[1]; j<=high_index[1]; ++j) {
            index[1] = j;
            for (k=low_index[2]; k<=high_index[2]; ++k) {
              index[2] = k;
              field[i][j][k] = CRBC_get_new_face_vals(upd, side, components[l], index);
            }
          }
        }

      } // end loop over sides
    } // end loop over sides

  } // end computeBoundary

writeExField function
---------------------

We use this function to output the :math:`E_x` field component to an ASCII VTK 
file format that can be opened with visualization software such as ParaView. ::

  void writeExField(yeeData *d, int id) {

    int i, j, k, l, n, cells;
    int extent[6] = {0, d->imax-2, 0, d->jmax-1, 0, d->kmax-1};

    char step[10];   
    char fileBaseName[] = "Ex_Field_";   
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
    fprintf(f, "DIMENSIONS %i %i %i\n", d->imax-1, d->jmax, d->kmax);

    // save the coordinates
    fprintf(f, "X_COORDINATES %i double \n", d->imax-1);
    for (i=0; i < d->imax-1; ++i)
      fprintf(f, "%f\n", d->dx/2.0 + i*d->dx);

    fprintf(f, "Y_COORDINATES %i double \n", d->jmax);
    for (j=0; j < d->jmax; ++j)
      fprintf(f, "%f\n", j*d->dy);

    fprintf(f, "Z_COORDINATES %i double \n", d->kmax);
    for (k=0; k < d->kmax; ++k)
      fprintf(f, "%f\n", k*d->dz);

    // set up a cell and field
    n = (d->imax-1) * d->jmax * d->kmax;
    cells = (d->imax-2) * (d->jmax-1) * (d->kmax-1);
    fprintf(f, "CELL_DATA %i\n", cells);
    fprintf(f, "POINT_DATA %i\n", n);
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "Ex 1 %i double\n", n);

    // now write out the data
    for (k=0; k < d->kmax; ++k) 
      for (j=0; j < d->jmax; ++j)
        for (i=0; i < d->imax-1; ++i)
          fprintf(f, "%f\n", d->Ex[i][j][k]);
  	

    // close the file
    fclose(f);
  }


Output
======

A visualization of the :math:`E_x` field
is available for reference. To better illustrate the behavior, we cut out a corner
and set the color scale to be reflective of the state roughly midway through the
simulation. 

.. raw:: html

  <div style="text-align: center">
    <iframe width="600" height="450" src="https://www.youtube.com/embed/4ngPbgkyyrY" frameborder="0" allowfullscreen></iframe>
  </div>

The files used to generate this movie can be generated by providing any command line
option at runtime, for instance ::

  ./yee.x -output

The following output is generated to the screen. ::

  The boundary updater was initialized with the following: 
    wave speed, c = 2.997925e+08 
    total time, T = 9.532874e-10 
    time step, dt = 1.906575e-12 
  The boundary updater calculated the following: 
    4 Ex faces 
    4 Ey faces 
    4 Ez faces 
    4 Ex edges 
    4 Ey edges 
    4 Ez edges 
    0 Ex corners 
    0 Ey corners 
    0 Ez corners 
  The faces were initialized with: 
    Left side in  x is using 5 recursions with a reflection coefficient of 1.831672e-05 
    Right side in x is using 5 recursions with a reflection coefficient of 1.831672e-05 
    Left side in  y is using 5 recursions with a reflection coefficient of 1.831672e-05 
    Right side in y is using 5 recursions with a reflection coefficient of 1.831672e-05 
    Left side in  z is using 5 recursions with a reflection coefficient of -1.000000e+00 
    Right side in z is using 5 recursions with a reflection coefficient of -1.000000e+00 
  The maximum reflection coefficient is 1.831672e-05 
  tstep = 10 	 |E| = 7.28158e-03 	 |H| = 5.36599e-06 	 l2 norm = 2.24865e-08 
  tstep = 20 	 |E| = 1.03198e-01 	 |H| = 6.66446e-05 	 l2 norm = 3.16034e-07 
  tstep = 30 	 |E| = 1.04193e+00 	 |H| = 6.10638e-04 	 l2 norm = 3.17502e-06 
  tstep = 40 	 |E| = 7.84034e+00 	 |H| = 4.10766e-03 	 l2 norm = 2.37798e-05 
  tstep = 50 	 |E| = 4.41031e+01 	 |H| = 2.01666e-02 	 l2 norm = 1.33166e-04 
  tstep = 60 	 |E| = 1.85500e+02 	 |H| = 7.15808e-02 	 l2 norm = 5.57775e-04 
  tstep = 70 	 |E| = 5.83387e+02 	 |H| = 1.80684e-01 	 l2 norm = 1.74770e-03 
  tstep = 80 	 |E| = 1.37182e+03 	 |H| = 3.13817e-01 	 l2 norm = 4.09711e-03 
  tstep = 90 	 |E| = 2.41184e+03 	 |H| = 3.45411e-01 	 l2 norm = 7.18710e-03 
  tstep = 100 	 |E| = 3.17028e+03 	 |H| = 1.75511e-01 	 l2 norm = 9.43554e-03 
  tstep = 110 	 |E| = 3.11557e+03 	 |H| = 1.67893e-01 	 l2 norm = 9.27259e-03 
  tstep = 120 	 |E| = 2.28914e+03 	 |H| = 3.51954e-01 	 l2 norm = 6.82297e-03 
  tstep = 130 	 |E| = 1.25778e+03 	 |H| = 3.37815e-01 	 l2 norm = 3.76175e-03 
  tstep = 140 	 |E| = 5.17905e+02 	 |H| = 2.21236e-01 	 l2 norm = 1.56090e-03 
  tstep = 150 	 |E| = 1.64024e+02 	 |H| = 1.39904e-01 	 l2 norm = 5.12647e-04 
  tstep = 160 	 |E| = 5.64689e+01 	 |H| = 1.15573e-01 	 l2 norm = 2.12176e-04 
  tstep = 170 	 |E| = 4.51007e+01 	 |H| = 1.05073e-01 	 l2 norm = 1.78560e-04 
  tstep = 180 	 |E| = 4.18910e+01 	 |H| = 1.08148e-01 	 l2 norm = 1.73883e-04 
  tstep = 190 	 |E| = 3.73272e+01 	 |H| = 1.13062e-01 	 l2 norm = 1.68524e-04 
  tstep = 200 	 |E| = 3.61365e+01 	 |H| = 1.02364e-01 	 l2 norm = 1.57257e-04 
  tstep = 210 	 |E| = 3.41130e+01 	 |H| = 8.90912e-02 	 l2 norm = 1.42400e-04 
  tstep = 220 	 |E| = 3.14542e+01 	 |H| = 7.98402e-02 	 l2 norm = 1.29501e-04 
  tstep = 230 	 |E| = 2.95496e+01 	 |H| = 7.30628e-02 	 l2 norm = 1.20164e-04 
  tstep = 240 	 |E| = 2.54773e+01 	 |H| = 7.28468e-02 	 l2 norm = 1.11426e-04 
  tstep = 250 	 |E| = 2.16452e+01 	 |H| = 7.15643e-02 	 l2 norm = 1.02879e-04 
  tstep = 260 	 |E| = 2.36933e+01 	 |H| = 5.73692e-02 	 l2 norm = 9.54272e-05 
  tstep = 270 	 |E| = 2.43793e+01 	 |H| = 4.58221e-02 	 l2 norm = 8.88876e-05 
  tstep = 280 	 |E| = 2.00561e+01 	 |H| = 5.20351e-02 	 l2 norm = 8.34511e-05 
  tstep = 290 	 |E| = 1.62060e+01 	 |H| = 5.54434e-02 	 l2 norm = 7.86657e-05 
  tstep = 300 	 |E| = 1.60938e+01 	 |H| = 5.10981e-02 	 l2 norm = 7.46621e-05 
  tstep = 310 	 |E| = 1.73192e+01 	 |H| = 4.40709e-02 	 l2 norm = 7.13901e-05 
  tstep = 320 	 |E| = 1.87057e+01 	 |H| = 3.48466e-02 	 l2 norm = 6.80002e-05 
  tstep = 330 	 |E| = 1.61399e+01 	 |H| = 3.79793e-02 	 l2 norm = 6.41802e-05 
  tstep = 340 	 |E| = 9.81647e+00 	 |H| = 4.83870e-02 	 l2 norm = 6.16067e-05 
  tstep = 350 	 |E| = 1.24888e+01 	 |H| = 4.18796e-02 	 l2 norm = 5.98749e-05 
  tstep = 360 	 |E| = 1.58622e+01 	 |H| = 2.88812e-02 	 l2 norm = 5.72361e-05 
  tstep = 370 	 |E| = 1.44103e+01 	 |H| = 3.02477e-02 	 l2 norm = 5.46660e-05 
  tstep = 380 	 |E| = 1.25167e+01 	 |H| = 3.31844e-02 	 l2 norm = 5.26400e-05 
  tstep = 390 	 |E| = 1.03637e+01 	 |H| = 3.58635e-02 	 l2 norm = 5.06683e-05 
  tstep = 400 	 |E| = 9.01898e+00 	 |H| = 3.69629e-02 	 l2 norm = 4.93670e-05 
  tstep = 410 	 |E| = 1.33182e+01 	 |H| = 2.43437e-02 	 l2 norm = 4.81166e-05 
  tstep = 420 	 |E| = 1.38390e+01 	 |H| = 1.76370e-02 	 l2 norm = 4.56796e-05 
  tstep = 430 	 |E| = 8.81529e+00 	 |H| = 3.16989e-02 	 l2 norm = 4.41672e-05 
  tstep = 440 	 |E| = 8.05623e+00 	 |H| = 3.23147e-02 	 l2 norm = 4.34384e-05 
  tstep = 450 	 |E| = 9.76612e+00 	 |H| = 2.72060e-02 	 l2 norm = 4.21260e-05 
  tstep = 460 	 |E| = 1.01791e+01 	 |H| = 2.45823e-02 	 l2 norm = 4.09485e-05 
  tstep = 470 	 |E| = 1.13321e+01 	 |H| = 1.88892e-02 	 l2 norm = 3.98171e-05 
  tstep = 480 	 |E| = 9.61705e+00 	 |H| = 2.25214e-02 	 l2 norm = 3.81614e-05 
  tstep = 490 	 |E| = 4.11735e+00 	 |H| = 3.16045e-02 	 l2 norm = 3.74872e-05 
  tstep = 500 	 |E| = 8.10059e+00 	 |H| = 2.52941e-02 	 l2 norm = 3.72155e-05 

.. rubric:: References

.. bibliography:: zcite.bib
   :encoding: UTF8
   :list: enumerated
   :filter: title % "Yee Cube"
