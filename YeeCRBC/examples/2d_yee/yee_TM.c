/**
    2-D FDTD code with CRBC/DAB boundary conditions

    This C code implements the finite difference time-domain solution of
    Maxwell's equations (curl formulation) for a transverse-magnetic problem.
    The grid is terminated using Double Absorbing Boundaries (DAB).

    This is meant to serve as an example of how one might use the CRBC/DAB 
    library. As such the Yee algorithm used is as simple as is practical to 
    achieve this purpose. Furthermore, when it is an option, we have opted
    for clarity rather than performance. We intend this to serve as a guide
    to using the CRBC/DAB library.
    
    This code is freely distributed with the goal of illustrating the CRBC/DAB
    library usage and research purposes. We provide no warranties of its 
    performance or suitability for any purpose.

*/ 

// for output, namely fprintf and printf
#include <stdio.h>   
// definitions of some useful functions: sqrt(), fabs(), etc
#include <math.h>   
// store and modify output filenames
#include <string.h>   
// standard libraries, we need this for malloc, free, etc.
#include <stdlib.h>   
// the 2D interface to the CRBC/DAB library
#include <2d_crbc_api.h>

#if USE_OPENMP
  #include <omp.h>
#endif

/* -----------------------------------------------------------------------------

                     Define a struct to hold data

------------------------------------------------------------------------------*/       

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
  // NOTE: Field Value Locations
  //
  // Ez - at (i, j)
  //  
  // Hx - at (i, j + 1/2)
  // Hy - at (i + 1/2, j)
  //
  // This means we have
  // (imax) * (jmax)   Ez field values 
  // (imax) * (jmax-1) Hx field values 
  // (imax-1) * (jmax) Hy field values 
  int imax;   
  int jmax;   

  // grid spacing in each direction
  double dx;   
  double dy;   

  // time step size, (we'll use 0.99*CFL)  
  double dt;  

  // boundary conditions
  // The type CRBC2d_Boundaries is an provided by the CRBC/DAB library and it is 
  // simply an enumeration of the supported boundary types. At this time, the 
  // supported boundary types are 
  //     CRBC2d_DIR  --- Homogeneous Dirichlet
  //     CRBC2d_NEUM --- Homogeneous Neumman
  //     CRBC2d_CRBC --- Complete Radiation Boundary Condition (implemented as a DAB)
  CRBC2d_Boundaries_t boundaries[4];

  // source parameters 
  // we will use an impulsive source that takes the form of a 
  // differentiated gaussian centered in the domain    
  double tw;      // pulse width   
  double t0;      // delay   
  double amp;     // Amplitude 

  // specify how often to generate output (in time steps) 
  int save_freq;   

  // H & E Field components  
  // NOTE: storing these as a 1D array and accessing appropriately is likely to 
  //       be more efficient, but we do this for better clarity. 
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
       

// function to allocate memory, generate parameters, etc.
void initialize(yeeData *d);  

// function to deallocate memory
void finalize(yeeData *d);  

// function to setup crbc parameters
void setup_crbc(yeeData *d, CrbcUpdater2d **upd); 

// function to compute the E-field updates
void computeE(yeeData *d, int *tstep);

// function to compute the H-field updates
void computeH(yeeData *d);

// function to compute the boundary updates (E-field only)
void computeBoundary(yeeData *d, CrbcUpdater2d *upd);
     
// function to write Ez field output
void writeEzField(yeeData *d, int id); 

       
/*-----------------------------------------------------------------------------*
 *                        
 *                              Main Routine
 * 
 * ----------------------------------------------------------------------------*/
int main(int argc, char *argv[]) {   

  int tstep, i, j;
  double norm, enorm, hnorm;
  int write_output_files = 0;

  // read in input to see if we should write output files. If output files are
  // enabled this program writes out the Ez field 450 times. Each file is
  // file is approximately 2.6 MB, so this generates about 1.2 GB of data.
  // By default, the file output is turned off.
  // There is only one option, so for simplicity we'll just assume that if we
  // recieve any command line option, then we should enable output instead of
  // actually parsing and identifying a specific option.
  if (argc > 1) {
    printf("The Ez field will be saved. \n");
    write_output_files = 1;
  }
  
  // declare a yee data structure
  yeeData data;

  // declare a boundary updater object.
  // This needs to be a pointer because the size is unknown at this time and will
  // be initialized later.
  CrbcUpdater2d *boundary_updater;

  /*----------------------------------------------------------------------------

                       Set the simulation parameter
 
  ----------------------------------------------------------------------------*/
  
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

  // boundary conditions
  // The type CRBC2d_Boundaries is an provided by the CRBC/DAB library and it is 
  // simply an enumeration of the supported boundary types. At this time, the 
  // supported boundary types are 
  //     CRBC2d_DIR  --- Homogeneous Dirichlet
  //     CRBC2d_NEUM --- Homogeneous Neumman
  //     CRBC2d_CRBC --- Complete Radiation Boundary Condition (implemented as a DAB)
  //
  // The library also provides enumerations that list the valid sides, they are
  //    CRBC2d_XLeft
  //    CRBC2d_XRight
  //    CRBC2d_YLeft
  //    CRBC2d_YRight
  //
  // Here, we'll set the boundary conditions so that we a waveguide with parallel
  // PEC plates with normals in the y-direction
  data.boundaries[CRBC2d_XLeft]  = CRBC2d_CRBC;
  data.boundaries[CRBC2d_XRight] = CRBC2d_CRBC;
  data.boundaries[CRBC2d_YLeft]  = CRBC2d_DIR;
  data.boundaries[CRBC2d_YRight] = CRBC2d_DIR;

  // source parameters 
  // we will use a so called impulsive source that takes the form of a 
  // differentiated gaussian centered in the domain    
  data.tw = 5E-11;           // pulse width   
  data.t0 = 4.0 * data.tw;   // delay (this needs to be >= 4*tw to be smooth)
  data.amp = 1000;           // Amplitude 

  // specify how often to generate output (in time steps) 
  data.save_freq = 10;   

  // allocate memory and compute the remaining parameters
  initialize(&data);   

  // Now, we need to setup the boundary updater. The updater has 3 basic 
  // parameters: 
  //    delta --- The minimum seperation of each boundary face from sources, 
  //              scatterers, and other inhomogenieties.
  //    T     --- The total time being simulated.
  //    P     --- The number of recursions to use to approximate the boundary.
  //              (This can be chosen by providing a tolerance)
  //
  // Additionally, it needs to know the wave speed, c, the boundary conditions
  // on the rest of the domain. Finally, it also needs to know some properties
  // of the discretization such as the time step size, mesh sizes, and the 
  // indexing for each CRBC face.
  //
  // Each of these is discussed in more detail in the setup_crbc(...) function. 
  //
  // We need to pass the reference to the boundary_updater pointer because it 
  // is currently unitialized.
  setup_crbc(&data, &boundary_updater); 
  
  /*----------------------------------------------------------------------------

                            Run the simulation
 
  ----------------------------------------------------------------------------*/

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
          enorm += data.Ez[i][j] * data.Ez[i][j];

      // loop over Hx
      for (i=0; i<data.imax; ++i)
        for (j=0; j<data.jmax-1; ++j)
          hnorm += data.Hx[i][j] * data.Hx[i][j];

      // loop over Hy
      for (i=0; i<data.imax-1; ++i)
        for (j=0; j<data.jmax; ++j)
          hnorm += data.Hy[i][j] * data.Hy[i][j];

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

  // free the memory associated with the solver
  finalize(&data);

  // free the memory associated with the boundary updater
  CRBC2d_delete_updater(boundary_updater);

  return 0;

} // end main


/*------------------------------------------------------------------------------

                           Initialization function

------------------------------------------------------------------------------*/
void initialize(yeeData *d) {

  int i, j, k;
  int imax = d->imax;
  int jmax = d->jmax;

  d->flag = 0;

  d->pi = 3.14159265358979;   
  d->C = 2.99792458e8;   

  // time step size, (we'll use 0.99*CFL)  
  d->dt = 0.99 / (d->C*sqrt(1.0/(d->dx*d->dx) + 1.0/(d->dy*d->dy))); 

  // calculate free space eps and mu
  d->mu0 = 4.0 * d->pi * 1.0E-7;   
  d->eps0 = 1.0 / (d->C * d->C * d->mu0);

  // calculate material epsilon and mu
  d->epsilon = d->epsR * d->eps0;
  d->mu = d->mu0;

  // calculate update coefficients
  d->CE = d->dt / d->epsilon;
  d->CH = d->dt / d->mu0;

  // allocate memory for the fields
  d->Ez = (double **) malloc((imax) * sizeof(double));
  for (i=0; i < d->imax; i++) {
    d->Ez[i] = (double *) malloc((jmax) * sizeof(double));
    for (j=0; j < d->jmax; j++) {
      d->Ez[i][j] = 0.0;
    }
  }

  d->Hx = (double **) malloc((imax) * sizeof(double));
  for (i=0; i < d->imax; i++) {
    d->Hx[i] = (double *) malloc((jmax-1) * sizeof(double));
    for (j=0; j < d->jmax-1; j++) {
      d->Hx[i][j] = 0.0;
    }
  }

  d->Hy = (double **) malloc((imax-1) * sizeof(double));
  for (i=0; i < d->imax-1; i++) {
    d->Hy[i] = (double *) malloc((jmax) * sizeof(double));
    for (j=0; j < d->jmax; j++) {
      d->Hy[i][j] = 0.0;
    }
  }

  d->flag = 1;
} 

/*------------------------------------------------------------------------------

                           Finalize function

------------------------------------------------------------------------------*/
void finalize(yeeData *d) {

  int i, j;

  if (d->flag != 0) {

    // delete memory for the fields
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

/*------------------------------------------------------------------------------

                      Function to initilize the CRBC updater

------------------------------------------------------------------------------*/
void setup_crbc(yeeData *d, CrbcUpdater2d **upd) {

  int i, l;
  double T;
  double h[2];
  double delta;
  int low_index[2];
  int high_index[2];
  int n;
  CRBC2d_Side_t side;
  

  // First we will create a new boundary updater. To do this, we need to provide
  // the total simulation time, the grid spacings, the time step size, the 
  // wave speed and the boundary conditions.
  h[0] = d->dx;
  h[1] = d->dy;
  T = d->dt * d->ntsteps; // total time is time step size * number of time steps
  *upd = CRBC2d_new_updater_tol(T, h, d->dt, d->C, d->boundaries, 8, 1e-3);

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

  // Before we set up the faces, print out the properties from the boundary 
  // updater to make sure they are correct
  printf("The boundary updater was initialized with the following: \n");
  printf("  wave speed, c = %e \n", CRBC2d_get_c(*upd));
  printf("  total time, T = %e \n", CRBC2d_get_T(*upd));
  printf("  time step, dt = %e \n", CRBC2d_get_dt(*upd));
  printf("The boundary updater calculated the following: \n");
  printf("  %i edges \n", CRBC2d_get_num_sides(*upd));
  printf("  %i corners \n", CRBC2d_get_num_corners(*upd));

  // Now set up the faces. Start by looping over all of the possible faces. We
  // follow the order given in CRBC2d_Side, so 
  //   l = 0 --> CRBC2d_XLeft
  //   l = 1 --> CRBC2d_XRight
  //   l = 2 --> CRBC2d_YLeft
  //   l = 3 --> CRBC2d_YRight
  for (l=0; l<4; ++l) {

    // First we need to calculate the minimum distance between the boundary and
    // source. Since we placed the source in the center of the domain, this is 
    // simple. If it is not possible to calculate this directly, using a lower
    // bound for the seperation is the safest thing to do, but it may result
    // in more work being done that is actually needed to achieve the desired 
    // accuracy.
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

    if (d->boundaries[side] == CRBC2d_CRBC) {

    // The CRBC updater library attempts to communicate with the solvers "native"
    // indexing scheme, so we need to tell it the upper and lower indexing extents.
    // These extents need to include all of the requested variables
    // on the boundary face as well the adjacent parallel plane of points in the 
    // interior of the domain.

      // Set the basic extents and then adjust them based on the side and field 
      // component. These should be inclusive extents;
      low_index[0] = 0;
      low_index[1] = 0;

      high_index[0] = d->imax-1; // minus 1 from 0-based indexing
      high_index[1] = d->jmax-1;

      switch (side) {

        case CRBC2d_XLeft:
          // For the left boundary in x, this means we need to include the data 
          // points of the requested component at i=0 and i=1 and all possible 
          // indices for j. Therefore, the lower indexing extents are 0 for
          // all indices. The boundary updater considers these extents to be 
          // inclusive.
          high_index[0] = 1;  
    
          break;

        case CRBC2d_XRight:
          // For the right boundary in x, this means we need to include the data
          // points of the requested component at i=imax-1 and i=imax-2 along with all of 
          // the possible indices for j. The boundary updater considers 
          // these extents to be inclusive.
          low_index[0] = d->imax-2;

          break;

        case CRBC2d_YLeft:
          // For the left boundary in y, this means we need to include the data 
          // points of the requested component at j=0 and j=1 along with all of 
          // the possible indices for i. The boundary updater considers 
          // these extents to be inclusive.
          high_index[1] = 1;
    
          break;  

        case CRBC2d_YRight:
          // For the right boundary in y, this means we need to include the data
          // points of the requested component at j=jmax-1 and j=jmax-2 along with all of 
          // the possible indices for i. The boundary updater considers 
          // extents to be inclusive.
          low_index[1] = d->jmax-2;

          break;

      } // end switch side

      // now initialize the face.
      if (CRBC2d_init_face(*upd, side, low_index, high_index, delta) != 0)
      {
        fprintf(stderr, "Error: init_face(...) failed \n");
        exit(-1);
      }
    } // end if
  } // end for over possible boundary faces (l)

  // now we'll print out some information about the recursions
  // NOTE  a reflection coefficient of -1 indicates that we are not performing
  // any updates on the face.
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

/*------------------------------------------------------------------------------

                          Compute E-field updates

------------------------------------------------------------------------------*/
void computeE(yeeData *d, int *tstep) {

  int i, j;

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

/*------------------------------------------------------------------------------

                          Compute H-field updates

------------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------------

                   Function to calculate the CRBC updates

------------------------------------------------------------------------------*/
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


  // now we can compute the updates to the boundaries
  if (CRBC2d_compute_updates(upd) != 0) {
    fprintf(stderr, "Error: compute_updates(...) failed \n");
    exit(-1);
  }

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

/*------------------------------------------------------------------------------
 
                                 Write output

------------------------------------------------------------------------------*/
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


  
