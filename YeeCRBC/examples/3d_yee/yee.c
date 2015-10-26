/**
    3-D FDTD code with CRBC/DAB boundary conditions

    This C code implements the finite difference time-domain solution of
    Maxwell's equations (curl formulation) over a 3-D Cartesian lattice.
    The grid is terminated using Double Absorbing Boundaries (DAB).

    This is meant to serve as an example of how one might use the CRBC/DAB 
    library. As such the Yee algorithm used is as simple as is practicle to 
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
// the interface to the CRBC/DAB library specifilized for the 3D Yee scheme
#include <3d_yee_crbc_api.h>

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
  // Ex - at (i + 1/2, j, k)
  // Ey - at (i, j + 1/2, k)
  // Ez - at (i, j, k + 1/2)
  //  
  // Hx - at (i, j + 1/2, k + 1/2)
  // Hy - at (i + 1/2, j, k + 1/2)
  // Hz - at (i + 1/2, j + 1/2, k) 
  //
  // This means we have
  // (imax-1) * (jmax) * (kmax)   Ex field values 
  // (imax) * (jmax-1) * (kmax)   Ey field values 
  // (imax) * (jmax) * (kmax-1)   Ez field values 
  // (imax) * (jmax-1) * (kmax-1) Hx field values 
  // (imax-1) * (jmax) * (kmax-1) Hy field values 
  // (imax-1) * (jmax-1) * (kmax) Hz field values 
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
  // The type CRBC_Boundaries is an provided by the CRBC/DAB library and it is 
  // simply an enumeration of the supported boundary types. At this time, the 
  // supported boundary types are 
  //     CRBC_PEC  --- Perfect Electric Conductor
  //     CRBC_CRBC --- Complete Radiation Boundary Condition (implemented as a DAB)
  CRBC_Boundaries_t boundaries[6];

  // source parameters 
  // we will use a so called impulsive source that takes the form of a 
  // differentiated gaussian centered in the domain    
  double tw;      // pulse width   
  double t0;      // delay   
  double amp;     // Amplitude 

  // specify how often to generate output (in time steps) 
  int save_freq;   

  // H & E Field components  
  // NOTE: storing these as a 1D array and accessing appropriately is likely to 
  //       be more efficient, but we do this for better clarity. 
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
       

// function to allocate memory, generate parameters, etc.
void initialize(yeeData *d);  

// function to deallocate memory
void finalize(yeeData *d);  

// function to setup crbc parameters
void setup_crbc(yeeData *d, YeeCrbcUpdater **upd); 

// function to compute the E-field updates
void computeE(yeeData *d, int *tstep);

// function to compute the H-field updates
void computeH(yeeData *d);

// function to compute the boundary updates (E-field only)
void computeBoundary(yeeData *d, YeeCrbcUpdater *upd);
     
// function to write Ex field output
void writeExField(yeeData *d, int id); 

       
/*-----------------------------------------------------------------------------*
 *                        
 *                              Main Routine
 * 
 * ----------------------------------------------------------------------------*/
int main(int argc, char *argv[]) {   

  int tstep, i, j, k;
  double norm, enorm, hnorm;
  int write_output_files = 0;

  // read in input to see if we should write output files. If output files are
  // enabled this program writes out the Ex field 50 times. this generates about 
  // 350 MB of data.
  // By default, the file output is turned off.
  // There is only one option, so for simplicity we'll just assume that if we
  // recieve any command line option, then we should enable output instead of
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

  /*----------------------------------------------------------------------------

                       Set the simulation parameter
 
  ----------------------------------------------------------------------------*/
  
  // relative permeability
  data.epsR = 1.0; // 1.0 corresponds to a vacuum  

  // number of time steps 
  data.ntsteps = 1200;   

  // grid size
  // NOTE: Field Value Locations
  //
  // Ex - at (i + 1/2, j, k)
  // Ey - at (i, j + 1/2, k)
  // Ez - at (i, j, k + 1/2)
  //  
  // Hx - at (i, j + 1/2, k + 1/2)
  // Hy - at (i + 1/2, j, k + 1/2)
  // Hz - at (i + 1/2, j + 1/2, k) 
  //
  // This means, e.g. we have
  // (imax-1) * (jmax) * (kmax) Ex field values    
  data.imax = 200;   
  data.jmax = 200;   
  data.kmax = 200;  

  // grid spacing in each direction
  data.dx = 1e-3;   
  data.dy = 1e-3;   
  data.dz = 1e-3;  

  // boundary conditions
  // The type CRBC_Boundaries is an provided by the CRBC/DAB library and it is 
  // simply an enumeration of the supported boundary types. At this time, the 
  // supported boundary types are 
  //     CRBC_PEC  --- Perfect Electric Conductor
  //     CRBC_CRBC --- Complete Radiation Boundary Condition (implemented as a DAB)
  //
  // The library also provides enumerations that list the valid sides, they are
  //    CRBC_XLeft
  //    CRBC_XRight
  //    CRBC_YLeft
  //    CRBC_YRight
  //    CRBC_ZLeft
  //    CRBC_ZRight 
  //
  // Here, we'll set the boundary conditions so that we have parallel PEC plates
  // with normals in the z-direction
  data.boundaries[CRBC_XLeft]  = CRBC_CRBC;
  data.boundaries[CRBC_XRight] = CRBC_CRBC;
  data.boundaries[CRBC_YLeft]  = CRBC_CRBC;
  data.boundaries[CRBC_YRight] = CRBC_CRBC;
  data.boundaries[CRBC_ZLeft]  = CRBC_CRBC;
  data.boundaries[CRBC_ZRight] = CRBC_CRBC;

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

    // test the restart capability (requires HDF5)
    #if USE_HDF5_RESTARTS
    if (tstep == data.ntsteps /2) {

      // This will write out the current state of the boundary library to 
      // several h5 files. In this case, it will save a basic parameter file
      // called hdf5_save_test.h5 and files for each of the three E-field
      // components in the files hdf5_save_test_ex.h5, hdf5_save_test_ey.h5, and
      // hdf5_save_test_ez.h5. Depending on the boundary configuration, it may
      // not save all of the boundary components.
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

  // free the memory associated with the solver
  finalize(&data);

  // free the memory associated with the boundary updater
  CRBC_delete_updater(boundary_updater);

  return 0;

} // end main


/*------------------------------------------------------------------------------

                           Initialization function

------------------------------------------------------------------------------*/
void initialize(yeeData *d) {

  int i, j, k;
  int imax = d->imax;
  int jmax = d->jmax;
  int kmax = d->kmax;

  d->flag = 0;

  d->pi = 3.14159265358979;   
  d->C = 2.99792458e8;   

  // time step size, (we'll use 0.99*CFL)  
  d->dt = 0.99 / (d->C*sqrt(1.0/(d->dx*d->dx) + 1.0/(d->dy*d->dy) + 1.0/(d->dz*d->dz))); 

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

/*------------------------------------------------------------------------------

                           Finalize function

------------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------------

                      Function to initilize the CRBC updater

------------------------------------------------------------------------------*/
void setup_crbc(yeeData *d, YeeCrbcUpdater **upd) {

  int i, l;
  double T;
  double h[3];
  double delta;
  int low_index[3];
  int high_index[3];
  int n;
  CRBC_Fields_t comps[3];
  CRBC_Side_t side;
  

  // First we will create a new boundary updater. To do this, we need to provide
  // the total simulation time, the grid spacings, the time step size, the 
  // wave speed and the boundary conditions.
  h[0] = d->dx;
  h[1] = d->dy;
  h[2] = d->dz;
  T = d->dt * d->ntsteps; // total time is time step size * number of time steps
  *upd = CRBC_new_updater(T, h, d->dt, d->C, d->boundaries);

  // alternatively one can call 
  // int P = 7;
  // upd = new_updater_p(T, h, d->dt, d->C, d->boundaries, P);
  // This does exactly the same thing but changes the number of recursions from
  // the default of 5 to 7.
  //
  // or 
  // int Pmax = 15;
  // double tol = 1e-3;
  // upd = new_updater_tol(T, h, d->dt, d->C, d->boundaries, Pmax, tol);
  // Will generate an updater that attempts to satsify the provided tolerance
  // and with fewer than Pmax recursions.

  // Before we set up the faces, print out the properties from the boundary 
  // updater to make sure they are correct
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

  // Now set up the faces. Start by looping over all of the possible faces. We
  // follow the order given in CRBC_Side, so 
  //   l = 0 --> CRBC_XLeft
  //   l = 1 --> CRBC_XRight
  //   l = 2 --> CRBC_YLeft
  //   l = 3 --> CRBC_YRight
  //   l = 4 --> CRBC_ZLeft
  //   l = 5 --> CRBC_ZRight
  for (l=0; l<6; ++l) {

    // First we need to calculate the minimum distance between the boundary and
    // source. Since we placed the source in the center of the domain, this is 
    // simple. If it is not possible to calculate this directly, using a lower
    // bound for the seperation is the safest thing to do, but it may result
    // in more work being done that is actually needed to achieve the desired 
    // accuracy.
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

    // Next, we'll get the field components that the updater expects. If this
    // is not a boundary that will be handled by the updater, this should 
    // return n = 0.
    CRBC_get_component_inputs(*upd, side, comps, &n);

    // Set up each component the updater has requested.
    for (i=0; i<n; ++i) {

      // The CRBC updater library attempts to communicate with the solvers "native"
      // indexing scheme, so we need to tell it the upper and lower indexing extents.
      // These extents need to include all of the requested component variables
      // on the boundary face (or just inside the boundary in the case of a normal 
      // component) as well the adjacent parallel plane of points in the interior 
      // of the domain.

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
          // For the left boundary in x, this means we need to include the data 
          // points of the requested component at i=0 and i=1 and all possible 
          // indices for j and k. Therefore, the lower indexing extents are 0 for
          // all indices. The upper extents vary by component. The boundary 
          // updater considers these extents to be inclusive.
          high_index[0] = 1;  
    
          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ey) {
            high_index[1]--;
          } else if (comps[i] == CRBC_Ez) {
            high_index[2]--;
          }

          break;

        case CRBC_XRight:
          // For the right boundary in x, this means we need to include the data
          // points of the requested component at i=imax-1 and i=imax-2 or 
          // i=imax-2 and i=imax-3 depending on the component along with all of 
          // the possible indices for j and k. The boundary updater considers 
          // these extents to be inclusive.
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
          // For the left boundary in y, this means we need to include the data 
          // points of the requested component at j=0 and j=1 along with all of 
          // the possible indices for i and k. The boundary updater considers 
          // these extents to be inclusive.
          high_index[1] = 1;
    
          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ex) {
            high_index[0]--;
          } else if (comps[i] == CRBC_Ez) {
            high_index[2]--;
          }
          
          break;  

        case CRBC_YRight:
          // For the right boundary in y, this means we need to include the data
          // points of the requested component at j=jmax-1 and j=jmax-2 or 
          // j=jmax-2 and j=jmax-3 depending on the component along with all of 
          // the possible indices for i and k. The boundary updater considers 
          // extents to be inclusive.
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
          // For the left boundary in z, this means we need to include the data 
          // points of the requested component at k=0 and k=1 along with all of 
          // the possible indices for j and k. The boundary updater considers 
          // these extents to be inclusive.
          high_index[2] = 1;
    
          // adjust the upper extent based on the field component
          if (comps[i] == CRBC_Ex) {
            high_index[0]--;
          } else if (comps[i] == CRBC_Ey) {
            high_index[1]--;
          }

          break;

        case CRBC_ZRight:
          // For the right boundary in z, this means we need to include the data
          // points of the requested component at k=kmax-1 and k=kmax-2 or 
          // k=kmax-2 and k=kmax-3 depending on the component along with all of 
          // the possible indices for i and j. The boundary updater considers 
          // these extents to be inclusive.
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

      // now initialize the face.
      if (CRBC_init_face(*upd, side, comps[i], low_index, high_index, delta) != 0)
      {
        fprintf(stderr, "Error: init_face(...) failed \n");
        exit(-1);
      }

    }  // end for over requested components (i)
  } // end for over possible boundary faces (l)

  // now we'll print out some information about the recursions
  // NOTE  a reflection coefficient of -1 indicates that we are not performing
  // any updates on the face.
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

/*------------------------------------------------------------------------------

                          Compute E-field updates

------------------------------------------------------------------------------*/
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
                                * ((*tstep*d->dt - d->t0) / d->tw) \
                                * exp(-pow(((*tstep*d->dt - d->t0) / d->tw), 2.0));

}

/*------------------------------------------------------------------------------

                          Compute H-field updates

------------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------------

                   Function to calculate the CRBC updates

------------------------------------------------------------------------------*/
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


  // now we can compute the updates to the boundaries
  if (CRBC_compute_updates(upd) != 0) {
    fprintf(stderr, "Error: compute_updates(...) failed \n");
    exit(-1);
  }

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

/*------------------------------------------------------------------------------
 
                                 Write output

------------------------------------------------------------------------------*/
void writeExField(yeeData *d, int id) {

  int i, j, k, n, cells;

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


  
