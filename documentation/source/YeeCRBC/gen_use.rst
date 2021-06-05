*************
General Usage
*************

Typical Work Flow
=================

Procedure for Calling the Library
---------------------------------

To use the library, the basic workflow is ::

                       -----------------------------
                       | initialize updater object |
                       -----------------------------
                                    |
                                    v
                       -----------------------------
                       | initialize boundary faces |
                       -----------------------------
                                    |
                                    v
                               -------------
           ------------------> | time step |--------------------
           ^                   -------------                   |
           |                                                   |
           |                                                   v 
  ---------------------      -------------------      ---------------------
  | copy new boundary |      |  updater object |      |  copy new values  |
  |  values into the  | <----|  computes new   | <----|  into the updater |
  |      solver       |      | boundary values |      |      object       |
  ---------------------      -------------------      ---------------------  

Initializing an Updater Object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To begin, a boundary updater object needs to be initialized. Currently, there
are three interfaces available: a 3D FDTD/Yee interface in C, a 2D interface in 
C, and the underlying C++ interface. In all three cases, a boundary object is
created in a similar manner. Each interface has three initializers: one to use
the default number of recursions (which we have set to 5), one to set the number
of recursions to use to a specific number, and finally one that chooses the
number of recursions to use based on a provided tolerance.

In the 3D FDTD/Yee interface these initializers are ::

  // create a new updater object with default settings
  YeeCrbcUpdater* CRBC_new_updater (const CoefType T,
                                    const CoefType h[3],  
                                    const CoefType dt,                  
                                    const CoefType c,
                                    const CRBC_Boundaries_t boundaries[6]);

  // create a new updater object and specify the number of recursions
  YeeCrbcUpdater* CRBC_new_updater_p (const CoefType T,
                                      const CoefType h[3],  
                                      const CoefType dt,                  
                                      const CoefType c,
                                      const CRBC_Boundaries_t boundaries[6],
                                      const int P);

  // create a new updater object and specify the a tolerance for the recursions
  YeeCrbcUpdater* CRBC_new_updater_tol (const CoefType T,
                                        const CoefType h[3],  
                                        const CoefType dt,                  
                                        const CoefType c,
                                        const CRBC_Boundaries_t boundaries[6],
                                        const int Pmax,
                                        const double tol);

Where the return type is provided by the library and is simply defined as ::

  struct YeeCrbcUpdater;
  typedef struct YeeCrbcUpdater YeeCrbcUpdater;

The ``CoefType`` is defined when the library is compiled and can either be
``double`` or ``float`` for the C interfaces and is a template parameter for the
C++ interfaces. Finally the ``CRBC_Boundaries_t`` is an enumeration 
type provided by the library defined ::

  typedef enum CRBC_Side {                                                            
    CRBC_XLeft,  /**< Left side in the x-coordinate direction */
    CRBC_XRight, /**< Right side in the x-coordinate direction */
    CRBC_YLeft,  /**< Left side in the y-coordinate direction */
    CRBC_YRight, /**< Right side in the y-coordinate direction */
    CRBC_ZLeft,  /**< Left side in the z-coordinate direction */
    CRBC_ZRight  /**< Right side in the z-coordinate direction */
  } CRBC_Side_t; 

The 2D interface is almost exactly the same: ::
                            
  // create a new updater object with default settings
  CrbcUpdater2d* CRBC2d_new_updater (const CoefType T,    // total time
                                     const CoefType h[2], // grid spacings  
                                     const CoefType dt,   // time step               
                                     const CoefType c,    // wave speed
                                     const CRBC2d_Boundaries_t boundaries[4]);

  // create a new updater object and specify the number of recursions
  CrbcUpdater2d* CRBC2d_new_updater_p (const CoefType T,
                                       const CoefType h[2],  
                                       const CoefType dt,                  
                                       const CoefType c,
                                       const CRBC2d_Boundaries_t boundaries[4],
                                       const int P);

  //  create a new updater object and specify the a tolerance for the recursions
  CrbcUpdater2d* CRBC2d_new_updater_tol (const CoefType T,
                                       const CoefType h[2],  
                                       const CoefType dt,                  
                                       const CoefType c,
                                       const CRBC2d_Boundaries_t boundaries[4],
                                       const int Pmax,
                                       const double tol);

Finally, the C++ interface differs somewhat. ::

    /// Constructor --- initialize only the parameters that are common to all of
    /// the faces in all circumstances.
    CrbcUpdates (const CoefType &T,
                 const CoefType h[DIM],  
                 const CoefType &dt,                  
                 const CoefType &c,
                 const Boundary boundaries[2*DIM]);

    /// Constructor --- initialize only the parameters that are common to all of
    /// the faces and set the number of recursions to be P in all directions
    CrbcUpdates (const CoefType &T,
                 const CoefType h[DIM],  
                 const CoefType &dt,                  
                 const CoefType &c,
                 const Boundary boundaries[2*DIM],
                 const int &P);

    /// Constructor --- initialize only the parameters that are common to all of
    /// the faces and set the number of recursions to be calculated based on the
    /// provided tolerance
    CrbcUpdates (const CoefType &T,
                 const CoefType h[DIM],  
                 const CoefType &dt,                  
                 const CoefType &c,
                 const Boundary boundaries[2*DIM],
                 const int &Pmax,
                 const double &tol);

Initializing Boundary Faces
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The library requires that each of the boundary faces (or sides in 2D) are initialized
individually. We do this for a number of reasons, but the most important is the
fact that each boundary face may pose a different degree of difficulty.

On each face, the boundary updater expects to get two things. The first is the
parameter, :math:`\delta`, which represents the separation from the boundary
face to any sources, scatterers, or other inhomogeneities. The second is 
the data extents that will serve as in the inputs an outputs of the library.
In order to function correctly, these extents need to include all of the points
on the boundary as well as the parallel plane (or line in 2D) of points immediately
interior to the boundary. 

For example, suppose we have a computational domain
with :math:`n_x` grid points in the x-direction, :math:`n_y` grid points in the
y-direction, and :math:`n_z` grid points in the z-direction.
Then, if we want to set up the boundary face on the right side in the x-direction,
we will require all of the points on the boundary, which would be the points at
:math:`i=n_x`, :math:`j=1,...,n_y`, and :math:`k=1,...,n_z` assuming the indexing
begins at one. Then we also need the parallel plane of points on the interior side
of this boundary, which are the points :math:`i=n_x-1`, :math:`j=1,...,n_y`, 
and :math:`k=1,...,n_z`.  
Therefore, our indexing extents should be :math:`[n_x-1, n_x] \times [1, n_y] \times
[1, n_z]`. However, the library expects the lower extents to be in one array and
the upper extents to be in a second array, so are inputs should be
:math:`[n_x-1, 1, 1]` and :math:`[n_x, n_y , n_z]` for the lower and upper 
extents, respectively. Note that the library considers these extents to be inclusive. 

For the 3D Yee interface this same procedure needs to be applied all of the components
of the **E** field required to generate the boundary conditions. Note that the 
library can provide a list of the components it expects by calling the function ::

  // get the components we expect to be inputted into the updater
  void CRBC_get_component_inputs (YeeCrbcUpdater *U, 
                           const CRBC_Side_t side, 
                           CRBC_Fields_t *comp, 
                           int *n);

The only odd thing here is the we sometimes require the normal **E** field component.
We note that this component is half of a grid space inside the boundary. In this 
case, we require the plane of points closest to the boundary and the adjacent 
parallel plane of points in the interior. These points are only needed as input
to update edges and corners.

In all three available interfaces, the function call to initialize the faces are
similar. The only practical difference is that in the 3D Yee interface, we also
need to know which component is being initialized. 

In the 3D Yee interface, the call is ::

  // initialize a component on a face
  int CRBC_init_face (YeeCrbcUpdater *U,
                      const CRBC_Side_t side,
                      const CRBC_Fields_t comp, 
                      const IndexType low_index[3],
                      const IndexType high_index[3],
                      const double delta);

In the 2D interface, the call is ::

  // initialize a boundary side
  int CRBC2d_init_face (CrbcUpdater2d *U,
                        const CRBC2d_Side_t side,
                        const IndexType low_index[2],
                        const IndexType high_index[2],
                        const double delta);

Finally, in the C++ interface, the call is ::

  void init_face(const int &side,
                 const IndexType low_index[DIM],
                 const IndexType high_index[DIM],
                 const double &delta);

Copying and Updating Values
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The remaining process of using the library is simply passing information back
and forth. After the solver used to update the interior values has taken a time
step, we need to load the set of points adjacent to the boundary that were updated
by the solver. To do this, we ask the boundary updater object to provide us with
the indexing extents it expects to receive as inputs. Depending on the interface,
this is one of the following calls ::

  void CRBC_get_input_extents (YeeCrbcUpdater *U, 
                       const CRBC_Side_t side, 
                       const CRBC_Fields_t comp, 
                       IndexType low[3], 
                       IndexType high[3]);

  void CRBC2d_get_input_extents (CrbcUpdater2d *U, 
                       const CRBC2d_Side_t side, 
                       IndexType low[2], 
                       IndexType high[2]);

  void get_input_extents(const int &side, 
                         IndexType low[DIM], 
                         IndexType high[DIM]) const;

Then, we simply copy the values from the solver into the boundary updater object
at these extents using one of the following calls, depending on the interface, ::

  void CRBC_load_face_data (YeeCrbcUpdater *U,
                            const CRBC_Side_t side,
                            const CRBC_Fields_t comp,
                            const IndexType P[3],
                            const ScalarType *val);

  void CRBC2d_load_face_data (CrbcUpdater2d *U,
                              const CRBC2d_Side_t side,
                              const IndexType P[2],
                              const ScalarType *val);

  void load_face_data(const int &side,
                      const IndexType P[DIM],
                      const DataType &val);

After doing this on all of the sides with a CRBC/DAB boundary condition, we can
ask the boundary updater object to compute the boundary objects with one of ::

  int CRBC_compute_updates (YeeCrbcUpdater *U);

  int CRBC2d_compute_updates (CrbcUpdater2d *U);

  void compute_updates();

Finally, we need to copy the updated boundary values into the solver. We do this
by asking the boundary updater object for extents it thinks it should output
data for and then copying the values. We note, of the 3D Yee interface, that it
is not necessary to copy the normal **E** field components into the solver because
they should already be correct (we needed the values to handle edges and corners).
This process involves first calling ::

  void CRBC_get_output_extents (YeeCrbcUpdater *U,
                                const CRBC_Side_t side, 
                                const CRBC_Fields_t comp, 
                                IndexType low[3], 
                                IndexType high[3]);

  void CRBC2d_get_output_extents (CrbcUpdater2d *U,
                                  const CRBC2d_Side_t side, 
                                  IndexType low[2], 
                                  IndexType high[2]);

  void get_output_extents(const int &side, 
                          IndexType low[DIM], 
                          IndexType high[DIM]) const;

and then copying the values using ::

  ScalarType CRBC_get_new_face_vals (YeeCrbcUpdater *U,
                                     const CRBC_Side_t side,
                                     const CRBC_Fields_t comp,
                                     const IndexType P[3]);

  ScalarType CRBC2d_get_new_face_vals (CrbcUpdater2d *U,
                                       const CRBC2d_Side_t side,
                                       const IndexType P[2]);

  DataType get_new_face_vals(const int &side,
                             const IndexType P[DIM]) const;

This process is repeated every time step.
  

More Details
^^^^^^^^^^^^

For more details, please refer to the Doxygen generated documentation <link>;
the documentation contained in the interface header files ``3d_yee_crbc_api.h``,
``2d_crbc_api.h``, and ``crbc_updates.hpp``; and the :doc:`examples`. 

Choice of Parameters
====================

The boundary conditions in this library have three basic parameters: the wave 
speed, :math:`c`, a separation parameter, :math:`\delta`, and a time parameter,
:math:`T`. We use these three parameters to form a single dimensionless parameter

.. math::

  \eta = \frac{\delta}{c T},

which we use as a measure of the difficulty of the simulation. The basic idea is
that larger values of :math:`\eta` are easier problems and smaller values are
harder problems. For instance, :math:`\eta \geq 1` implies that a boundary condition
is unnecessary because a wave cannot impact it in time :math:`T`. In the library,
we use values of :math:`1e-8 \leq \eta \leq 0.1`. If the calculated value falls
outside of this range, we simply choose :math:`\eta` to represent either an 
"easier" or "harder" problem as appropriate.

The parameter :math:`\delta` should be chosen as the minimum distance
from the boundary to any sources, scatterers, or other inhomogeneities in the 
interior of the domain. In general, if this parameter cannot be calculated directly,
a safe option is to compute a lower bound on this distance that is greater than
zero.

In general, the time parameter :math:`T` should be chosen to be the total amount
of time that the simulation is being run. However, in cases where all of the 
energy is expected to leave the domain quickly, particularly free space simulations,
it is sometimes reasonable to choose a value of :math:`T` that is on the order
of the time it will take for all of the energy to leave the domain.

Tolerance
---------

In general, we prefer the tolerance based system provided we have a way to choose
a reasonable tolerance. The tolerance is used to control the reflection coefficient
of the boundary. We have found that this typically provides a reasonable bound
on the relative error and is, therefore, generally better 
than simply guessing an appropriate number of recursions to use.

As a simple rule, it is best to choose a tolerance that is on the order of the
discretization/dispersion error for the length of the simulation. The reason for
this is that the accuracy of the boundary condition has little impact if it is 
more accurate than the discretization error in the interior. Conversely, if the 
boundary condition is less accurate than the interior, the error from the boundary
condition will dominate. That being said, it is often difficult to determine such
a tolerance before running simulations.

Number of Recursions
--------------------

With all else remaining equal, letting the number of recursions go to infinity 
will result in an analytically exact boundary condition. In computations, this 
is, of course, impossible. We have set the default number of recursions to be
5. We have chosen this primarily because this has seemed to be a good choice
in a majority of our experiments.  Additionally, it is roughly comparable in
computational work to a PML of 10 cells. For a typical simulation, using 5 recursions usually
gives a reflection coefficient that is :math:`\mathcal{O}(10^{-3})` or better (with variations
depending on the particular problem). Choosing fewer recursions will be less accurate
but the computational expense will be less. Choosing more recursions requires
more work, but results in greater accuracy. 

For reference, the work scales
linearly in the number of recursions on the faces, quadratically on the edges, 
and cubically on the corners. 
The net result is that work scales linearly on a sufficiently large computational
domain (e.g. there are many more points on a face than the number of recursions).

