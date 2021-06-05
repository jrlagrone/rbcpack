Examples
========

Simple Examples
---------------

To help demonstrate the usage of the CRBC/DAB library, we have constructed some
simple examples.

Transverse-magnetic Yee example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first example is a transverse-magnetic wave guide with a simple scatter. This
examples demonstrated the use of the 2D interface.

.. toctree::
   :maxdepth: 2

   yee_TM

3D Yee example
^^^^^^^^^^^^^^

Next we have an example using the full 3D Yee scheme. This example is set up to 
be a simulation between two parallel PEC plates with a Gaussian source in one of
the components. We note that this example is flexible in the sense that the 
boundaries can be changed to any configuration of PEC and radiating boundaries.
This example demonstrates the 3D Yee interface.

.. toctree::
   :maxdepth: 2

   3D_Yee_example_in_C


Advanced Examples
-----------------

We also present some more advanced examples that take advantage of some of the
more complex features of the underlying C++ code.

First we demonstrate running a scalar wave equation in 1D, 2D, and 3D taking 
advantage of the fact that the C++ interface is templated on dimension (note
that the underlying functions are all dimension independent, but the interface
only supports up to 3 dimensions).

.. toctree::
   :maxdepth: 2

   wave_eq_example_in_Cpp

Next we demonstrate a possible implementation of the 3D Yee scheme using the C++
interface and MPI (note that we are not offcially supporting MPI at this time)

.. toctree::
   :maxdepth: 2

   yee_mpi

