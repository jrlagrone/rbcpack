Install and Running the Included Examples
=========================================

To built the library we require CMake (version >= 2.8.12). 

To build the library, extract the library and set up a build directory. Then 
invoke cmake, make, and optionally make install. ::

  tar -zxvf yee_crbc_lib.tar.gz
  mkdir build
  cd build
  cmake ../yee_crbc_lib
  make
  make install

For some compilers, it may be necessary to add the flag "-DCMAKE_CXX_FLAGS:STRING=-std=c++11"
to the cmake command. For more detailed instructions and options for installation,
please refer to :doc:`../install`.

To test the build, the 2d Yee Scheme example (4500 iterations), 3d Yee (1200 iterations),
or wave equation examples can be run by navigating to the `example` directory in
the `build` folder, e.g.::

  cd build/examples/2d_yee
  ./yee_TM.x

For more details about these programs, please see :doc:`../examples`.

To build an example again (e.g. 2d_yee), it is only necessary to type make in 
the example directory, e.g. ::

  cd build/examples/2d_yee
  make

To manually compile and link the example (e.g. 2d_yee), note that one need to specify
the specify the location the library was installed (the default is `/usr/local/` ::

  gcc -L/usr/local/lib -I/usr/local/include/crbc yee_TM.c -o yee_TM.x -lyeecrbc -lstdc++ -lm

For more detailed information please refer to the following user guides:

.. toctree::
   :maxdepth: 2

   ../install
   ../supported_problems
   ../gen_use
   ../examples

