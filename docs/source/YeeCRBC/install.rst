************
Installation
************

CMake
=====

The simplest method to install the library and build the examples is to use 
CMake. For the default install, we suggest the following. First, extract the 
source ::

  tar -zxvf yee_crbc_lib.tar.gz

We suggest building the library in a separate directory. For this example we 
will call the directory 'build'. Then, we simply have to invoke CMake. ::

  mkdir build
  cd build
  cmake ../yee_crbc_lib
  make
  make install

This will install the library and build the examples.

The library uses some functionality of C++11. For most modern compilers, we have
found that it is not necessary to include any special compiler flags to use these
features (e.g. Intel 15+, GNU 4.8+). However, some compilers do require
we specifically tell it to use C++11. At this time it appears that CMake only
identifies these flags correctly for the GNU compilers, so it may be necessary to
tell CMake to include the appropriate flag. For example: ::

  cmake ../yee_crbc_lib -DCMAKE_CXX_FLAGS:STRING=-std=c++11 

The correct flag varies by compiler. Also note that if your compiler uses the GNU
STL implementation, it is necessary that the GNU library version be >=4.8.

Options
-------

There are additional options that can be enabled. 

Shared Library
^^^^^^^^^^^^^^

To build the shared version of the library use the command ::

    -DBUILD_SHARED="ON"

Install Path
^^^^^^^^^^^^

The default install places the library file in ``/usr/local/lib``
and places the header files in ``/usr/local/include/yeecrbc``, although this behavior
can vary from system to system. The install can be changed to place the 
library in ``/usr/lib`` and the headers in ``/usr/include/yeecrbc``, for example,
by adding the following CMake option. ::

  -DCMAKE_INSTALL_PREFIX:PATH=/usr 

Index and Data Types
^^^^^^^^^^^^^^^^^^^^

The defaults are to use the ``int`` type for indexing and the ``double`` type for
the field data type and coefficient data type. In the C interfaces, the indexing
type can be changed to the ``long int`` or ``long long`` data type by adding the 
following command line options to the CMake invocation ::

   -DINDEX_TYPE:STRING="long int"

or ::

   -DINDEX_TYPE:STRING="long long"

Similarly, the data type we use for a majority of the computations can be changed
from the default type of ``double`` to ``float`` by adding the following ::

   -DUSE_CRBC_SCALAR_FLOAT="ON"

and the coefficient data type can be changed from ``double`` to ``float`` with ::

   -DUSE_CRBC_COEF_FLOAT="ON"

Note that the underlying C++ code and interfaces is templated on the index and
data types, so these commands only affect the C interfaces.


Save and Restart
^^^^^^^^^^^^^^^^

The library includes the optional ability
to save the state of the boundary updater and restart from a previously
saved state. This functionality requires the use of the HDF5 library. To enable
this function from the command line add the following to the CMake invocation ::

  -DUSE_HDF5_RESTARTS=ON

Threading
^^^^^^^^^

Another option is to enable OpenMP threading. To do this, add ::

  -DUSE_OPENMP=ON

Debugging
^^^^^^^^^

To build in Debug mode, add ::

  -DCMAKE_BUILD_TYPE=Debug

Some additionaly error checking and messages can be enabled with the option ::

  -DENABLE_EXTRA_DEBUG="ON"

Note that builds with these options enabled significantly degrades performance.

Linking
-------

CMake Linking
^^^^^^^^^^^^^

In order to link with Cmake, a ``CMakeLists.txt`` file should include ::

  find_package(yeecrbc REQUIRED)
  target_link_libraries(foo yeecrbc)

Depending on the install location Cmake may or may not be able to locate the 
libraries configuration file. If it is not found, the variable ``yeecrbc_DIR`` 
needs to be set to the install location. Using the default install settings this
can be done with the following CMake command ::

  -Dyeecrbc_DIR=/usr/local/lib

Alternatively, this can be set in the ``CMakeLists.txt`` file by adding ::

  set(yeecrbc_DIR "/usr/local/lib")

Manual Linking
^^^^^^^^^^^^^^

To link to the library, add ``-lyeecrbc`` to the compile string. If the library
was not installed in a standard path, it is necessary to provide the path to
the library and the header files. For example, using the typical default install
settings ::

  -L/usr/local/lib -lyeecrbc -I/usr/local/include/yeecrbc

Since the majority of the library is written in C++, ``-lstdc++`` is often needed
to compile.

If the HDF5 restart capability is enabled, the HDF library also needs
to be included in the linking. Typically, this is done with ``-lhdf5``. If the
threading is enabled, the appropriate OpenMP compiler flag needs to be included.
For the GNU compilers, the flag is ``-fopenmp``; for Intel compilers, the flag
is ``-openmp``; and for Portland Group compilers, the flag is ``-mp``.

Running Dynamically Linked Executables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the shared version of the library, it may be necessary to provide the
path to the library when executing a program. There are a number of way to do
this, but the easiest way to do this is to specify the LD_LIBRARY_PATH. For
example (in bash) ::

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

Uninstalling
============

To uninstall the library, simply delete the file ``/usr/local/lib/libyeecrbc.a``
and the folder ``/usr/local/include/yeecrbc``. Note that the location of the include
folder and the library file may vary from system to system and depends on the
install configurations. In general the paths are ``<prefix>/lib/libyeecrbc.a``
and ``<prefix>/include/crbc``.
