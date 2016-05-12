# This file creates the CRBCYee library and optionally builds the examples

# We use some commands that require v. >= 3.0. This is available in the default
# ubuntu repositories or it can be compiled from source available at cmake.org
# (in our experience, it is straightforward to compile)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

# make sure a build type is set otherwise set to Release
if(NOT CMAKE_BUILD_TYPE) 
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE) 
endif(NOT CMAKE_BUILD_TYPE)

# find the math library (this is usually isn't needed, but it doesn't work under
# some system configurations without it
find_library(M_LIB m)
set(OPT_LIBS "${M_LIB}")

###
# Add the ability to compile optional features
# OPEN MP threading
option(USE_OPENMP "OpenMP threading" OFF)
option(ENABLE_EXTRA_DEBUG "Enable extra debugging options, significantly degrades performance" OFF)
option(USE_HDF5_RESTARTS "Enable HDF5 input/output" OFF)
option(USE_CRBC_SCALAR_FLOAT "Enable to use float for field data type, default is double" OFF)
option(USE_CRBC_COEF_FLOAT "Enable to use float for coefficient data type, default is double" OFF)
option(BUILD_SHARED "Build the library as a shared object" OFF)

set(INDEX_TYPE "INDEX_TYPE" CACHE STRING "Index type chosen by the user at CMake configure time")
set_property(CACHE INDEX_TYPE PROPERTY STRINGS "int" "long int" "long long")

if(INDEX_TYPE STREQUAL "long int")
  option(USE_CRBC_INDEX_LONG "Use long int for indexing" ON)
endif()

if(INDEX_TYPE STREQUAL "long long")
  option(USE_CRBC_INDEX_LONG_LONG "Use long int for indexing" ON)
endif()

# set options in the configuration folder
if(USE_OPENMP)
  set(USE_OPENMP_STRING 1)
endif()
if(USE_HDF5_RESTARTS)
  set(USE_HDF5_RESTARTS_STRING 1)
endif()
if(USE_CRBC_SCALAR_FLOAT)
  set(USE_OPENMP_STRING 1)
endif()
if(USE_CRBC_COEF_FLOAT)
  set(USE_CRBC_SCALAR_FLOAT_STRING 1)
endif()
if(USE_CRBC_INDEX_LONG)
  set(USE_CRBC_INDEX_LONG_STRING 1)
endif()
if(USE_CRBC_INDEX_LONG_LONG)
  set(USE_CRBC_INDEX_LONG_STRING 1)
endif()
if(ENABLE_EXTRA_DEBUG)
  add_definitions(-DCRBC_DEBUG)
endif()

# set hdf5 stuff if restarts are enabled
if(USE_HDF5_RESTARTS)
  
  # hdf5
  find_package(HDF5 REQUIRED)

  if(HDF5_FOUND)

      add_definitions(-DUSE_HDF5_RESTARTS)

      # with some compilers, the hdf5 include files aren't found when compiling
      # the c++ code. this seems to fix it ...
      include_directories(${HDF5_C_INCLUDE_DIR})

      set(OPT_LIBS ${OPT_LIBS} ${HDF5_LIBRARIES})

  endif()
endif()

# set OpenMP flags if threading enabled
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (yeecrbc_USE_OPENMP "${USE_OPENMP}")
  endif()
endif()

# build the library
ADD_SUBDIRECTORY(src)

# examples --- uncomment this to build the examples as part of the full build 
# process. The issue with this is that
# the source files we copy to the same location of the build are not the source
# files cmake uses to build the example so it might be confusing to users.
# ADD_SUBDIRECTORY(examples)

# so this is pretty much a hack to get the examples to be built in place so the
# source is next to the binaries. It's probably not a good way to do this, but it
# works better than the above.

# set the example install location so an user can change it if desired
set(example_dir "${CMAKE_CURRENT_BINARY_DIR}/examples" CACHE STRING "Set the install location for the examples.")

# copy the examples to build directory
file(COPY examples/3d_yee DESTINATION "${example_dir}")
file(COPY examples/2d_yee DESTINATION "${example_dir}")
file(COPY examples/wave_cxx DESTINATION "${example_dir}")
file(COPY examples/yee_mpi DESTINATION "${example_dir}")

# Create a CMakeList for the examples ...
set(yeecrbc_DIR "${CMAKE_INSTALL_PREFIX}/lib")
MESSAGE(STATUS "yee_crbcDIR" ${yeecrbc_DIR})
configure_file(examples/CMakeLists.txt.in ${example_dir}/CMakeLists.txt @ONLY)

# now add the examples directory in the build location ...
add_subdirectory(${example_dir} ${example_dir})


