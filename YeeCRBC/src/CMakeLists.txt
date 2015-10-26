# This file creates the CRBCYee library and optionally builds the examples

###
# set the project name and verions
project(yeecrbc)
set(yeecrbc_VERSION "1.0.0")

set(HDF5FILES "")

# set hdf5 stuff if restarts are enabled
if(USE_HDF5_RESTARTS)
  if(HDF5_FOUND)
    # define a list of files that depend on HDF5 routines
    set(HDF5FILES util/hdf5_interface.cpp)
    set(TEMPLATE_SOURCE util/hdf5_interface.hpp)
  endif()
endif()

# specify the templated header files --- these are source files, but we also want
# to include them in the install so they may be used for other applications than
# what we permit in the API
set(TEMPLATE_SOURCE
    ${TEMPLATE_SOURCE}
    grid/box.hpp
    grid/box_iterator.hpp
    grid/grid_array.hpp
    grid/point.hpp
    boundary_properties.hpp
    crbc_data.hpp
    crbc_updates.hpp
    recursions.hpp
    wave_equations.hpp)

# now set all of the sources
set(SOURCE 2d_crbc_api.cpp
           3d_yee_api.cpp 
           3d_yee_api_cxx.cpp
           optimal_cosines/optimal_cosines.c
           ${TEMPLATE_SOURCE}
           ${HDF5FILES})

# build the configuration header
configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h @ONLY)

# set the API header files
set(HEADERS config.h
            2d_crbc_api.h
            3d_yee_crbc_api.h 
            3d_yee_api.hpp 
            optimal_cosines/optimal_cosines.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/grid
                    ${CMAKE_CURRENT_SOURCE_DIR}/util
                    ${CMAKE_CURRENT_SOURCE_DIR}/optimal_cosines
                    ${CMAKE_CURRENT_BINARY_DIR}/src)

# now create the library targets
add_library(core OBJECT ${SOURCE} ${HEADERS})
set_target_properties(core PROPERTIES CXX_STANDARD 11 )
set_target_properties(core PROPERTIES CXX_STANDARD_REQUIRED ON )

if(BUILD_SHARED)
  set_property(TARGET core PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library(yeecrbc SHARED $<TARGET_OBJECTS:core>)
else()
  add_library(yeecrbc STATIC $<TARGET_OBJECTS:core>)
endif()



# set version
set_property(TARGET yeecrbc PROPERTY VERSION ${yeecrbc_VERSION})

# link optional libraries
target_link_libraries(yeecrbc ${OPT_LIBS})

# SET_TARGET_PROPERTIES(yeecrbc PROPERTIES LINKER_LANGUAGE CXX)

# Global includes. Used by all targets
# Note:
#   * header location in project: Foo/Source/foo/Bar.hpp
#   * header can be included by C++ code `#include <foo/Bar.hpp>`
#   * header location in project: ${CMAKE_CURRENT_BINARY_DIR}/bar_export.hpp
#   * header can be included by: `#include <bar_export.hpp>`
target_include_directories(
    yeecrbc PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/grid>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/util>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/optimal_cosines>"
)

####
# Install

# Layout.
#   * <prefix>/lib/cmake/<PROJECT-NAME>
#   * <prefix>/lib/
#   * <prefix>/include/
set(config_install_dir "lib/cmake/${PROJECT_NAME}")
set(include_install_dir "include/crbc")

# Configuration
set(version_config "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")


set(yeecrbc_INCLUDE_DIR $<INSTALL_INTERFACE:/include>/yeecrbc)
set(yeecrbc_LIB_DIR  $<INSTALL_INTERFACE:/lib>)

# Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Note: PROJECT_VERSION is used as a VERSION
write_basic_package_version_file("${version_config}"
                 VERSION "${yeecrbc_VERSION}" 
                 COMPATIBILITY SameMajorVersion
)

# Configure '<PROJECT-NAME>Config.cmake'
# Note: variable 'targets_export_name' used
configure_file("cmake/Config.cmake.in" "${project_config}" @ONLY)


#   * <prefix>/lib/libbar.a
#   * <prefix>/lib/libbaz.a
#   * header location after install: <prefix>/include/foo/Bar.hpp
#   * headers can be included by C++ code `#include <foo/Bar.hpp>`
install(TARGETS yeecrbc 
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib"
        RUNTIME DESTINATION "bin"
        INCLUDES DESTINATION "${include_install_dir}"
)

# install header files
install(FILES boundary_properties.hpp 
              crbc_data.hpp
              crbc_updates.hpp
              recursions.hpp
              wave_equations.hpp
              2d_crbc_api.h
              3d_yee_crbc_api.h 
              3d_yee_api.hpp 
              config.h
        DESTINATION 
              "${include_install_dir}"
        COMPONENT
               Devel
        )

install(DIRECTORY grid DESTINATION "${include_install_dir}"
        FILES_MATCHING PATTERN "*.hpp")

install(DIRECTORY util DESTINATION "${include_install_dir}"
        FILES_MATCHING PATTERN "*.hpp")

install(DIRECTORY optimal_cosines DESTINATION "${include_install_dir}"
        FILES_MATCHING PATTERN "*.h")

export(TARGETS yeecrbc 
       FILE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}LibraryDepends.cmake)

install(TARGETS yeecrbc EXPORT ${CMAKE_PROJECT_NAME}LibraryDepends 
        ARCHIVE DESTINATION "lib"
        LIBRARY DESTINATION "lib")

# Config

install(
    EXPORT "${targets_export_name}"
    NAMESPACE "${namespace}"
    DESTINATION "${config_install_dir}"
)

install(
    FILES "${project_config}" "${version_config}"
    DESTINATION "lib"
)

