cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(wave_eq)

# may be needed to find the yeecrbc config file
# set(yeecrbc_DIR "/path/to/yeecrbc")

# Find the yeecrbc library
if(NOT TARGET yeecrbc)

  find_package(yeecrbc QUIET
    HINTS ${yeecrbc_DIR} ../ ../../ $ENV{yeecrbc_DIR})

  IF(NOT ${yeecrbc_FOUND})
    MESSAGE(FATAL_ERROR "\n"
      "*** Could not locate the Yee CRBC library ***\n\n"
      "You may want to either pass a flag -Dyeecrbc_DIR=/path/to/yeecrbc to cmake\n"
      "or set an environment variable \"yeecrbc_DIR\" that contains this path."
    )
  ENDIF()
ENDIF()

add_executable(wave_eq.x wave_eq.cpp)

target_link_libraries(wave_eq.x yeecrbc)
