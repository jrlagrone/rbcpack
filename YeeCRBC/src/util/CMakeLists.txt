cmake_minimum_required(VERSION 2.8)

project(h5util)

add_library(h5util STATIC hdf5_interface.cpp hdf5_interface.hpp)
