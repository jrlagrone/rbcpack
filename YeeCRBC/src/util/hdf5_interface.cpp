/*
  Copyright 2015 John LaGrone

  This file is part of the Yee CRBC Library from RBCPACK.

  The Yee CRBC Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or (at your 
  option) any later version.

  The Yee CRBC Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the Yee CRBC Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hdf5_interface.hpp"

namespace util_hdf5
{

  // add explicit specialization for all required types:  
  template <> const hid_t DataMapH5::H5TypeMap<float>::h5_type = H5T_NATIVE_FLOAT;
  template <> const hid_t DataMapH5::H5TypeMap<double>::h5_type = H5T_NATIVE_DOUBLE;
  template <> const hid_t DataMapH5::H5TypeMap<long double>::h5_type = H5T_NATIVE_LDOUBLE;
  template <> const hid_t DataMapH5::H5TypeMap<short>::h5_type = H5T_NATIVE_SHORT;
  template <> const hid_t DataMapH5::H5TypeMap<unsigned short>::h5_type = H5T_NATIVE_USHORT;
  template <> const hid_t DataMapH5::H5TypeMap<int>::h5_type = H5T_NATIVE_INT;
  template <> const hid_t DataMapH5::H5TypeMap<unsigned int>::h5_type = H5T_NATIVE_UINT;
  template <> const hid_t DataMapH5::H5TypeMap<long>::h5_type = H5T_NATIVE_LONG;
  template <> const hid_t DataMapH5::H5TypeMap<unsigned long>::h5_type = H5T_NATIVE_ULONG;
  template <> const hid_t DataMapH5::H5TypeMap<long long>::h5_type = H5T_NATIVE_LLONG;
  template <> const hid_t DataMapH5::H5TypeMap<unsigned long long>::h5_type = H5T_NATIVE_ULLONG;

  /*------------------------------------------------------------------------------
                                 Output Defintions
------------------------------------------------------------------------------*/

// constructor
OutputHDF5::OutputHDF5(const std::string &fname, const bool compress)
{

  use_compression = compress;
  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  file_open = true;

}

// deconstructor
OutputHDF5::~OutputHDF5()
{
  close();
}

void OutputHDF5::close()
{
  if (file_open) {
    if (H5Fclose(file_id) < 0)
       std::cerr << "OutputHDF5::close: error closing file" << std::endl;
    file_open = false;
  }
}

void OutputHDF5::addGroup(const std::string &gname)
{
  hid_t  group_id;

  group_id = H5Gcreate2(file_id, gname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (H5Gclose(group_id) < 0)
    std::cerr << "OutputHDF5::addGroup: error closing group" << std::endl;
}

/*------------------------------------------------------------------------------
                                 Read Defintions
------------------------------------------------------------------------------*/
ReadHDF5::ReadHDF5(const std::string &fname)
{
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  file_open = true;
}

// deconstructor
ReadHDF5::~ReadHDF5()
{
  close();
}

void ReadHDF5::close()
{
  if (file_open) {
    if (H5Fclose(file_id) < 0)
       std::cerr << "ReadHDF5::close: error closing file" << std::endl;
    file_open = false;
  }
}


}
