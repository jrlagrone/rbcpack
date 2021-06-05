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

/*

  This provides a somewhat cleaner interface to output to HDF5 files. I'm not 
  sure this is compatible with the parallel HDF5, but it should be reasonable
  to modify ...

*/

#ifndef HDF5_INTERFACE_H_
#define HDF5_INTERFACE_H_

extern "C" {
#include "hdf5.h"
}
#include <string>
#include <iostream>

namespace util_hdf5
{

  /// class to figure out data mappings for HDF5
  class DataMapH5
  {
    public: 
      // define a structure to get the appropriate alias to the HDF datatypes
      template <typename T>
      struct H5TypeMap  {
        static hid_t const h5_type;
      };

      /// function to return the HDF5 type
      template<typename T>
      hid_t getH5Type() {return H5TypeMap<T>::h5_type;};
  };

  /// Class to handle file output to HDF5
  class OutputHDF5 : private DataMapH5
  {

    public:

      OutputHDF5() {};

      /// constructor creates an HDF5 file of the given name (this will
      /// overwrite the file if it exists
      /// 
      /// \param[in] fname file name
      /// \param[in] compress use compression
      OutputHDF5(const std::string &fname, const bool compress = true);

      /// destructor --- closes file 
      ~OutputHDF5();

      /// method to close the file
      void close();

      /// method to add a data set
      ///
      /// \param[in] dset_name name for the data set, remember to include leading
      ///            '/'. If creating a dataset that is not at the root level,
      ///             the "path" must have already been created
      /// \param[in] data array of data
      /// \param[in] n length of data array 
      template<typename DataType>
      void addDataSet(const std::string &dset_name, const DataType *data, const size_t &n);

      /// method to add a group
      /// \param[in] gname group name
      void addGroup(const std::string &gname);

      /// method to add an attribute to an existing data set
      ///
      /// \param[in] dset_name  name of the data set
      /// \param[in] attr_name  name for the attribute
      /// \param[in] attr_data  data array for attribute values
      /// \param[in] n          length of data array
      template<typename DataType>
      void addAttribute(const std::string &dset_name, 
                       const std::string &attr_name, 
                       const DataType *attr_data, 
                       const size_t &n);

      /// method to add an attribute to an existing data set
      ///
      /// \param[in] g_name     name of the group
      /// \param[in] attr_name  name for the attribute
      /// \param[in] attr_data  data array for attribute values
      /// \param[in] n          length of data array
      template<typename DataType>
      void addAttributeToGroup(const std::string &g_name, 
                       const std::string &attr_name, 
                       const DataType *attr_data, 
                       const size_t &n);

    private:

      bool file_open;
      bool use_compression;
      hid_t       file_id;

  }; // end class


  /// Class to handle file input from HDF5
  class ReadHDF5 : private DataMapH5
  {

    public:

      /// constructor opens an HDF5 file of the given name 
      /// 
      /// \param[in] fname file name
      ReadHDF5(const std::string &fname);

      /// destructor --- closes file 
      ~ReadHDF5();

      /// method to close the file
      void close();

      /// method to get a data set
      ///
      /// \param[in] dset_name name for the data set, remember to include leading
      ///            '/'. If creating a dataset that is not at the root level,
      ///             the "path" must have already been created
      /// \param[in] data array of data
      template<typename DataType>
      void getDataSet(const std::string &dset_name, DataType *data);

      /// method to get an attribute from data set
      template<typename DataType>
      void getAttribute(const std::string &dset_name, 
                       const std::string &attr_name, 
                       DataType *attr_data);

      /// method to get an attribute from group
      template<typename DataType>
      void getAttributeFromGroup(const std::string &dset_name, 
                       const std::string &attr_name, 
                       DataType *attr_data);

    private:

      bool file_open;
      hid_t       file_id;

  }; // end class


/*------------------------------------------------------------------------------
                                 Output Defintions
------------------------------------------------------------------------------*/
template<typename DataType>
void OutputHDF5::addDataSet(const std::string &dset_name, const DataType *data, const size_t &n)
{

  hid_t       dataset_id, dataspace_id, plist_id;
  hsize_t     dims[1];
  hid_t dtype = getH5Type<DataType>();

  dims[0] = n;

  dataspace_id = H5Screate_simple(1, dims, NULL);

  plist_id  = H5Pcreate(H5P_DATASET_CREATE);

  if (use_compression) {

    // chunk size ... not really sure what a good value is ...
    hsize_t cdims[1];
    cdims[0] = (n < 1000) ? n : 1000;

    if (H5Pset_chunk (plist_id, 1, cdims) < 0)
       std::cerr << "OutputHDF5::addDataSet: error setting chunk size" << std::endl;
    if (H5Pset_deflate (plist_id, 6) < 0)
       std::cerr << "OutputHDF5::addDataSet: error setting deflate value" << std::endl;

  }

  dataset_id = H5Dcreate2(file_id, dset_name.c_str(), dtype, dataspace_id, \
                          H5P_DEFAULT, plist_id, H5P_DEFAULT);

  if (data != NULL)
    if (H5Dwrite(dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
      std::cerr << "OutputHDF5::addDataSet: error writing data set" << std::endl;

  if (H5Sclose(dataspace_id) < 0)
    std::cerr << "OutputHDF5::addDataSet: error closing data space" << std::endl;
  if (H5Dclose(dataset_id) < 0)
    std::cerr << "OutputHDF5::addDataSet: error closing data set" << std::endl;
  if (H5Pclose(plist_id) < 0)
    std::cerr << "OutputHDF5::addDataSet: error closing property list" << std::endl;

}

template<typename DataType>
void OutputHDF5::addAttribute(const std::string &dset_name, 
                       const std::string &attr_name, 
                       const DataType *attr_data, 
                       const size_t &n)
{

  hid_t dataset_id, attribute_id, dataspace_id; 
  hsize_t dims = n;
  hid_t dtype = getH5Type<DataType>();
  
  // open a data set
  dataset_id = H5Dopen2(file_id, dset_name.c_str(), H5P_DEFAULT);

  // create dataspace for attribute
  dataspace_id = H5Screate_simple(1, &dims, NULL);

  // create attribute
  attribute_id = H5Acreate2 (dataset_id, attr_name.c_str(), dtype, dataspace_id, 
                             H5P_DEFAULT, H5P_DEFAULT);

  if (H5Awrite(attribute_id, dtype, attr_data) < 0)
    std::cerr << "OutputHDF5::addAttribute: error writing attribute" << std::endl;

  if (H5Aclose(attribute_id) < 0)
    std::cerr << "OutputHDF5::addAttribute: error closing attribute" << std::endl;
  if (H5Sclose(dataspace_id) < 0)
    std::cerr << "OutputHDF5::addAttribute: error closing data space" << std::endl;
  if (H5Dclose(dataset_id) < 0)
    std::cerr << "OutputHDF5::addAttribute: error closing data set" << std::endl;
}

template<typename DataType>
void OutputHDF5::addAttributeToGroup(const std::string &g_name, 
                       const std::string &attr_name, 
                       const DataType *attr_data, 
                       const size_t &n)
{

  hid_t group_id, attribute_id, dataspace_id; 
  hsize_t dims = n;
  hid_t dtype = getH5Type<DataType>();
  
  // open a group
  group_id = H5Gopen2(file_id, g_name.c_str(), H5P_DEFAULT);

  // create dataspace for attribute
  dataspace_id = H5Screate_simple(1, &dims, NULL);

  // create attribute
  attribute_id = H5Acreate2 (group_id, attr_name.c_str(), dtype, dataspace_id, 
                             H5P_DEFAULT, H5P_DEFAULT);

  if (H5Awrite(attribute_id, dtype, attr_data) < 0)
    std::cerr << "OutputHDF5::addAttributeToGroup: error writing attribute" << std::endl;

  if (H5Aclose(attribute_id) < 0)
    std::cerr << "OutputHDF5::addAttributeToGroup: error closing attribute" << std::endl;
  if (H5Sclose(dataspace_id) < 0)
    std::cerr << "OutputHDF5::addAttributeToGroup: error closing data space" << std::endl;
  if (H5Gclose(group_id) < 0)
    std::cerr << "OutputHDF5::addAttributeToGroup: error closing group" << std::endl;
}

/*------------------------------------------------------------------------------
                                 Read Defintions
------------------------------------------------------------------------------*/
template<typename DataType>
void ReadHDF5::getDataSet(const std::string &dset_name, DataType *data)
{

  hid_t       dataset_id;  /* identifiers */
  hid_t dtype = getH5Type<DataType>();

  dataset_id = H5Dopen2(file_id, dset_name.c_str(), H5P_DEFAULT);

  if (H5Dread(dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
    std::cerr << "ReadHDF5::getDataSet: error reading data set" << std::endl;

  if (H5Dclose(dataset_id) < 0)
    std::cerr << "ReadHDF5::getDataSet: error closing data set" << std::endl;

}

template<typename DataType>
void ReadHDF5::getAttribute(const std::string &dset_name, 
                       const std::string &attr_name, 
                       DataType *attr_data)
{
  hid_t dataset_id, attribute_id; 
  hid_t dtype = getH5Type<DataType>();

  dataset_id = H5Dopen2(file_id, dset_name.c_str(), H5P_DEFAULT);

  attribute_id = H5Aopen(dataset_id, attr_name.c_str(), H5P_DEFAULT);

  if (H5Aread(attribute_id, dtype, attr_data) < 0)
    std::cerr << "ReadHDF5::getAttribute: error reading attribute" << std::endl;

  if (H5Aclose(attribute_id) < 0)
    std::cerr << "ReadHDF5::getAttribute: error closing attribute" << std::endl;
  if (H5Dclose(dataset_id) < 0)
    std::cerr << "ReadHDF5::getAttribute: error closing data set" << std::endl;
}

template<typename DataType>
void ReadHDF5::getAttributeFromGroup(const std::string &dset_name, 
                       const std::string &attr_name, 
                       DataType *attr_data)
{
  hid_t group_id, attribute_id; 
  hid_t dtype = getH5Type<DataType>();

  group_id = H5Gopen2(file_id, dset_name.c_str(), H5P_DEFAULT);

  attribute_id = H5Aopen(group_id, attr_name.c_str(), H5P_DEFAULT);

  if (H5Aread(attribute_id, dtype, attr_data) < 0)
    std::cerr << "ReadHDF5::getAttributeFromGroup: error reading attribute" << std::endl;

  if (H5Aclose(attribute_id) < 0)
    std::cerr << "ReadHDF5::getAttributeFromGroup: error closing attribute" << std::endl;
  if (H5Gclose(group_id) < 0)
    std::cerr << "ReadHDF5::getAttributeFromGroup: error closing group" << std::endl;

}

} // end namespace

#endif // HDF5_INTERFACE_H_
