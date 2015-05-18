/*
 * HDFFile.h
 *
 *  Created on: Mar 30, 2015
 *      Author: philipp
 */

#ifndef SRC_HDFFILE_H_
#define SRC_HDFFILE_H_

#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include "H5Cpp.h"
#include <cstdint>
//#include <algorithm>

#include "Complex.hpp"
using boost::multi_array;
using boost::extents;

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

using namespace std;

typedef struct {
	double r;   /*real part */
	double i;   /*imaginary part */
} complex_t;

typedef multi_array<Complex, 3> ComplexArray3D;
typedef multi_array<Complex, 2> ComplexArray2D;

template<typename T> struct get_hdf5_data_type
{   static H5::PredType type()
    {
        //static_assert(false, "Unknown HDF5 data type");
        return H5::PredType::NATIVE_DOUBLE;
    }
};
template<> struct get_hdf5_data_type<char>                  {   H5::IntType type    {   H5::PredType::NATIVE_CHAR       };  };
//template<> struct get_hdf5_data_type<unsigned char>       {   H5::IntType type    {   H5::PredType::NATIVE_UCHAR      };  };
//template<> struct get_hdf5_data_type<short>               {   H5::IntType type    {   H5::PredType::NATIVE_SHORT      };  };
//template<> struct get_hdf5_data_type<unsigned short>      {   H5::IntType type    {   H5::PredType::NATIVE_USHORT     };  };
//template<> struct get_hdf5_data_type<int>                 {   H5::IntType type    {   H5::PredType::NATIVE_INT        };  };
//template<> struct get_hdf5_data_type<unsigned int>        {   H5::IntType type    {   H5::PredType::NATIVE_UINT       };  };
//template<> struct get_hdf5_data_type<long>                {   H5::IntType type    {   H5::PredType::NATIVE_LONG       };  };
//template<> struct get_hdf5_data_type<unsigned long>       {   H5::IntType type    {   H5::PredType::NATIVE_ULONG      };  };
template<> struct get_hdf5_data_type<long long>             {   H5::IntType type    {   H5::PredType::NATIVE_LLONG      };  };
template<> struct get_hdf5_data_type<unsigned long long>    {   H5::IntType type    {   H5::PredType::NATIVE_ULLONG     };  };
template<> struct get_hdf5_data_type<int8_t>                {   H5::IntType type    {   H5::PredType::NATIVE_INT8       };  };
template<> struct get_hdf5_data_type<uint8_t>               {   H5::IntType type    {   H5::PredType::NATIVE_UINT8      };  };
template<> struct get_hdf5_data_type<int16_t>               {   H5::IntType type    {   H5::PredType::NATIVE_INT16      };  };
template<> struct get_hdf5_data_type<uint16_t>              {   H5::IntType type    {   H5::PredType::NATIVE_UINT16     };  };
template<> struct get_hdf5_data_type<int32_t>               {   H5::IntType type    {   H5::PredType::NATIVE_INT32      };  };
template<> struct get_hdf5_data_type<uint32_t>              {   H5::IntType type    {   H5::PredType::NATIVE_UINT32     };  };
template<> struct get_hdf5_data_type<int64_t>               {   H5::IntType type    {   H5::PredType::NATIVE_INT64      };  };
template<> struct get_hdf5_data_type<uint64_t>              {   H5::IntType type    {   H5::PredType::NATIVE_UINT64     };  };

class HDFFile {
public:
	HDFFile();
	HDFFile(boost::filesystem::path p);
	virtual ~HDFFile();
	int SaveComplexArray2D(ComplexArray2D a, string dsetName);
	int SaveComplexArray3D(ComplexArray3D a, string dsetName);
	int SaveRealArray2D(multi_array<double,2> a, string dsetName);
	int SaveRealArray3D(multi_array<double,3> a, string dsetName);
	int SaveInt(int value,string dsetName);
	int SaveDouble(double value,string dsetName);
	int SaveString(string value,string dsetName);
	int SaveInt(hid_t locID, int value,string dsetName);
	int SaveDouble(hid_t locID,double value,string dsetName);
	int SaveString(hid_t locID,string value,string dsetName);
	int SaveAttributeInt(string dsetName, string attName, int data);
	int SaveAttributeDouble(string dsetName, string attName, double data);
	int SaveAttributeString(string dsetName, string attName, string data);
	Group* CreateGroup(string groupName);

protected:
	hid_t _c_id;
	H5File _file;
};

#endif /* SRC_HDFFILE_H_ */
