/*
 * HDFFile.cpp
 *
 *  Created on: Mar 30, 2015
 *      Author: philipp
 */

#include "HDFFile.hpp"
#include "hdf5_hl.h"

HDFFile::HDFFile() : _c_id(0){}
HDFFile::HDFFile(boost::filesystem::path p) {
	// TODO Auto-generated constructor stub
	/*used only to compute offsets */
	complex_t tmp;
	_c_id = H5Tcreate (H5T_COMPOUND, sizeof(tmp));
	H5Tinsert(_c_id, "r", 0, H5T_NATIVE_DOUBLE);
	H5Tinsert(_c_id, "i", sizeof(double), H5T_NATIVE_DOUBLE);

	_file = H5File(p.string(), H5F_ACC_TRUNC);
}

HDFFile::~HDFFile() {

	_file.close();
}
int HDFFile::SaveInt(int value,string dsetName){
	SaveInt(_file.getId(),value,dsetName);
}
int HDFFile::SaveDouble(double value,string dsetName){
	SaveDouble(_file.getId(),value,dsetName);
}
int HDFFile::SaveString(string value,string dsetName){
	SaveString(_file.getId(),value,dsetName);
}
int HDFFile::SaveInt(hid_t locID, int value,string dsetName){
	int RANK = 1;
	hsize_t     dims[RANK]={1};
	herr_t      status;
	return H5LTmake_dataset_int(locID,dsetName.c_str(),RANK,dims,&value);
}
int HDFFile::SaveDouble(hid_t locID,double value,string dsetName){
	int RANK = 1;
	hsize_t     dims[RANK]={1};
	herr_t      status;
	return H5LTmake_dataset_double(locID,dsetName.c_str(),RANK,dims,&value);
}
int HDFFile::SaveString(hid_t locID,string value,string dsetName){
	int RANK = 1;
	hsize_t     dims[RANK]={1};
	herr_t      status;
	return H5LTmake_dataset_string(locID,dsetName.c_str(),value.c_str());
}
int HDFFile::SaveAttributeInt(string dsetName, string attName, int data){
	int ATTR_SIZE = 1;
	return H5LTset_attribute_int(_file.getId(), dsetName.c_str(), attName.c_str(), &data, ATTR_SIZE);
}
int HDFFile::SaveAttributeDouble(string dsetName, string attName, double data){
	int ATTR_SIZE = 1;
	return H5LTset_attribute_double(_file.getId(), dsetName.c_str(), attName.c_str(), &data, ATTR_SIZE);
}
int HDFFile::SaveAttributeString(string dsetName, string attName, string data){
	int ATTR_SIZE = 1;
	return H5LTset_attribute_string(_file.getId(), dsetName.c_str(), attName.c_str(), data.c_str());
}
Group* HDFFile::CreateGroup(string groupName){
	return new Group( _file.createGroup(groupName.c_str()));
}
int HDFFile::SaveComplexArray2D(ComplexArray2D a, string dsetName){
	int RANK = 2;
	hsize_t     dims[RANK]={a.shape()[0],a.shape()[1]};
	herr_t      status;
	return H5LTmake_dataset(_file.getId(),dsetName.c_str(),RANK,dims,_c_id,reinterpret_cast<void *>(a.data()));
}
int HDFFile::SaveComplexArray3D(ComplexArray3D a, string dsetName){
	int RANK = 3;
	hsize_t     dims[RANK]={a.shape()[0],a.shape()[1],a.shape()[2]};
	herr_t      status;
	return H5LTmake_dataset(_file.getId(),dsetName.c_str(),RANK,dims,_c_id,(const void *)a.data());

}
int HDFFile::SaveRealArray2D(multi_array<double,2> a, string dsetName){
	int RANK = 2;
	hsize_t     dims[RANK]={a.shape()[0],a.shape()[1]};
	herr_t      status;
	return H5LTmake_dataset(_file.getId(),dsetName.c_str(),RANK,dims,H5T_NATIVE_DOUBLE,(const void *)a.data());
}
int HDFFile::SaveRealArray3D(multi_array<double,3> a, string dsetName){
	int RANK = 3;
	hsize_t     dims[RANK]={a.shape()[0],a.shape()[1],a.shape()[2]};
	herr_t      status;
	return H5LTmake_dataset(_file.getId(),dsetName.c_str(),RANK,dims,H5T_NATIVE_DOUBLE,(const void *)a.data());
}
