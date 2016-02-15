/*
 * HDFFile.cpp
 *
 *  Created on: Mar 30, 2015
 *      Author: philipp
 */

#include "HDFFile.hpp"
#include "hdf5_hl.h"
#include "H5Cpp.h"

HDFFile::HDFFile() :
		_complex_id(0) {
}
HDFFile::HDFFile(boost::filesystem::path p, bool saveComplexAsFloat2) :
		_saveComplexAsFloat2(saveComplexAsFloat2) {
	complex_t tmp;
	_complex_id = H5Tcreate(H5T_COMPOUND, sizeof(tmp));
	const H5std_string MEMBER1( "r" );
	const H5std_string MEMBER2( "i" );
	_complex_type = CompType(sizeof(complex_t));

#if FLOAT_PRECISION == 1
	H5Tinsert(_complex_id, "r", 0, H5T_NATIVE_FLOAT);
	H5Tinsert(_complex_id, "i", sizeof(float_tt), H5T_NATIVE_FLOAT);
	_complex_type.insertMember(MEMBER1, HOFFSET(complex_t, r), PredType::NATIVE_FLOAT);
	_complex_type.insertMember(MEMBER2, HOFFSET(complex_t, i), PredType::NATIVE_FLOAT);
#else
	H5Tinsert(_complex_id, "r", 0, H5T_NATIVE_DOUBLE);
	H5Tinsert(_complex_id, "i", sizeof(float_tt), H5T_NATIVE_DOUBLE);
	_complex_type.insertMember( MEMBER1, HOFFSET(complex_t, r), PredType::NATIVE_DOUBLE);
	_complex_type.insertMember( MEMBER2, HOFFSET(complex_t, i), PredType::NATIVE_DOUBLE);
#endif

	_file = H5File(p.string(), H5F_ACC_TRUNC);

#if FLOAT_PRECISION == 1
	_float_id = H5Tcopy(H5T_NATIVE_FLOAT);
#else
	_float_id = H5Tcopy(H5T_NATIVE_DOUBLE);
#endif
	H5Tlock(_float_id);
}
void HDFFile::CreateComplexDataSet(std::string path, std::vector<unsigned> dims) {
	CreateDataSet(path, _complex_type, dims);
}
void HDFFile::CreateFloatDataSet(std::string path, std::vector<unsigned> dims) {
	FloatType t(PredType::NATIVE_FLOAT);
	CreateDataSet(path, t, dims);
}
void HDFFile::SaveComplexSlice(std::string path, unsigned size_x, unsigned size_y, unsigned slice, multi_array<Complex, 2> a, bool float2) {
	if (float2){
		FloatType t(PredType::NATIVE_FLOAT);
		DataSlabIOFloat2(false, t, (void*) a.data(), path, size_x, size_y, slice);
	}
	else
		DataSlabIO(false, _complex_type, (void*) a.data(), path, size_x, size_y, slice);
}
void HDFFile::SaveFloatSlice(std::string path, unsigned size_x, unsigned size_y, unsigned slice, multi_array<float_tt, 2> a) {
	FloatType t(PredType::NATIVE_FLOAT);
	DataSlabIO(false, t, (void*) a.data(), path, size_x, size_y, slice);
}
void HDFFile::CreateDataSet(std::string path, DataType type, std::vector<unsigned> dims) {
	std::vector<hsize_t> d;
	for (int i = 0; i < dims.size(); i++)
		d.push_back((hsize_t) dims[i]);
	DataSpace space(d.size(), d.data());
	DSetCreatPropList plist;
	auto ds = _file.createDataSet(path.c_str(), type, space, plist);
	ds.close();
	space.close();
}
//start: a starting location for the hyperslab. In the example start is (0,1).
//stride: the number of elements to separate each element or block to be selected. In the example stride is (4,3).
//If the stride parameter is set to NULL, the stride size defaults to 1 in each dimension.
//count: the number of elements or blocks to select along each dimension. In the example, count is (2,4).
//block: the size of the block selected from the dataspace. In the example, block is (3,2). If the block parameter is set to
//NULL, the block size defaults to a single element in each dimension, as if the block array was set to all 1s.
void HDFFile::DataSlabIO(bool read, DataType datatype, void *arr, std::string path, unsigned size_x, unsigned size_y, unsigned slice) {
	std::vector<hsize_t> memory_dims = { size_x, size_y };
	std::vector<hsize_t> start = { slice, 0, 0 };
	std::vector<hsize_t> count = { 1, size_x, size_y };
	std::vector<hsize_t> stride = { size_y * size_x, size_x, 1 };
	hsize_t* block = NULL;
	DSetMemXferPropList props;

	auto set = _file.openDataSet(path.c_str());
	auto fspace = set.getSpace();
	auto memspace = DataSpace(memory_dims.size(), memory_dims.data());
	fspace.selectHyperslab(H5S_SELECT_SET, count.data(), start.data(), NULL, block);
	if (read) {
		set.read(arr, datatype, memspace, fspace, props);
	} else {
//		printf("Writing %s slice %d", path.c_str(), slice);
		set.write(arr, datatype, memspace, fspace, props);
	}
}
void HDFFile::DataSlabIOFloat2(bool read, DataType datatype, void *arr, std::string path, unsigned size_x, unsigned size_y, unsigned slice) {
	std::vector<hsize_t> memory_dims = { size_x, size_y, 2 };
	std::vector<hsize_t> start = { slice, 0, 0, 0 };
	std::vector<hsize_t> count = { 1, size_x, size_y, 2 };
	std::vector<hsize_t> stride = { size_y * size_x * 2, size_x * 2, 1 };
	hsize_t* block = NULL;
	DSetMemXferPropList props;

	auto set = _file.openDataSet(path.c_str());
	auto fspace = set.getSpace();
	auto memspace = DataSpace(memory_dims.size(), memory_dims.data());
	fspace.selectHyperslab(H5S_SELECT_SET, count.data(), start.data(), NULL, block);
	if (read) {
		set.read(arr, datatype, memspace, fspace, props);
	} else {
//		printf("Writing %s slice %d", path.c_str(), slice);
		set.write(arr, datatype, memspace, fspace, props);
	}
}
HDFFile::~HDFFile() {

	_file.close();
}
int HDFFile::SaveInt(int value, string dsetName) {
	SaveInt(_file.getId(), value, dsetName);
}
int HDFFile::SaveDouble(float_tt value, string dsetName) {
	SaveDouble(_file.getId(), value, dsetName);
}
int HDFFile::SaveString(string value, string dsetName) {
	SaveString(_file.getId(), value, dsetName);
}
int HDFFile::SaveInt(hid_t locID, int value, string dsetName) {
	int RANK = 1;
	hsize_t dims[1] = { 1 };
	herr_t status;
	return H5LTmake_dataset_int(locID, dsetName.c_str(), RANK, dims, &value);
}
int HDFFile::SaveDouble(hid_t locID, float_tt value, string dsetName) {
	int RANK = 1;
	hsize_t dims[1] = { 1 };
	herr_t status;
#if(FLOAT_PRECISION == 1)
	return H5LTmake_dataset_float(locID, dsetName.c_str(), RANK, dims, &value);
#else
	return H5LTmake_dataset_double(locID,dsetName.c_str(),RANK,dims,&value);
#endif
}
int HDFFile::SaveString(hid_t locID, string value, string dsetName) {
	int RANK = 1;
	hsize_t dims[1] = { 1 };
	herr_t status;
	return H5LTmake_dataset_string(locID, dsetName.c_str(), value.c_str());
}
int HDFFile::SaveAttributeInt(string dsetName, string attName, int data) {
	int ATTR_SIZE = 1;
	return H5LTset_attribute_int(_file.getId(), dsetName.c_str(), attName.c_str(), &data, ATTR_SIZE);
}
int HDFFile::SaveAttributeDouble(string dsetName, string attName, float_tt data) {
	int ATTR_SIZE = 1;
#if(FLOAT_PRECISION == 1)
	return H5LTset_attribute_float(_file.getId(), dsetName.c_str(), attName.c_str(), &data, ATTR_SIZE);
#else
	return H5LTset_attribute_double(_file.getId(), dsetName.c_str(), attName.c_str(), &data, ATTR_SIZE);
#endif
}
int HDFFile::SaveAttributeString(string dsetName, string attName, string data) {
	int ATTR_SIZE = 1;
	return H5LTset_attribute_string(_file.getId(), dsetName.c_str(), attName.c_str(), data.c_str());
}
Group* HDFFile::CreateGroup(string groupName) {
	return new Group(_file.createGroup(groupName.c_str()));
}
int HDFFile::SaveComplexArray2D(ComplexArray2D a, string dsetName) {
	int RANK = 2;
	hsize_t dims[2] = { a.shape()[0], a.shape()[1] };
	hsize_t dims2[3] = { a.shape()[0], a.shape()[1], 2 };
	herr_t status;
	if (_saveComplexAsFloat2)
		return H5LTmake_dataset(_file.getId(), dsetName.c_str(), RANK + 1, dims2, H5T_NATIVE_FLOAT, reinterpret_cast<const void *>(a.data()));
	else
		return H5LTmake_dataset(_file.getId(), dsetName.c_str(), RANK, dims, _complex_id, reinterpret_cast<const void *>(a.data()));
}
int HDFFile::SaveComplexArray3D(ComplexArray3D a, string dsetName) {
	int RANK = 3;
	hsize_t dims[3] = { a.shape()[0], a.shape()[1], a.shape()[2] };
	hsize_t dims2[4] = { a.shape()[0], a.shape()[1], a.shape()[2], 2 };
	herr_t status;
	if (_saveComplexAsFloat2)
		return H5LTmake_dataset(_file.getId(), dsetName.c_str(), RANK + 1, dims2, H5T_NATIVE_FLOAT, reinterpret_cast<const void *>(a.data()));
	else
		return H5LTmake_dataset(_file.getId(), dsetName.c_str(), RANK, dims, _complex_id, reinterpret_cast<const void *>(a.data()));
}
int HDFFile::SaveRealArray2D(multi_array<float_tt, 2> a, string dsetName) {
	int RANK = 2;
	hsize_t dims[2] = { a.shape()[0], a.shape()[1] };
	herr_t status;
#if(FLOAT_PRECISION == 1)
	return H5LTmake_dataset(_file.getId(), dsetName.c_str(), RANK, dims, H5T_NATIVE_FLOAT, (const void *) a.data());
#else
	return H5LTmake_dataset(_file.getId(),dsetName.c_str(),RANK,dims,H5T_NATIVE_DOUBLE,(const void *)a.data());
#endif
}
int HDFFile::SaveRealArray1D(multi_array<float_tt, 1> a, string dsetName) {
	int RANK = 1;
	hsize_t dims[RANK] = { a.shape()[0] };
	herr_t status;
#if(FLOAT_PRECISION == 1)
	return H5LTmake_dataset(_file.getId(), dsetName.c_str(), RANK, dims, H5T_NATIVE_FLOAT, (const void *) a.data());
#else
	return H5LTmake_dataset(_file.getId(),dsetName.c_str(),RANK,dims,H5T_NATIVE_DOUBLE,(const void *)a.data());
#endif
}
int HDFFile::SaveRealArray3D(multi_array<float_tt, 3> a, string dsetName) {
	int RANK = 3;
	hsize_t dims[3] = { a.shape()[0], a.shape()[1], a.shape()[2] };
	herr_t status;
#if(FLOAT_PRECISION == 1)
	return H5LTmake_dataset(_file.getId(), dsetName.c_str(), RANK, dims, H5T_NATIVE_FLOAT, (const void *) a.data());
#else
	return H5LTmake_dataset(_file.getId(),dsetName.c_str(),RANK,dims,H5T_NATIVE_DOUBLE,(const void *)a.data());
#endif
}
