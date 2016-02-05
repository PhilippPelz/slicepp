/*
 * main.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: philipp
 */

#include <HDFFile.hpp>
#include <algorithm>
#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;
int main(int argc, char *argv[]){
	// dataset dimensions set at run time
	int NX = 5,  NY = 6,  NZ = 7;

	// allocate array using the "extents" helper.
	// This makes it easier to see how big the array is
	multi_array<float_tt, 3>  float_data(extents[NX][NY][NZ]);
	ComplexArray3D complex_data(extents[NX][NY][NZ]);
	// use resize to change size when necessary
	// float_data.resize(extents[NX + 5][NY + 4][NZ + 3]);


	// This is how you would fill the entire array with a value (e.g. 3.0)
	std::fill_n(float_data.data(), float_data.num_elements(), 3.0);
	std::fill_n(complex_data.data(), complex_data.num_elements(), Complex(3.0,3.0));

	// initialise the array to some variables
	for (int ii = 0; ii != NX; ii++)
	    for (int jj = 0; jj != NY; jj++)
	        for (int kk = 0; kk != NZ; kk++)
	            float_data[ii][jj][kk]  = ii + jj + kk;

	// write to HDF5 format
	boost::filesystem::path p("test.h5");
	HDFFile f(p,false);
	Group* g = f.CreateGroup(string("config"));
	f.SaveRealArray3D(float_data,string("3D"));
	f.SaveComplexArray3D(complex_data,string("3Dcomplex"));
	f.SaveAttributeDouble(string("3D"),string("doubleAtt"),1.337);
	f.SaveAttributeInt(string("3D"),string("intAtt"),1337);
	f.SaveAttributeString(string("3D"),string("intAtt"),string("stringA"));
	f.SaveString(string("cool"),string("stringval"));
	f.SaveInt(1337,string("cool"));
	f.SaveDouble(1.337,string("cool1"));
	f.SaveString(g->getId(),string("cool"),string("stringval"));
	f.SaveInt(g->getId(),1337,string("cool"));
	f.SaveDouble(g->getId(),1.337,string("cool1"));
	delete g;
//	write_hdf5(file, "doubleArray", float_data );
}

