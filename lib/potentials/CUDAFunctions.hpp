/*
 * CUDAFunctions.hpp
 *
 *  Created on: Aug 3, 2015
 *      Author: wenxuan
 */

#ifndef CUDAFUNCTIONS_HPP_
#define CUDAFUNCTIONS_HPP_



#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <cufft.h>
#include <curand.h>
#include <cublas_v2.h>
#include "curand_kernel.h"
#include <vector>
#include <map>
#include "cuda_assert.hpp"
#include "cublas_assert.hpp"
#include "cufft_assert.hpp"
#include "projectedPotential.hpp"
#include "arrayfire.h"
#include "stemtypes_fftw3.hpp"
#include "config_IO/ConfigReader.hpp"
#include "af/cuda.h"

namespace slicepp {
class CUDAFunctions{
public:
	CUDAFunctions(superCellBoxPtr info, cModelConfPtr mc);
//	void GetPhaseGrating(cufftComplex* V_d, int s,  std::map<int, af::array> & atomPot);

	void GetSincAtomicPotential( cufftComplex* V, int Z);
	void GetAtomicPotential( cufftComplex* V, int Z);
	void GetAtomDeltaFunctions(cufftComplex* V, int Z,int slice);
	void PotentialToTransmission(cufftComplex* pot, cufftComplex* trans);
	void cmul(cufftComplex* a1, cufftComplex* a2);
	void SetComplex2D(cufftComplex* a, float real, float imag);
	void SetComplex3D(cufftComplex* a, float real, float imag);
	void printPotArray(cufftComplex* V_d, int nx, int ny);
	void printFloatArray(float_tt* f, int nx, int ny, int offset);
	void printIntArray(int *p, int size);
	void initPotArrays(int pix);
	void initArrays();
	void releaseArrays();
	void unlockArrays();
	void InitArrays();
protected:
	int myGSize(int size);
	int myBSize(int size);

	af::array _V_elem, _V_atom, _v_accum;
	cufftComplex *_V_elem_ptr, *_V_atom_ptr, *_V_accum_ptr;
	cModelConfPtr _mc;
	superCellBoxPtr _info;
	af::array xyzPos, occupancy, znums;
	float_tt *xyzPos_d, *occupancy_d;
	int *znums_d;
	int gS2D, slicePixels;
};
__global__ void putAtomDeltas (  cufftComplex* V, int nAt, int *Z, int Z0, float_tt *xyz, float_tt imPot, float_tt *occ, int s, int nx, int ny, int nslice, float_tt dx, float_tt dy, float_tt dz );
__global__ void divideBySinc ( cufftComplex* V, int nx, int ny, float_tt PI);
__global__ void multiplyWithProjectedPotential_d ( cufftComplex* V1, cufftComplex* V2, int nx, int ny);
__global__ void initialValues ( cuComplex* V, int size, float_tt initRe, float_tt initIm );
__global__ void potential2Transmission ( cufftComplex* t, cufftComplex* V, int size );
__device__ int sign( float x );
} /* namespace slicepp */
#endif /* CUDAFUNCTIONS_HPP_ */
