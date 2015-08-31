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
#include "float.h"
#include <math.h>
#include <string.h>
#include <cufft.h>
#include <curand.h>
#include <cublas_v2.h>
#include "curand_kernel.h"
#include <vector>
#include "cuda_assert.hpp"
#include "cublas_assert.hpp"
#include "cufft_assert.hpp"
#include "projectedPotential.hpp"
#include "arrayfire.h"

namespace QSTEM {
class CUDAFunctions{
public:
	CUDAFunctions();
	void phaseGrating(cufftComplex* V_d, int nAt, int nZ, float_tt* xyz_d, float_tt imPot, int* Z_d, int* Zlist, float_tt*  occ_d, int s, int nx, int ny, int nSlices, float_tt dx, float_tt dy, float_tt dz );
	void printPotArray(cufftComplex* V_d, int nx, int ny);
	void printFloatArray(float_tt* f, int nx, int ny, int offset);
	void printIntArray(int *p, int size);
	void initPotArrays(int slicePixels);
	void unlockArrays();
protected:
	int myGSize(int size);
	int myBSize(int size);

	af::array V1, V2, V3;
	cufftComplex *V1_d, *V2_d, *V3_d;
};
__global__ void squareAtoms_d (  cufftComplex* V, int nAt, int *Z, int Z0, float_tt *xyz, float_tt imPot, float_tt *occ, int s, int nx, int ny, int nslice, float_tt dx, float_tt dy, float_tt dz );
__global__ void divideBySinc ( cufftComplex* V, int nx, int ny, float_tt PI);
__global__ void multiplyWithProjectedPotential_d ( cufftComplex* V1, cufftComplex* V2, int nx, int ny);
__global__ void initialValues ( cuComplex* V, int size, float_tt initRe, float_tt initIm );
__global__ void potential2Transmission ( cufftComplex* t, cufftComplex* V, int size );
__device__ int mySignum_d( float x );
} /* namespace QSTEM */
#endif /* CUDAFUNCTIONS_HPP_ */
