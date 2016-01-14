/*
 * CUDA2DPotential.hpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */
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
#include "cuda_assert.hpp"
#include "cublas_assert.hpp"
#include "pot_base.hpp"
#include "CUDAFunctions.hpp"



#ifndef CUDA2DPOTENTIAL_HPP_
#define CUDA2DPOTENTIAL_HPP_
namespace slicepp {

class SliceBySliceCUDA2DPotential: public CPotential {

public:
	SliceBySliceCUDA2DPotential();
	SliceBySliceCUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc, cWaveConfPtr wc , PersistenceManagerPtr p);

	virtual ~SliceBySliceCUDA2DPotential();
	virtual void MakeSlices(superCellBoxPtr info);

protected:
//	virtual void SliceSetup();
	void progressCounter(int j, int jTot);
	void copyToDeviceInt(int *devdst, std::vector<int> src, int size);
	void copyToDeviceFloat(float_tt *devdst, std::vector<float_tt> src, int size);
	void copyDataToGPU(superCellBoxPtr info);

	virtual void AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ);
	virtual void ComputeAtomPotential(int znum){};
	virtual void SaveAtomicPotential(int znum);
	virtual af::array GetSlice(af::array t, unsigned idx);

	ComplexArray3D _t_host;
	CUDAFunctions *cf;
	cufftComplex *potential;
	float_tt *xyzPos_d, *occupancy_d;
	int *znums_d, *uniqueatoms;
	afcfloat *tafPtr;
	int slicePixels, numAtoms, numAtUnique;
	float_tt imPot;
	superCellBoxPtr _info;
	af::array xyzPos, occupancy, znums;
};
} /* namespace slicepp */

#endif /* CUDA2DPOTENTIAL_HPP_ */
