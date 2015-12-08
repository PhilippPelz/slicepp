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
namespace QSTEM {

class CUDA2DPotential: public CPotential {

public:
	CUDA2DPotential();
	CUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc , PersistenceManagerPtr p);

	virtual ~CUDA2DPotential();
	virtual void MakeSlices(superCellBoxPtr info);

protected:
//	virtual void SliceSetup();
	void progressCounter(int j, int jTot);
	void copyToDeviceInt(int *devdst, std::vector<int> src, int size);
	void copyToDeviceFloat(float_tt *devdst, std::vector<float_tt> src, int size);

	virtual void AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ);
	virtual void ComputeAtomPotential(int znum){};
	virtual void SaveAtomicPotential(int znum);
	virtual void SavePotential();
	virtual af::array GetSlice(af::array t, unsigned idx);
	void initPotArrays();
	void ComputeAtPot(superCellBoxPtr info);
	cufftComplex* toCxPtr(af::array a);

	//pointer to device transmission matrix
	cufftComplex *_t_d_ptr;
	CUDAFunctions *_cf;

	int slicePixels, numAtoms, numAtUnique;
	af::array _V_elem;
	af::array _V_accum;
	std::map<int,af::array> _Z_atomDeltas;

	cufftComplex *_V_elem_ptr, *_V_accum_ptr;
};
} /* namespace QSTEM */

#endif /* CUDA2DPOTENTIAL_HPP_ */
