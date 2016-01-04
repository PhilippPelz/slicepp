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

class CUDA2DPotential: public CPotential {

public:
	CUDA2DPotential();
	CUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc , PersistenceManagerPtr p);

	virtual ~CUDA2DPotential();
	virtual void MakeSlices(superCellBoxPtr info);

protected:
	void progressCounter(int j, int jTot);

	virtual void AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ);
	virtual void ComputeAtomPotential(int znum){};
	virtual void SaveAtomicPotential(int znum);
	virtual void SavePotential();
	virtual af::array GetSlice(af::array t, unsigned idx);
	void initPotArrays();
	void ComputeAtPot(superCellBoxPtr info);

	void fft(cufftComplex* V);
	void ifft(cufftComplex* V);
	void XplusequY(cufftComplex* X, cufftComplex* Y);

	//pointer to device transmission matrix
	cufftComplex *_t_d_ptr;
	cufftComplex *_V_elem_ptr;
	cufftComplex *_V_accum_ptr;

	CUDAFunctions *_cf;

	cufftHandle  _fftPlan;
	cublasHandle_t _cublasHandle;

	int _slicePixels, numAtoms, numAtUnique;

	std::map<int, cufftComplex*> _atomPot_d;
};
} /* namespace slicepp */

#endif /* CUDA2DPOTENTIAL_HPP_ */
