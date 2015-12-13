/*
 * CUDA2DPotential.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */

#include "CUDA2DPotential.hpp"
#include <stdio.h>
//#include <thrust/fill.h>
//#include <thrust/copy.h>

namespace QSTEM {
CUDA2DPotential::CUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) :
		CPotential(mc, oc, p) {
	cufft_assert(cufftPlan2d(&_fftPlan, _mc->nx, _mc->ny, CUFFT_C2C));
	cublas_assert(cublasCreate ( &_cublasHandle));
}

CUDA2DPotential::~CUDA2DPotential() {
	delete _cf;
	cublas_assert(cublasDestroy(_cublasHandle));
	cufft_assert(cufftDestroy(_fftPlan));
}
void CUDA2DPotential::initPotArrays() {
	_slicePixels = _mc->nx * _mc->ny;

	cuda_assert(cudaMalloc((void**) &_t_d_ptr, _mc->nSlices * _slicePixels * sizeof(cufftComplex)));
	cuda_assert(cudaMalloc((void**) &_V_elem_ptr, _slicePixels * sizeof(cufftComplex)));
	cuda_assert(cudaMalloc((void**) &_V_accum_ptr, _slicePixels * sizeof(cufftComplex)));

	_cf->SetComplex3D(_t_d_ptr, 0.f, 0.f);
//	cuda_assert(cudaDeviceSynchronize());
}
void CUDA2DPotential::ComputeAtPot(superCellBoxPtr info) {
	for (int Z : info->uniqueZ) {
		cufftComplex* _atomPot_d_ptr;
		cudaMalloc((void**)&_atomPot_d_ptr, _slicePixels*sizeof(cufftComplex));
		_atomPot_d[Z] = _atomPot_d_ptr;
//		BOOST_LOG_TRIVIAL(info)<< format("_atomPot_d_ptr: %d") % _atomPot_d_ptr;
		_cf->GetSincAtomicPotential(_atomPot_d_ptr, Z);

		if (_oc->SaveAtomicPotential) {
			cuda_assert(cudaDeviceSynchronize());
			SaveAtomicPotential(Z);
		}
	}
}
void CUDA2DPotential::fft(cufftComplex* V) {
	cufft_assert(cufftExecC2C( _fftPlan, V, V, CUFFT_FORWARD ));
}
void CUDA2DPotential::ifft(cufftComplex* V) {
	cufft_assert(cufftExecC2C(_fftPlan, V, V, CUFFT_INVERSE ));
}
void CUDA2DPotential::XplusequY(cufftComplex* X, cufftComplex* Y) {
	cufftComplex alpha;
	alpha.x = 1.0;
	alpha.y = 0.0;
	cublas_assert(cublasCaxpy(_cublasHandle, _slicePixels, &alpha, Y, 1, X, 1));
}
void CUDA2DPotential::MakeSlices(superCellBoxPtr info) {

	_cf = new CUDAFunctions(info, _mc);
	_cf->initArrays();
	initPotArrays();
	ComputeAtPot(info);
	af::timer time = af::timer::start();
	BOOST_LOG_TRIVIAL(info)<< "Calculating potential ...";

	for (int islice = 0; islice < _mc->nSlices; islice++) {
		progressCounter(islice, _mc->nSlices);
		_cf->SetComplex2D(_V_accum_ptr, 0.f, 0.f);
//		cuda_assert(cudaDeviceSynchronize());
		for (int& Z : info->uniqueZ) {
//			_V_elem.lock();
			_cf->SetComplex2D(_V_elem_ptr, 0.f, 0.f);
//			cuda_assert(cudaDeviceSynchronize());
			_cf->GetAtomDeltaFunctions(_V_elem_ptr, Z, islice);
//			cuda_assert(cudaDeviceSynchronize());

			if (_oc->SaveAtomDeltas) {
				_persist->SaveAtomDelta(_V_elem_ptr, islice, Z);
			}

			fft(_V_elem_ptr);
			_cf->cmul(_V_elem_ptr, _atomPot_d[Z]);
			ifft(_V_elem_ptr);

			if (_oc->SaveAtomConv) {
				_persist->SaveAtomConv(_V_elem_ptr, islice, Z);
			}

			XplusequY(_V_accum_ptr, _V_elem_ptr);
		}
		_cf->PotentialToTransmission(&_t_d_ptr[islice * _slicePixels], _V_accum_ptr);
//		cuda_assert(cudaDeviceSynchronize());
	}
	auto elapsed = af::timer::stop(time) * 1000;
	BOOST_LOG_TRIVIAL(info)<< format( "%g msec used for potential calculation (%g msec per atom)")
	% elapsed % (elapsed / info->atoms.size());

	_cf->releaseArrays();
	_t_d = af::array(_mc->ny, _mc->nx, _mc->nSlices,(afcfloat*)_t_d_ptr,afDevice);
	af::sync();
	if (_oc->SavePotential)
		SavePotential();
}

// don't forget to unlock the array once you are done with your kernel calls
cufftComplex* CUDA2DPotential::toCxPtr(af::array& a) {
	return (cufftComplex *) a.device<afcfloat>();
}
void CUDA2DPotential::copyToDeviceInt(int *devdst, std::vector<int> src, int size) {
	for (int i = 0; i < size; i++) {
		cuda_assert(cudaMemcpy(&devdst[i], &src[i], sizeof(int), cudaMemcpyHostToDevice));
	}
}

void CUDA2DPotential::copyToDeviceFloat(float_tt *devdst, std::vector<float_tt> src, int size) {
	for (int i = 0; i < size; i++) {
		cuda_assert(cudaMemcpy(&devdst[i], &src.data()[i], sizeof(float_tt), cudaMemcpyHostToDevice));
	}
}

af::array CUDA2DPotential::GetSlice(af::array t, unsigned idx) {
	return t(af::span, af::span, idx);
}

void CUDA2DPotential::AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ) {

}

void CUDA2DPotential::SaveAtomicPotential(int Z) {
	_atomPot[Z].resize(boost::extents[_mc->nx][_mc->ny]);
	cuda_assert(cudaMemcpy( _atomPot[Z].data(),_atomPot_d[Z], _mc->nx * _mc->ny* sizeof(cufftComplex), cudaMemcpyDeviceToHost));
	std::stringstream str;
	str << "atomicPotential_";
	str << Z;
	_persist->Save2DDataSet(_atomPot[Z], str.str());
}
void CUDA2DPotential::SavePotential() {
	_t.resize(boost::extents[_mc->nSlices][_mc->nx][_mc->ny]);
	_t_d.host(_t.data());
	CPotential::SavePotential();
}
void CUDA2DPotential::progressCounter(int j, int jTot) {
//	int interval = (jTot / 10);
//	if ((j % interval) == 0)
	loadbar(j, jTot);
}

} /* namespace QSTEM */

