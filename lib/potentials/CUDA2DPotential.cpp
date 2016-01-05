/*
 * CUDA2DPotential.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */

#include "CUDA2DPotential.hpp"
#include <stdio.h>
#include <thread>
#include <chrono>
//#include <thrust/fill.h>
//#include <thrust/copy.h>

namespace slicepp {
CUDA2DPotential::CUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) :
		CPotential(mc, oc, p) {
	cufft_assert(cufftPlan2d(&_fftPlan, _mc->n[0], _mc->n[1], CUFFT_C2C));
	cublas_assert(cublasCreate ( &_cublasHandle));
}

CUDA2DPotential::~CUDA2DPotential() {
	delete _cf;
	cublas_assert(cublasDestroy(_cublasHandle));
	cufft_assert(cufftDestroy(_fftPlan));
	for (auto& kv : _atomPot_d) {
		cuda_assert(cudaFree(kv.second));
	}
	cuda_assert(cudaFree(_t_d_ptr));
	cuda_assert(cudaFree(_V_elem_ptr));
	cuda_assert(cudaFree(_V_accum_ptr));
}

void CUDA2DPotential::printAfMem(){
	size_t alloc_bytes, alloc_buffers, lock_bytes, lock_buffers;

	af::deviceMemInfo (&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
	printf("AF GPU memory usage: alloc = %f, alloc_buffers = %f MB, lock = %f MB, lock_buffers = %f MB\n", (double) (alloc_bytes) / 1024.0 / 1024.0,
			(double) alloc_buffers / 1024.0 / 1024.0, (double) lock_bytes / 1024.0 / 1024.0, (double) lock_buffers / 1024.0 / 1024.0);
}
void CUDA2DPotential::printCudaMem(){
	size_t free_byte;
	size_t total_byte;

	cuda_assert(cudaMemGetInfo(&free_byte, &total_byte));
	printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", (double) (total_byte - free_byte) / 1024.0 / 1024.0,
			(double) free_byte / 1024.0 / 1024.0, (double) total_byte / 1024.0 / 1024.0);
}

void CUDA2DPotential::initPotArrays() {

	printAfMem();

	//	cuda_assert(cudaMalloc((void**)&_V_elem_ptr, _slicePixels * sizeof(cufftComplex)));
	//	cuda_assert(cudaMalloc((void**)&_t_d_ptr, _mc->n[2] * _slicePixels * sizeof(cufftComplex)));
	//	cuda_assert(cudaMalloc((void**)&_V_accum_ptr, _slicePixels * sizeof(cufftComplex)));

	_slicePixels = _mc->n[0] * _mc->n[1];
//	printf("pot elem: %d\n",_mc->n[2] * _slicePixels);

	auto t_d = af::array(_mc->n[2] * _slicePixels,c32);
	printAfMem();
	auto V_elem = af::array(_slicePixels,c32);
	printAfMem();
	auto V_accum = af::array(_slicePixels,c32);
	printAfMem();

	_t_d_ptr = (float2*)t_d.device<af::cfloat>();
	_V_elem_ptr = (float2*)V_elem.device<af::cfloat>();
	_V_accum_ptr = (float2*)V_accum.device<af::cfloat>();

	printf("_t_d_ptr:  %#08x\n", _t_d_ptr);
	printf("_V_elem_ptr:  %#08x\n", _V_elem_ptr);
	printf("_V_accum_ptr:  %#08x\n", _V_accum_ptr);

	cuda_assert(cudaDeviceSynchronize());
	std::this_thread::sleep_for(std::chrono::seconds(1));
//	_cf->SetComplex3D(_t_d_ptr, 0.f, 0.f);

}
void CUDA2DPotential::ComputeAtomicPotential(superCellBoxPtr info) {
	for (int Z : info->uniqueZ) {
		cufftComplex* _atomPot_d_ptr;
		cudaMalloc((void**) &_atomPot_d_ptr, _slicePixels * sizeof(cufftComplex));
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
	ComputeAtomicPotential(info);
	af::timer time = af::timer::start();
	BOOST_LOG_TRIVIAL(info) << "Calculating potential ...";

	for (int islice = 0; islice < _mc->n[2]; islice++) {
		progressCounter(islice, _mc->n[2]);
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
	BOOST_LOG_TRIVIAL(info) << format("%g msec used for potential calculation (%g msec per atom)") % elapsed % (elapsed / info->atoms.size());
	_cf->releaseArrays();
	_t_d = af::array(_mc->n[1], _mc->n[0], _mc->n[2], (afcfloat*)_t_d_ptr,afDevice);
	af::sync();
	if (_oc->SavePotential)
		SavePotential();
}

af::array CUDA2DPotential::GetSlice(af::array t, unsigned idx) {
	return t(af::span, af::span, idx);
}

void CUDA2DPotential::AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ) {

}

void CUDA2DPotential::SaveAtomicPotential(int Z) {
	_atomPot[Z].resize(boost::extents[_mc->n[0]][_mc->n[1]]);
	cuda_assert(cudaMemcpy(_atomPot[Z].data(), _atomPot_d[Z], _mc->n[0] * _mc->n[1] * sizeof(cufftComplex), cudaMemcpyDeviceToHost));
	std::stringstream str;
	str << "atomicPotential_";
	str << Z;
	_persist->SaveCx2DDataSet(_atomPot[Z], str.str());
}
void CUDA2DPotential::SavePotential() {
	_t.resize(boost::extents[_mc->n[2]][_mc->n[0]][_mc->n[1]]);
	_t_d.host(_t.data());
	CPotential::SavePotential();
}
void CUDA2DPotential::progressCounter(int j, int jTot) {
//	int interval = (jTot / 10);
//	if ((j % interval) == 0)
	loadbar(j, jTot);
}

} /* namespace slicepp */

