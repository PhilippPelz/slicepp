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

namespace slicepp {
CUDA2DPotential::CUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) :
		CPotential(mc, oc, p) {
	cufft_assert(cufftPlan2d(&_fftPlan, _mc->n[0], _mc->n[1], CUFFT_C2C));
	cufft_assert(cufftSetStream(_fftPlan, afcu::getStream(af::getDevice())));
	cublas_assert(cublasCreate ( &_cublasHandle));
	cublas_assert(cublasSetStream(_cublasHandle, afcu::getStream(af::getDevice())));
}

CUDA2DPotential::~CUDA2DPotential() {
	cublas_assert(cublasDestroy(_cublasHandle));
	cufft_assert(cufftDestroy(_fftPlan));
	for (auto& kv : _atomPot_d) {
		cuda_assert(cudaFree(kv.second));
	}
}

void CUDA2DPotential::printAfMem() {
	size_t alloc_bytes, alloc_buffers, lock_bytes, lock_buffers;

	af::deviceMemInfo(&alloc_bytes, &alloc_buffers, &lock_bytes, &lock_buffers);
	BOOST_LOG_TRIVIAL(info)<< format("GPU memory usage: allocated = %f MB") % ((double) (alloc_bytes) / 1024.0 / 1024.0);
}
void CUDA2DPotential::printCudaMem() {
	size_t free_byte;
	size_t total_byte;

	cuda_assert(cudaMemGetInfo(&free_byte, &total_byte));
	printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", (double) (total_byte - free_byte) / 1024.0 / 1024.0,
			(double) free_byte / 1024.0 / 1024.0, (double) total_byte / 1024.0 / 1024.0);
}

void CUDA2DPotential::initPotArrays(superCellBoxPtr info) {
	_slicePixels = _mc->n[0] * _mc->n[1];
	xyzPos = af::array(info->xyzPos.size(), (float_tt*) info->xyzPos.data(), afHost);
	occupancy = af::array(info->occupancy.size(), (float_tt*) info->occupancy.data(), afHost);
	znums = af::array(info->znums.size(), (int*) info->znums.data(), afHost);
	t_d = af::array(_mc->n[2] * _slicePixels, c32);
	V_elem = af::array(_slicePixels, c32);
	V_accum = af::array(_slicePixels, c32);

	_t_d_ptr = (float2*) t_d.device<af::cfloat>();
	_V_elem_ptr = (float2*) V_elem.device<af::cfloat>();
	_V_accum_ptr = (float2*) V_accum.device<af::cfloat>();
	printAfMem();
}
void CUDA2DPotential::ComputeAtomicPotential(superCellBoxPtr info) {
	for (int Z : info->uniqueZ) {
		cufftComplex* _atomPot_d_ptr;
		cudaMalloc((void**) &_atomPot_d_ptr, _slicePixels * sizeof(cufftComplex));
		_atomPot_d[Z] = _atomPot_d_ptr;
//		auto p = af::array(_slicePixels, c32);
//		_atomPot_d[Z] = (cufftComplex*) p.device<af::cfloat>();
//		printf("_atomPot_d[%d]:  %#08x\n", Z, _atomPot_d[Z]);
		_cf.GetSincAtomicPotential(_atomPot_d[Z], Z);
		if (_oc->SaveAtomicPotential) {
			sync();
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
void CUDA2DPotential::sync() {
	cuda_assert(cudaStreamSynchronize(afcu::getStream(af::getDevice())));
}
void CUDA2DPotential::MakeSlices(superCellBoxPtr info) {

	_cf = CUDAFunctions(info, _mc);
	initPotArrays(info);
	_cf.initArrays();
	ComputeAtomicPotential(info);
	af::timer time = af::timer::start();
	BOOST_LOG_TRIVIAL(info)<< "Calculating potential ...";
	xyzPos_d = xyzPos.device<float>();
	occupancy_d = occupancy.device<float>();
	znums_d = znums.device<int>();
//	printf("znums_d:  %#08x\n", znums_d);
//	printf("occupancy_d:  %#08x\n", occupancy_d);
//	printf("xyzPos_d:  %#08x\n", xyzPos_d);
	for (int islice = 0; islice < _mc->n[2]; islice++) {
		progressCounter(islice, _mc->n[2]);
		_cf.SetComplex2D(_V_accum_ptr, 0.f, 0.f);
		for (int& Z : info->uniqueZ) {
			_cf.SetComplex2D(_V_elem_ptr, 0.f, 0.f);
			_cf.GetAtomDeltaFunctions(_V_elem_ptr, Z, islice, xyzPos_d, occupancy_d, znums_d);
			if (_oc->SaveAtomDeltas) {
				_persist->SaveAtomDelta(_V_elem_ptr, islice, Z);
			}
			fft(_V_elem_ptr);
			_cf.cmul(_V_elem_ptr, _atomPot_d[Z]);
			ifft(_V_elem_ptr);

			if (_oc->SaveAtomConv) {
				_persist->SaveAtomConv(_V_elem_ptr, islice, Z);
			}

			XplusequY(_V_accum_ptr, _V_elem_ptr);
		}
		fft(_V_accum_ptr);
		_cf.limitBandwidth(_V_accum_ptr);
		ifft(_V_accum_ptr);
		_cf.PotentialToTransmission(&_t_d_ptr[islice * _slicePixels], _V_accum_ptr);
		sync();
	}
	auto elapsed = af::timer::stop(time) * 1000;
	BOOST_LOG_TRIVIAL(info)<< format("%g msec used for potential calculation (%g msec per atom)") % elapsed % (elapsed / info->atoms.size());
	_cf.releaseArrays();
	_t_d = af::array(_mc->n[1], _mc->n[0], _mc->n[2], (afcfloat*)_t_d_ptr,afDevice);

	if (_oc->SavePotential) {
		af::sync();

		SavePotential();
	}
}

af::array CUDA2DPotential::GetSlice(af::array t, unsigned idx) {
	return t(af::span, af::span, idx);
}

void CUDA2DPotential::AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ) {

}

void CUDA2DPotential::SaveAtomicPotential(int Z) {
	_atomPot[Z].resize(boost::extents[_mc->n[0]][_mc->n[1]]);
	cuda_assert(
			cudaMemcpyAsync(_atomPot[Z].data(), _atomPot_d[Z], _mc->n[0] * _mc->n[1] * sizeof(cufftComplex), cudaMemcpyDeviceToHost,
					afcu::getStream(af::getDevice())));
	std::stringstream str;
	str << "atomicPotential_";
	str << Z;
	_persist->SaveCx2DDataSet(_atomPot[Z], str.str());
}
void CUDA2DPotential::SavePotential() {
	_t.resize(boost::extents[_mc->n[2]][_mc->n[0]][_mc->n[1]]);
	_t_d.host(_t.data());
	_persist->SavePotential(_t);
	if (_oc->SaveProjectedPotential) {
		auto projectedPotential = af::sum(_t_d, 2);
		_persist->SaveProjectedPotential(projectedPotential);
	}
}
void CUDA2DPotential::progressCounter(int j, int jTot) {
	int interval = (jTot / 10);
	if ((j % interval) == 0)
		loadbar(j, jTot);
}

} /* namespace slicepp */

