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
#include <cufft.h>
#include "cublas_assert.hpp"
#include "cuda_assert.hpp"

namespace QSTEM {
CUDA2DPotential::CUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) :
		CPotential(mc, oc, p) {
	cufft_assert ( cufftPlan2d (_fftPlan, _mc->nx, _mc->ny, CUFFT_C2C ) );
	cublas_assert ( cublasCreate ( &_cublasHandle) );
}

CUDA2DPotential::~CUDA2DPotential() {
	delete _cf;
}
void CUDA2DPotential::initPotArrays() {
	_slicePixels = _mc->nx * _mc->ny;

    cudaMalloc((void**)&_t_d_ptr,  _mc->nSlices * _slicePixels*sizeof(afcfloat));
    cudaMalloc((void**)&_V_elem_ptr, _slicePixels*sizeof(afcfloat));
    cudaMalloc((void**)&_V_accum_ptr, _slicePixels*sizeof(afcfloat));

    _V_elem = af::array(_mc->nx, _mc->ny, (afcfloat*) _V_elem_ptr, afDevice);
	_V_accum = af::array(_mc->nx, _mc->ny,(afcfloat*) _V_accum_ptr, afDevice);

//	af::sync();
//	_t_d_ptr = toCxPtr(_t_d);
//	_V_elem_ptr = toCxPtr(_V_elem);
//	_V_accum_ptr = toCxPtr(_V_accum);
	_cf->SetComplex3D(_t_d_ptr, 0.f, 0.f);
	cuda_assert(cudaDeviceSynchronize());
	_t_d = af::array(_slicePixels * _mc->nSlices, (afcfloat*)_t_d_ptr, afDevice);
//	cuda_assert(cudaDeviceSynchronize());
}
void CUDA2DPotential::ComputeAtPot(superCellBoxPtr info) {
	for (int Z : info->uniqueZ) {
		afcfloat* _atomPot_d_ptr;
		cudaMalloc((void**)&_atomPot_d_ptr, _slicePixels*sizeof(afcfloat));
		_atomPot_d[Z] = af::array(_mc->ny * _mc->nx, _atomPot_d_ptr, afDevice);
//		_atomPot_d[Z] *= 1;
//		af::sync();
//		cufftComplex* pot = toCxPtr(_atomPot_d[Z]);
		BOOST_LOG_TRIVIAL(info)<< format("_atomPot_d_ptr: %d") % _atomPot_d_ptr;
		_cf->GetSincAtomicPotential((cufftComplex*)_atomPot_d_ptr, Z);
		cuda_assert(cudaDeviceSynchronize());
		_atomPot_d[Z].unlock();
//		af::sync();
//				cuda_assert(cudaDeviceSynchronize());
		_atomPot_d[Z] = af::moddims(_atomPot_d[Z], _mc->ny, _mc->nx);

		if (_oc->SaveAtomicPotential) {

			SaveAtomicPotential(Z);
		}
	}

}
void CUDA2DPotential::fft(cufftComplex* V){
	cufft_assert( cufftExecC2C( cufftPlanB, V, V, CUFFT_FORWARD ) );
}
void CUDA2DPotential::ifft(cufftComplex* V){
	cufft_assert( cufftExecC2C( cufftPlanB, V, V, CUFFT_INVERSE ) );
}
void CUDA2DPotential::XplusequY(cufftComplex* X, cufftComplex* Y){
	cublas_assert( cublasCaxpy( params->CU.cublasHandle, m12, &alpha, V1_d, 1, V_d, 1) );
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
		cuda_assert(cudaDeviceSynchronize());
		for (int& Z : info->uniqueZ) {
			_V_elem.lock();
			_cf->SetComplex2D(_V_elem_ptr, 0.f, 0.f);
			cuda_assert(cudaDeviceSynchronize());
			_cf->GetAtomDeltaFunctions(_V_elem_ptr, Z, islice);
			cuda_assert(cudaDeviceSynchronize());

			_V_elem.unlock();
			_V_accum.unlock();
//			if (_oc->SaveAtomDeltas) {
//				cuda_assert(cudaDeviceSynchronize());
//				float_tt rmin, rmax, aimin, aimax;
//				auto real = af::real(_V_elem);
//				auto imag = af::imag(_V_elem);
//				rmin = af::min<float_tt>(real);
//				rmax = af::max<float_tt>(real);
//				aimin = af::min<float_tt>(imag);
//				aimax = af::max<float_tt>(imag);
//				BOOST_LOG_TRIVIAL(info) << format("delta s=%-10d Z=%-3d (%-6.6g .. %-6.6g,i %-6.6g ... %-6.6g)") %islice%Z% rmin % rmax % aimin % aimax;
//				_persist->SaveAtomDelta(_V_elem, islice, Z);
//			}
//			cuda_assert(cudaDeviceSynchronize());


//			af::fft2InPlace(_V_elem);
//			_V_elem *= _atomPot_d[Z];
//			af::ifft2InPlace(_V_elem, _mc->nx * _mc->ny);

			fft(_V_elem_ptr);
			cf->cmul(_V_elem_ptr,_atomPot_d[Z]);
			ifft(_V_elem_ptr);

//			if (_oc->SaveAtomConv) {
//				cuda_assert(cudaDeviceSynchronize());
//				float_tt rmin, rmax, aimin, aimax;
//				auto real = af::real(_V_elem);
//				auto imag = af::imag(_V_elem);
//				rmin = af::min<float_tt>(real);
//				rmax = af::max<float_tt>(real);
//				aimin = af::min<float_tt>(imag);
//				aimax = af::max<float_tt>(imag);
//				BOOST_LOG_TRIVIAL(info)<< format("conv s=%-10d Z=%-3d (%-6.6g .. %-6.6g,i %-6.6g ... %-6.6g)") %islice%Z% rmin % rmax % aimin % aimax;
//				_persist->SaveAtomConv(_V_elem, islice, Z);
//			}

			_V_accum = _V_accum + _V_elem;
			af::sync();

//			_V_elem_ptr = toCxPtr(_V_elem);
//			BOOST_LOG_TRIVIAL(info)<< format("_V_elem_ptr: %d") % _V_elem_ptr;
		}

//		cuda_assert(cudaDeviceSynchronize());
		float_tt rmin, rmax, aimin, aimax;
		auto real = af::real(_V_accum);
		auto imag = af::imag(_V_accum);
		rmin = af::min<float_tt>(real);
		rmax = af::max<float_tt>(real);
		aimin = af::min<float_tt>(imag);
		aimax = af::max<float_tt>(imag);
		BOOST_LOG_TRIVIAL(info)<< format("accu s=%-16d (%-6.6g .. %-6.6g,i %-6.6g ... %-6.6g)") %islice% rmin % rmax % aimin % aimax;
//		_V_accum_ptr = toCxPtr(_V_accum);
		cufftComplex* V_slice_ptr = &_t_d_ptr[islice * _mc->nx * _mc->ny];
		_cf->PotentialToTransmission(V_slice_ptr, _V_accum_ptr);
		cuda_assert(cudaDeviceSynchronize());
		_V_accum.unlock();
		_V_accum_ptr = toCxPtr(_V_accum);
//		BOOST_LOG_TRIVIAL(info)<< format("_V_accum_ptr: %d") % _V_accum_ptr;
	}
	auto elapsed = af::timer::stop(time) * 1000;
	BOOST_LOG_TRIVIAL(info)<< format( "%g msec used for potential calculation (%g msec per atom)")
	% elapsed % (elapsed / info->atoms.size());

	_V_accum.unlock();
	_V_elem.unlock();

	_t_d = af::moddims(_t_d, _mc->ny, _mc->nx, _mc->nSlices);
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
	_atomPot_d[Z].host(_atomPot[Z].data());
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

