/*
 * CUDA2DPotential.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */

#include "CUDA2DPotential.hpp"
#include <stdio.h>


namespace QSTEM {
CUDA2DPotential::CUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) :
		CPotential(mc, oc, p) {
}

CUDA2DPotential::~CUDA2DPotential() {
	delete _cf;
}
void CUDA2DPotential::initPotArrays() {
	_V_elem = af::array(slicePixels, c32);
	_V_accum = af::array(slicePixels, c32);
	_t_d = af::array(_mc->nx * _mc->ny * _mc->nSlices, c32);
	af::sync();
	_t_d_ptr = toCxPtr(_t_d);
	_V_elem_ptr = (cufftComplex *) _V_elem.device<af::af_cfloat>();
	_V_accum_ptr = (cufftComplex *) _V_accum.device<af::af_cfloat>();
	af::moddims(_V_elem, _mc->ny, _mc->nx);
	af::moddims(_V_accum, _mc->ny, _mc->nx);
}
void CUDA2DPotential::ComputeAtPot(superCellBoxPtr info){
	for (int Z : info->uniqueZ) {
		auto pot = af::array(slicePixels, c32);
		auto atomdelta = af::array(slicePixels, c32);
		af::sync();
		_cf->GetSincAtomicPotential(toCxPtr(pot), Z);
		_atomPot_d[Z] = pot;
		if (_oc->SaveAtomicPotential)
			SaveAtomicPotential(Z);
	}
}
void CUDA2DPotential::MakeSlices(superCellBoxPtr info) {

	_cf = new CUDAFunctions(info, _mc);
	initPotArrays();

	ComputeAtPot(info);

	af::timer time = af::timer::start();
	BOOST_LOG_TRIVIAL(info)<< "Calculating potential ...";
	for (int islice = 0; islice < _mc->nSlices; islice++) {
		cufftComplex* V_slice = &_t_d_ptr[islice * _mc->nx * _mc->ny];
		progressCounter(islice, _mc->nSlices);

		_cf->SetComplex( V_slice, 0.f, 0.f);
		_cf->SetComplex( _V_accum_ptr, 0.f, 0.f);
		for (int& Z : info->uniqueZ) {
			_cf->SetComplex( _V_elem_ptr, 0.f, 0.f);

			_cf->GetAtomDeltaFunctions(_V_elem_ptr, Z, islice);

			af::sync();

			_V_elem.unlock();
			_V_accum.unlock();

			_persist->SaveAtomDelta(_V_elem,islice,Z);

			af::fft2InPlace(_V_elem);
			_V_elem *= _atomPot_d[Z];
			af::ifft2InPlace(_V_elem);

			_V_accum = _V_elem + _V_accum;
			af::sync();

			_t_d.lock();
			_V_accum.lock();
			_cf->PotentialToTransmission(V_slice, _V_accum_ptr);
			_t_d.unlock();
			_V_accum.unlock();
		}
		af::sync();
	}

	auto elapsed = af::timer::stop(time) * 1000;
	BOOST_LOG_TRIVIAL(info)<< format( "%g msec used for potential calculation (%g msec per atom)")
	% elapsed % (elapsed / info->atoms.size());

	_t_d = af::moddims(_t_d, _mc->ny, _mc->nx, _mc->nSlices);
	if (_oc->SavePotential)
		SavePotential();
}

cufftComplex* CUDA2DPotential::toCxPtr(af::array a) {
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
	ComplexArray2D tmp(boost::extents[_mc->nx][_mc->ny]);
	_atomPot_d[Z].host(tmp.data());
	_atomPot[Z] = tmp;
	CPotential::SaveAtomicPotential(Z);
}
void CUDA2DPotential::SavePotential() {
	_t.resize(boost::extents[_t_d.dims(2)][_t_d.dims(1)][_t_d.dims(0)]);
	_t_d.host(_t.data());
	CPotential::SavePotential();
}
void CUDA2DPotential::progressCounter(int j, int jTot) {
	int interval = (jTot / 10);
	if ((j % interval) == 0)
		loadbar(j, jTot);
}

} /* namespace QSTEM */

