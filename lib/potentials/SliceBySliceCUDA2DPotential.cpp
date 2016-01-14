/*
 * CUDA2DPotential.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */

#include "SliceBySliceCUDA2DPotential.hpp"
#include <stdio.h>

namespace slicepp {
SliceBySliceCUDA2DPotential::SliceBySliceCUDA2DPotential(cModelConfPtr mc, cOutputConfPtr oc, cWaveConfPtr wc, PersistenceManagerPtr p) :
		CPotential(mc, oc, wc, p) {
}

SliceBySliceCUDA2DPotential::~SliceBySliceCUDA2DPotential() {
}

void SliceBySliceCUDA2DPotential::MakeSlices(superCellBoxPtr info) {
	copyDataToGPU(info);
	_info = info;
	_t_d.unlock();
}

void SliceBySliceCUDA2DPotential::copyDataToGPU(superCellBoxPtr info) {
	cf = new CUDAFunctions(info, _mc,_wc);

	slicePixels = _mc->n[0] * _mc->n[1];
	numAtoms = info->atoms.size();
	numAtUnique = info->uniqueZ.size();

	_t_d = af::array(_mc->n[0] * _mc->n[1] * _mc->n[2], c32);
	af::sync();
	tafPtr = _t_d.device<afcfloat>();
	potential = (cufftComplex *) tafPtr;

	af::timer time = af::timer::start();

	// create the needed arrays on the gpu, with a few workarounds to get arrayfire working with cuda kernels
	xyzPos = af::array(info->xyzPos.size(), (float_tt *) info->xyzPos.data(), afHost);
	xyzPos *= 1;
	xyzPos_d = xyzPos.device<float_tt>();
	occupancy = af::array(info->occupancy.size(), (float_tt *) info->occupancy.data(), afHost);
	occupancy *= 1;
	occupancy_d = occupancy.device<float_tt>();
	znums = af::array(info->znums.size(), (int *) info->znums.data(), afHost);
	znums *= 1;
	znums_d = znums.device<int>();
	af::sync();
	uniqueatoms = info->uniqueZ.data();

	BOOST_LOG_TRIVIAL(info)<< format( "%g msec used copying data to gpu") % (af::timer::stop(time)*1000);

	if (_oc->SavePotential || _oc->ComputeFromProjectedPotential) {
		_t_host.resize(boost::extents[_mc->n[2]][_mc->n[0]][_mc->n[1]]);
	}
}
void SliceBySliceCUDA2DPotential::copyToDeviceInt(int *devdst, std::vector<int> src, int size) {
	for (int i = 0; i < size; i++) {
		cuda_assert(cudaMemcpyAsync(&devdst[i], &src[i], sizeof(int), cudaMemcpyHostToDevice, afcu::getStream(af::getDevice())));
	}
}

void SliceBySliceCUDA2DPotential::copyToDeviceFloat(float_tt *devdst, std::vector<float_tt> src, int size) {
	for (int i = 0; i < size; i++) {
		cuda_assert(cudaMemcpyAsync(&devdst[i], &src.data()[i], sizeof(float_tt), cudaMemcpyHostToDevice, afcu::getStream(af::getDevice())));
	}
}

af::array SliceBySliceCUDA2DPotential::GetSlice(af::array t, unsigned idx) {

	af::array slicepot = af::array(_mc->n[0], _mc->n[1], c32);
	af::sync();
	tafPtr = slicepot.device<afcfloat>();
	potential = (cufftComplex *) tafPtr;

	cf->initPotArrays(slicePixels);
//	cf->GetPhaseGrating(&potential[idx * _mc->n[0] * _mc->n[1]], idx, _atomPot_d);
	cf->unlockArrays();

	af::sync();
	slicepot.unlock();

	if (_oc->SavePotential) {
		_t_d(af::seq(idx * slicePixels, (idx + 1) * slicePixels - 1)) = slicepot(af::span);
	}

	if (idx == _mc->n[2] - 1) {
		xyzPos.unlock();
		occupancy.unlock();
		znums.unlock();
		_t_d = af::moddims(_t_d, _mc->n[0], _mc->n[1], _mc->n[2]);
		SavePotential();
	}
	return slicepot;

}

void SliceBySliceCUDA2DPotential::AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ) {

}
void SliceBySliceCUDA2DPotential::SaveAtomicPotential(int znum) {

}
void SliceBySliceCUDA2DPotential::progressCounter(int j, int jTot) {
	int interval = (jTot / 10);
	if ((j % interval) == 0)
		loadbar(j, jTot);
}

} /* namespace slicepp */

