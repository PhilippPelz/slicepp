/*
 * CUDA2DPotential.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */

#include "CUDA2DPotential.hpp"
#include <stdio.h>

namespace QSTEM {
CUDA2DPotential::CUDA2DPotential(const ConfigPtr& c, const PersistenceManagerPtr& persist) :
		CPotential(c, persist) {
}

CUDA2DPotential::~CUDA2DPotential() {
}

void CUDA2DPotential::MakeSlices(superCellBoxPtr info) {
	copyDataToGPU(info);
	if ((_c->Potential->CUDAOnTheFly == false) || (_c->ExperimentType == ExperimentType::PTYC)) {

		af::timer time = af::timer::start();
		cf->initPotArrays(slicePixels);
		BOOST_LOG_TRIVIAL(info)<< "Calculating potential ...";
		for (int islice = 0; islice < _c->Model->nSlices; islice++) {

			progressCounter(islice, _c->Model->nSlices);
			cf->phaseGrating(&potential[islice * slicePixels], numAtoms, numAtUnique, xyzPos_d, imPot, znums_d, info->uniqueatoms.data(), occupancy_d,
					islice, _c->Model->nx, _c->Model->ny, _c->Model->nSlices, _c->Model->dx, _c->Model->dy, _c->Model->dz);
			_t_device.unlock();
			af::sync();
			tafPtr = _t_device.device<afcfloat>();
			potential = (cufftComplex *) tafPtr;
		}
		cf->unlockArrays();

		BOOST_LOG_TRIVIAL(info)<< format( "%g msec used for real space potential calculation (%g msec per atom)")
		% (af::timer::stop(time)*1000) % ( (af::timer::stop(time)*1000) / info->atoms.size());
		_t_device.unlock();
		_t_device = af::moddims(_t_device, _c->Model->nx, _c->Model->ny, _c->Model->nSlices);

		SavePotential();

		xyzPos.unlock();
		occupancy.unlock();
		znums.unlock();
	} else {
		_t_device.unlock();
	}
}
void CUDA2DPotential::copyDataToGPU(superCellBoxPtr info) {
	imPot = _c->Wave->imPot;
	cf = new CUDAFunctions();

	slicePixels = _c->Model->nx * _c->Model->ny;
	numAtoms = info->atoms.size();
	numAtUnique = info->uniqueatoms.size();

	_t_device = af::array(_c->Model->nx * _c->Model->ny * _c->Model->nSlices, c32);
	af::sync();
	tafPtr = _t_device.device<afcfloat>();
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
	uniqueatoms = info->uniqueatoms.data();

	BOOST_LOG_TRIVIAL(info)<< format( "%g msec used copying data to gpu") % (af::timer::stop(time)*1000);

	if (_c->Output->SavePotential || _c->Output->ComputeFromProjectedPotential) {
		_t_host.resize(boost::extents[_c->Model->nSlices][_c->Model->nx][_c->Model->ny]);
	}
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
	if ((_c->Potential->CUDAOnTheFly == false) || (_c->ExperimentType == ExperimentType::PTYC)) {

		return t(af::span, af::span, idx);

	} else {

		af::array slicepot = af::array(_c->Model->nx, _c->Model->ny, c32);
		af::sync();
		tafPtr = slicepot.device<afcfloat>();
		potential = (cufftComplex *) tafPtr;

		cf->initPotArrays(slicePixels);
		cf->phaseGrating(potential, numAtoms, numAtUnique, xyzPos_d, imPot, znums_d, uniqueatoms, occupancy_d, idx, _c->Model->nx, _c->Model->ny,
				_c->Model->nSlices, _c->Model->dx, _c->Model->dy, _c->Model->dz);
		cf->unlockArrays();

		af::sync();
		slicepot.unlock();

		if (_c->Output->SavePotential) {
			_t_device(af::seq(idx * slicePixels, (idx + 1) * slicePixels - 1)) = slicepot(af::span);
		}

		if (idx == _c->Model->nSlices - 1) {
			xyzPos.unlock();
			occupancy.unlock();
			znums.unlock();
			_t_device = af::moddims(_t_device, _c->Model->nx, _c->Model->ny, _c->Model->nSlices);
			SavePotential();
		}
		return slicepot;

	}
}

void CUDA2DPotential::AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ) {

}
void CUDA2DPotential::SaveAtomicPotential(int znum) {

}
void CUDA2DPotential::progressCounter(int j, int jTot) {
	int interval = (jTot / 10);
	if ((j % interval) == 0)
		loadbar(j, jTot);
}

} /* namespace QSTEM */

