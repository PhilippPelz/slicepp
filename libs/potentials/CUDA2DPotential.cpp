/*
 * CUDA2DPotential.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */

#include "CUDA2DPotential.hpp"
#include <stdio.h>

namespace QSTEM {
CUDA2DPotential::CUDA2DPotential(const ConfigPtr& c,
		const PersistenceManagerPtr& persist) :
		CPotential(c, persist) {
	// TODO Auto-generated constructor stub

}

CUDA2DPotential::~CUDA2DPotential() {
	// TODO Auto-generated destructor stub
}

void CUDA2DPotential::MakeSlices(superCellBoxPtr info) {
	copyDataToGPU(info);
	if ((_c->Potential.CUDAOnTheFly == false) || (_c->ExperimentType == ExperimentType::PTYC)) {
		time_t time0, time1;
		cf->initPotArrays(slicePixels);
		BOOST_LOG_TRIVIAL(info)<< "Calculating potential ...";
		time(&time0);
		for (int islice = 0; islice < _c->Model.nSlices; islice++) {
			progressCounter(islice, _c->Model.nSlices);
			//cf->phaseGrating(cublasHandle, &potential[islice * slicePixels], numAtoms, numAtUnique, info->xyzPos, imPot, info->znums, info->uniqueatoms, info->occupancy, cufftPlanBatch, islice, _c->Model.nx, _c->Model.ny, _c->Model.nSlices, _c->Model.dx, _c->Model.dy, _c->Model.dz);
			cf->phaseGrating(&potential[islice * slicePixels], numAtoms,
					numAtUnique, xyzPos_d, imPot, znums_d, info->uniqueatoms.data(), occupancy_d, islice, _c->Model.nx, _c->Model.ny, _c->Model.nSlices, _c->Model.dx, _c->Model.dy, _c->Model.dz);
			_t_af.unlock();
			af::sync();
			tafPtr = _t_af.device<afcfloat>();
			potential = (cufftComplex *) tafPtr;
		}
		cf->unlockArrays();
		time(&time1);
		BOOST_LOG_TRIVIAL(info)<< format( "%g sec used for real space potential calculation (%g sec per atom)")
		% difftime(time1, time0)%( difftime(time1, time0) / info->atoms.size());
		_t_af.unlock();
		_t_af = af::moddims(_t_af, _c->Model.nx, _c->Model.ny, _c->Model.nSlices);

		SavePotential();

		xyzPos.unlock();
		occupancy.unlock();
		znums.unlock();
	}
	else{
		_t_af.unlock();
	}
}
void CUDA2DPotential::copyDataToGPU(superCellBoxPtr info) {
	time_t time0, time1;
	imPot = _c->Wave.imPot;
	cf = new CUDAFunctions();

	slicePixels = _c->Model.nx * _c->Model.ny;
	numAtoms = info->atoms.size();
	numAtUnique = info->uniqueatoms.size();

	_t_af = af::array(_c->Model.nx * _c->Model.ny * _c->Model.nSlices, c32);
	af::sync();
	tafPtr = _t_af.device<afcfloat>();
	potential = (cufftComplex *) tafPtr;

	time(&time0);
	xyzPos = af::array(info->xyzPos.size(), (float_tt *) info->xyzPos.data(), afHost);
	xyzPos *= 1;
	xyzPos_d = xyzPos.device<float_tt>();
	occupancy = af::array(info->occupancy.size(),(float_tt *) info->occupancy.data(), afHost);
	occupancy *= 1;
	occupancy_d = occupancy.device<float_tt>();
	znums = af::array(info->znums.size(), (int *) info->znums.data(), afHost);
	znums *= 1;
	znums_d = znums.device<int>();
	af::sync();
	uniqueatoms = info->uniqueatoms.data();
	time(&time1);
	BOOST_LOG_TRIVIAL(info)<< format( "%g sec used copying data to gpu")
	% difftime(time1, time0);
	if (_c->Output.SavePotential || _c->Output.ComputeFromProjectedPotential){
		_t.resize(boost::extents[_c->Model.nSlices][_c->Model.nx][_c->Model.ny]);
	}
}
void CUDA2DPotential::copyToDeviceInt(int *devdst, std::vector<int> src,
		int size) {
	for (int i = 0; i < size; i++) {
		cuda_assert(
				cudaMemcpy(&devdst[i], &src[i], sizeof(int),
						cudaMemcpyHostToDevice));
	}
}

void CUDA2DPotential::copyToDeviceFloat(float_tt *devdst,
		std::vector<float_tt> src, int size) {
	for (int i = 0; i < size; i++) {
		cuda_assert(
				cudaMemcpy(&devdst[i], &src.data()[i], sizeof(float_tt),
						cudaMemcpyHostToDevice));
	}
}

af::array CUDA2DPotential::GetSlice(af::array t, unsigned idx) {
	if ((_c->Potential.CUDAOnTheFly == false) || (_c->ExperimentType == ExperimentType::PTYC)) {
		return t(af::span, af::span, idx);
	}
	af::array slicepot = af::array(_c->Model.nx, _c->Model.ny, c32);
	af::sync();
	tafPtr = slicepot.device<afcfloat>();
	potential = (cufftComplex *) tafPtr;

	cf->initPotArrays(slicePixels);
	cf->phaseGrating(potential, numAtoms, numAtUnique, xyzPos_d, imPot, znums_d, uniqueatoms, occupancy_d, idx, _c->Model.nx,_c->Model.ny, _c->Model.nSlices, _c->Model.dx, _c->Model.dy,_c->Model.dz);
	cf->unlockArrays();

	af::sync();
	slicepot.unlock();
	if (_c->Output.SavePotential){
		_t_af(af::seq(idx * slicePixels, (idx + 1) * slicePixels - 1)) = slicepot(af::span);
	}
	if(idx == _c->Model.nSlices - 1){
		xyzPos.unlock();
		occupancy.unlock();
		znums.unlock();
		_t_af = af::moddims(_t_af, _c->Model.nx, _c->Model.ny, _c->Model.nSlices);
		SavePotential();
//		if (_c->Output.SavePotential){
//			_persist->SavePotential(_t);
//		}
//		if (_c->Output.SaveProjectedPotential) {
//			WriteProjectedPotential();
//		}
	}
	return slicepot;
}

void CUDA2DPotential::AddAtomToSlices(atom& atom, float_tt atomX,
		float_tt atomY, float_tt atomZ) {

}
void CUDA2DPotential::SaveAtomicPotential(int znum) {

}
void CUDA2DPotential::progressCounter(int j, int jTot) {
	int interval = (jTot / 10);
	if ((j % interval) == 0)
		loadbar(j, jTot);
}

} /* namespace QSTEM */

