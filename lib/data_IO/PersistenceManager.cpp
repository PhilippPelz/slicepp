/*
 * PersistenceManager.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: philipp
 */

#include "PersistenceManager.hpp"
#include <string>
#include "boost/format.hpp"

using boost::format;

namespace QSTEM {
PersistenceManager::PersistenceManager() :
		_waveSlicesAfterFT(ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterTransmit(ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterProp(
				ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterSlice(ComplexArray3D(boost::extents[1][1][1])), _probe(
				ComplexArray2D(boost::extents[1][1])), _projectedPotential(ComplexArray2D(boost::extents[1][1])), _info(NULL), _potSaved(false) {
}
PersistenceManager::PersistenceManager(const ConfigPtr c) :
		PersistenceManager() {
	_c = c;
	_file = HDFFile(_c->Output->savePath);
	_info = _file.CreateGroup("config");
}

PersistenceManager::~PersistenceManager() {
	delete _info;
}

void PersistenceManager::SaveProbe(ComplexArray2DPtr a) {
	_probe = ComplexArray2D(a);
}

void PersistenceManager::SaveProbe(af::array& wave) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);
	wave.host(a.data());
	SaveProbe(a);
}
void PersistenceManager::SaveAtomDelta(cuComplex* delta, int slice, int Z) {
	auto a = _atomDeltas[Z];
	if(a.get() == NULL){
		auto e3 = boost::extents[_c->Model->nSlices][_c->Model->nx][_c->Model->ny];
		_atomDeltas[Z] = boost::shared_ptr<ComplexArray3D>(new ComplexArray3D(e3));
		auto sh = _atomDeltas[Z]->shape();
//		BOOST_LOG_TRIVIAL(info)<< format("%s shape: [%d,%d,%d]") % "_atomDeltas" % sh[0] % sh[1] % sh[2];
	}
	cuda_assert(cudaMemcpy( _atomDeltas[Z]->operator [](slice).origin(),delta, _c->Model->nx * _c->Model->ny* sizeof(cuComplex), cudaMemcpyDeviceToHost));
}
void PersistenceManager::SaveAtomConv(cuComplex* delta, int slice, int Z) {
	auto a = _atomConv[Z];
	if(a.get() == NULL){
		auto e3 = boost::extents[_c->Model->nSlices][_c->Model->nx][_c->Model->ny];
		_atomConv[Z] = boost::shared_ptr<ComplexArray3D>(new ComplexArray3D(e3));
		auto sh = _atomConv[Z]->shape();
	}
	cuda_assert(cudaMemcpy( _atomConv[Z]->operator [](slice).origin(),delta, _c->Model->nx * _c->Model->ny* sizeof(cuComplex), cudaMemcpyDeviceToHost));
}
void PersistenceManager::SaveWaveAfterTransmit(ComplexArray2DPtr a, int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterTransmit[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterTransmit(af::array& wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);
	wave.host(a.data());
	SaveWaveAfterTransmit(a, slice);
}
void PersistenceManager::SaveWaveAfterTransform(ComplexArray2DPtr a, int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterFT[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterTransform(af::array& wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);
	wave.host(a.data());
	SaveWaveAfterTransform(a, slice);
}
void PersistenceManager::SaveWaveAfterPropagation(ComplexArray2DPtr a, int slice) {
#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterProp[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterPropagation(af::array& wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);

	wave.host(a.data());

	SaveWaveAfterPropagation(a, slice);
}
void PersistenceManager::SaveWaveAfterSlice(ComplexArray2DPtr a, int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterSlice[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterSlice(af::array& wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);

	wave.host(a.data());

	SaveWaveAfterSlice(a, slice);
}
void PersistenceManager::SavePotential(ComplexArray3D a) {
	_potSaved = true;
	_potential = a;
}

void PersistenceManager::SavePotential(af::array& data) {
	ComplexArray3D a(boost::extents[data.dims(2)][data.dims(0)][data.dims(1)]);
	data.host(a.data());
	SavePotential(a);
}

void PersistenceManager::SaveProjectedPotential(ComplexArray2DPtr a) {
	_projectedPotential = ComplexArray2D(a);
}

void PersistenceManager::SaveProjectedPotential(af::array& data) {
	ComplexArray2D a(boost::extents[data.dims(0)][data.dims(1)]);
	data.host(a.data());
	SaveProjectedPotential(a);
}

void PersistenceManager::Save3DDataSet(ComplexArray3DPtr a, string name) {
	_file.SaveComplexArray3D(ComplexArray3D(a), name);
}
void PersistenceManager::Save2DDataSet(ComplexArray2DPtr a, string name) {
	_file.SaveComplexArray2D(ComplexArray2D(a), name);
}

void PersistenceManager::Save2DDataSet(af::array& data, string name) {
	ComplexArray2D a(boost::extents[data.dims(0)][data.dims(1)]);
	data.host(a.data());
	Save2DDataSet(a, name);
}

void PersistenceManager::InitStorage() {
	auto e3 = boost::extents[_c->Model->nSlices][_c->Model->nx][_c->Model->ny];
	auto e2 = boost::extents[_c->Model->nx][_c->Model->ny];
	if (_c->Output->SavePotential || _c->Output->ComputeFromProjectedPotential)
		_potential.resize(e3);
	if (_c->Output->SaveWaveAfterSlice)
		_waveSlicesAfterSlice.resize(e3);
	if (_c->Output->SaveWaveAfterPropagation)
		_waveSlicesAfterProp.resize(e3);
	if (_c->Output->SaveWaveAfterTransform)
		_waveSlicesAfterFT.resize(e3);
	if (_c->Output->SaveWaveAfterTransmit)
		_waveSlicesAfterTransmit.resize(e3);
	if (_c->Output->saveProbe)
		_probe.resize(e2);
	if (_c->Output->SaveAtomDeltas)
		for (auto& kv : _atomDeltas)
			kv.second->resize(e3);
	if (_c->Output->SaveAtomConv)
		for (auto& kv : _atomConv)
			kv.second->resize(e3);
	if (_c->Output->SaveProjectedPotential || _c->Output->ComputeFromProjectedPotential)
		_projectedPotential.resize(e2);
}

void PersistenceManager::ResizeStorage(int xdim, int ydim) {
	auto e3 = boost::extents[_c->Model->nSlices][xdim][ydim];
	auto e2 = boost::extents[xdim][ydim];
	if (_c->Output->SavePotential || _c->Output->ComputeFromProjectedPotential)
		_potential.resize(boost::extents[_c->Model->nSlices][_c->Model->nx][_c->Model->ny]);
	if (_c->Output->SaveWaveAfterSlice)
		_waveSlicesAfterSlice.resize(e3);
	if (_c->Output->SaveWaveAfterPropagation)
		_waveSlicesAfterProp.resize(e3);
	if (_c->Output->SaveWaveAfterTransform)
		_waveSlicesAfterFT.resize(e3);
	if (_c->Output->SaveWaveAfterTransmit)
		_waveSlicesAfterTransmit.resize(e3);
	if (_c->Output->saveProbe)
		_probe.resize(e2);
	if (_c->Output->SaveAtomDeltas)
		for (auto& kv : _atomDeltas)
			kv.second->resize(e3);
	if (_c->Output->SaveAtomConv)
		for (auto& kv : _atomConv)
			kv.second->resize(e3);
	if (_c->Output->SaveProjectedPotential || _c->Output->ComputeFromProjectedPotential)
		_projectedPotential.resize(boost::extents[_c->Model->nx][_c->Model->ny]);
}

void PersistenceManager::StoreToDisc() {
	_file.SaveComplexArray2D(_probe, "probe");
	_file.SaveComplexArray2D(_projectedPotential, "projectedPotential");
	_file.SaveComplexArray3D(_waveSlicesAfterTransmit, "waveSlicesAfterTransmit");
	_file.SaveComplexArray3D(_waveSlicesAfterFT, "waveSlicesAfterFT");
	_file.SaveComplexArray3D(_waveSlicesAfterProp, "waveSlicesAfterProp");
	_file.SaveComplexArray3D(_waveSlicesAfterSlice, "waveSlicesAfterSlice");
	_file.SaveComplexArray3D(_potential, "potentialSlices");
	for (auto& kv : _atomDeltas) {
		string s = "atomDeltas_";
		s += std::to_string(kv.first);
		auto arr = kv.second;
		auto sh = arr->shape();
//		BOOST_LOG_TRIVIAL(info)<< format("%s shape: [%d,%d,%d]") % s % sh[0] % sh[1] % sh[2];
		_file.SaveComplexArray3D(*arr.get(), s);
	}
	for (auto& kv : _atomConv) {
		string s = "atomConv_";
		s += std::to_string(kv.first);
		auto arr = kv.second;
		auto sh = arr->shape();
//		BOOST_LOG_TRIVIAL(info)<< format("%s shape: [%d,%d,%d]") % s % sh[0] % sh[1] % sh[2];
		_file.SaveComplexArray3D(*arr.get(), s);
	}
}


void PersistenceManager::StoreToDiscMP(int pos, int x, int y) {
std::string info = "_" + std::to_string(pos) + "_(" + to_string(x) + ", " + to_string(y) + ")";
if (pos == 1) {
	_file.SaveComplexArray2D(_probe, "probe");
	_file.SaveComplexArray2D(_projectedPotential, "projectedPotential");
	_file.SaveComplexArray3D(_waveSlicesAfterTransmit, "waveSlicesAfterTransmit" + info);
	_file.SaveComplexArray3D(_waveSlicesAfterFT, "waveSlicesAfterFT" + info);
	_file.SaveComplexArray3D(_waveSlicesAfterProp, "waveSlicesAfterProp" + info);
	_file.SaveComplexArray3D(_waveSlicesAfterSlice, "waveSlicesAfterSlice" + info);
	_file.SaveComplexArray3D(_potential, "potentialSlices");
} else {
//		_file.SaveComplexArray3D(_waveSlicesAfterTransmit,"waveSlicesAfterTransmit" + info);
//		_file.SaveComplexArray3D(_waveSlicesAfterFT, "waveSlicesAfterFT");
//		_file.SaveComplexArray3D(_waveSlicesAfterProp, "waveSlicesAfterProp" + info);
	_file.SaveComplexArray3D(_waveSlicesAfterSlice, "waveSlicesAfterSlice" + info);
}
}

} /* namespace QSTEM */

