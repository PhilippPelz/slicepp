/*
 * PersistenceManager.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: philipp
 */

#include "PersistenceManager.hpp"

namespace QSTEM {
PersistenceManager::PersistenceManager() :
		_waveSlicesAfterFT(ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterTransmit(
				ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterProp(
				ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterSlice(
				ComplexArray3D(boost::extents[1][1][1])), _probe(
				ComplexArray2D(boost::extents[1][1])), _projectedPotential(
				ComplexArray2D(boost::extents[1][1])), _info(NULL), potSaved(false){
}
PersistenceManager::PersistenceManager(const ConfigPtr c) :
		PersistenceManager() {
	_c = c;
	_file = HDFFile(_c->Output.savePath);
	_info = _file.CreateGroup("config");
}

PersistenceManager::~PersistenceManager() {
	delete _info;
}
void PersistenceManager::SaveProbe(ComplexArray2DPtr a) {
	_probe = ComplexArray2D(a);
}

void PersistenceManager::SaveProbe(af::array wave) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);
	wave.host(a.data());
	SaveProbe(a);
}

void PersistenceManager::SaveWaveAfterTransmit(ComplexArray2DPtr a, int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements()/a.size(); j++) {
			_waveSlicesAfterTransmit[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterTransmit(af::array wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);
	wave.host(a.data());
	SaveWaveAfterTransmit(a, slice);
}
void PersistenceManager::SaveWaveAfterTransform(ComplexArray2DPtr a,
		int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements()/a.size(); j++) {
			_waveSlicesAfterFT[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterTransform(af::array wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);
	wave.host(a.data());
	SaveWaveAfterTransform(a, slice);
}
void PersistenceManager::SaveWaveAfterPropagation(ComplexArray2DPtr a,
		int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements()/a.size(); j++) {
			_waveSlicesAfterProp[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterPropagation(af::array wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);

	wave.host(a.data());

	SaveWaveAfterPropagation(a, slice);
}
void PersistenceManager::SaveWaveAfterSlice(ComplexArray2DPtr a, int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements()/a.size(); j++) {
			_waveSlicesAfterSlice[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterSlice(af::array wave, int slice) {
	ComplexArray2D a(boost::extents[wave.dims(0)][wave.dims(1)]);

	wave.host(a.data());

	SaveWaveAfterSlice(a, slice);
}
void PersistenceManager::SavePotential(ComplexArray3D a) {
	_potential = a;
}

void PersistenceManager::SavePotential(std::vector<af::array> data, int nslices) {
	potSaved = true;
	ComplexArray3D a(boost::extents[nslices][data[0].dims(0)][data[0].dims(1)]);
	for (int i = 0; i < nslices; i++){
		data[i].host(a[i].origin());
	}
	SavePotential(a);
}

void PersistenceManager::SaveProjectedPotential(ComplexArray2DPtr a) {
	_projectedPotential = ComplexArray2D(a);
}

void PersistenceManager::SaveProjectedPotential(af::array data) {
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

void PersistenceManager::Save2DDataSet(af::array data, string name) {
	ComplexArray2D a(boost::extents[data.dims(0)][data.dims(1)]);
	data.host(a.data());
	Save2DDataSet(a, name);
}

void PersistenceManager::InitStorage() {
	auto e3 = boost::extents[_c->Model.nSlices][_c->Model.nx][_c->Model.ny];
	auto e2 = boost::extents[_c->Model.nx][_c->Model.ny];
	_potential.resize(e3);
	_waveSlicesAfterSlice.resize(e3);
	_waveSlicesAfterProp.resize(e3);
	_waveSlicesAfterFT.resize(e3);
	_waveSlicesAfterTransmit.resize(e3);
	_probe.resize(e2);
	_projectedPotential.resize(e2);
}

void PersistenceManager::ResizeStorage(int xdim, int ydim) {
	auto e3 = boost::extents[_c->Model.nSlices][xdim][ydim];
	auto e2 = boost::extents[xdim][ydim];
	_potential.resize(e3);
	_waveSlicesAfterSlice.resize(e3);
	_waveSlicesAfterProp.resize(e3);
	_waveSlicesAfterFT.resize(e3);
	_waveSlicesAfterTransmit.resize(e3);
	_probe.resize(e2);
	_projectedPotential.resize(e2);
}



void PersistenceManager::StoreToDisc() {
	_file.SaveComplexArray2D(_probe, "probe");
	_file.SaveComplexArray2D(_projectedPotential, "projectedPotential");
	_file.SaveComplexArray3D(_waveSlicesAfterTransmit,"waveSlicesAfterTransmit");
	_file.SaveComplexArray3D(_waveSlicesAfterFT, "waveSlicesAfterFT");
	_file.SaveComplexArray3D(_waveSlicesAfterProp, "waveSlicesAfterProp");
	_file.SaveComplexArray3D(_waveSlicesAfterSlice, "waveSlicesAfterSlice");
	_file.SaveComplexArray3D(_potential, "potentialSlices");
}
} /* namespace QSTEM */

