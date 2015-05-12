/*
 * PersistenceManager.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: philipp
 */

#include "PersistenceManager.hpp"

namespace QSTEM {
PersistenceManager::PersistenceManager() {
}
PersistenceManager::PersistenceManager(const ConfigPtr c) : PersistenceManager() {
	_c = c;
	_file = HDFFile(_c->Output.savePath);
	_info = _file.CreateGroup("config");
	_waveSlicesAfterSlice.resize(extents[_c->Model.nSlices][_c->Model.nx][_c->Model.ny]);
	_waveSlicesAfterTransmit.resize(extents[_c->Model.nSlices][_c->Model.nx][_c->Model.ny]);
	_waveSlicesAfterFT.resize(extents[_c->Model.nSlices][_c->Model.nx][_c->Model.ny]);
	_waveSlicesAfterProp.resize(extents[_c->Model.nSlices][_c->Model.nx][_c->Model.ny]);
	_potential.resize(extents[_c->Model.nSlices][_c->Model.nx][_c->Model.ny]);
}

PersistenceManager::~PersistenceManager() {
	delete _info;
}
void PersistenceManager::SaveProbe(ComplexArray2DPtr a){
	_probe.resize(extents[_c->Model.nx][_c->Model.ny]);
	_probe = ComplexArray2D(a);
}
void PersistenceManager::SaveWaveAfterTransmit(ComplexArray2DPtr a,int slice){

#pragma omp parallel for
	for(int i=0;i<_c->Model.nx;i++)
		for(int j=0;j<_c->Model.ny;j++){
			_waveSlicesAfterTransmit[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterTransform(ComplexArray2DPtr a,int slice){

#pragma omp parallel for
	for(int i=0;i<_c->Model.nx;i++)
		for(int j=0;j<_c->Model.ny;j++){
			_waveSlicesAfterFT[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterPropagation(ComplexArray2DPtr a,int slice){

#pragma omp parallel for
	for(int i=0;i<_c->Model.nx;i++)
		for(int j=0;j<_c->Model.ny;j++){
			_waveSlicesAfterProp[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterSlice(ComplexArray2DPtr a,int slice){

#pragma omp parallel for
	for(int i=0;i<_c->Model.nx;i++)
		for(int j=0;j<_c->Model.ny;j++){
			_waveSlicesAfterSlice[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SavePotential(ComplexArray3D a){

//	BOOST_LOG_TRIVIAL(info) << a.shape()[0] << " " << _potential.shape()[0]<< a.shape()[1] << _potential.shape()[1]<< a.shape()[2] << _potential.shape()[2];
	_potential = ComplexArray3D(a);
}
void PersistenceManager::SaveProjectedPotential(ComplexArray2DPtr a){
	_projectedPotential.resize(extents[_c->Model.nx][_c->Model.ny]);
	_projectedPotential = ComplexArray2D(a);
}
void PersistenceManager::Save2DDataSet(ComplexArray2DPtr a, string name){
	_file.SaveComplexArray2D(ComplexArray2D(a),name);
}
void PersistenceManager::StoreToDisc(){
	_file.SaveComplexArray2D(_probe,"probe");
	_file.SaveComplexArray2D(_projectedPotential,"projectedPotential");
	_file.SaveComplexArray3D(_waveSlicesAfterTransmit,"waveSlicesAfterTransmit");
	_file.SaveComplexArray3D(_waveSlicesAfterFT,"waveSlicesAfterFT");
	_file.SaveComplexArray3D(_waveSlicesAfterProp,"waveSlicesAfterProp");
	_file.SaveComplexArray3D(_waveSlicesAfterSlice,"waveSlicesAfterSlice");
	_file.SaveComplexArray3D(_potential,"potentialSlices");
}
} /* namespace QSTEM */

