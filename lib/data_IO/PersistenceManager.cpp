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

namespace slicepp {
PersistenceManager::PersistenceManager() :
		_waveSlicesAfterFT(ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterTransmit(ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterProp(
				ComplexArray3D(boost::extents[1][1][1])), _waveSlicesAfterSlice(ComplexArray3D(boost::extents[1][1][1])), _probe(
				ComplexArray2D(boost::extents[1][1])), _projectedPotential(ComplexArray2D(boost::extents[1][1])), _info(NULL), _potSaved(false) {
}
PersistenceManager::PersistenceManager(const ConfigPtr c) :
		PersistenceManager() {
	_c = c;
	_file = HDFFile(_c->Output->SavePath);
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
	if (a.get() == NULL) {
		auto e3 = boost::extents[_c->Model->n[2]][_c->Model->n[0]][_c->Model->n[1]];
		_atomDeltas[Z] = boost::shared_ptr<ComplexArray3D>(new ComplexArray3D(e3));
		auto sh = _atomDeltas[Z]->shape();
//		BOOST_LOG_TRIVIAL(info)<< format("%s shape: [%d,%d,%d]") % "_atomDeltas" % sh[0] % sh[1] % sh[2];
	}
	cuda_assert(
			cudaMemcpy(_atomDeltas[Z]->operator [](slice).origin(), delta, _c->Model->n[0] * _c->Model->n[1] * sizeof(cuComplex),
					cudaMemcpyDeviceToHost));
}
void PersistenceManager::SaveAtomConv(cuComplex* delta, int slice, int Z) {
	auto a = _atomConv[Z];
	if (a.get() == NULL) {
		auto e3 = boost::extents[_c->Model->n[2]][_c->Model->n[0]][_c->Model->n[1]];
		_atomConv[Z] = boost::shared_ptr<ComplexArray3D>(new ComplexArray3D(e3));
		auto sh = _atomConv[Z]->shape();
	}
	cuda_assert(
			cudaMemcpy(_atomConv[Z]->operator [](slice).origin(), delta, _c->Model->n[0] * _c->Model->n[1] * sizeof(cuComplex),
					cudaMemcpyDeviceToHost));
}
void PersistenceManager::SaveMeasurement(af::array& m, int n) {
	m.host(_measurements[n].origin());
}
void PersistenceManager::SaveWaveAfterTransmit(ComplexArray2DPtr a, int slice) {
#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterTransmit[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterTransmit(af::array& wave, int slice) {
	wave.host(_waveSlicesAfterTransmit[slice].origin());
}
void PersistenceManager::SaveWaveAfterTransform(ComplexArray2DPtr a, int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterFT[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterTransform(af::array& wave, int slice) {
	wave.host(_waveSlicesAfterFT[slice].origin());
}
void PersistenceManager::SaveWaveAfterPropagation(ComplexArray2DPtr a, int slice) {
#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterProp[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterPropagation(af::array& wave, int slice) {
	wave.host(_waveSlicesAfterProp[slice].origin());
}
void PersistenceManager::SaveWaveAfterSlice(ComplexArray2DPtr a, int slice) {

#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a.num_elements() / a.size(); j++) {
			_waveSlicesAfterSlice[slice][i][j] = a[i][j];
		}
}
void PersistenceManager::SaveWaveAfterSlice(af::array& wave, int slice) {
	wave.host(_waveSlicesAfterSlice[slice].origin());
}
void PersistenceManager::SavePotential(ComplexArray3D a) {
	_potSaved = true;
	_potential = a;
}
void PersistenceManager::SavePositions(std::vector<int2> pos){
	FloatArray2D p(boost::extents[pos.size()][2]);
	for(int i = 0; i < pos.size();i++){
		p[i][0] = (float)pos[i].x;
		p[i][1] = (float)pos[i].y;
	}
	_file.SaveRealArray2D(p,"positions");
}
void PersistenceManager::SaveZnums(std::vector<int> Z){
	FloatArray2D p(boost::extents[Z.size()]);
	for(int i = 0; i < Z.size();i++){
		p[i] = (float)Z[i];
	}
	_file.SaveRealArray1D(p,"Znums");
}
void PersistenceManager::StoreMeasurements(){
	_file.SaveRealArray3D(_measurements, "measurements");
}
void PersistenceManager::SavePotential(af::array& data) {
	data.host(_potential.data());
}

void PersistenceManager::SaveProjectedPotential(ComplexArray2DPtr a) {
	_projectedPotential = ComplexArray2D(a);
}

void PersistenceManager::SaveProjectedPotential(af::array& data) {
	data.host(_projectedPotential.data());
}

void PersistenceManager::SaveCx3DDataSet(ComplexArray3DPtr a, string name) {
	_file.SaveComplexArray3D(ComplexArray3D(a), name);
}
void PersistenceManager::SaveF3DDataSet(FloatArray3DPtr a, string name) {
	_file.SaveRealArray3D(FloatArray3D(a), name);
}
void PersistenceManager::SaveCx2DDataSet(ComplexArray2DPtr a, string name) {
	_file.SaveComplexArray2D(ComplexArray2D(a), name);
}
void PersistenceManager::SaveCx2DDataSet(ComplexArray2DView2D a, string name) {
	_file.SaveComplexArray2D(ComplexArray2D(a), name);
}
void PersistenceManager::SaveF2DDataSet(FloatArray2DPtr a, string name) {
	_file.SaveRealArray2D(FloatArray2D(a), name);
}

void PersistenceManager::Save2DDataSet(af::array& d, string name) {
	printf("Saving %s: dims(%d,%d)\n", name.c_str(), d.dims(0), d.dims(1));
	if (d.type() == f32) {
		FloatArray2D a(boost::extents[d.dims(0)][d.dims(1)]);
		d.host(a.data());
		SaveF2DDataSet(a, name);
	} else if (d.type() == c32) {
		ComplexArray2D b(boost::extents[d.dims(0)][d.dims(1)]);
		d.host(b.data());
		SaveCx2DDataSet(b, name);
	}
}

void PersistenceManager::InitStorage() {
	auto e3 = boost::extents[_c->Model->n[2]][_c->Model->n[0]][_c->Model->n[1]];
	auto e2 = boost::extents[_c->Model->n[0]][_c->Model->n[1]];
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
	if (_c->Output->SaveProbe)
		_probe.resize(e2);
	if (_c->Output->SaveAtomDeltas)
		for (auto& kv : _atomDeltas)
			kv.second->resize(e3);
	if (_c->Output->SaveAtomConv)
		for (auto& kv : _atomConv)
			kv.second->resize(e3);
	if (_c->Output->SaveProjectedPotential || _c->Output->ComputeFromProjectedPotential)
		_projectedPotential.resize(e2);
	_measurements.resize(boost::extents[1][_c->Model->n[0]][_c->Model->n[1]]);
}

void PersistenceManager::ResizeStorage(int xdim, int ydim) {
	auto e3 = boost::extents[_c->Model->n[2]][xdim][ydim];
	auto e2 = boost::extents[xdim][ydim];
	if (_c->Output->SavePotential || _c->Output->ComputeFromProjectedPotential)
		_potential.resize(boost::extents[_c->Model->n[2]][_c->Model->n[0]][_c->Model->n[1]]);
	if (_c->Output->SaveWaveAfterSlice)
		_waveSlicesAfterSlice.resize(e3);
	if (_c->Output->SaveWaveAfterPropagation)
		_waveSlicesAfterProp.resize(e3);
	if (_c->Output->SaveWaveAfterTransform)
		_waveSlicesAfterFT.resize(e3);
	if (_c->Output->SaveWaveAfterTransmit)
		_waveSlicesAfterTransmit.resize(e3);
	if (_c->Output->SaveProbe)
		_probe.resize(e2);
	if (_c->Output->SaveAtomDeltas)
		for (auto& kv : _atomDeltas)
			kv.second->resize(e3);
	if (_c->Output->SaveAtomConv)
		for (auto& kv : _atomConv)
			kv.second->resize(e3);
	if (_c->Output->SaveProjectedPotential || _c->Output->ComputeFromProjectedPotential)
		_projectedPotential.resize(boost::extents[_c->Model->n[0]][_c->Model->n[1]]);
	_measurements.resize(boost::extents[_c->Scan->positions.size()][xdim][ydim]);
}

void PersistenceManager::StoreToDisc() {
	_file.SaveComplexArray2D(_probe, "probe");
	_file.SaveComplexArray2D(_projectedPotential, "projectedPotential");
	_file.SaveComplexArray3D(_waveSlicesAfterTransmit, "waveSlicesAfterTransmit");
	_file.SaveComplexArray3D(_waveSlicesAfterFT, "waveSlicesAfterFT");
	_file.SaveComplexArray3D(_waveSlicesAfterProp, "waveSlicesAfterProp");
	_file.SaveComplexArray3D(_waveSlicesAfterSlice, "waveSlicesAfterSlice");
	_file.SaveComplexArray3D(_potential, "potentialSlices");
	_file.SaveRealArray3D(_measurements, "measurements");
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
	} else {
//		_file.SaveComplexArray3D(_waveSlicesAfterTransmit,"waveSlicesAfterTransmit" + info);
//		_file.SaveComplexArray3D(_waveSlicesAfterFT, "waveSlicesAfterFT");
//		_file.SaveComplexArray3D(_waveSlicesAfterProp, "waveSlicesAfterProp" + info);
//		_file.SaveComplexArray3D(_waveSlicesAfterSlice, "waveSlicesAfterSlice" + info);
	}
}

} /* namespace slicepp */

