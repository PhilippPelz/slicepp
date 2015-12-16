/*
 * PersistenceManager.hpp
 *
 *  Created on: Mar 31, 2015
 *      Author: philipp
 */

#ifndef LIBS_DATA_IO_PERSISTENCEMANAGER_HPP_
#define LIBS_DATA_IO_PERSISTENCEMANAGER_HPP_

#include "HDFFile.hpp"
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include "config_IO/ConfigReader.hpp"
#include "cuda_assert.hpp"
#include "cublas_assert.hpp"
#include "cufft_assert.hpp"
#include <arrayfire.h>

namespace QSTEM {

class PersistenceManager {
public:
	PersistenceManager();
	PersistenceManager(const ConfigPtr c);
	void SaveProbe(ComplexArray2DPtr a);
	void SaveWaveAfterTransmit(ComplexArray2DPtr a, int slice);
	void SaveWaveAfterTransform(ComplexArray2DPtr a, int slice);
	void SaveWaveAfterSlice(ComplexArray2DPtr a, int slice);
	void SaveWaveAfterPropagation(ComplexArray2DPtr a, int slice);
	void SaveProbe(af::array& a);
	void SaveAtomDelta(cuComplex* delta, int slice, int Z);
	void SaveAtomConv(cuComplex* delta, int slice, int Z);
	void SaveWaveAfterTransmit(af::array& wave, int slice);
	void SaveWaveAfterTransform(af::array& wave, int slice);
	void SaveWaveAfterSlice(af::array& wave, int slice);
	void SaveWaveAfterPropagation(af::array& wave, int slice);
	void SavePotential(ComplexArray3D a);
	void SavePotential(af::array& data);
	void SaveProjectedPotential(ComplexArray2DPtr a);
	void SaveProjectedPotential(af::array& data);
	void Save2DDataSet(ComplexArray2DPtr a, string name);
	void Save3DDataSet(ComplexArray3DPtr a, string name);
	void Save2DDataSet(af::array& data, string name);
	void StoreToDisc();
	void StoreToDiscMP(int pos, int x, int y);
	void InitStorage();
	void ResizeStorage(int xdim, int ydim);
	virtual ~PersistenceManager();
	bool _potSaved;
protected:
	HDFFile _file;
	ConfigPtr _c;
	H5::Group* _info;
	ComplexArray3D _waveSlicesAfterTransmit;
	ComplexArray3D _waveSlicesAfterFT;
	ComplexArray3D _waveSlicesAfterProp, _waveSlicesAfterSlice;
	ComplexArray3D _potential;
	std::map<int,boost::shared_ptr<ComplexArray3D>> _atomDeltas, _atomConv;
	ComplexArray2D _projectedPotential;
	ComplexArray2D _probe;
};
typedef boost::shared_ptr<PersistenceManager> PersistenceManagerPtr;
} /* namespace QSTEM */

#endif /* LIBS_DATA_IO_PERSISTENCEMANAGER_HPP_ */
