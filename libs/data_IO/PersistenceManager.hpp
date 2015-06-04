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
#include "config_IO/read_qsc.hpp"
namespace QSTEM {

class PersistenceManager {
public:
	PersistenceManager();
	PersistenceManager(const ConfigPtr c);
	void SaveProbe(ComplexArray2DPtr a);
	void SaveWaveAfterTransmit(ComplexArray2DPtr a,int slice);
	void SaveWaveAfterTransform(ComplexArray2DPtr a,int slice);
	void SaveWaveAfterSlice(ComplexArray2DPtr a,int slice);
	void SaveWaveAfterPropagation(ComplexArray2DPtr a,int slice);
	void SavePotential(ComplexArray3D a);
	void SaveProjectedPotential(ComplexArray2DPtr a);
	void Save2DDataSet(ComplexArray2DPtr a, string name);
	void Save3DDataSet(ComplexArray3DPtr a, string name);
	void StoreToDisc();
	void InitStorage();
	virtual ~PersistenceManager();
protected:
	HDFFile _file;
	ConfigPtr _c;
	H5::Group* _info;
	ComplexArray3D _waveSlicesAfterTransmit;
	ComplexArray3D _waveSlicesAfterFT;
	ComplexArray3D _waveSlicesAfterProp,_waveSlicesAfterSlice;
	ComplexArray3D _potential;
	ComplexArray2D _projectedPotential;
	ComplexArray2D _probe;
};
typedef boost::shared_ptr<PersistenceManager> PersistenceManagerPtr;
} /* namespace QSTEM */



#endif /* LIBS_DATA_IO_PERSISTENCEMANAGER_HPP_ */
