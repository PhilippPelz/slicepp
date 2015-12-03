/*
 * detector_interface.hpp
 *
 *  Created on: Aug 17, 2015
 *      Author: wenxuan
 */

#ifndef DETECTOR_INTERFACE_HPP_
#define DETECTOR_INTERFACE_HPP_

#include "float.h"
#include "stemtypes_fftw3.hpp"
#include "config_IO/config_interface.hpp"
#include "config_IO/configs.hpp"
#include "data_IO/PersistenceManager.hpp"
#include "wavefunctions/wave_interface.hpp"
#include "arrayfire.h"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/function.hpp>

using boost::format;

namespace QSTEM{

class IDetector;
typedef boost::shared_ptr<IDetector> DetPtr;
typedef boost::function<IDetector*(const ConfigPtr& c,const PersistenceManagerPtr &persist)> detectorCreator;
typedef std::map<int,detectorCreator> DetectorFactory;

class IDetector
{
public:
	IDetector(const ConfigPtr& c,const PersistenceManagerPtr& persist) : _c(c), _persist(persist){};
	virtual void RecordImage(WavePtr w)=0;
	virtual void MultiplyMTF(af::array wave)=0;
	virtual ~IDetector(){};

protected:
	IDetector(){};
	ConfigPtr _c;
	PersistenceManagerPtr _persist;
};
}

#endif /* DETECTOR_INTERFACE_HPP_ */