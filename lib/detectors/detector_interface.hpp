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
#include "config_IO/configs.hpp"
#include "data_IO/PersistenceManager.hpp"
#include "wavefunctions/wave_interface.hpp"
#include "arrayfire.h"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/function.hpp>

using boost::format;

namespace slicepp{

class IDetector;
typedef boost::shared_ptr<IDetector> DetPtr;
typedef boost::function<IDetector*(cDetectorConfPtr dc,PersistenceManagerPtr p)> detectorCreator;
typedef std::map<DetectorType,detectorCreator> DetectorFactory;

class IDetector
{
public:
	IDetector(cDetectorConfPtr dc,PersistenceManagerPtr p) : _dc(dc), _p(p){};
	virtual void RecordImage(af::array& wave)=0;
	virtual ~IDetector(){};

protected:
	IDetector(){};
	cDetectorConfPtr _dc;
	PersistenceManagerPtr _p;
};
}

#endif /* DETECTOR_INTERFACE_HPP_ */
