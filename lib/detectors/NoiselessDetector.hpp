/*
 * NoiselessDetector.hpp
 *
 *  Created on: Jan 14, 2016
 *      Author: philipp
 */

#ifndef NOISELESSDETECTOR_HPP_
#define NOISELESSDETECTOR_HPP_
#include "detector_interface.hpp"
namespace slicepp {

class NoiselessDetector: public IDetector {
public:
	NoiselessDetector(cDetectorConfPtr dc,PersistenceManagerPtr p);
	virtual ~NoiselessDetector();
	void RecordImage(af::array& wave);
};

} /* namespace slicepp */

#endif /* NOISELESSDETECTOR_HPP_ */
