/*
 * DirectDetector.hpp
 *
 *  Created on: Jan 14, 2016
 *      Author: philipp
 */

#ifndef DIRECTDETECTOR_HPP_
#define DIRECTDETECTOR_HPP_
#include "detector_interface.hpp"
namespace slicepp {

class DirectDetector: public IDetector {
public:
	DirectDetector(cDetectorConfPtr dc,PersistenceManagerPtr p);
	virtual ~DirectDetector();
	void RecordImage(af::array& wave);
protected:
	af::array Noise(af::array&);
};

} /* namespace slicepp */

#endif /* DIRECTDETECTOR_HPP_ */
