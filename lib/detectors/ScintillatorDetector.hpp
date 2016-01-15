/*
 * FlatAreaDetector.hpp
 *
 *  Created on: Aug 17, 2015
 *      Author: wenxuan
 */

#ifndef FLATAREADETECTOR_HPP_
#define FLATAREADETECTOR_HPP_

#include "detector_interface.hpp"
#include <string>

namespace slicepp{
class DLL_EXPORT ScintillatorDetector: public IDetector {
public:
	ScintillatorDetector(cDetectorConfPtr dc,PersistenceManagerPtr p);
	void RecordImage(af::array& wave);
	virtual ~ScintillatorDetector();

protected:
	ScintillatorDetector();
	af::array anscombeNoise(af::array&);
	af::array MultiplyMTF(af::array&);

};
}
#endif /* FLATAREADETECTOR_HPP_ */
