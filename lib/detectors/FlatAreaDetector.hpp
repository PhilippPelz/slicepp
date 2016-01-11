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
class DLL_EXPORT FlatAreaDetector: public IDetector {
public:
	FlatAreaDetector(cDetectorConfPtr dc,PersistenceManagerPtr p);
	void RecordImage(af::array& wave);
	virtual ~FlatAreaDetector();

protected:
	FlatAreaDetector();
	af::array anscombeNoise(af::array&);
	af::array MultiplyMTF(af::array&);
	int _numSaved;
};
}
#endif /* FLATAREADETECTOR_HPP_ */
