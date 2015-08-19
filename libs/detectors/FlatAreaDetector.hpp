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

namespace QSTEM{
class DLL_EXPORT FlatAreaDetector: public IDetector {
public:
	FlatAreaDetector(const ConfigPtr& c,const PersistenceManagerPtr& persist);
	void RecordImage(WavePtr w);
	void MultiplyMTF(af::array wave);
	virtual ~FlatAreaDetector();


protected:
	FlatAreaDetector();
	void anscombeNoise(af::array slice, float_tt dose);
	int _nx, _ny, _numSaved;
	af::array _image;
};
}
#endif /* FLATAREADETECTOR_HPP_ */
