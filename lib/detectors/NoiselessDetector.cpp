/*
 * NoiselessDetector.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: philipp
 */

#include "NoiselessDetector.hpp"

namespace slicepp {

NoiselessDetector::NoiselessDetector(cDetectorConfPtr dc,PersistenceManagerPtr p):IDetector(dc,p) {
	// TODO Auto-generated constructor stub

}
void NoiselessDetector::RecordImage(af::array& wave){

}
NoiselessDetector::~NoiselessDetector() {
	// TODO Auto-generated destructor stub
}

} /* namespace slicepp */
