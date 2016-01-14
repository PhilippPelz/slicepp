/*
 * DirectDetector.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: philipp
 */

#include "DirectDetector.hpp"

namespace slicepp {

DirectDetector::DirectDetector(cDetectorConfPtr dc,PersistenceManagerPtr p):IDetector(dc,p) {
	// TODO Auto-generated constructor stub

}
void DirectDetector::RecordImage(af::array& wave){

}
DirectDetector::~DirectDetector() {
	// TODO Auto-generated destructor stub
}

} /* namespace slicepp */
