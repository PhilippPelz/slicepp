/*
 * CoherentCBED.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: philiipp
 */

#include "CoherentCBED.hpp"

namespace QSTEM {

CoherentCBED::CoherentCBED(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p,const PersistenceManagerPtr& pers) :
				CoherentSinglePositionExperiment(c,s,w,p,pers){
}

CoherentCBED::~CoherentCBED() {
}
void CoherentCBED::PostSpecimenProcess(){
	_wave->ToFourierSpace();
	_wave->fftShift();
}
} /* namespace QSTEM */
