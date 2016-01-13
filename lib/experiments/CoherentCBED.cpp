/*
 * CoherentCBED.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: philiipp
 */

#include "CoherentCBED.hpp"
#include "afhelpers.hpp"

namespace slicepp {

CoherentCBED::CoherentCBED(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p, const DetPtr& d, const PersistenceManagerPtr& pers) :
				CoherentSinglePositionExperiment(c,s,w,p,d,pers){
}

CoherentCBED::~CoherentCBED() {
}
void CoherentCBED::PostSpecimenProcess(){
	_wave->ToFourierSpace();
	auto dp = fftShift(_wave->GetWaveAF());
	_det->RecordImage(dp);
}
} /* namespace slicepp */
