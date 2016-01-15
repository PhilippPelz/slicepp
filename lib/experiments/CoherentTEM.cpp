/*
 * CoherentTEM.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: philiipp
 */

#include "CoherentTEM.hpp"
#include "afhelpers.hpp"

namespace slicepp {

CoherentTEM::CoherentTEM(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p, const DetPtr& d, const PersistenceManagerPtr& pers) :
		_aberration(c->Model->wavelength*1e10,c->Model->OLaberrations,"Objective Lens Aberrations:"),
		CoherentSinglePositionExperiment(c,s,w,p,d,pers){
}

CoherentTEM::~CoherentTEM() {
}
void CoherentTEM::PostSpecimenProcess(){
	_wave->ApplyCTF();
	_det->RecordImage(_wave->GetWaveAF());
}
} /* namespace slicepp */
