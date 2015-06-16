/*
 * CoherentTEM.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: philiipp
 */

#include "CoherentTEM.hpp"

namespace QSTEM {

CoherentTEM::CoherentTEM(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p,const PersistenceManagerPtr& pers) :
		CoherentSinglePositionExperiment(c,s,w,p,pers){
	// TODO Auto-generated constructor stub
}

CoherentTEM::~CoherentTEM() {
	// TODO Auto-generated destructor stub
}
void CoherentTEM::PostSpecimenProcess(){
	_wave->ApplyTransferFunction();
}
} /* namespace QSTEM */
