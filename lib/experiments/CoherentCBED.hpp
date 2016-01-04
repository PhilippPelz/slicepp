/*
 * CoherentCBED.hpp
 *
 *  Created on: Jun 12, 2015
 *      Author: philiipp
 */

#ifndef LIBS_EXPERIMENTS_COHERENTCBED_HPP_
#define LIBS_EXPERIMENTS_COHERENTCBED_HPP_

#include "CoherentSinglePositionExperiment.hpp"

namespace slicepp {

class CoherentCBED: public CoherentSinglePositionExperiment {
public:
	CoherentCBED(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p, const DetPtr& d, const PersistenceManagerPtr& pers);
	virtual ~CoherentCBED();
	//virtual void PostSpecimenProcess();
};

} /* namespace slicepp */

#endif /* LIBS_EXPERIMENTS_COHERENTCBED_HPP_ */
