/*
 * CoherentTEM.hpp
 *
 *  Created on: Jun 12, 2015
 *      Author: philiipp
 */

#ifndef LIBS_EXPERIMENTS_COHERENTTEM_HPP_
#define LIBS_EXPERIMENTS_COHERENTTEM_HPP_

#include "CoherentSinglePositionExperiment.hpp"
#include "Aberrations.hpp"

namespace slicepp {

class CoherentTEM: public CoherentSinglePositionExperiment {
public:
	CoherentTEM(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p, const DetPtr& d, const PersistenceManagerPtr& pers);
	virtual ~CoherentTEM();
	virtual void PostSpecimenProcess();
protected:
	Aberrations _aberration;
};

} /* namespace slicepp */

#endif /* LIBS_EXPERIMENTS_COHERENTTEM_HPP_ */
