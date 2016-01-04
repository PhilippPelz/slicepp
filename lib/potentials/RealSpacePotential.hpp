/*
 * RealSpacePotential.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: philipp
 */

#ifndef REALSPACEPOTENTIAL_HPP_
#define REALSPACEPOTENTIAL_HPP_

#include "pot_base.hpp"
#include <vector>


namespace slicepp {


class RealSpacePotential : public CPotential
{
public :
	inline RealSpacePotential(cModelConfPtr mc, cOutputConfPtr oc , PersistenceManagerPtr p): CPotential(mc,oc,p){};
	virtual ~RealSpacePotential();
	void AtomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
protected:
	void AddAtomRealSpace(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ);
	virtual void _AddAtomRealSpace(atom& atom,float_tt atomX, unsigned int ix,
			float_tt atomY, unsigned int iy,float_tt atomZ, unsigned int iatomZ)=0;
	virtual void SaveAtomicPotential(int znum)=0;
	std::map<int, atomBoxPtr> m_atomBoxes;
};

}
#endif /* REALSPACEPOTENTIAL_HPP_ */
