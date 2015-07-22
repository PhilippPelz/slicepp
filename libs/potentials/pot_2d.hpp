#ifndef POTENTIAL_2D_H
#define POTENTIAL_2D_H

#include "RealSpacePotential.hpp"

namespace QSTEM {

class C2DPotential: public RealSpacePotential {
public:
	C2DPotential(const ConfigPtr& c,const PersistenceManagerPtr& persist) ;
	virtual void DisplayParams();
	virtual void AtomBoxLookUp(complex_tt &val, int Znum, float_tt x,float_tt y, float_tt z, float_tt B);
	virtual void AddAtomToSlices(atom& atom );
	virtual void _AddAtomRealSpace(atom& atom, float_tt atomBoxX,unsigned int ix, float_tt atomBoxY,
			unsigned int iy, float_tt atomZ,unsigned int iatomZ);

	bool CheckAtomZInBounds(float_tt atomZ);
	void CenterAtomZ(atom& atom, float_tt &z);
protected:
	virtual void ComputeAtomPotential(int znum);
	virtual void SliceSetup();
	virtual void SaveAtomicPotential(int znum);
};

}

#endif
