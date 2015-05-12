#ifndef POTENTIAL_2D_H
#define POTENTIAL_2D_H

#include "RealSpacePotential.hpp"

namespace QSTEM {

class C2DPotential: public RealSpacePotential {
public:
	C2DPotential(const ConfigPtr c, PersistenceManagerPtr p) ;
	virtual void DisplayParams();
	virtual void AtomBoxLookUp(complex_tt &val, int Znum, float_tt x,float_tt y, float_tt z, float_tt B);
	bool CheckAtomZInBounds(float_tt atomZ);
	void CenterAtomZ(atom& atom, float_tt &z);
	virtual void AddAtomToSlices(atom& atom,float_tt atomX, float_tt atomY, float_tt atomZ);
	virtual void _AddAtomRealSpace(atom& atom, float_tt atomBoxX,unsigned int ix, float_tt atomBoxY,
			unsigned int iy, float_tt atomZ,unsigned int iatomZ);
protected:
	virtual void ComputeAtomPotential(int znum);
	virtual void SliceSetup();
	friend class CPotFactory;
	// Create an instance of this class, wrapped in a shared ptr
	//     This should not be inherited - any subclass needs its own implementation.
//	static PotPtr Create(const ConfigPtr configReader) {
//		return PotPtr(new C2DPotential(configReader));
//	}
};

}

#endif
