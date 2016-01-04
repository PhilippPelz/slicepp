#ifndef POTENTIAL_3D_H
#define POTENTIAL_3D_H

#include "RealSpacePotential.hpp"

namespace slicepp
{

class C3DPotential : public RealSpacePotential
{
public:
  C3DPotential(cModelConfPtr mc, cOutputConfPtr oc , PersistenceManagerPtr p);
  virtual void DisplayParams();
  virtual void AtomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
  //virtual void makeSlices(int nlayer, char *fileName, atom *center);
  void CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z);
  bool CheckAtomZInBounds(float_tt atomZ);
  virtual void AddAtomToSlices(atom& atom );
  void _AddAtomRealSpace(atom &atom,
                         float_tt atomBoxX, unsigned ix,
                         float_tt atomBoxY, unsigned iy,
                         float_tt atomZ, unsigned iAtomZ);
protected:
  unsigned m_sliceStep;  // number of float_tt to advance to next slice (2*m_nx*m_ny)
  virtual void ComputeAtomPotential(int Znum);
  virtual void SaveAtomicPotential(int znum);
  virtual void SliceSetup();

  float_tt _ddz;

  friend class CPotFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
//  static PotPtr Create(const ConfigPtr configReader){return PotPtr(new C3DPotential(configReader));};
};

}

#endif
