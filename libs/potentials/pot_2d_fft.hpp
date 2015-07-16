/*
  QSTEM - image simulation for TEM/STEM/CBED
  Copyright (C) 2000-2010  Christoph Koch
  Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "pot_base.hpp"
#include "fftw++.hpp"

namespace QSTEM
{

class C2DFFTPotential : public CPotential
{
public:
	C2DFFTPotential();
	C2DFFTPotential(const ConfigPtr& c,const PersistenceManagerPtr& persist);
	virtual void DisplayParams();


	virtual void MakeSlices(superCellBoxPtr info);
	virtual void AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ);
protected:
	virtual void AddAtomPeriodic(atom& atom,
			float_tt atomBoxX, int ix,
			float_tt atomBoxY, int iy,
			float_tt atomZ);
	virtual void AddAtomNonPeriodic(atom& atom,
			float_tt atomBoxX, int ix,
			float_tt atomBoxY, int iy,
			float_tt atomZ);
	complex_tt *GetAtomPotential2D(int Znum, float_tt B);
	virtual void SliceSetup();
	virtual void ComputeAtomPotential(int Znum) ;
	virtual void SaveAtomicPotential(int znum);
	virtual void CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z);
	std::map<int, af::array> _atomPot;

	friend class CPotFactory;
	// Create an instance of this class, wrapped in a shared ptr
	//     This should not be inherited - any subclass needs its own implementation.
	//  static PotPtr Create(const ConfigPtr configReader){return PotPtr(new C2DFFTPotential(configReader));}
};

}
