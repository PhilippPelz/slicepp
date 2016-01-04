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
#include <map>


namespace slicepp {

class C3DFFTPotential: public CPotential {
public:
	C3DFFTPotential(cModelConfPtr mc, cOutputConfPtr oc , PersistenceManagerPtr p);
	virtual void DisplayParams();
	virtual void AddAtomToSlices(atom& atom );

protected:
	virtual void AddAtomPeriodic(atom& atom, unsigned ix, unsigned iy);
	virtual void AddAtomNonPeriodic(atom& atom, unsigned int ix, unsigned int iy);
	virtual void ComputeAtomPotential(int Znum);
	virtual void SliceSetup();
	virtual void SaveAtomicPotential(int znum);
	virtual void CleanUp();

	FloatArray2DView getSlicedAtomicPotentialAt(ratio r,int Znum);
	void GetAtomPotentialOffset3D(unsigned Znum, float_tt B, unsigned &nzSub,
			unsigned &Nr, unsigned &Nz_lut, float_tt q, ComplexVector &output);

	// collection of atom potentials
	std::map<unsigned, ComplexArray2D> _atPot;
	std::map<unsigned, ComplexArray2D> _atPotOffset;
	// integrated potential for each fraction of the total slice;
	std::map<ratio,FloatArray2D> _intAtPot;

	// step width of one atomic potential slice;
	float_tt _dz;
	// number of atomic potential slices per total potential slice
	unsigned _nzPerSlice;
	// atom radius in units of the atomic potential sampling
	unsigned _nrAtomZ;
	// atom diameter in units of the atomic potential sampling
	unsigned _ndiaAtomZ;
	// current Znum for which _intAtPots are stored intAtPot is refreshed once a new Znum is found
	unsigned _currentZnum;


	int getOffsetZ(float_tt d);
};

}
