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

#include "stemtypes_fftw3.hpp"
#include "structure_IO/crystal.hpp"
#include "scatfactsRez.hpp"
#include "data_IO/PersistenceManager.hpp"
#include "pot_interface.hpp"

//#include <thrust/device_ptr.h>
//#include <thrust/device_vector.h>
//#include <thrust/complex.h>

#include <map>
#include "arrayfire.h"

#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#define OVERSAMPLING 2
#define OVERSAMPLINGZ 3*OVERSAMPLING

namespace slicepp
{
class CPotential : public IPotential
{
public:
	CPotential(cModelConfPtr mc, cOutputConfPtr oc , PersistenceManagerPtr persist);

	void GetSizePixels(unsigned &nx, unsigned &ny) const;
	void WriteSlice(unsigned idx, std::string prefix);
	void WriteProjectedPotential();
	void SetNSlices(unsigned slices);
	complex_tt GetSlicePixel(unsigned iz, unsigned ix, unsigned iy);

	virtual ~CPotential();
	virtual void CenterAtomZ(atom& atom, float_tt &z);
	virtual void DisplayParams();
	virtual void MakeSlices(superCellBoxPtr info);
	virtual void ReadPotential(std::string &fileName, unsigned subSlabIdx);

	virtual af::array GetSlice(af::array& t, unsigned idx);
	virtual af::array GetSubPotential(int startx, int starty, int nx, int ny);
	virtual af::array& GetPotential();
protected:
	CPotential();

	void ResizeSlices();
	void MakePhaseGratings();
	void ReadSlice(const std::string &fileName, ComplexArray2DView slice, unsigned idx);
	void SetScatteringFactors(float_tt kmax);

	virtual void SavePotential();
	virtual void SliceSetup();
	virtual void AddAtomToSlices(atom& atom){};
	virtual void ComputeAtomPotential(int znum){};
	virtual void SaveAtomicPotential(int znum);
	virtual void CleanUp(){};

	//transmission array
	ComplexArray3D _t;
	//transmission array on GPU
	af::array _t_d;

	//map of atomic potential arrays
	std::map<int, ComplexArray2D> _atomPot;

	int _nx,_ny;
	// oversampled resolutions
	float_tt _ddx, _ddy;
	float_tt _dkx,_dky, _dkz,_kmax,_kmax2;
	// total potential thickness
	float_tt _totalThickness;
	// step width in which radial V(r,z) is defined
	float_tt _dr;
	float_tt _atomRadius2;
	float_tt _offsetX, _offsetY;
	// atom radius in units of dx
	int _nRadX;
	// atom radius in units of dy
	int _nRadY;
	// atom radius in units of dz
	int _nRadZ;
	// transverse squared atom radius
	int _nRad2Trans;
	// number of total potential sampling points in x
	int _ndiaAtomX;
	// number of total potential sampling points in y
	int _ndiaAtomY;
	int _boxNx, _boxNy, m_boxNz;
	// atom radius in units of the potential sampling
	unsigned _nrAtomTrans ;
	std::vector<float_tt> _sliceThicknesses;

	float_tt sfLUT(float_tt s,int atKind);
	void splinh( float_tt x[], float_tt y[],  std::vector<float_tt>& b,std::vector<float_tt>& c,std::vector<float_tt>& d, int n);
	float_tt seval( float_tt *x, float_tt *y,std::vector<float_tt>& b,std::vector<float_tt>& c, std::vector<float_tt>& d, int n, float_tt x0 );
};
}

#endif
