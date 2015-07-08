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

#include "../stemtypes_fftw3.hpp"
#include "structure_IO/crystal.hpp"
#include "scatfactsRez.hpp"
#include <data_IO/PersistenceManager.hpp>
#include "pot_interface.hpp"
#include <map>
#include "arrayfire.h"

#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#define OVERSAMPLING 2
#define OVERSAMPLINGZ 3*OVERSAMPLING

namespace QSTEM
{
class CPotential : public IPotential
{
public:
	CPotential(const ConfigPtr& c,const PersistenceManagerPtr& persist);

	void GetSizePixels(unsigned &nx, unsigned &ny) const;
	void WriteSlice(unsigned idx, std::string prefix);
	void WriteProjectedPotential();
	void SetSliceThickness(float_tt thickness_Angstroms);
	void SetSliceThickness(std::vector<float_tt> thickness_Angstroms);
	void SetNSlices(unsigned slices);
	complex_tt GetSlicePixel(unsigned iz, unsigned ix, unsigned iy);

	virtual ~CPotential();
	virtual void CenterAtomZ(atom& atom, float_tt &z);
	virtual void DisplayParams();
	virtual void MakeSlices(superCellBoxPtr info);
	virtual void Refresh();
	virtual void ReadPotential(std::string &fileName, unsigned subSlabIdx);
	virtual void SetStructure(StructurePtr structure);

//	inline ComplexArray2DView GetSlice(unsigned idx){return _t[boost::indices[idx][range(0,_c->Model.nx)][range(0,_c->Model.ny)]];}
	af::array GetSlice(unsigned idx);
protected:
	CPotential();
	virtual void SliceSetup();
	void ResizeSlices();
	void MakePhaseGratings();
	void BandlimitTransmissionFunction();
	void ReadSlice(const std::string &fileName, ComplexArray2DView slice, unsigned idx);

	virtual void AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ)=0;
	virtual void ComputeAtomPotential(int znum)=0;
	virtual void SaveAtomicPotential(int znum)=0;

	StructureBuilderPtr _sb;
	ComplexArray3D _t;
	af::array _t_af;

	float_tt _ddx, _ddy, _ddz;   // oversampled resolutions
	float_tt _dkx,_dky, _dkz,_kmax2;
	float_tt _totalThickness;
	float_tt _dr, _atomRadius2;
	float_tt _offsetX, _offsetY;
	int _iRadX, _iRadY, _iRadZ, _iRad2;
	int _nx,_ny;
	int _boxNx, _boxNy, m_boxNz;
	int m_scatFactor = 0 ;  // TODO: NOT USED YET The scattering factor type.  One of: 0 (Doyle-Turner); 1 (Wieck-Kohl); 2 (Custom)

	std::vector<float_tt> _sliceThicknesses;
	std::vector<float_tt> _slicePos;

	float_tt sfLUT(float_tt s,int atKind);
	void splinh( float_tt x[], float_tt y[],  std::vector<float_tt>& b,std::vector<float_tt>& c,std::vector<float_tt>& d, int n);
	float_tt seval( float_tt *x, float_tt *y,std::vector<float_tt>& b,std::vector<float_tt>& c,
			std::vector<float_tt>& d, int n, float_tt x0 );
};
}

#endif
