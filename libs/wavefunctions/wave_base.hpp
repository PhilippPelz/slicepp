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

#ifndef WAVE_BASE_H
#define WAVE_BASE_H

#include "stemtypes_fftw3.hpp"
#include "wave_interface.hpp"
#include "omp.h"
#include "boost/format.hpp"
#include <boost/log/trivial.hpp>
#include "data_IO/PersistenceManager.hpp"
#include "fftw++.hpp"

namespace QSTEM
{

class QSTEM_HELPER_DLL_EXPORT CBaseWave : public IWave
{
public:
	CBaseWave(const ConfigPtr& c,const PersistenceManagerPtr& p);
	CBaseWave( const CBaseWave& other );

	void DisplayParams();
	void ToRealSpace();
	void ToFourierSpace();
	bool IsRealSpace(){return m_realSpace;}
	void GetSizePixels(int &x, int &y) const ;
	int GetTotalPixels() const {return _nx*_ny;}
	void GetResolution(float_tt &x, float_tt &y) const ;
	void GetPositionOffset(int &x, int &y) const ;
	float_tt GetK2(int ix, int iy) const ;
	inline float_tt GetKX2(int ix) const {return m_kx2[ix];}
	inline float_tt GetKY2(int iy) const {return m_ky2[iy];}
	inline float_tt GetK2Max() const {return m_k2max;}
	inline float_tt GetVoltage()  const {return m_v0;}
	inline float_tt GetWavelength()  const {return m_wavlen;}
	inline ComplexArray2D GetWave() const {return  _wave;}
	inline float_tt GetPixelIntensity(int i) const {
		return abs2(_wave.data()[i]);
	}
	inline float_tt GetPixelIntensity(int x, int y) const  {return GetPixelIntensity(x+_nx*y);}

	void WriteBeams(int absoluteSlice);
	float_tt GetIntegratedIntensity() const ;

	virtual void ApplyTransferFunction(std::vector<complex_tt> &wave);
	virtual ~CBaseWave();
	virtual WavePtr Clone()=0;
	virtual void FormProbe();
	virtual void GetExtents(int& nx, int& ny) const;
	virtual void Transmit(ComplexArray2DView t);
	virtual void PropagateToNextSlice();
	virtual void InitializePropagators();
	virtual void ShiftTo(float_tt x, float_tt y);

protected:
	ComplexArray2D _prop, _wave;
	std::vector<float_tt> m_propxr, m_propxi, m_propyr, m_propyi;
	bool m_realSpace;
	int m_detPosX, m_detPosY;
	int _nx, _ny;		      /* size of wavefunc and diffpat arrays */
	int m_Scherzer;
	int m_printLevel;
	float_tt m_dx, m_dy;  // physical pixel size of wavefunction array
	float_tt m_k2max;
	float_tt m_v0;
	float_tt m_wavlen;

	std::vector<int> m_position;
	RealVector m_kx2,m_ky2,m_kx,m_ky;

	//float_tt **m_avgArray;
	//float_tt m_thickness;
	//float_tt m_intIntensity;
	//float_tt m_electronScale;
	//float_tt m_beamCurrent;
	//float_tt m_dwellTime;

	fftwpp::fft2d _forward, _backward;

	void Initialize(std::string input_ext, std::string output_ext);
	void InitializeKVectors();
	float_tt Wavelength(float_tt keV);

};
}

#endif
