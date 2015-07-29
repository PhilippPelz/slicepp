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

#include <stdio.h>
#include <cstdlib>
#include "stemtypes_fftw3.hpp"
#include "wave_interface.hpp"
#include "omp.h"
#include "boost/format.hpp"
#include <boost/log/trivial.hpp>
#include "data_IO/PersistenceManager.hpp"
#include "fftw++.hpp"
#include <arrayfire.h>
#include <af/util.h>

namespace QSTEM
{

class DLL_EXPORT CBaseWave : public IWave
{
public:
	CBaseWave(const ConfigPtr& c,const PersistenceManagerPtr& p);
	CBaseWave( const CBaseWave& other );
	bool IsRealSpace(){return m_realSpace;}
	int GetTotalPixels() const {return _nx*_ny;}
	void DisplayParams();
	void ToRealSpace();
	void ToFourierSpace();
	void GetSizePixels(int &x, int &y) const ;
	void GetResolution(float_tt &x, float_tt &y) const ;
	void GetPositionOffset(int &x, int &y) const;
	void SetPositionOffset(int x, int y);
	virtual void GetK2();
	inline float_tt GetVoltage()  const {return m_v0;}
	inline float_tt GetWavelength()  const {return m_wavlen;}
	inline float_tt GetPixelIntensity(int i) const { return abs2(_wave.data()[i]); }
	inline float_tt GetPixelIntensity(int x, int y) const  {return GetPixelIntensity(x+_nx*y);}
	inline ComplexArray2D GetWave() const {return  _wave;}

	inline af::array GetPixelIntensity() const {return af::real(_wave_af)*af::real(_wave_af) + af::imag(_wave_af)*af::imag(_wave_af);}
	inline af::array GetWaveAF() const {return  _wave_af;}
	inline af::array GetProbe() const {return  _probe;}

	void WriteBeams(int absoluteSlice);
	float_tt GetIntegratedIntensity() const ;

	virtual ~CBaseWave();
	virtual WavePtr Clone()=0;
	virtual void FormProbe();
	virtual void ResetProbe();
	virtual void GetExtents(int& nx, int& ny) const;
	virtual void Transmit(af::array t);
	virtual void PropagateToNextSlice();
	virtual void InitializePropagators();
	virtual void ShiftTo(float_tt x, float_tt y);
	virtual af::array fftShift(af::array _wave);
	virtual af::array ifftShift(af::array _wave);
	virtual void ApplyTransferFunction();
protected:
	ComplexArray2D _wave;
	af::array _prop, _wave_af, _probe;
	af::array _condition, _zero; //for propagation
	std::vector<float_tt> m_propxr, m_propxi, m_propyr, m_propyi;
	bool m_realSpace;
	int m_detPosX, m_detPosY;
	int _nx, _ny;		      /* size of wavefunc and diffpat arrays */
	int m_Scherzer;
	int m_printLevel;
	int slice_c;
	float_tt m_dx, m_dy;  // physical pixel size of wavefunction array
	float_tt m_k2max;
	float_tt m_v0;
	float_tt m_wavlen;

	std::vector<int> m_position;
	af::array m_kx2, m_ky2, m_kx, m_ky, m_k2, _k2max;
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
