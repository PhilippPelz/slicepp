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

#include "omp.h"
#include "boost/format.hpp"
#include <boost/log/trivial.hpp>
#include <arrayfire.h>
#include <af/util.h>

#include "stemtypes_fftw3.hpp"
#include "wave_interface.hpp"
#include "data_IO/PersistenceManager.hpp"
namespace slicepp
{


class DLL_EXPORT CBaseWave : public IWave
{
public:
	CBaseWave(cWaveConfPtr wc, cModelConfPtr mc, PersistenceManagerPtr p);
	bool IsRealSpace(){return _realSpace;}
	int GetTotalPixels() const {return _wc->n[0]*_wc->n[1];}
	void DisplayParams();
	void ToRealSpace();
	void ToFourierSpace();
	inline float_tt GetPixelIntensity(int i) const { return abs2(_wave.data()[i]); }
	inline float_tt GetPixelIntensity(int x, int y) const  {return GetPixelIntensity(x+_wc->n[0]*y);}
	inline ComplexArray2D GetWave() const {return  _wave;}

	inline af::array GetPixelIntensity() const {return af::real(_wave_af)*af::real(_wave_af) + af::imag(_wave_af)*af::imag(_wave_af);}
	inline af::array& GetWaveAF() {return  _wave_af;}
	inline af::array& GetProbe() {return  _probe;}

	void WriteBeams(int absoluteSlice);
	float_tt GetIntegratedIntensity() const ;

	virtual ~CBaseWave();
	virtual WavePtr Clone()=0;
	virtual void FormProbe();
	virtual void ResetProbe();
	virtual void Transmit(af::array& t);
	virtual void PropagateToNextSlice();
	virtual void InitializePropagators();
	virtual af::array fftShift(af::array& _wave);
	virtual af::array ifftShift(af::array& _wave);
	virtual void ApplyTransferFunction();

protected:
	ComplexArray2D _wave;
	af::array _prop;
	af::array _wave_af;
	af::array _probe;
	af::array _condition;
	bool _realSpace;
	float_tt _k2max;

	af::array _kx;
	af::array _ky;
	af::array _k2;

	void Initialize(std::string input_ext, std::string output_ext);
	void InitializeKVectors();

	af::array GetKx();
	af::array GetKy();
	af::array GetK2();
};
}

#endif
