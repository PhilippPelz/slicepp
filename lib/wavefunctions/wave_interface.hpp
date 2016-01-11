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

#ifndef WAVE_INTERFACE_H
#define WAVE_INTERFACE_H

#include <data_IO/PersistenceManager.hpp>
#include "config_IO/configs.hpp"
#include "stemtypes_fftw3.hpp"

#include <map>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/function.hpp>
#include <arrayfire.h>

namespace slicepp {
class IWave;
typedef boost::shared_ptr<IWave> WavePtr;
typedef boost::function<IWave*(cWaveConfPtr wc, cModelConfPtr mc, PersistenceManagerPtr p)> waveCreator;
typedef std::map<int, waveCreator> WaveFactory;
class IWave {
public:
	IWave(cWaveConfPtr wc, cModelConfPtr mc, PersistenceManagerPtr p) {
		_persist = p;
		_wc = wc;
		_mc = mc;
	}
	virtual ~IWave() {
	}
	;
	virtual void FormProbe()=0;
	virtual void DisplayParams()=0;
	virtual void ToRealSpace()=0;
	virtual void ToFourierSpace()=0;
	virtual bool IsRealSpace()=0;
	virtual void ApplyTransferFunction()=0;
	virtual WavePtr Clone()=0;
	virtual ComplexArray2D GetWave() const = 0;
	virtual af::array& GetWaveAF() = 0;
	virtual af::array& GetProbe() = 0;
	virtual float_tt GetPixelIntensity(int i) const =0;
	virtual float_tt GetPixelIntensity(int x, int y) const =0;
	virtual float_tt GetIntegratedIntensity() const =0;
	virtual void Transmit(af::array& t)=0;
	virtual void PropagateToNextSlice()=0;
	virtual void InitializePropagators()=0;
	virtual void ResetProbe()=0;
	virtual af::array fftShift(af::array& _wave) =0;
	virtual af::array ifftShift(af::array& _wave) =0;
protected:
	PersistenceManagerPtr _persist;
	cWaveConfPtr _wc;
	cModelConfPtr _mc;
};

}
#endif
