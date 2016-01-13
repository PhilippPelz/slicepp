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

#ifndef CONVERGENT_WAVE_H
#define CONVERGENT_WAVE_H

#include "wave_base.hpp"
#include "Aberrations.hpp"
namespace slicepp {

class DLL_EXPORT CConvergentWave: public CBaseWave {
public:
	CConvergentWave(cWaveConfPtr wc, cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p);

	virtual void DisplayParams();
	virtual ~CConvergentWave() {
	}
	;
	virtual WavePtr Clone();
	virtual void FormProbe();
	void ConstructWave();

protected:
	float_tt _alpha_max; /* convergence angle */
	bool _smoothen; /* smoothen the probe wave function */
	bool _isGaussian;
	float_tt _sigma;
	float_tt _CLA;

	Aberrations _aberration;
};
}

#endif
