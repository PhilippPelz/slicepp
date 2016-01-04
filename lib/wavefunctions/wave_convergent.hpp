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
namespace slicepp
{

class DLL_EXPORT CConvergentWave : public CBaseWave
{
public:
	CConvergentWave(cWaveConfPtr wc, cModelConfPtr mc, PersistenceManagerPtr p);

	virtual void DisplayParams();
	virtual ~CConvergentWave(){};
	virtual WavePtr Clone();
	virtual void FormProbe();
	void ConstructWave();
protected:
	// Coefficients to aberration function:
	float_tt m_a33, m_a31;
	float_tt m_a44, m_a42;
	float_tt m_a55, m_a53, m_a51;
	float_tt m_a66, m_a64, m_a62;
	float_tt m_phi33, m_phi31;
	float_tt m_phi44, m_phi42;
	float_tt m_phi55, m_phi53, m_phi51;
	float_tt m_phi66, m_phi64, m_phi62;

	float_tt m_astigMag;				/* astigmatism*/
	float_tt m_astigAngle;				/* angle of astigmatism */
	float_tt m_C5;
	float_tt m_dE_E;  // energy spread of emitted electrons
	float_tt m_dV_V;  // acc. voltage fluctuations
	float_tt m_dI_I;  // lens current fluctuations
	float_tt m_alpha;   /* convergence angle */
	float_tt m_Cc;      /* chromatic aberration */
	float_tt m_Cs;      /* spherical aberration */
	float_tt m_df0;				/* defocus */

	bool _smoothen;                          /* smoothen the probe wave function */
	bool _isGaussian;
	float_tt _gaussScale;

	float_tt _CLA, m_rmin, m_rmax, m_aimin, m_aimax;

	std::string m_Scherzer;
private:
	friend class CWaveFactory;
	// Create an instance of this class, wrapped in a shared ptr
	//     This should not be inherited - any subclass needs its own implementation.
	//  static WavePtr Create(const ConfigPtr reader){
	//    return WavePtr(new CConvergentWave(reader));
	//  }
};

}

#endif
