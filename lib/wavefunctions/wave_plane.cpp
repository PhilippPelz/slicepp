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

#include "wave_plane.hpp"

static std::string imageFilePrefix="image";
using boost::format;

namespace slicepp
{

CPlaneWave::CPlaneWave(cWaveConfPtr wc, cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) : CBaseWave(wc,mc,oc,p)
{
}

WavePtr CPlaneWave::Clone()
{
	return WavePtr(new CPlaneWave(*this));
}

void CPlaneWave::DisplayParams()
{
	CBaseWave::DisplayParams();

	// TODO: transmit tiltBack status somehow (it's an experiment parameter, not a wave parameter.
	BOOST_LOG_TRIVIAL(info) << format("* Beam tilt:            x=%g deg, y=%g deg")
			% (_wc->tilt[0]*RAD2DEG) % (_wc->tilt[1] *RAD2DEG);
}

void CPlaneWave::FormProbe()
{
	CBaseWave::FormProbe();
	float_tt scale = 1/sqrt((float_tt)(_wc->n[0]*_wc->n[1]));
	if ((_wc->tilt[0] == 0) && (_wc->tilt[1] == 0)) {
		af::array theta = af::constant(0, _wc->n[0], _wc->n[1]);
		_wave_af = scale*af::complex(cos(theta), sin(theta));
	}
	else {
		TiltBeam();
	}
	_probe = _wave_af.copy();
}



void CPlaneWave::TiltBeam(bool tiltBack)
{
	if ((_wc->tilt[0] != 0) || (_wc->tilt[1] != 0))
	{
		float_tt scale = 1/sqrt(_wc->n[0]*_wc->n[1]);
		int direction = tiltBack ? -1 : 1;

		// produce a tilted wave function (btiltx,btilty):
		float_tt ktx = direction*2.0*M_PI*sin(_wc->tilt[0])/_mc->wavelength;
		float_tt kty = direction*2.0*M_PI*sin(_wc->tilt[1])/_mc->wavelength;
		unsigned px=_wc->n[0]*_wc->n[1];
		af::array x = _mc->d[0]*(af::range(_wc->n[0]) - _wc->n[0]/2);
		af::array y = _mc->d[1]*(af::range(_wc->n[1]) - _wc->n[1]/2);
		af::array x2D = af::tile(x.T(), 1, _wc->n[1]);
		af::array y2D = af::tile(y, _wc->n[0]);
		_wave_af = af::complex(cos(ktx*x2D + kty*y2D)*scale, sin(ktx*x2D + kty*y2D)*scale);
	}
}

void CPlaneWave::TiltBack()
{
	TiltBeam(true);
}


} // end namespace slicepp
