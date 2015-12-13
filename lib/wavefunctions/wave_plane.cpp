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

namespace QSTEM
{

CPlaneWave::CPlaneWave(const boost::shared_ptr<WaveConfig> wc, const boost::shared_ptr<ModelConfig> mc,	const PersistenceManagerPtr& p) : CBaseWave(wc,mc,p)
{
}

/** Copy constructor - used to copy wave just before dispatching multiple threads for STEM simulations */
CPlaneWave::CPlaneWave(const CPlaneWave& other) : CBaseWave(other)
{
	// TODO: need to copy arrays and anything pointed to - anything that needs to be thread-local
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
			% (_wc->tiltX*RAD2DEG) % (_wc->tiltY *RAD2DEG);
}

void CPlaneWave::FormProbe()
{
	CBaseWave::FormProbe();
	float_tt scale = 1/sqrt((float_tt)(_wc->nx*_wc->ny));
	if ((_wc->tiltX == 0) && (_wc->tiltY == 0)) {
		af::array theta = af::constant(0, _wc->nx, _wc->ny);
		_wave_af = scale*af::complex(cos(theta), sin(theta));
	}
	else {
		TiltBeam();
	}
	_probe = af::array(_wave_af);
}



void CPlaneWave::TiltBeam(bool tiltBack)
{
	if ((_wc->tiltX != 0) || (_wc->tiltY != 0))
	{
		float_tt scale = 1/sqrt(_wc->nx*_wc->ny);
		int direction = tiltBack ? -1 : 1;

		// produce a tilted wave function (btiltx,btilty):
		float_tt ktx = direction*2.0*M_PI*sin(_wc->tiltX)/GetWavelength();
		float_tt kty = direction*2.0*M_PI*sin(_wc->tiltY)/GetWavelength();
		unsigned px=_wc->nx*_wc->ny;
		af::array x = _mc->dx*(af::range(_wc->nx) - _wc->nx/2);
		af::array y = _mc->dy*(af::range(_wc->ny) - _wc->ny/2);
		af::array x2D = af::tile(x.T(), 1, _wc->ny);
		af::array y2D = af::tile(y, _wc->nx);
		_wave_af = af::complex(cos(ktx*x2D + kty*y2D)*scale, sin(ktx*x2D + kty*y2D)*scale);
//		for (unsigned i=0; i<_wc->nx; i++)
//			for (unsigned j=0; j<_wc->ny; j++)
//			{
//				float_tt x = m_dx*(i-_wc->nx/2);
//				float_tt y = m_dy*(j-_wc->ny/2);
//				_wave[i][j]=complex_tt((float_tt)cos(ktx*x+kty*y)*scale,(float_tt)sin(ktx*x+kty*y)*scale);
//			}
	}
}

void CPlaneWave::TiltBack()
{
	TiltBeam(true);
}


} // end namespace QSTEM
