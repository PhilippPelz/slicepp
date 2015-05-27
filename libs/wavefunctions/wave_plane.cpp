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

CPlaneWave::CPlaneWave(const ConfigPtr& c,const PersistenceManagerPtr& p) : CBaseWave(c,p)
{
	_nx = c->Model.nx;
	_ny = c->Model.ny;
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
			% (_config->Wave.tiltX*RAD2DEG) % (_config->Wave.tiltY *RAD2DEG);
}

void CPlaneWave::FormProbe()
{
	float_tt scale = 1/sqrt((float_tt)_nx*_ny);
	if ((_config->Wave.tiltX == 0) && (_config->Wave.tiltY == 0)) {
		for (unsigned i=0; i<_nx; i++)
			for (unsigned j=0; j<_ny; j++)
			{
				_wave[i][j]=complex_tt(scale,scale);

			}
	}
	else {
		TiltBeam();
	}
}

void CPlaneWave::TiltBeam(bool tiltBack)
{
	if ((_config->Wave.tiltX != 0) || (_config->Wave.tiltY != 0))
	{
		float_tt scale = 1/sqrt(_nx*_ny);
		int direction = tiltBack ? -1 : 1;

		// produce a tilted wave function (btiltx,btilty):
		float_tt ktx = direction*2.0*M_PI*sin(_config->Wave.tiltX)/GetWavelength();
		float_tt kty = direction*2.0*M_PI*sin(_config->Wave.tiltY)/GetWavelength();
		unsigned px=_nx*_ny;
		for (unsigned i=0; i<_nx; i++)
			for (unsigned j=0; j<_ny; j++)
			{
				float_tt x = m_dx*(i-_nx/2);
				float_tt y = m_dy*(j-_ny/2);
				_wave[i][j]=complex_tt((float_tt)cos(ktx*x+kty*y)*scale,(float_tt)sin(ktx*x+kty*y)*scale);
			}
	}
}

void CPlaneWave::TiltBack()
{
	TiltBeam(true);
}


} // end namespace QSTEM
