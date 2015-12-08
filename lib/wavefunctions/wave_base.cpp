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
#include <arrayfire.h>
#include <af/util.h>
#include "wave_base.hpp"
using boost::format;

namespace QSTEM {

CBaseWave::CBaseWave(const boost::shared_ptr<WaveConfig> wc, const boost::shared_ptr<ModelConfig> mc, const PersistenceManagerPtr p) :
		IWave(wc, mc, p) {
	_nx = _wc->nx;
	_ny = _wc->ny;
	_dx = _mc->dx;
	_dy = _mc->dy;
	_E = _mc->EnergykeV;
	m_realSpace = true;
	_pixelDose = _wc->pixelDose;
	// TODO: where does this belong?
	//m_electronScale = m_beamCurrent*m_dwellTime*MILLISEC_PICOAMP;
	_wavlen = Wavelength(_E);
}
void CBaseWave::InitializePropagators() {
	float_tt scale = _mc->dz * PI * GetWavelength();
	//t = exp(-i pi lam k^2 dz)
	af::array s, kx2D, ky2D;
// Tile the arrays to create 2D versions of the k vectors
	kx2D = af::tile(m_kx2, 1, _ny);
	ky2D = af::tile(m_ky2.T(), _nx);
	s = scale * (kx2D + ky2D);
	_prop = af::complex(af::cos(s), af::sin(s));

//	TODO:
	_prop = fftShift(_prop);
	_persist->Save2DDataSet(_prop, "Propagator");
	//_prop = fftShift(_prop);
}

af::array CBaseWave::fftShift(af::array wave) {
	return af::shift(wave, wave.dims(0) / 2, wave.dims(1) / 2);
}

af::array CBaseWave::ifftShift(af::array wave) {
	return af::shift(wave, (wave.dims(0) + 1) / 2, (wave.dims(1) + 1) / 2);
}

void CBaseWave::ApplyTransferFunction() {
	// TODO: transfer function should be passed as a 1D vector that is half the size of the wavefunc.
	//       It should be applied by a radial lookup table (with interpolation?)
	//       Alternatively, is it easier to just use a 2D CTF?
	//       Whatever you do, use m_transferFunction as the storage for it.
	int px = GetTotalPixels();

	// multiply wave (in rec. space) with transfer function and write result to imagewave
	ToFourierSpace();
	for (int i = 0; i < _nx; i++)
		for (int j = 0; j < _ny; j++) {
			// here, we apply the CTF:
			// 20140110 - MCS - I think this is where Christoph wanted to apply the CTF - nothing is done ATM.

			// TODO: use these for calculating a radius (to get the CTF value from)
			//ix=i%m_nx;
			//iy=i/m_ny;

			//		wave[i][j] = m_wave[i][j];
		}
	ToRealSpace();
}
void CBaseWave::PropagateToNextSlice() {
	//_wave_af = fftShift(_wave_af);
	_wave_af = _condition * (_wave_af * _prop);
	//_wave_af = ifftShift(_wave_af);
} /* end propagate */

void CBaseWave::Transmit(af::array t_af) {
	_wave_af *= t_af;
} /* end transmit() */
/** Copy constructor - make sure arrays are deep-copied */
CBaseWave::CBaseWave(const CBaseWave &other) :
		CBaseWave(other._wc, other._mc, other._persist) {
	// TODO: make sure arrays are deep copied
	other.GetSizePixels(_nx, _ny);
	other.GetResolution(_dx, _dy);
	_E = other.GetVoltage();
}

CBaseWave::~CBaseWave() {
}

void CBaseWave::InitializeKVectors() {
	float_tt ax = _mc->dx * _nx;
	float_tt by = _mc->dy * _ny;

	m_kx = (af::range(_nx) - _nx / 2) / ax;
	m_kx2 = m_kx * m_kx;

	m_ky = (af::range(_ny) - _ny / 2) / ax;
	m_ky2 = m_ky * m_ky;

	m_k2max = _nx / (2.0F * ax);
	if (_ny / (2.0F * by) < m_k2max)
		m_k2max = _ny / (2.0F * by);
	m_k2max = 2.0 / 3.0 * m_k2max;
	m_k2max = m_k2max * m_k2max;

	GetK2();
	_k2max = af::constant(m_k2max, _nx, _ny);
	_condition = (m_k2 < _k2max);
	_condition = fftShift(_condition);
}

void CBaseWave::GetExtents(int& nx, int& ny) const {
	nx = _nx;
	ny = _ny;
}
void CBaseWave::FormProbe() {

	_zero = af::constant(0, _nx, _ny);
	_wave_af = af::complex(_zero, _zero);
	_wave.resize(boost::extents[_nx][_ny]);
	InitializeKVectors();
}

void CBaseWave::ResetProbe() {
	_wave_af = af::array(_probe);
}

void CBaseWave::DisplayParams() {
	BOOST_LOG_TRIVIAL(info)<<
	"*****************************  Wave  Parameters **************************************************";
	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info)<<format("* Real space res.:      %gA (=%gmrad)")% (1.0/m_k2max)%(GetWavelength()*m_k2max*1000.0);
	BOOST_LOG_TRIVIAL(info)<<format("* Reciprocal space res: dkx=%g, dky=%g (1/A)")% (1.0/(_nx*_mc->dx))%(1.0/(_ny*_mc->dy));

	BOOST_LOG_TRIVIAL(info)<<format("* Beams:                %d x %d ")%_nx%_ny;

	BOOST_LOG_TRIVIAL(info)<<format("* Acc. voltage:         %g (lambda=%gA)")%_E%(Wavelength(_E));

	BOOST_LOG_TRIVIAL(info)<<format("* Probe array:          %d x %d pixels")%_nx%_ny;
	BOOST_LOG_TRIVIAL(info)<<format("*                       %g x %gA")%(_nx*_mc->dx)%(_ny*_mc->dy);
}

/*
 //TODO: where does this belong?
 inline void CBaseWave::GetElectronScale(float_tt &electronScale)
 {
 electronScale=m_electronScale;
 }
 */
void CBaseWave::GetSizePixels(int &x, int &y) const {
	x = _nx;
	y = _ny;
}

void CBaseWave::GetResolution(float_tt &x, float_tt &y) const {
	x = _dx;
	y = _dy;
}

void CBaseWave::GetK2() {
	af::array kx2D, ky2D;
	kx2D = af::tile(m_kx2, 1, _ny);
	ky2D = af::tile(m_ky2.T(), _nx);
	m_k2 = kx2D + ky2D;
}

float_tt CBaseWave::GetIntegratedIntensity() const {
	float_tt intensity;
	intensity = af::sum<float_tt>(af::sqrt(af::real(_wave_af) * af::real(_wave_af) + af::imag(_wave_af) * af::imag(_wave_af)));
	// TODO: divide by px or not?
	return (intensity) / (_nx * _ny);
}

/*--------------------- wavelength() -----------------------------------*/
/*
 return the electron wavelength (in Angstroms)
 keep this is one place so I don't have to keep typing in these
 constants (that I can never remember anyhow)

 ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
 (The American Institute of Physics, New York) 1989
 page 4.

 kev = electron energy in keV

 */
float_tt CBaseWave::Wavelength(float_tt kev) {
	double w;
	const double emass = 510.99906; /* electron rest mass in keV */
	const double hc = 12.3984244; /* Planck's const x speed of light*/

	/* electron wavelength in Angstroms */
	return hc / sqrt(kev * (2 * emass + kev));
} /* end wavelength() */

// FFT to Fourier space, but only if we're current in real space
void CBaseWave::ToFourierSpace() {
	if (IsRealSpace()) {
		m_realSpace = false;
		_wave_af = af::fft2(_wave_af);

	}
}

// FFT back to realspace, but only if we're currently in Fourier space
void CBaseWave::ToRealSpace() {
	if (!IsRealSpace()) {
		m_realSpace = true;
		_wave_af = af::ifft2(_wave_af);
	}
}

} //end namespace QSTEM
