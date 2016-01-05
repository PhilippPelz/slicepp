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

namespace slicepp {

CBaseWave::CBaseWave(cWaveConfPtr wc, cModelConfPtr mc, PersistenceManagerPtr p) :
		IWave(wc, mc, p) {
	_realSpace = true;
}
void CBaseWave::InitializePropagators() {
	float_tt scale = _mc->d[2] * PI * _mc->wavelength;
	//t = exp(-i pi lam k^2 dz)

	auto s = scale * _k2;
	_prop = af::complex(af::cos(s), af::sin(s));

	af::array ck2max = af::constant(_k2max, _wc->n[0], _wc->n[1]);
	af::array bandwidthLimit = (_k2 < ck2max);
	_prop *= bandwidthLimit;
	_persist->Save2DDataSet(_prop, "Propagator");
}

af::array CBaseWave::fftShift(af::array& wave) {
	return af::shift(wave, wave.dims(0) / 2, wave.dims(1) / 2);
}

af::array CBaseWave::ifftShift(af::array& wave) {
	return af::shift(wave, (wave.dims(0) + 1) / 2, (wave.dims(1) + 1) / 2);
}

void CBaseWave::ApplyTransferFunction() {
	int px = GetTotalPixels();

	ToFourierSpace();
	for (int i = 0; i < _wc->n[0]; i++)
		for (int j = 0; j < _wc->n[1]; j++) {
		}
	ToRealSpace();
}
void CBaseWave::PropagateToNextSlice() {
	_wave_af = _wave_af * _prop;
} /* end propagate */

void CBaseWave::Transmit(af::array& t_af) {
	_wave_af *= t_af;
} /* end transmit() */

CBaseWave::~CBaseWave() {
}

void CBaseWave::InitializeKVectors() {
	float_tt ax = _mc->d[0] * _wc->n[0];
	float_tt by = _mc->d[1] * _wc->n[1];

//	printf("dx: %f , lambda:%g", _mc->d[0], _mc->wavelength);

	auto kx = (af::range(_wc->n[0]) - _wc->n[0] / 2) * (_mc->wavelength * 1e10) / ax;
	kx = af::shift(kx, kx.dims(0) / 2);
	_kx = af::tile(kx, 1, _wc->n[1]);

	auto ky = (af::range(_wc->n[1]) - _wc->n[1] / 2) * (_mc->wavelength * 1e10) / ax;
	ky = af::shift(ky, ky.dims(0) / 2);
	_ky = af::tile(ky.T(), _wc->n[0], 1);

	_k2max = _wc->n[0] / (2.0F * ax);
	if (_wc->n[1] / (2.0F * by) < _k2max)
		_k2max = _wc->n[1] / (2.0F * by);
	_k2max = 2.0 / 3.0 * _k2max;
	_k2max = _k2max * _k2max;

	_k2 = _kx * _kx + _ky * _ky;
	_persist->Save2DDataSet(_kx, "_kx");
	_persist->Save2DDataSet(_ky, "_ky");
	_persist->Save2DDataSet(_k2, "_k2");
}

void CBaseWave::FormProbe() {
	_wave.resize(boost::extents[_wc->n[0]][_wc->n[1]]);
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
	BOOST_LOG_TRIVIAL(info)<<format("* Real space res.:      %gA (=%gmrad)")% (1.0/_k2max)%(_mc->wavelength*_k2max*1000.0);
	BOOST_LOG_TRIVIAL(info)<<format("* Reciprocal space res: dkx=%g, dky=%g (1/A)")% (1.0/(_wc->n[0]*_mc->d[0]))%(1.0/(_wc->n[1]*_mc->d[1]));

	BOOST_LOG_TRIVIAL(info)<<format("* Beams:                %d x %d ")%_wc->n[0]%_wc->n[1];

	BOOST_LOG_TRIVIAL(info)<<format("* Acc. voltage:         %g kV (lambda=%gA)")%_mc->EnergykeV%_mc->wavelength;

	BOOST_LOG_TRIVIAL(info)<<format("* Probe array:          %d x %d pixels")%_wc->n[0]%_wc->n[1];
	BOOST_LOG_TRIVIAL(info)<<format("*                       %g x %gA")%(_wc->n[0]*_mc->d[0])%(_wc->n[1]*_mc->d[1]);
}

float_tt CBaseWave::GetIntegratedIntensity() const {
	float_tt intensity;
	intensity = af::sum<float_tt>(af::sqrt(af::real(_wave_af) * af::real(_wave_af) + af::imag(_wave_af) * af::imag(_wave_af)));
	// TODO: divide by px or not?
	return (intensity) / (_wc->n[0] * _wc->n[1]);
}

// FFT to Fourier space, but only if we're current in real space
void CBaseWave::ToFourierSpace() {
	if (_realSpace) {
		_realSpace = false;
		af::fft2InPlace(_wave_af);
	}
}

// FFT back to realspace, but only if we're currently in Fourier space
void CBaseWave::ToRealSpace() {
	if (!_realSpace) {
		_realSpace = true;
		af::ifft2InPlace(_wave_af);
	}
}

} //end namespace slicepp
