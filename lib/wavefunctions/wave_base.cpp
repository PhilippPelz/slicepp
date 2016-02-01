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

CBaseWave::CBaseWave(cWaveConfPtr wc, cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) :
		_OLaberrations(mc->wavelength*1e10,mc->OLaberrations,"Objective Lens Aberrations"),
		IWave(wc, mc, oc,p) {
	_realSpace = true;
}
void CBaseWave::InitializePropagators() {

}
void CBaseWave::PropagateToNextSlice() {
	_wave_af *= _prop;
}
void CBaseWave::Transmit(af::array& t_af) {
	_wave_af *= t_af;
}
void CBaseWave::ApplyCTF(){
	auto phaseplate = _OLaberrations.getPhasePlate(_kabs,_kx,_ky,_persist);
	ToFourierSpace();
	_wave_af *= phaseplate;
	ToRealSpace();
}
CBaseWave::~CBaseWave() {
}
void CBaseWave::FormProbe() {
	_wave.resize(boost::extents[_wc->n[0]][_wc->n[1]]);

	float_tt ax = _mc->d[0] * _wc->n[0];
	float_tt by = _mc->d[1] * _wc->n[1];

	auto kx = (af::range(_wc->n[0]) - _wc->n[0] / 2) / ax;
	kx = af::shift(kx, kx.dims(0) / 2);
	_kx = af::tile(kx, 1, _wc->n[1]);

	auto ky = (af::range(_wc->n[1]) - _wc->n[1] / 2) / ax;
	ky = af::shift(ky, ky.dims(0) / 2);
	_ky = af::tile(ky.T(), _wc->n[0], 1);

	_k2max = _wc->n[0] / (2.0F * ax);
	if (_wc->n[1] / (2.0F * by) < _k2max)
		_k2max = _wc->n[1] / (2.0F * by);
	_k2max = 2.0 / 3.0 * _k2max;
	_k2max = _k2max * _k2max;

	auto _k2 = _kx * _kx + _ky * _ky;

	float_tt scale = - _mc->d[2] * PI * (_mc->wavelength*1e10);
	//t = exp(-i pi lam k^2 dz)

	auto s = scale * _k2;
	_prop = af::complex(af::cos(s), af::sin(s));

	af::array bandwidthLimit = (_k2 > _k2max);
	_prop(bandwidthLimit) = 0;

	if(_oc->SaveProbe)
		_persist->Save2DDataSet(_prop, "propagator");

	_kabs = af::sqrt(_k2);
}

void CBaseWave::ResetProbe() {
	_wave_af = _probe.copy();
	_realSpace = true;
}

void CBaseWave::DisplayParams() {
	BOOST_LOG_TRIVIAL(info)<<
	"*****************************  Wave  Parameters **************************************************";
	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info)<<format("* Real space res.:      %gA (=%gmrad)")% (1.0/_k2max)%(_mc->wavelength*_k2max*1000.0);
	BOOST_LOG_TRIVIAL(info)<<format("* Reciprocal space res: dkx=%g, dky=%g (1/A)")% (1.0/(_wc->n[0]*_mc->d[0]))%(1.0/(_wc->n[1]*_mc->d[1]));
	BOOST_LOG_TRIVIAL(info)<<format("* Acc. voltage:         %g kV (lambda=%gA)")%_mc->EnergykeV%_mc->wavelength;
	BOOST_LOG_TRIVIAL(info)<<format("* Probe array:          %d x %d pixels")%_wc->n[0]%_wc->n[1];
	BOOST_LOG_TRIVIAL(info)<<format("*                       %g x %gA")%(_wc->n[0]*_mc->d[0])%(_wc->n[1]*_mc->d[1]);
	BOOST_LOG_TRIVIAL(info) <<
		"**************************************************************************************************";
	_OLaberrations.DisplayParams();
}

float_tt CBaseWave::GetIntegratedIntensity() const {
	return af::sum<float_tt>(GetIntensity());
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
