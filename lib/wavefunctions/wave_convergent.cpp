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
 along with this program.  If not, see <http://wwwc->gnu.org/licenses/>.
 */

#include "wave_convergent.hpp"
#include "boost/log/trivial.hpp"
#include "afhelpers.hpp"
using boost::format;

namespace slicepp {

CConvergentWave::CConvergentWave(cWaveConfPtr wc, cModelConfPtr mc, cOutputConfPtr oc, PersistenceManagerPtr p) :
		_aberration(mc->wavelength*1e10,wc->aberrations,"Convergent Beam"),
		CBaseWave(wc, mc,oc, p) {
	_smoothen = wc->IsSmooth;
	_sigma = wc->gaussFWHM / 2.35482;
	_isGaussian = wc->IsGaussian;
	_alpha_max = wc->alpha * 0.001;
	_CLA = wc->AISaperture;
}

WavePtr CConvergentWave::Clone() {
	return WavePtr(new CConvergentWave(*this));
}

void CConvergentWave::DisplayParams() {
	CBaseWave::DisplayParams();

	BOOST_LOG_TRIVIAL(info)<< format("* Aperture half angle:  %g mrad") % (_alpha_max/1000);
	BOOST_LOG_TRIVIAL(info)<< format("* AIS aperture:");
	if (_CLA > 0)
		BOOST_LOG_TRIVIAL(info)<< format("*                  %g A") % (_CLA);
		else BOOST_LOG_TRIVIAL(info)<< format("                        none");
	_aberration.DisplayParams();
}

void CConvergentWave::FormProbe() {
	CBaseWave::FormProbe();
	ConstructWave();
}

void CConvergentWave::ConstructWave() {
	_realSpace = false;
	auto ones = af::constant(1, _wc->n[0], _wc->n[1], f32);
	auto zeros = af::constant(0, _wc->n[0], _wc->n[1], f32);
	auto weights = af::complex(ones, zeros);
	float lambdaA = 1e10 * _mc->wavelength;

	auto kmax = sin(_alpha_max) / lambdaA;
	auto alpha = af::asin(_kabs * lambdaA);
	weights(alpha > _alpha_max) = 0;

	if (_smoothen) {
		auto dk = min(1 / (_wc->n[0] * _mc->d[0]), 1 / (_wc->n[1] * _mc->d[1]));
		auto dEdge = 2 / (kmax / dk);
		auto alpha_norm = alpha / _alpha_max;
		auto ind = (alpha_norm > (1 - dEdge)) && (alpha_norm < (1 + dEdge));
		auto tmp = (0.5 * (1 - af::sin(PI / (2 * dEdge) * (alpha_norm - 1))));
		tmp = af::complex(tmp,zeros);
		weights(ind) = tmp(ind);
	}
	_persist->Save2DDataSet(weights,"weights");

	auto phase = _aberration.getPhasePlate(_kabs,_kx,_ky,_persist);
	_persist->Save2DDataSet(phase,"phaseplate");
	_wave_af = weights * phase;

	/* Fourier transform into real space */
	ToRealSpace();
	_wave_af = fftShift(_wave_af);
	if (_isGaussian) {
		auto ry = (af::range(_wc->n[1]) - _wc->n[1] / 2);
//		ry = af::shift(ry, ry.dims(0) / 2);
		ry = af::tile(ry.T(), _wc->n[0], 1);
		auto rx = (af::range(_wc->n[0]) - _wc->n[0] / 2);
//		rx = af::shift(rx, rx.dims(0) / 2);
		rx = af::tile(rx.T(), _wc->n[1], 1);
		auto r2 = rx*rx+ry*ry;

		auto r = af::exp(-r2/2/(_sigma*_sigma));
		_wave_af *= r.as(c32);
	}

	float_tt sum = af::sum<float_tt>(GetIntensity());
	float_tt scale_s = 1.0 / sum;
	_wave_af *= (float_tt) sqrt(scale_s);

	auto real = af::real(_wave_af);
	auto imag = af::imag(_wave_af);
	float_tt rmin = af::min<float_tt>(real);
	float_tt rmax = af::max<float_tt>(real);
	float_tt aimin = af::min<float_tt>(imag);
	float_tt aimax = af::max<float_tt>(imag);

	_probe = _wave_af.copy();

	BOOST_LOG_TRIVIAL(info)<<format("wave value range (%f .. %f,i %f ... %f)") % rmin % rmax % aimin % aimax;
}

}
 // end namespace slicepp
