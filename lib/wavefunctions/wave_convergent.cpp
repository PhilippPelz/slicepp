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
using boost::format;

namespace slicepp {

CConvergentWave::CConvergentWave(cWaveConfPtr wc, cModelConfPtr mc, PersistenceManagerPtr p) :
		CBaseWave(wc, mc, p) {
	m_Cs = wc->Cs;
	m_C5 = wc->C5;
	m_Cc = wc->Cc;
	m_df0 = wc->Defocus;
	m_Scherzer = "";
	m_astigMag = wc->Astigmatism;
	m_a33 = wc->a_33;
	m_a31 = wc->a_31;
	m_a44 = wc->a_44;
	m_a42 = wc->a_42;
	m_a55 = wc->a_55;
	m_a53 = wc->a_53;
	m_a51 = wc->a_51;
	m_a66 = wc->a_66;
	m_a64 = wc->a_64;
	m_a62 = wc->a_62;
	m_astigAngle = wc->AstigmatismAngle;
	m_phi33 = wc->phi_33;
	m_phi31 = wc->phi_31;
	m_phi44 = wc->phi_44;
	m_phi42 = wc->phi_42;
	m_phi55 = wc->phi_55;
	m_phi53 = wc->phi_53;
	m_phi51 = wc->phi_51;
	m_phi66 = wc->phi_66;
	m_phi64 = wc->phi_64;
	m_phi62 = wc->phi_62;
	_smoothen = wc->IsSmooth;
	_sigma = wc->gaussFWHM / 2.35482;
	_isGaussian = wc->IsGaussian;
	m_dE_E = wc->dE_E;
	m_dI_I = wc->dI_I;
	m_dV_V = wc->dV_V;
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

	BOOST_LOG_TRIVIAL(info)<< format("* C_3 (C_s):            %g mm") % (m_Cs*1e-7);
	BOOST_LOG_TRIVIAL(info)<< format("* C_1 (Defocus):        %g nm%s") % (0.1*m_df0) % ((m_Scherzer == "Scherzer") ? " (Scherzer)" : (m_Scherzer=="") ? " (opt.)":"");
	BOOST_LOG_TRIVIAL(info)<< format("* Astigmatism:          %g nm, %g deg") % (0.1*m_astigMag) % (RAD2DEG*m_astigAngle);

	if (m_a33 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_3,3:                %g nm, phi=%g deg") %(m_a33*1e-1)%(m_phi33*RAD2DEG);
	if (m_a31 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_3,1:                %g nm, phi=%g deg") %(m_a31*1e-1)%(m_phi31*RAD2DEG);

	if (m_a44 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_4,4:                %g um, phi=%g deg") %(m_a44*1e-4)%(m_phi44*RAD2DEG);
	if (m_a42 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_4,2:                %g um, phi=%g deg") %(m_a42*1e-4)%(m_phi42*RAD2DEG);

	if (m_a55 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_5,5:                %g um, phi=%g deg") %(m_a55*1e-4)%(m_phi55*RAD2DEG);
	if (m_a53 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_5,3:                %g um, phi=%g deg") %(m_a53*1e-4)%(m_phi53*RAD2DEG);
	if (m_a51 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_5,1:                %g um, phi=%g deg") %(m_a51*1e-4)%(m_phi51*RAD2DEG);

	if (m_a66 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_6,6:                %g um, phi=%g deg") %(m_a66*1e-7)%(m_phi66*RAD2DEG);
	if (m_a64 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_6,4:                %g um, phi=%g deg") %(m_a64*1e-7)%(m_phi64*RAD2DEG);
	if (m_a62 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_6,2:                %g um, phi=%g deg") %(m_a62*1e-7)%(m_phi62*RAD2DEG);
	if (m_C5 != 0)
		BOOST_LOG_TRIVIAL(info)<< format("* C_5:                  %g mm") % (m_C5*1e-7);

	BOOST_LOG_TRIVIAL(info)<< format("* C_c:                  %g mm") %(m_Cc*1e-7);
	BOOST_LOG_TRIVIAL(info)<< format("* Damping dE/E: %g / %g ") % (sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I)*_mc->EnergykeV*1e3) %(_mc->EnergykeV*1e3);
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
	auto k = af::sqrt(_k2) * lambdaA;
	auto phi = af::atan2(_ky * lambdaA, _kx * lambdaA);
	auto alpha = af::asin(k);
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
	auto chi = zeros;
	chi += 0.5 * af::pow(k, 2) * (m_df0 + m_Cc * m_dE_E + (m_astigMag * af::cos(2.0 * (phi - m_astigAngle))));
	if ((m_a33 > 0) || (m_a31 > 0))
		chi += af::pow(k, 3) * (m_a33 * af::cos(3.0 * (phi - m_phi33)) + m_a31 * af::cos(phi - m_phi31)) / 3.0;
	if ((m_a44 > 0) || (m_a42 > 0) || (m_Cs != 0))
		chi += af::pow(k, 4) * (m_a44 * af::cos(4.0 * (phi - m_phi44)) + m_a42 * af::cos(2.0 * (phi - m_phi42)) + m_Cs) / 4.0;
	if ((m_a55 > 0) || (m_a53 > 0) || (m_a51 > 0))
		chi += af::pow(k, 5) * (m_a55 * af::cos(5.0 * (phi - m_phi55)) + m_a53 * af::cos(3.0 * (phi - m_phi53)) + m_a51 * af::cos(phi - m_phi51))
				/ 5.0;
	if ((m_a66 > 0) || (m_a64 > 0) || (m_a62 = 0) || (m_C5 != 0))
		chi += af::pow(k, 6)
				* (m_a66 * af::cos(6.0 * (phi - m_phi66)) + m_a64 * af::cos(4.0 * (phi - m_phi64)) + m_a62 * af::cos(2.0 * (phi - m_phi62)) + m_C5)
				/ 6.0;

	chi *= 2 * PI / lambdaA;

	_wave_af = weights * af::complex(af::cos(chi), -af::sin(chi));

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
		auto gauss = af::complex(r,0);
		_wave_af *= r;
	}

	float_tt sum = af::sum<float_tt>(af::real(_wave_af) * af::real(_wave_af) + af::imag(_wave_af) * af::imag(_wave_af));
	float_tt scale_s = 1.0 / sum;
	_wave_af *= (float_tt) sqrt(scale_s);

	sum = 0.0;
	auto real = af::real(_wave_af);
	auto imag = af::imag(_wave_af);
	float_tt rmin = af::min<float_tt>(real);
	float_tt rmax = af::max<float_tt>(real);
	float_tt aimin = af::min<float_tt>(imag);
	float_tt aimax = af::max<float_tt>(imag);

	m_rmin = rmin;
	m_rmax = rmax;
	m_aimin = aimin;
	m_aimax = aimax;
	_probe = _wave_af.copy();

	BOOST_LOG_TRIVIAL(info)<<format("wave value range (%f .. %f,i %f ... %f)") % rmin % rmax % aimin % aimax;
}

}
 // end namespace slicepp
