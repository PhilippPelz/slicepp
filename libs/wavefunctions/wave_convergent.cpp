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

#include "wave_convergent.hpp"
#include "boost/log/trivial.hpp"
using boost::format;

namespace QSTEM
{

CConvergentWave::CConvergentWave(const ConfigPtr& c,const PersistenceManagerPtr& p) : CBaseWave(c,p)
{
	// TODO: where does beam current belong?
	//configReader->ReadDoseParameters(m_beamCurrent, m_dwellTime);
	const WaveConfig& w = c->Wave;
	m_Cs = w.Cs;
	m_C5= w.C5;
	m_Cc= w.Cc;
	m_df0= w.Defocus;
	m_Scherzer= "";
	m_astigMag= w.Astigmatism;
	m_a33= w.a_33;
	m_a31= w.a_31;
	m_a44= w.a_44;
	m_a42= w.a_42;
	m_a55= w.a_55;
	m_a53= w.a_53;
	m_a51= w.a_51;
	m_a66= w.a_66;
	m_a64= w.a_64;
	m_a62= w.a_62;
	m_astigAngle= w.AstigmatismAngle;
	m_phi33= w.phi_33;
	m_phi31= w.phi_31;
	m_phi44= w.phi_44;
	m_phi42= w.phi_42;
	m_phi55= w.phi_55;
	m_phi53= w.phi_53;
	m_phi51= w.phi_51;
	m_phi66= w.phi_66;
	m_phi64= w.phi_64;
	m_phi62= w.phi_62;
	_smoothen= w.Smooth;
	_gaussScale= w.gaussScale;
	_isGaussian= w.Gaussian;
	m_dE_E= w.dE_E;
	m_dI_I= w.dI_I;
	m_dV_V= w.dV_V;
	m_alpha= w.alpha;
	_CLA= w.AISaperture;
	m_printLevel = c->Output.LogLevel;
}

/** Copy constructor - used to copy wave just before dispatching multiple threads for STEM simulations */
CConvergentWave::CConvergentWave(const CConvergentWave& other) : CBaseWave(other)
{
	// TODO: need to copy arrays and anything pointed to - anything that needs to be thread-local
}


WavePtr CConvergentWave::Clone()
{
	return WavePtr(new CConvergentWave(*this));
}

void CConvergentWave::DisplayParams()
{
	CBaseWave::DisplayParams();

	BOOST_LOG_TRIVIAL(info)<< format("* Aperture half angle:  %g mrad") % m_alpha;
	BOOST_LOG_TRIVIAL(info)<< format("* AIS aperture:");
	if (_CLA > 0)
		BOOST_LOG_TRIVIAL(info)<< format("*                  %g A") % (_CLA);
	else BOOST_LOG_TRIVIAL(info)<< format("                        none");

	BOOST_LOG_TRIVIAL(info)<< format("* C_3 (C_s):            %g mm") % (m_Cs*1e-7);
	BOOST_LOG_TRIVIAL(info)<< format("* C_1 (Defocus):        %g nm%s") % (0.1*m_df0) % ((m_Scherzer == "Scherzer") ? " (Scherzer)" : (m_Scherzer=="") ? " (opt.)":"");
	BOOST_LOG_TRIVIAL(info)<< format("* Astigmatism:          %g nm, %g deg") % (0.1*m_astigMag) % (RAD2DEG*m_astigAngle);

	if (m_a33 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_3,3:                %g nm, phi=%g deg") %(m_a33*1e-1)%(m_phi33*RAD2DEG);
	if (m_a31 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_3,1:                %g nm, phi=%g deg") %(m_a31*1e-1)%(m_phi31*RAD2DEG);

	if (m_a44 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_4,4:                %g um, phi=%g deg") %(m_a44*1e-4)%(m_phi44*RAD2DEG);
	if (m_a42 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_4,2:                %g um, phi=%g deg") %(m_a42*1e-4)%(m_phi42*RAD2DEG);

	if (m_a55 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_5,5:                %g um, phi=%g deg") %(m_a55*1e-4)%(m_phi55*RAD2DEG);
	if (m_a53 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_5,3:                %g um, phi=%g deg") %(m_a53*1e-4)%(m_phi53*RAD2DEG);
	if (m_a51 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_5,1:                %g um, phi=%g deg") %(m_a51*1e-4)%(m_phi51*RAD2DEG);

	if (m_a66 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_6,6:                %g um, phi=%g deg") %(m_a66*1e-7)%(m_phi66*RAD2DEG);
	if (m_a64 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_6,4:                %g um, phi=%g deg") %(m_a64*1e-7)%(m_phi64*RAD2DEG);
	if (m_a62 > 0)BOOST_LOG_TRIVIAL(info)<< format("* a_6,2:                %g um, phi=%g deg") %(m_a62*1e-7)%(m_phi62*RAD2DEG);
	if (m_C5 != 0)BOOST_LOG_TRIVIAL(info)<< format("* C_5:                  %g mm") % (m_C5*1e-7);

	BOOST_LOG_TRIVIAL(info)<< format("* C_c:                  %g mm") %(m_Cc*1e-7);
	BOOST_LOG_TRIVIAL(info)<< format("* Damping dE/E: %g / %g ") % (sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I)*m_v0*1e3) %(m_v0*1e3);

	/*
    // TODO: where does beam current belong?
  printf("* beam current:         %g pA\n",m_beamCurrent);
  printf("* dwell time:           %g msec (%g electrons)\n",
         m_dwellTime,m_electronScale);
	 */
}

/**********************************************
 * This function creates a incident STEM probe
 * at position (dx,dy)
 * with parameters given in muls
 *
 * The following Abberation functions are being used:
 * 1) ddf = Cc*dE/E + Cc2*(dE/E)^2,
 *    Cc, Cc2 = chrom. Abber. (1st, 2nd order) [1]
 * 2) chi(qx,qy) = (2*pi/lambda)*{0.5*C1*(qx^2+qy^2)+
 *                 0.5*C12a*(qx^2-qy^2)+
 *                 C12b*qx*qy+
 *                 C21a/3*qx*(qx^2+qy^2)+
 *                 ...
 *                 +0.5*C3*(qx^2+qy^2)^2
 *                 +0.125*C5*(qx^2+qy^2)^3
 *                 ... (need to finish)
 *
 *
 *    qx = acos(kx/K), qy = acos(ky/K)
 *
 * References:
 * [1] J. Zach, M. Haider,
 *    "Correction of spherical and Chromatic Abberation
 *     in a low Voltage SEM", Optik 98 (3), 112-118 (1995)
 * [2] O.L. Krivanek, N. Delby, A.R. Lupini,
 *    "Towards sub-Angstroem Electron Beams",
 *    Ultramicroscopy 78, 1-11 (1999)
 *
 *********************************************/
#define SMOOTH_EDGE 5 // make a smooth edge on AIS aperture over +/-SMOOTH_EDGE pixels
void CConvergentWave::FormProbe()
{
	CBaseWave::FormProbe();
	ConstructWave();
}  /* end probe() */


void CConvergentWave::ConstructWave(){
	unsigned iy, ixmid, iymid;
	float_tt rmin, rmax, aimin, aimax;
	float_tt k2max, alpha;
	float_tt ax = _nx*m_dx;
	float_tt by = _ny*m_dy;
	float_tt dx = ax-_nx/2*m_dx;
	float_tt dy = by-_ny/2*m_dy;
	float_tt avgRes = sqrt(0.5*(m_dx*m_dx+m_dy*m_dy));
	float_tt edge = SMOOTH_EDGE*avgRes;

	float_tt sum = 0.0;
	bool condition;
	af::array kx, ky, kx2, ky2, k2, k2max_af, ktheta, ktheta2, phi, chi, real, imag, scale;
	af::array x, y, pixel, rx, ry, rx2, ry2;
	af::array delta;
	/********************************************************
	 * formulas from:
	 * http://cimesg1.epfl.ch/CIOL/asu94/ICT_8.html
	 *
	 * dE_E = dE/E = energy spread of emitted electrons
	 * dV_V = dV/V = acc. voltage fluctuations
	 * dI_I = dI/I = lens current fluctuations
	 * delta defocus in Angstroem (Cc in A)
	 *******************************************************/
	float_tt delta_c = m_Cc*m_dE_E;
	if (m_printLevel > 2) printf("defocus offset: %g nm (Cc = %g)\n",delta_c,m_Cc);
	delta = af::constant(delta_c, _nx, _ny);

	/**********************************************************
	 *  Calculate misc constants
	 *********************************************************/
//
//	float_tt rx = 1.0/ax;
//	float_tt rx2 = rx * rx;
//	float_tt ry = 1.0/by;
//	float_tt ry2 = ry * ry;
	rx = af::constant(1.0/ax, _nx, _ny);
	rx2 = rx*rx;
	ry = af::constant(1.0/by, _nx, _ny);
	ry2 = ry*ry;

	ixmid = _nx/2;
	iymid = _ny/2;

	// Start in Fourier space
	ToFourierSpace();
	/* convert convergence angle from mrad to rad */
	alpha = 0.001*m_alpha;
	k2max = sin(alpha)/m_wavlen;  /* = K0*sin(alpha) */
	k2max = k2max * k2max;

	/*   Calculate MTF
       NOTE zero freg is in the bottom left corner and
       expandes into all other corners - not in the center
       this is required for FFT

       PIXEL = diagonal width of pixel squared
       if a pixel is on the apertur boundary give it a weight
       of 1/2 otherwise 1 or 0
	 */
	pixel = ( rx2 + ry2 );
	scale = af::constant(1.0/sqrt((float_tt)_nx*(float_tt)_ny), _nx, _ny);

	kx = af::range(_nx);
	ky = af::range(_ny);
	kx(af::seq(ixmid + 1, af::end)) -= _nx;
	kx = af::tile(kx, 1, _ny);
	ky(af::seq(iymid + 1, af::end)) -= _ny;
	ky = af::tile(ky.T(), _nx);
	ky2 = ky*ky*ry2;
	kx2 = kx*kx*rx2;
	k2 = kx2 + ky2;
	ktheta2 = k2*(m_wavlen*m_wavlen);
	ktheta = af::sqrt(ktheta2);
	phi = af::atan2(ry*ky, rx*kx);
	chi = ktheta2*(m_df0 +delta + m_astigMag*af::cos(2.0*(phi-m_astigAngle)))/2.0;

	ktheta2 *= ktheta;  //ktheta^3
	if ((m_a33 > 0) || (m_a31 > 0))
		chi += ktheta2*(m_a33*af::cos(3.0*(phi-m_phi33)) + m_a31*af::cos(phi-m_phi31))/3.0;

	ktheta2 *= ktheta;  //ktheta^4
	if ((m_a44 > 0) || (m_a42 > 0) || (m_Cs != 0))
		chi += ktheta2*(m_a44*af::cos(4.0*(phi-m_phi44)) + m_a42*af::cos(2.0*(phi-m_phi42)) + m_Cs)/4.0;

	ktheta2 *= ktheta;  //ktheta^5
	if ((m_a55 > 0) || (m_a53 > 0) || (m_a51 > 0))
		chi += ktheta2*(m_a55*af::cos(5.0*(phi-m_phi55)) + m_a53*af::cos(3.0*(phi-m_phi53)) + m_a51*af::cos(phi-m_phi51))/5.0;

	ktheta2 *= ktheta;  //ktheta^6
	if ((m_a66 > 0) || (m_a64 > 0) || (m_a62 = 0) || (m_C5 != 0))
		chi += ktheta2*(m_a66*af::cos(6.0*(phi-m_phi66))+m_a64*af::cos(4.0*(phi-m_phi64))+m_a62*af::cos(2.0*(phi-m_phi62))+m_C5)/6.0;

	af::array condition2;

	chi *= 2*PI/m_wavlen;
	chi -= 2.0*PI*( (kx*dx/ax) + (ky*dy/by) );
	k2max_af = af::constant(k2max, _nx, _ny);
	real = af::constant(0, _nx, _ny);
	imag = af::constant(0, _nx, _ny);
	condition2 = (k2 <= k2max_af);
	real += condition2*scale * af::cos(chi);
	imag += condition2*scale * af::sin(chi);

	condition2 = ( _smoothen != 0)*( af::abs(k2-k2max_af) <= pixel);
	real -= af::constant(0.5, _nx, _ny)*condition2*scale * af::cos(chi);
	imag -= af::constant(1.5, _nx, _ny)*condition2*scale * af::sin(chi);
	_wave_af = af::complex(real, imag);

	/* Fourier transform into real space */
	ToRealSpace();
	/**********************************************************
	 * display cross section of probe intensity
	 */
	/* multiply with gaussian in Real Space in order to avoid artifacts */
	af::array r;
	if (_isGaussian) {
		af::array kx, ky;
		kx = (af::range(_nx) - _nx/2);
		ky = (af::range(_ny) - _ny/2);
		kx = af::tile(kx, 1, _ny);
		ky = af::tile(ky.T(), _nx);
		r = af::exp((-1)*(kx*kx + ky*ky)/(_nx*_ny*_gaussScale));
		_wave_af *= r;
	}
	/* Apply AIS aperture in Real Space */
	// printf("center: %g,%g\n",dx,dy);
	if (_CLA > 0) {
		x = af::range(_nx)*m_dx-dx;
		y = af::range(_ny)*m_dy-dy;
		x = af::tile(x, 1, _ny);
		y = af::tile(y.T(), _nx);
		r = af::sqrt(x*x+y*y);
		delta = r - af::constant(0.5*_CLA, _nx, _ny);
		delta += af::constant(edge, _nx, _ny);
		condition2 = !(delta >_zero);
		_wave_af *= condition;
		condition2 = !(delta < af::constant(-edge, _nx, _ny));
		scale = (!condition* scale) + (condition*0.5*(af::constant(1, _nx, _ny)-cos(delta*PI/edge)));
		_wave_af = (!condition* _wave_af) + (condition*scale*_wave_af);
	}
	/*  Normalize probe intensity to unity  */



 	sum = af::sum<float_tt>(af::real(_wave_af)*af::real(_wave_af) + af::imag(_wave_af)*af::imag(_wave_af));
	float_tt scale_s = 1.0 / sum;
	_wave_af *= (float_tt)sqrt(scale_s);

	//sum = af::sum<float_tt>(af::sqrt(af::real(_wave_af)*af::real(_wave_af) + af::imag(_wave_af)*af::imag(_wave_af)));
	sum = 0.0;
	real = af::real(_wave_af);
	imag = af::imag(_wave_af);
	rmin = af::min<float_tt>(real);
	rmax = af::max<float_tt>(real);
	aimin = af::min<float_tt>(imag);
	aimax = af::max<float_tt>(imag);

	m_rmin = rmin;
	m_rmax = rmax;
	m_aimin = aimin;
	m_aimax = aimax;
	_probe = af::array(_wave_af);

	BOOST_LOG_TRIVIAL(info) << format("wave value range (%f .. %f,i %f ... %f)") % rmin % rmax % aimin % aimax;
	/**********************************************************/
}

} // end namespace QSTEM
