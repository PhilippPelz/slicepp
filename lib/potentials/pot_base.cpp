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

#include "pot_base.hpp"
#include "matrixlib.hpp"

#include <omp.h>
#include <algorithm>
#include <boost/format.hpp>

using boost::format;

const int BUF_LEN = 256;

namespace slicepp {

CPotential::CPotential(cModelConfPtr mc, cOutputConfPtr oc , PersistenceManagerPtr persist) :
		_ddx(0), _ddy(0) ,_nrAtomTrans(0),
		_dkx(0), _dky(0), _dkz(0), _kmax(0), _kmax2(0), _totalThickness(0),
		_dr(0), _atomRadius2(0), _offsetX(0), _offsetY(0), _nRadX(0),
		_nRadY(0), _nRadZ(0), _nRad2Trans(0), _ndiaAtomX(0), _ndiaAtomY(0), _boxNx(0),
		_boxNy(0), m_boxNz(0),
		IPotential(mc,oc,persist) {

}
CPotential::~CPotential() {
}

void CPotential::DisplayParams() {
	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info)<<format("* Log level:            %d") % _oc->LogLevel;
	BOOST_LOG_TRIVIAL(info)<<format("* Model Sampling:       %g x %g x %g A") % _mc->d[0]% _mc->d[1]% _mc->d[2];
//	BOOST_LOG_TRIVIAL(info)<<format("* Pot. array offset:    (%g,%g,%g) A") % _c->Structure->xOffset%_c->Structure->yOffset% _c->Structure->zOffset;
	BOOST_LOG_TRIVIAL(info)<<format("* Potential periodic:   (x,y): %s ") %((_mc->periodicXY) ? "yes" : "no");
	BOOST_LOG_TRIVIAL(info)<<format("* Potential array:      %d x %d")% _mc->n[0]% _mc->n[1];
	BOOST_LOG_TRIVIAL(info)<<format("*                       %g x %gA") % (_mc->n[0] * _mc->d[0])% (_mc->n[1] * _mc->d[1]);
	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************";
}

void CPotential::ReadPotential(std::string &fileName, unsigned subSlabIdx) {
}

void CPotential::ReadSlice(const std::string &fileName,
		ComplexArray2DView slice, unsigned idx) {
}

void CPotential::SliceSetup() {
	_ddx = _mc->d[0] / (double) OVERSAMPLING;
	_ddy = _mc->d[1] / (double) OVERSAMPLING;

	_ndiaAtomX = 2 * OVERSAMPLING * (int) ceil(_mc->ratom / _mc->d[0] );
	_ndiaAtomY = 2 * OVERSAMPLING * (int) ceil(_mc->ratom / _mc->d[1] );

	_dkx = 0.5 * OVERSAMPLING / ((_ndiaAtomX) * _mc->d[0]);
	_dky = 0.5 * OVERSAMPLING / ((_ndiaAtomY) * _mc->d[1]);

	_kmax = 0.5 * _ndiaAtomX * _dkx / (double) OVERSAMPLING; // largest k that we'll admit
	_kmax2 = _kmax * _kmax;

	/* For now we don't care, if the box has only small prime factors, because we will not fourier transform it  especially not very often. */
	_boxNx = (int) (_mc->ratom / _ddx + 2.0);
	_boxNy = (int) (_mc->ratom / _ddy + 2.0);

	_totalThickness = _mc->d[2] * _mc->n[2];
	_dr = min(_mc->d[0],_mc->d[1]) / OVERSAMPLING;
	_nrAtomTrans =  (int) ceil(_mc->ratom / _dr + 0.5);

	_nRadX = (int) ceil(_mc->ratom / _mc->d[0]);
	_nRadY = (int) ceil(_mc->ratom / _mc->d[1]);
	_nRadZ = (int) ceil(_mc->ratom / _mc->d[2]) ;
	_nRad2Trans = _nRadX * _nRadX + _nRadY * _nRadY;

	_atomRadius2 = _mc->ratom * _mc->ratom;
	_sliceThicknesses.resize(_mc->n[2]);

	_sliceThicknesses[0] = _mc->d[2];

	_t.resize(boost::extents[_mc->n[2]][_mc->n[0]][_mc->n[1]]);
	std::fill(_t.data(), _t.data() + _t.size(), complex_tt(0, 0));
}

complex_tt CPotential::GetSlicePixel(unsigned iz, unsigned ix, unsigned iy) {
	return _t[iz][ix][iy];
}

void CPotential::SetScatteringFactors(float_tt kmax) {
	scatPar[0][N_SF - 1] = 1.2 * kmax;
	scatPar[0][N_SF - 2] = 1.1 * kmax;
	scatPar[0][N_SF - 3] = kmax;

	if (scatPar[0][N_SF - 4] > scatPar[0][N_SF - 3]) {
		int ix = 0;
		if (1) {	// TODO: why if(1)? dont use all factors?
			// set additional scattering parameters to zero:
			for (; ix < N_SF - 10; ix++) {
				if (scatPar[0][N_SF - 4 - ix]
						< scatPar[0][N_SF - 3] - 0.001 * (ix + 1))
					break;
				scatPar[0][N_SF - 4 - ix] = scatPar[0][N_SF - 3]
						- 0.001 * (ix + 1);
				for (unsigned iy = 1; iy < N_ELEM; iy++)
					scatPar[iy][N_SF - 4 - ix] = 0;
			}
		} else {
			for (; ix < 20; ix++) {
				if (scatPar[0][N_SF - 4 - ix] < scatPar[0][N_SF - 3])
					break;
				scatPar[0][N_SF - 4 - ix] = scatPar[0][N_SF - 3];
				for (unsigned iy = 1; iy < N_ELEM; iy++)
					scatPar[iy][N_SF - 4 - ix] = scatPar[iy][N_SF - 3];
			}
		}
		if (_oc->LogLevel < 2)
				BOOST_LOG_TRIVIAL(info)<< format("* Reduced angular range of scattering factor to %g/Å!") % scatParOffs[0][N_SF - 4 - ix];
		}
	}

void CPotential::CenterAtomZ(atom& atom, float_tt &z) {
}

void CPotential::MakeSlices(superCellBoxPtr info) {
	time_t time0, time1;
	SliceSetup();

	for (std::vector<int>::iterator a = info->uniqueZ.begin(); a < info->uniqueZ.end(); a = a + 1) {
		ComputeAtomPotential(*a);
		if (_oc->SaveAtomicPotential)
			SaveAtomicPotential(*a);
	}

	time(&time0);
	int atomsAdded = 0;

	BOOST_LOG_TRIVIAL(info)<< "Adding atoms to slices ...";

#pragma omp parallel for shared(atomsAdded,info)
	for (std::vector<atom>::iterator a = info->atoms.begin(); a < info->atoms.end(); a = a + 1) {
		atom atom(a);
		if (atom.Znum == 0) continue;
		AddAtomToSlices(atom);
		BOOST_LOG_TRIVIAL(trace)<< format("Adding atom : (%3.3f, %3.3f, %3.3f) Z=%d") % atom.r[0] % atom.r[1] % atom.r[2] % atom.Znum;

		atomsAdded++;

		int interval = (info->atoms.size() / 20) == 0 ? 1 : (info->atoms.size() / 20);
		if ((atomsAdded % interval) == 0 && omp_get_thread_num() == 1)
			loadbar(atomsAdded + 1, info->atoms.size());
	}

	CleanUp();
	_t_d = af::array(_mc->n[0], _mc->n[1], _mc->n[2], (afcfloat *)_t.data());
	MakePhaseGratings();
	time(&time1);
	BOOST_LOG_TRIVIAL(info)<< format( "%g sec used for real space potential calculation (%g sec per atom)")
	% difftime(time1, time0)%( difftime(time1, time0) / info->atoms.size());

	SavePotential();
}

void CPotential::MakePhaseGratings() {
	float_tt mm0 = 1.0F + _mc->EnergykeV / 511.0F; // relativistic corr. factor gamma
	float_tt scale = mm0 * _mc->wavelength;

	BOOST_LOG_TRIVIAL(info)<<format("Making phase gratings for %d layers (scale=%g rad/VA, gamma=%g, sigma=%g) ... ")
	%_t.shape()[0]%scale%mm0%_mc->sigma;

	float_tt minph = 3.1, maxph = 0, minabs = 100, maxabs = 0;
		_t_d = scale * af::real(_t_d);
		maxph = af::max<float_tt>(_t_d);
		minph = af::min<float_tt>(_t_d);
		_t_d = af::complex(af::cos(_t_d), af::sin(_t_d));
	BOOST_LOG_TRIVIAL(info)<<format("Phase values %g ... %g")%minph%maxph;
}

void CPotential::WriteSlice(unsigned idx, std::string prefix) {
	char buf[255];
	std::map<std::string, float_tt> params;
	params["Thickness"] = _mc->d[2];
	params["dx"] = _mc->d[0];
	params["dy"] = _mc->d[1];
	sprintf(buf, "Projected Potential (slice %d)", idx);
	std::string comment = buf;
	std::stringstream filename;
	filename << prefix << idx;
	// TODO use persist
	//	m_imageIO->WriteImage(m_trans1[boost::indices[idx][range(0,_config->Model.nx)][range(0,_config->Model.ny)]], filename.str().c_str(), params, comment);
}

void CPotential::WriteProjectedPotential() {
	std::map<std::string, float_tt> params;
	char buf[255];
	ComplexArray2D sum(extents[_mc->n[0]][_mc->n[1]]);
	float_tt potVal = 0;

	for (unsigned iz = 0; iz < _mc->n[2]; iz++){
		for (unsigned ix = 0; ix < _mc->n[0]; ix++) {
			for (unsigned iy = 0; iy < _mc->n[1]; iy++) {
				sum[ix][iy] += _t[iz][ix][iy];
			}
		}
	}
	_persist->SaveProjectedPotential(sum);
}

void CPotential::SaveAtomicPotential(int znum) {
	std::stringstream str;
	str << "atomicPotential_";
	str << znum;
	_persist->SaveCx2DDataSet(_atomPot[znum], str.str());
}

/*------------------ splinh() -----------------------------*/
/*
 fit a quasi-Hermite cubic spline

 [1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
 'A New Method of Interpolation and Smooth
 Curve Fitting Based on Local Procedures'

 [2] H.Akima, Comm. ACM, 15(1972)p.914-918

 E. Kirkland 4-JUL-85
 changed zero test to be a small nonzero number 8-jul-85 ejk
 converted to C 24-jun-1995 ejk

 The inputs are:
 x[n] = array of x values in ascending order, each X(I) must
 be unique
 y[n] = array of y values corresponding to X(N)
 n = number of data points must be 2 or greater

 The outputs are (with z=x-x(i)):
 b[n] = array of spline coeficients for (x-x[i])
 c[n] = array of spline coeficients for (x-x[i])**2
 d[n] = array of spline coeficients for (x-x[i])**3
 ( x[i] <= x <= x[i+1] )
 To interpolate y(x) = yi + bi*z + c*z*z + d*z*z*z

 The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
 interval. NOTE that the last set of coefficients,
 b[n-1], c[n-1], d[n-1] are meaningless.
 */
void CPotential::splinh(float_tt x[], float_tt y[], std::vector<float_tt>& b,
		std::vector<float_tt>& c, std::vector<float_tt>& d, int n) {
#define SMALL 1.0e-25

	int i, nm1, nm4;
	float_tt m1, m2, m3, m4, m5, t1, t2, m54, m43, m32, m21, x43;

	if (n < 4)
		return;

	/* Do the first end point (special case), and get starting values */

	m5 = (y[3] - y[2]) / (x[3] - x[2]); /* mx = slope at pt x */
	m4 = (y[2] - y[1]) / (x[2] - x[1]);
	m3 = (y[1] - y[0]) / (x[1] - x[0]);

	m2 = m3 + m3 - m4; /* eq. (9) of reference [1] */
	m1 = m2 + m2 - m3;

	m54 = fabs(m5 - m4);
	m43 = fabs(m4 - m3);
	m32 = fabs(m3 - m2);
	m21 = fabs(m2 - m1);

	if ((m43 + m21) > SMALL)
		t1 = (m43 * m2 + m21 * m3) / (m43 + m21);
	else
		t1 = 0.5 * (m2 + m3);

	/* Do everything up to the last end points */

	nm1 = n - 1;
	nm4 = n - 4;

	for (i = 0; i < nm1; i++) {

		if ((m54 + m32) > SMALL)
			t2 = (m54 * m3 + m32 * m4) / (m54 + m32);
		else
			t2 = 0.5 * (m3 + m4);

		x43 = x[i + 1] - x[i];
		b[i] = t1;
		c[i] = (3.0 * m3 - t1 - t1 - t2) / x43;
		d[i] = (t1 + t2 - m3 - m3) / (x43 * x43);

		m1 = m2;
		m2 = m3;
		m3 = m4;
		m4 = m5;
		if (i < nm4) {
			m5 = (y[i + 4] - y[i + 3]) / (x[i + 4] - x[i + 3]);
		} else {
			m5 = m4 + m4 - m3;
		}

		m21 = m32;
		m32 = m43;
		m43 = m54;
		m54 = fabs(m5 - m4);
		t1 = t2;
	}

	return;

} /* end splinh() */

/*----------------------- seval() ----------------------*/
/*
 Interpolate from cubic spline coefficients

 E. Kirkland 4-JUL-85
 modified to do a binary search for efficiency 13-Oct-1994 ejk
 converted to C 26-jun-1995 ejk
 fixed problem on end-of-range 16-July-1995 ejk

 The inputs are:
 x[n] = array of x values in ascending order, each x[i] must
 be unique
 y[n] = array of y values corresponding to x[n]
 b[n] = array of spline coeficients for (x-x[i])
 c[n] = array of spline coeficients for (x-x[i])**2
 d[n] = array of spline coeficients for (x-x[i])**3
 n = number of data points
 x0 = the x value to interpolate at
 (x[i] <= x <= x[i+1]) and all inputs remain unchanged

 The value returned is the interpolated y value.

 The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
 interval. NOTE that the last set of coefficients,
 b[n-1], c[n-1], d[n-1] are meaningless.
 */
float_tt CPotential::seval(float_tt *x, float_tt *y, std::vector<float_tt>& b,
		std::vector<float_tt>& c, std::vector<float_tt>& d, int n,
		float_tt x0) {
	int i, j, k;
	float_tt z, seval1;

	/* exit if x0 is outside the spline range */
	if (x0 <= x[0])
		i = 0;
	else if (x0 >= x[n - 2])
		i = n - 2;
	else {
		i = 0;
		j = n;
		do {
			k = (i + j) / 2;
			if (x0 < x[k])
				j = k;
			else if (x0 >= x[k])
				i = k;
		} while ((j - i) > 1);
	}

	z = x0 - x[i];
	seval1 = y[i] + (b[i] + (c[i] + d[i] * z) * z) * z;

	return (seval1);

} /* end seval() */

af::array CPotential::GetSubPotential(int startx, int starty, int nx, int ny){
	return _t_d(af::seq(startx, startx + nx -1), af::seq(starty, starty + ny -1), af::span);
}

af::array& CPotential::GetPotential(){
	return _t_d;
}

void CPotential::GetSizePixels(unsigned int &nx, unsigned int &ny) const {
	nx = _mc->n[0];
	ny = _mc->n[1];
}

af:: array CPotential::GetSlice(af::array& t, unsigned idx){
	return t(af::span, af::span, idx);
}

void CPotential::SavePotential(){
	_persist->SavePotential(_t);
	if (_oc->SaveProjectedPotential || _oc->ComputeFromProjectedPotential){
//		if(!_persist->_potSaved){
//			_t.resize(boost::extents[_t_d.dims(2)][_t_d.dims(1)][_t_d.dims(0)]);
//			_t_d.host(_t.data());
//		}
		WriteProjectedPotential();
	}

}
void CPotential::ResizeSlices() {
}

} // end namespace slicepp
