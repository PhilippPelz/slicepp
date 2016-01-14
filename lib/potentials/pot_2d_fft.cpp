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

#include "pot_2d_fft.hpp"
#include "fftw3.h"
const int k_fftMeasureFlag = FFTW_ESTIMATE;
#include <boost/format.hpp>

using boost::format;

namespace slicepp {

C2DFFTPotential::C2DFFTPotential(cModelConfPtr mc, cOutputConfPtr oc, cWaveConfPtr wc, PersistenceManagerPtr p) :
		CPotential(mc, oc,wc, p) {
	_atomPot = std::map<int, ComplexArray2D>();
}

void C2DFFTPotential::DisplayParams() {
	CPotential::DisplayParams();
	BOOST_LOG_TRIVIAL(info)<<format("* Potential calculation: 2D (FFT method)");
}

void C2DFFTPotential::MakeSlices(superCellBoxPtr info) {
	CPotential::MakeSlices(info);
	/* check whether we have constant slice thickness */
	for (unsigned i = 0; i < _mc->n[2]; i++) {
		if (_sliceThicknesses[0] != _sliceThicknesses[i]) {
			printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n", i);
			break;
		}
	}
}

void C2DFFTPotential::AddAtomToSlices(atom& atom) {
	unsigned iAtomX = (int) floor(atom.r[0] / _mc->d[0]);
	unsigned iAtomY = (int) floor(atom.r[1] / _mc->d[1]);

	if (_mc->periodicXY) {
		AddAtomPeriodic(atom, atom.r[0], iAtomX, atom.r[1], iAtomY, atom.r[2]);
	} else {
		AddAtomNonPeriodic(atom, atom.r[0], iAtomX, atom.r[1], iAtomY, atom.r[2]);
	}
}
void C2DFFTPotential::CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z) {

}
void C2DFFTPotential::AddAtomNonPeriodic(atom& atom, float_tt atomBoxX, int iAtomX, float_tt atomBoxY, int iAtomY, float_tt atomZ) {
	int iAtomZ = (int) floor(atomZ / _mc->d[2]);
	int iax0, iay0, potentialOffsetX = 0, potentialOffsetY = 0;
	int iOffsetX = rint(_mc->offset[0] / _mc->d[0]);
	int iOffsetY = rint(_mc->offset[1] / _mc->d[1]);
	if (iAtomX - _nRadX + iOffsetX < 0) {
		iax0 = 0;
		potentialOffsetX = abs(iAtomX - _nRadX + iOffsetX) * OVERSAMPLING;
	} else {
		iax0 = iAtomX - _nRadX + iOffsetX;
	}
	if (iAtomY - _nRadY + iOffsetY < 0) {
		iay0 = 0;
		potentialOffsetY = abs(iAtomY - _nRadY + iOffsetY) * OVERSAMPLING;
	} else {
		iay0 = iAtomY - _nRadY + iOffsetY;
	}
	int iax1 = iAtomX + _nRadX + iOffsetX >= _mc->n[0] ? _mc->n[0] - 1 : iAtomX + _nRadX + iOffsetX;
	int iay1 = iAtomY + _nRadY + iOffsetY >= _mc->n[1] ? _mc->n[1] - 1 : iAtomY + _nRadY + iOffsetY;
	float_tt ddx = (atomBoxX / _mc->d[0] - iAtomX);
	float_tt ddy = (atomBoxY / _mc->d[1] - iAtomY);
	int iOffsX = (int) floor(ddx);
	int iOffsY = (int) floor(ddy);
	ddx -= (float_tt) iOffsX;
	ddy -= (float_tt) iOffsY;
	float_tt s11 = (1 - ddx) * (1 - ddy);
	float_tt s12 = (1 - ddx) * ddy;
	float_tt s21 = ddx * (1 - ddy);
	float_tt s22 = ddx * ddy;
	ComplexArray2D pot = _atomPot[atom.Znum];
	complex_tt added = complex_tt(0, 0);
	BOOST_LOG_TRIVIAL(trace)<< format("atom xyz (%-02.3f,%-02.3f,%-02.3f) Ixyz (%-3d,%-3d,%-3d) iax (%-3d .. %-3d) iay (%-3d .. %-3d)")
	% atom.r[0] % atom.r[1] % atom.r[2] % iAtomX % iAtomY % iAtomZ % iax0 %iax1%iay0%iay1;

	for (int iax = iax0; iax < iax1; iax++) {
		for (int iay = iay0; iay < iay1; iay++) {
			int xindex = iOffsX + OVERSAMPLING * (iax - iax0) + potentialOffsetX;
			int yindex = iOffsY + OVERSAMPLING * (iay - iay0) + potentialOffsetY;
			float_tt vz = (s11 * pot[xindex][yindex] + s12 * pot[xindex + 1][yindex] + s21 * pot[xindex][yindex + 1]
					+ s22 * pot[xindex + 1][yindex + 1]).real();
			_t[iAtomZ][iax][iay] += complex_tt(vz, 0);

			added += vz;
		}
	}
}

void C2DFFTPotential::AddAtomPeriodic(atom& atom, float_tt atomBoxX, int iAtomX, float_tt atomBoxY, int iAtomY, float_tt atomZ) {
	unsigned iAtomZ = (int) floor(atomZ / _mc->d[2]);
	unsigned iax0 = iAtomX - _nRadX + _mc->n[0];
	unsigned iax1 = iAtomX + _nRadX + _mc->n[0];
	unsigned iay0 = iAtomY - _nRadY + _mc->n[1];
	unsigned iay1 = iAtomY + _nRadY + _mc->n[1];
//	_nx = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx);
//	m_iRadX = (int) ceil( _config->Potential.AtomRadiusAngstrom / _config->Model.dx);
	float_tt ddx = (atomBoxX / _mc->d[0] - iAtomX);
	float_tt ddy = (atomBoxY / _mc->d[1] - iAtomY);
	unsigned iOffsX = (int) floor(ddx);
	unsigned iOffsY = (int) floor(ddy);
	ddx -= (float_tt) iOffsX;
	ddy -= (float_tt) iOffsY;

	float_tt s22 = (1 - ddx) * (1 - ddy);
	float_tt s21 = (1 - ddx) * ddy;
	float_tt s12 = ddx * (1 - ddy);
	float_tt s11 = ddx * ddy;

	ComplexArray2DPtr pot = _atomPot[atom.Znum];
	complex_tt added = complex_tt(0, 0);
	for (int iax = iax0; iax < iax1; iax++) { // TODO: should use ix += OVERSAMPLING
		for (int iay = iay0; iay < iay1; iay++) {
			int xindex = (iOffsX + OVERSAMPLING * (iax - iax0));
			int yindex = iOffsY + OVERSAMPLING * (iay - iay0);
			float_tt vz = (s11 * pot[xindex][yindex] + s12 * pot[xindex + 1][yindex] + s21 * pot[xindex][yindex + 1]
					+ s22 * pot[xindex + 1][yindex + 1]).real();
			_t[iAtomZ][iax % _mc->n[0]][iay % _mc->n[1]] += complex_tt(vz, 0);
			added += vz;
		}
	}
}
void C2DFFTPotential::SliceSetup() {
	CPotential::SliceSetup();
	if (_atomPot.size() == 0) {
		_nx = 2 * OVERSAMPLING * (int) ceil(_mc->ratom / _mc->d[0]);
		_ny = 2 * OVERSAMPLING * (int) ceil(_mc->ratom / _mc->d[1]);
		_dkx = 0.5 * OVERSAMPLING / ((_nx) * _mc->d[0]);
		_dky = 0.5 * OVERSAMPLING / ((_ny) * _mc->d[1]);
		_kmax2 = 0.5 * _nx * _dkx / (float_tt) OVERSAMPLING; // largest k that we'll admit

		BOOST_LOG_TRIVIAL(info)<< format("Cutoff scattering angle: kmax = %g (1/Å)") % _kmax2;
		scatPar[0][N_SF - 1] = 1.2 * _kmax2;
		scatPar[0][N_SF - 2] = 1.1 * _kmax2;
		scatPar[0][N_SF - 3] = _kmax2;
		if (scatPar[0][N_SF - 4] > scatPar[0][N_SF - 3]) {
			unsigned ix = 0;
			if (1) {
				// set additional scattering parameters to zero:
				for (ix; ix < N_SF - 10; ix++) {
					if (scatPar[0][N_SF - 4 - ix] < scatPar[0][N_SF - 3] - 0.001 * (ix + 1))
						break;
					scatPar[0][N_SF - 4 - ix] = scatPar[0][N_SF - 3] - 0.001 * (ix + 1);
					for (unsigned iy = 1; iy < N_ELEM; iy++)
						scatPar[iy][N_SF - 4 - ix] = 0;
				}
			}

			if (_oc->LogLevel < 2)
				BOOST_LOG_TRIVIAL(info)<< format("Reduced angular range of scattering factor to %g/Å!") % scatParOffs[0][N_SF - 4 - ix];
			}
		_kmax2 *= _kmax2;
	}
}
void C2DFFTPotential::ComputeAtomPotential(int Znum) {

	std::vector<float_tt> splinb(N_SF, 0), splinc(N_SF, 0), splind(N_SF, 0);
	float_tt B = 0; // TODO what does this do   _config->Model.UseTDS ? 0 : atom->dw;
	if (_atomPot.count(Znum) == 0) {
		// setup cubic spline interpolation:
		splinh(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF);
		_atomPot.insert(std::pair<int, ComplexArray2D>(Znum, ComplexArray2D(boost::extents[_ndiaAtomX][_ndiaAtomY])));
		for (unsigned ix = 0; ix < _ndiaAtomX; ix++) {
			float_tt kx = _dkx * (ix < _ndiaAtomX / 2 ? ix : _ndiaAtomX - ix);
			for (unsigned iy = 0; iy < _ndiaAtomY; iy++) {
				float_tt ky = _dky * (iy < _ndiaAtomY / 2 ? iy : _ndiaAtomY - iy);
				float_tt s2 = (kx * kx + ky * ky);
				// if this is within the allowed circle:
				if (s2 < _kmax2) {
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					float_tt f = seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(s2)) * exp(-s2 * B * 0.25);
					float_tt phase = PI * (kx * _mc->d[0] * _nx + ky * _mc->d[1] * _ny);
					_atomPot[Znum][ix][iy] = complex_tt(f * cos(phase), f * sin(phase));
				}
			}
		}

#if FLOAT_PRECISION == 1
		fftwf_complex *ptr = reinterpret_cast<fftwf_complex*>(_atomPot[Znum].data());
		fftwf_plan plan = fftwf_plan_dft_2d(_nx, _ny, ptr, ptr, FFTW_BACKWARD,
		FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
#else
		fftw_complex *ptr=(fftw_complex *)(_atomPot[Znum].data());
		fftw_plan plan = fftw_plan_dft_2d(_nx,_ny,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
#endif

		for (unsigned ix = 0; ix < _ndiaAtomX; ix++)
			for (unsigned iy = 0; iy < _ndiaAtomY; iy++) {
				_atomPot[Znum][ix][iy] *= _dkx * _dky * (OVERSAMPLING * OVERSAMPLING);
			}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
		BOOST_LOG_TRIVIAL(info)<< format("Created 2D %d x %d potential array for Z=%d (B=%g A^2)") % _nx % _ny% Znum% B;
	}

}


#define PHI_SCALE 47.87658
// #define SHOW_SINGLE_POTENTIAL 0
////////////////////////////////////////////////////////////////////////////
// This function should be used yet, because it computes the projected
// potential wrongly, since it doe not yet perform the projection!!!
complex_tt *C2DFFTPotential::GetAtomPotential2D(int Znum, float_tt B) {

}
#undef PHI_SCALE
#undef SHOW_SINGLE_POTENTIAL

} // end namespace slicepp
