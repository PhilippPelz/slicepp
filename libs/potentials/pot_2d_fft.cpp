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
#include <boost/format.hpp>
using boost::format;
namespace QSTEM {

C2DFFTPotential::C2DFFTPotential(const ConfigPtr& c,
		const PersistenceManagerPtr& p) :
		CPotential(c, p) {
	_atomPot = std::map<int, af::array>();
}

void C2DFFTPotential::DisplayParams() {
	CPotential::DisplayParams();
	BOOST_LOG_TRIVIAL(info)<<format("* Potential calculation: 2D (FFT method)");
}

void C2DFFTPotential::MakeSlices(superCellBoxPtr info) {
	CPotential::MakeSlices(info);
	/* check whether we have constant slice thickness */
	for (unsigned i = 0; i < _c->Model.nSlices; i++) {
		if (_sliceThicknesses[0] != _sliceThicknesses[i]) {
			printf(
					"Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",
					i);
			break;
		}
	}
}

void C2DFFTPotential::AddAtomToSlices(atom& atom, float_tt atomX,
		float_tt atomY, float_tt atomZ) {
	unsigned iAtomX = (int) floor(atomX / _c->Model.dx);
	unsigned iAtomY = (int) floor(atomY / _c->Model.dy);
	if (_c->Potential.periodicXY) {
		AddAtomPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	} else {
		AddAtomNonPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
	}
}
void C2DFFTPotential::CenterAtomZ(std::vector<atom>::iterator &atom,
		float_tt &z) {

}
void C2DFFTPotential::AddAtomNonPeriodic(atom& atom, float_tt atomBoxX,
		int iAtomX, float_tt atomBoxY, int iAtomY, float_tt atomZ) {
	int iAtomZ = (int) floor(atomZ / _c->Model.dz);
	int iax0, iay0, potentialOffsetX = 0, potentialOffsetY = 0;
	int iOffsetX = rint(_c->Structure.xOffset / _c->Model.dx);
	int iOffsetY = rint(_c->Structure.yOffset / _c->Model.dy);
	if (iAtomX - _iRadX + iOffsetX < 0) {
		iax0 = 0;
		potentialOffsetX = abs(iAtomX - _iRadX + iOffsetX) * OVERSAMPLING;
	} else {
		iax0 = iAtomX - _iRadX + iOffsetX;
	}
	if (iAtomY - _iRadY + iOffsetY < 0) {
		iay0 = 0;
		potentialOffsetY = abs(iAtomY - _iRadY + iOffsetY) * OVERSAMPLING;
	} else {
		iay0 = iAtomY - _iRadY + iOffsetY;
	}
	int iax1 =
			iAtomX + _iRadX + iOffsetX >= _c->Model.nx ?
					_c->Model.nx - 1 : iAtomX + _iRadX + iOffsetX;
	int iay1 =
			iAtomY + _iRadY + iOffsetY >= _c->Model.ny ?
					_c->Model.ny - 1 : iAtomY + _iRadY + iOffsetY;
	float_tt ddx = (atomBoxX / _c->Model.dx - iAtomX);
	float_tt ddy = (atomBoxY / _c->Model.dy - iAtomY);
	int iOffsX = (int) floor(ddx);
	int iOffsY = (int) floor(ddy);
	ddx -= (float_tt) iOffsX;
	ddy -= (float_tt) iOffsY;
	float_tt s11 = (1 - ddx) * (1 - ddy);
	float_tt s12 = (1 - ddx) * ddy;
	float_tt s21 = ddx * (1 - ddy);
	float_tt s22 = ddx * ddy;

	complex_tt added = complex_tt(0, 0);

	af::array iax = af::seq(iax0, iax1 - 1, 1);
	af::array iay = af::seq(iay0, iay1 - 1, 1);
	af::array xindex = (iOffsX + OVERSAMPLING * (iax - iax0));
	af::array yindex = (iOffsY + OVERSAMPLING * (iay - iay0));

	af::array vz;

	// shift the potentials, the cost is reasonable
	af::array pot12 = af::shift(_atomPot[atom.Znum], -1, 0);
	af::array pot21 = af::shift(_atomPot[atom.Znum], 0, -1);
	af::array pot22 = af::shift(_atomPot[atom.Znum], -1, -1);

	vz = s11 * _atomPot[atom.Znum] + s12 * pot12 + s21 * pot21 + s22 * pot22;
	added = af::sum<float_tt>(af::real(vz(yindex, xindex)));
	vz = af::real(vz);
	iax = (iax) % (_c->Model.nx);
	af::array xc = (iax < 0);
	iax += _c->Model.nx * xc;
	iay = (iay) % (_c->Model.ny);
	af::array yc = (iay < 0);
	iay += _c->Model.ny * yc;
	_t_af(iay, iax, iAtomZ) += vz(yindex, xindex);

//ComplexArray2D pot = _atomPot[atom.Znum];
	BOOST_LOG_TRIVIAL(trace)<< format("atom xyz (%-02.3f,%-02.3f,%-02.3f) Ixyz (%-3d,%-3d,%-3d) iax (%-3d .. %-3d) iay (%-3d .. %-3d)")
	% atom.r[0] % atom.r[1] % atom.r[2] % iAtomX % iAtomY % iAtomZ % iax0 %iax1%iay0%iay1;
}

void C2DFFTPotential::AddAtomPeriodic(atom& atom, float_tt atomBoxX, int iAtomX,
		float_tt atomBoxY, int iAtomY, float_tt atomZ) {
	unsigned iAtomZ = (int) floor(atomZ / _c->Model.dz);
	unsigned iax0 = iAtomX - _iRadX + _c->Model.nx;
	unsigned iax1 = iAtomX + _iRadX + _c->Model.nx;
	unsigned iay0 = iAtomY - _iRadY + _c->Model.ny;
	unsigned iay1 = iAtomY + _iRadY + _c->Model.ny;
//	_nx = 2 * OVERSAMPLING * (int) ceil(_config->Potential.AtomRadiusAngstrom / _config->Model.dx);
//	m_iRadX = (int) ceil( _config->Potential.AtomRadiusAngstrom / _config->Model.dx);

	float_tt ddx = (atomBoxX / _c->Model.dx - iAtomX);
	float_tt ddy = (atomBoxY / _c->Model.dy - iAtomY);
	unsigned iOffsX = (int) floor(ddx);
	unsigned iOffsY = (int) floor(ddy);
	ddx -= (float_tt) iOffsX;
	ddy -= (float_tt) iOffsY;

	float_tt s22 = (1 - ddx) * (1 - ddy);
	float_tt s21 = (1 - ddx) * ddy;
	float_tt s12 = ddx * (1 - ddy);
	float_tt s11 = ddx * ddy;

	complex_tt added = complex_tt(0, 0);

	af::array iax = af::seq(iax0, iax1 - 1, 1);
	af::array iay = af::seq(iay0, iay1 - 1, 1);
	af::array xindex = (iOffsX + OVERSAMPLING * (iax - iax0));
	af::array yindex = (iOffsY + OVERSAMPLING * (iay - iay0));

	af::array vz;

	// shift the potentials, the cost is reasonable
	af::array pot12 = af::shift(_atomPot[atom.Znum], -1, 0);
	af::array pot21 = af::shift(_atomPot[atom.Znum], 0, -1);
	af::array pot22 = af::shift(_atomPot[atom.Znum], -1, -1);

	vz = s11 * _atomPot[atom.Znum] + s12 * pot12 + s21 * pot21 + s22 * pot22;
	added = af::sum<float_tt>(af::real(vz(yindex, xindex)));
	vz = af::real(vz);
	iax = (iax) % (_c->Model.nx);
	af::array xc = (iax < 0);
	iax += _c->Model.nx * xc;
	iay = (iay) % (_c->Model.ny);
	af::array yc = (iay < 0);
	iay += _c->Model.ny * yc;
	_t_af(iay, iax, iAtomZ) += vz(yindex, xindex);
	// flattened array computations are faster and do not produce weird errors
}
void C2DFFTPotential::SliceSetup() {
	CPotential::SliceSetup();
	if (_atomPot.size() == 0) {
		_nx = 2 * OVERSAMPLING
				* (int) ceil(_c->Potential.AtomRadiusAngstrom / _c->Model.dx);
		_ny = 2 * OVERSAMPLING
				* (int) ceil(_c->Potential.AtomRadiusAngstrom / _c->Model.dy);
		_dkx = 0.5 * OVERSAMPLING / ((_nx) * _c->Model.dx);
		_dky = 0.5 * OVERSAMPLING / ((_ny) * _c->Model.dy);
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
					if (scatPar[0][N_SF - 4 - ix]
							< scatPar[0][N_SF - 3] - 0.001 * (ix + 1))
						break;
					scatPar[0][N_SF - 4 - ix] = scatPar[0][N_SF - 3]
							- 0.001 * (ix + 1);
					for (unsigned iy = 1; iy < N_ELEM; iy++)
						scatPar[iy][N_SF - 4 - ix] = 0;
				}
			}

			if (_c->Output.LogLevel < 2)
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
		_atomPot.insert(std::pair<int, af::array>(Znum, af::constant(0, _nx, _ny)));
		ComplexArray2D atomPot(boost::extents[_nx][_ny]);
		for (int ix = 0; ix < _nx; ix++) {
			float_tt kx = _dkx * (ix < _nx / 2 ? ix : _nx - ix);
			for (int iy = 0; iy < _ny; iy++) {
				float_tt ky = _dky * (iy < _ny / 2 ? iy : _ny - iy);
				float_tt s2 = (kx * kx + ky * ky);
				// if this is within the allowed circle:
				if (s2 < _kmax2) {
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					float_tt f = seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(s2))* exp(-s2 * B * 0.25);

					float_tt phase = PI * (kx * _c->Model.dx * _nx + ky * _c->Model.dy * _ny);
					atomPot[ix][iy] = complex_tt(f * cos(phase),f * sin(phase));
				}else{
					atomPot[ix][iy] = complex_tt(0, 0);
				}
			}
		}

#if FLOAT_PRECISION == 1
		fftwf_complex *ptr = reinterpret_cast<fftwf_complex*>(atomPot.data());
		fftwf_plan plan = fftwf_plan_dft_2d(_nx, _ny, ptr, ptr, FFTW_BACKWARD,
		FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
#else
		fftw_complex *ptr=(fftw_complex *)(atomPot.data());
		fftw_plan plan = fftw_plan_dft_2d(_nx,_ny,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
#endif
//		_backward.fft(_atomPot[Znum].data());

		for (unsigned ix = 0; ix < _nx; ix++)
			for (unsigned iy = 0; iy < _ny; iy++) {
				atomPot[ix][iy] *= _dkx * _dky
						* (OVERSAMPLING * OVERSAMPLING);
//				cout<<_atomPot[Znum][ix][iy]<<endl;
			}

		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
		_atomPot[Znum] = af::array(_nx, _ny, (afcfloat *)atomPot.data(), afHost);
		BOOST_LOG_TRIVIAL(info)<< format("Created 2D %d x %d potential array for Z=%d (B=%g A^2)") % _nx % _ny% Znum% B;
	}

}
void C2DFFTPotential::SaveAtomicPotential(int znum) {
	std::stringstream str;
	str << "atomicPotential_";
	str << znum;
	_persist->Save2DDataSet(_atomPot[znum], str.str());
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

} // end namespace QSTEM
