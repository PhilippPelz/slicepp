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

#include "pot_3d_fft.hpp"

#include <cstring>
#include <string>
#include <boost/format.hpp>

using boost::format;

namespace QSTEM {

C3DFFTPotential::C3DFFTPotential(const ConfigPtr& c, const PersistenceManagerPtr& p) :
		_ndiaAtomZ(0),
		_nzPerSlice(0),
		_nrAtomZ(0),
		_currentZnum(0),
		CPotential(c, p) {
}

void C3DFFTPotential::DisplayParams() {
	CPotential::DisplayParams();
	BOOST_LOG_TRIVIAL(info)<< format("* Potential calculation: 3D (FFT method)");
	BOOST_LOG_TRIVIAL(info)<< format("* Will use %d sampling points per slice, total nz=%d (%d)") % _nzPerSlice % _ndiaAtomZ % (_nzPerSlice >> 1);
	BOOST_LOG_TRIVIAL(info)<< format("* dkx = %g, nx = %d, kmax2 = %g") % _dkx % _ndiaAtomX % _kmax;
	BOOST_LOG_TRIVIAL(info)<< format("* Cutoff scattering angle: kmax=%g (1/A), dk=(%g,%g %g)") % _kmax % _dkx % _dkx % _dkz;
	BOOST_LOG_TRIVIAL(info)<< format("* At. Pot. arrays:      %d x %d (r x z)") % (2*_nrAtomTrans) % _ndiaAtomZ;
	BOOST_LOG_TRIVIAL(info)<< format("* Creating up to %d intermediate arrays of size %d x %d to store sliced potentials")
			% _nzPerSlice % (2*_nrAtomTrans) % _ndiaAtomZ;
	BOOST_LOG_TRIVIAL(info)<< format("* Potential offset: will use %d sampling points per slice, total nz=%d (%d)")
			% _nzPerSlice % _c->Model->nSlices % (_nzPerSlice >> 1);
}
/********************************************************************************
 * Clean up after the total potential has been created
 ********************************************************************************/
void C3DFFTPotential::CleanUp(){
	_intAtPot.clear();
	_atPot.clear();
	_atPotOffset.clear();
}
//				// 2*pi*e*a0 = 47.88;
//				// with e = 14.4 Volt-Angstrom
//				// a0 = 0.5292 Angstrom
/********************************************************************************
 * Computes and saves a properly integrated, sliced version of V(r,z) corresponding to the ratio
 * that the atom center's position has in its current slice
 *
 * 2*pi*e*a0 = 47.88;
 * with e = 14.4 Volt-Angstrom
 * a0 = 0.5292 Angstrom
 *
 * kirkland 2nd ed. p. 253
 ********************************************************************************/
FloatArray2DView C3DFFTPotential::getSlicedAtomicPotentialAt(ratio rat, int Znum) {
	if(_currentZnum == 0 || _currentZnum != Znum){
		_intAtPot.clear();
		_currentZnum = Znum;
	}
	ComplexArray2DView atomPot = _atPot[Znum][boost::indices[range(0, 2 * _nrAtomTrans)][range(0, _ndiaAtomZ)]];
	if (_intAtPot.count(rat) == 0) {
		float_tt s =  47.88 * _dz;
		_intAtPot[rat] = FloatArray2D();
		_intAtPot[rat].resize(boost::extents[2 * _nrAtomTrans][2 * _nRadZ+1]);
		FloatArray2D& ia = _intAtPot[rat];
		std::fill(_intAtPot[rat].data(), _intAtPot[rat].data() + _intAtPot[rat].size(), 0);
		for (int iax = -_nRadX; iax <= _nRadX; iax++) {
			float_tt x2 = (iax  ) * _c->Model->dx;
			x2 *= x2;
			for (int iay = -_nRadY; iay <= _nRadY; iay++) {
				float_tt y2 = (iay  ) * _c->Model->dy;
				y2 *= y2;
				float_tt r = sqrt(x2 + y2); // TODO: LUT for sqrt or precompute r
				int ir = (int) floor(r / _dr);
				float_tt ddr = r / _dr - ir;
				if (ir < _nrAtomTrans - 1) {
					// first slice
					int i = 0;
					unsigned end = (1 - rat).numerator();
					float_tt tmp = 0;
					for (unsigned iaz = 0; iaz < end; iaz++) {
						tmp += atomPot[ir][iaz].real();
					}
					ia[ir][i++] = tmp * s * end;
					// middle slices
					while (end + _nzPerSlice < _ndiaAtomZ) {
						tmp = 0;
						for (unsigned iaz = end; iaz < end + _nzPerSlice;
								iaz++) {
							tmp += atomPot[ir][iaz].real();
						}
						ia[ir][i] = tmp * s * _nzPerSlice;
						i++;
						end += _nzPerSlice;
					}
					tmp = 0;
					// last slice
					for (int iaz = end; iaz < _ndiaAtomZ; iaz++) {
						tmp += atomPot[ir][iaz].real();
					}
					ia[ir][i] = tmp * s * (_ndiaAtomZ - end);
				}
			}
		}
	}
	return _intAtPot[rat][boost::indices[range(0, _nrAtomTrans)][range(0, 2 * _nRadZ+1)]];
}
void C3DFFTPotential::AddAtomNonPeriodic(atom& atom, unsigned int iAtomX, unsigned int iAtomY) {
	ComplexVector atPotOffsPtr;
	float_tt fatomZ = atom.r[2] / _c->Model->dz;
	unsigned iAtomZ = (int) ceil(fatomZ);
	float_tt dAtomZ = abs(fatomZ - iAtomZ + 1);
	unsigned ndAtomZ = (unsigned) floor(dAtomZ / _dz);
	FloatArray2DView atPot = getSlicedAtomicPotentialAt(ratio(ndAtomZ,_nzPerSlice),atom.Znum);

	//we took care that we the atoms stay within the slices in z while building the crystal/superstructure
	int ixstart = (int) iAtomX - (int) _nRadX < 0 ? 0 : iAtomX - _nRadX;
	int iystart = (int) iAtomY - (int) _nRadY < 0 ? 0 : iAtomY - _nRadY;
	int izstart = iAtomZ - _nRadZ;

	int ixend = (int) iAtomX + (int) _nRadX >= _c->Model->nx ? _c->Model->nx - 1 : iAtomX + _nRadX;
	int iyend = (int) iAtomY + (int) _nRadY >= _c->Model->ny ? _c->Model->ny - 1 : iAtomY + _nRadY;
	int izend = iAtomZ + _nRadZ;

	if (_c->Potential->UseQPotentialOffsets) {
//			GetAtomPotentialOffset3D(atom->Znum,muls,m_tds ? 0 : atoms[iatom].dw,nzPerSlice,nRadius,Nz_lut,atom->q, atPotOffsPtr);
	}

	for (int iax = ixstart; iax <= ixend; iax++) {
		float_tt x2 = (iax + 0.5) * _c->Model->dx - atom.r[0];
		x2 *= x2;
		for (int iay = iystart; iay <= iyend; iay++) {
			float_tt y2 = (iay + 0.5) * _c->Model->dy - atom.r[1];
			y2 *= y2;
			float_tt r = sqrt(x2 + y2); // TODO: LUT for sqrt or precompute r
			int ir = (int) floor(r / _dr);
			float_tt ddr = r / _dr - ir;
			if (ir < _nrAtomTrans - 2) {
				for (int iaz = izstart; iaz <= izend; iaz++) {
					float_tt potVal = (1-ddr) * atPot[ir+1][iaz-izstart]
									+ ddr * atPot[ir][iaz-izstart] ;
					if (atPotOffsPtr.size() != 0) { // TODO
					}
					_t[iaz][iay][iax] += potVal;
				}
			}
		}
	}
}

void C3DFFTPotential::AddAtomToSlices(atom& atom) {
	unsigned iAtomX = (int) floor(atom.r[0] / _c->Model->dx);
	unsigned iAtomY = (int) floor(atom.r[1] / _c->Model->dy);

	if (_c->Potential->periodicXY) {
		AddAtomPeriodic(atom,  iAtomX,  iAtomY);
	} else {
		AddAtomNonPeriodic(atom,  iAtomX,  iAtomY);
	}
}

void C3DFFTPotential::AddAtomPeriodic(atom& atom, unsigned iAtomX, unsigned iAtomY) {
	int iAtomZ = (int) rint(atom.r[2] / _c->Model->dz) + _nRadZ;
	int iax0 = iAtomX - _nRadX;
	int iax1 = iAtomX + _nRadX;
	int iay0 = iAtomY - _nRadY;
	int iay1 = iAtomY + _nRadY;

	//	int iaz0 = iAtomZ + atomRadiusZOffset - m_iRadZ < 0 ? -(iAtomZ+atomRadiusZOffset) : -m_iRadZ;
	//	int iaz1 =	iAtomZ + m_iRadZ >= _config->Model.nSlices ? _config->Model.nSlices - (iAtomZ+atomRadiusZOffset) - 1 : m_iRadZ;

	// if (iatom < 2) printf("iatomZ: %d, %d cz=%g, %g: %d, %d\n",iAtomZ,iaz0,_config->Model.dz,atomZ,(int)(-2.5-atomZ),(int)(atomZ+2.5));
	// printf("%d: iatomZ: %d, %d cz=%g, %g\n",iatom,iAtomZ,iaz0,_config->Model.dz,atomZ);
	//	bool isAtomZWithinSlices = (iAtomZ + iaz0 < _config->Model.nSlices) && (iAtomZ + iaz1 >= 0);
	//	if (isAtomZWithinSlices) {
	ComplexArray2D pot = _atPot[atom.Znum];

#if USE_Q_POT_OFFSETS
	// retrieve the pointer to the array of charge-dependent potential offset
	// This function will return NULL; if the charge of this atom is zero:
	ComplexVector atPotOffsPtr;
	BOOST_LOG_TRIVIAL(fatal) << "getAtomPotentialOffset3D not implemented yet";
	getAtomPotentialOffset3D(atoms[iatom].Znum,muls,m_tds ? 0 : atoms[iatom].dw,&_nzPerSlice,&_ndiaAtomX/2,&_ndiaAtomZ/2,atoms[iatom].q, atPotOffsPtr);
#endif // USE_Q_POT_OFFSETS

	int iOffsLimHi = _ndiaAtomX / 2 * (_ndiaAtomZ / 2 - 1);
	int iOffsLimLo = -_ndiaAtomX / 2 * (_ndiaAtomZ / 2 - 1);
	int iOffsStep = _nzPerSlice * _ndiaAtomX / 2;

	// Slices around the slice that this atom is located in must be affected by this atom:
	// iaz must be relative to the first slice of the atom potential box.
	for (int iax = iax0; iax < iax1; iax++) {
		for (int iay = iay0; iay < iay1; iay++) {
			float_tt x2 = iax * _c->Model->dx - atom.r[0];
			x2 *= x2;
			float_tt y2 = iay * _c->Model->dy - atom.r[1];
			y2 *= y2;
			float_tt ddr = sqrt(x2 + y2) / _dr;
			int iCurrentRadius = (int) floor(ddr);
			// add in different slices, once r has been defined
			if (iCurrentRadius < _ndiaAtomX / 2 - 1) {
				float fractionAboveCurrentRadius = ddr - (double) iCurrentRadius;
				// Include interpolation in z-direction as well (may do it in a very smart way here !):

				float_tt dOffsZ = (iAtomZ - _nRadZ - atom.r[2] / _c->Model->dz) * _nzPerSlice;
#if Z_INTERPOLATION
				int iOffsZ = (int)dOffsZ;
				float_tt ddz = fabs(dOffsZ - (double)iOffsZ);
#else // Z_INTERPOLATION
				int iOffsZ = (int) (dOffsZ + 0.5);
#endif // Z_INTERPOLATION
				iOffsZ *= _ndiaAtomX / 2;

				for (int iaz = -_nRadZ; iaz <= _nRadZ; iaz++) {
					float_tt potVal = 0;
					// iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/_config->Model.dz)*nzSub+0.5);
					if (iOffsZ < 0) {
						if (iOffsZ > iOffsLimLo) {
							// TODO fix this to use ComplexArray2D, fix the confusing indices
							//do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
							potVal = (1-ddz)*((1-fractionAboveCurrentRadius)*pot[iCurrentRadius-iOffsZ+_ndiaAtomX/2][0]
									+fractionAboveCurrentRadius*pot[iCurrentRadius+1-iOffsZ+_ndiaAtomX/2][0])+
							ddz *((1-fractionAboveCurrentRadius)*pot[iCurrentRadius-iOffsZ][0]+
									fractionAboveCurrentRadius*pot[iCurrentRadius+1-iOffsZ][0]);
#if USE_Q_POT_OFFSETS
							// add the charge-dependent potential offset
							if (atPotOffsPtr != NULL) {
								potVal += atoms[iatom].q*((1-ddz)*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius-iOffsZ+_ndiaAtomX/2][0]+
												fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1-iOffsZ+_ndiaAtomX/2][0])+
										ddz *((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius-iOffsZ ][0]+
												fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1-iOffsZ ][0]));
							}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
							//FIXME
							potVal =
									(1 - fractionAboveCurrentRadius)
											* pot[iCurrentRadius - iOffsZ
													+ _ndiaAtomX / 2][1000].real()
											+ fractionAboveCurrentRadius
													* pot[iCurrentRadius + 1
															- iOffsZ + _ndiaAtomX / 2][1000].real();
#if USE_Q_POT_OFFSETS
							// add the charge-dependent potential offset
							if (atPotOffsPtr != NULL) {
								potVal += atoms[iatom].q*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius-iOffsZ+_ndiaAtomX/2][0]+
										fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1-iOffsZ+_ndiaAtomX/2][0]);
							}
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
						}
					} else {
						// select the pointer to the right layer in the lookup box
						// printf("%4d: iOffsZ: %d, iaz: %d (slice: %d, pos: %g [%d .. %d])\n",iatom,iOffsZ,iaz,iAtomZ+iaz,atomZ,iaz0,iaz1);
						if (iOffsZ < iOffsLimHi) {
							// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
							potVal = (1-ddz)*((1-fractionAboveCurrentRadius)*pot[iCurrentRadius+iOffsZ][0]+
									fractionAboveCurrentRadius*pot[iCurrentRadius+1+iOffsZ][0])+
							ddz *((1-fractionAboveCurrentRadius)*pot[iCurrentRadius+iOffsZ+_ndiaAtomX][0]+
									fractionAboveCurrentRadius*pot[iCurrentRadius+1+iOffsZ+_ndiaAtomX][0]);
#if USE_Q_POT_OFFSETS
							// add the charge-dependent potential offset
							if (atPotOffsPtr != NULL) {
								potVal += atoms[iatom].q*((1-ddz)*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius+iOffsZ][0]+
												fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1+iOffsZ][0])+
										ddz *((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius+iOffsZ+_ndiaAtomX/2][0]+
												fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1+iOffsZ+_ndiaAtomX/2][0]));
							}
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
							//FIXME
							potVal =
									(1 - fractionAboveCurrentRadius)
											* pot[iCurrentRadius + iOffsZ][1000].real()
											+ fractionAboveCurrentRadius
													* pot[iCurrentRadius + 1
															+ iOffsZ][1000].real();
#if USE_Q_POT_OFFSETS
							// add the charge-dependent potential offset
							if (atPotOffsPtr != NULL) {
								potVal += atoms[iatom].q*((1-fractionAboveCurrentRadius)*atPotOffsPtr[iCurrentRadius+iOffsZ][0]+
										fractionAboveCurrentRadius*atPotOffsPtr[iCurrentRadius+1+iOffsZ][0]);
							}
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
						}
					}
					int ix = (iAtomX + iax + _c->Model->nx) % _c->Model->nx;
					int iy = (iAtomY + iay + _c->Model->ny) % _c->Model->ny;
					_t[iAtomZ + iaz][iy][ix] += potVal;

					BOOST_LOG_TRIVIAL(trace)<< boost::format("_trans1[%d ][%d][%d] += %f\n")
					% (iAtomZ+iaz)%((iAtomY+iay+_c->Model->ny) % _c->Model->ny)% ((iAtomX+iax+_c->Model->nx) % _c->Model->nx)%potVal;
					iOffsZ += iOffsStep;
				} // for iaz
			} // if ir < Nr-1

			// advance pointer to next complex potential point:
			//				potPtr += 2;
			// make imaginary part zero for now
			// *potPtr = 0;
			// potPtr++;

			// wrap around when end of y-line is reached:
			//				if (++iay % m_ny == 0)
			//					potPtr -= 2 * m_ny;
			//					iay = 0;
		}
	}
	//	} // iaz0+iAtomZ < m_slices
}

/************************************************************************************//**
 * Sets up all member variables used for the potential calculation
 *
 * The FFT-resolution in the z-direction must be high enough to avoid
 * artifacts due to premature cutoff of the rec. space scattering factor
 * we will therefore make it roughly the same as the x-resolution
 * However, we will make sure that a single slice contains an even number
 * of sampling points.
 *
 ***************************************************************************************/
void C3DFFTPotential::SliceSetup() {
	CPotential::SliceSetup();

	for (unsigned i = 0; i < _c->Model->nSlices; i++) {
		if (_sliceThicknesses[0] != _sliceThicknesses[i]) {
			BOOST_LOG_TRIVIAL(warning)<< format("Warning: slice thickness not constant, will give wrong results (iz=%d)!") % i;
		}
	}

	if (_atPot.size() == 0) {
		_nzPerSlice = (int) floor(OVERSAMPLING * _c->Model->dz / min(_c->Model->dx,_c->Model->dy));
		if (_nzPerSlice % 2 != 0) _nzPerSlice++;
		_dz = _c->Model->dz / _nzPerSlice;
		_nrAtomZ =  (int) ceil(_c->Potential->ratom /  _c->Model->dz) * _nzPerSlice;
		_ndiaAtomZ = (2 * _nrAtomZ);
		_dkz = _nzPerSlice / (_ndiaAtomZ * _c->Model->dz);
		SetScatteringFactors(_kmax);
	}
}

/********************************************************************************
 * Create Lookup table for 3D potential due to neutral atoms
 *
 * What this look-up procedure will do is to supply V(r,z) computed from fe(q).
 * Since V(r,z) is rotationally symmetric we might as well compute
 * V(x,y,z) at y=0, i.e. V(x,z).
 * In order to do the proper 3D inverse FT without having to do a complete 3D FFT
 * we will pre-compute the qy-integral for y=0.
 *
 * We also make sure that the potential touches zero at least somewhere. This will avoid
 * sharp edges that could produce ringing artifacts.
 * It is certainly debatable whether this is a good approach, or not.
 ********************************************************************************/
#define PHI_SCALE 47.87658
void C3DFFTPotential::ComputeAtomPotential(int Znum) {
	float_tt B = 0; // TODO how to compute this   _config->Model.UseTDS ? 0 : atom->dw; multiply scattering factor with Debye-Waller factor:
	std::vector<float_tt> splinb(N_SF), splinc(N_SF), splind(N_SF);
	if (_atPot.count(Znum) == 0) {
		// setup cubic spline interpolation:
		splinh((scatPar[0]), (scatPar[Znum]), splinb, splinc, splind, N_SF);
		_atPot[Znum] = ComplexArray2D();
		_atPot[Znum].resize(boost::extents[2*_nrAtomTrans][_ndiaAtomZ]);
		std::fill(_atPot[Znum].data(), _atPot[Znum].data() + _atPot[Znum].size(), complex_tt(0, 0));

		// define x-and z-position of atom center: The atom should sit in the top-left corner,
		// however (nzPerSlice+1)/2 above zero in z-direction
		int izOffset = (_nzPerSlice - 1);
		float_tt xPos = -2.0 * PI * 0.0; // or m_dx*nx/(OVERSAMPLING), if in center
		float_tt zPos = -2.0 * PI * (_c->Model->dz / _nzPerSlice * (izOffset/2));

		for (int iz = 0; iz < _ndiaAtomZ; iz++) {
			float_tt kz = _dkz * (iz < _nrAtomZ ? iz : iz - _ndiaAtomZ);
			for (int ix = 0; ix < 2*_nrAtomTrans; ix++) {
				float_tt kx = _dkx * (ix < _nrAtomTrans ? ix : ix - 2*_nrAtomTrans);
				float_tt kxz2 = (kx * kx + kz * kz);
				if (kxz2 < _kmax2) {
					float_tt f = seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(kxz2)) * exp(-kxz2 * B * 0.25);
					// perform the qy-integration for qy <> 0:
					for (int iy = 0; iy < _ndiaAtomY; iy++) {
						float_tt kxyz2 = _dky * iy;
						kxyz2 = kxyz2 * kxyz2 + kxz2;
						if (kxyz2 < _kmax2) {
							float_tt ss = seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(kxyz2)) * exp(-kxyz2 * B * 0.25);
							f += 2 * ss;
						} else
							break;
					}
					f *= _dky;
					// note that the factor 2 is missing in the phase (2pi k*r) this places the atoms in the center of the box.
					float_tt phase = kx * xPos + kz * zPos;
					_atPot[Znum][ix][iz] = complex_tt(f * cos(phase), f * sin(phase)); // *zScale
				}
			}
		}

		_persist->Save2DDataSet(_atPot[Znum], "atompotK_" + std::to_string(Znum));

		// This converts the 2D kx-kz map of the scattering factor to a 2D real space map.
		time_t time0, time1;
#if FLOAT_PRECISION ==1
		fftwf_complex *ptr=reinterpret_cast<fftwf_complex*>(_atPot[Znum].data());
		fftwf_plan plan = fftwf_plan_dft_2d(2*_nrAtomTrans,_ndiaAtomZ,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
#else
		fftw_complex *ptr = reinterpret_cast<fftw_complex*>(_atPot[Znum].data());
		fftw_plan plan = fftw_plan_dft_2d(2*_nrAtomTrans, _ndiaAtomZ, ptr, ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
#endif

		_persist->Save2DDataSet(_atPot[Znum], "atompotR_" + std::to_string(Znum));


//		for (unsigned ix = 0; ix < _nrAtomTrans; ix++) {
//			for (unsigned iz = 0; iz < _nrAtomZ; iz++) {
//				float_tt zScale = 0;
//				// Integrate over nzPerSlice neighboring layers here:::::::::
//				for (int iiz = 0; iiz <= izOffset; iiz++) {
//					unsigned idx = iz + iiz;
//					if (idx < _nrAtomZ) {
//						zScale += tmp[ix][idx].real();
//					}
//				}
//				if (zScale < 0)
//					zScale = 0;
//
//				// assign the iz-th slice the sum of the 3 other slices:
//				// and divide by unit cell volume (since this is in 3D):
//				// Where does the '1/2' come from??? OVERSAMPLING*OVERSAMP_Y/8 = 1/2
//				// remember, that s=0.5*k;
//
//				_atPot[Znum][ix][iz] = 47.88 * _dkx * _dkz * zScale;
//				// normalization of fftw
//				_atPot[Znum][ix][iz] /= (_ndiaAtomZ*2*_nrAtomTrans);
//				// 2*pi*e*a0 = 47.88;
//				// with e = 14.4 Volt-Angstrom
//				// a0 = 0.5292 Angstrom
//
//			}
//		}
		BOOST_LOG_TRIVIAL(info)<< format("Created 3D (r-z) %d x %d potential array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)")
		% (_ndiaAtomX / 2) % (_ndiaAtomZ / 2) % Znum % B % _dkx % _dkx % _dkz % izOffset;


	}
}
void C3DFFTPotential::SaveAtomicPotential(int znum) {
	std::stringstream str;
	str << "atomicPotential_";
	str << znum;
	_persist->Save2DDataSet(_atPot[znum], str.str());
}
/********************************************************************************
 * TODO: not implemented yet
 *
 * Lookup function for 3D potential offset due to charged atoms (ions)
 ********************************************************************************/
void C3DFFTPotential::GetAtomPotentialOffset3D(unsigned Znum, float_tt B,
		unsigned &nzSub, unsigned &Nr, unsigned &Nz_lut, float_tt q,
		ComplexVector &output) {

	// if there is no charge to this atom, return NULL:
//	if (q == 0)
//		return;
//
//	unsigned nx, ny, nz, nzPerSlice;
//	float_tt dkx, dky, dkz;
//	std::vector<float_tt> splinb(N_SF), splinc(N_SF), splind(N_SF);
//	float_tt kmax2;
//
//	ComplexVector temp;
//
//
//
//	// initialize this atom, if it has not been done yet:
//	if (_atPotOffset.count(Znum) == 0) {
//		// setup cubic spline interpolation:
//		splinh(scatParOffs[0], scatParOffs[Znum], splinb, splinc, splind, N_SF);
//
//		_atPotOffset[Znum] = ComplexVector(nx * nz / 4);
//		memset((void*) &temp[0], 0, nx * nz * sizeof(complex_tt));
//
//		//complex_tt *temp=(complex_tt *)fftw_malloc(nx*nz*sizeof(complex_tt));
//		//memset(temp, 0, nx*nz*sizeof(complex_tt));
//
//		float_tt kzmax = dkz * nz / 2.0;
//		// define x-and z-position of atom center:
//		// The atom should sit in the top-left corner,
//		// however (nzPerSlice+1)/2 above zero in z-direction
//		float_tt xPos = -2.0 * PI * 0.0; // or m_dx*nx/(OVERSAMPLING), if in center
//		unsigned izOffset = (nzPerSlice - 1) / 2;
//		float_tt zPos = -2.0 * PI * (_c->Model->dz / nzPerSlice * (izOffset));
//
//		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
//		for (unsigned iz = 0; iz < nz; iz++) {
//			float_tt kz = dkz * (iz < nz / 2 ? iz : iz - nz);
//			// We also need to taper off the potential in z-direction
//			// in order to avoid cutoff artifacts.
//			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
//			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
//			for (unsigned ix = 0; ix < nx; ix++) {
//				float_tt kx = dkx * (ix < nx / 2 ? ix : ix - nx);
//				float_tt s2 = (kx * kx + kz * kz);
//				// if this is within the allowed circle:
//				if (s2 < kmax2) {
//					unsigned ind3d = ix + iz * nx;
//					// f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
//					// multiply scattering factor with Debye-Waller factor:
//					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
//					float_tt f = seval(scatParOffs[0], scatParOffs[Znum],
//							splinb, splinc, splind, N_SF, sqrt(s2))
//							* exp(-s2 * B * 0.25);
//					// perform the qy-integration for qy <> 0:
//					for (unsigned iy = 1; iy < nx; iy++) {
//						float_tt s3 = dky * iy;
//						s3 = s3 * s3 + s2;
//						if (s3 < kmax2) {
//							f += 2 * seval(scatPar[0], scatPar[Znum], splinb, splinc, splind, N_SF, sqrt(s3))
//									* exp(-s3 * B * 0.25);
//						} else
//							break;
//					}
//					f *= dkx;
//					// note that the factor 2 is missing in the phase (2pi k*r)
//					// this places the atoms in the center of the box.
//					float_tt phase = kx * xPos + kz * zPos;
//					temp[ind3d] = complex_tt(f * cos(phase), f * sin(phase)); // *zScale
//					//          temp[ind3d][1] = f*sin(phase);        // *zScale
//					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
//				}
//			}
//		} // for iz....
//
//
//#if FLOAT_PRECISION ==1
//		fftwf_complex *ptr = (fftwf_complex *) &temp[0];
//		fftwf_plan plan = fftwf_plan_dft_2d(nz, nx, ptr, ptr, FFTW_BACKWARD,FFTW_ESTIMATE);
//		fftwf_execute(plan);
//		fftwf_destroy_plan(plan);
//#else
//		fftw_complex *ptr = (fftw_complex *) &temp[0];
//		fftw_plan plan = fftw_plan_dft_2d(nz, nx, ptr, ptr, FFTW_BACKWARD,
//				FFTW_ESTIMATE);
//		fftw_execute(plan);
//		fftw_destroy_plan(plan);
//#endif
//
//		// We also make sure that the potential touches zero at least somewhere. This will avoid
//		// sharp edges that could produce ringing artifacts.
//		// It is certainly debatable whether this is a good apprach, or not.
//		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
//		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
//		for (unsigned ix = 0; ix < nx / 2; ix++)
//			for (unsigned iz = 0; iz < nz / 2; iz++) {
//				unsigned ind3d = ix + iz * nx / 2;
//				float_tt zScale = 0;
//				// Integrate over nzPerSlice neighboring layers here:::::::::
//				for (int iiz = -izOffset; iiz <= izOffset; iiz++) {
//					if (iz + izOffset + iiz < nz / 2)
//						zScale += temp[ix + (iz + izOffset + iiz) * nx].real();
//				}
//				if (zScale < 0)
//					zScale = 0;
//				// assign the iz-th slice the sum of the 3 other slices:
//				// and divide by unit cell volume (since this is in 3D):
//				// Where does the '1/2' come from??? OVERSAMPLING*OVERSAMP_Y/8 = 1/2
//				// if nothing has changed, then OVERSAMPLING=2 OVERSAMP_Z=18.
//				// remember, that s=0.5*k;
//				// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
//				//    implicitly sets imaginary part to 0.
//				_atPotOffset[Znum][ind3d] = 47.8658 * dkx * dkz / (nz) * zScale;
//			}
//		if (_c->Output->LogLevel < 2)
//			BOOST_LOG_TRIVIAL(info)<<
//			format("Created 3D (r-z) %d x %d potential offset array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)")
//			%(nx / 2)%( nz / 2)% Znum% B% dkx% dky% dkz% izOffset;
//		} // end of creating atom if not exist...
//	Nz_lut = nz / 2;
//	nzSub = nzPerSlice;
//	Nr = nx / 2;
//	output = _atPotOffset[Znum];
}

}// end namespace QSTEM

