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

#include "crystal.hpp"
#include "stemtypes_fftw3.hpp"
#include "matrixlib.hpp"

#include "random.hpp"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <nlopt.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/array.hpp>

using boost::format;
using namespace nlopt;

#define PI180 1.7453292519943e-2
#define THZ_AMU_HBAR 0.15745702964189 /* A°^2*THz*amu/(hbar) */
// 4.46677327584453 /* 1e10/sqrt(THz*amu/(pi*hbar)) */
#define THZ_HBAR_KB 1.90963567802059 /* THz*hbar/kB */

namespace QSTEM {

static const float_tt k_wobScale = 1.0 / (8 * M_PI * M_PI);
static const float_tt k_sq3 = 1.0 / sqrt(3.0); /* sq3 is an additional needed factor which stems from
 * int_-infty^infty exp(-x^2/2) x^2 dx = sqrt(pi)
 * introduced in order to match the wobble factor with <u^2>
 */

CrystalBuilder::CrystalBuilder(StructureReaderPtr r, cStructureConfPtr sc, cModelConfPtr mc, cOutputConfPtr oc) :
		_minX(0), _maxX(0), _minY(0), _maxY(0), _minZ(0), _maxZ(0), _offsetX(0), _offsetY(0), m_adjustCubeSize(false), m_phononFile(
				boost::filesystem::path()), m_ax(0), _superCellBox(new superCellBox()), _baseAtoms(std::vector<atom>()), _atoms(std::vector<atom>()), _Mm(
				FloatArray2D(boost::extents[3][3])), m_MmInv(FloatArray2D(boost::extents[3][3])), IStructureBuilder(r, sc, mc, oc) {
	m_wobble_temp_scale = sqrt(_sc->temperatureK / 300.0);
	_tiltX = sc->crystalTiltX;
	_tiltY = sc->crystalTiltY;
	_tiltZ = sc->crystalTiltZ;
}
//(const std::vector<float_tt> &x, std::vector<float_tt> &grad, void* f_data)
double rotationCostFunc(const std::vector<double>& rotAngles, std::vector<double> &grad, void* f_data) {
	zoneAxisOptParams* p = reinterpret_cast<zoneAxisOptParams*>(f_data);
	std::vector<int> zone = p->zone;
	std::vector<int> refZone = p->refZone;
	armamat M = p->M;
	//	std::vector<float_tt> zone,std::vector<float_tt> refZone
	float_tt phi_x = rotAngles[0] * PI / 180;
	float_tt phi_y = rotAngles[1] * PI / 180;
	float_tt phi_z = rotAngles[2] * PI / 180;
	float_tt cx = cos(phi_x), sx = sin(phi_x), cy = cos(phi_y), sy = sin(phi_y), cz = cos(phi_z), sz = sin(phi_z);
	armamat Mx = {1,0,0,0,cx,-sx,0,sx,cx};
	armamat My = {cy,0,-sy,0,1,0,sy,0,cy};
	armamat Mz = {cz,sz,0,-sz,cz,0,0,0,1};
	BOOST_LOG_TRIVIAL(info)<< format("angles %f %f %f") % rotAngles[0]% rotAngles[1]% rotAngles[2];

	Mx.reshape(3, 3);
	My.reshape(3, 3);
	Mz.reshape(3, 3);
	//	std::cout << "Mx: "<< endl<< Mx << std::endl;
	//	std::cout << "My: "<< endl<< My << std::endl;
	//	std::cout << "Mz: "<< endl<< Mz << std::endl;
	//	std::cout << "M: "<< endl<< M << std::endl;

	armamat Mrot = Mx * My * Mz;
	//	std::cout << "Mrot: "<< endl<< Mrot << std::endl;
	armavec zone1 = {(float_tt)zone[0],(float_tt)zone[1],(float_tt)zone[2]};
	armavec refZone1 = {(float_tt)refZone[0],(float_tt)refZone[1],(float_tt)refZone[2]};
	armavec v = Mrot * M * zone1;

	v = v / sqrt(arma::sum(v % v));
	//	std::cout << "refZone1: "<< endl<< refZone1 << std::endl;
	//	std::cout << "v: "<< endl<< v << std::endl;
	auto tmp = v - refZone1;
	//	std::cout << "tmp: "<< endl<< tmp << endl;
	auto tmp2 = tmp % tmp;
	//	std::cout << "tmp2: "<< endl<< tmp2 << endl;
	float_tt chi2 = arma::sum(tmp2);
	//	std::cout << "chi2: "<< chi2 << endl;
	if (refZone[1] == 0) {
		armavec vtmp = {1,0,0};
		armavec vx = Mrot*M*vtmp;
		chi2 += vx[2]*vx[2];
	}
	//	std::cout << "chi2: "<< chi2 << endl;
	return chi2;
	//	Mx[0] ={1,0,0};
	//	    cx = cos(phi_x);    sx = sin(phi_x);
	//	    cy = cos(phi_y);    sy = sin(phi_y);
	//	    cz = cos(phi_z);    sz = sin(phi_z);
	//
	//	    Mx = [1 0 0;0 cx -sx;0 sx cx];
	//	    My = [cy 0 sy;0 1 0; -sy 0 cy];
	//	    Mz = [cz -sz 0;sz cz 0; 0 0 1];
	//	    Mrot = Mz*My*Mx;
	//	v = Mrot*M*zone;
	//	v = v/sqrt(sum(v.^2));
	//	% [v refZone]
	//	chi2 = sum((v-refZone).^2);
	//	% check also the direction of the x-vector:
	//	% It should point to the right.
	//	if (refZone(1) == 0)
	//	    vx = Mrot*M*[1;0;0];
	//	    chi2 = chi2+vx(2)^2;
	//	end
}

void CrystalBuilder::CalculateTiltAnglesFromZoneAxis() {
	auto refZone = _sc->zoneAxis;
	std::vector<int> refZone2(3), zone2(3);
	zone2[0] = 0;
	zone2[1] = 0;
	zone2[2] = 1;

	int len = sqrt(refZone[0] * refZone[0] + refZone[1] * refZone[1] + refZone[2] * refZone[2]);
	for (int i = 0; i < 3; i++) {
		refZone2[i] = refZone[i] / len;
	}

	armamat M = {_Mm[0][0],_Mm[1][0],_Mm[2][0],
		_Mm[0][1],_Mm[1][1],_Mm[2][1],
		_Mm[0][2],_Mm[1][2],_Mm[2][2]};
	M.reshape(3, 3);

	zoneAxisOptParams data(M, zone2, refZone2);

	opt o = opt(LN_NELDERMEAD, 3);
	o.set_lower_bounds(-180);
	o.set_upper_bounds(180);
	o.set_min_objective(&rotationCostFunc, (void*) &data);
	o.set_stopval(1e-15);
	double f_final;
	std::vector<double> x(3);
	x[0] = x[1] = x[2] = 0;

	o.optimize(x, f_final);

	BOOST_LOG_TRIVIAL(info)<< format("Rotating (%g,%g,%g) degrees from zone axis [0,0,1] to [%d,%d,%d]")
	% x[0]%x[1]%x[2]%refZone[0]%refZone[1]%refZone[2];

	_tiltX = x[0];
	_tiltY = x[1];
	_tiltZ = x[2];

	//	if nargin < 3
	//	    refZone = [0; 0; 1];
	//	end
	//	if nargin < 2
	//	    M = eye(3);
	//	end
	//	if isempty(M)
	//	    M = eye(3);
	//	end
	//
	//
	//	% The rotations are in the following sequence: X-Y-Z
	//	% The rotation about the z-axis will be last.
	//	% The rotation angles must tilt the vector 'zone' to refZone
	//	[Ny,Nx] = size(zone);
	//	if (Ny == 1), zone = zone.'; end
	//	% zone = M*zone;
	//	zone2 = zone/sqrt(sum(zone.^2));
	//
	//	[Ny,Nx] = size(refZone);
	//	if (Ny == 1), refZone = refZone.'; end
	//	% refZone = M*refZone;
	//	refZone2 = refZone/sqrt(sum(refZone.^2));
	//
	//
	//
	//	rotAngles = 20*rand(1,3);
	//
	//	[rotAngles,chi2] = fminsearch(@(rotAngles) rotationCostFunc(rotAngles,M,zone2,refZone2), rotAngles);
	//	% [rotAngles,chi2] = ga(@(rotAngles) rotationCostFunc(rotAngles,zone),3);
	//	% rotAngles(3) = 0;
	//	[chi2,Mrot] = rotationCostFunc(rotAngles,M,zone2,refZone2);
	//	fprintf('chi2 = %f\n',chi2);
	//	% Mrot*M
	//
	//	function [chi2,Mrot] = rotationCostFunc(rotAngles,M,zone,refZone)
	//	% Generate the rotation matrix with which one can rotate any atom position:
	//	Mrot = zeros(3);
	//	phi_x = rotAngles(1)*pi/180;
	//	phi_y = rotAngles(2)*pi/180;
	//	phi_z = rotAngles(3)*pi/180;
	//
	//	if (1)
	//	    cx = cos(phi_x);    sx = sin(phi_x);
	//	    cy = cos(phi_y);    sy = sin(phi_y);
	//	    cz = cos(phi_z);    sz = sin(phi_z);
	//
	//	    Mx = [1 0 0;0 cx -sx;0 sx cx];
	//	    My = [cy 0 sy;0 1 0; -sy 0 cy];
	//	    Mz = [cz -sz 0;sz cz 0; 0 0 1];
	//	    Mrot = Mz*My*Mx;
	//	else
	//	    Mrot(1,1) = cos(phi_z)*cos(phi_y);
	//	    Mrot(1,2) = cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x);
	//	    Mrot(1,3) = cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x);
	//
	//	    Mrot(2,1) = sin(phi_z)*cos(phi_y);
	//	    Mrot(2,2) = sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x);
	//	    Mrot(2,3) = sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x);
	//
	//	    Mrot(3,1) = -sin(phi_y);
	//	    Mrot(3,2) = cos(phi_y)*sin(phi_x);
	//	    Mrot(3,3) = cos(phi_y)*cos(phi_x);
	//	end
	//
	//	v = Mrot*M*zone;
	//	v = v/sqrt(sum(v.^2));
	//	% [v refZone]
	//	chi2 = sum((v-refZone).^2);
	//	% check also the direction of the x-vector:
	//	% It should point to the right.
	//	if (refZone(1) == 0)
	//	    vx = Mrot*M*[1;0;0];
	//	    chi2 = chi2+vx(2)^2;
	//	end
}
// Calculates only cell parameters
superCellBoxPtr CrystalBuilder::DisplaceAtoms() {
	//	std::vector<atom>::iterator at = _atoms.begin(), end = _atoms.end();
	//	FloatArray2D u(boost::extents[1][3]);
	//
	//	for (at; at != end; ++at) {
	//		EinsteinDisplacement(u, (*at));
	//		// Add obtained u to atomic position
	//		(*at).r[0] += u[0][0];
	//		(*at).r[1] += u[0][1];
	//		(*at).r[2] += u[0][2];
	//	}
	return _superCellBox;
}
superCellBoxPtr CrystalBuilder::Build() {
	bool handleVacancies = true;
	int jz, i2, j, ix, iy, iz = 0;
	float_tt boxXmin = 0, boxXmax = 0, boxYmin = 0, boxYmax = 0, boxZmin = 0, boxZmax = 0;
	float_tt boxCenterX, boxCenterY, boxCenterZ, boxCenterXrot, boxCenterYrot, boxCenterZrot, bcX, bcY, bcZ;
	float_tt totOcc, choice, lastOcc;
	float_tt u[3];
	static int ncoord_old = 0;
	int ncx = _sc->nCellX;
	int ncy = _sc->nCellY;
	int ncz = _sc->nCellZ;
	bool crystalIsTilted = (_tiltX != 0) || (_tiltY != 0) || (_tiltX != 0);

	ReadFromFile();

	if (handleVacancies) {
		qsort(&_baseAtoms[0], _baseAtoms.size(), sizeof(atom), &CrystalBuilder::AtomCompareZYX);
	}
	for (int i = _baseAtoms.size() - 1; i >= 0; i--) {
		if (_mc->UseTDS) {
			m_u2[_baseAtoms[i].Znum] = 0;
		}
	}

	if (_sc->rotateToZoneAxis) {
		CalculateTiltAnglesFromZoneAxis();
	}

	if (_sc->isBoxed) {
		TiltBoxed(_baseAtoms.size(), handleVacancies);
		boxCenterX = _sc->boxX;
		boxCenterY = _sc->boxY;
		boxCenterZ = _sc->boxZ;
		boxXmax = _sc->boxX;
		boxYmax = _sc->boxY;
		boxZmax = _sc->boxZ;
	} else {
		ReplicateUnitCell(handleVacancies);

		BOOST_LOG_TRIVIAL(trace)<< "Atoms after replication of unit cells";
		BOOST_LOG_TRIVIAL(trace)<< format("size: %d") % _atoms.size();
#pragma omp parallel for
		for (int j = 0; j < _atoms.size(); j++) {
			BOOST_LOG_TRIVIAL(trace)<< format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % _atoms[j].r[0] % _atoms[j].r[1] % _atoms[j].r[2];
			// This converts also to cartesian coordinates
			float_tt x = _Mm[0][0] * _atoms[j].r[0] + _Mm[1][0] * _atoms[j].r[1] + _Mm[2][0] * _atoms[j].r[2];
			float_tt y = _Mm[0][1] * _atoms[j].r[0] + _Mm[1][1] * _atoms[j].r[1] + _Mm[2][1] * _atoms[j].r[2];
			float_tt z = _Mm[0][2] * _atoms[j].r[0] + _Mm[1][2] * _atoms[j].r[1] + _Mm[2][2] * _atoms[j].r[2];

			_atoms[j].r[0] = x;
			_atoms[j].r[1] = y;
			_atoms[j].r[2] = z;

			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % _atoms[j].r[0] % _atoms[j].r[1] % _atoms[j].r[2];
		}

		//  tilt around the center of the full crystal
		bcX = ncx / 2.0;
		bcY = ncy / 2.0;
		bcZ = ncz / 2.0;
		u[0] = _Mm[0][0] * bcX + _Mm[1][0] * bcY + _Mm[2][0] * bcZ;
		u[1] = _Mm[0][1] * bcX + _Mm[1][1] * bcY + _Mm[2][1] * bcZ;
		u[2] = _Mm[0][2] * bcX + _Mm[1][2] * bcY + _Mm[2][2] * bcZ;
		boxCenterX = u[0];
		boxCenterY = u[1];
		boxCenterZ = u[2];
		// Determine the size of the (rotated) super cell
		for (int icx = 0; icx <= ncx; icx += ncx) {
			for (int icy = 0; icy <= ncy; icy += ncy) {
				for (int icz = 0; icz <= ncz; icz += ncz) {
					u[0] = _Mm[0][0] * (icx - bcX) + _Mm[1][0] * (icy - bcY) + _Mm[2][0] * (icz - bcZ);
					u[1] = _Mm[0][1] * (icx - bcX) + _Mm[1][1] * (icy - bcY) + _Mm[2][1] * (icz - bcZ);
					u[2] = _Mm[0][2] * (icx - bcX) + _Mm[1][2] * (icy - bcY) + _Mm[2][2] * (icz - bcZ);
					if (crystalIsTilted)
						RotateVect(u, u, _tiltX, _tiltY, _tiltZ); // simply applies rotation matrix
					// x = u[0]+boxCenterXrot; y = u[1]+boxCenterYrot; z = u[2]+boxCenterZrot;
					float_tt x = u[0] + boxCenterX;
					float_tt y = u[1] + boxCenterY;
					float_tt z = u[2] + boxCenterZ;
					if ((icx == 0) && (icy == 0) && (icz == 0)) {
						boxXmin = boxXmax = x;
						boxYmin = boxYmax = y;
						boxZmin = boxZmax = z;
					} else {
						boxXmin = boxXmin > x ? x : boxXmin;
						boxXmax = boxXmax < x ? x : boxXmax;
						boxYmin = boxYmin > y ? y : boxYmin;
						boxYmax = boxYmax < y ? y : boxYmax;
						boxZmin = boxZmin > z ? z : boxZmin;
						boxZmax = boxZmax < z ? z : boxZmax;
					}
				}
			}
		}

#pragma omp parallel
		for (int j = 0; j < _atoms.size(); j++) {
			// Tilt
			if (crystalIsTilted) {
				u[0] = _atoms[j].r[0] - boxCenterX;
				u[1] = _atoms[j].r[1] - boxCenterY;
				u[2] = _atoms[j].r[2] - boxCenterZ;
				RotateVect(u, u, _tiltX, _tiltY, _tiltZ); // simply applies rotation matrix
				u[0] += boxCenterX;
				u[1] += boxCenterY;
				u[2] += boxCenterZ;
				_atoms[j].r[0] = u[0];
				_atoms[j].r[1] = u[1];
				_atoms[j].r[2] = u[2];
			} /* if tilts != 0 ... */

			// Rebase to topleft
			_atoms[j].r[0] -= boxXmin;
			_atoms[j].r[1] -= boxYmin;
			_atoms[j].r[2] -= boxZmin;
			BOOST_LOG_TRIVIAL(trace)<< format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % _atoms[j].r[0] % _atoms[j].r[1] % _atoms[j].r[2];

			//apply offsets
			if (_mc->HasOffset()) {
				_atoms[j].r[0] += _mc->xOffset;
				_atoms[j].r[1] += _mc->yOffset;
				_atoms[j].r[2] += _mc->zOffset;
			}
		}
	}

	for (int j = 0; j < _atoms.size(); j++) {
		_xyzPos.push_back(_atoms[j].r[0]);
		_xyzPos.push_back(_atoms[j].r[1]);
		_xyzPos.push_back(_atoms[j].r[2]);
		_znums.push_back(_atoms[j].Znum);
		_occupancy.push_back(_atoms[j].occ);
	}

	_sizeX = boxXmax - boxXmin;
	_sizeY = boxYmax - boxYmin;
	_sizeZ = boxZmax - boxZmin;
	CalculateCrystalBoundaries();
	_superCellBox->atoms = _atoms;
	_superCellBox->ax = _sizeX;
	_superCellBox->by = _sizeY;
	_superCellBox->cz = _sizeZ;
	_superCellBox->cmx = boxCenterX;
	_superCellBox->cmy = boxCenterY;
	_superCellBox->cmz = boxCenterZ;
	_superCellBox->uniqueZ = _uniqueAtoms;
	_superCellBox->xyzPos = _xyzPos;
	_superCellBox->znums = _znums;
	_superCellBox->occupancy = _occupancy;
	return _superCellBox;
}

CrystalBuilder::~CrystalBuilder() {
}
void CrystalBuilder::SetFillUnitCell(bool val) {
	_fillUnitCell = val;
}
void CrystalBuilder::ReadFromFile() {
	_r->ReadCellParams(_Mm);
	_r->ReadAtoms(_baseAtoms, _uniqueAtoms, true);
	CalculateCellDimensions();
	int s = _baseAtoms.size();
//	BOOST_LOG_TRIVIAL(info)<< format("Read %d atoms, tds: %d") % s % (_mc->UseTDS);
}
void CrystalBuilder::CalculateCrystalBoundaries() {
	int largestZindex = -1;
#pragma omp parallel for
	for (int i = 0; i < _atoms.size(); i++) {
#pragma omp critical
		{
			if (_atoms[i].r[0] < _minX)
				_minX = _atoms[i].r[0];
		}
#pragma omp critical
		{
			if (_atoms[i].r[0] > _maxX)
				_maxX = _atoms[i].r[0];
		}
#pragma omp critical
		{
			if (_atoms[i].r[1] < _minY)
				_minY = _atoms[i].r[1];
		}
#pragma omp critical
		{
			if (_atoms[i].r[1] > _maxY)
				_maxY = _atoms[i].r[1];
		}
#pragma omp critical
		{
			if (_atoms[i].r[2] < _minZ)
				_minZ = _atoms[i].r[2];
		}
#pragma omp critical
		{
			if (_atoms[i].r[2] > _maxZ) {
				_maxZ = _atoms[i].r[2];
				largestZindex = i;
			}
		}
	}
	BOOST_LOG_TRIVIAL(trace)<<format("largest Z at index %d : (%f,%f,%f)") %largestZindex%_atoms[largestZindex].r[0]%_atoms[largestZindex].r[1]%_atoms[largestZindex].r[2];
}

void CrystalBuilder::Init(unsigned run_number) {
	BOOST_LOG_TRIVIAL(trace)<<format("range of thermally displaced atoms (%d atoms): ") %_atoms.size();
	BOOST_LOG_TRIVIAL(trace)<<format("X: %g .. %g")% _minX% _maxX;
	BOOST_LOG_TRIVIAL(trace)<<format("Y: %g .. %g")% _minY% _maxY;
	BOOST_LOG_TRIVIAL(trace)<<format("Z: %g .. %g")% _minZ% _maxZ;

	qsort(&_atoms[0], _atoms.size(), sizeof(atom),&CrystalBuilder::AtomCompareZnum);
	WriteStructure(run_number);
}

void CrystalBuilder::DisplayParams() {
	BOOST_LOG_TRIVIAL(info)<<
	"***************************** Structure Parameters ***********************************************";
	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info)<<format("* Input file:           %s") % _sc->structureFilename.string().c_str();

	if (_sc->isBoxed == false)
	BOOST_LOG_TRIVIAL(info)<<format(" Unit cell:            ax=%g by=%g cz=%g")% m_ax% m_by% m_cz;
	else {
		BOOST_LOG_TRIVIAL(info)<<format(" Size of Cube:         ax=%g by=%g cz=%g")
		% _sc->boxX% _sc->boxY% _sc->boxZ;
		BOOST_LOG_TRIVIAL(info)<<format(" Cube size adjusted:   %s")%( m_adjustCubeSize ? "yes" : "no");
	}
	BOOST_LOG_TRIVIAL(info)<<format(" Super cell:           %d x %d x %d unit cells")
	%_sc->nCellX% _sc->nCellY% _sc->nCellZ;
	BOOST_LOG_TRIVIAL(info)<<format(" Number of atoms:      %d (super cell)")% _atoms.size();
	BOOST_LOG_TRIVIAL(info)<<format(" Crystal tilt:         x=%g deg, y=%g deg, z=%g deg")
	%(_tiltX)% (_tiltY )%( _tiltZ );
	BOOST_LOG_TRIVIAL(info)<<format(" Model dimensions:     ax=%gA, by=%gA, cz=%gA (after tilt)")
	%m_ax% m_by% m_cz;
	BOOST_LOG_TRIVIAL(info)<<format(" Temperature:          %gK") %_sc->temperatureK;
	if (_mc->UseTDS)
	BOOST_LOG_TRIVIAL(info)<<format(" TDS:                  yes");
	else
	BOOST_LOG_TRIVIAL(info)<<format(" TDS:                  no");
}
void CrystalBuilder::SetCellParameters(float_tt ax, float_tt by, float_tt cz) {
	m_ax = ax;
	m_by = by;
	m_cz = cz;
}

void CrystalBuilder::GetCellAngles(float_tt &alpha, float_tt &beta, float_tt &gamma) {
	alpha = m_cAlpha;
	beta = m_cBeta;
	gamma = m_cGamma;
}

void CrystalBuilder::GetCellParameters(float_tt &ax, float_tt &by, float_tt &cz) {
	ax = m_ax;
	by = m_by;
	cz = m_cz;
}

void CrystalBuilder::OffsetCenter(atom &center) {
	/*
	 define the center of our unit cell by moving the atom specified
	 by "center" at position (0.5,0.5,0.0)
	 */

	float_tt dx = m_ax / 2.0f - center.r[0];
	float_tt dy = m_by / 2.0f - center.r[1];
	float_tt dz = -center.r[2];
	for (size_t i = 0; i < _atoms.size(); i++) {
		_atoms[i].r[0] += dx;
		if (_atoms[i].r[0] < 0.0f)
			_atoms[i].r[0] += m_ax;
		else if (_atoms[i].r[0] > m_ax)
			_atoms[i].r[0] -= m_ax;
		_atoms[i].r[1] += dy;
		if (_atoms[i].r[1] < 0.0f)
			_atoms[i].r[1] += m_by;
		else if (_atoms[i].r[1] > m_by)
			_atoms[i].r[1] -= m_by;
		_atoms[i].r[2] += dz;
		if (_atoms[i].r[2] < 0.0f)
			_atoms[i].r[2] += m_cz;
		else if (_atoms[i].r[2] > m_cz)
			_atoms[i].r[2] -= m_cz;
	}
}

// Uses m_Mm to calculate ax, by, cz, and alpha, beta, gamma
void CrystalBuilder::CalculateCellDimensions() {
	m_ax = sqrt(_Mm[0][0] * _Mm[0][0] + _Mm[0][1] * _Mm[0][1] + _Mm[0][2] * _Mm[0][2]);
	m_by = sqrt(_Mm[1][0] * _Mm[1][0] + _Mm[1][1] * _Mm[1][1] + _Mm[1][2] * _Mm[1][2]);
	m_cz = sqrt(_Mm[2][0] * _Mm[2][0] + _Mm[2][1] * _Mm[2][1] + _Mm[2][2] * _Mm[2][2]);
	m_cGamma = atan2(_Mm[1][1], _Mm[1][0]);
	m_cBeta = acos(_Mm[2][0] / m_cz);
	m_cAlpha = acos(_Mm[2][1] * sin(m_cGamma) / m_cz + cos(m_cBeta) * cos(m_cGamma));
	m_cGamma /= (float) PI180;
	m_cBeta /= (float) PI180;
	m_cAlpha /= (float) PI180;
}

// TODO: old version return a pointer to new atom positions.  Can we do this in place?
//		If not, just set atoms to the new vector.
void CrystalBuilder::TiltBoxed(int ncoord, bool handleVacancies) {
	int atomKinds = 0, i, jVac, jequal, jChoice, i2, ix, iy, iz, atomSize;
	static float_tt *axCell, *byCell, *czCell = NULL;
	FloatArray2D a(boost::extents[3][1]), aOrig(boost::extents[3][1]), b(boost::extents[3][1]);
	float_tt x, y, z, dx = 0, dy = 0, dz = 0;
	float_tt totOcc, lastOcc, choice;

	int Ncells = _sc->nCellX * _sc->nCellY * _sc->nCellZ;
	FloatArray2D u(boost::extents[1][3]), uf(boost::extents[1][3]);
	FloatArray2D Mm(boost::extents[3][3]), Mminv(boost::extents[3][3]), MpRed(boost::extents[3][3]), // conversion lattice to obtain red. prim. coords/ from reduced cubic rect.
	MpRedInv(boost::extents[3][3]), //conversion lattice to obtain red. cub. coords from reduced primitive lattice coords
	MbPrim(boost::extents[3][3]), // float_tt version of primitive lattice basis
	MbPrimInv(boost::extents[3][3]), // float_tt version of inverse primitive lattice basis
	MmOrig(boost::extents[3][3]), MmOrigInv(boost::extents[3][3]);
	dx = _offsetX;
	dy = _offsetY;
	// We need to copy the transpose of m_Mm to Mm.
	for (ix = 0; ix < 3; ix++)
		for (iy = 0; iy < 3; iy++)
			Mm[ix][iy] = _Mm[iy][ix];

	memcpy(MmOrig.data(), Mm.data(), 3 * 3 * sizeof(float_tt));
	Inverse_3x3(MmOrigInv, MmOrig);
	/* remember that the angles are in rad: */
	RotateMatrix(Mm, Mm, _tiltX / PI180, _tiltY / PI180, _tiltZ / PI180);
	Inverse_3x3(Mminv, Mm);  // computes Mminv from Mm!
	/* find out how far we will have to go in unit of unit cell vectors.  when creating the supercell by
	 * checking the number of unit cell vectors necessary to reach every corner of the supercell box.
	 */
	MatrixProduct(Mminv, 3, 3, a, 3, 1, b);
	showMatrix(_Mm, 3, 3, "_M");
	showMatrix(Mm, 3, 3, "M");
	// showMatrix(Mminv,3,3,"M");
	unsigned nxmax, nymax, nzmax;
	unsigned nxmin = nxmax = (int) floor(b[0][0] - dx);
	unsigned nymin = nymax = (int) floor(b[1][0] - dy);
	unsigned nzmin = nzmax = (int) floor(b[2][0] - dz);
	for (ix = 0; ix <= 1; ix++)
		for (iy = 0; iy <= 1; iy++)
			for (iz = 0; iz <= 1; iz++) {
				a[0][0] = ix * _sc->boxX - dx;
				a[1][0] = iy * _sc->boxY - dy;
				a[2][0] = iz * _sc->boxZ - dz;

				MatrixProduct(Mminv, 3, 3, a, 3, 1, b);

				if (nxmin > (int) floor(b[0][0]))
					nxmin = (int) floor(b[0][0]);
				if (nxmax < (int) ceil(b[0][0]))
					nxmax = (int) ceil(b[0][0]);
				if (nymin > (int) floor(b[1][0]))
					nymin = (int) floor(b[1][0]);
				if (nymax < (int) ceil(b[1][0]))
					nymax = (int) ceil(b[1][0]);
				if (nzmin > (int) floor(b[2][0]))
					nzmin = (int) floor(b[2][0]);
				if (nzmax < (int) ceil(b[2][0]))
					nzmax = (int) ceil(b[2][0]);
			}
	jVac = 0;
	for (i = 0; i < ncoord;) {
		atom at(_baseAtoms[i]);
		if ((handleVacancies) && (at.Znum > 0)) {
			totOcc = at.occ;
			for (jequal = i + 1; jequal < ncoord; jequal++) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase the total occupany and the counter jequal.
				if ((fabs(at.r[0] - _baseAtoms[jequal].r[0]) < 1e-6) && (fabs(at.r[1] - _baseAtoms[jequal].r[1]) < 1e-6)
						&& (fabs(at.r[2] - _baseAtoms[jequal].r[2]) < 1e-6)) {
					totOcc += _baseAtoms[jequal].occ;
				} else
					break;
			}
		} else {
			jequal = i + 1;
			totOcc = 1;
		}
		for (ix = nxmin; ix <= nxmax; ix++) {
			for (iy = nymin; iy <= nymax; iy++) {
				for (iz = nzmin; iz <= nzmax; iz++) {
					atom newAtom(_baseAtoms[i]);
					u[0][0] = u[0][1] = u[0][2] = 0;
					// atom position in cubic reduced coordinates:
					aOrig[0][0] = ix + at.r[0];
					aOrig[1][0] = iy + at.r[1];
					aOrig[2][0] = iz + at.r[2];

//					BOOST_LOG_TRIVIAL(trace) << format("aOrig=(%3.3f,%3.3f,%3.3f)") % aOrig[0][0] % aOrig[1][0] % aOrig[2][0];

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					// All we need to decide is whether to include the atom at all (if totOcc < 1
					// of which of the atoms at equal positions to include
					jChoice = i;  // This will be the atom we wil use.
					if ((totOcc < 1) || (jequal > i + 1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values.
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totOcc < 1.0)
							choice = ran1();
						else
							choice = totOcc * ran1();
						lastOcc = 0;
						for (i2 = i; i2 < jequal; i2++) {
							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc)
							if ((choice < lastOcc) || (choice >= lastOcc + _baseAtoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								// atoms[atomCount].Znum =  0;  // vacancy
								jVac++;
							} else {
								jChoice = i2;
							}
							lastOcc += _baseAtoms[i2].occ;
						}
					}
					if (_mc->UseTDS)
						switch (_mc->displacementType) {
							case DisplacementType::Einstein:
								EinsteinDisplacement(u, _baseAtoms[i]);
								break;
							case DisplacementType::Phonon:
								PhononDisplacement(u, jChoice, ix, iy, iz, _atoms[jChoice], false);
								break;
							case DisplacementType::None:
							default:
								break;
						}

					a[0][0] = aOrig[0][0] + u[0][0];
					a[1][0] = aOrig[1][0] + u[0][1];
					a[2][0] = aOrig[2][0] + u[0][2];

					MatrixProduct(Mm, 3, 3, aOrig, 3, 1, b);

					// b now contains atom positions in cartesian coordinates */
					x = b[0][0] + dx;
					y = b[1][0] + dy;
					z = b[2][0] + dz;

					bool atomIsInBox = (x >= 0) && (x <= _sc->boxX) && (y >= 0) && (y <= _sc->boxY) && (z >= 0) && (z <= _sc->boxZ);
					if (atomIsInBox) {
						newAtom.r =
						armavec( {(float_tt)x,(float_tt)y,(float_tt)z});
						newAtom.dw = _baseAtoms[jChoice].dw;
						newAtom.occ = _baseAtoms[jChoice].occ;
						newAtom.q = _baseAtoms[jChoice].q;
						newAtom.Znum = _baseAtoms[jChoice].Znum;
						_atoms.push_back(newAtom);
						BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % _atoms.size() % newAtom.r[0] % newAtom.r[1] % newAtom.r[2];
					}
				}
			}
		}
		i = jequal;
	}
	BOOST_LOG_TRIVIAL(trace)<<format("Removed %d atoms because of multiple occupancy or occupancy < 1") % jVac;
	m_ax = _sc->boxX;
	m_by = _sc->boxY;
	m_cz = _sc->boxZ;
}
////////////////////////////////////////////////////////////////////////
// replicateUnitCell
// 
// Replicates the unit cell NcellX x NCellY x NCellZ times
// applies phonon displacement and removes vacancies and atoms appearing
// on same position:
// ncoord is the number of atom positions that has already been read.
// memory for the whole atom-array of size natom has already been allocated
// but the sites beyond natom are still empty.
void CrystalBuilder::ReplicateUnitCell(int handleVacancies) {
	int i, j, i2, jChoice, ncx, ncy, ncz, icx, icy, icz;
	int jequal; // Number of atoms that share a position (mixed position if > 1)
	int jVac;  // Number of vacancies total in replicated supercell
	int jCell; // Offset (in number of coordinates) from origin cell
	int atomKinds = 0;
	float_tt totalOccupancy;
	float_tt choice, lastOcc;
	FloatArray2D u(boost::extents[1][3]);
	ncx = _sc->nCellX;
	ncy = _sc->nCellY;
	ncz = _sc->nCellZ;

	_atoms.resize(ncx * ncy * ncz * _baseAtoms.size());

	//////////////////////////////////////////////////////////////////////////////
	// Look for atoms which share the same position:
	jVac = 0;  // no atoms have been removed yet
	for (i = _baseAtoms.size() - 1; i >= 0;) {

		////////////////
		if ((handleVacancies) && (_atoms[i].Znum > 0)) {
			totalOccupancy = _atoms[i].occ;
			for (jequal = i - 1; jequal >= 0; jequal--) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase
				// the total occupany and the counter jequal.
				if ((fabs(_atoms[i].r[0] - _atoms[jequal].r[0]) < 1e-6) && (fabs(_atoms[i].r[1] - _atoms[jequal].r[1]) < 1e-6)
						&& (fabs(_atoms[i].r[2] - _atoms[jequal].r[2]) < 1e-6)) {
					totalOccupancy += _atoms[jequal].occ;

				} else
					break;

			} // jequal-loop
		} else {
			jequal = i - 1;
			totalOccupancy = 1;
		}

		////////////////
		/* replicate unit cell ncx,y,z times: */
		/* We have to start with the last atoms first, because once we added the displacements
		 * to the original unit cell (icx=icy=icz=0), we cannot use those positions
		 * as unit cell coordinates for the other atoms anymore
		 */
		for (icx = ncx - 1; icx >= 0; icx--) {
			for (icy = ncy - 1; icy >= 0; icy--) {
				for (icz = ncz - 1; icz >= 0; icz--) {
					u[0][0] = u[0][1] = u[0][2] = 0;
					jCell = (icz + icy * ncz + icx * ncy * ncz) * _baseAtoms.size();
					j = jCell + i;
					/* We will also add the phonon displacement to the atomic positions now: */
					_atoms[j].dw = _baseAtoms[i].dw;
					_atoms[j].occ = _baseAtoms[i].occ;
					_atoms[j].q = _baseAtoms[i].q;
					_atoms[j].Znum = _baseAtoms[i].Znum;

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					jChoice = i;
					if ((totalOccupancy < 1) || (jequal < i - 1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1
						//
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totalOccupancy < 1.0)
							choice = ran1();
						else
							choice = totalOccupancy * ran1();
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2 = i; i2 > jequal; i2--) {
							_atoms[jCell + i2].dw = _baseAtoms[i2].dw;
							_atoms[jCell + i2].occ = _baseAtoms[i2].occ;
							_atoms[jCell + i2].q = _baseAtoms[i2].q;
							_atoms[jCell + i2].Znum = _baseAtoms[i2].Znum;

							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc)
							if ((choice < lastOcc) || (choice >= lastOcc + _baseAtoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								_atoms[jCell + i2].Znum = 0;  // vacancy
								jVac++;
							} else {
								jChoice = i2;
							}
							lastOcc += _baseAtoms[i2].occ;
						}
					}

					if (_mc->UseTDS)
						switch (_mc->displacementType) {
							case DisplacementType::Einstein:
								EinsteinDisplacement(u, _baseAtoms[i]);
								break;
							case DisplacementType::Phonon:
								PhononDisplacement(u, jChoice, icx, icy, icz, _atoms[jChoice], false);
								break;
							case DisplacementType::None:
							default:
								break;
						}
					for (i2 = i; i2 > jequal; i2--) {
						float_tt x = _baseAtoms[i2].r[0] + icx + u[0][0];
						float_tt y = _baseAtoms[i2].r[1] + icy + u[0][1];
						float_tt z = _baseAtoms[i2].r[2] + icz + u[0][2];
						_atoms[jCell + i2].r = armavec( {x, y, z});
					}
				}
			}
		}
		i = jequal;
	}
	if ((jVac > 0))
		BOOST_LOG_TRIVIAL(trace)<<format("Removed %d atoms because of occupancies < 1 or multiple atoms in the same place") %jVac;
	}

void CrystalBuilder::EinsteinDisplacement(FloatArray2D& u, atom &_atom) {
	/* convert the Debye-Waller factor to sqrt(<u^2>) */
	float_tt wobble = m_wobble_temp_scale * sqrt(_atom.dw * k_wobScale);
	u[0][0] = (wobble * k_sq3 * gasdev());
	u[0][1] = (wobble * k_sq3 * gasdev());
	u[0][2] = (wobble * k_sq3 * gasdev());
	///////////////////////////////////////////////////////////////////////
	// Book keeping:
	m_u2[_atom.Znum] += u[0][0] * u[0][0] + u[0][1] * u[0][1] + u[0][2] * u[0][2];
	//ux += u[0]; uy += u[1]; uz += u[2];
	m_u2Count[_atom.Znum]++;

	/* Finally we must convert the displacement for this atom back into its fractional
	 * coordinates so that we can add it to the current position in vector a
	 */
	//	std::vector<float_tt> uf(3, 0);
	FloatArray2D uf(boost::extents[1][3]);
	MatrixProduct(u, 1, 3, m_MmInv, 3, 3, uf);
	u = uf;
	memcpy(u.data(), uf.data(), 3 * sizeof(float_tt));
}

/*******************************************************************************
 * int phononDisplacement:
 * This function will calculate the phonon displacement for a given atom i of the
 * unit cell, which has been replicated to the larger cell (icx,icy,icz)
 * The phonon displacement is either defined by the phonon data file, or,
 * the Einstein model, if the appropriate flags in the muls struct are set
 * The displacement will be given in fractional coordinates of a single unit cell.
 *
 * Input parameters:
 * Einstein-mode:
 * need only: dw, Znum, atomCount
 * atomCount: give statistics report, if 0, important only for non-Einstein mode
 * maxAtom: total number of atoms (will be called first, i.e. atomCount=maxAtoms-1:-1:0)
 *
 * Phonon-file mode:
 * ...
 *
 ********************************************************************************/
//  phononDisplacement(u,muls,jChoice,icx,icy,icz,j,atoms[jChoice].dw,*natom,jz);
//  j == atomCount
void CrystalBuilder::PhononDisplacement(FloatArray2D &u, int id, int icx, int icy, int icz, atom &atom, bool printReport) {
	int ix, iy, idd; // iz;
	static FILE *fpPhonon = NULL;
	static int Nk, Ns;  // number of k-vectors and atoms per primitive unit cell
	static float_tt *massPrim;   // masses for every atom in primitive basis
	static FloatArray2D omega;     // array of eigenvalues for every k-vector
	static complex_tt ***eigVecs;  // array of eigenvectors for every k-vector
	static FloatArray2D kVecs(boost::extents[Nk][3]); // array for Nk 3-dim k-vectors
	static FloatArray2D q1(boost::extents[3 * Ns][Nk]), q2(boost::extents[3 * Ns][Nk]);
	int ik, lambda, icoord; // Ncells, nkomega;
	float_tt kR, kRi, kRr, wobble;
	static float_tt *u2T, ux = 0, uy = 0, uz = 0; // u2Collect=0; // Ttotal=0;
	std::map<unsigned, float_tt> u2;
	std::map<unsigned, unsigned> u2Count;
	// static float_tt uxCollect=0,uyCollect=0,uzCollect=0;
	static int *u2CountT, runCount = 1, u2Size = -1;
	static long iseed = 0;
	FloatArray2D Mm(boost::extents[3][3]), MmInv(boost::extents[3][3]);
	// static float_tt **MmOrig=NULL,**MmOrigInv=NULL;
	static float_tt *axCell, *byCell, *czCell;
	static float_tt wobScale = 0, sq3, scale = 0;
	std::vector<float_tt> uf(3), b(3);

	float_tt dw = atom.dw;

	if (_mc->UseTDS == 0)
		return;

	// We need to copy the transpose of m_Mm to Mm.
	// we therefore cannot use the following command:
	for (ix = 0; ix < 3; ix++)
		for (iy = 0; iy < 3; iy++)
			Mm[ix][iy] = _Mm[iy][ix];

	Inverse_3x3(MmInv, Mm);

	/***************************************************************************
	 * Thermal Diffuse Scattering according to accurate phonon-dispersion
	 *
	 * Information in the phonon file will be stored in binary form as follows:
	 * Nk (number of k-points: 32-bit integer)
	 * Ns (number of atomic species 32-bit integer)
	 * M_1 M_2 ... M_Ns  (32-bit floats)
	 * kx(1) ky(1) kz(1) (32-bit floats)
	 * w_1(1) q_11 q_21 ... q_(3*Ns)1    (32-bit floats)
	 * w_2(1) q_12 q_22 ... q_(3*Ns)2
	 * :
	 * w_(3*Ns)(1) q_1Ns q_2Ns ... q_(3*Ns)Ns
	 * kx(2) ky(2) kz(2)
	 * :
	 * kx(Nk) ky(Nk) kz(Nk)
	 * :
	 *
	 *
	 * Note: only k-vectors in half of the Brillouin zone must be given, since
	 * w(k) = w(-k)
	 * also: 2D arrays will be read slowly varying index = first index (i*m+j)
	 **************************************************************************/

	if (wobScale == 0) {
		wobScale = 1.0 / (8 * M_PI * M_PI);
		sq3 = 1.0 / sqrt(3.0); /* sq3 is an additional needed factor which stems from
		 * int_-infty^infty exp(-x^2/2) x^2 dx = sqrt(pi)
		 * introduced in order to match the wobble factor with <u^2>
		 */
		scale = (float) sqrt(_sc->temperatureK / 300.0);
	}

	if (fpPhonon == NULL) {
		//		if ((fpPhonon = fopen(m_phononFile.string().c_str(), "r")) != NULL) {
		//			if (2 * sizeof(float) != sizeof(fftwf_complex)) {
		//				printf(
		//						"phononDisplacement: data type mismatch: fftw_complex != 2*float!\n");
		//				exit(0);
		//			}
		//			fread(&Nk, sizeof(int), 1, fpPhonon);
		//			fread(&Ns, sizeof(int), 1, fpPhonon);
		//			massPrim = (float_tt *) malloc(Ns * sizeof(float)); // masses for every atom in primitive basis
		//			fread(massPrim, sizeof(float), Ns, fpPhonon);
		//			kVecs.resize(boost::extents[Nk][3]);
		//			omega.resize(boost::extents[Nk][3*Ns]);
		//			 * omega is given in THz, but the 2pi-factor
		//			 * is still there, i.e. f=omega/2pi
		//			 */
		//			eigVecs = complex3D(Nk, 3 * Ns, 3 * Ns, "eigVecs"); // array of eigenvectors for every k-vector
		//			for (ix = 0; ix < Nk; ix++) {
		//				fread(kVecs[ix], sizeof(float), 3, fpPhonon);  // k-vector
		//				for (iy = 0; iy < 3 * Ns; iy++) {
		//					fread(omega[ix] + iy, sizeof(float), 1, fpPhonon);
		//					fread(eigVecs[ix][iy], 2 * sizeof(float), 3 * Ns, fpPhonon);
		//				}
		//			}
		//
		//			/* convert omega into q scaling factors, since we need those, instead of true omega:    */
		//			/* The 1/sqrt(2) term is from the dimensionality ((q1,q2) -> d=2)of the random numbers */
		//			for (ix = 0; ix < Nk; ix++) {
		//				for (idd = 0; idd < Ns; idd++)
		//					for (iy = 0; iy < 3; iy++) {
		//						// quantize the energy distribution:
		//						// tanh and exp give different results will therefore use exp
		//						// nkomega = (int)(1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/_sc->temperatureK));
		//						// wobble  =      (1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/_sc->temperatureK)-0.5);
		//						// nkomega = (int)(1.0/(exp(THZ_HBAR_KB*omega[ix][iy+3*id]/_sc->temperatureK)-1)+0.5);
		//						if (omega[ix][iy + 3 * idd] > 1e-4) {
		//							wobble =
		//									_sc->temperatureK > 0 ?
		//											(1.0
		//													/ (exp(
		//															THZ_HBAR_KB
		//															* omega[ix][iy
		//																		+ 3
		//																		* idd]
		//																		/ _sc->temperatureK)
		//															- 1)) :
		//															0;
		//							// if (ix == 0) printf("%g: %d %g\n",omega[ix][iy+3*id],nkomega,wobble);
		//							wobble = sqrt(
		//									(wobble + 0.5)
		//									/ (2 * M_PI * Nk * 2 * massPrim[idd]
		//																	* omega[ix][iy + 3 * idd]
		//																				* THZ_AMU_HBAR));
		//						} else
		//							wobble = 0;
		//						/* Ttotal += 0.25*massPrim[id]*((wobble*wobble)/(2*Ns))*
		//						 omega[ix][iy+3*id]*omega[ix][iy+3*id]*AMU_THZ2_A2_KB;
		//						 */
		//						omega[ix][iy + 3 * idd] = wobble;
		//					}  // idd
		//				// if (ix == 0) printf("\n");
		//			}
		//			// printf("Temperature: %g K\n",Ttotal);
		//			// printf("%d %d %d\n",(int)(0.4*(double)Nk/11.0),(int)(0.6*(double)Nk),Nk);
		//
		//		}
		//		fclose(fpPhonon);
	}  // end of if phononfile

	//
	// in the previous bracket: the phonon file is only read once.
	/////////////////////////////////////////////////////////////////////////////////////
	if (Nk > 800)
		printf("Will create phonon displacements for %d k-vectors - please wait ...\n", Nk);
	for (lambda = 0; lambda < 3 * Ns; lambda++)
		for (ik = 0; ik < Nk; ik++) {
			q1[lambda][ik] = (omega[ik][lambda] * gasdev());
			q2[lambda][ik] = (omega[ik][lambda] * gasdev());
		}
	// printf("Q: %g %g %g\n",q1[0][0],q1[5][8],q1[0][3]);
	// id seems to be the index of the correct atom, i.e. ranges from 0 .. Natom
	printf("created phonon displacements %d, %d, %d %d (eigVecs: %d %d %d)!\n", atom.Znum, Ns, Nk, id, Nk, 3 * Ns, 3 * Ns);
	/* loop over k and lambda:  */
	memset(u.data(), 0, 3 * sizeof(float_tt));
	for (lambda = 0; lambda < 3 * Ns; lambda++)
		for (ik = 0; ik < Nk; ik++) {
			// if (kVecs[ik][2] == 0){
			kR = 2 * M_PI * (icx * kVecs[ik][0] + icy * kVecs[ik][1] + icz * kVecs[ik][2]);
			//  kR = 2*M_PI*(blat[0][0]*kVecs[ik][0]+blat[0][1]*kVecs[ik][1]+blat[0][2]*kVecs[ik][2]);
			kRr = cos(kR);
			kRi = sin(kR);
			for (icoord = 0; icoord < 3; icoord++) {
				u[0][icoord] += q1[lambda][ik]
						* (eigVecs[ik][lambda][icoord + 3 * id].real() * kRr - eigVecs[ik][lambda][icoord + 3 * id].imag() * kRi)
						- q2[lambda][ik] * (eigVecs[ik][lambda][icoord + 3 * id].real() * kRi + eigVecs[ik][lambda][icoord + 3 * id].imag() * kRr);
			}
		}
	// printf("u: %g %g %g\n",u[0],u[1],u[2]);
	/* Convert the cartesian displacements back to reduced coordinates
	 */
	///////////////////////////////////////////////////////////////////////
	// Book keeping:
	u2[atom.Znum] += u[0][0] * u[0][0] + u[0][1] * u[0][1] + u[0][2] * u[0][2];
	ux += u[0][0];
	uy += u[0][1];
	uz += u[0][2];
	u[0][0] /= m_ax;
	u[0][1] /= m_by;
	u[0][2] /= m_cz;
	u2Count[atom.Znum]++;

	return;
}

int CrystalBuilder::AtomCompareZnum(const void *atPtr1, const void *atPtr2) {
	int comp = 0;
	atom *atom1, *atom2;

	atom1 = (atom *) atPtr1;
	atom2 = (atom *) atPtr2;
	/* Use the fact that z is the first element in the atom struct */
	comp = (atom1->Znum == atom2->Znum) ? 0 : ((atom1->Znum > atom2->Znum) ? -1 : 1);
	return comp;
}

int CrystalBuilder::AtomCompareZYX(const void *atPtr1, const void *atPtr2) {
	int comp = 0;
	atom *atom1, *atom2;

	atom1 = (atom *) atPtr1;
	atom2 = (atom *) atPtr2;
	/* Use the fact that z is the first element in the atom struct */
	comp = (atom1->r[2] == atom2->r[2]) ? 0 : ((atom1->r[2] > atom2->r[2]) ? 1 : -1);
	if (comp == 0) {
		comp = (atom1->r[1] == atom2->r[1]) ? 0 : ((atom1->r[1] > atom2->r[1]) ? 1 : -1);
		if (comp == 0) {
			comp = (atom1->r[0] == atom2->r[0]) ? 0 : ((atom1->r[0] > atom2->r[0]) ? 1 : -1);
		}
	}
	return comp;
}

void CrystalBuilder::WriteStructure(unsigned run_number) {
}

void CrystalBuilder::GetCrystalBoundaries(float_tt &min_x, float_tt &max_x, float_tt &min_y, float_tt &max_y, float_tt &min_z, float_tt &max_z) {
	if (_mc->periodicXY) {
		min_x = 0;
		max_x = _sc->nCellX * m_ax;
		min_y = 0;
		max_y = _sc->nCellX * m_by;
		min_z = _minZ;
		max_z = _maxZ;
	} else {
		min_x = _minX;
		max_x = _maxX;
		min_y = _minY;
		max_y = _maxY;
		min_z = _minZ;
		max_z = _maxZ;
	}
}
std::vector<int> CrystalBuilder::GetUniqueAtoms() {
	return _uniqueAtoms;
}
} // end namespace QSTEM

// *******************  Matrix manipulation ***********************
//    For now, use our own internal routines as has been done always.
//    For the future, consider using a linear algebra library instead - Eigen, BLAS/ATLAS/MKL/GOTO

#include "matrixlib.hpp"

namespace QSTEM {

void CrystalBuilder::Inverse_3x3(FloatArray2D& res, const FloatArray2D& a) {
// use function from matrixlib for now
	return inverse_3x3(res, a);
}

void CrystalBuilder::RotateVect(float_tt *vectIn, float_tt *vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z) {
	return rotateVect(vectIn, vectOut, phi_x, phi_y, phi_z);
}

void CrystalBuilder::MatrixProduct(const FloatArray2D& a, int Nxa, int Nya, const FloatArray2D& b, int Nxb, int Nyb, FloatArray2D& c) {
	return matrixProduct(a, Nxa, Nya, b, Nxb, Nyb, c);
}

void CrystalBuilder::RotateMatrix(const FloatArray2D& matrixIn, FloatArray2D& matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z) {
	return rotateMatrix(matrixIn, matrixOut, phi_x, phi_y, phi_z);
}

// ******************  end matrix manipulation

}// end namespace QSTEM
