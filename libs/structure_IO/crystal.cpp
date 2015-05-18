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


#include <stdlib.h>
#include "crystal.hpp"
#include "memory_fftw3.hpp"
#include "stemtypes_fftw3.hpp"
#include "random.hpp"
#include "structure_factories.hpp"
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

CrystalBuilder::CrystalBuilder() :
												m_minX(0),
												m_maxX(0),
												m_minY(0),
												m_maxY(0),
												m_minZ(0),
												m_maxZ(0),
												m_phononFile(boost::filesystem::path()),
												m_ax(0),
												_superCellBox(new superCellBox()),
												m_Mm(FloatArray2D(boost::extents[3][3])),
												m_MmInv(FloatArray2D(boost::extents[3][3]))
{}
CrystalBuilder::CrystalBuilder(boost::filesystem::path p): CrystalBuilder(){

}
CrystalBuilder::CrystalBuilder(StructureReaderFactory sfac,const ConfigPtr& c) : CrystalBuilder()  {
	_config = c;
	m_wobble_temp_scale = sqrt(_config->Structure.temperatureK / 300.0);
	boost::filesystem::path p(c->Structure.structureFilename);
//	BOOST_LOG_TRIVIAL(info) << p.extension().string();
	_structureReader = StructureReaderPtr(sfac[p.extension().string()](p));
	_structureWriter = CStructureWriterFactory::Get()->GetWriter(c->Structure.structureFilename,m_ax, m_by, m_cz);
}
//(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
double rotationCostFunc(const std::vector<double>& rotAngles,  std::vector<double> &grad, void* f_data){
	zoneAxisOptParams* p = reinterpret_cast<zoneAxisOptParams*>(f_data);
	std::vector<int> zone = p->zone;
	std::vector<int> refZone = p->refZone;
	arma::mat M = p->M;
//	std::vector<float_tt> zone,std::vector<float_tt> refZone
	float_tt	phi_x = rotAngles[0]*PI/180;
	float_tt	phi_y = rotAngles[1]*PI/180;
	float_tt	phi_z = rotAngles[2]*PI/180;
	float_tt cx = cos(phi_x),sx = sin(phi_x),
			cy = cos(phi_y),    sy = sin(phi_y),
			cz = cos(phi_z)  ,  sz = sin(phi_z);
	arma::mat Mx = {1,0,0,0,cx,-sx,0,sx,cx};
	arma::mat My = {cy,0,-sy,0,1,0,sy,0,cy};
	arma::mat Mz = {cz,sz,0,-sz,cz,0,0,0,1};
	BOOST_LOG_TRIVIAL(info)<< format("angles %f %f %f") % rotAngles[0]% rotAngles[1]% rotAngles[2];

	Mx.reshape(3,3);
	My.reshape(3,3);
	Mz.reshape(3,3);
//	std::cout << "Mx: "<< endl<< Mx << std::endl;
//	std::cout << "My: "<< endl<< My << std::endl;
//	std::cout << "Mz: "<< endl<< Mz << std::endl;
//	std::cout << "M: "<< endl<< M << std::endl;

	arma::mat Mrot = Mx*My*Mz;
//	std::cout << "Mrot: "<< endl<< Mrot << std::endl;
	arma::vec zone1 = {(double)zone[0],(double)zone[1],(double)zone[2]};
	arma::vec refZone1 = {(double)refZone[0],(double)refZone[1],(double)refZone[2]};
	arma::vec v = Mrot*M*zone1;

	v = v/sqrt(arma::sum(v%v));
//	std::cout << "refZone1: "<< endl<< refZone1 << std::endl;
//	std::cout << "v: "<< endl<< v << std::endl;
	auto tmp = v-refZone1;
//	std::cout << "tmp: "<< endl<< tmp << endl;
	auto tmp2 = tmp%tmp;
//	std::cout << "tmp2: "<< endl<< tmp2 << endl;
	double chi2 = arma::sum(tmp2);
//	std::cout << "chi2: "<< chi2 << endl;
	if(refZone[1]==0){
		arma::vec vtmp = {1,0,0};
		arma::vec vx = Mrot*M*vtmp;
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

void CrystalBuilder::CalculateTiltAnglesFromZoneAxis(){
	auto refZone = _config->Structure.zoneAxis;
	std::vector<int> refZone2(3),zone2(3) ;
	zone2[0]=0;
	zone2[1]=0;
	zone2[2]=1;

	int len = sqrt(refZone[0]*refZone[0]+refZone[1]*refZone[1]+refZone[2]*refZone[2]);
	for(int i=0;i<3;i++){
		refZone2[i] = refZone[i]/len;
	}

	arma::mat M = {m_Mm[0][0],m_Mm[1][0],m_Mm[2][0],
			m_Mm[0][1],m_Mm[1][1],m_Mm[2][1],
			m_Mm[0][2],m_Mm[1][2],m_Mm[2][2]};
	M.reshape(3,3);

	zoneAxisOptParams data(M,zone2,refZone2);

	opt o = opt(LN_NELDERMEAD, 3);
	o.set_lower_bounds(-180);
	o.set_upper_bounds(180);
	o.set_min_objective(&rotationCostFunc,(void*)&data);
	o.set_stopval(1e-15);
	double f_final;
	std::vector<double> x(3);
	x[0]=x[1]=x[2]=0;

	o.optimize(x,f_final);

	BOOST_LOG_TRIVIAL(info) << format("Rotating (%g,%g,%g) degrees from zone axis [0,0,1] to [%d,%d,%d]")
		% x[0]%x[1]%x[2]%refZone[0]%refZone[1]%refZone[2];

	_config->Structure.crystalTiltX = x[0];
	_config->Structure.crystalTiltY = x[1];
	_config->Structure.crystalTiltZ = x[2];


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
superCellBoxPtr CrystalBuilder::Build(){
	bool handleVacancies = true;
	int ncoord = _baseAtoms.size(), ncx, ncy, ncz, icx, icy, icz, jz;
	int i2, j, ix, iy, iz = 0;
	float_tt boxXmin = 0, boxXmax = 0, boxYmin = 0, boxYmax = 0, boxZmin = 0, boxZmax = 0;
	float_tt boxCenterX, boxCenterY, boxCenterZ, boxCenterXrot, boxCenterYrot, boxCenterZrot, bcX, bcY, bcZ;
	float_tt totOcc, choice, lastOcc;
	float_tt *u = NULL;
	FloatArray2D Mm = m_Mm;
	static int ncoord_old = 0;
	u = (float_tt *) malloc(3 * sizeof(float_tt));
	ncx = _config->Structure.nCellX;
	ncy = _config->Structure.nCellY;
	ncz = _config->Structure.nCellZ;

	ReadFromFile();

	unsigned natom = _baseAtoms.size() * ncx * ncy * ncz;
	_atoms.resize(natom);

	if (handleVacancies) {
		qsort(&_baseAtoms[0], ncoord, sizeof(atom),&CrystalBuilder::AtomCompareZYX);
	}
	/***********************************************************
	 * Read actual Data
	 ***********************************************************/
	for (int i = _baseAtoms.size() - 1; i >= 0; i--) {
		if (_config->Model.UseTDS) {
			m_u2[_baseAtoms[i].Znum] = 0;
		}
	}

	if(_config->Structure.rotateToZoneAxis){
		CalculateTiltAnglesFromZoneAxis();
	}

	if (_config->Structure.isBoxed) {
		/* at this point the atoms should have fractional coordinates */
		TiltBoxed(ncoord, handleVacancies);
		boxCenterX = _config->Structure.boxX;
		boxCenterY = _config->Structure.boxY;
		boxCenterZ = _config->Structure.boxZ;
	} else {  // work in NCell mode
		// atoms are in fractional coordinates so far, we need to convert them to
		// add the phonon displacement in this condition, because there we can
		// actually do the correct Eigenmode treatment.
		// but we will probably just do Einstein vibrations anyway:
		ReplicateUnitCell(handleVacancies);
		natom = ncoord * ncx * ncy * ncz;

		BOOST_LOG_TRIVIAL(trace) << "Atoms after replication of unit cells";

		for (int j = 0; j < natom; j++) {
			//			LOG(INFO) << format("atom %d: (%3.3f, %3.3f, %3.3f)\n") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
			// This converts also to cartesian coordinates
			float_tt x = Mm[0][0] * _atoms[j].x + Mm[1][0] * _atoms[j].y + Mm[2][0] * _atoms[j].z;
			float_tt y = Mm[0][1] * _atoms[j].x + Mm[1][1] * _atoms[j].y + Mm[2][1] * _atoms[j].z;
			float_tt z = Mm[0][2] * _atoms[j].x + Mm[1][2] * _atoms[j].y + Mm[2][2] * _atoms[j].z;

			_atoms[j].x = x;
			_atoms[j].y = y;
			_atoms[j].z = z;

			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % _atoms[j].x % _atoms[j].y % _atoms[j].z;
		}

		/***************************************************************
		 * Now let us tilt around the center of the full crystal
		 */

		bcX = ncx / 2.0;
		bcY = ncy / 2.0;
		bcZ = ncz / 2.0;
		u[0] = Mm[0][0] * bcX + Mm[1][0] * bcY + Mm[2][0] * bcZ;
		u[1] = Mm[0][1] * bcX + Mm[1][1] * bcY + Mm[2][1] * bcZ;
		u[2] = Mm[0][2] * bcX + Mm[1][2] * bcY + Mm[2][2] * bcZ;
		boxCenterX = u[0];
		boxCenterY = u[1];
		boxCenterZ = u[2];
		// Determine the size of the (rotated) super cell
		for (icx = 0; icx <= ncx; icx += ncx)
			for (icy = 0; icy <= ncy; icy += ncy)
				for (icz = 0; icz <= ncz; icz += ncz) {
					u[0] = Mm[0][0] * (icx - bcX) + Mm[1][0] * (icy - bcY)+ Mm[2][0] * (icz - bcZ);
					u[1] = Mm[0][1] * (icx - bcX) + Mm[1][1] * (icy - bcY)+ Mm[2][1] * (icz - bcZ);
					u[2] = Mm[0][2] * (icx - bcX) + Mm[1][2] * (icy - bcY)+ Mm[2][2] * (icz - bcZ);
					if ((_config->Structure.crystalTiltX != 0) || (_config->Structure.crystalTiltY != 0) || (_config->Structure.crystalTiltZ != 0))
						RotateVect(u, u, _config->Structure.crystalTiltX, _config->Structure.crystalTiltY, _config->Structure.crystalTiltZ); // simply applies rotation matrix
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

		//		printf("(%f, %f, %f): %f .. %f, %f .. %f, %f .. %f\n",m_ax,m_by,m_cz,boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax);
		if ((_config->Structure.crystalTiltX != 0)
				|| (_config->Structure.crystalTiltY != 0)
				|| (_config->Structure.crystalTiltZ != 0)) {
			for (int j =  0; j < (natom); j++) {
				u[0] = _atoms[j].x - boxCenterX;
				u[1] = _atoms[j].y - boxCenterY;
				u[2] = _atoms[j].z - boxCenterZ;
				RotateVect(u, u, _config->Structure.crystalTiltX, _config->Structure.crystalTiltY, _config->Structure.crystalTiltZ); // simply applies rotation matrix
				u[0] += boxCenterX;
				u[1] += boxCenterY;
				u[2] += boxCenterZ;
				_atoms[j].x = u[0];
				_atoms[j].y = u[1];
				_atoms[j].z = u[2];
			}
		} /* if tilts != 0 ... */

		// rebase to some atom at (0,0,0)
		BOOST_LOG_TRIVIAL(trace) << "Atoms after tilting and going to cartesian";
		for (int j = 0; j < natom; j++) {
			//			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)\n") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
			_atoms[j].x -= boxXmin;
			_atoms[j].y -= boxYmin;
			_atoms[j].z -= boxZmin;
			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % _atoms[j].x % _atoms[j].y % _atoms[j].z;
		}

		// Offset the atoms in x- and y-directions:
		// Do this after the rotation!
		if ((_config->Structure.xOffset != 0) || (_config->Structure.yOffset != 0)) {
			for (int j = 0 ; j < natom; j++) {
				_atoms[j].x += _config->Structure.xOffset;
				_atoms[j].y += _config->Structure.yOffset;
			}
		}
	} // end of Ncell mode conversion to cartesian coords and tilting.

	_sizeX = boxXmax - boxXmin;
	_sizeY = boxYmax - boxYmin;
	_sizeZ = boxZmax - boxZmin;
	CalculateCrystalBoundaries();
	free(u);
	_superCellBox->atoms = _atoms;
	_superCellBox->ax = _sizeX;
	_superCellBox->by = _sizeY;
	_superCellBox->cz = _sizeZ;
	_superCellBox->cmx = boxCenterX;
	_superCellBox->cmy = boxCenterY;
	_superCellBox->cmz = boxCenterZ;
	_superCellBox->uniqueatoms = _uniqueAtoms;
	return _superCellBox;
}


CrystalBuilder::~CrystalBuilder() {}
void CrystalBuilder::SetFillUnitCell(bool val){
	_fillUnitCell = val;
}
void CrystalBuilder::ReadFromFile(){
	_structureReader->ReadCellParams(m_Mm);
	_structureReader->ReadAtoms(_baseAtoms,_uniqueAtoms,true);
	CalculateCellDimensions();
	int s = _baseAtoms.size();
//	BOOST_LOG_TRIVIAL(info)<< format("Read %d atoms, tds: %d") % s % (_config->Model.UseTDS);
}
void CrystalBuilder::CalculateCrystalBoundaries() {
	int largestZindex = -1;
#pragma omp parallel for
	for (int i = 0; i < _atoms.size(); i++) {
#pragma omp critical
		{
			if (_atoms[i].x < m_minX)
				m_minX = _atoms[i].x;
		}
#pragma omp critical
		{
			if (_atoms[i].x > m_maxX)
				m_maxX = _atoms[i].x;
		}
#pragma omp critical
		{
			if (_atoms[i].y < m_minY)
				m_minY = _atoms[i].y;
		}
#pragma omp critical
		{
			if (_atoms[i].y > m_maxY)
				m_maxY = _atoms[i].y;
		}
#pragma omp critical
		{
			if (_atoms[i].z < m_minZ)
				m_minZ = _atoms[i].z;
		}
#pragma omp critical
		{
			if (_atoms[i].z > m_maxZ){
				m_maxZ = _atoms[i].z;
				largestZindex = i;
			}
		}
	}
	BOOST_LOG_TRIVIAL(trace)<<format("largest Z at index %d : (%f,%f,%f)")
																						%largestZindex%_atoms[largestZindex].x%_atoms[largestZindex].y%_atoms[largestZindex].z;
}

void CrystalBuilder::Init(unsigned run_number) {
	BOOST_LOG_TRIVIAL(trace)<<format("range of thermally displaced atoms (%d atoms): ") %_atoms.size();
	BOOST_LOG_TRIVIAL(trace)<<format("X: %g .. %g")% m_minX% m_maxX;
	BOOST_LOG_TRIVIAL(trace)<<format("Y: %g .. %g")% m_minY% m_maxY;
	BOOST_LOG_TRIVIAL(trace)<<format("Z: %g .. %g")% m_minZ% m_maxZ;

	qsort(&_atoms[0], _atoms.size(), sizeof(atom),&CrystalBuilder::AtomCompareZnum);
	WriteStructure(run_number);
}

void CrystalBuilder::DisplayParams() {
	BOOST_LOG_TRIVIAL(info) <<
			"***************************** Structure Parameters ***********************************************";
	BOOST_LOG_TRIVIAL(info) <<
			"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info)<<format("* Input file:           %s") %  _config->Structure.structureFilename.c_str();

	if (_config->Structure.isBoxed == false)
		BOOST_LOG_TRIVIAL(info)<<format(" Unit cell:            ax=%g by=%g cz=%g")% m_ax% m_by% m_cz;
	else {
		BOOST_LOG_TRIVIAL(info)<<format(" Size of Cube:         ax=%g by=%g cz=%g")
						% _config->Structure.boxX% _config->Structure.boxY%	_config->Structure.boxZ;
		BOOST_LOG_TRIVIAL(info)<<format(" Cube size adjusted:   %s")%( m_adjustCubeSize ? "yes" : "no");
	}
	BOOST_LOG_TRIVIAL(info)<<format(" Super cell:           %d x %d x %d unit cells")
													%_config->Structure.nCellX%		_config->Structure.nCellY% _config->Structure.nCellZ;
	BOOST_LOG_TRIVIAL(info)<<format(" Number of atoms:      %d (super cell)")% _atoms.size();
	BOOST_LOG_TRIVIAL(info)<<format(" Crystal tilt:         x=%g deg, y=%g deg, z=%g deg")
													%(_config->Structure.crystalTiltX)% (_config->Structure.crystalTiltY )%( _config->Structure.crystalTiltZ );
	BOOST_LOG_TRIVIAL(info)<<format(" Model dimensions:     ax=%gA, by=%gA, cz=%gA (after tilt)")
													%m_ax% m_by% m_cz;
	BOOST_LOG_TRIVIAL(info)<<format(" Temperature:          %gK") %_config->Structure.temperatureK;
	if (_config->Model.UseTDS)
		BOOST_LOG_TRIVIAL(info)<<format(" TDS:                  yes");
	else
		BOOST_LOG_TRIVIAL(info)<<format(" TDS:                  no");
}
void CrystalBuilder::SetSliceThickness(ModelConfig& mc){
	float_tt max_x, min_x, max_y, min_y, max_z, min_z, zTotal;
	auto box = Build(); // TODO find a way not to have to do this
	GetCrystalBoundaries(min_x, max_x, min_y, max_y, min_z, max_z);
	zTotal = max_z - min_z;

	switch (mc.SliceThicknessCalculation) {
	case SliceThicknessCalculation::Auto:
		mc.sliceThicknessAngstrom = (zTotal/((int)zTotal))+0.01*(zTotal/((int)zTotal));
		mc.nSlices = (int)zTotal;
		break;
	case SliceThicknessCalculation::NumberOfSlices:
		mc.sliceThicknessAngstrom = (zTotal/mc.nSlices)+0.01*(zTotal/mc.nSlices);
		break;
	case SliceThicknessCalculation::Thickness:
		mc.nSlices = (int)(zTotal / mc.sliceThicknessAngstrom);
		break;
	default:
		break;
	}
}
void CrystalBuilder::SetResolution(ModelConfig& mc, const PotentialConfig pc){
	float_tt max_x, min_x, max_y, min_y, max_z, min_z, zTotal;
	GetCrystalBoundaries(min_x, max_x, min_y, max_y, min_z, max_z);

	switch(mc.ResolutionCalculation) {
	case ResolutionCalculation::FILLN:
		mc.dx = (max_x - min_x)/mc.nx;
		mc.dy = (max_y - min_y)/mc.ny;
		_config->Structure.xOffset = 0;
		_config->Structure.yOffset = 0;
		break;
	case ResolutionCalculation::FILLRES:
		mc.nx = ceil((max_x - min_x) / mc.dx) ;
		mc.ny = ceil((max_y - min_y) / mc.dy) ;
		mc.dx = (max_x - min_x)/mc.nx;
		mc.dy = (max_y - min_y)/mc.ny;
		_config->Structure.xOffset = 0;
		_config->Structure.yOffset = 0;
		break;
	case ResolutionCalculation::BOXRES:
		mc.nx =  _config->Structure.boxX / mc.dx ;
		mc.ny =  _config->Structure.boxY / mc.dy ;
		if(mc.CenterSample) {
			_config->Structure.xOffset = _config->Structure.boxX/2 - (max_x-min_x)/2;
			_config->Structure.yOffset = _config->Structure.boxY/2 - (max_y-min_y)/2;
		} else {
			_config->Structure.xOffset = 0;
			_config->Structure.yOffset = 0;
		}
		break;
	case ResolutionCalculation::BOXN:
		mc.dx =  _config->Structure.boxX / mc.nx ;
		mc.dy =  _config->Structure.boxY / mc.ny ;
		if(mc.CenterSample) {
			_config->Structure.xOffset = _config->Structure.boxX/2 - (max_x-min_x)/2;
			_config->Structure.yOffset = _config->Structure.boxY/2 - (max_y-min_y)/2;
		} else {
			_config->Structure.xOffset = 0;
			_config->Structure.yOffset = 0;
		}
		break;

	}
}
void CrystalBuilder::SetCellParameters(float_tt ax, float_tt by, float_tt cz) {
	m_ax = ax;
	m_by = by;
	m_cz = cz;
}

void CrystalBuilder::SetNCells(unsigned nx, unsigned ny, unsigned nz) {
	_config->Structure.nCellX = nx;
	_config->Structure.nCellY = ny;
	_config->Structure.nCellZ = nz;
	// TODO: should this ever be false?
	bool handleVacancies = true;
	ReplicateUnitCell(handleVacancies);
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

	float_tt dx = m_ax / 2.0f - center.x;
	float_tt dy = m_by / 2.0f - center.y;
	float_tt dz = -center.z;
	for (size_t i = 0; i < _atoms.size(); i++) {
		_atoms[i].x += dx;
		if (_atoms[i].x < 0.0f)
			_atoms[i].x += m_ax;
		else if (_atoms[i].x > m_ax)
			_atoms[i].x -= m_ax;
		_atoms[i].y += dy;
		if (_atoms[i].y < 0.0f)
			_atoms[i].y += m_by;
		else if (_atoms[i].y > m_by)
			_atoms[i].y -= m_by;
		_atoms[i].z += dz;
		if (_atoms[i].z < 0.0f)
			_atoms[i].z += m_cz;
		else if (_atoms[i].z > m_cz)
			_atoms[i].z -= m_cz;
	}
}

// This function reads the atomic positions from fileName and also adds 
// Thermal displacements to their positions, if _config->Model.UseTDS is turned on.
void CrystalBuilder::MakeCrystal(bool handleVacancies) {
	int ncoord = _baseAtoms.size(), ncx, ncy, ncz, icx, icy, icz, jz;
	int i2, j, ix, iy, iz = 0;
	float_tt boxXmin = 0, boxXmax = 0, boxYmin = 0, boxYmax = 0, boxZmin = 0, boxZmax = 0;
	float_tt boxCenterX, boxCenterY, boxCenterZ, boxCenterXrot, boxCenterYrot, boxCenterZrot, bcX, bcY, bcZ;
	float_tt totOcc, choice, lastOcc;
	float_tt *u = NULL;
	FloatArray2D Mm = m_Mm;
	static int ncoord_old = 0;
	u = (float_tt *) malloc(3 * sizeof(float_tt));
	ncx = _config->Structure.nCellX;
	ncy = _config->Structure.nCellY;
	ncz = _config->Structure.nCellZ;

	unsigned natom = _baseAtoms.size() * ncx * ncy * ncz;
	_atoms.resize(natom);

	if (handleVacancies) {
		qsort(&_baseAtoms[0], ncoord, sizeof(atom),&CrystalBuilder::AtomCompareZYX);
	}
	/***********************************************************
	 * Read actual Data
	 ***********************************************************/
	for (int i = _baseAtoms.size() - 1; i >= 0; i--) {
		if (_config->Model.UseTDS) {
			m_u2[_baseAtoms[i].Znum] = 0;
		}
	}
	if (_config->Structure.isBoxed) {
		/* at this point the atoms should have fractional coordinates */
		TiltBoxed(ncoord, handleVacancies);
	} else {  // work in NCell mode
		// atoms are in fractional coordinates so far, we need to convert them to
		// add the phonon displacement in this condition, because there we can
		// actually do the correct Eigenmode treatment.
		// but we will probably just do Einstein vibrations anyway:
		ReplicateUnitCell(handleVacancies);
		natom = ncoord * ncx * ncy * ncz;

		BOOST_LOG_TRIVIAL(trace) << "Atoms after replication of unit cells";

		for (int j = 0; j < natom; j++) {
			//			LOG(INFO) << format("atom %d: (%3.3f, %3.3f, %3.3f)\n") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
			// This converts also to cartesian coordinates
			float_tt x = Mm[0][0] * _atoms[j].x + Mm[1][0] * _atoms[j].y + Mm[2][0] * _atoms[j].z;
			float_tt y = Mm[0][1] * _atoms[j].x + Mm[1][1] * _atoms[j].y + Mm[2][1] * _atoms[j].z;
			float_tt z = Mm[0][2] * _atoms[j].x + Mm[1][2] * _atoms[j].y + Mm[2][2] * _atoms[j].z;

			_atoms[j].x = x;
			_atoms[j].y = y;
			_atoms[j].z = z;

			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % _atoms[j].x % _atoms[j].y % _atoms[j].z;
		}

		/***************************************************************
		 * Now let us tilt around the center of the full crystal
		 */

		bcX = ncx / 2.0;
		bcY = ncy / 2.0;
		bcZ = ncz / 2.0;
		u[0] = Mm[0][0] * bcX + Mm[1][0] * bcY + Mm[2][0] * bcZ;
		u[1] = Mm[0][1] * bcX + Mm[1][1] * bcY + Mm[2][1] * bcZ;
		u[2] = Mm[0][2] * bcX + Mm[1][2] * bcY + Mm[2][2] * bcZ;
		boxCenterX = u[0];
		boxCenterY = u[1];
		boxCenterZ = u[2];
		// Determine the size of the (rotated) super cell
		for (icx = 0; icx <= ncx; icx += ncx)
			for (icy = 0; icy <= ncy; icy += ncy)
				for (icz = 0; icz <= ncz; icz += ncz) {
					u[0] = Mm[0][0] * (icx - bcX) + Mm[1][0] * (icy - bcY)+ Mm[2][0] * (icz - bcZ);
					u[1] = Mm[0][1] * (icx - bcX) + Mm[1][1] * (icy - bcY)+ Mm[2][1] * (icz - bcZ);
					u[2] = Mm[0][2] * (icx - bcX) + Mm[1][2] * (icy - bcY)+ Mm[2][2] * (icz - bcZ);
					if ((_config->Structure.crystalTiltX != 0) || (_config->Structure.crystalTiltY != 0) || (_config->Structure.crystalTiltZ != 0))
						RotateVect(u, u, _config->Structure.crystalTiltX, _config->Structure.crystalTiltY, _config->Structure.crystalTiltZ); // simply applies rotation matrix
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

		//		printf("(%f, %f, %f): %f .. %f, %f .. %f, %f .. %f\n",m_ax,m_by,m_cz,boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax);
		if ((_config->Structure.crystalTiltX != 0) || (_config->Structure.crystalTiltY != 0) || (_config->Structure.crystalTiltZ != 0)) {
			for (int j =  0; j < (natom); j++) {
				u[0] = _atoms[j].x - boxCenterX;
				u[1] = _atoms[j].y - boxCenterY;
				u[2] = _atoms[j].z - boxCenterZ;
				RotateVect(u, u, _config->Structure.crystalTiltX, _config->Structure.crystalTiltY, _config->Structure.crystalTiltZ); // simply applies rotation matrix
				u[0] += boxCenterX;
				u[1] += boxCenterY;
				u[2] += boxCenterZ;
				_atoms[j].x = u[0];
				_atoms[j].y = u[1];
				_atoms[j].z = u[2];
			}
		} /* if tilts != 0 ... */

		// rebase to some atom at (0,0,0)
		BOOST_LOG_TRIVIAL(trace) << "Atoms after tilting and going to cartesian";
		for (int j = 0; j < natom; j++) {
			//			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)\n") % j % m_atoms[j].x % m_atoms[j].y % m_atoms[j].z;
			_atoms[j].x -= boxXmin;
			_atoms[j].y -= boxYmin;
			_atoms[j].z -= boxZmin;
			BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f)") % j % _atoms[j].x % _atoms[j].y % _atoms[j].z;
		}

		// Offset the atoms in x- and y-directions:
		// Do this after the rotation!
		if ((_config->Structure.xOffset != 0) || (_config->Structure.yOffset != 0)) {
			for (int j = 0 ; j < natom; j++) {
				_atoms[j].x += _config->Structure.xOffset;
				_atoms[j].y += _config->Structure.yOffset;
			}
		}
	} // end of Ncell mode conversion to cartesian coords and tilting.
	// end of loop over atoms

	_sizeX = boxXmax - boxXmin;
	_sizeY = boxYmax - boxYmin;
	_sizeZ = boxZmax - boxZmin;
	CalculateCrystalBoundaries();
	free(u);
}

// Uses m_Mm to calculate ax, by, cz, and alpha, beta, gamma
void CrystalBuilder::CalculateCellDimensions() {
	m_ax = sqrt(m_Mm[0][0] * m_Mm[0][0] + m_Mm[0][1] * m_Mm[0][1]+ m_Mm[0][2] * m_Mm[0][2]);
	m_by = sqrt(m_Mm[1][0] * m_Mm[1][0] + m_Mm[1][1] * m_Mm[1][1]+ m_Mm[1][2] * m_Mm[1][2]);
	m_cz = sqrt(m_Mm[2][0] * m_Mm[2][0] + m_Mm[2][1] * m_Mm[2][1]+ m_Mm[2][2] * m_Mm[2][2]);
	m_cGamma = atan2(m_Mm[1][1], m_Mm[1][0]);
	m_cBeta = acos(m_Mm[2][0] / m_cz);
	m_cAlpha = acos(m_Mm[2][1] * sin(m_cGamma) / m_cz + cos(m_cBeta) * cos(m_cGamma));
	m_cGamma /= (float) PI180;
	m_cBeta /= (float) PI180;
	m_cAlpha /= (float) PI180;
}

// TODO: old version return a pointer to new atom positions.  Can we do this in place?
//		If not, just set atoms to the new vector.
void CrystalBuilder::TiltBoxed(int ncoord, bool handleVacancies) {
	int atomKinds = 0;
	int iatom, jVac, jequal, jChoice, i2, ix, iy, iz, atomCount = 0, atomSize;
	static float_tt *axCell, *byCell, *czCell = NULL;
	//	static float_tt **Mm = NULL, **Mminv = NULL, **MpRed = NULL, **MpRedInv =
	//			NULL;
	//	static float_tt **MbPrim = NULL, **MbPrimInv = NULL, **MmOrig = NULL,
	//			**MmOrigInv = NULL;
	//static float_tt **a = NULL,**aOrig = NULL,**b= NULL,**bfloor=NULL,**blat=NULL;
	//	std::vector<float_tt> a(3, 0), aOrig(3, 0), b(3, 0), bfloor(3, 0), blat(3, 0);
	//static float_tt *uf;
	FloatArray2D a(boost::extents[3][1]),
			aOrig(boost::extents[3][1]),
			b(boost::extents[3][1]) ;
	static int oldAtomSize = 0;
	double x, y, z, dx, dy, dz;
	double totOcc, lastOcc, choice;
	atom newAtom;
	int Ncells;
	FloatArray2D u(boost::extents[1][3]),uf(boost::extents[1][3]) ;

	FloatArray2D Mm(boost::extents[3][3]),
			Mminv(boost::extents[3][3]),
			MpRed(boost::extents[3][3]),// conversion lattice to obtain red. prim. coords/ from reduced cubic rect.
			MpRedInv(boost::extents[3][3]),//conversion lattice to obtain red. cub. coords from reduced primitive lattice coords
			MbPrim(boost::extents[3][3]),// double version of primitive lattice basis
			MbPrimInv(boost::extents[3][3]),// double version of inverse primitive lattice basis
			MmOrig(boost::extents[3][3]),
			MmOrigInv(boost::extents[3][3]);

	//unsigned jz;

	Ncells = _config->Structure.nCellX * _config->Structure.nCellY * _config->Structure.nCellZ;


	dx = 0;
	dy = 0;
	dz = 0;
	dx = m_offsetX;
	dy = m_offsetY;
	// We need to copy the transpose of m_Mm to Mm.
	// we therefore cannot use the following command:
	// memcpy(Mm[0],m_Mm[0],3*3*sizeof(double));
	for (ix = 0; ix < 3; ix++)
		for (iy = 0; iy < 3; iy++)
			Mm[ix][iy] = m_Mm[iy][ix];

	memcpy(MmOrig.data(), Mm.data(), 3 * 3 * sizeof(double));
	Inverse_3x3(MmOrigInv, MmOrig);
	/* remember that the angles are in rad: */
	RotateMatrix(Mm, Mm, _config->Structure.crystalTiltX/PI180, _config->Structure.crystalTiltY/PI180, _config->Structure.crystalTiltZ/PI180);
	Inverse_3x3(Mminv, Mm);  // computes Mminv from Mm!
	/* find out how far we will have to go in unit of unit cell vectors.  when creating the supercell by
	 * checking the number of unit cell vectors necessary to reach every corner of the supercell box.
	 */
	// matrixProduct(a,1,3,Mminv,3,3,b);
	MatrixProduct(Mminv, 3, 3, a, 3, 1, b);
	// showMatrix(Mm,3,3,"M");
	// showMatrix(Mminv,3,3,"M");
	unsigned nxmax, nymax, nzmax;
	unsigned nxmin = nxmax = (int) floor(b[0][0] - dx);
	unsigned nymin = nymax = (int) floor(b[1][0] - dy);
	unsigned nzmin = nzmax = (int) floor(b[2][0] - dz);
	for (ix = 0; ix <= 1; ix++)
		for (iy = 0; iy <= 1; iy++)
			for (iz = 0; iz <= 1; iz++) {
				a[0][0] = ix * _config->Structure.boxX - dx;
				a[1][0] = iy * _config->Structure.boxY - dy;
				a[2][0] = iz * _config->Structure.boxZ - dz;

				// matrixProduct(a,1,3,Mminv,3,3,b);
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

	atomSize =	(1 + (nxmax - nxmin) * (nymax - nymin) * (nzmax - nzmin) * ncoord);
	if (atomSize != oldAtomSize) {
		_atoms.resize(atomSize);
		oldAtomSize = atomSize;
	}
	// showMatrix(Mm,3,3,"Mm");
	// showMatrix(Mminv,3,3,"Mminv");
	// printf("Range: (%d..%d, %d..%d, %d..%d)\n",
	// nxmin,nxmax,nymin,nymax,nzmin,nzmax);

	atomCount = 0;
	jVac = 0;
	for (iatom = 0; iatom < ncoord;) {
		// printf("%d: (%g %g %g) %d\n",iatom,unitAtoms[iatom].x,unitAtoms[iatom].y,
		//   unitAtoms[iatom].z,unitAtoms[iatom].Znum);
		memcpy(&newAtom, &_baseAtoms[iatom], sizeof(atom));
		//for (jz=0;jz<m_Znums.size();jz++)	if (m_Znums[jz] == newAtom.Znum) break;
		// allocate more memory, if there is a new element
		/*
		 if (jz == atomKinds) {
		 atomKinds++;
		 if (atomKinds > m_atomKinds) {
		 m_Znums = (int *)realloc(m_Znums,atomKinds*sizeof(int));
		 m_atomKinds = atomKinds;
		 // printf("%d kinds (%d)\n",atomKinds,atoms[i].Znum);
		 }
		 m_Znums[jz] = newAtom.Znum;
		 }
		 */
		/////////////////////////////////////////////////////
		// look for atoms at equal position
		if ((handleVacancies) && (newAtom.Znum > 0)) {
			totOcc = newAtom.occ;
			for (jequal = iatom + 1; jequal < ncoord; jequal++) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase
				// the total occupany and the counter jequal.
				if ((fabs(newAtom.x - _baseAtoms[jequal].x) < 1e-6)
						&& (fabs(newAtom.y - _baseAtoms[jequal].y) < 1e-6)
						&& (fabs(newAtom.z - _baseAtoms[jequal].z) < 1e-6)) {
					totOcc += _baseAtoms[jequal].occ;
				} else
					break;
			} // jequal-loop
		} else {
			jequal = iatom + 1;
			totOcc = 1;
		}

		unsigned atomCount = 0;

		// printf("%d: %d\n",atomCount,jz);
		for (ix = nxmin; ix <= nxmax; ix++) {
			for (iy = nymin; iy <= nymax; iy++) {
				for (iz = nzmin; iz <= nzmax; iz++) {
					// atom position in cubic reduced coordinates:
					aOrig[0][0] = ix + newAtom.x;
					aOrig[1][0] = iy + newAtom.y;
					aOrig[2][0] = iz + newAtom.z;

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					// All we need to decide is whether to include the atom at all (if totOcc < 1
					// of which of the atoms at equal positions to include
					jChoice = iatom;  // This will be the atom we wil use.
					if ((totOcc < 1) || (jequal > iatom + 1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values.
						//
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totOcc < 1.0)
							choice = ran1();
						else
							choice = totOcc * ran1();
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2 = iatom; i2 < jequal; i2++) {
							// atoms[atomCount].Znum = unitAtoms[i2].Znum;
							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc)
							if ((choice < lastOcc)
									|| (choice >= lastOcc + _baseAtoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								// atoms[atomCount].Znum =  0;  // vacancy
								jVac++;
							} else {
								jChoice = i2;
							}
							lastOcc += _baseAtoms[i2].occ;
						}
					}
					// here we select the index of our chosen atom
					//if (jChoice != iatom) {
					//std::vector<unsigned>::iterator item = std::find(m_Znums.begin(), m_Znums.end(), m_baseAtoms[jChoice].Znum);
					//if (item!=m_Znums.end())
					//jz=item-m_Znums.begin();
					//for (jz=0;jz<m_Znums.size();jz++)	if (m_Znums[jz] == m_baseAtoms[jChoice].Znum) break;
					//}

					// here we need to call phononDisplacement:
					// phononDisplacement(u,muls,iatom,ix,iy,iz,atomCount,atoms[i].dw,*natom,atoms[i].Znum);
					if (_config->Model.displacementType == DisplacementType::Einstein) {
						if (_config->Model.UseTDS) {
							EinsteinDisplacement(u,newAtom);
							a[0][0] = aOrig[0][0] + u[0][0];
							a[1][0] = aOrig[1][0] + u[0][1];
							a[2][0] = aOrig[2][0] + u[0][2];
						} else {
							a[0][0] = aOrig[0][0];
							a[1][0] = aOrig[1][0];
							a[2][0] = aOrig[2][0];
						}
					} else {
						BOOST_LOG_TRIVIAL(fatal)<<format("Cannot handle phonon-distribution mode for boxed sample yet. Exiting...");
						exit(0);
					}
					// matrixProduct(aOrig,1,3,Mm,3,3,b);
					MatrixProduct(Mm, 3, 3, aOrig, 3, 1,b);

					// if (atomCount < 2) {showMatrix(a,1,3,"a");showMatrix(b,1,3,"b");}
					// b now contains atom positions in cartesian coordinates */
					x = b[0][0] + dx;
					y = b[1][0] + dy;
					z = b[2][0] + dz;

					// include atoms that are within the box
					if ((x >= 0) &&
							(x <= _config->Structure.boxX) && (y >= 0) &&
							(y <= _config->Structure.boxY) && (z >= 0) &&
							(z <= _config->Structure.boxZ)) {
						// matrixProduct(a,1,3,Mm,3,3,b);
						MatrixProduct(Mm, 3, 3, a, 3, 1, b);
						_atoms[atomCount].x = b[0][0] + dx;
						_atoms[atomCount].y = b[1][0] + dy;
						_atoms[atomCount].z = b[2][0] + dz;
						_atoms[atomCount].dw = _baseAtoms[jChoice].dw;
						_atoms[atomCount].occ = _baseAtoms[jChoice].occ;
						_atoms[atomCount].q = _baseAtoms[jChoice].q;
						_atoms[atomCount].Znum = _baseAtoms[jChoice].Znum;


						atomCount++;
						/*
						 if (m_baseAtoms[jChoice].Znum > 22)
						 printf("Atomcount: %d, Z = %d\n",atomCount,m_baseAtoms[jChoice].Znum);
						 */
					}
				} /* iz ... */
			} /* iy ... */
		} /* ix ... */
		iatom = jequal;
	} /* iatom ... */
	BOOST_LOG_TRIVIAL(trace)<<format("Removed %d atoms because of multiple occupancy or occupancy < 1") % jVac;
	m_ax = _config->Structure.boxX;
	m_by = _config->Structure.boxY;
	m_cz = _config->Structure.boxZ;
	// call phononDisplacement again to update displacement data:
	// Not supported yet
	//	PhononDisplacement(u, iatom, ix, iy, iz, newAtom, false);
}  // end of 'tiltBoxed(...)'

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
	double totalOccupancy;
	double choice, lastOcc;

	ncx = _config->Structure.nCellX;
	ncy = _config->Structure.nCellY;
	ncz = _config->Structure.nCellZ;
	FloatArray2D u(boost::extents[1][3]);

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
				if ((fabs(_atoms[i].x - _atoms[jequal].x) < 1e-6)
						&& (fabs(_atoms[i].y - _atoms[jequal].y) < 1e-6)
						&& (fabs(_atoms[i].z - _atoms[jequal].z) < 1e-6)) {
					totalOccupancy += _atoms[jequal].occ;

				} else break;

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
					if ((totalOccupancy < 1) || (jequal < i-1)) { // found atoms at equal positions or an occupancy less than 1!
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

					if (_config->Model.UseTDS)
						switch(_config->Model.displacementType){
						case DisplacementType::Einstein:
							EinsteinDisplacement(u,_baseAtoms[i]);
							break;
						case DisplacementType::Phonon:
							PhononDisplacement(u, jChoice, icx, icy, icz, _atoms[jChoice], false);
							break;
						case DisplacementType::None:
						default:
							break;
						}
					for (i2 = i; i2 > jequal; i2--) {
						_atoms[jCell + i2].x = _baseAtoms[i2].x + icx + u[0][0];
						_atoms[jCell + i2].y = _baseAtoms[i2].y + icy + u[0][1];
						_atoms[jCell + i2].z = _baseAtoms[i2].z + icz + u[0][2];
						BOOST_LOG_TRIVIAL(trace) << format("atom %d: (%3.3f, %3.3f, %3.3f) Z=%d")
												% (jCell + i2) % _atoms[jCell + i2].x % _atoms[jCell + i2].y
												% _atoms[jCell + i2].z %  _atoms[jCell + i2].Znum ;
					}
				}  // for (icz=ncz-1;icz>=0;icz--)
			} // for (icy=ncy-1;icy>=0;icy--)
		} // for (icx=ncx-1;icx>=0;icx--)
		i = jequal;
	} // for (i=ncoord-1;i>=0;)
	if ((jVac > 0) )
		BOOST_LOG_TRIVIAL(trace)<<format("Removed %d atoms because of occupancies < 1 or multiple atoms in the same place") %jVac;
}

void CrystalBuilder::DisplaceAtoms() {
	std::vector<atom>::iterator at = _atoms.begin(), end = _atoms.end();
	FloatArray2D u(boost::extents[1][3]);

	for (at; at != end; ++at) {
		EinsteinDisplacement(u, (*at));
		// Add obtained u to atomic position
		(*at).x += u[0][0];
		(*at).y += u[0][1];
		(*at).z += u[0][2];
	}
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
	MatrixProduct(u, 1, 3, m_MmInv, 3, 3,uf);
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
void CrystalBuilder::PhononDisplacement(FloatArray2D &u, int id, int icx,
		int icy, int icz, atom &atom, bool printReport) {
	int ix, iy, idd; // iz;
	static FILE *fpPhonon = NULL;
	static int Nk, Ns;  // number of k-vectors and atoms per primitive unit cell
	static float_tt *massPrim;   // masses for every atom in primitive basis
	static float_tt **omega;     // array of eigenvalues for every k-vector
	static complex_tt ***eigVecs;  // array of eigenvectors for every k-vector
	static float_tt **kVecs;     // array for Nk 3-dim k-vectors
	static float_tt **q1 = NULL, **q2 = NULL;
	int ik, lambda, icoord; // Ncells, nkomega;
	double kR, kRi, kRr, wobble;
	static float_tt *u2T, ux = 0, uy = 0, uz = 0; // u2Collect=0; // Ttotal=0;
	std::map<unsigned, float_tt> u2;
	std::map<unsigned, unsigned> u2Count;
	// static double uxCollect=0,uyCollect=0,uzCollect=0;
	static int *u2CountT, runCount = 1, u2Size = -1;
	static long iseed = 0;
	FloatArray2D Mm(boost::extents[3][3]),
			MmInv(boost::extents[3][3]);
	// static float_tt **MmOrig=NULL,**MmOrigInv=NULL;
	static float_tt *axCell, *byCell, *czCell;
	static float_tt wobScale = 0, sq3, scale = 0;
	std::vector<float_tt> uf(3), b(3);

	float_tt dw = atom.dw;

	if (_config->Model.UseTDS == 0)
		return;

	// We need to copy the transpose of m_Mm to Mm.
	// we therefore cannot use the following command:
	for (ix = 0; ix < 3; ix++)
		for (iy = 0; iy < 3; iy++)
			Mm[ix][iy] = m_Mm[iy][ix];

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
		scale = (float) sqrt(_config->Structure.temperatureK / 300.0);
	}

	if (fpPhonon == NULL) {
		if ((fpPhonon = fopen(m_phononFile.string().c_str(), "r")) != NULL) {
			if (2 * sizeof(float) != sizeof(fftwf_complex)) {
				printf(
						"phononDisplacement: data type mismatch: fftw_complex != 2*float!\n");
				exit(0);
			}
			fread(&Nk, sizeof(int), 1, fpPhonon);
			fread(&Ns, sizeof(int), 1, fpPhonon);
			massPrim = (float_tt *) malloc(Ns * sizeof(float)); // masses for every atom in primitive basis
			fread(massPrim, sizeof(float), Ns, fpPhonon);
			kVecs = float2D(Nk, 3, "kVecs");
			omega = float2D(Nk, 3 * Ns, "omega"); /* array of eigenvalues for every k-vector
			 * omega is given in THz, but the 2pi-factor
			 * is still there, i.e. f=omega/2pi
			 */
			eigVecs = complex3D(Nk, 3 * Ns, 3 * Ns, "eigVecs"); // array of eigenvectors for every k-vector
			for (ix = 0; ix < Nk; ix++) {
				fread(kVecs[ix], sizeof(float), 3, fpPhonon);  // k-vector
				for (iy = 0; iy < 3 * Ns; iy++) {
					fread(omega[ix] + iy, sizeof(float), 1, fpPhonon);
					fread(eigVecs[ix][iy], 2 * sizeof(float), 3 * Ns, fpPhonon);
				}
			}

			/* convert omega into q scaling factors, since we need those, instead of true omega:    */
			/* The 1/sqrt(2) term is from the dimensionality ((q1,q2) -> d=2)of the random numbers */
			for (ix = 0; ix < Nk; ix++) {
				for (idd = 0; idd < Ns; idd++)
					for (iy = 0; iy < 3; iy++) {
						// quantize the energy distribution:
						// tanh and exp give different results will therefore use exp
						// nkomega = (int)(1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/_config->Structure.temperatureK));
						// wobble  =      (1.0/tanh(THZ_HBAR_2KB*omega[ix][iy+3*id]/_config->Structure.temperatureK)-0.5);
						// nkomega = (int)(1.0/(exp(THZ_HBAR_KB*omega[ix][iy+3*id]/_config->Structure.temperatureK)-1)+0.5);
						if (omega[ix][iy + 3 * idd] > 1e-4) {
							wobble =
									_config->Structure.temperatureK > 0 ?
											(1.0
													/ (exp(
															THZ_HBAR_KB
															* omega[ix][iy
																		+ 3
																		* idd]
																		/ _config->Structure.temperatureK)
															- 1)) :
															0;
							// if (ix == 0) printf("%g: %d %g\n",omega[ix][iy+3*id],nkomega,wobble);
							wobble = sqrt(
									(wobble + 0.5)
									/ (2 * M_PI * Nk * 2 * massPrim[idd]
																	* omega[ix][iy + 3 * idd]
																				* THZ_AMU_HBAR));
						} else
							wobble = 0;
						/* Ttotal += 0.25*massPrim[id]*((wobble*wobble)/(2*Ns))*
						 omega[ix][iy+3*id]*omega[ix][iy+3*id]*AMU_THZ2_A2_KB;
						 */
						omega[ix][iy + 3 * idd] = wobble;
					}  // idd
				// if (ix == 0) printf("\n");
			}
			// printf("Temperature: %g K\n",Ttotal);
			// printf("%d %d %d\n",(int)(0.4*(double)Nk/11.0),(int)(0.6*(double)Nk),Nk);
			q1 = float2D(3 * Ns, Nk, "q1");
			q2 = float2D(3 * Ns, Nk, "q2");

		}
		fclose(fpPhonon);
	}  // end of if phononfile

	//
	// in the previous bracket: the phonon file is only read once.
	/////////////////////////////////////////////////////////////////////////////////////
	if (Nk > 800)
		printf("Will create phonon displacements for %d k-vectors - please wait ...\n",Nk);
	for (lambda = 0; lambda < 3 * Ns; lambda++)
		for (ik = 0; ik < Nk; ik++) {
			q1[lambda][ik] = (omega[ik][lambda] * gasdev());
			q2[lambda][ik] = (omega[ik][lambda] * gasdev());
		}
	// printf("Q: %g %g %g\n",q1[0][0],q1[5][8],q1[0][3]);
	// id seems to be the index of the correct atom, i.e. ranges from 0 .. Natom
	printf("created phonon displacements %d, %d, %d %d (eigVecs: %d %d %d)!\n",
			atom.Znum, Ns, Nk, id, Nk, 3 * Ns, 3 * Ns);
	/* loop over k and lambda:  */
	memset(u.data(), 0, 3 * sizeof(double));
	for (lambda = 0; lambda < 3 * Ns; lambda++)
		for (ik = 0; ik < Nk; ik++) {
			// if (kVecs[ik][2] == 0){
			kR = 2 * M_PI
					* (icx * kVecs[ik][0] + icy * kVecs[ik][1]
															+ icz * kVecs[ik][2]);
			//  kR = 2*M_PI*(blat[0][0]*kVecs[ik][0]+blat[0][1]*kVecs[ik][1]+blat[0][2]*kVecs[ik][2]);
			kRr = cos(kR);
			kRi = sin(kR);
			for (icoord = 0; icoord < 3; icoord++) {
				u[0][icoord] +=
						q1[lambda][ik]
								   * (eigVecs[ik][lambda][icoord + 3 * id].real() * kRr
										   - eigVecs[ik][lambda][icoord + 3 * id].imag()
										   * kRi)
										   - q2[lambda][ik]
														* (eigVecs[ik][lambda][icoord + 3 * id].real()
																* kRi
																+ eigVecs[ik][lambda][icoord
																					  + 3 * id].imag() * kRr);
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
	comp = (atom1->Znum == atom2->Znum) ?
			0 : ((atom1->Znum > atom2->Znum) ? -1 : 1);
	return comp;
}

int CrystalBuilder::AtomCompareZYX(const void *atPtr1, const void *atPtr2) {
	int comp = 0;
	atom *atom1, *atom2;

	atom1 = (atom *) atPtr1;
	atom2 = (atom *) atPtr2;
	/* Use the fact that z is the first element in the atom struct */
	comp = (atom1->z == atom2->z) ? 0 : ((atom1->z > atom2->z) ? 1 : -1);
	if (comp == 0) {
		comp = (atom1->y == atom2->y) ? 0 : ((atom1->y > atom2->y) ? 1 : -1);
		if (comp == 0) {
			comp = (atom1->x == atom2->x) ?
					0 : ((atom1->x > atom2->x) ? 1 : -1);
		}
	}
	return comp;
}

void CrystalBuilder::WriteStructure(unsigned run_number) {
	// to_string is C++0x - may not work on older compilers
	_structureWriter->Write(_atoms, std::to_string(run_number));
	/*
	 if (m_cfgFile != NULL) {
	 sprintf(buf,"%s/%s",m_folder,m_cfgFile);
	 // append the TDS run number
	 if (strcmp(buf+strlen(buf)-4,".cfg") == 0) *(buf+strlen(buf)-4) = '\0';
	 if (_config->Model.UseTDS) sprintf(buf+strlen(buf),"_%d.cfg",m_avgCount);
	 else sprintf(buf+strlen(buf),".cfg");

	 // printf("Will write CFG file <%s> (%d)\n",buf,_config->Model.UseTDS)
	 writeCFG(atoms,natom,buf,muls);
	 }
	 */
}

void CrystalBuilder::GetCrystalBoundaries(float_tt &min_x, float_tt &max_x, float_tt &min_y, float_tt &max_y, float_tt &min_z, float_tt &max_z) {
	if(_config->Potential.periodicXY){
		min_x = 0;
		max_x = _config->Structure.nCellX * m_ax;
		min_y = 0;
		max_y = _config->Structure.nCellX * m_by;
		min_z = m_minZ;
		max_z = m_maxZ;
	} else {
		min_x = m_minX;
		max_x = m_maxX;
		min_y = m_minY;
		max_y = m_maxY;
		min_z = m_minZ;
		max_z = m_maxZ;
	}
}
std::vector<int> CrystalBuilder::GetUniqueAtoms(){
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

void CrystalBuilder::RotateVect(float_tt *vectIn, float_tt *vectOut, float_tt phi_x,
		float_tt phi_y, float_tt phi_z) {
	return rotateVect(vectIn, vectOut, phi_x, phi_y, phi_z);
}

void CrystalBuilder::MatrixProduct(const FloatArray2D& a,int Nxa, int Nya, const FloatArray2D& b,int Nxb, int Nyb, FloatArray2D& c) {
	return matrixProduct(a, Nxa, Nya, b, Nxb, Nyb, c);
}

void CrystalBuilder::RotateMatrix(const FloatArray2D& matrixIn, FloatArray2D& matrixOut,
		float_tt phi_x, float_tt phi_y, float_tt phi_z) {
	return rotateMatrix(matrixIn, matrixOut, phi_x, phi_y, phi_z);
}

// ******************  end matrix manipulation

}// end namespace QSTEM
