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


#include "pot_3d.hpp"

namespace QSTEM
{

C3DPotential::C3DPotential(const ConfigPtr c, PersistenceManagerPtr p ) : RealSpacePotential(c,p )
{
	m_sliceStep = 2*_config->Model.nx*_config->Model.ny;
	m_boxNz = (int)(_config->Potential.AtomRadiusAngstrom/m_ddz+2.0);
	m_c += 2*_config->Potential.AtomRadiusAngstrom;
}
void C3DPotential::ComputeAtomPotential(int Znum){

}
void C3DPotential::SliceSetup(){
	CPotential::SliceSetup();
}
void C3DPotential::DisplayParams()
{
	CPotential::DisplayParams();
	printf("* Potential calculation: 3D (real space method)");
}

void C3DPotential::AtomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B)
{
	float_tt dx, dy, dz;
	int ix, iy, iz;

	// does the atom box lookup or calculation
	CPotential::AtomBoxLookUp(val, Znum, x, y, z, B);

	/***************************************************************
	 * Do the trilinear interpolation
	 */
	val = complex_tt(0,0);
	if (x*x+y*y+z*z > m_atomRadius2) {
		return;
	}
	x = fabs(x);
	y = fabs(y);
	z = fabs(z);
	ix = (int)(x/m_ddx);
	iy = (int)(y/m_ddy);
	iz = (int)(z/m_ddz);
	dx = x-(float_tt)ix*m_ddx;
	dy = y-(float_tt)iy*m_ddy;
	dz = z-(float_tt)iz*m_ddz;
	if ((dx < 0) || (dy<0) || (dz<0)) {
		/* printf("Warning, dx(%g), dy(%g), dz(%g) < 0, (x=%g, y=%g, z=%g)\n",dx,dy,dz,x,y,z);
		 */
		if (dx < 0) dx = 0.0;
		if (dy < 0) dy = 0.0;
		if (dz < 0) dz = 0.0;
	}

	if (m_atomBoxes[Znum]->B > 0) {
		float_tt real = (1.0-dz)*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy].real()+
				dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy].real())+
				dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy+1].real()+
						dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy+1].real()))+
						dz*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy].real()+
								dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy].real())+
								dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy+1].real()+
										dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy+1].real()));
		float_tt imag = (1.0-dz)*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy].imag()+
				dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy].imag())+
				dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy+1].imag()+
						dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy+1].imag()))+
						dz*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy].imag()+
								dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy].imag())+
								dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy+1].imag()+
										dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy+1].imag()));
		val = complex_tt(real,imag);
	}
	else {
		val = (1.0-dz)*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz][ix][iy]+
				dx*m_atomBoxes[Znum]->rpotential[iz][ix+1][iy])+
				dy*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz][ix][iy+1]+
						dx*m_atomBoxes[Znum]->rpotential[iz][ix+1][iy+1]))+
						dz*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz+1][ix][iy]+
								dx*m_atomBoxes[Znum]->rpotential[iz+1][ix+1][iy])+
								dy*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz+1][ix][iy+1]+
										dx*m_atomBoxes[Znum]->rpotential[iz+1][ix+1][iy+1]));
	}
}

bool C3DPotential::CheckAtomZInBounds(float_tt atomZ)
{
	/*
	 * c = the thickness of the current slab.
	 *
	 * if the z-position of this atom is outside the potential slab
	 * we won't consider it and skip to the next
	 */
	return ((atomZ + _config->Potential.AtomRadiusAngstrom < m_c) && (atomZ - _config->Potential.AtomRadiusAngstrom + _config->Model.sliceThicknessAngstrom >= 0));
}

void C3DPotential::AddAtomToSlices(atom& atom,
		float_tt atomX, float_tt atomY, float_tt atomZ)
{
	if (!_config->Potential.periodicZ && CheckAtomZInBounds(atomZ))
	{
		AddAtomRealSpace(atom, atomX, atomY, atomZ);
	}
}

void C3DPotential::_AddAtomRealSpace(atom &atom,
		float_tt atomBoxX, unsigned ix,
		float_tt atomBoxY, unsigned iy,
		float_tt atomBoxZ, unsigned iAtomZ)
{
	unsigned iz;
	complex_tt dPot;

	/* calculate the range which we have left to cover with z-variation */
	/* iRadZ is the number of slices (rounded up) that this atom
	 * will contribute to, given its current x,y-radius
	 */
	float_tt r2sqr = atomBoxX*atomBoxX + atomBoxY*atomBoxY;
	// TODO: iRadZ is also calculated in the base class, one level up.  Which is correct?
	unsigned iRadZ = (unsigned)(sqrt(m_atomRadius2-r2sqr)/m_sliceThicknesses[0]+1.0);
	/* loop through the slices that this atoms contributes to */
	for (int iaz=-m_iRadZ;iaz <=m_iRadZ;iaz++) {
		if (!_config->Potential.periodicZ) {
			if (iaz+iAtomZ < 0) {
				if (-iAtomZ <= m_iRadZ) iaz = -iAtomZ;
				else break;
				if (abs(iaz)>_config->Model.nSlices) break;
			}
			if (iaz+iAtomZ >= _config->Model.nSlices)        break;
		}
		atomBoxZ = (double)(iAtomZ+iaz+0.5)*m_sliceThicknesses[0]-atomBoxZ;
		/* shift into the positive range */
		iz = (iaz+iAtomZ+32*_config->Model.nSlices) % _config->Model.nSlices;
		/* x,y,z is the true vector from the atom center
		 * We can look up the proj potential at that spot
		 * using trilinear extrapolation.
		 */
		AtomBoxLookUp(dPot,atom.Znum,atomBoxX,atomBoxY,atomBoxZ, _config->Model.UseTDS ? 0 : atom.dw);
		// printf("access: %d %d %d\n",iz,ix,iy);
		//    unsigned idx=ix*m_ny+iy;
		m_trans1[iz][iy][ix] += dPot;
		//    m_trans[iz][idx]+=dPot;
		//m_trans[iz][ix][iy][0] += dPot[0];
		//m_trans[iz][ix][iy][1] += dPot[1];
	} /* end of for iaz=-iRadZ .. iRadZ */
}


void C3DPotential::CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z)
{
	z = atom->z + _config->Structure.zOffset;
}


} // end namespace QSTEM