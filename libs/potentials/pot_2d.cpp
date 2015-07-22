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

#include "pot_2d.hpp"

namespace QSTEM
{

C2DPotential::C2DPotential(const ConfigPtr& c,const PersistenceManagerPtr& p)  : RealSpacePotential(c, p)
{
	m_boxNz = 1;
}

void C2DPotential::SliceSetup() {
CPotential::SliceSetup();
}
void  C2DPotential::ComputeAtomPotential(int znum){

}
void C2DPotential::DisplayParams()
{
	CPotential::DisplayParams();
  printf("* Potential calculation: 2D (real space method)");
}

void C2DPotential::AtomBoxLookUp(complex_tt &sum, int Znum, float_tt x, float_tt y, float_tt z, float_tt B) 
{
  float_tt dx, dy;
  int ix, iy;
  RealSpacePotential::AtomBoxLookUp(sum, Znum, x, y, z, B);

  /***************************************************************
   * Do the trilinear interpolation
   */
  sum = 0.0;
//  sum[1] = 0.0;
  if (x*x+y*y+z*z > _atomRadius2) {
    return;
  }
  x = fabs(x);
  y = fabs(y);
  ix = (int)(x/_ddx);
  iy = (int)(y/_ddy);
  dx = x-(float_tt)ix*_ddx;
  dy = y-(float_tt)iy*_ddy;
  
  if ((dx < 0) || (dy<0) ) {
    /* printf("Warning, dx(%g), dy(%g), dz(%g) < 0, (x=%g, y=%g, z=%g)\n",dx,dy,dz,x,y,z);
     */
    if (dx < 0) dx = 0.0;
    if (dy < 0) dy = 0.0;
  }
  
  if (m_atomBoxes[Znum]->B > 0) {
    float_tt real = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy].real()+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy].real())+
                       dy*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy+1].real()+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy+1].real());
    float_tt imag = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy].imag()+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy].imag())+
                       dy*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy+1].imag()+
                       dx*m_atomBoxes[Znum]->potential[0][ix+1][iy+1].imag());
    sum = complex_tt(real,imag);
  }
  else {
    sum = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->rpotential[0][ix][iy]+
                       dx*m_atomBoxes[Znum]->rpotential[0][ix+1][iy])+
                       dy*((1.0-dx)*m_atomBoxes[Znum]->rpotential[0][ix][iy+1]+
		       dx*m_atomBoxes[Znum]->rpotential[0][ix+1][iy+1]);
  }
}
void C2DPotential::SaveAtomicPotential(int znum){

}
bool C2DPotential::CheckAtomZInBounds(float_tt atomZ)
{
  /*
   * c = the thickness of the current slab.
   *
   * if the z-position of this atom is outside the potential slab
   * we won't consider it and skip to the next
   */
  return ((atomZ<_totalThickness) && (atomZ>=0));
}


void C2DPotential::AddAtomToSlices(atom& atom )
{
  // Note that if you override this method, you should do the following check to make sure the atom is in bounds.
  // skip atoms that are beyond the cell's boundaries
	// TODO z periodic potenital
//  if (!_config->Potential.periodicZ)
//    {
//      if (atomZ > _config->Structure. || atomZ<0) return;
//    }

  // Calls parent class method, which in turn calls method below after computing ix, iy, iAtomZ
  AddAtomRealSpace(atom, atom.r[0], atom.r[1], atom.r[2]);
}

void C2DPotential::_AddAtomRealSpace(atom& atom,
                                     float_tt atomBoxX, unsigned int ix, 
                                     float_tt atomBoxY, unsigned int iy, 
                                     float_tt atomZ, unsigned int iAtomZ)
{
  complex_tt dPot;
  // Coordinates in the total potential space
  unsigned iz;
  

  if (!_c->Potential.periodicZ) {
    if (iAtomZ < 0) return;
    if (iAtomZ >= _c->Model.nSlices) return;
  }                
  iz = (iAtomZ+32*_c->Model.nSlices) % _c->Model.nSlices;         /* shift into the positive range */
  // x, y are the coordinates in the space of the atom box
  AtomBoxLookUp(dPot,atom.Znum,atomBoxX,atomBoxY,0, _c->Model.UseTDS ? 0 : atom.dw);
  float_tt atomBoxZ = (double)(iAtomZ+1)*_sliceThicknesses[0]-atomZ;

  unsigned idx=ix*_c->Model.ny+iy;

  /* split the atom if it is close to the top edge of the slice */
  if ((atomBoxZ<0.15*_sliceThicknesses[0]) && (iz >0)) {
//TODO use trans1    m_trans[iz][idx] += complex_tt(0.5*dPot.real(),0.5*dPot.imag());
//TODO use trans1    m_trans[iz-1][idx] += complex_tt(0.5*dPot.real(),0.5*dPot.imag());
  }
  /* split the atom if it is close to the bottom edge of the slice */
  else {
    if ((atomBoxZ>0.85*_sliceThicknesses[0]) && (iz < _c->Model.nSlices-1)) {
// TODO use trans1     m_trans[iz][idx] += complex_tt(0.5*dPot.real(),0.5*dPot.imag());
// TODO use trans1     m_trans[iz+1][idx] += complex_tt(0.5*dPot.real(),0.5*dPot.imag());
    }
    else {
// TODO use trans1      m_trans[iz][idx] += dPot;
    }
  }
}

void C2DPotential::CenterAtomZ(atom& atom, float_tt &z)
{
  CPotential::CenterAtomZ(atom, z);
  z += 0.5*_c->Model.dz;
}

} // end namespace QSTEM
