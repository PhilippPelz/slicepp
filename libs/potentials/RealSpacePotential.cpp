/*
 * RealSpacePotential.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: philipp
 */

#include "RealSpacePotential.hpp"
using boost::format;
namespace QSTEM {

void RealSpacePotential::AddAtomRealSpace(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ) {
	// TODO add iatom in function parameters
//	unsigned iatom = atom - _structureBuilder->m_atoms.begin();
	unsigned iatom = 0;
	//	CenterAtomZ(atom, atomZ);

	/* Warning: will assume constant slice thickness ! */
	/* do not round here: atomX=0..dx -> iAtomX=0 */
	unsigned iAtomX = (int) floor(atomX / _config->Model.dx);
	unsigned iAtomY = (int) floor(atomY / _config->Model.dy);
	unsigned iAtomZ = (int) floor(atomZ / m_sliceThicknesses[0]);

	for (int iax = -m_iRadX; iax <= m_iRadX; iax++) {
		if (!_config->Potential.periodicXY) {
			if (iax + iAtomX < 0) {
				iax = -iAtomX;
				if (abs(iax) > m_iRadX)
					break;
			}
			if (iax + iAtomX >= _config->Model.nx)
				break;
		}
		float_tt x = (iAtomX + iax) * _config->Model.dx - atomX;
		unsigned ix = (iax + iAtomX + 16 * _config->Model.nx) % _config->Model.nx; /* shift into the positive range */
		for (int iay = -m_iRadY; iay <= m_iRadY; iay++) {
			if (!_config->Potential.periodicXY) {
				if (iay + iAtomY < 0) {
					iay = -iAtomY;
					if (abs(iay) > m_iRadY)
						break;
				}
				if (iay + iAtomY >= _config->Model.ny)
					break;
			}
			float_tt y = (iAtomY + iay) * _config->Model.dy - atomY;
			unsigned iy = (iay + iAtomY + 16 * _config->Model.ny) % _config->Model.ny; /* shift into the positive range */
			float_tt r2sqr = x * x + y * y;
			if (r2sqr <= m_atomRadius2) {
				// This (virtual) method is meant to be implemented by subclasses,
				//    for specific functionality varying by dimensionality.
				_AddAtomRealSpace(atom, x, ix, y, iy, atomZ, iAtomZ);
			}
		}
	}
}

RealSpacePotential::~RealSpacePotential() {
	// TODO Auto-generated destructor stub
}

} /* namespace QSTEM */
