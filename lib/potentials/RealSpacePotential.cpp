/*
 * RealSpacePotential.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: philipp
 */

#include "RealSpacePotential.hpp"
using boost::format;
namespace slicepp {
/****************************************************************************
 * function: atomBoxLookUp - looks up potential at position x, y, z, relative to atom center
 *
 * Znum = element
 * x,y,z = real space position (in A)
 * B = Debye-Waller factor, B=8 pi^2 <u^2>
 ***************************************************************************/
void RealSpacePotential::AtomBoxLookUp(complex_tt &val, int Znum, float_tt x,
		float_tt y, float_tt z, float_tt B) {
	int boxNx, boxNy, boxNz;
	float_tt dx, dy, dz, ddx, ddy, ddz;
	float_tt maxRadius2;
	char fileName[256], systStr[256];
	int tZ, tnx, tny, tnz, tzOversample;
	float_tt tdx, tdy, tdz, tv0, tB;
	FILE *fp;
	int numRead = 0, dummy;

	/* initialize all the atoms to non-used */
	if (!m_atomBoxes.count(Znum)) {
		m_atomBoxes[Znum] = atomBoxPtr(new atomBox());
		m_atomBoxes[Znum]->B = -1.0;

		BOOST_LOG_TRIVIAL(trace)<<format("Atombox has real space resolution of %g x %g x %gA (%d x %d x %d pixels)")
				%ddx% ddy% ddz% boxNx% boxNy% boxNz;
	}
	// printf("Debugging: %d %g %g: %g",Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B));

	/* Creating/Reading a atombox for every new kind of atom, but only as needed */
	if (fabs(m_atomBoxes[Znum]->B - B) > 1e-6) {
		//  printf("Debugging 1 (%d: %.7g-%.7g= %.7g), %d",
		//	   Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B),fabs(m_atomBoxes[Znum]->B - B) > 1e-6);
		m_atomBoxes[Znum]->B = B;
		/* Open the file with the projected potential for this particular element
		 */
		sprintf(fileName, "potential_%d_B%d.prj", Znum, (int) (100.0 * B));
		if ((fp = fopen(fileName, "r")) == NULL) {
			sprintf(systStr, "scatpot %s %d %g %d %d %d %g %g %g %d %g",
					fileName, Znum, B, boxNx, boxNy, boxNz, ddx, ddy, ddz,
					OVERSAMPLINGZ, _mc->EnergykeV);
			BOOST_LOG_TRIVIAL(info)<<format("Could not find precalculated potential for Z=%d, will calculate now.")% Znum;
			BOOST_LOG_TRIVIAL(info)<<format("Calling: %s") %systStr;
			system(systStr);
			for (dummy = 0; dummy < 10000; dummy++)
				;
			if ((fp = fopen(fileName, "r")) == NULL) {
				BOOST_LOG_TRIVIAL(error)<<format("cannot calculate projected potential using scatpot - exit!");
				exit(0);
			}

		}
		fgets(systStr, 250, fp);
		sscanf(systStr, "%d %f %d %d %d %f %f %f %d %f", &tZ, &tB, &tnx, &tny,
				&tnz, &tdx, &tdy, &tdz, &tzOversample, &tv0);
		/* If the parameters in the file don't match the current ones,
		 * we need to create a new potential file
		 */
		if ((tZ != Znum) || (fabs(tB - B) > 1e-6) || (tnx != boxNx)
				|| (tny != boxNy) || (tnz != boxNz) || (fabs(tdx - ddx) > 1e-5)
				|| (fabs(tdy - ddy) > 1e-5) || (fabs(tdz - ddz) > 1e-5)
				|| (tzOversample != OVERSAMPLINGZ)
				|| (tv0 != _mc->EnergykeV)) {
			BOOST_LOG_TRIVIAL(info)<<format("Potential input file %s has the wrong parameters") %fileName;
			BOOST_LOG_TRIVIAL(info) << format( "Parameters:");
			BOOST_LOG_TRIVIAL(info) << format("file:    Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g")
					%tZ% tB% tnx% tny% tnz% tdx% tdy% tdz% tzOversample% tv0;
			BOOST_LOG_TRIVIAL(info) << format("program: Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g")
					%Znum% B% boxNx% boxNy% boxNz% ddx% ddy% ddz%(OVERSAMPLINGZ)% _mc->EnergykeV;
			BOOST_LOG_TRIVIAL(info) << format("will create new potential file, please wait ...");

			/* Close the old file, Create a new potential file now
			 */
			fclose(fp);
			sprintf(systStr, "scatpot %s %d %g %d %d %d %g %g %g %d %g",
					fileName, Znum, B, boxNx, boxNy, boxNz, ddx, ddy, ddz,
					OVERSAMPLINGZ, _mc->EnergykeV);
			system(systStr);
			if ((fp = fopen(fileName, "r")) == NULL) {
				BOOST_LOG_TRIVIAL(error)<<format("cannot calculate projected potential using scatpot - exit!");
				exit(0);
			}
			fgets(systStr, 250, fp);
		}

		/* Finally we can read in the projected potential
		 */
		if (B == 0) {
			m_atomBoxes[Znum]->rpotential.resize(boost::extents[boxNz][boxNx][boxNy]);
			numRead = fread(m_atomBoxes[Znum]->rpotential.data(),
					sizeof(float_tt), (size_t) (boxNx * boxNy * boxNz), fp);
		} else {
			m_atomBoxes[Znum]->potential.resize(boost::extents[boxNz][boxNx][boxNy]);
			//			= complex3D(boxNz, boxNx, boxNy,
			//					"atomBox");
			numRead = fread(m_atomBoxes[Znum]->potential.data(),
					sizeof(complex_tt), (size_t) (boxNx * boxNy * boxNz), fp);
		}

		/* writeImage_old(m_atomBoxes[Znum]->potential[0],boxNx,boxNy, 0.0,"potential.img");
		 system("showimage potential.img");
		 */
		fclose(fp);

		if (numRead == boxNx * boxNy * boxNz) {
			BOOST_LOG_TRIVIAL(info)<<format("Sucessfully read in the projected potential");
		} else {
			BOOST_LOG_TRIVIAL(error)<<format("error while reading potential file %s: read %d of %d values")
					%fileName% numRead%( boxNx * boxNy * boxNz);
			exit(0);
		}
	}
}
void RealSpacePotential::AddAtomRealSpace(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ) {
	// TODO add iatom in function parameters
//	unsigned iatom = atom - _structureBuilder->m_atoms.begin();
	unsigned iatom = 0;
	//	CenterAtomZ(atom, atomZ);

	/* Warning: will assume constant slice thickness ! */
	/* do not round here: atomX=0..dx -> iAtomX=0 */
	unsigned iAtomX = (int) floor(atomX / _mc->d[0]);
	unsigned iAtomY = (int) floor(atomY / _mc->d[1]);
	unsigned iAtomZ = (int) floor(atomZ / _sliceThicknesses[0]);

	for (int iax = -_nRadX; iax <= _nRadX; iax++) {
		if (!_mc->periodicXY) {
			if (iax + iAtomX < 0) {
				iax = -iAtomX;
				if (abs(iax) > _nRadX)
					break;
			}
			if (iax + iAtomX >= _mc->n[0])
				break;
		}
		float_tt x = (iAtomX + iax) * _mc->d[0] - atomX;
		unsigned ix = (iax + iAtomX + 16 * _mc->n[0]) % _mc->n[0]; /* shift into the positive range */
		for (int iay = -_nRadY; iay <= _nRadY; iay++) {
			if (!_mc->periodicXY) {
				if (iay + iAtomY < 0) {
					iay = -iAtomY;
					if (abs(iay) > _nRadY)
						break;
				}
				if (iay + iAtomY >= _mc->n[1])
					break;
			}
			float_tt y = (iAtomY + iay) * _mc->d[1] - atomY;
			unsigned iy = (iay + iAtomY + 16 * _mc->n[1]) % _mc->n[1]; /* shift into the positive range */
			float_tt r2sqr = x * x + y * y;
			if (r2sqr <= _atomRadius2) {
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

} /* namespace slicepp */
