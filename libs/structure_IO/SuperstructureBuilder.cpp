/*
 * SuperstructureCreator.cpp
 *
 *  Created on: Apr 20, 2015
 *      Author: philipp
 */

#include "SuperstructureBuilder.hpp"
#include <boost/format.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <string>
using boost::format;
using boost::algorithm::trim;
namespace QSTEM {

SuperstructureBuilder::SuperstructureBuilder(StructureReaderFactory sfac,const ConfigPtr c) : IStructureBuilder(sfac,c) {
	atRad = atRadf;
	covRad = covRadf;
	_nGrains = 0;
	_sfac = sfac;
	_filePath = boost::filesystem::path(c->Structure.structureFilename);
	_superCell = superCellBoxPtr(new superCellBox());
}

SuperstructureBuilder::~SuperstructureBuilder() {
	// TODO Auto-generated destructor stub
}
superCellBoxPtr SuperstructureBuilder::Build(){
	char outFileName[64],datFileName[64],moldyName[128];
	char *str;
	// atom *atomPtr;
	int g,nstart,j,ak,count;
	FILE *mainfp;
	double charge;
	int moldyFlag = 0;     // suppress creation of ..._moldy.in file
	int distPlotFlag = 0;  // suppress creation of disList.dat


	/* Let's set the radii of certain elements by hand:*/
	atRad[39-1]=2.27; covRad[39-1]=2.62;  // Y
	atRad[57-1]=2.74; covRad[57-1]=1.69;  // La
	atRad[71-1]=2.25; covRad[71-1]=1.56;  // Lu


	_c->Structure.nCellX = 1;
	_c->Structure.nCellY = 1;
	_c->Structure.nCellZ = 1;

	BOOST_LOG_TRIVIAL(info) << format("Starting to build superstructure...");

	if (!readParams(_c->Structure.structureFilename.c_str()))
		exit(0);
	if (_nGrains > 0) {
		// loop, until we find crystalline grains:
		for (g=0; g<_nGrains; g++) if (grains[g].amorphFlag == CRYSTALLINE) break;
		/* make the crystalline part of the super cell */
		if (g<_nGrains) makeSuperCell();


		/* if there is also an amorphous part, then add it now, and also write a
		 * MD input file, so that the amorphous phase atoms can be relaxed
		 */
		for (g=0; g<_nGrains; g++) if (grains[g].amorphFlag != CRYSTALLINE) break;
		if (g<_nGrains) {
			if (moldyFlag) {
				sprintf(moldyName,"%s",datFileName);
				moldyName[strlen(datFileName)-4] = '\0';
				strcat(moldyName,"_moldy.in");
				if (( mainfp = fopen(moldyName,"w")) == NULL) {
					BOOST_LOG_TRIVIAL(info) << format("Could not open moldy input file %s!")%moldyName;
					exit(0);
				}
			}
			else {
				mainfp = NULL;
			}
			if (_nGrains > 1) {
				//					writeFrameWork(mainfp,superCell);
				computeCenterofMass();
				nstart = _superCell->atoms.size();
				switch (grains[g].amorphFlag) {
				case 1: makeAmorphous();
				break;
				case 2: makeSpecial(distPlotFlag);
				break;
				}
				//					writeAmorphous(mainfp,superCell,nstart,_superCell->atoms.size());
			}
			else {
				switch (grains[g].amorphFlag) {
				case 1: makeAmorphous();
				break;
				case 2: makeSpecial(distPlotFlag);
				break;
				}
				//					writeAmorphous( mainfp, superCell,0,_superCell->atoms.size());
			}
			if (moldyFlag)	fclose( mainfp );


		}
	}
	if (0) { // (moldyFlag) {
		///////////////////////////////////////////////////////////////
		// write Moldy input file, without presence of amorphous phase:
		sprintf(moldyName,"%s",datFileName);
		moldyName[strlen(datFileName)-4] = '\0';
		strcat(moldyName,"_moldy.in");
		if (( mainfp = fopen(moldyName,"w")) == NULL) {
			BOOST_LOG_TRIVIAL(info) << format("Could not open moldy input file %s!")%moldyName;
			exit(0);
		}
		// writeFrameWork(fp,superCell);
		// computeCenterofMass();
		// superCell2Moldy(fp,superCell);
		fclose(mainfp);
	} // end of: if moldyFlag ...
	strcpy(outFileName,datFileName);

	// atomPtr = readUnitCell(&natoms,fileName,&muls);
	// writePDB(atomPtr,nat  /* reset the input file and advance to the next crystal row */

	str = strchr(outFileName,'.');
	if (str == NULL) str=outFileName+strlen(outFileName);
	sprintf(str,".cfg");
	//		muls->ax = (float)superCell.ax;
	//		muls->by = (float)superCell.by;
	//		muls->c	= (float)superCell.cz;

	removeVacancies(_superCell->atoms);

	//		printf("will write cfg file to %s\n",outFileName);
	//		writeCFG(superCell.atoms, superCell.natoms, outFileName, muls);
	//		printf("wrote cfg file to %s\n",outFileName);

	/**************************************************************
	 * find the charge for the Y-atoms, in order to remain neutral:
	 */
	charge = 0.0;
	if (0) {
		//			for (ak=0;ak<muls->atomKinds;ak++) {
		//				count =0;
		//				for (j=0;j<superCell.natoms;j++) {
		//					if (muls->Znums[ak] == superCell.atoms[j].Znum) count++;
		//				}
		//				printf("Z=%3d: %d\n",muls->Znums[ak],count);
		//				switch (muls->Znums[ak]) {
		//				case  7: charge += count*(-3.0); break;
		//				case  8: charge += count*(-2.0);  break;
		//				case  38: charge += count*(2.0);  break;
		//				case  22: charge += count*(4.0);  break;
		//				case 14: charge += count*  4.0;  break;
		//				}
		//			}
	}
	else {
		for (j=0;j<_superCell->atoms.size();j++) {
			charge += _superCell->atoms[j].q*_superCell->atoms[j].occ;
		}
	}
	// BOOST_LOG_TRIVIAL(info) << format("Total charge: %g, i.e. %g %s\n",charge,charge,(charge > 0) ? "holes" : "electrons");
	BOOST_LOG_TRIVIAL(info) << format("Total charge: %g")%charge;
	if (charge > 0) BOOST_LOG_TRIVIAL(info) << format(", i.e. %g holes")%charge;
	if (charge < 0) BOOST_LOG_TRIVIAL(info) << format(", i.e. %g electrons")%-charge;

	return _superCell;
}
int SuperstructureBuilder::removeVacancies(std::vector<atom>& atoms){
	int natomsFinal;
	int i,i2,j,jz;
	float_tt totOcc,lastOcc, choice;
	long idum = -(long)time(NULL);
	int printLevel = 1;

	// printf("Time: %d\n",time(NULL));
	natomsFinal = atoms.size();
	qsort((void *)&atoms[0],atoms.size(),sizeof(atom),&CrystalBuilder::AtomCompareZYX);

	for(jz = 0,i=0;i<atoms.size();i++) {
		if (atoms[i].Znum > 0) {
			totOcc = atoms[i].occ;
			for (j=i+1;j<atoms.size();j++) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase
				// the total occupany and the counter j.
				if ((fabs(atoms[i].x-atoms[j].x) < 0.1) && (fabs(atoms[i].y-atoms[j].y) < 0.1) && (fabs(atoms[i].z-atoms[j].z) < 0.1)) {
					totOcc += atoms[j].occ;
				}
				else break;
			} // j-loop
			// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
			// do something about it:
			// printf("%d: tocc = %g\n",i,totOcc);
			if ((totOcc < 1) || (j > i+1)) { // found atoms at equal positions!
				// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values.
				//
				// if the total occupancy is less than 1 -> make sure we keep this
				// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
				if (totOcc < 1.0) choice = ran1();
				else choice = totOcc*ran1();
				lastOcc = 0;
				for (i2=i;i2<j;i2++) {
					// if choice does not match the current atom:
					// choice will never be 0 or 1(*totOcc)
					if ((choice <lastOcc) || (choice >=lastOcc+atoms[i2].occ)) {
						atoms[i2].Znum =  0;  // vacancy
						jz++;
					}
					lastOcc += atoms[i2].occ;
					// if we keep this atom!
					atoms[i2].occ =  1.0;  // to avoid this atom being reduced in occupancy again when reading the cfg file
				}
			}
			i = j-1;
		}
	}
	if ((jz > 0 ) &&(printLevel)) printf("Removed %d atoms because of occupancies < 1 or multiple atoms in the same place\n",jz);

	// end of managing atoms at equal position
	//////////////////////////////////////////////////////////////////////////////

	natomsFinal = atoms.size()-jz;
	// printf("%d atoms left\n",atoms.size()Final);
	// We now need to move all the atoms with zero Znum to the very end.
	if (jz > 0) {
		qsort(&atoms[0],atoms.size(),sizeof(atom),&CrystalBuilder::AtomCompareZnum);
		// for (i=0;i<atoms.size()Final;i++) printf("%3d: %d (%.1f %.1f %.1f): occ=%.1f\n",i,atoms[i].Znum,atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].occ);
	}
	return natomsFinal;
}
void SuperstructureBuilder::computeCenterofMass(){
	int i;
	double cmx=0.0,cmy=0.0,cmz=0.0;

	for (i=0;i<_superCell->atoms.size();i++) {
		cmx += _superCell->atoms[i].x;
		cmy += _superCell->atoms[i].y;
		cmz += _superCell->atoms[i].z;
	}
	_superCell->cmx = cmx/(_superCell->atoms.size()*_superCell->ax);
	_superCell->cmy = cmy/(_superCell->atoms.size()*_superCell->by);
	_superCell->cmz = cmz/(_superCell->atoms.size()*_superCell->cz);
	printf("Center of mass of frame work: (%g, %g, %g)\n",_superCell->cmx,_superCell->cmy,_superCell->cmz);
}
int SuperstructureBuilder::parOpen( const char *fileName )
{
	int i;
	BOOST_LOG_TRIVIAL(trace) << format("parOpen operating on: %s ") % fileName ;
	if ( fpStack.size() == 0 )
	{
		fpStack.resize(STACK_SIZE);
		for ( i = 0; i < STACK_SIZE; i++ )
		{
			fpStack[i] = NULL;
		}

		fpStack[0] = _fpParam;
	}

	if ( _fpParam != NULL )
	{
		printf( "DEBUG: RAM, fpParam is not NULL, stack not cleared correctly\n" );
		fclose( _fpParam );
	}

	_fpParam = fopen( fileName, "r" );
	return ( _fpParam != NULL );
}

int SuperstructureBuilder::readParams(const char *datFileName){
	int printFlag = 1;
	char title[128], parStr1[128],*str;
	std::string parStr;
	int gCount,i,Nkind;
	char unitCellFile[64];
	atom *tempCell;
	boost::filesystem::path unitCellFilePath;


	if (!parOpen(datFileName)) {
		BOOST_LOG_TRIVIAL(info) << format("Could not open data input file %s") %datFileName;
		return 0;
	}
	resetParamFile(_fpParam);
	while(readparam(_fpParam, "crystal:",parStr,0)) _nGrains++;
	resetParamFile(_fpParam);
	while(readparam(_fpParam, "amorph:",parStr,0)) _nGrains++;
	resetParamFile(_fpParam);
	while(readparam(_fpParam, "special:",parStr,0)) _nGrains++;
	BOOST_LOG_TRIVIAL(info) << format("Found data for %d grain(s) (crystalline frame work and amorphous)")%_nGrains;
	if (_nGrains == 0) return 0;

	grains.resize(_nGrains);
	/* Now we will loop through all the grains and lok for the necessary
	 * data for each grain
	 */
	if (readparam(_fpParam, "box:",parStr,1)) {
//		BOOST_LOG_TRIVIAL(info) << parStr;
		sscanf(parStr.c_str(),"%lf %lf %lf",&(_superCell->ax),&(_superCell->by),&(_superCell->cz));
	}
	else {
		BOOST_LOG_TRIVIAL(info) << format("Size of super cell box not defined - exit!");
		exit(0);
	}

	/* reset the input file and advance to the next crystal row */
	resetParamFile(_fpParam);
	gCount = -1;
	/* We will look for the following tokens:
	 * tilt: tiltx,tilty,tiltz;
	 * translation: shiftx shifty shiftz
	 * plane: vectX vectY vectZ pointX pointY pointZ
	 */
	while (readNextParam(_fpParam, title,parStr1)) {
		// printf("%s\n",parStr);
		/* if we found a new crystal ... */
		if (strncmp(title,"crystal:",8) == 0) {
			gCount++;
			grains[gCount].amorphFlag = 0;
			grains[gCount].density = 0;
			grains[gCount].rmin = 0;
			grains[gCount].rFactor = 1.0;
			grains[gCount].sphereRadius = 0;
			/* first we want to extract crystal name and file name for
			 * unit cell data file from this one line
			 */
			grains[gCount].name = string(parStr1);
			str = strnext(parStr1,(char*)" \t");
			if (str != NULL) {
				// printf("%s, name: %s\n",parStr,grains[gCount].name);
				//str = strnext(str," \t");
				//if (str != NULL) {
				strcpy(unitCellFile,str);
				if ((str = strchr(unitCellFile,' ')) != NULL)
					*str = '\0';
				if ((str = strchr(unitCellFile,'\t')) != NULL)
					*str = '\0';
			}
			else {
				BOOST_LOG_TRIVIAL(error) << format("Error: no unit cell data file specified for crystal %s") %
						grains[gCount].name;
				return 0;
			}

			grains[gCount].nplanes = 0;
			grains[gCount].planes = NULL;
			boost::filesystem::path p(unitCellFile);
			_c->Structure.structureFilename = p.string();
			CrystalBuilder cryst(_sfac,_c);
			cryst.ReadFromFile();
			auto atoms = cryst.GetUnitCellAtoms();
			auto unique = cryst.GetUniqueAtoms();
			for(int Znum : unique){
				bool thisZexists = (std::find(_superCell->uniqueatoms.begin(), _superCell->uniqueatoms.end(), Znum)) != _superCell->uniqueatoms.end();
				if(false == thisZexists) _superCell->uniqueatoms.push_back(Znum);
			}
			std::copy(atoms.begin(), atoms.end(), std::back_inserter(grains[gCount].unitCell));
			cryst.GetCellAngles(grains[gCount].alpha,grains[gCount].beta,grains[gCount].gamma);
			cryst.GetCellParameters(grains[gCount].ax, grains[gCount].by, grains[gCount].cz);
			grains[gCount].M = cryst.GetCellMatrix();
		}
		/***************************************************
		 * amorphous stuff
		 */
		else if (strncmp(title,"amorph:",7) == 0) {
			gCount++;
			grains[gCount].amorphFlag = AMORPHOUS;
			grains[gCount].density = 0;
			grains[gCount].rmin = 1000;
			grains[gCount].rFactor = 1.2;
			grains[gCount].sphereRadius = 0;
			/* first we want to extract crystal name and file name for
			 * unit cell data file from this one line
			 */
			grains[gCount].name = string(parStr1);
			str = strnext(parStr1,reinterpret_cast<const char*>(" \t"));
			if (str != NULL) {
				// printf("%s, name: %s\n",parStr,grains[gCount].name);
				//str = strnext(str," \t");
				//if (str != NULL) {
				strcpy(unitCellFile,str);
				if ((str = strchr(unitCellFile,' ')) != NULL)
					*str = '\0';
				if ((str = strchr(unitCellFile,'\t')) != NULL)
					*str = '\0';
			}
			else {
				BOOST_LOG_TRIVIAL(info) << format("Error: no unit cell data file specified for crystal %s") % grains[gCount].name;
				return 0;
			}

			// sscanf(parStr,"%s %s",grains[gCount].name,unitCellFile);
			grains[gCount].nplanes = 0;
			grains[gCount].planes = NULL;

			boost::filesystem::path p(unitCellFile);
			_c->Structure.structureFilename = p.string();
			CrystalBuilder cryst(_sfac,_c);
			cryst.SetFillUnitCell(true);
			cryst.ReadFromFile();
			auto atoms = cryst.GetUnitCellAtoms();
			auto unique = cryst.GetUniqueAtoms();
			for(int Znum : unique){
				bool thisZexists = (std::find(_superCell->uniqueatoms.begin(), _superCell->uniqueatoms.end(), Znum)) != _superCell->uniqueatoms.end();
				if(false == thisZexists) _superCell->uniqueatoms.push_back(Znum);
			}
			std::copy(atoms.begin(), atoms.end(), std::back_inserter(grains[gCount].unitCell));
			grains[gCount].alpha = 0;
			grains[gCount].beta  = 0;
			grains[gCount].gamma = 0;
			grains[gCount].ax = 0;
			grains[gCount].by = 0;
			grains[gCount].cz = 0;
		}
		/***************************************************
		 * code for specially distributed amorphous stuff
		 */
		else if (strncmp(title,"special:",8) == 0) {
			gCount++;
			grains[gCount].amorphFlag = SPECIAL_GRAIN;
			grains[gCount].density = 0;
			grains[gCount].rmin = 1000;
			grains[gCount].rFactor = 1.2;
			grains[gCount].sphereRadius = 0;
			/* first we want to extract crystal name and file name for
			 * unit cell data file from this one line
			 */
			grains[gCount].name = string(parStr1);
			grains[gCount].unitCell = std::vector<atom>();
			grains[gCount].nplanes = 0;
			grains[gCount].planes = NULL;
			grains[gCount].alpha = 0;
			grains[gCount].beta  = 0;
			grains[gCount].gamma = 0;
			grains[gCount].ax = 0;
			grains[gCount].by = 0;
			grains[gCount].cz = 0;
		}  // end of "special"
		/* if we found tilt data */
		else if (gCount >= 0) {
			if (strncmp(title,"tilt:",5) == 0) {
				sscanf(parStr1,"%lf %lf %lf",&(grains[gCount].tiltx),
						&(grains[gCount].tilty),&(grains[gCount].tiltz));
				if (strstr(parStr1,"degree") != NULL) {
					grains[gCount].tiltx *= PI180;
					grains[gCount].tilty *= PI180;
					grains[gCount].tiltz *= PI180;
				}
			}
			/* assign density */
			else if (strncmp(title,"density:",8) == 0) {
				sscanf(parStr1,"%lf",&(grains[gCount].density));
				grains[gCount].rmin = pow(sqrt(15.0/144.0)/grains[gCount].density,1.0/3.0);
			}
			/* assign density factor */
			else if (strncmp(title,"rmin:",5) == 0) {
				sscanf(parStr1,"%lf",&(grains[gCount].rmin));
				grains[gCount].density = sqrt(15.0/144)*pow(grains[gCount].rmin,3);
			}
			else if (strncmp(title,"r-factor:",9) == 0) {
				sscanf(parStr1,"%lf",&(grains[gCount].rFactor));
			}
			/* if we found shift data */
			else if (strncmp(title,"translation:",12) == 0) {
				sscanf(parStr1,"%lf %lf %lf",&(grains[gCount].shiftx),
						&(grains[gCount].shifty),&(grains[gCount].shiftz));
			}
			else if (strncmp(title,"sphere:",7) == 0) {
				sscanf(parStr1,"%lf %lf %lf %lf",&(grains[gCount].sphereRadius),
						&(grains[gCount].sphereX),&(grains[gCount].sphereY),&(grains[gCount].sphereZ));
			}
			/* if we found a new plane for this crystal */
			else if (strncmp(title,"plane:",6) == 0) {
				grains[gCount].nplanes++;
				grains[gCount].planes = (plane *) realloc(grains[gCount].planes,
						grains[gCount].nplanes*sizeof(plane));
				if (grains[gCount].planes == NULL) {
					printf("Sorry, could not allocate memory for new plane\n");
					return 0;
				}
				sscanf(parStr1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
						&(grains[gCount].planes[grains[gCount].nplanes-1].pointX),
						&(grains[gCount].planes[grains[gCount].nplanes-1].pointY),
						&(grains[gCount].planes[grains[gCount].nplanes-1].pointZ),
						&(grains[gCount].planes[grains[gCount].nplanes-1].vect1X),
						&(grains[gCount].planes[grains[gCount].nplanes-1].vect1Y),
						&(grains[gCount].planes[grains[gCount].nplanes-1].vect1Z),
						&(grains[gCount].planes[grains[gCount].nplanes-1].vect2X),
						&(grains[gCount].planes[grains[gCount].nplanes-1].vect2Y),
						&(grains[gCount].planes[grains[gCount].nplanes-1].vect2Z));
				crossProduct(&(grains[gCount].planes[grains[gCount].nplanes-1].vect1X),
						&(grains[gCount].planes[grains[gCount].nplanes-1].vect2X),
						&(grains[gCount].planes[grains[gCount].nplanes-1].normX));
			} /* end of if plane ... */
			else if (strncmp(title,"atom:",5) == 0) {
				//				grains[gCount].natoms++;
				//				grains[gCount].unitCell = (atom *)realloc(grains[gCount].unitCell,
				//					grains[gCount].natoms*sizeof(atom));
				//				// Now read Znum, r (->z), and count (->y)
				//				sscanf(parStr1,"%d %f %d %f",&(grains[gCount].unitCell[grains[gCount].natoms-1].Znum),
				//					&(grains[gCount].unitCell[grains[gCount].natoms-1].z),&Nkind,
				//					&(grains[gCount].unitCell[grains[gCount].natoms-1].y));
				//
				//				// assign number of atoms directly, if specified in input file
				//				if (Nkind > 0) grains[gCount].unitCell[grains[gCount].natoms-1].y = (float)Nkind;
				//				for (i=0;i<superCrystal->GetNumberOfAtomTypes();i++)
				//					if (superCrystal->GetAtomTypes()[i] == grains[gCount].unitCell[grains[gCount].natoms-1].Znum) break;
				//				if (i == superCrystal->GetNumberOfAtomTypes()) {
				//					superCrystal->GetNumberOfAtomTypes()++;
				//					superCrystal->GetAtomTypes() = (int *) realloc(superCrystal->GetAtomTypes(),superCrystal->GetNumberOfAtomTypes()*sizeof(int));
				//					superCrystal->GetAtomTypes()[i] = grains[gCount].unitCell[grains[gCount].natoms-1].Znum;
				//				}
			} /* end of if "atom:" */
		} /* end of if gCount >=0 */
	}
	parClose();
	if (printFlag) DisplayParams();
	return 1;
}
void SuperstructureBuilder::parClose()
{
	if ( _fpParam != NULL )
		fclose( _fpParam );
	_fpParam = NULL;
}

void SuperstructureBuilder::makeSuperCell(){
	int g,p,iatom,ix,iy,iz,atomCount = 0;
	static std::vector<float_tt> axCell(3),byCell(3),czCell(3);
	FloatArray2D Mm(boost::extents[3][3]),
			Mminv(boost::extents[3][3]),
			Mrot(boost::extents[3][3]),
			Mr(boost::extents[3][3]),
			Mr2(boost::extents[3][3]);
	arma::rowvec a(3),b(3);
	float_tt maxLength,dx,dy,dz,d,dxs,dys,dzs;
	atom newAtom;
	// float_tt xpos,ypos,zpos;
	int nxmin = 20000,nxmax=0,nymin=20000,nymax=0,nzmin=20000,nzmax=0;

	/* claculate maximum length in supercell box, which is naturally
	 * the room diagonal:
	 */
	//	maxLength = vectLength(&(_superCell->ax));
	maxLength = sqrt(_superCell->ax*_superCell->ax+
			_superCell->by*_superCell->by+
			_superCell->cz*_superCell->cz);


	atomCount = _superCell->atoms.size();
	for (g=0;g<_nGrains;g++) {
		/********************************************************
		 * if this grain is a crystalline one ...
		 */
		if (grains[g].amorphFlag == 0) {
			dx = grains[g].shiftx/_superCell->ax;
			dy = grains[g].shifty/_superCell->by;
			dz = grains[g].shiftz/_superCell->cz;
			/* find the rotated unit cell vectors .. */
			makeCellVect(grains[g], axCell, byCell, czCell);
			// showMatrix(Mm,3,3,"M");
			///////////////////////////////////////////////////////////////////
			float_tt	phi_x = grains[g].tiltx;
			float_tt	phi_y = grains[g].tilty;
			float_tt	phi_z = grains[g].tiltz;
			float_tt cx = cos(phi_x),sx = sin(phi_x),
					cy = cos(phi_y),    sy = sin(phi_y),
					cz = cos(phi_z)  ,  sz = sin(phi_z);
			arma::mat Mx = {1,0,0,0,cx,-sx,0,sx,cx};
			arma::mat My = {cy,0,-sy,0,1,0,sy,0,cy};
			arma::mat Mz = {cz,sz,0,-sz,cz,0,0,0,1};
			Mx.reshape(3,3);
			My.reshape(3,3);
			Mz.reshape(3,3);
			arma::mat Mrot = Mx*My*Mz;
			arma::mat M = grains[g].M;
			arma::mat Minv = arma::inv(M);
			arma::mat Mnew = M*Mrot;
			arma::mat MnewInv = arma::inv(Mnew);

//			std::cout << "Mrot: "<< endl<< Mrot << std::endl;
//			std::cout << "M: "<< endl<< M << std::endl;
//			std::cout << "Minv: "<< endl<< Minv << std::endl;
//			std::cout << "Mnew: "<< endl<< Mnew << std::endl;
//			std::cout << "MnewInv: "<< endl<< MnewInv << std::endl;

			for (ix=0;ix<=1;ix++) for (iy=0;iy<=1;iy++) for (iz=0;iz<=1;iz++) {
				a[0]=ix*_superCell->ax;
				a[1]=iy*_superCell->by;
				a[2]=iz*_superCell->cz;

				b = a*MnewInv;
//				std::cout << "b: "<< endl<< b << std::endl;
				if (nxmin > (int)floor(b[0]-dx)) nxmin=(int)floor(b[0]-dx);
				if (nxmax < (int)ceil( b[0]-dx)) nxmax=(int)ceil( b[0]-dx);
				if (nymin > (int)floor(b[1]-dy)) nymin=(int)floor(b[1]-dy);
				if (nymax < (int)ceil( b[1]-dy)) nymax=(int)ceil( b[1]-dy);
				if (nzmin > (int)floor(b[2]-dz)) nzmin=(int)floor(b[2]-dz);
				if (nzmax < (int)ceil( b[2]-dz)) nzmax=(int)ceil( b[2]-dz);
			}
			// nxmin--;nxmax++;nymin--;nymax++;nzmin--;nzmax++;
			//			superCell.atoms = (atom *)realloc(superCell.atoms,(superCell.natoms+1+
			//				(nxmax-nxmin)*(nymax-nymin)*
			//				(nzmax-nzmin)*grains[g].natoms)*
			//				sizeof(atom));
			_superCell->atoms.resize(_superCell->atoms.size()+1+ (nxmax-nxmin)*(nymax-nymin)* (nzmax-nzmin)*grains[g].unitCell.size());
			// showMatrix(Mm,3,3,"Mm");
			// showMatrix(Mminv,3,3,"Mminv");
			printf("Grain %d: range: (%d..%d, %d..%d, %d..%d)\n", g,nxmin,nxmax,nymin,nymax,nzmin,nzmax);

			dx = grains[g].shiftx;
			dy = grains[g].shifty;
			dz = grains[g].shiftz;
			for (iatom=0;iatom<grains[g].unitCell.size();iatom++) {
				memcpy(&newAtom,&(grains[g].unitCell[iatom]),sizeof(atom));
				// We need to convert the cartesian coordinates of this atom
				// to fractional ones:
				b[0] = newAtom.x;
				b[1] = newAtom.y;
				b[2] = newAtom.z;
				a = b*Minv;
				newAtom.x = a[0];
				newAtom.y = a[1];
				newAtom.z = a[2];
				//printf("%2d: %d (%g,%g,%g) (%g,%g,%g)\n",iatom,newAtom.Znum,
				//		 b[0][0],b[0][1],b[0][2],a[0][0],a[0][1],a[0][2]);
				for (ix=nxmin;ix<=nxmax;ix++) {
					for (iy=nymin;iy<=nymax;iy++) {
						for (iz=nzmin;iz<=nzmax;iz++) {
							/* atom position in reduced coordinates: */
							// a[0][0] = ix+newAtom.x; a[0][1] = iy+newAtom.y; a[0][2] = iz+newAtom.z;
							a[0] = newAtom.x+ix;
							a[1] = newAtom.y+iy;
							a[2] = newAtom.z+iz;
//							matrixProduct(a,1,3,Mm,3,3,b);
							b = a*M;
							_superCell->atoms[atomCount].x  = b[0]+dx;
							_superCell->atoms[atomCount].y  = b[1]+dy;
							_superCell->atoms[atomCount].z  = b[2]+dz;
							if ((_superCell->atoms[atomCount].x >= 0) &&
									(_superCell->atoms[atomCount].x < _superCell->ax) &&
									(_superCell->atoms[atomCount].y >= 0) &&
									(_superCell->atoms[atomCount].y < _superCell->by) &&
									(_superCell->atoms[atomCount].z >= 0) &&
									(_superCell->atoms[atomCount].z < _superCell->cz)) {
								// If this is a sphere:
								if (grains[g].sphereRadius > 0) {
									dxs = _superCell->atoms[atomCount].x - grains[g].sphereX;
									dys = _superCell->atoms[atomCount].y - grains[g].sphereY;
									dzs = _superCell->atoms[atomCount].z - grains[g].sphereZ;
									if (dxs*dxs+dys*dys+dzs*dzs < grains[g].sphereRadius*grains[g].sphereRadius) {
										_superCell->atoms[atomCount].dw = newAtom.dw;
										_superCell->atoms[atomCount].occ = newAtom.occ;
										_superCell->atoms[atomCount].q = newAtom.q;
										_superCell->atoms[atomCount].Znum = newAtom.Znum;
										atomCount++;
									}
								}
								// If this is a straight-edged grain
								else {
									for (p=0;p<grains[g].nplanes;p++) {
										/*
											printf("hello %d (%g %g %g)\n",g,
											_superCell->atoms[atomCount].x,_superCell->atoms[atomCount].y,
											_superCell->atoms[atomCount].z);
										 */
										d = findLambda(grains[g].planes+p,&(_superCell->atoms[atomCount].z),-1);

										/*
											printf("%3d lambda: %g (%g %g %g), (%g %g %g), %d\n",atomCount,d,
											newAtom.x,newAtom.y,newAtom.z,
											_superCell->atoms[atomCount].x,_superCell->atoms[atomCount].y,
											_superCell->atoms[atomCount].z,grains[g].nplanes);
										 */
										if (d < 0)
											break;
									}
									/* if all the previous test have been successful, this atom is IN,
									 * which means that we also need to copy the other data elements
									 * for this atom.
									 */
									if (p == grains[g].nplanes) {
										_superCell->atoms[atomCount].q = newAtom.q;
										_superCell->atoms[atomCount].dw = newAtom.dw;
										_superCell->atoms[atomCount].occ = newAtom.occ;
										_superCell->atoms[atomCount].Znum = newAtom.Znum;
										atomCount++;
									}
								} // if this is a sphere or not ...
							}
						} /* iz ... */
					} /* iy ... */
				} /* ix ... */
			} /* iatom ... */
			//			_superCell->natoms = atomCount;
			_superCell->atoms.resize(atomCount);
		} /* end of if !amorph,i.e. crystalline */
	} /* g=0..nGrains .. */
	BOOST_LOG_TRIVIAL(info) << format("Supercell contains %d atoms.") % _superCell->atoms.size();
}
std::vector<int> SuperstructureBuilder::GetUniqueAtoms(){

}
void SuperstructureBuilder::DisplayParams(){
	int g,p;

	printf("Supercell dimensions: %g x %g x %g A\n",
			_superCell->ax,_superCell->by,_superCell->cz);

	for (g=0;g<_nGrains;g++) {
		if (grains[g].amorphFlag == 0) {
			BOOST_LOG_TRIVIAL(info) << format("%d: Grain %s (ax=%g, by=%g, cz=%g, alpha=%g, beta=%g, gamma=%g)") %
					g%grains[g].name%grains[g].ax%grains[g].by%grains[g].cz%
					grains[g].alpha%grains[g].beta%grains[g].gamma;
			BOOST_LOG_TRIVIAL(info) << format("%d atoms in unit cell")%grains[g].unitCell.size();
			BOOST_LOG_TRIVIAL(info) << format("tilt=(%g, %g, %g)rad shift=(%g, %g, %g)A") %
					grains[g].tiltx%grains[g].tilty%grains[g].tiltz%
					grains[g].shiftx%grains[g].shifty%grains[g].shiftz;
		}
		else {
			BOOST_LOG_TRIVIAL(info) << format("%d: Grain %s (density=%g, rmin=%g, r-factor=%g")%
					g%grains[g].name%grains[g].density%grains[g].rmin%grains[g].rFactor;
		}
		printf("planes:\n");
		for (p=0;p<grains[g].nplanes;p++)
			printf("vect1=(%g %g %g) vect2=(%g %g %g) point=(%g %g %g) normal=(%g %g %g)\n",
					grains[g].planes[p].vect1X,grains[g].planes[p].vect1Y,
					grains[g].planes[p].vect1Z,
					grains[g].planes[p].vect2X,grains[g].planes[p].vect2Y,
					grains[g].planes[p].vect2Z,
					grains[g].planes[p].pointX,grains[g].planes[p].pointY,
					grains[g].planes[p].pointZ,
					grains[g].planes[p].normX,grains[g].planes[p].normY,
					grains[g].planes[p].normZ);
	} // for g=0 ...
}
void SuperstructureBuilder::SetSliceThickness(ModelConfig& mc){

}
void SuperstructureBuilder::SetResolution(ModelConfig& mc, const PotentialConfig pc){

}
void SuperstructureBuilder::makeAmorphous(){
	int g,p,iatom,ix,iy,iz,ic,atomCount = 0,amorphSites,amorphAtoms,randCount;
	static double *axCell,*byCell,*czCell=NULL;
	static double **Mm = NULL;
	double rCellx,rCelly,rCellz;
	double d,r;
	std::vector<atom> amorphCell;
	// atom newAtom;
	// double xpos,ypos,zpos;
	int nx,ny,nz;
	std::vector<int> randArray;


	if (Mm == NULL) {
		Mm = double2D(3,3,"Mm");
		axCell=Mm[0]; byCell=Mm[1]; czCell=Mm[2];
	}

	for (g=0;g<_nGrains;g++) {
		/********************************************************
		 * if this grain is an amorphous one ...
		 */
		if (grains[g].amorphFlag == AMORPHOUS) {
			r = grains[g].rmin/grains[g].rFactor;
			/* create an hexagonally closed packed unit cell for initial amorphous structure
			 * The length of each vector is now 1
			 */
			axCell[0] = r;     axCell[1] = 0;               axCell[2] = 0;
			byCell[0] = 0.5*r; byCell[1] = 0.5*sqrt(3.0)*r; byCell[2] = 0;
			czCell[0] = 0.5*r; czCell[1] = 0.5/sqrt(3.0)*r; czCell[2] = sqrt(5.0)/6*r;
			/* size of rectangular cell containing 4 atoms: */
			rCellx = r; rCelly = sqrt(3.0)*r; rCellz = (sqrt(5.0)*r)/3.0;

			/* determine number of unit cells in super cell */
			nx = (int)(_superCell->ax/rCellx);
			ny = (int)(_superCell->by/rCelly);
			nz = (int)(_superCell->cz/rCellz);
			amorphSites = 4*nx*ny*nz;
			amorphCell.resize(amorphSites+1);

			atomCount = 0;
			for (ix=0;ix<=nx;ix++) {
				for (iy=0;iy<=ny;iy++) {
					for (iz=0;iz<=nz;iz++) {
						for (ic=0;ic<4;ic++) {
							/* check if this atom and any of the 4 atoms per rect. unit cell lie within the super cell
							 */
							amorphCell[atomCount].x  = ix*rCellx-(ic==3)*axCell[0]+(ic % 2)*byCell[0]+(ic>1)*czCell[0];
							amorphCell[atomCount].y  = iy*rCelly-(ic==3)*axCell[1]+(ic % 2)*byCell[1]+(ic>1)*czCell[1];
							amorphCell[atomCount].z  = iz*rCellz-(ic==3)*axCell[2]+(ic % 2)*byCell[2]+(ic>1)*czCell[2];
							if ((amorphCell[atomCount].x >= 0) &&
									(amorphCell[atomCount].x < _superCell->ax) &&
									(amorphCell[atomCount].y >= 0) &&
									(amorphCell[atomCount].y < _superCell->by) &&
									(amorphCell[atomCount].z >= 0) &&
									(amorphCell[atomCount].z < _superCell->cz)) {
								for (p=0;p<grains[g].nplanes;p++) {
									d = findLambda(grains[g].planes+p,&(amorphCell[atomCount].z),-1);
									if (d < 0)
										break;
								}
								/* if all the previous test have been successful, this atom is IN */
								if (p == grains[g].nplanes) atomCount++;
							}
						} /* ic ... */
					} /* iz ... */
				} /* iy ... */
			} /* ix ... */
			amorphSites = atomCount;
			/****************************************************************************************
			 * Now we have all the sites within the bounding planes on which we can put atoms
			 */
			/* the true number of amorphous atoms is # of sites / rFactor^3 */
			amorphAtoms = (int)floor((double)amorphSites/pow(grains[g].rFactor,3.0));
			if (amorphAtoms > amorphSites) amorphAtoms = amorphSites;
			randCount = amorphSites;

			randArray.resize(amorphSites);
			for (ix=0;ix<amorphSites;ix++) randArray[ix] = ix;
			/*
			printf("memory allocation: sC.atoms: %d .. %d, rArray: %d .. %d\n",
			(int)_superCell->atoms,(int)_superCell->atoms+(atomCount+amorphAtoms+1)*sizeof(atom),
			(int)randArray,(int)randArray+amorphSites*sizeof(int));
			 */

			for (ix=amorphAtoms;ix>0;ix--) {
				do {
					iy = (int)((double)rand()*(double)(randCount-1)/(double)(RAND_MAX));
					if (iy >= randCount) iy = randCount-1;
					if (iy < 0) iy = 0;
					//	printf("%5d, iy: %d, sites: %d, atoms: %d  ",ix,iy,randCount,amorphAtoms);
					iz = randArray[iy];
					if (iz > amorphSites) {
						printf("%5d, iy: %d, sites: %d, atoms: %d  ",ix,iy,randCount,amorphAtoms);
						printf("iz: %d (%d)\n",iz,(int)_superCell->atoms[atomCount].z);
						printf("makeAmorphous: Error because of overlapping memory areas!\n");
						for (iz=0;iz<=amorphAtoms;iz++)
							printf("iz=%d: %d\n",iz,randArray[iz]);
						exit(0);
					}
				} while (iz > amorphSites);

				/* replace already choosen sites with unused ones, so that we don't occupy
				 * any site twice
				 */
				if (iy == randCount-1)  randCount--;
				else
					randArray[iy] = randArray[--randCount];


				iatom = ix % grains[g].unitCell.size();

				atom a;

				a.q = grains[g].unitCell[iatom].q;
				a.dw = grains[g].unitCell[iatom].dw;
				a.occ = grains[g].unitCell[iatom].occ;
				a.Znum = grains[g].unitCell[iatom].Znum;
				a.x = amorphCell[iz].x;
				a.y = amorphCell[iz].y;
				a.z = amorphCell[iz].z;
				_superCell->atoms.push_back(a);
			}
		} /* end of if amorph,i.e. crystalline */
	} /* g=0..nGrains .. */
}
#define SQR(x) ((x)*(x))
#define TRIAL_COUNT 1000000
void SuperstructureBuilder::makeSpecial(int distPlotFlag) {
	int p,g,iatom,atomCount = 0,amorphAtoms;
	double d,r,x,y,z,dist,volume;
	int i,j,Znum,count;
	long seed;
	float_tt pos[3],center[3],grainBound[6];
	int trials = 0,type;
	char *ptr;

	seed = -(long)(time(NULL));  // initialize random number generator.

	for (g=0;g<_nGrains;g++) {
		/********************************************************
		 * if this grain is a special one ...
		 */
		if (grains[g].amorphFlag == SPECIAL_GRAIN) {
			//  type = atoi(grains[g].name);
			trim(grains[g].name);
			type = std::stoi(grains[g].name);
			// printf("Distribution type: %d \n",type);
			BOOST_LOG_TRIVIAL(info) << format("%s: distribution type: %d (%s)")%grains[g].name,type%
					(type == 2) ? "double gaussian" : (type == 1 ? "single gaussian" : "random");
			/**************************************************************
			 * We would like to calculate the Volume of this grain.
			 * We do this by brute force, by simply checking, if randomly
			 * distributed points are within the grain, or not.
			 *************************************************************/
			grainBound[0] = _superCell->ax;  grainBound[1] = 0;
			grainBound[2] = _superCell->by;  grainBound[3] = 0;
			grainBound[4] = _superCell->cz;  grainBound[5] = 0;
			for (count=0, i=0; i<TRIAL_COUNT;i++) {
				// remember that we have to create a vector with
				// z,y,x, because that is the way the atom struct is
				pos[2] = _superCell->ax*ran1( );
				pos[1] = _superCell->by*ran1( );
				pos[0] = _superCell->cz*ran1( );
				for (p=0;p<grains[g].nplanes;p++) {
					d = findLambda((plane*)(grains[g].planes+p),pos,-1);
					if (d < 0) break;
				}
				// if all the previous tests have been successful, this atom is IN
				if (p == grains[g].nplanes) {
					count++;
					// center of this grain
					center[0] += pos[2]; center[1] += pos[1]; center[2] += pos[0];
					// boundaries in X-direction of this grain
					if (grainBound[0] > pos[2]) grainBound[0] = pos[2]; // xmin
					if (grainBound[1] < pos[2]) grainBound[1] = pos[2]; // xmax
					if (grainBound[2] > pos[1]) grainBound[2] = pos[1]; // ymin
					if (grainBound[3] < pos[1]) grainBound[3] = pos[1]; // ymax
					if (grainBound[4] > pos[0]) grainBound[4] = pos[0]; // zmin
					if (grainBound[5] < pos[0]) grainBound[5] = pos[0]; // zmax
				}
			}
			center[0] /= (double)count;
			center[1] /= (double)count;
			center[2] /= (double)count;
			volume = _superCell->ax*_superCell->by*_superCell->cz*(double)count/(double)TRIAL_COUNT;
			printf("Volume: %gA^3, %g %%\n",volume,(double)(100*count)/(double)TRIAL_COUNT);
			printf("boundaries: x: %g..%g, y: %g..%g, z: %g..%g\n",
					grainBound[0],grainBound[1],grainBound[2],
					grainBound[3],grainBound[4],grainBound[5]);

			// First we need to find out how much memory we need to reserve for this grain
			amorphAtoms = 0;
			for (iatom=0;iatom<grains[g].unitCell.size();iatom++) {
				if (grains[g].unitCell[iatom].y < 1.0) {  // if this is a concentration, and no count
					grains[g].unitCell[iatom].y *= volume;  // then convert it to number of atoms
				}
				amorphAtoms += (int)(grains[g].unitCell[iatom].y);
			}

			//			superCell.atoms = (atom *)realloc(superCell.atoms,(amorphAtoms+superCell.natoms+1)*
			//				sizeof(atom));
			//			if (superCell.atoms == NULL) {
			//				printf("makeAmorphous: Could not allocate memory for additional atoms!\n");
			//				exit(0);
			//			}
			atomCount = _superCell->atoms.size();  // start adding amorphous atoms, where we left off.

			// Now we can loop through and add these atoms randomly to the grain
			for (iatom=0;iatom<grains[g].unitCell.size();iatom++) {
				r = grains[g].unitCell[iatom].z;             // radius of this atom
				count = (int)(grains[g].unitCell[iatom].y);  // number of atoms of this kind
				Znum = grains[g].unitCell[iatom].Znum;
				covRad[Znum-1] = r;                          // set radius of other atoms also
				for (j=0;j<count;j++) {
					do { // make it lie within the grain bounding planes
						do { // make the atoms not touch eachother
							// z = superCell.cz*ran1();
							// y = superCell.by*ran1();
							z = grainBound[4]+ran1()*(grainBound[5]-grainBound[4]);
							y = grainBound[2]+ran1()*(grainBound[3]-grainBound[2]);
							if (fabs(_superCell->cz-z) < 2e-5) z = 0.0;
							if (fabs(_superCell->by-y) < 2e-5) y = 0.0;
							if (iatom > 0) {
								x = grainBound[0]+ran1()*(grainBound[1]-grainBound[0]);
							}
							else {
								switch (type) {
								case 0:
									x = grainBound[0]+ran1()*(grainBound[1]-grainBound[0]);
									break;
								case 1:
									x = xDistrFun1(center[0],0.08*(grainBound[1]-grainBound[0]));
									break;
								case 2:
									x = xDistrFun2(center[0],0.80*(grainBound[1]-grainBound[0]),
											0.08*(grainBound[1]-grainBound[0]));
									break;
								default:
									x = grainBound[0]+ran1()*(grainBound[1]-grainBound[0]);
								}
							}
							if (fabs(_superCell->ax-x) < 2e-5) x = 0.0;
							// Now we must check, whether atoms overlap
							for (i=0;i<atomCount;i++) {
								for (p=-1;p<=1;p++) {
									dist = sqrt(SQR(x-_superCell->atoms[i].x)+
											SQR(y-_superCell->atoms[i].y+p*_superCell->by)+
											SQR(z-_superCell->atoms[i].z));
									if (dist < r+covRad[_superCell->atoms[i].Znum-1]) break;
								}
								if (p < 2) break;
							}
							trials++;
							if (trials % amorphAtoms == 0)
								printf("Average trials per atom: %d times, success: %g %%\n",
										trials/amorphAtoms,100.0*(atomCount-_superCell->atoms.size())/
										(double)amorphAtoms);
						} while (i < atomCount);
						// try until we find one that does not touch any other

						// _superCell->atoms[atomCount].dw = 0.0;
						atom a;

						a.dw = 0.45*28.0/(double)(2.0*Znum);
						a.occ  = 1.0;
						a.q    = 0;
						a.Znum = Znum;
						a.x    = x;
						a.y    = y;
						a.z    = z;

						_superCell->atoms.push_back(a);

						for (p=0;p<grains[g].nplanes;p++) {
							d = findLambda(grains[g].planes+p,&(_superCell->atoms[atomCount].z),-1);
							if (d < 0) break;
						}
						// if all the previous tests have been successful, this atom is IN
					} while(p < grains[g].nplanes);
					atomCount++;
				} // for j=0..count
				printf("%d (%d): %d \n",iatom,Znum,count);
			} // for iatom = 0..natoms
			printf("\n%d amorphous atoms, volume: %gA^3 (%g%%), center: %g, width: %g\n",
					_superCell->atoms.size(),
					volume,100.0*volume/(_superCell->ax*_superCell->by*_superCell->cz),center[0],
					grainBound[1]-grainBound[0]);
			switch (type) {
			case 2:
				xDistrFun2(0.0,0.0,0.0);
				break;
			case 1:
				xDistrFun1(0.0,0.0);
				break;
			}
		} // if special_grain
	} // g=0..nGrains ..
	/*******************************************************************
	 * Now we must produce distribution plots of the different atom kinds
	 */
	// TODO implement makeDistrPlot whenever needed
	//	if (distPlotFlag) makeDistrPlot(superCell.atoms,superCell.natoms,superCell.ax);
}
/* This function creates the single gaussian distruted x-values for the dopand atoms
 */
double SuperstructureBuilder::xDistrFun1(double xcenter,double width) {
	static long idum = 0;
	static int count = 0;
	static double x2=0,xavg=0;
	double x;

	if (idum == 0) idum  = -(long)(time(NULL));

	if (width >0) {
		count ++;
		x = xcenter+width*gasdev();
		xavg += x;
		x2 += SQR(x-xcenter);
		return x;
	}
	else {
		printf("Statistics (%d): xavg=%g, x2dev=%g\n",count,xavg/(double)count,sqrt(x2/(double)count));
		xavg = 0;
		x2 = 0;
		return 0;
	}
}
/* This function creates the double gaussian distruted x-values for the dopand atoms
 * width1 i the distance between the 2 peaks
 * width 2 is the width of a single peak.
 */
double SuperstructureBuilder::xDistrFun2(double xcenter,double width1,double width2) {
	static long idum = 0;
	static long seed = 0;
	static int count = 0;
	static double x2=0,xavg=0;
	double x,dx;

	if (idum == 0) idum  = -(long)(time(NULL));
	if (seed == 0) seed  = -(long)(time(NULL));

	if (width2 >0) {
		count ++;
		dx = width2*gasdev();
		x = xcenter+width1*((ran1() > 0.5)-0.5)+dx;
		xavg += x;
		x2 += SQR(dx);
		return x;
	}
	else {
		printf("Statistics (%d): xavg=%g, x2dev=%g\n",count,xavg/(double)count,sqrt(x2/(double)count));
		xavg = 0;
		x2 = 0;
		return 0;
	}
}
// one can run "xmgr -nxy disList.dat &" to view the data produced by this function
#define DR 1.1
void makeDistrPlot(atom *atoms,int natoms,double ax) {
	//	int j,i,count,ind;
	//	int **list;
	//	FILE *pltfp;
	//
	//	printf("Atom kinds: %d: ",muls->atomKinds);
	//	for (i=0;i<muls->atomKinds;i++) printf(" %3d ",muls->Znums[i]);
	//	printf("\n");
	//
	//	count = (int)(ax/DR+1);
	//	list = int2D(muls->atomKinds,count,"list");
	//	memset(list[0],0,count*muls->atomKinds*sizeof(int));
	//	for (j=0;j<natoms;j++) {
	//		ind = (int)(atoms[j].x/DR);
	//		if (ind < 0) ind = 0;
	//		if (ind >= count) ind = count;
	//		for (i=0;i<muls->atomKinds;i++) if (muls->Znums[i] == atoms[j].Znum) break;
	//		if (i==muls->atomKinds) {
	//			// printf("Error: wrong Z (%d)\n",atoms[j].Znum);
	//		}
	//		else list[i][ind]++;
	//	}
	//	pltfp = fopen("disList.dat","w");
	//	for (j=0;j<count;j++) {
	//		fprintf( pltfp,"%.3f ",j*DR);
	//		for (i=0;i<muls->atomKinds;i++) fprintf( pltfp,"%d ",list[i][j]);
	//		fprintf( pltfp,"\n");
	//	}
	//	fclose(pltfp);
}
}
