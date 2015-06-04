/*
 * SuperstructureCreator.h
 *
 *  Created on: Apr 20, 2015
 *      Author: philipp
 */

#ifndef LIBS_STRUCTURE_IO_SUPERSTRUCTURECREATOR_H_
#define LIBS_STRUCTURE_IO_SUPERSTRUCTURECREATOR_H_

#include "matrixlib.hpp"
#include "stemtypes_fftw3.hpp"
#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <boost/filesystem.hpp>
#include <structure_IO/IStructureBuilder.hpp>
#include "random.hpp"
#include "readparams.hpp"
#include "crystal.hpp"

namespace QSTEM {
#define NAME_BUF_LEN 64
#define CRYSTALLINE 0
#define AMORPHOUS 1
#define SPECIAL_GRAIN 2
#define STACK_SIZE 5

class SuperstructureBuilder : public IStructureBuilder {
public:
	SuperstructureBuilder(StructureReaderPtr& r,const ConfigPtr& c);
	virtual ~SuperstructureBuilder();
	int removeVacancies(std::vector<atom>& atoms);
	void computeCenterofMass();
	int readParams(const char *datFileName);
	void makeSuperCell();
	void makeAmorphous();
	void makeSpecial(int distPlotFlag);
	double xDistrFun1(double xcenter,double width);
	double xDistrFun2(double xcenter,double width1,double width2);
	virtual superCellBoxPtr Build();
	virtual superCellBoxPtr DisplaceAtoms();
	virtual std::map<unsigned, float_tt> GetU2(){//TODO
		std::map<unsigned, float_tt> a;
		return 	a;};
	virtual std::vector<int> GetUniqueAtoms();
	virtual void DisplayParams();
	virtual void SetSliceThickness(ModelConfig& mc);
	virtual void SetResolution(ModelConfig& mc, const PotentialConfig pc);

	superCellBoxPtr _superCell;
	boost::filesystem::path _filePath;
	std::vector<float_tt> _atRadf,_covRadf,_atRad, _covRad;
	std::vector<grainBox> grains;
protected:
	FILE *_fpParam=NULL;
	std::vector<FILE*> fpStack;
	int parOpen( const char *fileName );
	void parClose();
	int _nGrains = 0;
};

} /* namespace QSTEM */

#endif /* LIBS_STRUCTURE_IO_SUPERSTRUCTURECREATOR_H_ */
