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
#include "memory_fftw3.hpp"
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
	SuperstructureBuilder(StructureReaderFactory& sfac,const ConfigPtr& c);
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
	virtual std::map<unsigned, float_tt> GetU2(){//TODO
		std::map<unsigned, float_tt> a;
		return 	a;};
	virtual std::vector<int> GetUniqueAtoms();
	virtual void DisplayParams();
	virtual void SetSliceThickness(ModelConfig& mc);
	virtual void SetResolution(ModelConfig& mc, const PotentialConfig pc);

	superCellBoxPtr _superCell;
	boost::filesystem::path _filePath;
	std::vector<float_tt> atRadf = {0.79,0.49, 2.05,1.40,1.17,0.91,0.75,0.65,0.57,0.51, 2.23,1.72,0.00,1.46,0.00,0.00,0.97,0.88, 2.77,2.23,2.09,2.00,1.92,1.85,1.79,1.72,1.67,1.62,1.57,1.53,1.81,1.52,1.33,1.22,1.12,1.03};
	std::vector<float_tt> covRadf = {0.32,0.93, 1.23,0.90,0.82,0.77,0.75,0.73,0.72,0.71, 1.54,1.36,0.00,1.90,0.00,0.00,0.99,0.98, 2.03,1.74,1.44,1.32,1.63,1.18,1.17,1.17,1.16,1.15,1.17,1.25,1.26,1.22,1.20,1.16,1.14,1.12};
	std::vector<float_tt> atRad, covRad;
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
