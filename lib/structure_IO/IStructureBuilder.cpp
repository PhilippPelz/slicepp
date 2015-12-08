/*
 * IStructureBuilder.cpp
 *
 *  Created on: Apr 21, 2015
 *      Author: philipp
 */

#include "IStructureBuilder.hpp"

namespace QSTEM {

IStructureBuilder::IStructureBuilder(StructureReaderPtr r,cStructureConfPtr sc, cModelConfPtr mc, cOutputConfPtr oc) {
	_r = r;
	_sc = sc;
	_mc = mc;
	_oc = oc;
	_r->SetFilePath(sc->structureFilename);
}

IStructureBuilder::~IStructureBuilder() {
	// TODO Auto-generated destructor stub
}

} /* namespace QSTEM */
