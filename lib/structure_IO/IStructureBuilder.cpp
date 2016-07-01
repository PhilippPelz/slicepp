/*
 * IStructureBuilder.cpp
 *
 *  Created on: Apr 21, 2015
 *      Author: philipp
 */

#include "IStructureBuilder.hpp"

namespace slicepp {

IStructureBuilder::IStructureBuilder(StructureReaderFactory f,cStructureConfPtr sc, cModelConfPtr mc, cOutputConfPtr oc) {
	_f = f;
	_sc = sc;
	_mc = mc;
	_oc = oc;
//	_r->SetFilePath(sc->StructureFilename);
}

IStructureBuilder::~IStructureBuilder() {
	// TODO Auto-generated destructor stub
}

} /* namespace slicepp */
