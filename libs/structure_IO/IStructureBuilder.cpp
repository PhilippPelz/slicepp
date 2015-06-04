/*
 * IStructureBuilder.cpp
 *
 *  Created on: Apr 21, 2015
 *      Author: philipp
 */

#include "IStructureBuilder.hpp"

namespace QSTEM {

IStructureBuilder::IStructureBuilder(StructureReaderPtr& r,const ConfigPtr& c) {
	_r = r;
	_c = c;
	_r->SetFilePath(_c->Structure.structureFilename);
}

IStructureBuilder::~IStructureBuilder() {
	// TODO Auto-generated destructor stub
}

} /* namespace QSTEM */
