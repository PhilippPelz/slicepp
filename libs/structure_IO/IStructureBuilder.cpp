/*
 * IStructureBuilder.cpp
 *
 *  Created on: Apr 21, 2015
 *      Author: philipp
 */

#include "IStructureBuilder.hpp"

namespace QSTEM {

IStructureBuilder::IStructureBuilder(StructureReaderFactory& sfac,const ConfigPtr& c) {
	_sfac = sfac;
	_c = c;
}

IStructureBuilder::~IStructureBuilder() {
	// TODO Auto-generated destructor stub
}

} /* namespace QSTEM */
