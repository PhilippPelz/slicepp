/*
 * IStructureBuilder.hpp
 *
 *  Created on: Apr 21, 2015
 *      Author: philipp
 */

#ifndef LIBS_STRUCTURE_IO_ISTRUCTUREBUILDER_HPP_
#define LIBS_STRUCTURE_IO_ISTRUCTUREBUILDER_HPP_

#include "config_IO/read_qsc.hpp"
#include "structureInterface.hpp"

namespace QSTEM {



class IStructureBuilder {
public:
	IStructureBuilder(StructureReaderFactory& sfac,const ConfigPtr& c);
	virtual ~IStructureBuilder();
	virtual superCellBoxPtr Build()=0;
	virtual std::vector<int> GetUniqueAtoms()=0;
	virtual std::map<unsigned, float_tt> GetU2()=0;
	virtual void DisplayParams()=0;
	virtual void SetSliceThickness(ModelConfig& mc)=0;
	virtual void SetResolution(ModelConfig& mc, const PotentialConfig pc)=0;
protected:
	IStructureBuilder(){};
	StructureReaderFactory _sfac;
	ConfigPtr _c;
};
typedef boost::shared_ptr<IStructureBuilder> StructureBuilderPtr;
typedef boost::function<IStructureBuilder*(StructureReaderFactory& sfac,const ConfigPtr& c) > structureBuilderCreator;
typedef std::map<string,structureBuilderCreator> StructureBuilderFactory;
} /* namespace QSTEM */

#endif /* LIBS_STRUCTURE_IO_ISTRUCTUREBUILDER_HPP_ */
