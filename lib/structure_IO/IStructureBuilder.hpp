/*
 * IStructureBuilder.hpp
 *
 *  Created on: Apr 21, 2015
 *      Author: philipp
 */

#ifndef LIBS_STRUCTURE_IO_ISTRUCTUREBUILDER_HPP_
#define LIBS_STRUCTURE_IO_ISTRUCTUREBUILDER_HPP_

#include "config_IO/ConfigReader.hpp"
#include "structureInterface.hpp"

#include<list>

namespace slicepp {

class IStructureBuilder {
public:
	IStructureBuilder(StructureReaderFactory f,cStructureConfPtr sc, cModelConfPtr mc, cOutputConfPtr oc);
	virtual ~IStructureBuilder();
	virtual superCellBoxPtr Build()=0;
	virtual superCellBoxPtr DisplaceAtoms()=0;
	virtual std::vector<int> GetUniqueAtoms()=0;
	virtual std::map<unsigned, float_tt> GetU2()=0;
	virtual void DisplayParams()=0;
	// Register a function that is called on every atom when it has its final position in the sample
	void RegisterAtomVisitor(boost::function<void(const atom& a)>);
protected:
	IStructureBuilder(){};
	StructureReaderFactory _f;
	cStructureConfPtr _sc;
	cModelConfPtr _mc;
	cOutputConfPtr _oc;
};

typedef boost::shared_ptr<IStructureBuilder> StructureBuilderPtr;
typedef boost::function<IStructureBuilder*(StructureReaderFactory f,cStructureConfPtr sc, cModelConfPtr mc, cOutputConfPtr oc)> structureBuilderCreator;
typedef std::map<string,structureBuilderCreator> StructureBuilderFactory;

} /* namespace slicepp */

#endif /* LIBS_STRUCTURE_IO_ISTRUCTUREBUILDER_HPP_ */
