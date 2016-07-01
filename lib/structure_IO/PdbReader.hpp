/*
 * PdbReader.h
 *
 *  Created on: Feb 16, 2015
 *      Author: philipp
 */

#ifndef LIBS_STRUCTURE_IO_PdbREADER_HPP_
#define LIBS_STRUCTURE_IO_PdbREADER_HPP_

#include "structureInterface.hpp"
#include "stemtypes_fftw3.hpp"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

using namespace OpenBabel;

namespace slicepp {

class PdbReader: public IStructureReader {
public:
	PdbReader(const boost::filesystem::path &file);
	virtual ~PdbReader();

	virtual int ReadCellParams(FloatArray2D& Mm);
	virtual int ReadAtoms(std::vector<atom> &atoms,
			std::vector<int> &uniqueAtoms,bool fillUnitCell);
protected:
	int _atomsAfterFillUC, _atomsBeforeFillUC;

private:
  friend class CStructureReaderFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static StructureReaderPtr Create(const boost::filesystem::path &file){
    return StructureReaderPtr(new PdbReader(file));}
};
}
#endif /* LIBS_STRUCTURE_IO_PdbREADER_HPP_ */
