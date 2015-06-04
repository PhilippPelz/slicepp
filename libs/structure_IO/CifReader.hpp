/*
 * CifReader.h
 *
 *  Created on: Feb 16, 2015
 *      Author: philipp
 */

#ifndef LIBS_STRUCTURE_IO_CIFREADER_HPP_
#define LIBS_STRUCTURE_IO_CIFREADER_HPP_

#include "structureInterface.hpp"
#include "stemtypes_fftw3.hpp"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

using namespace OpenBabel;

namespace QSTEM {

class CifReader: public IStructureReader {
public:
	CifReader(const boost::filesystem::path &file);
	virtual ~CifReader();

	virtual int ReadCellParams(FloatArray2D& Mm);
	virtual int ReadAtoms(std::vector<atom> &atoms,
			std::vector<int> &uniqueAtoms,bool fillUnitCell);
protected:
//	OBUnitCell* _unitCell;
//	OBConversion _conversion;
//	OBMol _mol;
	int _atomsAfterFillUC, _atomsBeforeFillUC;

private:
  friend class CStructureReaderFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static StructureReaderPtr Create(const boost::filesystem::path &file){
    return StructureReaderPtr(new CifReader(file));}
};
}
#endif /* LIBS_STRUCTURE_IO_CIFREADER_HPP_ */
