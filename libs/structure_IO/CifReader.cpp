/*
 * CifReader.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: philipp
 */

#include "CifReader.hpp"
#include <openbabel/math/vector3.h>
#include <openbabel/op.h>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
using namespace std;
using namespace OpenBabel;
using boost::format;

namespace QSTEM {
CifReader::CifReader(const boost::filesystem::path &file) {
	_filepath = file.string();

}

CifReader::~CifReader() {
	// TODO Auto-generated destructor stub
}

int CifReader::ReadCellParams(FloatArray2D& Mm) {
	ifstream in(_filepath);

	stringstream out;
	OBConversion _conversion(&in, &out);
	OBOp* pOp = OBOp::FindType("fillUC");

	if (_conversion.SetInFormat("CIF")) {
		OBMol _mol;

		_conversion.ReadFile(&_mol, _filepath);
		OBUnitCell* _unitCell = (OBUnitCell*)_mol.GetData(OBGenericDataType::UnitCell);
		_atomsBeforeFillUC = _mol.NumAtoms();
		int i = 1;
		if (!pOp)
			BOOST_LOG_TRIVIAL(error)<< format("Could not fill the unit cell, because the operation is missing from OpenBabel");
		pOp->Do((OBBase*) &_mol, "strict", NULL, NULL);
		_atomsAfterFillUC = _mol.NumAtoms();

//		BOOST_LOG_TRIVIAL(info)<< format("Asymmetric unit filled with %d atoms.") % _atomsBeforeFillUC;
//		BOOST_LOG_TRIVIAL(info)<< format("Unit cell       filled with %d atoms.") % _atomsAfterFillUC;

		matrix3x3 m = _unitCell->GetCellMatrix();
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				Mm[i][j] = m.Get(i,j);
	}
}
int CifReader::ReadAtoms(std::vector<atom> &atoms,
		std::vector<int> &uniqueAtoms,bool fillUnitCell) {
	ifstream in(_filepath);

	stringstream out;
	OBConversion _conversion(&in, &out);
	OBOp* pOp = OBOp::FindType("fillUC");
	if (_conversion.SetInFormat("CIF")) {
		OBMol _mol;

		_conversion.ReadFile(&_mol, _filepath);
		OBUnitCell* _unitCell = (OBUnitCell*)_mol.GetData(OBGenericDataType::UnitCell);
		_atomsBeforeFillUC = _mol.NumAtoms();
		int i = 0;
		if(fillUnitCell){
			if (!pOp)
				BOOST_LOG_TRIVIAL(error)<< format("Could not fill the unit cell, because the operation is missing from OpenBabel");
			pOp->Do((OBBase*) &_mol, "strict", NULL, NULL);
			_atomsAfterFillUC = _mol.NumAtoms();
		}

		std::list<int> uniqueZ;
		atoms.resize(_atomsAfterFillUC);

		OBAtom* at = nullptr;
		while (at = _mol.GetAtom(++i)) {
//			BOOST_LOG_TRIVIAL(info)<< format("Atom %d ") % (i-1);
			vector3 v = at->GetVector();
			vector3 v1 = _unitCell->CartesianToFractional(v);
			atom* a = &atoms[i-1];

			a->Znum = at->GetAtomicNum();
			a->mass = at->GetAtomicMass();
			a->dw = 0.45*28.0/(double)(2.0*a->Znum);;
			a->q = at->GetFormalCharge();
			a->x = v1.x();
			a->y = v1.y();
			a->z = v1.z();
			a->occ = 1;

			bool thisZexists = (std::find(uniqueZ.begin(), uniqueZ.end(),
					atoms[i-1].Znum)) != uniqueZ.end();
			if (!thisZexists) {
				uniqueZ.push_back(atoms[i-1].Znum);
				uniqueAtoms.push_back(atoms[i-1].Znum);
			}
		}
		return 1;
	}
}
}
