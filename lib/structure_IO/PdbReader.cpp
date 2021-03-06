/*
 * PdbReader.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: philipp
 */

#include "PdbReader.hpp"
#include <openbabel/math/vector3.h>
#include <openbabel/op.h>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <math.h>
using namespace std;
using namespace OpenBabel;
using boost::format;

namespace slicepp {
PdbReader::PdbReader(const boost::filesystem::path &file) : _atomsAfterFillUC(0), _atomsBeforeFillUC(0), IStructureReader(file){
}

PdbReader::~PdbReader() {
	// TODO Auto-generated destructor stub
}

int PdbReader::ReadCellParams(FloatArray2D& Mm) {
	ifstream in(_path.string());

	stringstream out;
	OBConversion _conversion(&in, &out);
	OBOp* pOp = OBOp::FindType("fillUC");

	if (_conversion.SetInFormat("Pdb")) {
		OBMol _mol;
		string tmp = _path.string();
		tmp.erase(std::remove(tmp.begin(), tmp.end(), '\n'), tmp.end());
		if(_conversion.ReadFile(&_mol, tmp)){
			OBUnitCell* _unitCell = (OBUnitCell*)_mol.GetData(OBGenericDataType::UnitCell);
			_atomsBeforeFillUC = _mol.NumAtoms();
			int i = 1;
			//			_unitCell->FillUnitCell(&_mol);
			if (!pOp)
				BOOST_LOG_TRIVIAL(error)<< format("Could not fill the unit cell, because the operation is missing from OpenBabel");
			pOp->Do((OBBase*) &_mol, "strict", NULL, NULL);
			_atomsAfterFillUC = _mol.NumAtoms();
			matrix3x3 m = _unitCell->GetCellMatrix();
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
					if( fabs(m.Get(i,j)) < 1e-6 )
						Mm[i][j] = 0;
					else
						Mm[i][j] = m.Get(i,j);
		}
		else{
			BOOST_LOG_TRIVIAL(error)<< format("Could not read the Pdb file %s.") % _path.string();
		}
	}
}
int PdbReader::ReadAtoms(std::vector<atom> &atoms,
		std::vector<int> &uniqueAtoms,bool fillUnitCell) {
	ifstream in(_path.string());

	stringstream out;
	OBConversion _conversion(&in, &out);
	OBOp* pOp = OBOp::FindType("fillUC");
	if (_conversion.SetInFormat("PDB")) {
		OBMol _mol;

		string tmp = _path.string();
		tmp.erase(std::remove(tmp.begin(), tmp.end(), '\n'), tmp.end());
		if(_conversion.ReadFile(&_mol, tmp)){
			int i = 0;
			_atomsAfterFillUC = _mol.NumAtoms();
			std::list<int> uniqueZ;
			atoms.resize(_atomsAfterFillUC);

			OBAtom* at = nullptr;
			while (at = _mol.GetAtom(++i)) {
				//			BOOST_LOG_TRIVIAL(info)<< format("Atom %d ") % (i-1);
				vector3 v = at->GetVector();
//				vector3 v1 = _unitCell->CartesianToFractional(v);
				atom* a = &atoms[i-1];

				a->Znum = at->GetAtomicNum();
				a->mass = at->GetAtomicMass();
				a->dw = 0.45*28.0/(float_tt)(2.0*a->Znum);;
				a->q = at->GetFormalCharge();
				a->occ = 1;
				a->r = {(float_tt)v.x(),(float_tt)v.y(),(float_tt)v.z()};

				BOOST_LOG_TRIVIAL(trace)<< format("Base atom %d: (%3.3f,%3.3f,%3.3f)") % i % v.x()%v.y()%v.z();

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
}
