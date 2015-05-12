/*
  QSTEM - image simulation for TEM/STEM/CBED
  Copyright (C) 2000-2010  Christoph Koch
  Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CRYSTAL_H
#define CRYSTAL_H

#include "stemtypes_fftw3.hpp"
#include "config_IO/config_reader_factory.hpp"
#include "config_IO/read_qsc.hpp"
#include "structure_readers.hpp"
#include <string>
#include <map>
#include "structure_IO/IStructureBuilder.hpp"
#include <armadillo>
#include <boost/filesystem.hpp>

namespace QSTEM
{

struct zoneAxisOptParams{
	arma::mat M;
	std::vector<int> zone;
	std::vector<int> refZone;

	zoneAxisOptParams(arma::mat M, std::vector<int> zone, std::vector<int> refZone){
		this->M = M;
		this->zone = zone;
		this->refZone = refZone;
	}
};
class CrystalBuilder;
typedef boost::shared_ptr<CrystalBuilder> StructurePtr;

class  CrystalBuilder : public IStructureBuilder{
public:

	CrystalBuilder(boost::filesystem::path p);
	CrystalBuilder(StructureReaderFactory sfac,const ConfigPtr c);
	~CrystalBuilder();

	void ReadFromFile();
	void Init(unsigned run_number);
	void MakeCrystal(bool handleVacancies);
	void TiltBoxed(int ncoord,bool handleVacancies);
	void EinsteinDisplacement(FloatArray2D& u, atom &_atom);
	void PhononDisplacement(FloatArray2D &u,int id,int icx,int icy, int icz,atom &atom,bool printReport);
	void ReplicateUnitCell(int handleVacancies);
	void WriteStructure(unsigned run_number);

	void DisplaceAtoms();
	void CalculateCrystalBoundaries();
	void CalculateTiltAnglesFromZoneAxis();

	void SetFillUnitCell(bool val);
	void SetCellParameters(float_tt ax, float_tt by, float_tt cz);
	void SetNCells(unsigned nx, unsigned ny, unsigned nz);

	float_tt GetCZ(){return m_cz;}
	void GetCellAngles(float_tt &alpha, float_tt &beta, float_tt &gamma);
	void GetCellParameters(float_tt &ax, float_tt &by, float_tt &cz);
	void GetCrystalBoundaries(float_tt &min_x, float_tt &max_x, float_tt &min_y, float_tt &max_y,
			float_tt &min_z, float_tt &max_z);
	inline arma::mat GetCellMatrix(){
		arma::mat M = {m_Mm[0][0],m_Mm[1][0],m_Mm[2][0],
			m_Mm[0][1],m_Mm[1][1],m_Mm[2][1],
			m_Mm[0][2],m_Mm[1][2],m_Mm[2][2]};
		M.reshape(3,3);
			return M;}
	inline unsigned GetNumberOfAtomTypes(){return m_u2.size();}
	inline std::map<unsigned, float_tt> GetU2(){return m_u2;}
	inline float_tt GetU2avg(unsigned znum){return m_u2[znum]/m_u2Count[znum];} //returns Mean Squared displacement
	inline void GetAtom(unsigned idx, atom &_atom){_atom=_atoms[idx];}
	inline unsigned GetNumberOfCellAtoms(){return _baseAtoms.size();}
	inline unsigned GetNumberOfAtoms(){return _atoms.size();}
	inline virtual float_tt GetU2(unsigned znum){return m_u2[znum];}
	inline const std::vector<atom> GetUnitCellAtoms(){return _baseAtoms;}

	virtual superCellBoxPtr Build();
	virtual void DisplayParams();
	virtual std::vector<int> GetUniqueAtoms();
	virtual void SetSliceThickness(ModelConfig& mc);
	virtual void SetResolution(ModelConfig& mc, const PotentialConfig pc);

	static int AtomCompareZnum(const void *atPtr1,const void *atPtr2);
	static int AtomCompareZYX(const void *atPtr1,const void *atPtr2);

	std::vector<atom> _baseAtoms; // The atoms read directly from the input file (no alteration)
	std::vector<atom> _atoms; // The atoms after duplication, tilt, and phonon shaking
	std::vector<int> _uniqueAtoms;
	float_tt _sizeX,_sizeY,_sizeZ;
protected:

	CrystalBuilder();
	boost::filesystem::path m_structureFile;
	bool _fillUnitCell;
	FloatArray2D m_Mm;                     /* metric matrix Mm(ax,by,cz,alpha,beta,gamma).  Used to go from fractional
                                          coordinates to physical cartesian coordinates.  */
	FloatArray2D m_MmInv;                  /* inverse of metric matrix.  Used to go from physical, cartesian coordinates
                                          back to fractional coordinates. */
	float_tt m_ax, m_by, m_cz;           /* lattice parameters */
	float_tt m_cAlpha, m_cBeta, m_cGamma;
	float_tt m_maxX, m_minX;
	float_tt m_maxY, m_minY;
	float_tt m_maxZ, m_minZ;

	bool m_adjustCubeSize;
	float_tt m_offsetX, m_offsetY;
	float_tt m_ctiltx, m_ctilty, m_ctiltz;  /* crystal tilt in mrad */
	unsigned m_nCellX, m_nCellY, m_nCellZ;  /* number of unit cells in x-y-z dir*/
	float_tt m_wobble_temp_scale;
	std::map<unsigned, float_tt> m_u2;  /* (current) rms displacement of atoms */
	std::map<unsigned, unsigned> m_u2Count; /* count of atoms that have been displaced.  Used for computing m_u2avg. */

	boost::filesystem::path m_phononFile;

	StructureReaderPtr _structureReader;
	StructureWriterPtr _structureWriter;
	ConfigPtr _config;
	superCellBoxPtr _superCellBox;

	void CalculateCellDimensions();
	void OffsetCenter(atom &center);
	void Inverse_3x3(FloatArray2D& res, const FloatArray2D& a);
	void RotateVect(float_tt *vectIn,float_tt *vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
	void MatrixProduct(const FloatArray2D& a,int Nxa, int Nya, const FloatArray2D& b,int Nxb, int Nyb, FloatArray2D& c);
	void RotateMatrix(const FloatArray2D& matrixIn, FloatArray2D& matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);

};
}
#endif









