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

#ifndef CONFIG_READER_INTERFACE_H
#define CONFIG_READER_INTERFACE_H

#include "stemtypes_fftw3.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <locale>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/log/trivial.hpp>
using namespace std;
using boost::property_tree::ptree;

namespace QSTEM {

enum class ExperimentType {
	STEM = 3, CBED = 1, TEM = 4, NBED = 2, PTYC = 5
};
enum class SliceThicknessCalculation {
	Auto = 1, SliceThickness = 2, NumberOfSlices = 3
};
enum class StructureFactorType {
	WeickKohl = 1, Rez = 2
};
enum class SaveLevel {
	Everything = 1, Something = 2, Results = 3
};
enum class ResolutionCalculation {
	FILLRES = 1, FILLN = 2, BOXRES = 3, BOXN = 4
};
enum class DisplacementType {
	Einstein = 1, Phonon = 2, None = 3
};

class IPropertyTreeReader {
public:
	virtual void Read(ptree& t) = 0;
protected:
	virtual ~IPropertyTreeReader() {
	}
	;
};

}

#endif
