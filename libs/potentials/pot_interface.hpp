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

#include "stemtypes_fftw3.hpp"
#include "config_IO/config_interface.hpp"
#include "structure_IO/crystal.hpp"
#include "config_IO/configs.hpp"
#include "data_IO/PersistenceManager.hpp"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/function.hpp>

using boost::format;

#ifndef POTENTIAL_INTERFACE_H
#define POTENTIAL_INTERFACE_H

namespace QSTEM
{

class IPotential;
typedef boost::shared_ptr<IPotential> PotPtr;
typedef boost::function<IPotential*(const ConfigPtr& c,const PersistenceManagerPtr &persist) > potentialCreator;
typedef std::map<string,potentialCreator> PotentialFactory;
typedef PotPtr (*CreatePotentialFn)(ConfigPtr);


class IPotential
{
public:
	IPotential(const ConfigPtr& c,const PersistenceManagerPtr& persist) : _c(c), _persist(persist){} ;
	virtual void DisplayParams(){};
	virtual void MakeSlices(superCellBoxPtr info)=0;
	virtual void ReadPotential(std::string &fileName, unsigned subSlabIdx){};

	virtual void GetSizePixels(unsigned &nx, unsigned &ny) const =0;
	virtual void WriteSlice(unsigned idx, std::string prefix)=0;
	virtual void WriteProjectedPotential()=0;
//	virtual ComplexArray2DView GetSlice(unsigned idx)=0;
	virtual af::array GetSlice(unsigned idx) =0;
	virtual complex_tt GetSlicePixel(unsigned iz, unsigned ix, unsigned iy)=0;
	virtual af::array GetSubPotential(int startx, int starty, int nx, int ny){};
	virtual af::array GetPotential(){};
	virtual ~IPotential(){};
protected:
	IPotential(){};
	ConfigPtr _c;
	PersistenceManagerPtr _persist;
};

}

#endif
