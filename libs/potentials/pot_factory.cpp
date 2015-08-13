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


#include "pot_factory.hpp"
#include "CUDA2DPotential.hpp"
#include "pot_2d.hpp"
#include "pot_2d_fft.hpp"
#include "pot_3d.hpp"
#include "pot_3d_fft.hpp"

namespace QSTEM
{

CPotFactory::CPotFactory()
{
//	Register("2D",    &C2DPotential::Create);
//	Register("2DFFT", &C2DFFTPotential::Create);
//	Register("3D",    &C3DPotential::Create);
//	Register("3DFFT", &C3DFFTPotential::Create);
}

void CPotFactory::Register(const std::string &potName, CreatePotentialFn pfnCreate)
{
	m_FactoryMap[potName] = pfnCreate;
}

PotPtr CPotFactory::GetPotential(const std::string &name,const ConfigPtr c)
{
//	std::cout << name << std::endl;
	FactoryMap::iterator it = m_FactoryMap.find(name);
	if( it != m_FactoryMap.end() )
		return it->second(c);
	return PotPtr();
}

PotPtr CPotFactory::GetPotential(bool _3D, bool fft,const ConfigPtr c)
{
	std::stringstream str;
	str << (_3D ? "3D" : "2D");
	str << (fft ? "FFT" : "");
	str << std::ends;
	return GetPotential(str.str(),c);
}


PotPtr CPotFactory::GetPotential(const ConfigPtr c)
{
  bool _3D = c->Potential.Use3D, fft = c->Potential.UseFFT;
  PotPtr pot=GetPotential(_3D, fft,c);

  return pot;
}

} // end namespace QSTEM
