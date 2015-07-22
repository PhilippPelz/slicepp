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

#ifndef PLANE_WAVE_H
#define PLANE_WAVE_H

#include "wave_base.hpp"

namespace QSTEM
{

class DLL_EXPORT CPlaneWave : public CBaseWave
{
public:
  CPlaneWave(const ConfigPtr& c,const PersistenceManagerPtr& p);
  CPlaneWave(const CPlaneWave& other);
  virtual void FormProbe();
  virtual void FormProbe(int &scanx, int &scany);
  void TiltBeam(bool tiltBack=false);
  void TiltBack();
  virtual void DisplayParams();

  WavePtr Clone();

  // ReadImage is for TEM mode
  void ReadImage();
  void WriteImage();
protected:
  RealVector m_image;               /* Real-space image output */
private:
  friend class CWaveFactory;
//  static WavePtr Create(const ConfigPtr reader){  return WavePtr(new CPlaneWave(reader));  }
};

typedef boost::shared_ptr<CPlaneWave> PlaneWavePtr;

}

#endif
