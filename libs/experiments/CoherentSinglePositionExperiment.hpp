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

#ifndef EXPERIMENT_CBED_H
#define EXPERIMENT_CBED_H

#include "base.hpp"

namespace QSTEM
{

class CoherentSinglePositionExperiment : public BaseExperiment
{
public:
    CoherentSinglePositionExperiment(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p,const PersistenceManagerPtr& pers);

    void Run();

    virtual void DisplayParams();
    virtual void WriteBeams(unsigned absoluteSlice);
    virtual ~CoherentSinglePositionExperiment(){};

protected:
    virtual void PostSpecimenProcess()=0;
    void PostSliceProcess();
    void CollectIntensity(unsigned absoluteSlice);
    void SaveImages();

    unsigned _nbout = 1;				/* number of recorded beams */
    float_tt **_pendelloesung;
    bool _lbeams;			/* flag indicating whether to record beams */
    float_tt _scanXStart, _scanYStart;     /* The beam position on the sample */
    bool _showProbe;            /* if true, saves a plot of the probe */
    bool m_storeSeries;
};
}
#endif


