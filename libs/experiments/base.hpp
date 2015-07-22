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

#ifndef EXPERIMENT_BASE_H
#define EXPERIMENT_BASE_H

#include "experiment_interface.hpp"
#include "wavefunctions/wave_interface.hpp"
#include "structure_IO/crystal.hpp"
#include "data_IO/PersistenceManager.hpp"
#include "potentials/pot_interface.hpp"
#include "fftw++.hpp"

namespace QSTEM
{
class DLL_EXPORT BaseExperiment : public IExperiment
{
public:
    BaseExperiment(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p,const PersistenceManagerPtr& pers);
    virtual ~BaseExperiment(){};
    virtual void DisplayProgress(int flag);
    virtual void DisplayParams();
    virtual void Run()=0;
    virtual void SaveImages()=0;
    void SetResolution(superCellBoxPtr);
    void SetSliceThickness(superCellBoxPtr);

protected:
    virtual void PostSliceProcess(unsigned absoluteSlice) {}; // Called in RunMuls after a slice is transmitted/propagated through.  Override as desired.
    virtual void CollectIntensity(unsigned absoluteSlice)=0;
    virtual int RunMultislice();
    virtual void Transmit(WavePtr wave, unsigned sliceIdx);
    virtual void Propagate(WavePtr wave, float_tt dz);
    virtual void AddDPToAvgArray(const WavePtr &wave);

    void ReadAvgArray();
    void ReadAvgArray(unsigned navg);
    void ReadAvgArray(unsigned posX, unsigned posY);

    void InitializePropagators(WavePtr wave);

    void _WriteAvgArray(std::string &fileName, std::string &comment,
                        std::map<std::string, double> &params,
                        std::vector<unsigned> &position);

    inline void WriteAvgArray(std::string comment="Average Array",
                              std::map<std::string, double>params = std::map<std::string, double>())
    {
    }
    inline void WriteAvgArray(unsigned navg, std::string comment="Average Array",
                              std::map<std::string, double>params = std::map<std::string, double>())
    {
    }
    inline void WriteAvgArray(unsigned posX, unsigned posY, std::string comment="Average Array",
                              std::map<std::string, double>params = std::map<std::string, double>())
    {
    }

    void fft_normalize(WavePtr wave);

    bool m_tds;
    unsigned m_avgRuns, _runCount;  // number of runs to average; runs currently averaged
    unsigned m_printLevel;

    ConfigPtr _c;
    // the electron wave
    WavePtr _wave;
    // the sample potential
    PotPtr _pot;
    // the structure builder
    StructureBuilderPtr _structureBuilder;
    // saving stuff
    PersistenceManagerPtr _persist;

    float_tt m_intIntensity;  // Integrated intensity from experiment - if too low,
    // your wave array is too small, and the beam is being scattered beyond it.

    unsigned m_cellDiv;		// The number of sub-slabs that the supercell is divided into
    bool m_equalDivs;			// Whether or not all sub-slabs are the same size
    unsigned m_outputInterval;  // The number of slices between saving intermediate output files
    unsigned m_totalSliceCount; // The total number of slices that we've progressed through (all sub-slabs included)
    float_tt m_thickness;       // The total thickness of the sample at the current slice
    float_tt m_dz;
    RealVector m_chisq;
    QSTEM::ExperimentType m_mode;      // String representing the multislice mode (e.g. TEM, STEM, etc.)
    RealVector m_avgArray;   // The averaged diffraction pattern (m_avgCount says how many are averaged currently)
               /* integer offset for positioning probe within potential array */

    std::vector<float_tt> m_propxr, m_propxi, m_propyr, m_propyi;
};

typedef boost::function<BaseExperiment*(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p,const PersistenceManagerPtr& pers)> experimentCreator;
typedef std::map<ExperimentType,experimentCreator> ExperimentFactory;

}
#endif
