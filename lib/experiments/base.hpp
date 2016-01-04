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
#include "detectors/detector_interface.hpp"
#include "wavefunctions/wave_interface.hpp"
#include "structure_IO/crystal.hpp"
#include "data_IO/PersistenceManager.hpp"
#include "potentials/pot_interface.hpp"

namespace slicepp {
class DLL_EXPORT BaseExperiment: public IExperiment {
public:
	BaseExperiment(ConfigPtr c, StructureBuilderPtr s, WavePtr w, PotPtr p, DetPtr d, PersistenceManagerPtr pers);
	virtual ~BaseExperiment() {
	}
	;
	virtual void DisplayProgress(int flag);
	virtual void DisplayParams();
	virtual void Run()=0;
	virtual void SaveImages()=0;
	void SetResolution(superCellBoxPtr);
	void SetSliceThickness(superCellBoxPtr);

protected:
	// Called in RunMuls after a slice is transmitted/propagated through.
	virtual void PostSliceProcess(unsigned absoluteSlice) {
	}
	;
	virtual void CollectIntensity(unsigned absoluteSlice)=0;
	// Called after the wave has propagated through all slices
	virtual void PostSpecimenProcess();
	// Run the multislice algorithm on the given potential array
	virtual int RunMultislice(af::array t_af);

	unsigned _runCount;  // number of runs to average; runs currently averaged

	ConfigPtr _c;
	// the electron wavefunction
	WavePtr _wave;
	// the sample potential
	PotPtr _pot;
	// the structure builder
	StructureBuilderPtr _structureBuilder;
	// saving stuff
	PersistenceManagerPtr _persist;
	// detector
	DetPtr _det;

	float_tt m_intIntensity; // Integrated intensity from experiment - if too low, your wave array is too small, and the beam is being scattered beyond it.
	RealVector m_chisq;
};

typedef boost::function<
		BaseExperiment*(ConfigPtr c, StructureBuilderPtr s, WavePtr w, PotPtr p, DetPtr d, PersistenceManagerPtr pers)> experimentCreator;
typedef std::map<int, experimentCreator> ExperimentFactory;

}
#endif
