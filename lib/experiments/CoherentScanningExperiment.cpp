/*
 * Ptychograph.cpp
 *
 *  Created on: Jul 20, 2015
 *      Author: wenxuan
 */

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

#include "CoherentScanningExperiment.hpp"
#include "afhelpers.hpp"

#include "random.hpp"
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include "math.h"
using boost::format;

namespace slicepp {

CoherentScanningExperiment::CoherentScanningExperiment(const ConfigPtr& c, const StructureBuilderPtr& s, const WavePtr& w, const PotPtr& p,
		const DetPtr& d, const PersistenceManagerPtr& pers) :
		BaseExperiment(c, s, w, p, d, pers) {
}

void CoherentScanningExperiment::DisplayParams() {
}

void CoherentScanningExperiment::Run() {
	DisplayProgress(-1);

	auto time = af::timer::start();
	auto box = _structureBuilder->Build();
	SetResolution(box);
	SetSliceThickness(box);
	auto elapsed = af::timer::stop(time) * 1000;

	BOOST_LOG_TRIVIAL(info)<< format( "%g msec used for building structure.") % elapsed;

	time = af::timer::start();
	_persist->ResizeStorage(_c->Wave->n[0], _c->Wave->n[1]);
	_wave->FormProbe();
	_wave->InitializePropagators();
	_wave->DisplayParams();
	_pot->DisplayParams();
	auto scanPositions = _c->Scan->positions;
	elapsed = af::timer::stop(time) * 1000;

	auto image = af::array(_c->Wave->n[1], _c->Wave->n[0],f32);

	BOOST_LOG_TRIVIAL(info)<< format( "%g msec used for initialization.") % elapsed;

	for (_runCount = 0; _runCount < _c->Model->TDSRuns; _runCount++) {
		auto box = _structureBuilder->DisplaceAtoms();
		_pot->MakeSlices(box);
		if (_c->Output->SaveProbe)
			_persist->SaveProbe(_wave->GetProbe());
		BOOST_LOG_TRIVIAL(info)<< format("Calculating for %.1f position(s)") % scanPositions.size();
		for (auto i = 0; i != scanPositions.size(); i++) {
			int xp = scanPositions[i].x, yp = scanPositions[i].y;
			BOOST_LOG_TRIVIAL(info)<< format("==== Position %.lf (%.lf, %.lf) ====") % (i + 1) % xp % yp;
			RunMultislice(_pot->GetSubPotential(xp, yp, _c->Wave->n[0], _c->Wave->n[1]));

			PostSpecimenProcess();

			_persist->StoreToDiscMP(i + 1, xp, yp);
			_wave->ResetProbe();
		}

		_persist->SavePositions(scanPositions);
		DisplayProgress(1);
		BOOST_LOG_TRIVIAL(info)<< "Finished saving...";
	}
}
void CoherentScanningExperiment::PostSpecimenProcess() {
	_wave->ToFourierSpace();
	auto dp = fftShift(_wave->GetWaveAF());
	_det->RecordImage(dp);
}

void CoherentScanningExperiment::WriteBeams(unsigned int absoluteSlice) {

}

void CoherentScanningExperiment::CollectIntensity(unsigned absoluteSlice) {
	WriteBeams(absoluteSlice);
}

void CoherentScanningExperiment::PostSliceProcess() {
}

void CoherentScanningExperiment::SaveImages() {

}
} // end namespace slicepp

