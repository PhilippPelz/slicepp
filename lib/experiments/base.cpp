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

#include "base.hpp"
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
using boost::format;

namespace slicepp {

BaseExperiment::BaseExperiment(ConfigPtr c, StructureBuilderPtr s, WavePtr w, PotPtr p, DetPtr d, PersistenceManagerPtr pers) :
		IExperiment() {
	_c = c;
	_persist = pers;
	_structureBuilder = s;
	_wave = w;
	_pot = p;
	_det = d;

	DisplayParams();
}

void BaseExperiment::DisplayParams() {
	FILE* fpDir;
	char systStr[64];
	double k2max, temp;
	int i, j;
	static char Date[16], Time[16];
	time_t caltime;
	struct tm* mytime;
	caltime = time(NULL);
	mytime = localtime(&caltime);
	strftime(Date, 12, "%Y:%m:%d", mytime);
	strftime(Time, 9, "%H:%M:%S", mytime);

	BOOST_LOG_TRIVIAL(info)<<
	"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info)<< format("* Running program slice++ (version %s) in %d mode") % VERSION %
	static_cast<int>(_c->ExpType);
	BOOST_LOG_TRIVIAL(info)<< format("* Date: %s, Time: %s") % Date % Time;
	BOOST_LOG_TRIVIAL(info)<< format("* Output file/folder:          %s") % _c->Output->SavePath;
	BOOST_LOG_TRIVIAL(info)<< format("* TDS:                         %d runs") % _c->Model->TDSRuns;
	BOOST_LOG_TRIVIAL(info)<<
	"**************************************************************************************************\n";
}
void BaseExperiment::SetResolution(superCellBoxPtr b) {
	float_tt max_x = b->ax, max_y = b->by, max_z = b->cz, zTotal;
	auto mc = _c->Model;
	auto wc = _c->Wave;
	auto sc = _c->Structure;
	auto dc = _c->Detector;

	//TODO: nx,ny = integer multiple of # unit cells
	//TODO: super cell size = N * unit cell size

	switch (mc->ResCalcType) {
		case ResolutionCalculation::FILLN:
			if (sc->isBoxed) {

			} else {

			}
			mc->n[0] = (mc->n[0] % 2 != 0) ? mc->n[0] + 1 : mc->n[0];
			mc->n[1] = (mc->n[1] % 2 != 0) ? mc->n[1] + 1 : mc->n[1];
			mc->d[0] = (max_x - 0) / mc->n[0];
			mc->d[1] = (max_y - 0) / mc->n[1];
			mc->offset[0] = 0;
			mc->offset[1] = 0;
			break;
		case ResolutionCalculation::FILLRES:
			mc->n[0] = ceil((max_x - 0) / mc->d[0]);
			mc->n[1] = ceil((max_y - 0) / mc->d[1]);
			mc->n[0] = (mc->n[0] % 2 != 0) ? mc->n[0] + 1 : mc->n[0];
			mc->n[1] = (mc->n[1] % 2 != 0) ? mc->n[1] + 1 : mc->n[1];
			mc->d[0] = (max_x - 0) / mc->n[0];
			mc->d[1] = (max_y - 0) / mc->n[1];
			mc->offset[0] = 0;
			mc->offset[1] = 0;
			break;
		case ResolutionCalculation::BOXRES:
			mc->n[0] = sc->box[0] / mc->d[0];
			mc->n[1] = sc->box[1] / mc->d[1];
			mc->n[0] = (mc->n[0] % 2 != 0) ? mc->n[0] + 1 : mc->n[0];
			mc->n[1] = (mc->n[1] % 2 != 0) ? mc->n[1] + 1 : mc->n[1];
			if (mc->CenterSample) {
				mc->offset[0] = sc->box[0] / 2 - (max_x - 0) / 2;
				mc->offset[1] = sc->box[1] / 2 - (max_y - 0) / 2;
			} else {
				mc->offset[0] = 0;
				mc->offset[1] = 0;
			}
			break;
		case ResolutionCalculation::BOXN:
			mc->n[0] = (mc->n[0] % 2 != 0) ? mc->n[0] + 1 : mc->n[0];
			mc->n[1] = (mc->n[1] % 2 != 0) ? mc->n[1] + 1 : mc->n[1];
			mc->d[0] = sc->box[0] / mc->n[0];
			mc->d[1] = sc->box[1] / mc->n[1];
			if (mc->CenterSample) {
				mc->offset[0] = sc->box[0] / 2 - (max_x - 0) / 2;
				mc->offset[1] = sc->box[1] / 2 - (max_y - 0) / 2;
			} else {
				mc->offset[0] = 0;
				mc->offset[1] = 0;
			}
			break;

	}
	if (_c->ExpType != PTYCHO) {
		wc->n[0] = mc->n[0];
		wc->n[1] = mc->n[1];
		dc->n[0] = mc->n[0];
		dc->n[1] = mc->n[1];
	}
}
void BaseExperiment::SetSliceThickness(superCellBoxPtr b) {
	float_tt max_x = b->ax, max_y = b->by, max_z = b->cz, zTotal = b->cz;
	auto mc = _c->Model;

	switch (mc->SliceCalcType) {
		case SliceThicknessCalculation::Auto:
			mc->d[2] = (zTotal / ((int) zTotal)) + 0.01 * (zTotal / ((int) zTotal));
			mc->n[2] = (int) zTotal + 1;
			break;
		case SliceThicknessCalculation::NumberOfSlices:
			mc->d[2] = (zTotal / mc->n[2]) + 0.01 * (zTotal / mc->n[2]);
			break;
		case SliceThicknessCalculation::SliceThickness:
			mc->n[2] = (int) (zTotal / mc->d[2]);
			break;
		default:
			break;
	}
	int atomRadiusSlices = ceil(mc->ratom / mc->d[2]);
}
void BaseExperiment::DisplayProgress(int flag) {
	// static double timer;
	static double timeAvg = 0;
	static double intensityAvg = 0;
	static time_t time0, time1;
	double curTime;
	int jz;

	if (flag < 0) {
		time(&time0);
		// timer = cputim();
		return;
	}
	time(&time1);
	curTime = difftime(time1, time0);
	/*   curTime = cputim()-timer;
	 if (curTime < 0) {
	 printf("timer: %g, curr. time: %g, diff: %g\n",timer,cputim(),curTime);
	 }
	 */
	if (_c->Output->LogLevel > 0) {
		if (_c->Model->UseTDS) {
			timeAvg = ((_runCount) * timeAvg + curTime) / (_runCount + 1);
			intensityAvg = ((_runCount) * intensityAvg + m_intIntensity) / (_runCount + 1);
			BOOST_LOG_TRIVIAL(info)<< format("********************** run %3d ************************") %
			(_runCount + 1);

			std::map<unsigned, float_tt> displacements(_structureBuilder->GetU2());
			std::map<unsigned, float_tt>::iterator disp = displacements.begin(), end = displacements.end();

			BOOST_LOG_TRIVIAL(info)<< format("* <u>: %3d |") % (*disp++).first;
			while (disp != end)
				BOOST_LOG_TRIVIAL(info)<< format(" %8d |") % (*disp++).first;

			BOOST_LOG_TRIVIAL(info)<< " intensity | time(sec) |    chi^2  |";
			BOOST_LOG_TRIVIAL(info)<< "*";

			// ComputeAverageU2();

			disp = displacements.begin();
			while (disp != end)
				BOOST_LOG_TRIVIAL(info)<< format(" %8f |") % (*disp++).second;
			BOOST_LOG_TRIVIAL(info)<< format(" %9f | %9f | %9f |") % m_intIntensity % curTime %
			(_runCount > 0 ? m_chisq[_runCount - 1] : 0);
			BOOST_LOG_TRIVIAL(info)<< format("*");

			BOOST_LOG_TRIVIAL(info)<< format(" %9f | %9f ") % intensityAvg % timeAvg;
		} else {
			BOOST_LOG_TRIVIAL(info) << format("Finished calculations after %.1f s") % curTime;
		}
	} // end of printLevel check.

	time(&time0);
	//  timer = cputim();
}

int BaseExperiment::RunMultislice(af::array t) {
	BOOST_LOG_TRIVIAL(info)<< "Propagating through slices ...";
	af::timer time = af::timer::start();
	for (int i = 0; i < _c->Model->n[2]; i++) {

		af::array slice = t(af::span, af::span, i);
		_wave->Transmit(slice);
		if (_c->Output->SaveWaveAfterTransmit)
			_persist->SaveWaveAfterTransmit(_wave->GetWaveAF(), i);

		_wave->ToFourierSpace();
		if (_c->Output->SaveWaveAfterTransform)
			_persist->SaveWaveAfterTransform(_wave->GetWaveAF(), i);

		_wave->PropagateToNextSlice();
		if (_c->Output->SaveWaveAfterPropagation)
			_persist->SaveWaveAfterPropagation(_wave->GetWaveAF(), i);

		CollectIntensity(i);

		_wave->ToRealSpace();

		if (_c->Output->SaveWaveAfterSlice && i % _c->Output->SaveWaveIterations == 0)
			_persist->SaveWaveAfterSlice(_wave->GetWaveAF(), i);
		PostSliceProcess(i);

		if (_c->Output->LogLevel <= 2) { ///info
			if (i % (int) ceil(_c->Model->n[2] / 10.0) == 0)
				loadbar(i + 1, _c->Model->n[2]);
			auto psi = _wave->GetIntegratedIntensity();
			BOOST_LOG_TRIVIAL(info)<< format("slice %-3d I=%-3.3f") % i % (psi);
		}
	}
	BOOST_LOG_TRIVIAL(info)<< format( "%g ms used for wave propagation (%g us per slice)")
	% (af::timer::stop(time)*1000) %( af::timer::stop(time)*1e6 / _c->Model->n[2]);
	return 0;
}

void BaseExperiment::PostSpecimenProcess() {
	_det->RecordImage(_wave->GetWaveAF());
}

} // end namespace slicepp
