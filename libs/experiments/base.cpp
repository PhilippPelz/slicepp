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
#include "fftw++.hpp"
using boost::format;

namespace QSTEM
{

BaseExperiment::BaseExperiment(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p,const PersistenceManagerPtr& pers)
: IExperiment()
{
	_c = c;
	_persist = pers;
	_structureBuilder = s;
	_wave = w;
	_pot = p;
	m_equalDivs = true;
	m_avgArray = RealVector();
	m_saveLevel = static_cast<unsigned>(c->Output.SaveLevel);

	DisplayParams();
}

void BaseExperiment::DisplayParams()
{
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

	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info) << format("* Running program slice++ (version %s) in %d mode") % VERSION %
			static_cast<int>(_c->ExperimentType);
	BOOST_LOG_TRIVIAL(info) << format("* Date: %s, Time: %s") % Date % Time;
	BOOST_LOG_TRIVIAL(info) << format("* Output file/folder:          %s") %
			_c->Output.savePath.string().c_str();
			(m_equalDivs ? "equal" : "non-equal");
	BOOST_LOG_TRIVIAL(info) << format("* TDS:                         %d runs") % _c->Model.TDSRuns;
	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************\n";
}
void BaseExperiment::SetResolution(superCellBoxPtr b){
	float_tt max_x = b->ax, max_y=b->by, max_z=b->cz, zTotal;
	ModelConfig& mc = _c->Model;

	//TODO: nx,ny = integer multiple of # unit cells
	//TODO: super cell size = N * unit cell size

	switch(mc.ResolutionCalculation) {
	case ResolutionCalculation::FILLN:
		if(_c->Structure.isBoxed){

		} else {

		}
		mc.nx = (mc.nx % 2 != 0) ? mc.nx+1 : mc.nx;
		mc.ny = (mc.ny % 2 != 0) ? mc.ny+1 : mc.ny;
		mc.dx = (max_x - 0)/mc.nx;
		mc.dy = (max_y - 0)/mc.ny;
		_c->Structure.xOffset = 0;
		_c->Structure.yOffset = 0;
		break;
	case ResolutionCalculation::FILLRES:
		mc.nx = ceil((max_x - 0) / mc.dx) ;
		mc.ny = ceil((max_y - 0) / mc.dy) ;
		mc.nx = (mc.nx % 2 != 0) ? mc.nx+1 : mc.nx;
		mc.ny = (mc.ny % 2 != 0) ? mc.ny+1 : mc.ny;
		mc.dx = (max_x - 0)/mc.nx;
		mc.dy = (max_y - 0)/mc.ny;
		_c->Structure.xOffset = 0;
		_c->Structure.yOffset = 0;
		break;
	case ResolutionCalculation::BOXRES:
		mc.nx =  _c->Structure.boxX / mc.dx ;
		mc.ny =  _c->Structure.boxY / mc.dy ;
		mc.nx = (mc.nx % 2 != 0) ? mc.nx+1 : mc.nx;
		mc.ny = (mc.ny % 2 != 0) ? mc.ny+1 : mc.ny;
		if(mc.CenterSample) {
			_c->Structure.xOffset = _c->Structure.boxX/2 - (max_x-0)/2;
			_c->Structure.yOffset = _c->Structure.boxY/2 - (max_y-0)/2;
		} else {
			_c->Structure.xOffset = 0;
			_c->Structure.yOffset = 0;
		}
		break;
	case ResolutionCalculation::BOXN:
		mc.nx = (mc.nx % 2 != 0) ? mc.nx+1 : mc.nx;
		mc.ny = (mc.ny % 2 != 0) ? mc.ny+1 : mc.ny;
		mc.dx =  _c->Structure.boxX / mc.nx ;
		mc.dy =  _c->Structure.boxY / mc.ny ;
		if(mc.CenterSample) {
			_c->Structure.xOffset = _c->Structure.boxX/2 - (max_x-0)/2;
			_c->Structure.yOffset = _c->Structure.boxY/2 - (max_y-0)/2;
		} else {
			_c->Structure.xOffset = 0;
			_c->Structure.yOffset = 0;
		}
		break;

	}
}
void BaseExperiment::SetSliceThickness(superCellBoxPtr b){
	float_tt max_x = b->ax, max_y=b->by, max_z=b->cz, zTotal=b->cz;
	ModelConfig& mc = _c->Model;

	switch (mc.SliceThicknessCalculation) {
	case SliceThicknessCalculation::Auto:
		mc.dz = (zTotal/((int)zTotal))+0.01*(zTotal/((int)zTotal));
		mc.nSlices = (int)zTotal+1;
		break;
	case SliceThicknessCalculation::NumberOfSlices:
		mc.dz = (zTotal/mc.nSlices)+0.01*(zTotal/mc.nSlices);
		break;
	case SliceThicknessCalculation::SliceThickness:
		mc.nSlices = (int)(zTotal / mc.dz);
		break;
	default:
		break;
	}
	int atomRadiusSlices = ceil(_c->Potential.AtomRadiusAngstrom / _c->Model.dz);
	if(_c->Potential.Use3D)
		_c->Model.nSlices += 2 * atomRadiusSlices;
}
void BaseExperiment::DisplayProgress(int flag)
{
	// static double timer;
	static double timeAvg = 0;
	static double intensityAvg = 0;
	static time_t time0, time1;
	double curTime;
	int jz;

	if(flag < 0) {
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
	if(_c->Output.LogLevel > 0) {
		if(_c->Model.UseTDS) {
			timeAvg = ((_runCount)*timeAvg + curTime) / (_runCount + 1);
			intensityAvg = ((_runCount)*intensityAvg + m_intIntensity) / (_runCount + 1);
			BOOST_LOG_TRIVIAL(info) << format("********************** run %3d ************************") %
					(_runCount + 1);

			std::map<unsigned, float_tt> displacements(_structureBuilder->GetU2());
			std::map<unsigned, float_tt>::iterator disp = displacements.begin(), end = displacements.end();

			BOOST_LOG_TRIVIAL(info) << format("* <u>: %3d |") % (*disp++).first;
			while(disp != end)
				BOOST_LOG_TRIVIAL(info) << format(" %8d |") % (*disp++).first;

			BOOST_LOG_TRIVIAL(info) << " intensity | time(sec) |    chi^2  |";
			BOOST_LOG_TRIVIAL(info) << "*";

			// ComputeAverageU2();

			disp = displacements.begin();
			while(disp != end)
				BOOST_LOG_TRIVIAL(info) << format(" %8f |") % (*disp++).second;
			BOOST_LOG_TRIVIAL(info) << format(" %9f | %9f | %9f |") % m_intIntensity % curTime %
					(_runCount > 0 ? m_chisq[_runCount - 1] : 0);
			BOOST_LOG_TRIVIAL(info) << format("*");

			/*
	    // TODO: averaging should be handled on this class, not on lower level crystal class.
	    atom = atomTypes.begin();
	    while (atom!=end) printf(" %8f |",(float)(m_crystal->GetU2avg((*atom++))));
			 */
			BOOST_LOG_TRIVIAL(info) << format(" %9f | %9f ") % intensityAvg % timeAvg;
		} else {
			BOOST_LOG_TRIVIAL(info) << format("Finished calculations after %.1f s") % curTime;
		}
	} // end of printLevel check.

	time(&time0);
	//  timer = cputim();
}

/*
// TODO: this is very broken.  displaced here from Crystal because I want to handle averages outside of the
//     lower level classes.
void CExperimentBase::ComputeAverageU2()
{

  (*z)->second /= u2Count[(*z)->first];
  if (runCount > 0)
    m_u2avg[(*z)] = sqrt(((runCount-1)*(m_u2avg[(*z)]*m_u2avg[(*z)])+u2[(*z)])/runCount);
  else
    m_u2avg[(*z)] = sqrt(m_u2[(*z)]);
}
 */


int BaseExperiment::RunMultislice()
{
	int printFlag = 0;
	int showEverySlice = 1;
	int islice, i, ix, iy, mRepeat;
	float_tt cztot = 0.0;
	float_tt wavlen, sum = 0.0; //,zsum=0.0
	float_tt x, y;
	int absolute_slice;
	char outStr[64];
	double fftScale;
	int nx, ny;

	_wave->GetSizePixels(nx, ny);

	printFlag = (_c->Output.LogLevel > 3);
	fftScale = 1.0 / (nx * ny);

	wavlen = _wave->GetWavelength();

	int nx1, ny1;
	_wave->GetExtents(nx1, ny1);
	m_avgArray.resize(nx1 * ny1);

	/*  calculate the total specimen thickness and echo */
	cztot = 0.0;

	if(printFlag) {
		for(islice = 0; islice < _c->Model.nSlices; islice++) {
			cztot += _c->Model.dz;
		}
		BOOST_LOG_TRIVIAL(info) << format("Specimen thickness: %g Angstroms\n") % cztot;
	}

	BOOST_LOG_TRIVIAL(info) << "Propagating through slices ...";
	for(islice = 0; islice < _c->Model.nSlices; islice++) {

		_wave->Transmit(_pot->GetSlice(islice));
//		Transmit(wave, islice);
		if(_c->Output.SaveWaveAfterTransmit) _persist->SaveWaveAfterTransmit(_wave->GetWaveAF(), islice);

		_wave->ToFourierSpace();

		if(_c->Output.SaveWaveAfterTransform) _persist->SaveWaveAfterTransform(_wave->GetWaveAF(), islice);

		_wave->PropagateToNextSlice();
//		Propagate(wave, islice);

		if(_c->Output.SaveWaveAfterPropagation) _persist->SaveWaveAfterPropagation(_wave->GetWaveAF(), islice);

		CollectIntensity(islice);

		_wave->ToRealSpace();

		if(_c->Output.SaveWaveAfterSlice && islice % _c->Output.SaveWaveIterations == 0) _persist->SaveWaveAfterSlice(_wave->GetWaveAF(),islice);
		PostSliceProcess(islice);

		if(islice % (int)ceil(_c->Model.nSlices / 10.0) == 0)
			loadbar(islice + 1, _c->Model.nSlices);
	} /* end for(islice...) */
	return 0;
} // end of runMulsSTEM


void BaseExperiment::AddDPToAvgArray(const WavePtr& wave)
{
	int nx, ny;
	wave->GetExtents(nx, ny);
	unsigned px = nx * ny;
	// get the pointer to the first data element, and do 1D addressing (it's faster)
	float_tt chisq;

//	const float_tt* dp = wave->GetDPPointer();
// TODO replace this
	for(unsigned i = 0; i < px; i++) {
//		float_tt t = m_avgArray[i] * _runCount + dp[i] / (_runCount + 1);
//		chisq += (m_avgArray[i] - t) * (m_avgArray[i] - t);
//		m_avgArray[i] = t;
		// TODO see what this does and fix it
	}
#pragma omp atomic
	m_chisq[_runCount] += chisq / px;
}

void BaseExperiment::_WriteAvgArray(std::string& fileName,
		std::string& comment,
		std::map<std::string, double>& params,
		std::vector<unsigned>& position)
{
	// params["dx"]=1.0/(m_nx*m_dx);
	// params["dy"]=1.0/(m_ny*m_dy);
	params["Thickness"] = m_thickness;
	// TODO use persist _WriteAvgArray
	//		m_imageIO->WriteImage(m_avgArray, fileName, params, comment, position);
}

void BaseExperiment::ReadAvgArray()
{
	std::vector<unsigned> position;
	// TODO use persist
	//		m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void BaseExperiment::ReadAvgArray(unsigned navg)
{
	std::vector<unsigned> position(1);
	position[0] = navg;
	// TODO use persist
	//		m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void BaseExperiment::ReadAvgArray(unsigned positionx, unsigned positiony)
{
	std::vector<unsigned> position(2);
	position[0] = positionx;
	position[1] = positiony;
	// TODO use persist
	//		m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void BaseExperiment::fft_normalize(WavePtr wave)
{
	ComplexArray2DPtr w = wave->GetWave();
	int nx, ny;
	wave->GetExtents(nx, ny);

	float_tt fftScale = 1.0 / (nx * ny);
	for(unsigned i = 0; i < nx; i++)
		for(unsigned j = 0; j < ny; j++) {
			w[i][j] = complex_tt(w[i][j].real() * fftScale, w[i][j].imag() * fftScale);
		}
}

} // end namespace QSTEM
