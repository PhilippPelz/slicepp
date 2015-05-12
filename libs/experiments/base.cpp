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

CExperimentBase::CExperimentBase(ConfigPtr c, StructureBuilderPtr s, WavePtr w, PotPtr p, PersistenceManagerPtr pers)
: IExperiment()
{
	_config = c;
	_persist = pers;
	_structureBuilder = s;
	_wave = w;
	_pot = p;
	m_equalDivs = true;
	m_saveLevel = static_cast<unsigned>(c->Output.SaveLevel);
	int atomRadiusSlices = ceil(_config->Potential.AtomRadiusAngstrom / c->Model.sliceThicknessAngstrom);
	if(_config->Potential.Use3D)
		c->Model.nSlices += 2 * atomRadiusSlices;

	DisplayParams();
	_structureBuilder->DisplayParams();
	_pot->DisplayParams();
	_wave->DisplayParams();
	m_avgArray = RealVector();
}

void CExperimentBase::DisplayParams()
{
	FILE* fpDir;
	char systStr[64];
	double k2max, temp;
	int i, j;
	static char Date[16], Time[16];
	time_t caltime;
	struct tm* mytime;
	const double pi = 3.1415926535897;

	/*
    if (wave->printLevel < 1) {
    if ((fpDir = fopen(muls.folder.c_str(),"r"))) {
      fclose(fpDir);
      // printf(" (already exists)\n");
    }
    else {
      sprintf(systStr,"mkdir %s",muls.folder.c_str());
      system(systStr);
      // printf(" (created)\n");
    }
    return;
    }
	 */
	caltime = time(NULL);
	mytime = localtime(&caltime);
	strftime(Date, 12, "%Y:%m:%d", mytime);
	strftime(Time, 9, "%H:%M:%S", mytime);

	BOOST_LOG_TRIVIAL(info) <<
	"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info) << format("* Running program slice++ (version %s) in %d mode") % VERSION %
			static_cast<int>(_config->ExperimentType);
	BOOST_LOG_TRIVIAL(info) << format("* Date: %s, Time: %s") % Date % Time;
	BOOST_LOG_TRIVIAL(info) << format("* Output file/folder:          ./%s/ ") %
			_config->Output.savePath.string().c_str();
			(m_equalDivs ? "equal" : "non-equal");
	BOOST_LOG_TRIVIAL(info) << format("* TDS:                  %d runs)") % _config->Model.TDSRuns;
}

void CExperimentBase::DisplayProgress(int flag)
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
	if(_config->Output.LogLevel > 0) {
		if(_config->Model.UseTDS) {
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
			BOOST_LOG_TRIVIAL(info) << format("**************** finished after %.1f sec ******************") % curTime;
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


void CExperimentBase::InitializePropagators(WavePtr wave)
{
	int nx, ny;

	wave->GetSizePixels(nx, ny);
	m_propxr.resize(nx);
	m_propxi.resize(nx);
	m_propyr.resize(ny);
	m_propyi.resize(ny);

	float_tt scale = _config->Model.sliceThicknessAngstrom * PI;

	//		BOOST_LOG_TRIVIAL(debug) << format("* InitializePropagators") ;

#pragma omp parallel for
	for(int ixa = 0; ixa < nx; ixa++) {
		float_tt t = scale * (wave->GetKX2(ixa) * wave->GetWavelength());
		m_propxr[ixa] = (float_tt)cos(t);
		m_propxi[ixa] = (float_tt)-sin(t);
		//			BOOST_LOG_TRIVIAL(debug) << format("* m_propx[%d] = (%g,%g)") % ixa % m_propxr[ixa] %
		// m_propxi[ixa];
	}
#pragma omp parallel for
	for(int iya = 0; iya < ny; iya++) {
		float_tt t = scale * (wave->GetKY2(iya) * wave->GetWavelength());
		m_propyr[iya] = (float_tt)cos(t);
		m_propyi[iya] = (float_tt)-sin(t);
	}
}

int CExperimentBase::RunMultislice(WavePtr wave)
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

	wave->GetSizePixels(nx, ny);

	printFlag = (_config->Output.LogLevel > 3);
	fftScale = 1.0 / (nx * ny);

	wavlen = wave->GetWavelength();

	int nx1, ny1;
	wave->GetExtents(nx1, ny1);
	m_avgArray.resize(nx1 * ny1);

	/*  calculate the total specimen thickness and echo */
	cztot = 0.0;

	if(printFlag) {
		for(islice = 0; islice < _config->Model.nSlices; islice++) {
			cztot += _config->Model.sliceThicknessAngstrom;
		}
		BOOST_LOG_TRIVIAL(info) << format("Specimen thickness: %g Angstroms\n") % cztot;
	}

	_wave->InitializePropagators();
//	InitializePropagators(wave);
	BOOST_LOG_TRIVIAL(info) << "Propagating through slices ...";
	for(islice = 0; islice < _config->Model.nSlices; islice++) {
		absolute_slice = (m_totalSliceCount + islice);

		_wave->Transmit(_pot->GetSlice(islice));
//		Transmit(wave, islice);

		if(_config->Output.SaveWaveAfterTransmit) _persist->SaveWaveAfterTransmit(wave->GetWave(), islice);

		_wave->ToFourierSpace();

		if(_config->Output.SaveWaveAfterTransform) _persist->SaveWaveAfterTransform(wave->GetWave(), islice);

		_wave->PropagateToNextSlice();
//		Propagate(wave, islice);

		if(_config->Output.SaveWaveAfterPropagation) _persist->SaveWaveAfterPropagation(wave->GetWave(), islice);

		CollectIntensity(absolute_slice);

		wave->ToRealSpace();

		if(_config->Output.SaveWaveAfterSlice) _persist->SaveWaveAfterSlice(_wave->GetWave(),islice);
		PostSliceProcess(absolute_slice);

		if(islice % (int)ceil(_config->Model.nSlices / 10.0) == 0)
			loadbar(islice + 1, _config->Model.nSlices);
	} /* end for(islice...) */
	return 0;
} // end of runMulsSTEM

void CExperimentBase::Propagate(WavePtr wave, float_tt dz)
{
	float_tt wr, wi, tr, ti;
	float_tt scale, t;
	float_tt dzs = 0;

	float_tt dx, dy;
	int nx, ny;

	wave->GetResolution(dx, dy);
	wave->GetSizePixels(nx, ny);

	ComplexArray2DPtr w = wave->GetWave();

#pragma omp parallel for private(wr, wi, tr, ti)
	for(int i = 0; i < nx; i++)
		for(int j = 0; i < ny; i++) {
			try {
				if((wave->GetKX2(i) + wave->GetKY2(j)) < wave->GetK2Max()) {
					wr = w[i][j].real();
					wi = w[i][j].imag();
					tr = wr * m_propyr[j] - wi * m_propyi[j];
					ti = wr * m_propyi[j] + wi * m_propyr[j];
					w[i][j] = complex_tt(tr * m_propxr[i] - ti * m_propxi[i], tr * m_propxi[i] + ti * m_propxr[i]);
				} else {
					w[i][j] = 0.0F;
				}
			} catch(const std::exception& e) {
				std::cerr << e.what();
			}
		} /* end for(ix..) */
} /* end propagate */

/*------------------------ transmit() ------------------------*/
/*
transmit the wavefunction thru one layer
(simply multiply wave by transmission function)

waver,i[ix][iy]  = real and imaginary parts of wavefunction
transr,i[ix][iy] = real and imag parts of transmission functions

nx, ny = size of array

on entrance waver,i and transr,i are in real space

only waver,i will be changed by this routine
 */
void CExperimentBase::Transmit(WavePtr wave, unsigned sliceIdx)
{
	double wr, wi, tr, ti;
	int nx, ny;
	ComplexArray2DPtr w = wave->GetWave();
	wave->GetSizePixels(nx, ny);

	/*  trans += posx; */
	for(unsigned ix = 0; ix < nx; ix++) {
		for(unsigned iy = 0; iy < ny; iy++) {
			complex_tt t = _pot->GetSlicePixel(sliceIdx, ix , iy );
			wr = w[ix][iy].real();
			wi = w[ix][iy].imag();
			tr = t.real();
			ti = t.imag();
			w[ix][iy] *= t;
			BOOST_LOG_TRIVIAL(trace) << boost::format("w=(%g,%g) t=(%2.3f,%2.3f) w*t=(%g,%g)") % wr % wi % tr % ti %
					w[ix][iy].real() % w[ix][iy].imag();
		} /* end for(iy.. ix .) */
	}
} /* end transmit() */

void CExperimentBase::AddDPToAvgArray(const WavePtr& wave)
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

void CExperimentBase::_WriteAvgArray(std::string& fileName,
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

void CExperimentBase::ReadAvgArray()
{
	std::vector<unsigned> position;
	// TODO use persist
	//		m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::ReadAvgArray(unsigned navg)
{
	std::vector<unsigned> position(1);
	position[0] = navg;
	// TODO use persist
	//		m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::ReadAvgArray(unsigned positionx, unsigned positiony)
{
	std::vector<unsigned> position(2);
	position[0] = positionx;
	position[1] = positiony;
	// TODO use persist
	//		m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::fft_normalize(WavePtr wave)
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
