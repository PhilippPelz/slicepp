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

#include "Ptychograph.hpp"
#include "random.hpp"
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include "math.h"
using boost::format;

namespace QSTEM
{

Ptychograph::Ptychograph(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p, const DetPtr& d, const PersistenceManagerPtr& pers):BaseExperiment(c,s,w,p,d,pers)
{
	m_mode=ExperimentType::PTYC;
	_scan = ScanPtr(new Scan(c));
	_lbeams = false;
}

void Ptychograph::DisplayParams()
{
}

void Ptychograph::Run()
{
	int ix,iy,i,pCount,result;
	double timer,timerTot ;
	float_tt t=0;
	FloatArray2D avgPendelloesung(boost::extents[_nbout][_c->Model.nSlices]);
	std::map<std::string, double> params;
	std::vector<unsigned> position(1);         // Used to indicate the number of averages
	std::vector<std::pair<int, int>> scanPositions;

	m_chisq.resize(_c->Model.TDSRuns);
	timerTot = 0; /* cputim();*/
	DisplayProgress(-1);

	auto box = _structureBuilder->Build();

	SetResolution(box);
	SetSliceThickness(box);

	_wave->FormProbe();
	scanPositions = _scan->GetScanPositions();
//	_persist->InitStorage();
	_persist->ResizeStorage(_c->Wave.nx, _c->Wave.ny);
	_wave->InitializePropagators();
	_wave->DisplayParams();
	_pot->DisplayParams();
	for (_runCount = 0;_runCount < _c->Model.TDSRuns;_runCount++) {
		auto box = _structureBuilder->DisplaceAtoms();
		//PotPtr p  _pot ;
		//p ->MakeSlices(box);

		_pot->MakeSlices(box);
		if (_c->Output.saveProbe) _persist->SaveProbe(_wave->GetProbe());
		BOOST_LOG_TRIVIAL(info) << format("Calculating for %.1f position(s)") % scanPositions.size();
		for(std::vector<std::pair<int, int>>::size_type i = 0; i != scanPositions.size(); i++) {
			_wave->SetPositionOffset(scanPositions[i].first, scanPositions[i].second);

			int xp = scanPositions[i].first, yp = scanPositions[i].second;
			//  _pot = p->GetSubPotential(xp,yp,size)
			// auto e = CoherentSinglePositionExperiment(p,...);
			// e.Run();
			BOOST_LOG_TRIVIAL(info) << format("\n==== Position %.lf (%.lf, %.lf) ====") % (i + 1) % xp % yp;
			//RunMultislice(p->GetSubPotential(xp,yp,size));
			RunMultislice(_pot->GetSubPotential(xp, yp, _c->Wave.nx, _c->Wave.ny));
			if (_runCount == 0) {
				if (_lbeams) {
					for (iy=0;iy<_c->Model.nSlices ;iy++) {
						for (ix=0;ix<_nbout;ix++) {
							avgPendelloesung[ix][iy] = _pendelloesung[ix][iy];
						}
					}
				}
			}
			else {
				m_storeSeries = 1;
				if (m_saveLevel == 0)	m_storeSeries = 0;
				else if (_runCount % m_saveLevel != 0) m_storeSeries = 0;
				if (m_storeSeries)
				{
					params["1/Wavelength"] = 1.0/_wave->GetWavelength();
					WriteAvgArray(_runCount+1, "Averaged Diffraction pattern, unit: 1/A", params);
				}
			}
			PostSpecimenProcess();
			_persist->StoreToDiscMP(i + 1, xp, yp);
			_wave->ResetProbe();
		}
	DisplayProgress(1);
	BOOST_LOG_TRIVIAL(info) << "Saving to disc...";
	BOOST_LOG_TRIVIAL(info) << "Finished saving...";

	}

}


void Ptychograph::WriteBeams(unsigned int absoluteSlice)
{

}

void Ptychograph::CollectIntensity(unsigned absoluteSlice)
{
	WriteBeams(absoluteSlice);
}

void Ptychograph::PostSliceProcess( )
{
}

void Ptychograph::SaveImages()
{

}
} // end namespace QSTEM



