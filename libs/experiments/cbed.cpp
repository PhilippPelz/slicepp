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

#include "cbed.hpp"
#include "random.hpp"
#include "memory_fftw3.hpp"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
using boost::format;

namespace QSTEM
{

CExperimentCBED::CExperimentCBED(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p,const PersistenceManagerPtr& pers) : CExperimentBase(c,s,w,p,pers)
{
	m_mode=ExperimentType::CBED;
	_lbeams = false;
	_config = c;

}

void CExperimentCBED::DisplayParams()
{
}

void CExperimentCBED::Run()
{
	int ix,iy,i,pCount,result;
	FILE *avgFp,*fp,*fpPos=0;
	double timer,timerTot;
	double probeCenterX,probeCenterY,probeOffsetX,probeOffsetY;
	float_tt t=0;
	float_tt **avgPendelloesung = NULL;
	int nx, ny;
	_wave->GetSizePixels(nx, ny);
	std::map<std::string, double> params;
	std::vector<unsigned> position(1);         // Used to indicate the number of averages

	m_chisq.resize(_config->Model.TDSRuns);

	if (_lbeams) {
		_pendelloesung = NULL;
		if (avgPendelloesung == NULL) {
			avgPendelloesung = float2D(_nbout,	_config->Model.nSlices*_config->Potential.NSubSlabs,"pendelloesung");
		}
	}
	probeCenterX = _scanXStart;
	probeCenterY = _scanYStart;
	timerTot = 0; /* cputim();*/
	DisplayProgress(-1);
	_wave->FormProbe();

	for (_runCount = 0;_runCount < _config->Model.TDSRuns;_runCount++) {

		auto box = _structureBuilder->Build();
		_pot->MakeSlices(box);

		m_totalSliceCount = 0;
		pCount = 0;

		if (_config->Output.saveProbe)
			_persist->SaveProbe(_wave->GetWave());

		if (_lbeams){
			if (_pendelloesung != NULL)
				free(_pendelloesung[0]);
			free(avgPendelloesung[0]);
			_pendelloesung = NULL;
			avgPendelloesung = float2D(_nbout,_config->Model.nSlices*_config->Potential.NSubSlabs,"pendelloesung");
		}

		RunMultislice(_wave);

		if (_runCount == 0) {
			/* move the averaged (raw data) file to the target directory as well */
			// TODO: make sure that DP average gets created properly
			//sprintf(avgName,"%s/diffAvg_%d.img",m_folder.c_str(),_runCount+1);
			//sprintf(systStr,"mv %s/diff.img %s",m_folder.c_str(),avgName);
			//system(systStr);
			if (_lbeams) {
				for (iy=0;iy<_config->Model.nSlices*_config->Potential.NSubSlabs;iy++) {
					for (ix=0;ix<_nbout;ix++) {
						avgPendelloesung[ix][iy] = _pendelloesung[ix][iy];
					}
				}
			}
		} // of if m_avgCount == 0 ...
		else {

			m_storeSeries = 1;
			if (m_saveLevel == 0)	m_storeSeries = 0;
			else if (_runCount % m_saveLevel != 0) m_storeSeries = 0;

			if (m_storeSeries)
			{
				params["1/Wavelength"] = 1.0/_wave->GetWavelength();
				WriteAvgArray(_runCount+1, "Averaged Diffraction pattern, unit: 1/A", params);
			}

			/*************************************************************/

			/***********************************************************
			 * Average over the pendelloesung plot as well
			 */
			if (_lbeams) {
				for (iy=0;iy<_config->Model.nSlices*_config->Potential.NSubSlabs;iy++) {
					for (ix=0;ix<_nbout;ix++) {
						avgPendelloesung[ix][iy] =((float_tt)_runCount*avgPendelloesung[ix][iy]+_pendelloesung[ix][iy])/(float_tt)(_runCount+1);
					}
				}
			}
		} /* else ... if avgCount was greater than 0 */

		if (_lbeams) {
			/**************************************************************
			 * The diffraction spot intensities of the selected
			 * diffraction spots are now stored in the 2 dimensional array
			 * m_pendelloesung[beam][slice].
			 * We can write the array to a file and display it, just for
			 * demonstration purposes
			 *************************************************************/
			char systStr[255];
			sprintf(systStr,"%s/pendelloesung.dat",m_outputLocation.c_str());
			if ((fp=fopen(systStr,"w")) !=NULL) {
				BOOST_LOG_TRIVIAL(info) << "Writing Pendelloesung data";
				for (iy=0;iy<_config->Model.nSlices*_config->Potential.NSubSlabs;iy++) {
					/* write the thicknes in the first column of the file */
					fprintf(fp,"%g",iy*_config->Model.sliceThicknessAngstrom);//((float)(m_potential->GetNSlices()*_config->Potential.NSubSlabs)));
					/* write the beam intensities in the following columns */
					for (ix=0;ix<_nbout;ix++) {
						fprintf(fp,"\t%g",avgPendelloesung[ix][iy]);
					}
					fprintf(fp,"\n");
				}
				fclose(fp);
			}
			else {
				BOOST_LOG_TRIVIAL(error) << "Could not open file for pendelloesung plot";
			}
		} /* end of if lbeams ... */
		DisplayProgress(1);
	} /* end of for m_avgCount=0.. */
	_persist->StoreToDisc();
}

void CExperimentCBED::CollectIntensity(unsigned absoluteSlice)
{
	WriteBeams(absoluteSlice);
}

void CExperimentCBED::SaveImages()
{
	// This is called at the designated save interval, and one final time at the end of each run.
	//   TODO: define the images that should be saved here.
}

void CExperimentCBED::WriteBeams(unsigned int absoluteSlice)
{
	// TODO: this needs to be reconsidered in terms of not using static variables

	/*
  if ((fp1 == NULL) || (fpAmpl == NULL) || (fpPhase == NULL)) {
    scale = 1.0F / ( ((float_tt)m_nx) * ((float_tt)m_ny) );
    hbeam = (*muls).hbeams;
    kbeam = (*muls).kbeams;
    if ((hbeam.empty()) || (kbeam.empty())) {
      printf("ERROR: hbeam or kbeam == NULL!\n");
      exit(0);
    }

    sprintf(fileAmpl,"%s/beams_amp.dat",(*muls).folder.c_str());
    sprintf(filePhase,"%s/beams_phase.dat",(*muls).folder.c_str());
    sprintf(fileBeam,"%s/beams_all.dat",(*muls).folder.c_str());
    fp1 = fopen(fileBeam, "w" );
    fpAmpl = fopen( fileAmpl, "w" );
    fpPhase = fopen( filePhase, "w" );
    if(fp1==NULL) {
      printf("can't open file %s\n", fileBeam);
      exit(0);
    }
    if(fpAmpl==NULL) {
      printf("can't open amplitude file %s\n",fileAmpl);
      exit(0);
    }
    if(fpPhase==NULL) {
      printf("can't open phase file %s\n", filePhase);
      exit(0);
    }
    fprintf(fp1, " (h,k) = ");
    for(ib=0; ib<(*muls).nbout; ib++) {
      fprintf(fp1," (%d,%d)", muls->hbeams[ib],  muls->kbeams[ib]);
    }
    fprintf( fp1, "\n" );
    fprintf( fp1, "nslice, (real,imag) (real,imag) ...\n\n");
    for( ib=0; ib<muls->nbout; ib++)
      {
        // printf("beam: %d [%d,%d]",ib,hbeam[ib],kbeam[ib]);			
        if(hbeam[ib] < 0 ) hbeam[ib] = muls->nx + hbeam[ib];
        if(kbeam[ib] < 0 ) kbeam[ib] = muls->ny + kbeam[ib];
        if(hbeam[ib] < 0 ) hbeam[ib] = 0;
        if(kbeam[ib] < 0 ) kbeam[ib] = 0;
        if(hbeam[ib] > muls->nx-1 ) hbeam[ib] = muls->nx-1;
        if(kbeam[ib] > muls->ny-1 ) kbeam[ib] = muls->ny-1;
        // printf(" => [%d,%d] %d %d\n",hbeam[ib],kbeam[ib],muls->nx,muls->ny);			
      }

    // setup of beam files, include the t=0 information 
    fprintf( fpAmpl, "%g",0.0);
    fprintf( fpPhase, "%g",0.0);
    for( ib=0; ib<muls->nbout; ib++) {
      ampl = 0.0;
      if ((hbeam[ib] == 0) && (kbeam[ib]==0))
        ampl = 1.0;
      fprintf(fpAmpl,"\t%g",ampl);
      fprintf(fpPhase,"\t%g",0.0);
    }
    fprintf( fpAmpl, "\n");
    fprintf( fpPhase, "\n");
  } // end of if fp1 == NULL ... i.e. setup 


  zsum += (*muls).cz[ilayer];

  fprintf( fp1, "%g", zsum);
  fprintf( fpAmpl, "%g",zsum);
  fprintf( fpPhase, "%g",zsum);
  for( ib=0; ib<(*muls).nbout; ib++) {
    fprintf(fp1, "\t%g\t%g",
            rPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][0],
            iPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][1]);
    ampl = (float_tt)sqrt(rPart*rPart+iPart*iPart);
    phase = (float_tt)atan2(iPart,rPart);	
    fprintf(fpAmpl,"\t%g",ampl);
    fprintf(fpPhase,"\t%g",phase);
  }
  fprintf( fp1, "\n");
  fprintf( fpAmpl, "\n");
  fprintf( fpPhase, "\n");
	 */
}

void CExperimentCBED::PostSliceProcess(unsigned absoluteSlice)
{
	//	InterimWave(absoluteSlice);
	//	if (_config->Output.LogLevel < 2) {
	//		ComplexArray2D w = m_wave->GetWave();
	//		float_tt potVal = w[0][0].real();
	//		float_tt ddx = potVal;
	//		float_tt ddy = potVal;
	//		for (unsigned ix = 0; ix <  _config->Model.nx; ix++)
	//			for (unsigned iy = 0; iy < _config->Model.ny ; iy++){
	//				potVal = w[ix][iy].real();
	//				if (ddy < potVal)
	//					ddy = potVal;
	//				if (ddx > potVal)
	//					ddx = potVal;
	//			}
	//		BOOST_LOG_TRIVIAL(info)<<format("Saving (complex) wave layer %d to file (r: %g..%g)")%absoluteSlice% ddx% ddy;
	//	}
}

} // end namespace QSTEM
