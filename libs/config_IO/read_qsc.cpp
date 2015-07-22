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

#include "read_qsc.hpp"
#include "readparams.hpp"
#include <boost/algorithm/string.hpp>
#include "string.h"
#include <string>
// used for splitting strings where necessary
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
using boost::format;

namespace QSTEM
{
void StructureConfig::Read(ptree& t){
	structureFilename=boost::filesystem::path(t.get<string>("structure.structure_filename"));
	nCellX=t.get<int>("structure.ncellx");
	nCellY=t.get<int>("structure.ncelly");
	nCellZ=t.get<int>("structure.ncellz");
	temperatureK=t.get<float_tt>("structure.temperatureK");
	crystalTiltX=t.get<float_tt>("structure.crystalTiltX");
	crystalTiltY=t.get<float_tt>("structure.crystalTiltY");
	crystalTiltZ=t.get<float_tt>("structure.crystalTiltZ");
	xOffset=t.get<float_tt>("structure.xOffset");
	yOffset=t.get<float_tt>("structure.yOffset");
	boxX=t.get<float_tt>("structure.boxX");
	boxY=t.get<float_tt>("structure.boxY");
	boxZ=t.get<float_tt>("structure.boxZ");
	isBoxed=t.get<bool>("structure.isBoxed");
	rotateToZoneAxis=t.get<bool>("structure.rotateToZoneAxis");
	zoneAxis.resize(3);
	string zoneAxisStr=t.get<string>("structure.zoneAxis");
	std::vector<std::string> strs;
	boost::split(strs, zoneAxisStr, boost::is_any_of(","));
	if(strs.size() < 3){
		BOOST_LOG_TRIVIAL(error) << format("Zone Axis vector must have 3 dimensions");
	} else {
		std::string::size_type sz;
		for(int i=0;i<3;i++){
			zoneAxis[i] = std::stoi(strs[i],&sz);
		}
	}
	xOffset=t.get<float_tt>("structure.xOffset");
	yOffset=t.get<float_tt>("structure.yOffset");
	zOffset=t.get<float_tt>("structure.zOffset");
}
void ModelConfig::Read(ptree& t){
	UseTDS=t.get<bool>("model.tds");
	displacementType = static_cast<QSTEM::DisplacementType>(t.get<int>("model.displacementType"));
	TiltBack=t.get<bool>("model.tiltBack");
	CenterSlices=t.get<bool>("model.centerSlices");
	CenterSample=t.get<bool>("model.centerSample");
	TDSRuns=t.get<int>("model.tdsRuns");
	nx=t.get<int>("model.nx");
	ny=t.get<int>("model.ny");
	nSlices =t.get<int>("model.slices");
	SliceThicknessCalculation = static_cast<QSTEM::SliceThicknessCalculation>(t.get<int>("model.sliceThicknessCalculation"));
	ResolutionCalculation = static_cast<QSTEM::ResolutionCalculation>(t.get<int>("model.resolutionCalculation"));
	dz=t.get<float_tt>("model.sliceThicknessAngstrom");
	dx=t.get<float_tt>("model.resolutionXAngstrom");
	dy=t.get<float_tt>("model.resolutionYAngstrom");
	beamTiltX=t.get<float_tt>("model.beamTiltX");
	beamTiltY=t.get<float_tt>("model.beamTiltY");
}
void PotentialConfig::Read(ptree& t){
	PlotVrr=t.get<bool>("model.potential.plotVr_r");
	Use3D=t.get<bool>("model.potential.3D");
	UseFFT=t.get<bool>("model.potential.FFT");
	periodicXY=t.get<bool>("model.potential.periodicXY");
	periodicZ=t.get<bool>("model.potential.periodicZ");
	StructureFactorType  = static_cast<QSTEM::StructureFactorType>(t.get<int>("model.potential.structureFactors"));
	AtomRadiusAngstrom=t.get<float_tt>("model.potential.atomRadiusAngstrom");
}
void OutputConfig::Read(ptree& t){
	LogLevel=t.get<int>("output.loglevel");
	SaveWaveIterations=t.get<int>("output.SaveWaveAfterNSlices");
	ShowProbe=t.get<bool>("output.showProbe");
	PendelloesungPlot=t.get<bool>("output.pendelloesungPlot");
	SavePotential=t.get<bool>("output.savePotential");
	SaveWaveAfterTransmit=t.get<bool>("output.SaveWaveAfterTransmit");
	SaveWaveAfterTransform=t.get<bool>("output.SaveWaveAfterTransform");
	SaveWaveAfterSlice=t.get<bool>("output.SaveWaveAfterSlice");
	SaveWaveAfterPropagation=t.get<bool>("output.SaveWaveAfterPropagation");
	saveProbe = t.get<bool>("output.saveProbe");
	SaveAtomicPotential=t.get<bool>("output.SaveAtomicPotential");
	SaveProjectedPotential=t.get<bool>("output.saveProjectedPotential");
	readPotential=t.get<bool>("output.readPotential");
	string p = t.get<string>("output.savePath");
	savePath = boost::filesystem::path(p);
	string p1 = t.get<string>("output.logFileName");
	LogFileName = boost::filesystem::path(p1);
	WriteLogFile = t.get<bool>("output.writeLogFile");
}
void WaveConfig::Read(ptree& t){
	type=t.get<int>("wave.type");
	Cs=t.get<float_tt>("wave.Cs");
	C5=t.get<float_tt>("wave.C5");
	Cc=t.get<float_tt>("wave.Cc");
	dV_V=t.get<float_tt>("wave.dV/V");
	alpha=t.get<float_tt>("wave.alpha");
	Defocus=t.get<float_tt>("wave.defocus");
	Astigmatism=t.get<float_tt>("wave.astigmatism");
	AstigmatismAngle=t.get<float_tt>("wave.astigmatismAngle");
	Smooth=t.get<bool>("wave.smooth");
	Gaussian=t.get<bool>("wave.gaussian");
	a_33=t.get<float_tt>("wave.a_33");
	a_31=t.get<float_tt>("wave.a_31");
	a_44=t.get<float_tt>("wave.a_44");
	a_42=t.get<float_tt>("wave.a_42");
	a_55=t.get<float_tt>("wave.a_55");
	a_53=t.get<float_tt>("wave.a_53");
	a_51=t.get<float_tt>("wave.a_51");
	a_66=t.get<float_tt>("wave.a_66");
	a_64=t.get<float_tt>("wave.a_64");
	a_62=t.get<float_tt>("wave.a_62");
	phi_33=t.get<float_tt>("wave.phi_33");
	phi_31=t.get<float_tt>("wave.phi_31");
	phi_44=t.get<float_tt>("wave.phi_44");
	phi_42=t.get<float_tt>("wave.phi_42");
	phi_55=t.get<float_tt>("wave.phi_55");
	phi_53=t.get<float_tt>("wave.phi_53");
	phi_51=t.get<float_tt>("wave.phi_51");
	phi_66=t.get<float_tt>("wave.phi_66");
	phi_64=t.get<float_tt>("wave.phi_64");
	phi_62=t.get<float_tt>("wave.phi_62");
	gaussScale=t.get<float_tt>("wave.gaussScale");
	dI_I=t.get<float_tt>("wave.dI/I");
	dE_E=t.get<float_tt>("wave.dE/E");
	AISaperture=t.get<float_tt>("wave.AISaperture");
	tiltX=t.get<float_tt>("wave.tiltX");
	tiltY=t.get<float_tt>("wave.tiltY");
//	posX=t.get<float_tt>("wave.posX");
//	posY=t.get<float_tt>("wave.posY");
}
void BeamConfig::Read(ptree& t){
	EnergykeV =t.get<float_tt>("beam.energy_keV");
	SourceDiameterAngstrom=t.get<float_tt>("beam.sourceDiameterAngstrom");
	BeamCurrentpA=t.get<float_tt>("beam.beamCurrentpA");
	DwellTimeMsec=t.get<float_tt>("beam.dwellTimeMsec");

	float_tt w;
	const float_tt emass=510.99906; /* electron rest mass in keV */
	const float_tt hc=12.3984244; /* Planck's const x speed of light*/

	/* electron wavelength in Angstroms */
	wavelength = hc/sqrt( EnergykeV * ( 2*emass + EnergykeV ) );

	float_tt s, pi, x;
	x = ( emass + EnergykeV ) / ( 2.0*emass + EnergykeV);
	pi = 4.0 * atan( 1.0 );
	sigma = 2.0 * pi * x / (wavelength*EnergykeV);  // 2*pi*kz*(1+kev/emaxx)/(2*emass+kev)

}

void ScanConfig::Read(ptree& t){
	Scanx = t.get<int>("scan.scanx");
	Scany = t.get<int>("scan.scany");
	xPos = t.get<int>("scan.x Start Position");
	yPos = t.get<int>("scan.y Start Position");
	xStep = t.get<int>("scan.xStep");
	yStep = t.get<int>("scan.yStep");
	scanType = t.get<int>("scan.scanType");
}

Config::Config(ptree& t){
	ExperimentType =static_cast<QSTEM::ExperimentType>(t.get<int>("mode"));
	nThreads = t.get<int>("nthreads");
	Structure =StructureConfig();
	Model =ModelConfig();
	Potential =PotentialConfig();
	Output =OutputConfig();
	Wave =WaveConfig();
	Beam = BeamConfig();
	Scan = ScanConfig();


	Structure.Read(t);
	Model.Read(t);
	Potential.Read(t);
	Output.Read(t);
	Wave.Read(t);
	Beam.Read(t);
	Scan.Read(t);
}


CQscReader::CQscReader(boost::filesystem::path &filename) : IConfigReader()
{
	m_fp = fopen(filename.string().c_str(), "r" );
	if( m_fp == NULL ) {
		BOOST_LOG_TRIVIAL(error) << format("Cannot open file %s\n") % filename.string().c_str();
		m_isValid=false;
	}
	else
	{
		m_isValid=true;
	}
}

CQscReader::~CQscReader()
{
	// make sure the file is closed
	if(m_fp!=NULL)
		fclose(m_fp);
}

void CQscReader::ReadMode(std::string &mode)
{
	readparam(m_fp,"mode:",mode,1);
}

void CQscReader::ReadPrintLevel(unsigned int &printLevel)
{
	if (readparam(m_fp,"print level:",m_buf,1))
		printLevel = atoi(m_buf.c_str());
}

void CQscReader::ReadSaveLevel(unsigned int &saveLevel)
{
	if (readparam(m_fp,"save level:",m_buf,1))
		saveLevel = atoi(m_buf.c_str());

}

void CQscReader::ReadPotentialOutputInterval(unsigned &displayPotCalcInterval)
{
	if (readparam(m_fp,"potential progress interval:",m_buf,1))
		displayPotCalcInterval = atoi(m_buf.c_str());
}

void CQscReader::ReadSTEMProgressInterval(unsigned int &displayProgInterval)
{
	if (readparam(m_fp,"propagation progress interval:",m_buf,1))
		displayProgInterval = atoi(m_buf.c_str());
}

// TODO: this should remove quotes from the output directory/filename
void CQscReader::ReadStructureFileName(boost::filesystem::path &structure_path)
{
	if (!readparam(m_fp,"filename:",m_buf,1)) exit(0);
	structure_path = boost::filesystem::path(m_buf);
}

void CQscReader::ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ)
{
	if (readparam(m_fp,"NCELLX:",m_buf,1)) nCellX=atoi(m_buf.c_str());
	if (readparam(m_fp,"NCELLY:",m_buf,1)) nCellY=atoi(m_buf.c_str());

	// NCellZ is a bit of a special case - the number of sub-slabs is stored in the same line.
	//  Need to split that line, and take only the first part here.
	if (readparam(m_fp,"NCELLZ:",m_buf,1)) {
		std::vector<std::string> values;
		boost::split(values, m_buf, boost::is_any_of(" "));
		nCellZ=atoi(values[0].c_str());
		//sscanf(m_buf,"%s",answer);
		//nCellZ = atoi(answer);
	}
}

void CQscReader::ReadNSubSlabs(unsigned &cellDiv)
{
	cellDiv=0;
	// This is stored on the same line as nCellZ, hence the duplication.
	//    We don't store nCellZ again here.
	if (readparam(m_fp, "NCELLZ:",m_buf,1)) {
		std::vector<std::string> values;
		boost::split(values, m_buf, boost::is_any_of(" "));
		if (values.size()>1)
			cellDiv=atoi(values[1].c_str());
		// cellDiv=1 when cell is undivided
		if (cellDiv==0)
			cellDiv=1;
		//sscanf(m_buf,"%s",answer);
		//if ((strPtr = strchr(answer,'/')) != NULL) {
		//strPtr[0] = '\0';
		//  cellDiv = atoi(strPtr+1);
		//}
	}
}

/*************************************************
 * Read the beam tilt parameters
 */
void CQscReader::ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool &tiltBack)
{
	std::vector<std::string> values;

	if (readparam(m_fp, "Beam tilt X:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		btiltx=atof(values[0].c_str());
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			btiltx *= PI/180.0;
	}
	if (readparam(m_fp, "Beam tilt Y:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		btilty=atof(values[0].c_str());
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			btilty *= PI/180.0;
	}
	if (readparam(m_fp, "Tilt back:",m_buf,1)) {
		tiltBack = IsBufferYes(m_buf);
	}

}
void CQscReader::ReadCrystalCubeAndTilt(float_tt &tiltx, float_tt &tilty, float_tt &tiltz, 
		float_tt &cubex, float_tt &cubey, float_tt &cubez,
		bool &adjustCubeSize)
{
	std::vector<std::string> values;
	/*************************************************
	 * Read the crystal tilt parameters
	 */
	answer[0] = '\0';
	if (readparam(m_fp, "Crystal tilt X:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		tiltx=atof(values[0].c_str());
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			tiltx *= PI/180.0;
	}
	answer[0] = '\0';
	if (readparam(m_fp, "Crystal tilt Y:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		tilty=atof(values[0].c_str());
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			tilty *= PI/180.0;
	}
	answer[0] = '\0';
	if (readparam(m_fp, "Crystal tilt Z:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		tiltz=atof(values[0].c_str());
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			tiltz *= PI/180.0;
	}
	cubex = 0; cubey = 0; cubez = 0; /* in A */
	if (readparam(m_fp, "Cube:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		cubex=atof(values[0].c_str());
		cubey=atof(values[1].c_str());
		cubez=atof(values[2].c_str());
	}

	adjustCubeSize = false;
	if (readparam(m_fp, "Adjust cube size with tilt:",m_buf,1)) {
		adjustCubeSize = IsBufferYes(m_buf);
	}
}

void CQscReader::ReadTemperatureData(bool &doTDS, float_tt &tdsTemperature, boost::filesystem::path &phononFile, 
		bool &useEinstein)
{
	if (readparam(m_fp, "tds:",m_buf,1)) {
		doTDS = IsBufferYes(m_buf);
	}
	else doTDS = false;

	if (readparam(m_fp, "temperature:",m_buf,1)) tdsTemperature=atof(m_buf.c_str());
	else tdsTemperature = 300.0;

	useEinstein = true;  //phononFile = NULL;
	std::string readout;
	if (readparam(m_fp, "phonon-File:",readout,1)) {
		useEinstein = false;
		// convert the read out string into a boost::filesystem::path
		phononFile=readout;
	}
}

void CQscReader::ReadSliceOffset(float_tt &xOffset, float_tt &yOffset)
{
	xOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam(m_fp, "xOffset:",m_buf,1)) xOffset=atof(m_buf.c_str());
	yOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam(m_fp, "yOffset:",m_buf,1)) yOffset=atof(m_buf.c_str());
	// printf("Reading Offset: %f, %f\n",xOffset,yOffset);
}

void CQscReader::ReadProbeArraySize(unsigned &nx, unsigned &ny)
{
	// Required parameter
	if (!readparam(m_fp, "nx:",m_buf,1)) exit(0);
	nx=atoi(m_buf.c_str());
	if (readparam(m_fp, "ny:",m_buf,1)) ny=atoi(m_buf.c_str());
	else ny = nx;
}

void CQscReader::ReadResolution(float_tt &resolutionX, float_tt &resolutionY)
{
	if (readparam(m_fp, "resolutionX:",m_buf,1)) resolutionX=atof(m_buf.c_str());
	if (readparam(m_fp, "resolutionY:",m_buf,1)) resolutionY=atof(m_buf.c_str());
}

void CQscReader::ReadVoltage(float_tt &voltage)
{
	// Required parameter
	if (!readparam(m_fp, "v0:",m_buf,1)) exit(0); voltage=atof(m_buf.c_str());
}

void CQscReader::ReadSliceParameters(bool &centerSlices, float_tt &sliceThickness, 
		unsigned &nslices, unsigned &outputInterval,
		float_tt &zOffset)
{
	if (readparam(m_fp, "center slices:",m_buf,1)) {
		centerSlices = IsBufferYes(m_buf);
	}
	// just in case the answer was not exactly 1 or 0:
	// centerSlices = (centerSlices) ? 1 : 0;

	if (readparam(m_fp, "slice-thickness:",m_buf,1)) sliceThickness=atof(m_buf.c_str());
	if (readparam(m_fp, "slices:",m_buf,1)) nslices=atoi(m_buf.c_str());

	// read the output interval:
	if (readparam(m_fp, "slices between outputs:",m_buf,1)) outputInterval=atoi(m_buf.c_str());
	if (readparam(m_fp, "zOffset:",m_buf,1)) zOffset=atof(m_buf.c_str());
}

void CQscReader::ReadPeriodicParameters(bool &periodicXY, bool &periodicZ)
{
	if (readparam(m_fp, "periodicXY:",m_buf,1)) {
		periodicXY = IsBufferYes(m_buf);
	}
	if (readparam(m_fp, "periodicZ:",m_buf,1)) {
		periodicZ = IsBufferYes(m_buf);
	}
}

void CQscReader::ReadBandLimitTrans(bool &bandlimittrans)
{
	if (readparam(m_fp, "bandlimit f_trans:",m_buf,1)) {
		bandlimittrans = IsBufferYes(m_buf);
	}
}

void CQscReader::ReadLoadPotential(bool &readPotential, boost::filesystem::path &filename)
{
	if (readparam(m_fp, "read potential:",m_buf,1)) {
		readPotential = IsBufferYes(m_buf);
	}
	if (readPotential)
	{
		// TODO: what is the format for reading the potential filename?
	}
}

void CQscReader::ReadPotentialOutputParameters(bool &savePotential, bool &saveTotalPotential, 
		bool &plotPotential)
{
	if (readparam(m_fp, "save potential:",m_buf,1)) {
		savePotential = IsBufferYes(m_buf);
	}
	if (readparam(m_fp, "save projected potential:",m_buf,1)) {
		saveTotalPotential = IsBufferYes(m_buf);
	}
	if (readparam(m_fp, "plot V(r)*r:",m_buf,1)) {
		plotPotential = IsBufferYes(m_buf);
	}
}

void CQscReader::ReadPotentialCalculationParameters(bool &fftPotential, bool &potential3D)
{
	if (readparam(m_fp, "one time integration:",m_buf,1)) {
		fftPotential = IsBufferYes(m_buf);
	}
	if (readparam(m_fp, "potential3D:",m_buf,1)) {
		potential3D = IsBufferYes(m_buf);
	}
}

void CQscReader::ReadAverageParameters(unsigned &avgRuns, bool &storeSeries)
{
	if (readparam(m_fp, "Runs for averaging:",m_buf,1))
		avgRuns=atoi(m_buf.c_str());

	if (readparam(m_fp, "Store TDS diffr. patt. series:",m_buf,1)) {
		storeSeries = IsBufferYes(m_buf);
	}
}

void CQscReader::ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
		float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN)
{
	if (readparam(m_fp, "scan_x_start:",m_buf,1))  scanXStart=atof(m_buf.c_str());
	if (readparam(m_fp, "scan_x_stop:",m_buf,1))   scanXStop=atof(m_buf.c_str());
	if (readparam(m_fp, "scan_x_pixels:",m_buf,1)) scanXN=atoi(m_buf.c_str());
	if (readparam(m_fp, "scan_y_start:",m_buf,1))  scanYStart=atof(m_buf.c_str());
	if (readparam(m_fp, "scan_y_stop:",m_buf,1))   scanYStop=atof(m_buf.c_str());
	if (readparam(m_fp, "scan_y_pixels:",m_buf,1)) scanYN=atoi(m_buf.c_str());
}


void CQscReader::ReadOutputName(boost::filesystem::path &fileOrFolderName)
{
	/**********************************************************************
	 * Parameters for image display and directories, etc.
	 */
	fileOrFolderName = "data";
	//sprintf(fileOrFolderName,"data");
	std::string readout;
	readparam(m_fp, "Folder:",readout,1);
	// convert the read back file or folder name into a boost::filesystem::path
	fileOrFolderName=readout;
}

void CQscReader::ReadAtomRadius(float_tt &radius)
{
	/*******************************************************************
	 * Read in parameters related to the calculation of the projected
	 * Potential
	 *******************************************************************/
	if (readparam(m_fp, "atom radius:",m_buf,1))
		radius=atof(m_buf.c_str()); /* in A */
}

void CQscReader::ReadStructureFactorType(std::string &scatFactor)
{
	readparam(m_fp, "Structure Factors:",scatFactor,1);
}

void CQscReader::ReadPendelloesungParameters(std::vector<int> &hbeams, std::vector<int> &kbeams, 
		bool &lbeams, unsigned &nbout)
{
	/*************************************************************
	 * read in the beams we want to plot in the pendeloesung plot
	 * Only possible if not in STEM or CBED mode
	 *************************************************************/
	resetParamFile(m_fp);
	if (readparam(m_fp, "Pendelloesung plot:",m_buf,1)) {
		lbeams = IsBufferYes(m_buf);
	}
	if (lbeams) {
		std::vector<std::string> values;

		while (readparam(m_fp, "beam:",m_buf,0)) nbout++;
		BOOST_LOG_TRIVIAL(trace) << format("will record %d beams") % nbout;
		hbeams.resize(nbout);
		kbeams.resize(nbout);
		/* now read in the list of detectors: */
		resetParamFile(m_fp);
		for (unsigned i=0;i<nbout;i++) {
			if (!readparam(m_fp, "beam:",m_buf,0)) break;
			boost::split(values, m_buf, boost::is_any_of(" "));
			hbeams[i] = atoi(values[0].c_str());
			kbeams[i] = atoi(values[0].c_str());
			//sscanf(m_buf,"%d %d",&(hbeams[i]),&(kbeams[i]));
		}
	}
}

void CQscReader::ReadNumberOfDetectors(int &numDetectors)
{
	numDetectors=0;
	resetParamFile(m_fp);
	while (readparam(m_fp, "detector:",m_buf, 0)) ++numDetectors;
}

void CQscReader::ReadDetectorParameters(int det_idx, float_tt &rInside, float_tt &rOutside, std::string &name, 
		float_tt &shiftX, float_tt &shiftY)
{
	int file_det_idx=0;
	std::vector<std::string> values;

	resetParamFile(m_fp);
	while (readparam(m_fp, "detector:",m_buf,0)) {
		if (det_idx==file_det_idx)
		{
			boost::split(values, m_buf, boost::is_any_of(" "));
			rInside=atof(values[0].c_str());
			rOutside=atof(values[1].c_str());
			name=values[2];
			shiftX=atof(values[3].c_str());
			shiftY=atof(values[4].c_str());
			return;
		}
		++file_det_idx;
	}
}

void CQscReader::ReadDoseParameters(float_tt &beamCurrent, float_tt &dwellTimeMs)
{
	///// read beam current and dwell time ///////////////////////////////
	if (readparam(m_fp, "beam current:",m_buf,1)) {
		beamCurrent=atof(m_buf.c_str()); /* in pA */
	}
	if (readparam(m_fp, "dwell time:",m_buf,1)) {
		dwellTimeMs=atof(m_buf.c_str()); /* in msec */
	}
}

void CQscReader::ReadProbeParameters(float_tt &dE_E, float_tt &dI_I, float_tt &dV_V, float_tt &alpha, float_tt &aAIS)
{
	/**********************************************************************
	 * Read STEM/CBED probe parameters
	 */
	if (readparam(m_fp, "dE/E:",m_buf,1))
		dE_E = atof(m_buf.c_str());
	if (readparam(m_fp, "dI/I:",m_buf,1))
		dI_I = atof(m_buf.c_str());
	if (readparam(m_fp, "dV/V:",m_buf,1))
		dV_V = atof(m_buf.c_str());

	// required parameter
	if (!readparam(m_fp, "alpha:",m_buf,1)) exit(0);
	alpha=atof(m_buf.c_str()); /* in mrad */

	if (readparam(m_fp, "AIS aperture:",m_buf,1))
		aAIS=atof(m_buf.c_str()); /* in A */
}
void CQscReader::ReadSourceRadius(float_tt &sourceRadius)
{
	if (readparam(m_fp, "Source Size (diameter):",m_buf,1))
		sourceRadius = atof(m_buf.c_str())/2.0;
}

void CQscReader::ReadSmoothingParameters(bool &smooth, float_tt &gaussScale, bool &gaussFlag)
{
	if (readparam(m_fp, "smooth:",m_buf,1)) smooth = IsBufferYes(m_buf);
	if (readparam(m_fp, "gaussian:",m_buf,1)) {
		std::vector<std::string> values;
		boost::split(values, m_buf, boost::is_any_of(" "));
		gaussFlag = IsBufferYes(values[0]);
		if (gaussFlag && values.size()>1)
			gaussScale=atof(values[1].c_str());
	}
}

void CQscReader::ReadTomoParameters(float_tt &tomoTilt, float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
		float_tt &zoomFactor)
{
	std::vector<std::string> values;
	// in case this file has been written by the tomography function, read the current tilt:
	if (readparam(m_fp, "tomo tilt:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		tomoTilt=atof(values[0].c_str());  /* in mrad */
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			tomoTilt *= 1000*PI/180.0;
	}

	if (readparam(m_fp, "tomo start:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		tomoStart=atof(values[0].c_str());  /* in mrad */
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			tomoStart *= 1000*PI/180.0;
	}
	if (readparam(m_fp, "tomo step:",m_buf,1)) {
		boost::split(values, m_buf, boost::is_any_of(" "));
		tomoStep=atof(values[0].c_str());  /* in mrad */
		if (values.size()>1 && ('d'==tolower(values[1][0])))
			tomoStep *= 1000*PI/180.0;
	}

	if (readparam(m_fp, "tomo count:",m_buf,1))
		tomoCount = atoi(m_buf.c_str());
	if (readparam(m_fp, "zoom factor:",m_buf,1))
		zoomFactor=atof(m_buf.c_str());
	if ((tomoStep == 0) && (tomoStep > 1))
		tomoStep = -2.0*tomoStart/(float_tt)(tomoCount - 1);
}

void CQscReader::ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5, float_tt &Cc,
		float_tt &df0, std::string &Scherzer, float_tt &astig,
		float_tt &a33, float_tt &a31,
		float_tt &a44, float_tt &a42,
		float_tt &a55, float_tt &a53, float_tt &a51,
		float_tt &a66, float_tt &a64, float_tt &a62)
{
	// require that people specify Cs
	if (!readparam(m_fp, "Cs:",m_buf,1))  exit(0);
	Cs=atof(m_buf.c_str()); /* in mm */

	if (readparam(m_fp, "C5:",m_buf,1)) {
		C5=atof(m_buf.c_str()); /* in mm */
	}
	if (readparam(m_fp, "Cc:",m_buf,1))
		Cc = atof(m_buf.c_str());

	if (readparam(m_fp, "defocus:",m_buf,1)) {
		if (isalpha(m_buf[0]))
			Scherzer=m_buf;
		else {
			df0=atof(m_buf.c_str());
			Scherzer = "No";
		}
	}
	// Astigmatism:
	if (readparam(m_fp, "astigmatism:",m_buf,1)) astig=atof(m_buf.c_str());

	if (readparam(m_fp, "a_33:",m_buf,1)) {a33=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_31:",m_buf,1)) {a31=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_44:",m_buf,1)) {a44=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_42:",m_buf,1)) {a42=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_55:",m_buf,1)) {a55=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_53:",m_buf,1)) {a53=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_51:",m_buf,1)) {a51=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_66:",m_buf,1)) {a66=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_64:",m_buf,1)) {a64=atof(m_buf.c_str()); }
	if (readparam(m_fp, "a_62:",m_buf,1)) {a62=atof(m_buf.c_str()); }
}

void CQscReader::ReadAberrationAngles(float_tt &astig,
		float_tt &phi33, float_tt &phi31,
		float_tt &phi44, float_tt &phi42,
		float_tt &phi55, float_tt &phi53, float_tt &phi51,
		float_tt &phi66, float_tt &phi64, float_tt &phi62)
{
	if (readparam(m_fp, "astigmatism angle:",m_buf,1)) astig=atof(m_buf.c_str());

	if (readparam(m_fp, "phi_33:",m_buf,1)) {phi33=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_31:",m_buf,1)) {phi31=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_44:",m_buf,1)) {phi44=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_42:",m_buf,1)) {phi42=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_55:",m_buf,1)) {phi55=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_53:",m_buf,1)) {phi53=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_51:",m_buf,1)) {phi51=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_66:",m_buf,1)) {phi66=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_64:",m_buf,1)) {phi64=atof(m_buf.c_str()); }
	if (readparam(m_fp, "phi_62:",m_buf,1)) {phi62=atof(m_buf.c_str()); }
}

bool CQscReader::IsBufferYes(std::string &buf)
{
	boost::algorithm::to_lower(buf);
	return ("yes"==m_buf);
}

} // end namespace QSTEM
