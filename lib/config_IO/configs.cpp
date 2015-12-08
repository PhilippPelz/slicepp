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

#include "configs.hpp"
#include "readparams.hpp"
#include <boost/algorithm/string.hpp>
#include "string.h"
#include <string>
// used for splitting strings where necessary
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
using boost::format;

namespace QSTEM {

void StructureConfig::Read(ptree& t) {
	structureFilename = boost::filesystem::path(t.get<string>("structure.structure_filename"));
	nCellX = t.get<int>("structure.ncellx");
	nCellY = t.get<int>("structure.ncelly");
	nCellZ = t.get<int>("structure.ncellz");
	temperatureK = t.get<float_tt>("structure.temperatureK");
	crystalTiltX = t.get<float_tt>("structure.crystalTiltX");
	crystalTiltY = t.get<float_tt>("structure.crystalTiltY");
	crystalTiltZ = t.get<float_tt>("structure.crystalTiltZ");
	boxX = t.get<float_tt>("structure.boxX");
	boxY = t.get<float_tt>("structure.boxY");
	boxZ = t.get<float_tt>("structure.boxZ");
	isBoxed = t.get<bool>("structure.isBoxed");
	rotateToZoneAxis = t.get<bool>("structure.rotateToZoneAxis");
	zoneAxis.resize(3);
	string zoneAxisStr = t.get<string>("structure.zoneAxis");
	std::vector<std::string> strs;
	boost::split(strs, zoneAxisStr, boost::is_any_of(","));
	if (strs.size() < 3) {
		BOOST_LOG_TRIVIAL(error)<< format("Zone Axis vector must have 3 dimensions");
	} else {
		std::string::size_type sz;
		for(int i=0;i<3;i++) {
			zoneAxis[i] = std::stoi(strs[i],&sz);
		}
	}

}
StructureConfPtr StructureConfig::Clone() const {
	return StructureConfPtr(new StructureConfig(*this));
}
void ModelConfig::Read(ptree& t) {
	UseTDS = t.get<bool>("model.tds");
	displacementType = static_cast<QSTEM::DisplacementType>(t.get<int>("model.displacementType"));
	TiltBack = t.get<bool>("model.tiltBack");
	CenterSlices = t.get<bool>("model.centerSlices");
	CenterSample = t.get<bool>("model.centerSample");
	TDSRuns = t.get<int>("model.tdsRuns");
	nx = t.get<int>("model.nx");
	ny = t.get<int>("model.ny");
	nSlices = t.get<int>("model.slices");
	SliceThicknessCalculation = static_cast<QSTEM::SliceThicknessCalculation>(t.get<int>("model.sliceThicknessCalculation"));
	ResolutionCalculation = static_cast<QSTEM::ResolutionCalculation>(t.get<int>("model.resolutionCalculation"));
	dz = t.get<float_tt>("model.sliceThicknessAngstrom");
	dx = t.get<float_tt>("model.resolutionXAngstrom");
	dy = t.get<float_tt>("model.resolutionYAngstrom");
	beamTiltX = t.get<float_tt>("model.beamTiltX");
	beamTiltY = t.get<float_tt>("model.beamTiltY");
	SourceDiameterAngstrom = t.get<float_tt>("beam.sourceDiameterAngstrom");
	BeamCurrentpA = t.get<float_tt>("beam.beamCurrentpA");
	PlotVrr = t.get<bool>("model.plotVr_r");
	periodicXY = t.get<bool>("model.periodicXY");
	periodicZ = t.get<bool>("model.periodicZ");
	StructureFactorType = static_cast<QSTEM::StructureFactorType>(t.get<int>("model.structureFactors"));
	CUDAOnTheFly = t.get<bool>("model.CUDAOnTheFly");
	PotentialType = t.get<std::string>("model.type");
	ratom = t.get<float_tt>("model.atomRadiusAngstrom");
	DoZInterpolation = t.get<bool>("model.DoZInterpolation");
	UseQPotentialOffsets = t.get<bool>("model.UseQPotentialOffsets");
	imPot = t.get<float>("wave.imaginary potential factor");
	xOffset = t.get<float_tt>("structure.xOffset");
	yOffset = t.get<float_tt>("structure.yOffset");
	zOffset = t.get<float_tt>("structure.zOffset");

	EnergykeV = t.get<float_tt>("beam.energy_keV");
	float_tt w;
	const float_tt emass = 510.99906; /* electron rest mass in keV */
	const float_tt hc = 12.3984244; /* Planck's const x speed of light*/

	/* electron wavelength in Angstroms */
	wavelength = hc / sqrt(EnergykeV * (2 * emass + EnergykeV));

	float_tt s, pi, x;
	x = (emass + EnergykeV) / (2.0 * emass + EnergykeV);
	pi = 4.0 * atan(1.0);
	sigma = 2.0 * pi * x / (wavelength * EnergykeV);  // 2*pi*kz*(1+kev/emaxx)/(2*emass+kev)
}

void OutputConfig::Read(ptree& t) {
	LogLevel = t.get<int>("output.loglevel");
	SaveWaveIterations = t.get<int>("output.SaveWaveAfterNSlices");
	PendelloesungPlot = t.get<bool>("output.pendelloesungPlot");
	SavePotential = t.get<bool>("output.savePotential");
	SaveWaveAfterTransmit = t.get<bool>("output.SaveWaveAfterTransmit");
	SaveWaveAfterTransform = t.get<bool>("output.SaveWaveAfterTransform");
	SaveWaveAfterSlice = t.get<bool>("output.SaveWaveAfterSlice");
	SaveWaveAfterPropagation = t.get<bool>("output.SaveWaveAfterPropagation");
	saveProbe = t.get<bool>("output.saveProbe");
	SaveAtomicPotential = t.get<bool>("output.SaveAtomicPotential");
	SaveProjectedPotential = t.get<bool>("output.saveProjectedPotential");
	SaveAtomDeltas = t.get<bool>("output.SaveAtomDeltas");
	readPotential = t.get<bool>("output.readPotential");
	string p = t.get<string>("output.savePath");
	savePath = boost::filesystem::path(p);
	string p1 = t.get<string>("output.logFileName");
	LogFileName = boost::filesystem::path(p1);
	WriteLogFile = t.get<bool>("output.writeLogFile");
	ComputeFromProjectedPotential = t.get<bool>("output.ComputeFromProjectedPotential");
}
void WaveConfig::Read(ptree& t) {
	type = t.get<int>("wave.type");
	Cs = t.get<float_tt>("wave.Cs");
	C5 = t.get<float_tt>("wave.C5");
	Cc = t.get<float_tt>("wave.Cc");
	dV_V = t.get<float_tt>("wave.dV/V");
	alpha = t.get<float_tt>("wave.alpha");
	Defocus = t.get<float_tt>("wave.defocus");
	Astigmatism = t.get<float_tt>("wave.astigmatism");
	AstigmatismAngle = t.get<float_tt>("wave.astigmatismAngle");
	Smooth = t.get<bool>("wave.smooth");
	Gaussian = t.get<bool>("wave.gaussian");
	a_33 = t.get<float_tt>("wave.a_33");
	a_31 = t.get<float_tt>("wave.a_31");
	a_44 = t.get<float_tt>("wave.a_44");
	a_42 = t.get<float_tt>("wave.a_42");
	a_55 = t.get<float_tt>("wave.a_55");
	a_53 = t.get<float_tt>("wave.a_53");
	a_51 = t.get<float_tt>("wave.a_51");
	a_66 = t.get<float_tt>("wave.a_66");
	a_64 = t.get<float_tt>("wave.a_64");
	a_62 = t.get<float_tt>("wave.a_62");
	phi_33 = t.get<float_tt>("wave.phi_33");
	phi_31 = t.get<float_tt>("wave.phi_31");
	phi_44 = t.get<float_tt>("wave.phi_44");
	phi_42 = t.get<float_tt>("wave.phi_42");
	phi_55 = t.get<float_tt>("wave.phi_55");
	phi_53 = t.get<float_tt>("wave.phi_53");
	phi_51 = t.get<float_tt>("wave.phi_51");
	phi_66 = t.get<float_tt>("wave.phi_66");
	phi_64 = t.get<float_tt>("wave.phi_64");
	phi_62 = t.get<float_tt>("wave.phi_62");
	gaussScale = t.get<float_tt>("wave.gaussScale");
	dI_I = t.get<float_tt>("wave.dI/I");
	dE_E = t.get<float_tt>("wave.dE/E");
	AISaperture = t.get<float_tt>("wave.AISaperture");
	tiltX = t.get<float_tt>("wave.tiltX");
	tiltY = t.get<float_tt>("wave.tiltY");
	nx = t.get<int>("wave.nx");
	ny = t.get<int>("wave.ny");

	pixelDose = t.get<float_tt>("wave.pixel dose");

}

void ScanConfig::Read(ptree& t) {

	xPos = t.get<int>("scan.x Start Position");
	yPos = t.get<int>("scan.y Start Position");
	xStep = t.get<int>("scan.xStep");
	yStep = t.get<int>("scan.yStep");
	scanType = t.get<int>("scan.scanType");
}
void DetectorConfig::Read(ptree& t) {
	type = t.get<int>("detector.type");
	mtfA = t.get<float_tt>("detector.mtfA");
	mtfB = t.get<float_tt>("detector.mtfB");
	mtfC = t.get<float_tt>("detector.mtfC");
	mtfD = t.get<float_tt>("detector.mtfD");
	DwellTimeMsec = t.get<float_tt>("beam.dwellTimeMsec");
	nx = t.get<int>("wave.nx");
	ny = t.get<int>("wave.ny");
}
Config::Config(ptree& t, boost::filesystem::path configPath) {
	ExperimentType = static_cast<QSTEM::ExperimentType>(t.get<int>("mode"));
	nThreads = t.get<int>("nthreads");

	Structure = boost::shared_ptr<StructureConfig>(new StructureConfig());
	Model = boost::shared_ptr<ModelConfig>(new ModelConfig());
	Output = boost::shared_ptr<OutputConfig>(new OutputConfig());
	Wave = boost::shared_ptr<WaveConfig>(new WaveConfig());
	Scan = boost::shared_ptr<ScanConfig>(new ScanConfig());
	Detector = boost::shared_ptr<DetectorConfig>(new DetectorConfig());

	Output->configPath = configPath;
	Structure->Read(t);
	Model->Read(t);
	Output->Read(t);
	Wave->Read(t);
	Scan->Read(t);
	Detector->Read(t);
}

} // end namespace QSTEM
