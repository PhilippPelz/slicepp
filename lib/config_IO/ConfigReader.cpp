/*
 * ConfigReader.cpp
 *
 *  Created on: Dec 14, 2015
 *      Author: philipp
 */

#include "ConfigReader.hpp"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

namespace QSTEM {

ConfigReader::ConfigReader() {
	// TODO Auto-generated constructor stub

}

ConfigReader::~ConfigReader() {
	// TODO Auto-generated destructor stub
}

ConfigPtr ConfigReader::Read(boost::filesystem::path configFile){
	ConfigPtr c = ConfigPtr(new Config());

	ptree t;
	bpt::json_parser::read_json(configFile.string().c_str(),pt);

	c->Output->configPath = boost::filesystem::path(_configFile).parent_path();
	ExperimentType = static_cast<QSTEM::ExperimentType>(t.get<int>("mode"));
	nThreads = t.get<int>("nthreads");

	c->Detector->type = t.get<int>("detector.type");
	c->Detector->mtfA = t.get<float_tt>("detector.mtfA");
	c->Detector->mtfB = t.get<float_tt>("detector.mtfB");
	c->Detector->mtfC = t.get<float_tt>("detector.mtfC");
	c->Detector->mtfD = t.get<float_tt>("detector.mtfD");
	c->Detector->DwellTimeMsec = t.get<float_tt>("beam.dwellTimeMsec");
	c->Detector->nx = t.get<int>("wave.nx");
	c->Detector->ny = t.get<int>("wave.ny");

	c->Scan->xPos = t.get<int>("scan.x Start Position");
	c->Scan->yPos = t.get<int>("scan.y Start Position");
	c->Scan->xStep = t.get<int>("scan.xStep");
	c->Scan->yStep = t.get<int>("scan.yStep");
	c->Scan->scanType = t.get<int>("scan.scanType");

	c->Wave->type = t.get<int>("wave.type");
	c->Wave->Cs = t.get<float_tt>("wave.Cs");
	c->Wave->C5 = t.get<float_tt>("wave.C5");
	c->Wave->Cc = t.get<float_tt>("wave.Cc");
	c->Wave->dV_V = t.get<float_tt>("wave.dV/V");
	c->Wave->alpha = t.get<float_tt>("wave.alpha");
	c->Wave->Defocus = t.get<float_tt>("wave.defocus");
	c->Wave->Astigmatism = t.get<float_tt>("wave.astigmatism");
	c->Wave->AstigmatismAngle = t.get<float_tt>("wave.astigmatismAngle");
	c->Wave->Smooth = t.get<bool>("wave.smooth");
	c->Wave->Gaussian = t.get<bool>("wave.gaussian");
	c->Wave->a_33 = t.get<float_tt>("wave.a_33");
	c->Wave->a_31 = t.get<float_tt>("wave.a_31");
	c->Wave->a_44 = t.get<float_tt>("wave.a_44");
	c->Wave->a_42 = t.get<float_tt>("wave.a_42");
	c->Wave->a_55 = t.get<float_tt>("wave.a_55");
	c->Wave->a_53 = t.get<float_tt>("wave.a_53");
	c->Wave->a_51 = t.get<float_tt>("wave.a_51");
	c->Wave->a_66 = t.get<float_tt>("wave.a_66");
	c->Wave->a_64 = t.get<float_tt>("wave.a_64");
	c->Wave->a_62 = t.get<float_tt>("wave.a_62");
	c->Wave->phi_33 = t.get<float_tt>("wave.phi_33");
	c->Wave->phi_31 = t.get<float_tt>("wave.phi_31");
	c->Wave->phi_44 = t.get<float_tt>("wave.phi_44");
	c->Wave->phi_42 = t.get<float_tt>("wave.phi_42");
	c->Wave->phi_55 = t.get<float_tt>("wave.phi_55");
	c->Wave->phi_53 = t.get<float_tt>("wave.phi_53");
	c->Wave->phi_51 = t.get<float_tt>("wave.phi_51");
	c->Wave->phi_66 = t.get<float_tt>("wave.phi_66");
	c->Wave->phi_64 = t.get<float_tt>("wave.phi_64");
	c->Wave->phi_62 = t.get<float_tt>("wave.phi_62");
	c->Wave->gaussScale = t.get<float_tt>("wave.gaussScale");
	c->Wave->dI_I = t.get<float_tt>("wave.dI/I");
	c->Wave->dE_E = t.get<float_tt>("wave.dE/E");
	c->Wave->AISaperture = t.get<float_tt>("wave.AISaperture");
	c->Wave->tiltX = t.get<float_tt>("wave.tiltX");
	c->Wave->tiltY = t.get<float_tt>("wave.tiltY");
	c->Wave->nx = t.get<int>("wave.nx");
	c->Wave->ny = t.get<int>("wave.ny");
	c->Wave->pixelDose = t.get<float_tt>("wave.pixel dose");

	c->Output->LogLevel = t.get<int>("output.loglevel");
	c->Output->SaveWaveIterations = t.get<int>("output.SaveWaveAfterNSlices");
	c->Output->PendelloesungPlot = t.get<bool>("output.pendelloesungPlot");
	c->Output->SavePotential = t.get<bool>("output.savePotential");
	c->Output->SaveWaveAfterTransmit = t.get<bool>("output.SaveWaveAfterTransmit");
	c->Output->SaveWaveAfterTransform = t.get<bool>("output.SaveWaveAfterTransform");
	c->Output->SaveWaveAfterSlice = t.get<bool>("output.SaveWaveAfterSlice");
	c->Output->SaveWaveAfterPropagation = t.get<bool>("output.SaveWaveAfterPropagation");
	c->Output->saveProbe = t.get<bool>("output.saveProbe");
	c->Output->SaveAtomicPotential = t.get<bool>("output.SaveAtomicPotential");
	c->Output->SaveProjectedPotential = t.get<bool>("output.saveProjectedPotential");
	c->Output->SaveAtomDeltas = t.get<bool>("output.SaveAtomDeltas");
	c->Output->SaveAtomConv = t.get<bool>("output.SaveAtomConv");
	c->Output->readPotential = t.get<bool>("output.readPotential");
	string p = t.get<string>("output.savePath");
	c->Output->savePath = boost::filesystem::path(p);
	string p1 = t.get<string>("output.logFileName");
	c->Output->LogFileName = boost::filesystem::path(p1);
	c->Output->WriteLogFile = t.get<bool>("output.writeLogFile");
	c->Output->ComputeFromProjectedPotential = t.get<bool>("output.ComputeFromProjectedPotential");

	c->Model->UseTDS = t.get<bool>("model.tds");
	c->Model->displacementType = static_cast<QSTEM::DisplacementType>(t.get<int>("model.displacementType"));
	c->Model->TiltBack = t.get<bool>("model.tiltBack");
	c->Model->CenterSlices = t.get<bool>("model.centerSlices");
	c->Model->CenterSample = t.get<bool>("model.centerSample");
	c->Model->TDSRuns = t.get<int>("model.tdsRuns");
	c->Model->nx = t.get<int>("model.nx");
	c->Model->ny = t.get<int>("model.ny");
	c->Model->nSlices = t.get<int>("model.slices");
	c->Model->SliceThicknessCalculation = static_cast<QSTEM::SliceThicknessCalculation>(t.get<int>("model.sliceThicknessCalculation"));
	c->Model->ResolutionCalculation = static_cast<QSTEM::ResolutionCalculation>(t.get<int>("model.resolutionCalculation"));
	c->Model->dz = t.get<float_tt>("model.sliceThicknessAngstrom");
	c->Model->dx = t.get<float_tt>("model.resolutionXAngstrom");
	c->Model->dy = t.get<float_tt>("model.resolutionYAngstrom");
	c->Model->beamTiltX = t.get<float_tt>("model.beamTiltX");
	c->Model->beamTiltY = t.get<float_tt>("model.beamTiltY");
	c->Model->SourceDiameterAngstrom = t.get<float_tt>("beam.sourceDiameterAngstrom");
	c->Model->BeamCurrentpA = t.get<float_tt>("beam.beamCurrentpA");
	c->Model->PlotVrr = t.get<bool>("model.plotVr_r");
	c->Model->periodicXY = t.get<bool>("model.periodicXY");
	c->Model->periodicZ = t.get<bool>("model.periodicZ");
	c->Model->StructureFactorType = static_cast<QSTEM::StructureFactorType>(t.get<int>("model.structureFactors"));
	c->Model->CUDAOnTheFly = t.get<bool>("model.CUDAOnTheFly");
	c->Model->PotentialType = t.get<std::string>("model.type");
	c->Model->ratom = t.get<float_tt>("model.atomRadiusAngstrom");
	c->Model->DoZInterpolation = t.get<bool>("model.DoZInterpolation");
	c->Model->UseQPotentialOffsets = t.get<bool>("model.UseQPotentialOffsets");
	c->Model->ImagPot = t.get<float>("wave.imaginary potential factor");
	c->Model->xOffset = t.get<float_tt>("structure.xOffset");
	c->Model->yOffset = t.get<float_tt>("structure.yOffset");
	c->Model->zOffset = t.get<float_tt>("structure.zOffset");

	c->Model->EnergykeV = t.get<float_tt>("beam.energy_keV");

    const float E0 = EnergykeV*1e3;
    const float m0 = 9.1093822f;
    const float c0  = 2.9979246f;
    const float e  = 1.6021766f;
    const float h  = 6.6260696f;
    const float pi = 3.1415927f;

    c->Model->_gamma = 1.f + E0 * e / m0 / c0 / c0 * 1e-4f;
    c->Model->wavelength = h / sqrt ( 2.f * m0 * e ) * 1e-9f / sqrt ( E0 * ( 1.f + E0 * e / 2.f / m0 / c0 / c0 * 1e-4f ) ); // electron wavelength (m), see De Graef p. 92
    c->Model->sigma =  2.f * pi * _gamma * wavelength * m0 * e / h / h * 1e18f; // interaction constant (1/(Vm))

    c->Structure->structureFilename = boost::filesystem::path(t.get<string>("structure.structure_filename"));
    c->Structure->nCellX = t.get<int>("structure.ncellx");
    c->Structure->nCellY = t.get<int>("structure.ncelly");
    c->Structure->nCellZ = t.get<int>("structure.ncellz");
    c->Structure->temperatureK = t.get<float_tt>("structure.temperatureK");
    c->Structure->crystalTiltX = t.get<float_tt>("structure.crystalTiltX");
    c->Structure->crystalTiltY = t.get<float_tt>("structure.crystalTiltY");
    c->Structure->crystalTiltZ = t.get<float_tt>("structure.crystalTiltZ");
    c->Structure->boxX = t.get<float_tt>("structure.boxX");
    c->Structure->boxY = t.get<float_tt>("structure.boxY");
    c->Structure->boxZ = t.get<float_tt>("structure.boxZ");
    c->Structure->isBoxed = t.get<bool>("structure.isBoxed");
    c->Model->rotateToZoneAxis = t.get<bool>("structure.rotateToZoneAxis");
    c->Model->zoneAxis.resize(3);
	string zoneAxisStr = t.get<string>("structure.zoneAxis");
	std::vector<std::string> strs;
	boost::split(strs, zoneAxisStr, boost::is_any_of(","));
	if (strs.size() < 3) {
		BOOST_LOG_TRIVIAL(error)<< format("Zone Axis vector must have 3 dimensions");
	} else {
		std::string::size_type sz;
		for(int i=0;i<3;i++) {
			c->Model->zoneAxis[i] = std::stoi(strs[i],&sz);
		}
	}
}

} /* namespace QSTEM */
