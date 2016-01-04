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
#include <boost/log/trivial.hpp>


namespace slicepp {

using boost::property_tree::ptree;


Config::Config() {
	Structure = boost::shared_ptr<StructureConfig>(StructureConfig_new());
	Model = boost::shared_ptr<ModelConfig>(ModelConfig_new());
	Output = boost::shared_ptr<OutputConfig>(OutputConfig_new());
	Wave = boost::shared_ptr<WaveConfig>(WaveConfig_new());
	Scan = boost::shared_ptr<ScanConfig>(ScanConfig_new());
	Detector = boost::shared_ptr<DetectorConfig>(DetectorConfig_new());
}

Config::Config(c_Config& c) {
	Structure = boost::shared_ptr<StructureConfig>(c.Structure);
	Model = boost::shared_ptr<ModelConfig>(c.Model);
	Output = boost::shared_ptr<OutputConfig>(c.Output);
	Wave = boost::shared_ptr<WaveConfig>(c.Wave);
	Scan = boost::shared_ptr<ScanConfig>(c.Scan);
	Detector = boost::shared_ptr<DetectorConfig>(c.Detector);

	SavePath = boost::filesystem::path(c.Output->SavePath);
	StructureFilename = boost::filesystem::path(c.Structure->structureFilename);
	ConfigPath = boost::filesystem::path(c.Output->ConfigPath);
}

ConfigReader::ConfigReader() {
	// TODO Auto-generated constructor stub

}


ConfigPtr ConfigReader::Read(boost::filesystem::path configFile){
	ConfigPtr c = ConfigPtr(new Config());

	ptree t;
	boost::property_tree::json_parser::read_json(configFile.string().c_str(),t);

	c->Output->ConfigPath = configFile.parent_path().string().c_str();
	c->ExperimentType = t.get<int>("mode");
	c->nThreads = t.get<int>("nthreads");

	c->Detector->type = (DetectorType)t.get<int>("detector.type");
	c->Detector->mtfA = t.get<float_tt>("detector.mtfA");
	c->Detector->mtfB = t.get<float_tt>("detector.mtfB");
	c->Detector->mtfC = t.get<float_tt>("detector.mtfC");
	c->Detector->mtfD = t.get<float_tt>("detector.mtfD");
	c->Detector->DwellTimeMsec = t.get<float_tt>("beam.dwellTimeMsec");
	c->Detector->n[0] = t.get<int>("wave.nx");
	c->Detector->n[1] = t.get<int>("wave.ny");
	c->Detector->MaxElectronCounts = t.get<float_tt>("wave.pixel dose");

	c->Scan->xPos = t.get<int>("scan.x Start Position");
	c->Scan->yPos = t.get<int>("scan.y Start Position");
	c->Scan->xStep = t.get<int>("scan.xStep");
	c->Scan->yStep = t.get<int>("scan.yStep");
	c->Scan->scanType = t.get<int>("scan.scanType");

	c->Wave->type = (WaveType)t.get<int>("wave.type");
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


	c->Output->LogLevel = t.get<int>("output.loglevel");
	c->Output->SaveWaveIterations = t.get<int>("output.SaveWaveAfterNSlices");
	c->Output->PendelloesungPlot = t.get<bool>("output.pendelloesungPlot");
	c->Output->SavePotential = t.get<bool>("output.savePotential");
	c->Output->SaveWaveAfterTransmit = t.get<bool>("output.SaveWaveAfterTransmit");
	c->Output->SaveWaveAfterTransform = t.get<bool>("output.SaveWaveAfterTransform");
	c->Output->SaveWaveAfterSlice = t.get<bool>("output.SaveWaveAfterSlice");
	c->Output->SaveWaveAfterPropagation = t.get<bool>("output.SaveWaveAfterPropagation");
	c->Output->SaveProbe = t.get<bool>("output.saveProbe");
	c->Output->SaveAtomicPotential = t.get<bool>("output.SaveAtomicPotential");
	c->Output->SaveProjectedPotential = t.get<bool>("output.saveProjectedPotential");
	c->Output->SaveAtomDeltas = t.get<bool>("output.SaveAtomDeltas");
	c->Output->SaveAtomConv = t.get<bool>("output.SaveAtomConv");
	c->Output->readPotential = t.get<bool>("output.readPotential");
//	string p = t.get<string>("output.savePath").c_str();
//	c->Output->savePath = boost::filesystem::path(p);
	c->Output->SavePath = t.get<string>("output.savePath").c_str();
//	string p1 = t.get<string>("output.logFileName");
//	c->Output->LogFileName = boost::filesystem::path(p1);
	c->Output->LogFileName = t.get<string>("output.logFileName").c_str();;
	c->Output->WriteLogFile = t.get<bool>("output.writeLogFile");
	c->Output->ComputeFromProjectedPotential = t.get<bool>("output.ComputeFromProjectedPotential");

	c->Model->UseTDS = t.get<bool>("model.tds");
	c->Model->DisplaceType = (DisplacementType)t.get<int>("model.displacementType");
	c->Model->SliceCalcType = (SliceThicknessCalculation)t.get<int>("model.sliceThicknessCalculation");
	c->Model->ResCalcType = (ResolutionCalculation)t.get<int>("model.resolutionCalculation");
	c->Model->StructFactorType = (StructureFactorType)t.get<int>("model.structureFactors");

	c->Model->TiltBack = t.get<bool>("model.tiltBack");
	c->Model->CenterSample = t.get<bool>("model.centerSample");
	c->Model->TDSRuns = t.get<int>("model.tdsRuns");
	c->Model->n[0] = t.get<int>("model.nx");
	c->Model->n[1] = t.get<int>("model.ny");
	c->Model->n[2] = t.get<int>("model.slices");

	c->Model->d[2] = t.get<float_tt>("model.sliceThicknessAngstrom");
	c->Model->d[0] = t.get<float_tt>("model.resolutionXAngstrom");
	c->Model->d[1] = t.get<float_tt>("model.resolutionYAngstrom");
	c->Model->beamTiltX = t.get<float_tt>("model.beamTiltX");
	c->Model->beamTiltY = t.get<float_tt>("model.beamTiltY");
	c->Model->SourceDiameterAngstrom = t.get<float_tt>("beam.sourceDiameterAngstrom");
	c->Model->BeamCurrentpA = t.get<float_tt>("beam.beamCurrentpA");
	c->Model->PlotVrr = t.get<bool>("model.plotVr_r");
	c->Model->periodicXY = t.get<bool>("model.periodicXY");
	c->Model->periodicZ = t.get<bool>("model.periodicZ");

	c->Model->PotentialType = t.get<std::string>("model.type").c_str();
	c->Model->ratom = t.get<float_tt>("model.atomRadiusAngstrom");
	c->Model->DoZInterpolation = t.get<bool>("model.DoZInterpolation");
	c->Model->UseQPotentialOffsets = t.get<bool>("model.UseQPotentialOffsets");
	c->Model->ImagPot = t.get<float>("model.ImagPot");
	c->Model->offset[0] = t.get<float_tt>("structure.xOffset");
	c->Model->offset[1] = t.get<float_tt>("structure.yOffset");
	c->Model->offset[2] = t.get<float_tt>("structure.zOffset");

	c->Model->EnergykeV = t.get<float_tt>("beam.energy_keV");

    const float E0 = c->Model->EnergykeV*1e3;
    const float m0 = 9.1093822f;
    const float c0  = 2.9979246f;
    const float e  = 1.6021766f;
    const float h  = 6.6260696f;
    const float pi = 3.1415927f;

    c->Model->_gamma = 1.f + E0 * e / m0 / c0 / c0 * 1e-4f;
    c->Model->wavelength = h / sqrt ( 2.f * m0 * e ) * 1e-9f / sqrt ( E0 * ( 1.f + E0 * e / 2.f / m0 / c0 / c0 * 1e-4f ) ); // electron wavelength (m), see De Graef p. 92
    c->Model->sigma =  2.f * pi * c->Model->_gamma * c->Model->wavelength * m0 * e / h / h * 1e18f; // interaction constant (1/(Vm))

    c->Structure->structureFilename = t.get<string>("structure.structure_filename").c_str();
    c->Structure->nCells[0] = t.get<int>("structure.ncellx");
    c->Structure->nCells[1] = t.get<int>("structure.ncelly");
    c->Structure->nCells[2] = t.get<int>("structure.ncellz");
    c->Structure->T_Kelvin = t.get<float_tt>("structure.temperatureK");
    c->Structure->crystalTilt[0] = t.get<float_tt>("structure.crystalTiltX");
    c->Structure->crystalTilt[1] = t.get<float_tt>("structure.crystalTiltY");
    c->Structure->crystalTilt[2] = t.get<float_tt>("structure.crystalTiltZ");
    c->Structure->box[0] = t.get<float_tt>("structure.boxX");
    c->Structure->box[1] = t.get<float_tt>("structure.boxY");
    c->Structure->box[1] = t.get<float_tt>("structure.boxZ");
    c->Structure->isBoxed = t.get<bool>("structure.isBoxed");
    c->Structure->rotateToZoneAxis = t.get<bool>("structure.rotateToZoneAxis");
	string zoneAxisStr = t.get<string>("structure.zoneAxis");
	std::vector<std::string> strs;
	boost::split(strs, zoneAxisStr, boost::is_any_of(","));
	if (strs.size() < 3) {
//		BOOST_LOG_TRIVIAL(error)<< format("Zone Axis vector must have 3 dimensions");
	} else {
		std::string::size_type sz;
		for(int i=0;i<3;i++) {
			c->Structure->zoneAxis[i] = std::stoi(strs[i],&sz);
		}
	}
	auto sp = boost::filesystem::path(c->Output->SavePath);
	auto lf = boost::filesystem::path(c->Output->LogFileName);
	auto sf = boost::filesystem::path(c->Structure->structureFilename);
	if(sp.has_relative_path()){
		c->Output->SavePath = (configFile.parent_path() / c->Output->SavePath).string().c_str();
	}
	if(lf.has_relative_path()){
		c->Output->LogFileName = (configFile.parent_path() / c->Output->LogFileName).string().c_str();
	}
	if(sf.has_relative_path()){
		c->Structure->structureFilename = (configFile.parent_path() / c->Structure->structureFilename).string().c_str();
	}

	return c;
}

} /* namespace slicepp */
