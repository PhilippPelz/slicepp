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

#ifndef READ_QSC_H
#define READ_QSC_H

#include "config_interface.hpp"
#include <boost/filesystem.hpp>

namespace QSTEM {

class StructureConfig;
typedef boost::shared_ptr<const StructureConfig> cStructureConfPtr;
typedef boost::shared_ptr<StructureConfig> StructureConfPtr;

class DLL_EXPORT StructureConfig: IPropertyTreeReader {
public:
	boost::filesystem::path structureFilename;
	std::vector<int> zoneAxis;
	int nCellX = 8;
	int nCellY = 8;
	int nCellZ = 8;
	bool rotateToZoneAxis;
	float_tt temperatureK;
	float_tt crystalTiltX;
	float_tt crystalTiltY;
	float_tt crystalTiltZ;
	float_tt boxX;
	float_tt boxY;
	float_tt boxZ;
	bool isBoxed;
	virtual void Read(ptree& t);
	StructureConfPtr Clone() const;
};

class DLL_EXPORT ModelConfig: IPropertyTreeReader {
public:
	bool UseTDS;
	bool TiltBack;
	bool CenterSlices;
	bool CenterSample;
	bool rotateToZoneAxis;
	int TDSRuns, nx, ny, nSlices;
	QSTEM::SliceThicknessCalculation SliceThicknessCalculation;
	QSTEM::ResolutionCalculation ResolutionCalculation;
	float_tt dz, dx, dy, beamTiltX, beamTiltY, SourceDiameterAngstrom, BeamCurrentpA, xOffset, yOffset, zOffset;
	std::vector<int> zoneAxis;
	DisplacementType displacementType;
	bool Use3D, UseFFT, CUDAOnTheFly, PlotVrr, periodicXY, periodicZ, DoZInterpolation, UseQPotentialOffsets;
	QSTEM::StructureFactorType StructureFactorType;
	//atom radius in angstrom
	float_tt ratom;
	std::string PotentialType;
	float_tt EnergykeV, wavelength, sigma, _gamma;
	float_tt imPot;

	virtual void Read(ptree& t);
	bool HasOffset() const {
		return xOffset != 0 || yOffset != 0 || zOffset != 0;
	}
};

class DLL_EXPORT WaveConfig: IPropertyTreeReader {
public:
	float_tt Cs, C5, Cc, dV_V, alpha, Defocus, Astigmatism, AstigmatismAngle, a_33, a_31, a_44, a_42, a_55, a_53, a_51, a_66, a_64, a_62, phi_33,
			phi_31, phi_44, phi_42, phi_55, phi_53, phi_51, phi_66, phi_64, phi_62, gaussScale, dI_I, dE_E, AISaperture, tiltX, tiltY, posX, posY,
			pixelDose;

	bool Smooth, Gaussian;
	int type, nx, ny;

	virtual void Read(ptree& t);
};

class DLL_EXPORT OutputConfig: IPropertyTreeReader {
public:
	int LogLevel, SaveWaveIterations;
	bool SavePotential, SaveProjectedPotential, WriteLogFile, saveProbe, SaveWaveAfterTransmit, SaveWaveAfterTransform, SaveWaveAfterPropagation,
			SaveWaveAfterSlice, SaveAtomicPotential, ComputeFromProjectedPotential,
			SaveAtomDeltas,SaveAtomConv;
	boost::filesystem::path savePath, LogFileName;
	boost::filesystem::path configPath;

	// TODO: deprecated
	bool PendelloesungPlot, readPotential;

	virtual void Read(ptree& t);
};
class DLL_EXPORT DetectorConfig: IPropertyTreeReader {
public:
	float_tt mtfA, mtfB, mtfC, mtfD, DwellTimeMsec;
	int type;
	int nx, ny;
	virtual void Read(ptree& t);
};
class DLL_EXPORT ScanConfig: IPropertyTreeReader {
public:
	int xPos, yPos, xStep, yStep, scanType;
	virtual void Read(ptree& t);
};

class DLL_EXPORT Config {
public:
	Config() :
			nThreads(1), ExperimentType(QSTEM::ExperimentType::CBED) {
	}
	;
	Config(ptree& t, boost::filesystem::path configPath);
	int nThreads;


	QSTEM::ExperimentType ExperimentType;

	boost::shared_ptr<StructureConfig> Structure;
	boost::shared_ptr<ModelConfig> Model;
	boost::shared_ptr<OutputConfig> Output;
	boost::shared_ptr<WaveConfig> Wave;
	boost::shared_ptr<ScanConfig> Scan;
	boost::shared_ptr<DetectorConfig> Detector;
};


typedef boost::shared_ptr<const Config> cConfigPtr;
typedef boost::shared_ptr<Config> ConfigPtr;

typedef boost::shared_ptr<const ModelConfig> cModelConfPtr;
typedef boost::shared_ptr<ModelConfig> ModelConfPtr;

typedef boost::shared_ptr<const OutputConfig> cOutputConfPtr;
typedef boost::shared_ptr< OutputConfig> OutputConfPtr;

typedef boost::shared_ptr<const WaveConfig> cWaveConfPtr;
typedef boost::shared_ptr< WaveConfig> WaveConfPtr;

typedef boost::shared_ptr<const ScanConfig> cScanConfPtr;
typedef boost::shared_ptr< ScanConfig> ScanConfPtr;

typedef boost::shared_ptr<const DetectorConfig> cDetectorConfPtr;
typedef boost::shared_ptr< DetectorConfig> DetectorConfPtr;

} // end namespace QSTEM

#endif
