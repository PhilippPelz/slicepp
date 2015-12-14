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

struct StructureConfig;
typedef boost::shared_ptr<const StructureConfig> cStructureConfPtr;
typedef boost::shared_ptr<StructureConfig> StructureConfPtr;

typedef struct DLL_EXPORT StructureConfig{
public:
	boost::filesystem::path structureFilename;
	std::vector<int> zoneAxis;
	int nCellX;
	int nCellY;
	int nCellZ;
	bool rotateToZoneAxis;
	float_tt temperatureK;
	float_tt crystalTiltX;
	float_tt crystalTiltY;
	float_tt crystalTiltZ;
	float_tt boxX;
	float_tt boxY;
	float_tt boxZ;
	bool isBoxed;

	StructureConfPtr Clone() const;
} StructureConfig;

typedef struct DLL_EXPORT ModelConfig{
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
	bool Use3D;
	bool UseFFT;
	bool CUDAOnTheFly;
	bool PlotVrr;
	bool periodicXY;
	bool periodicZ;
	bool DoZInterpolation;
	bool UseQPotentialOffsets;
	QSTEM::StructureFactorType StructureFactorType;
	//atom radius in angstrom
	float_tt ratom;
	std::string PotentialType;
	float_tt EnergykeV, wavelength, sigma, _gamma;
	float_tt ImagPot;


	bool HasOffset() const {
		return xOffset != 0 || yOffset != 0 || zOffset != 0;
	}
} ModelConfig;

typedef struct DLL_EXPORT WaveConfig{
public:
	float_tt Cs;
	float_tt C5;
	float_tt Cc;
	float_tt dV_V;
	float_tt alpha;
	float_tt Defocus;
	float_tt Astigmatism;
	float_tt AstigmatismAngle;
	float_tt a_33;
	float_tt a_31;
	float_tt a_44;
	float_tt a_42;
	float_tt a_55;
	float_tt a_53;
	float_tt a_51;
	float_tt a_66;
	float_tt a_64;
	float_tt a_62;
	float_tt phi_33;
	float_tt phi_31;
	float_tt phi_44;
	float_tt phi_42;
	float_tt phi_55;
	float_tt phi_53;
	float_tt phi_51;
	float_tt phi_66;
	float_tt phi_64;
	float_tt phi_62;
	float_tt gaussScale;
	float_tt dI_I;
	float_tt dE_E;
	float_tt AISaperture;
	float_tt tiltX;
	float_tt tiltY;
	float_tt posX;
	float_tt posY, pixelDose;

	bool Smooth;
	bool Gaussian;
	int type;
	int nx;
	int ny;


} WaveConfig;

typedef struct DLL_EXPORT OutputConfig{
public:
	int LogLevel;
	int SaveWaveIterations;
	bool SavePotential;
	bool SaveProjectedPotential;
	bool WriteLogFile;
	bool saveProbe;
	bool SaveWaveAfterTransmit;
	bool SaveWaveAfterTransform;
	bool SaveWaveAfterPropagation;
	bool SaveWaveAfterSlice;
	bool SaveAtomicPotential;
	bool ComputeFromProjectedPotential, SaveAtomDeltas;
	bool SaveAtomConv;
	boost::filesystem::path savePath;
	bool LogFileName;
	boost::filesystem::path configPath;

	// TODO: deprecated
	bool PendelloesungPlot;
	bool readPotential;


} OutputConfig;

typedef struct DLL_EXPORT DetectorConfig{
public:
	float_tt mtfA;
	float_tt mtfB;
	float_tt mtfC;
	float_tt mtfD;
	float_tt DwellTimeMsec;
	int type;
	int nx;
	int ny;

} DetectorConfig;

typedef struct DLL_EXPORT ScanConfig{
public:
	int xPos;
	int yPos;
	int xStep;
	int yStep;
	int scanType;

} ScanConfig;

typedef struct DLL_EXPORT Config {
public:
	Config();

	int nThreads;
	ExperimentType ExperimentType;

	boost::shared_ptr<StructureConfig> Structure;
	boost::shared_ptr<ModelConfig> Model;
	boost::shared_ptr<OutputConfig> Output;
	boost::shared_ptr<WaveConfig> Wave;
	boost::shared_ptr<ScanConfig> Scan;
	boost::shared_ptr<DetectorConfig> Detector;
} Config;

typedef boost::shared_ptr<const Config> cConfigPtr;
typedef boost::shared_ptr<Config> ConfigPtr;

typedef boost::shared_ptr<const ModelConfig> cModelConfPtr;
typedef boost::shared_ptr<ModelConfig> ModelConfPtr;

typedef boost::shared_ptr<const OutputConfig> cOutputConfPtr;
typedef boost::shared_ptr<OutputConfig> OutputConfPtr;

typedef boost::shared_ptr<const WaveConfig> cWaveConfPtr;
typedef boost::shared_ptr<WaveConfig> WaveConfPtr;

typedef boost::shared_ptr<const ScanConfig> cScanConfPtr;
typedef boost::shared_ptr<ScanConfig> ScanConfPtr;

typedef boost::shared_ptr<const DetectorConfig> cDetectorConfPtr;
typedef boost::shared_ptr<DetectorConfig> DetectorConfPtr;

} // end namespace QSTEM

#endif
