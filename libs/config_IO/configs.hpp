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

class DLL_EXPORT StructureConfig: IPropertyTreeReader {
public:
	boost::filesystem::path structureFilename;
	std::vector<int> zoneAxis;
	int nCellX = 8;
	int nCellY = 8;
	int nCellZ = 8;
	bool isBoxed, rotateToZoneAxis;
	float_tt temperatureK, crystalTiltX, crystalTiltY, crystalTiltZ, xOffset,
			yOffset, zOffset, boxX, boxY, boxZ;

	virtual void Read(ptree& t);
	bool HasOffset() {
		return xOffset != 0 || yOffset != 0 || zOffset != 0;
	}
	;
};

class DLL_EXPORT ModelConfig: IPropertyTreeReader {
public:
	bool UseTDS, TiltBack, CenterSlices, CenterSample, rotateToZoneAxis;
	int TDSRuns, nx, ny, nSlices;
	QSTEM::SliceThicknessCalculation SliceThicknessCalculation;
	QSTEM::ResolutionCalculation ResolutionCalculation;
	float_tt dz, dx, dy, beamTiltX, beamTiltY,SourceDiameterAngstrom,BeamCurrentpA;
	std::vector<int> zoneAxis;
	DisplacementType displacementType;
	virtual void Read(ptree& t);
};

class DLL_EXPORT PotentialConfig: IPropertyTreeReader {
public:
	bool Use3D, UseFFT, CUDAOnTheFly, PlotVrr,periodicXY,periodicZ,DoZInterpolation,UseQPotentialOffsets;
	QSTEM::StructureFactorType StructureFactorType;
	//atom radius in angstrom
	float_tt ratom;
	std::string PotentialType;
	virtual void Read(ptree& t);
};

class DLL_EXPORT WaveConfig: IPropertyTreeReader {
public:
	float_tt Cs, C5, Cc, dV_V, alpha, Defocus, Astigmatism, AstigmatismAngle,
			a_33, a_31, a_44, a_42, a_55, a_53, a_51, a_66, a_64, a_62, phi_33,
			phi_31, phi_44, phi_42, phi_55, phi_53, phi_51, phi_66, phi_64,
			phi_62, gaussScale, dI_I, dE_E, AISaperture, tiltX, tiltY, posX,
			posY, imPot, pixelDose;
	float_tt EnergykeV, wavelength, sigma;
	bool Smooth, Gaussian;
	int type, nx, ny;

	virtual void Read(ptree& t);
};

class DLL_EXPORT OutputConfig: IPropertyTreeReader {
public:
	int LogLevel, SaveWaveIterations;
	bool SavePotential, SaveProjectedPotential, WriteLogFile, saveProbe,
			SaveWaveAfterTransmit, SaveWaveAfterTransform,
			SaveWaveAfterPropagation, SaveWaveAfterSlice, SaveAtomicPotential, ComputeFromProjectedPotential;
	boost::filesystem::path savePath, LogFileName;

	// TODO: deprecated
	bool PendelloesungPlot, readPotential;

	virtual void Read(ptree& t);
};
class DLL_EXPORT DetectorConfig : IPropertyTreeReader{
public:
	float_tt mtfA, mtfB, mtfC, mtfD,DwellTimeMsec;
	int type;
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

	boost::filesystem::path configPath;
	QSTEM::ExperimentType ExperimentType;

	boost::shared_ptr<StructureConfig> Structure;
	boost::shared_ptr<ModelConfig> Model;
	boost::shared_ptr<PotentialConfig> Potential;
	boost::shared_ptr<OutputConfig> Output;
	boost::shared_ptr<WaveConfig> Wave;
	boost::shared_ptr<ScanConfig> Scan;
	boost::shared_ptr<DetectorConfig> Detector;
};
typedef boost::shared_ptr<Config> ConfigPtr;
} // end namespace QSTEM

#endif
