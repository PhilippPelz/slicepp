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

#ifndef configs_hpp
#define configs_hpp

typedef enum {
	STEM = 3, CBED = 1, TEM = 4, NBED = 2, PTYCHO = 5
} ExperimentType;

typedef enum  {
	Auto = 1, SliceThickness = 2, NumberOfSlices = 3
} SliceThicknessCalculation;

typedef enum  {
	WeickKohl = 1, Rez = 2
} StructureFactorType;

typedef enum  {
	FILLRES = 1, FILLN = 2, BOXRES = 3, BOXN = 4
} ResolutionCalculation;

typedef enum  {
	Einstein = 1, Phonon = 2, None = 3
} DisplacementType;

typedef enum  {
	Plane = 1, Convergent =2
} WaveType;

typedef enum  {
	Scintillator = 1, Direct = 2
} DetectorType;

typedef enum  {
	FFT2D = 1,
	FFT3D = 2,
	Real2D = 3,
	Real3D = 4,
	CUDA2D = 5
} PotentialType;

typedef struct StructureConfig {
	char StructureFilename[1000];

	float T_Kelvin;

	float crystalTilt[3];
	float box[3];
	int zoneAxis[3];
	int nCells[3];

	bool isBoxed;
	bool rotateToZoneAxis;

} StructureConfig;

typedef struct ModelConfig {
	bool TiltBack;
	bool CenterSample;

	bool UseTDS;
	int TDSRuns;

	int n[3];
	float d[3];
	float offset[3];

	float beamTiltX;
	float beamTiltY;
	float SourceDiameterAngstrom;
	float BeamCurrentpA;

	bool PlotVrr;
	bool periodicXY;
	bool periodicZ;
	bool DoZInterpolation;
	bool UseQPotentialOffsets;

	StructureFactorType StructFactorType;
	SliceThicknessCalculation SliceCalcType;
	ResolutionCalculation ResCalcType;
	DisplacementType DisplaceType;

	//atom radius in angstrom
	float ratom;
	PotentialType PotType;
	float EnergykeV;
	float wavelength;
	float sigma;
	float _gamma;
	float ImagPot;
} ModelConfig;

typedef struct WaveConfig {

	float Cs;
	float C5;
	float Cc;

	float alpha;
	float Defocus;
	float Astigmatism;
	float AstigmatismAngle;

	float a_33;
	float a_31;
	float a_44;
	float a_42;
	float a_55;
	float a_53;
	float a_51;
	float a_66;
	float a_64;
	float a_62;
	float phi_33;
	float phi_31;
	float phi_44;
	float phi_42;
	float phi_55;
	float phi_53;
	float phi_51;
	float phi_66;
	float phi_64;
	float phi_62;
	float gaussScale;

	float dI_I;
	float dE_E;
	float dV_V;

	float AISaperture;
	float tiltX;
	float tiltY;

	bool Smooth;
	bool Gaussian;

	WaveType type;
	int nx;
	int ny;

} WaveConfig;

typedef struct OutputConfig {

	int LogLevel;
	int SaveWaveIterations;
	bool SavePotential;
	bool SaveProjectedPotential;
	bool WriteLogFile;
	bool SaveProbe;
	bool SaveWaveAfterTransmit;
	bool SaveWaveAfterTransform;
	bool SaveWaveAfterPropagation;
	bool SaveWaveAfterSlice;
	bool SaveAtomicPotential;
	bool ComputeFromProjectedPotential, SaveAtomDeltas;
	bool SaveAtomConv;

	char LogFileName[1000];
	char SavePath[1000];
	char ConfigPath[1000];

	// TODO: deprecated
	bool PendelloesungPlot;
	bool readPotential;

} OutputConfig;

typedef struct DetectorConfig {
	float mtfA;
	float mtfB;
	float mtfC;
	float mtfD;
	float DwellTimeMsec;
	DetectorType type;
	int n[2];
	float MaxElectronCounts;
} DetectorConfig;

typedef struct ScanConfig {
	int xPos;
	int yPos;
	int xStep;
	int yStep;
	int scanType;
} ScanConfig;

typedef struct c_Config {
	int nThreads;
	ExperimentType ExpType;

	StructureConfig* Structure;
	ModelConfig* Model;
	OutputConfig* Output;
	WaveConfig* Wave;
	ScanConfig* Scan;
	DetectorConfig* Detector;
} c_Config;

StructureConfig* StructureConfig_clone(const StructureConfig* cloneMe);

StructureConfig* StructureConfig_new();
ModelConfig* ModelConfig_new();
WaveConfig* WaveConfig_new();
OutputConfig* OutputConfig_new();
DetectorConfig* DetectorConfig_new();
ScanConfig* ScanConfig_new();
c_Config* c_Config_new();

void StructureConfig_delete(StructureConfig* c);
void ModelConfig_delete(ModelConfig* c);
void WaveConfig_delete(WaveConfig);
void OutputConfig_delete(OutputConfig* c);
void DetectorConfig_delete(DetectorConfig* c);
void ScanConfig_delete(ScanConfig* c);
void c_Config_delete(c_Config* c);

#endif
