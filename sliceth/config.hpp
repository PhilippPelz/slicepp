#ifndef __config_hpp__
#define __config_hpp__

typedef struct structure{
	char* structureFilename;
	int zoneAxis[3];
	int nCellX;
	int nCellY;
	int nCellZ;
	bool isBoxed, rotateToZoneAxis;
	float temperatureK;
	float crystalTiltX;
	float crystalTiltY;
	float crystalTiltZ;
	float xOffset;
	float yOffset;
	float zOffset;
	float boxX;
	float boxY;
	float boxZ;
} StructureConfig;

typedef struct beam {
	float E;
	float SourceDiameterAngstrom;
	float BeamCurrentpA;
	float DwellTimeMsec;
	float wavelength;
} BeamConfig;

typedef struct m {
	bool UseTDS, TiltBack, CenterSlices, CenterSample, rotateToZoneAxis;
	int TDSRuns, nx, ny, nSlices;
	int SliceThicknessCalculation;
	int ResolutionCalculation;
	float dz;
	float dx;
	float dy;
	float beamTiltX;
	float beamTiltY;
	int displacementType;
} ModelConfig;

typedef struct potential {

	bool Use3D, UseFFT, CUDAOnTheFly, PlotVrr,periodicXY,periodicZ,DoZInterpolation,UseQPotentialOffsets;
	int StructureFactorType;
	float ratom;
	char* PotentialType;
} PotentialConfig;

typedef struct wave {
	float Cs, C5, Cc, dV_V, alpha, Defocus, Astigmatism, AstigmatismAngle,
			a_33, a_31, a_44, a_42, a_55, a_53, a_51, a_66, a_64, a_62, phi_33,
			phi_31, phi_44, phi_42, phi_55, phi_53, phi_51, phi_66, phi_64,
			phi_62, gaussScale, dI_I, dE_E, AISaperture, tiltX, tiltY, posX,
			posY, imPot, pixelDose;
	bool Smooth, Gaussian;
	int type, nx, ny;
} WaveConfig;

typedef struct out {
	int LogLevel, SaveWaveIterations;
	bool SavePotential, SaveProjectedPotential, WriteLogFile, saveProbe,
			SaveWaveAfterTransmit, SaveWaveAfterTransform,
			SaveWaveAfterPropagation, SaveWaveAfterSlice, SaveAtomicPotential, ComputeFromProjectedPotential;
	char* savePath, LogFileName;
} OutputConfig;

typedef struct d{
	float mtfA, mtfB, mtfC, mtfD;
	int type;
} DetectorConfig;

typedef struct scan {
	int xPos, yPos, xStep, yStep, scanType;
} ScanConfig;

typedef struct conf{
	int nThreads;
	char* configPath;
	int ExperimentType;
	StructureConfig Structure;
	ModelConfig Model;
	PotentialConfig Potential;
	OutputConfig Output;
	WaveConfig Wave;
	BeamConfig Beam;
	ScanConfig Scan;
  DetectorConfig Detector;
} Config;
#endif
