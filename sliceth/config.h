#ifndef __config_hpp__
#define __config_hpp__

typedef struct structure_config{
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
} structure_config;

typedef struct beam_config {
	float E;
	float SourceDiameterAngstrom;
	float BeamCurrentpA;
	float DwellTimeMsec;
	float wavelength;
} beam_config;

typedef struct model_config {
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
} model_config;

typedef struct potential_config {

	bool Use3D, UseFFT, CUDAOnTheFly, PlotVrr,periodicXY,periodicZ,DoZInterpolation,UseQPotentialOffsets;
	int StructureFactorType;
	float ratom;
	char* PotentialType;
} potential_config;

typedef struct wave_config {
	float Cs, C5, Cc, dV_V, alpha, Defocus, Astigmatism, AstigmatismAngle,
			a_33, a_31, a_44, a_42, a_55, a_53, a_51, a_66, a_64, a_62, phi_33,
			phi_31, phi_44, phi_42, phi_55, phi_53, phi_51, phi_66, phi_64,
			phi_62, gaussScale, dI_I, dE_E, AISaperture, tiltX, tiltY, posX,
			posY, imPot, pixelDose;
	bool Smooth, Gaussian;
	int type, nx, ny;
} wave_config;

typedef struct out_config {
	int LogLevel, SaveWaveIterations;
	bool SavePotential, SaveProjectedPotential, WriteLogFile, saveProbe,
			SaveWaveAfterTransmit, SaveWaveAfterTransform,
			SaveWaveAfterPropagation, SaveWaveAfterSlice, SaveAtomicPotential, ComputeFromProjectedPotential;
	char* savePath, LogFileName;
} out_config;

typedef struct det_config{
	float mtfA, mtfB, mtfC, mtfD;
	int type;
} det_config;

typedef struct scan_config {
	int xPos, yPos, xStep, yStep, scanType;
} scan_config;

typedef struct config{
	int nThreads;
	char* configPath;
	int ExperimentType;
	structure_config Structure;
	model_config Model;
	potential_config Potential;
	out_config Output;
	wave_config Wave;
	beam_config Beam;
	scan_config Scan;
  	det_config Detector;
} config;
#endif
