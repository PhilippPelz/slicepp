local ffi = require 'ffi'

ffi.cdef([[
typedef struct StructureConfig {
	const char* structureFilename;
	int nCellX;
	int nCellY;
	int nCellZ;
	float temperatureK;
	float crystalTiltX;
	float crystalTiltY;
	float crystalTiltZ;
	float boxX;
	float boxY;
	float boxZ;
	bool isBoxed;
	int zoneAxis[3];
	bool rotateToZoneAxis;
} StructureConfig;

typedef struct ModelConfig {

	bool UseTDS;
	bool TiltBack;
	bool CenterSlices;
	bool CenterSample;

	int TDSRuns, nx, ny, nSlices;

	float dz;
	float dx;
	float dy;
	float beamTiltX;
	float beamTiltY;
	float SourceDiameterAngstrom;
	float BeamCurrentpA;
	float xOffset;
	float yOffset;
	float zOffset;


	bool Use3D;
	bool UseFFT;
	bool CUDAOnTheFly;
	bool PlotVrr;
	bool periodicXY;
	bool periodicZ;
	bool DoZInterpolation;
	bool UseQPotentialOffsets;

	int StructureFactorType;
	int SliceThicknessCalculation;
	int ResolutionCalculation;
	int displacementType;

	//atom radius in angstrom
	float ratom;
	const char* PotentialType;
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
	float dV_V;
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
	float AISaperture;
	float tiltX;
	float tiltY;
	float posX;
	float posY, pixelDose;

	bool Smooth;
	bool Gaussian;
	int type;
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

	const char* LogFileName;
	const char* SavePath;
	const char* ConfigPath;

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
	int type;
	int nx;
	int ny;
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
	int ExperimentType;

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

void run_simulation(c_Config* conf);
]])

local C = ffi.load('slice')
return C
