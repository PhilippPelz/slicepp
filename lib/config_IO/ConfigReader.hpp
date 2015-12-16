/*
 * ConfigReader.hpp
 *
 *  Created on: Dec 14, 2015
 *      Author: philipp
 */

#ifndef CONFIGREADER_HPP_
#define CONFIGREADER_HPP_

#include "configs.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include "stemtypes_fftw3.hpp"

namespace QSTEM {



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

typedef boost::shared_ptr<const StructureConfig> cStructureConfPtr;
typedef boost::shared_ptr<StructureConfig> StructureConfPtr;

class DLL_EXPORT Config {
public:
	Config();
	Config(c_Config& c);

	int nThreads;
	int ExperimentType;

	boost::filesystem::path StructureFilename;
	boost::filesystem::path SavePath;
	boost::filesystem::path ConfigPath;

	StructureConfPtr Structure;
	ModelConfPtr Model;
	OutputConfPtr Output;
	WaveConfPtr Wave;
	ScanConfPtr Scan;
	DetectorConfPtr Detector;
};

typedef boost::shared_ptr<const Config> cConfigPtr;
typedef boost::shared_ptr<Config> ConfigPtr;

typedef struct ConfigReader {
public:
	ConfigReader();
	ConfigPtr Read(boost::filesystem::path configPath);
} ConfigReader;

} /* namespace QSTEM */

#endif /* CONFIGREADER_HPP_ */
