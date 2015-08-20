/*
 * Bootstrapper.h
 *
 *  Created on: Mar 31, 2015
 *      Author: philipp
 */

#ifndef STEM3_BOOTSTRAPPER_H_
#define STEM3_BOOTSTRAPPER_H_

#include "config_IO/configs.hpp"
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/functional/factory.hpp>
#include <boost/function.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;
using boost::property_tree::info_parser::read_info;
namespace bpt = boost::property_tree;
//using namespace std;

#include <wavefunctions/wave_base.hpp>
#include <wavefunctions/wave_interface.hpp>
#include <wavefunctions/wave_convergent.hpp>
#include <wavefunctions/wave_plane.hpp>
#include "experiments/base.hpp"
#include "experiments/CoherentCBED.hpp"
#include "experiments/CoherentTEM.hpp"
#include "experiments/Ptychograph.hpp"
#include "potentials/pot_interface.hpp"
#include "potentials/pot_2d.hpp"
#include "potentials/pot_2d_fft.hpp"
#include "potentials/pot_3d.hpp"
#include "potentials/pot_3d_fft.hpp"
#include "potentials/CUDA2DPotential.hpp"
#include "detectors/detector_interface.hpp"
#include "detectors/FlatAreaDetector.hpp"
#include "stemtypes_fftw3.hpp"
#include "structure_IO/structureInterface.hpp"
#include "structure_IO/IStructureBuilder.hpp"
#include "structure_IO/SuperstructureBuilder.hpp"
#include "structure_IO/CifReader.hpp"

//#define BOOST_FUNCTIONAL_FACTORY_SUPPORT_NONE_T

namespace QSTEM {

class Bootstrapper {
public:
	Bootstrapper(int argc, char *argv[]);
	virtual ~Bootstrapper();
	void Initialize();
	ExperimentPtr GetExperiment();
protected:
	ConfigPtr _c;
	ExperimentPtr _e;
	void RegisterWaveTypes();
	void RegisterExperimentTypes();
	void RegisterPotentialTypes();
	void RegisterDetectorTypes();
	void RegisterStructureTypes();
	void RegisterStructureBuilders();
	void RegisterStructureReaders();
	void SetSliceThickness(StructurePtr s);
	void SetResolution(StructurePtr s);
	DetectorFactory _detectorFactory;
	WaveFactory _waveFactory;
	ExperimentFactory _experimentFactory;
	PotentialFactory _potentialFactory;
	StructureReaderFactory _structureReaderFactory;
	StructureBuilderFactory _structureBuilderFactory;
	std::string _configFile;
	boost::filesystem::path _configPath;
};

} /* namespace QSTEM */

#endif /* STEM3_BOOTSTRAPPER_H_ */
