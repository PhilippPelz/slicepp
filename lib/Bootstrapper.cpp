/*
 * Bootstrapper.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: philipp
 */

#include "Bootstrapper.hpp"
#include "fftw3.h"
#include "arrayfire.h"
namespace slicepp {

//Initial sampling suggestions for calculating images of thin specimens
//1. The transmission function t (x ) (5.25) should be symmetrically
//bandwidth limited.
//2. The real space pixel size Δx and Δy should be less than about
//d0/3 to d0/4 where d0 is the resolution of the instrument.
//3. The reciprocal space pixel size Δkx and Δky should be less than
//about k
//ob j/10 where λkob j is the maximum semiangle in the
//objective aperture.
//4. The maximum angle in reciprocal space λkx−max and λky −max
//should be about twice the maximum angle in the objective aperture (for CTEM) or slightly bigger than the maximum angle on
//the ADF-STEM detector (for STEM)

Bootstrapper::Bootstrapper(int argc, char *argv[]) {
	string config;
	if (argc < 2)   config = "config.json";
	else    config=argv[1];

	_configFile = boost::filesystem::path(config);

	logging::add_console_log(std::cout, keywords::format = ">> %Message%");

}
Bootstrapper::Bootstrapper(){}
Bootstrapper::~Bootstrapper() {}

void Bootstrapper::Initialize(c_Config* conf){
	ConfigPtr c = ConfigPtr(new Config(*conf));

	fftw_init_threads();
	fftw_plan_with_nthreads(c->nThreads);
	omp_set_num_threads(c->nThreads);

	af::info();

	RegisterWaveTypes();
	RegisterPotentialTypes();
	RegisterStructureTypes();
	RegisterExperimentTypes();
	RegisterStructureReaders();
	RegisterStructureBuilders();
	RegisterDetectorTypes();

	logging::core::get()->set_filter
	(
	    logging::trivial::severity >= static_cast<logging::trivial::severity_level>(c->Output->LogLevel)
	);

	if(c->Output->WriteLogFile)
		logging::add_file_log
			(
				keywords::file_name = std::string(c->Output->LogFileName),
				// This makes the sink to write log records that look like this:
				// 1: <normal> A normal severity message
				// 2: <error> An error severity message
				keywords::format =
				(
					expr::stream
						<< expr::attr< unsigned int >("LineID")
						<< ": <" << logging::trivial::severity
						<< "> " << expr::smessage
				),
				keywords::auto_flush = true
			);

	boost::filesystem::path p(c->Structure->structureFilename);

	std::string str = c->Model->PotentialType;
 	std::transform(str.begin(), str.end(),str.begin(), ::toupper);

	auto sreader = StructureReaderPtr(_structureReaderFactory[".cif"](p));
	auto structureBuilder = StructureBuilderPtr(_structureBuilderFactory[p.extension().string()](sreader,c->Structure,c->Model,c->Output));
	auto persist = PersistenceManagerPtr(new PersistenceManager(c));
	auto w = WaveConfPtr(c->Wave);
	auto m = ModelConfPtr(c->Model);
	auto wave = WavePtr(_waveFactory[c->Wave->type](w,m,persist));
	auto detector = DetPtr(_detectorFactory[c->Detector->type](c->Detector,persist));
	auto potential = PotPtr(_potentialFactory[str](c->Model,c->Output,persist));
	_e = ExperimentPtr( _experimentFactory[c->ExperimentType](c,structureBuilder,wave,potential,detector, persist));
}
void Bootstrapper::Initialize(){

	auto cr = ConfigReader();
	auto c = cr.Read(_configFile);

	fftw_init_threads();
	fftw_plan_with_nthreads(c->nThreads);
	omp_set_num_threads(c->nThreads);

	af::info();

	RegisterWaveTypes();
	RegisterPotentialTypes();
	RegisterStructureTypes();
	RegisterExperimentTypes();
	RegisterStructureReaders();
	RegisterStructureBuilders();
	RegisterDetectorTypes();

	logging::core::get()->set_filter
	(
	    logging::trivial::severity >= static_cast<logging::trivial::severity_level>(c->Output->LogLevel)
	);

	if(c->Output->WriteLogFile)
		logging::add_file_log
			(
				keywords::file_name = std::string(c->Output->LogFileName),
				// This makes the sink to write log records that look like this:
				// 1: <normal> A normal severity message
				// 2: <error> An error severity message
				keywords::format =
				(
					expr::stream
						<< expr::attr< unsigned int >("LineID")
						<< ": <" << logging::trivial::severity
						<< "> " << expr::smessage
				),
				keywords::auto_flush = true
			);

	boost::filesystem::path p(c->Structure->structureFilename);

	std::string str = c->Model->PotentialType;
 	std::transform(str.begin(), str.end(),str.begin(), ::toupper);

	auto sreader = StructureReaderPtr(_structureReaderFactory[".cif"](p));
	auto structureBuilder = StructureBuilderPtr(_structureBuilderFactory[p.extension().string()](sreader,c->Structure,c->Model,c->Output));
	auto persist = PersistenceManagerPtr(new PersistenceManager(c));
	auto w = WaveConfPtr(c->Wave);
	auto m = ModelConfPtr(c->Model);
	auto wave = WavePtr(_waveFactory[c->Wave->type](w,m,persist));
	auto detector = DetPtr(_detectorFactory[c->Detector->type](c->Detector,persist));
	auto potential = PotPtr(_potentialFactory[str](c->Model,c->Output,persist));
	_e = ExperimentPtr( _experimentFactory[c->ExperimentType](c,structureBuilder,wave,potential,detector, persist));
}

ExperimentPtr Bootstrapper::GetExperiment(){
	return _e;
}
void Bootstrapper::RegisterStructureBuilders(){
	_structureBuilderFactory[".gbm"] = boost::factory<SuperstructureBuilder*>();
	_structureBuilderFactory[".cif"] = boost::factory<CrystalBuilder*>();
}
void Bootstrapper::RegisterStructureReaders(){
	_structureReaderFactory[".cif"] = boost::factory<CifReader*>();
}

void Bootstrapper::RegisterWaveTypes(){
	_waveFactory[1] = boost::factory<CPlaneWave*>();
	_waveFactory[2] = boost::factory<CConvergentWave*>();
}
void Bootstrapper::RegisterExperimentTypes(){
	_experimentFactory[CBED] = boost::factory<CoherentCBED*>();
	_experimentFactory[TEM] = boost::factory<CoherentTEM*>();
	_experimentFactory[PTYCHO] = boost::factory<Ptychograph*>();
}
void Bootstrapper::RegisterPotentialTypes(){
	_potentialFactory["3DFFT"]= boost::factory<C3DFFTPotential*>();
	_potentialFactory["2DFFT"]= boost::factory<C2DFFTPotential*>();
	_potentialFactory["3D"]= boost::factory<C3DPotential*>();
	_potentialFactory["2D"]= boost::factory<C2DPotential*>();
	_potentialFactory["CUDA"]= boost::factory<CUDA2DPotential*>();
}

void Bootstrapper::RegisterDetectorTypes(){
	_detectorFactory[1]= boost::factory<FlatAreaDetector*>();
}

void Bootstrapper::RegisterStructureTypes(){

}

} /* namespace slicepp */
