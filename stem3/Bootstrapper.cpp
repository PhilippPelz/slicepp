/*
 * Bootstrapper.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: philipp
 */

#include "Bootstrapper.hpp"

namespace QSTEM {

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
	std::string fileName;
	if (argc < 2)   fileName = "config.json";
	else    fileName=argv[1];

	ptree pt;
	bpt::json_parser::read_json(fileName,pt);
//	read_info(fileName,pt);
//	bpt::json_parser::write_json(std::string("config.json"),pt);
	_c = ConfigPtr(new Config(pt));

	logging::core::get()->set_filter
	(
	    logging::trivial::severity >= static_cast<logging::trivial::severity_level>(_c->Output.LogLevel)
	);
	logging::add_console_log(std::cout, keywords::format = ">> %Message%");
	if(_c->Output.WriteLogFile)
		logging::add_file_log
			(
				keywords::file_name = _c->Output.LogFileName.string(),
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
}

Bootstrapper::~Bootstrapper() {
	// TODO Auto-generated destructor stub
}

void Bootstrapper::Initialize(){
	fftw_init_threads();
	fftw_plan_with_nthreads(_c->nThreads);
	fftwpp::fftw::maxthreads = _c->nThreads;
	omp_set_num_threads(_c->nThreads);

	RegisterWaveTypes();
	RegisterPotentialTypes();
	RegisterStructureTypes();
	RegisterExperimentTypes();
	RegisterStructureReaders();
	RegisterStructureBuilders();
	boost::filesystem::path p(_c->Structure.structureFilename);
	auto sreader = StructureReaderPtr(_structureReaderFactory[".cif"](p));
	auto structureBuilder = StructureBuilderPtr(
			_structureBuilderFactory[p.extension().string()](sreader,_c));
//	structureBuilder->SetSliceThickness(_c->Model);
//	structureBuilder->SetResolution(_c->Model,_c->Potential);
	auto persist = PersistenceManagerPtr(new PersistenceManager(_c));

	std::stringstream str;
	str << (_c->Potential.Use3D ? "3D" : "2D");
	str << (_c->Potential.UseFFT ? "FFT" : "");

	auto wave = WavePtr(_waveFactory[_c->Wave.type](_c,persist));
	auto potential = PotPtr(_potentialFactory[str.str()](_c,persist));
	_e = ExperimentPtr( _experimentFactory[_c->ExperimentType](_c,structureBuilder,wave,potential,persist));
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
	_experimentFactory[ExperimentType::CBED] = boost::factory<CExperimentCBED*>();
	_experimentFactory[ExperimentType::STEM] = boost::factory<CExperimentSTEM*>();
	_experimentFactory[ExperimentType::TEM] = boost::factory<CExperimentTEM*>();
}
void Bootstrapper::RegisterPotentialTypes(){
	_potentialFactory["3DFFT"]= boost::factory<C3DFFTPotential*>();
	_potentialFactory["2DFFT"]= boost::factory<C2DFFTPotential*>();
	_potentialFactory["3D"]= boost::factory<C3DPotential*>();
	_potentialFactory["2D"]= boost::factory<C2DPotential*>();
}
void Bootstrapper::RegisterStructureTypes(){

}

} /* namespace QSTEM */
