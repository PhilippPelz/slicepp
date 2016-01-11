#ifndef EXPERIMENT_PTYC_H
#define EXPERIMENT_PTYC_H

#include "base.hpp"
namespace slicepp {
class Ptychograph: public BaseExperiment {
public:
	Ptychograph(const ConfigPtr& c, const StructureBuilderPtr& s, const WavePtr& w, const PotPtr& p, const DetPtr& d,
			const PersistenceManagerPtr& pers);
	void Run();

	virtual void DisplayParams();
	virtual void WriteBeams(unsigned absoluteSlice);
	virtual ~Ptychograph() {
	}
	;
	void PostSpecimenProcess();

protected:
	void CollectIntensity(unsigned absoluteSlice);
	void PostSliceProcess();
	void SaveImages();
};
}
#endif

