#ifndef EXPERIMENT_PTYC_H
#define EXPERIMENT_PTYC_H

#include "base.hpp"
#include "Scan.hpp"
namespace slicepp
{
class Ptychograph : public BaseExperiment
{
public:
	Ptychograph(const ConfigPtr& c,const StructureBuilderPtr& s,const WavePtr& w,const PotPtr& p, const DetPtr& d, const PersistenceManagerPtr& pers);
    void Run();

    virtual void DisplayParams();
    virtual void WriteBeams(unsigned absoluteSlice);
    virtual ~Ptychograph(){};

protected:
    void CollectIntensity(unsigned absoluteSlice);
    void PostSliceProcess();
    void SaveImages();
    unsigned _nbout = 1;
    float_tt **_pendelloesung;
    ScanPtr _scan;
    bool _lbeams;			/* flag indicating whether to record beams */
    int scanx, scany;     /* The beam position on the sample */
    bool _showProbe;            /* if true, saves a plot of the probe */
    bool m_storeSeries;
};
}
#endif


