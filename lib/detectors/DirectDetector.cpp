/*
 * DirectDetector.cpp
 *
 *  Created on: Jan 14, 2016
 *      Author: philipp
 */

#include "DirectDetector.hpp"

namespace slicepp {

DirectDetector::DirectDetector(cDetectorConfPtr dc,PersistenceManagerPtr p):IDetector(dc,p) {
	// TODO Auto-generated constructor stub

}
void DirectDetector::RecordImage(af::array& w){
	auto abs = af::abs(w);
	auto I = abs*abs;
//	float_tt max = af::max<float_tt>(I);
//	I *= _dc->MaxElectronCounts / max;
//	I = af::round(Noise(I));
//	I = af::round((I));
//	_p->Save2DDataSet(I, "Iround");
	_p->SaveMeasurement(I,_numSaved);
	_numSaved++;
}
af::array DirectDetector::Noise(af::array& I){
//	y = std.*randn(1000,1) + mean;
	af::array r = af::randn(I.dims());
//	_p->Save2DDataSet(r, "r");
//	_p->Save2DDataSet(I, "I");
	r = af::sqrt(I) * r + I;
//	_p->Save2DDataSet(r, "gaussian");
//	auto poisson = 0.25*af::pow2(r) + 0.306186 * af::pow(r,-1) - 11.0/8*af::pow(r,-2) + 0.765465 * af::pow(r,-3) - 0.125;
//	_p->Save2DDataSet(poisson, "poisson");
//	poisson /= _dc->MaxElectronCounts;
//	r(r < FLT_MIN) = 0;
	return r;
}
DirectDetector::~DirectDetector() {
	// TODO Auto-generated destructor stub
}

} /* namespace slicepp */
