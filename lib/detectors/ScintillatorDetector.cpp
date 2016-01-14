/*
 * FlatAreaDetector.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: wenxuan
 */

#include "ScintillatorDetector.hpp"
#include "afhelpers.hpp"

namespace slicepp{
ScintillatorDetector::ScintillatorDetector(cDetectorConfPtr dc,PersistenceManagerPtr p):IDetector(dc,p) {
	_numSaved = 0;
}

ScintillatorDetector::~ScintillatorDetector() {
}

void ScintillatorDetector::RecordImage(af::array& wave){
	af::array image = wave;
//	_p->Save2DDataSet(image, "wave before noise");
//	if (_dc->MaxElectronCounts > FLT_EPSILON ){
//		image = anscombeNoise(image);
//	}
//	_p->Save2DDataSet(image, "wave after noise");
	af::fft2InPlace(image);
	image = MultiplyMTF(image);
	af::ifft2InPlace(image);
//	_p->Save2DDataSet(image, "wave after mtf");
	image = af::pow2(af::abs(image));
	_p->SaveMeasurement(image,_numSaved);
	_numSaved++;
}
//https://en.wikipedia.org/wiki/Anscombe_transform#cite_note-3
af::array ScintillatorDetector::anscombeNoise(af::array& wave){
	af::array wavem = af::real(wave);
	af::sync();
	_p->Save2DDataSet(wavem, "wavem");
	wavem *= _dc->MaxElectronCounts;

//	af::array filter = wavem > 1e-2f;
	af::array y = af::randn(wave.dims(0), wave.dims(1)) * wavem;
	auto poisson = 0.25*af::pow2(y) + 0.306186 * af::pow(y,-1) - 11.0/8*af::pow(y,-2) + 0.765465 * af::pow(y,-3) - 0.125;
	_p->Save2DDataSet(poisson, "poisson");
	poisson /= _dc->MaxElectronCounts;
	poisson(poisson < FLT_MIN) = 0;
	return af::complex(poisson, af::imag(wave));
}
af::array ScintillatorDetector::MultiplyMTF(af::array& wave){
	auto ix = (af::range(af::dim4(_dc->n[0],_dc->n[1]),1) - _dc->n[0]/2)/_dc->n[0];
	auto iy = (af::range(af::dim4(_dc->n[0],_dc->n[1]),0) - _dc->n[1]/2)/_dc->n[1];
	ix = fftShift(ix);
	iy = fftShift(iy);
	auto mtf = af::sqrt(ix *ix + iy*iy);
	mtf = (_dc->mtfA * af::exp((-1) * _dc->mtfC * mtf) + _dc->mtfB * af::exp((-1) * _dc->mtfD * mtf * mtf));
	ix *= PI;
	iy *= PI;
	mtf *= ((af::sin ( ix ) + FLT_EPSILON ) / ( ix + FLT_EPSILON ) ) * ( (af::sin ( iy ) + FLT_EPSILON ) / ( iy + FLT_EPSILON ) );
//	_p->Save2DDataSet(mtf, "mtf");
	return wave * mtf.as(c32);
}
}
