/*
 * FlatAreaDetector.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: wenxuan
 */

#include "FlatAreaDetector.hpp"

namespace slicepp{
FlatAreaDetector::FlatAreaDetector(cDetectorConfPtr dc,PersistenceManagerPtr p):IDetector(dc,p) {
	_numSaved = 0;
}

FlatAreaDetector::~FlatAreaDetector() {
}

void FlatAreaDetector::RecordImage(af::array& wave){
	af::array image = wave;
	_p->Save2DDataSet(image, "wave before noise");
	if (_dc->MaxElectronCounts > FLT_EPSILON ){
		image = anscombeNoise(image);
	}
	_p->Save2DDataSet(image, "wave after noise");
	af::fft2InPlace(image);
	image = MultiplyMTF(image);
	af::ifft2InPlace(image);
	_p->Save2DDataSet(image, "wave after mtf");
	_numSaved++;
	image = af::pow2(af::abs(image));
	_p->Save2DDataSet(image, "Detected Image_" + std::to_string(_numSaved));
}
//https://en.wikipedia.org/wiki/Anscombe_transform#cite_note-3
af::array FlatAreaDetector::anscombeNoise(af::array& wave){
	af::array wavem = af::real(wave) * _dc->MaxElectronCounts;
//	af::array filter = wavem > 1e-2f;
	af::array y = af::randn(wave.dims(0), wave.dims(1)) * wavem;
	auto poisson = 0.25*af::pow2(y) + 0.306186 * af::pow(y,-1) - 11.0/8*af::pow(y,-2) + 0.765465 * af::pow(y,-3) - 0.125;
	poisson /= _dc->MaxElectronCounts;
//	x = x * filter * af::sqrt( 1 - af::exp( (-1) * wavem / 0.777134f)) + wavem * (!filter);
//	x += filter * 2 * af::sqrt( wavem + 0.375f ) - 0.25 / af::sqrt(wavem);
//	x = af::round( 0.25f * x * x  - 0.375f );
//	poisson = af::round(poisson);
	poisson(poisson < FLT_MIN) = 0;
	return af::complex(poisson, af::imag(wave));
}
af::array FlatAreaDetector::MultiplyMTF(af::array& wave){
	af::array ix, iy, temp, mtf;
	ix = af::seq(_dc->n[0]);
	iy = af::seq(_dc->n[1]);
	temp = ix > _dc->n[0]/2;
	ix -= temp * (_dc->n[0]/2);
	temp = iy > _dc->n[1]/2;
	iy -= temp * (_dc->n[1]/2);
	ix = af::tile(ix/_dc->n[0], 1, _dc->n[1]);
	iy = af::tile(iy.T()/_dc->n[1], _dc->n[0]);
	mtf = af::sqrt(ix *ix + iy*iy);
	mtf = (_dc->mtfA * af::exp((-1) * _dc->mtfC * mtf) + _dc->mtfB * af::exp((-1) * _dc->mtfD * mtf * mtf));
	ix *= PI;
	iy *= PI;
	mtf *= ((af::sin ( ix ) + FLT_EPSILON ) / ( ix + FLT_EPSILON ) ) * ( (af::sin ( iy ) + FLT_EPSILON ) / ( iy + FLT_EPSILON ) );
	_p->Save2DDataSet(mtf, "mtf");
	return wave * mtf;
}
}
