/*
 * FlatAreaDetector.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: wenxuan
 */

#include "FlatAreaDetector.hpp"

namespace QSTEM{
FlatAreaDetector::FlatAreaDetector(const ConfigPtr& c, const PersistenceManagerPtr& persist):IDetector(c, persist) {
	_numSaved = 0;
}

FlatAreaDetector::~FlatAreaDetector() {
}

void FlatAreaDetector::RecordImage(WavePtr w){
	_nx = _c->Wave.nx;
	_ny = _c->Wave.ny;
	_image = af::array(w->GetWaveAF());
	float_tt scale = 1/(float_tt)(_nx *_ny);
	af::fft2(_image);
	if ( w->GetPixelDose() > FLT_EPSILON ){
		af::ifft2(_image);
		anscombeNoise(_image, w->GetPixelDose());
		af::fft2(_image);
	}
	MultiplyMTF(_image);
	af::ifft2(_image);
	_numSaved++;
	_persist->Save2DDataSet(_image, "Detected Image_" + std::to_string(_numSaved));
}

void FlatAreaDetector::anscombeNoise(af::array wave, float_tt dose){
	af::array wavem = af::real(wave) * dose;
	af::array filter = wavem > 1e-2f;
	af::array x = af::randn(wave.dims(0), wave.dims(1));
	x = x * filter * af::sqrt( 1 - af::exp( (-1) * wavem / 0.777134f)) + x * (!filter);
	x += filter * 2 * af::sqrt( wavem + 0.375f ) - 0.25 / af::sqrt( wavem);
	x = af::round( 0.25f * x * x  - 0.375f );
	filter = x < FLT_MIN;
	x *= filter;
	wave = af::complex(x/dose, af::imag(wave));
}
void FlatAreaDetector::MultiplyMTF(af::array wave ){
	af::array ix, iy, temp, mtf;
	ix = af::seq(_nx);
	iy = af::seq(_ny);
	temp = ix > _nx/2;
	ix -= temp * (_nx/2);
	temp = iy > _ny/2;
	iy -= temp * (_ny/2);
	ix = af::tile(ix/_nx, 1, _ny);
	iy = af::tile(iy.T()/_ny, _nx);
	mtf = af::sqrt(ix *ix + iy*iy);
	mtf = (_c->Detector.mtfA * af::exp((-1) * _c->Detector.mtfC * mtf) + _c->Detector.mtfB * af::exp((-1) * _c->Detector.mtfC * mtf * mtf));
	ix *= PI;
	iy *= PI;
	mtf *= ((af::sin ( ix ) + FLT_EPSILON ) / ( ix + FLT_EPSILON ) ) * ( (af::sin ( iy ) + FLT_EPSILON ) / ( iy + FLT_EPSILON ) );
	wave *= mtf;
}
}
