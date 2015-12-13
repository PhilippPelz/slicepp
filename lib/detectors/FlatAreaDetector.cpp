/*
 * FlatAreaDetector.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: wenxuan
 */

#include "FlatAreaDetector.hpp"

namespace QSTEM{
FlatAreaDetector::FlatAreaDetector(cDetectorConfPtr dc,PersistenceManagerPtr p):IDetector(dc,p) {
	_numSaved = 0;
}

FlatAreaDetector::~FlatAreaDetector() {
}

void FlatAreaDetector::RecordImage(WavePtr w){
	_image = af::array(w->GetWaveAF());
	float_tt scale = 1/(float_tt)(_dc->nx *_dc->ny);
	af::fft2(_image);
	if ( w->GetPixelDose() > FLT_EPSILON ){
		af::ifft2(_image);
		anscombeNoise(_image, w->GetPixelDose());
		af::fft2(_image);
	}
	MultiplyMTF(_image);
	af::ifft2(_image);
	_numSaved++;
	_p->Save2DDataSet(_image, "Detected Image_" + std::to_string(_numSaved));
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
	ix = af::seq(_dc->nx);
	iy = af::seq(_dc->ny);
	temp = ix > _dc->nx/2;
	ix -= temp * (_dc->nx/2);
	temp = iy > _dc->ny/2;
	iy -= temp * (_dc->ny/2);
	ix = af::tile(ix/_dc->nx, 1, _dc->ny);
	iy = af::tile(iy.T()/_dc->ny, _dc->nx);
	mtf = af::sqrt(ix *ix + iy*iy);
	mtf = (_dc->mtfA * af::exp((-1) * _dc->mtfC * mtf) + _dc->mtfB * af::exp((-1) * _dc->mtfC * mtf * mtf));
	ix *= PI;
	iy *= PI;
	mtf *= ((af::sin ( ix ) + FLT_EPSILON ) / ( ix + FLT_EPSILON ) ) * ( (af::sin ( iy ) + FLT_EPSILON ) / ( iy + FLT_EPSILON ) );
	wave *= mtf;
}
}
