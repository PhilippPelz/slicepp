/*
 * afhelpers.cpp
 *
 *  Created on: Jan 13, 2016
 *      Author: philipp
 */

#include "afhelpers.hpp"

af::array fftShift(af::array& _wave){
	return af::shift(_wave, _wave.dims(0) / 2, _wave.dims(1) / 2);
}
af::array ifftShift(af::array& _wave){
	return af::shift(_wave, (_wave.dims(0) + 1) / 2, (_wave.dims(1) + 1) / 2);
}
