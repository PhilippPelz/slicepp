/*
 * afhelpers.hpp
 *
 *  Created on: Jan 13, 2016
 *      Author: philipp
 */

#ifndef AFHELPERS_HPP_
#define AFHELPERS_HPP_
#include <arrayfire.h>
#include <af/util.h>

af::array fftShift(af::array& _wave);
af::array ifftShift(af::array& _wave);

#endif /* AFHELPERS_HPP_ */
