/*
 * ScannedParameter.hpp
 *
 *  Created on: Aug 26, 2015
 *      Author: wenxuan
 */
#include "stemtypes_fftw3.hpp"
#include <complex>
#include <float.h>
#ifndef ScannedParameter_HPP_
#define ScannedParameter_HPP_

namespace QSTEM {

class ScannedParameter {
public:
	ScannedParameter(float_tt value);
	ScannedParameter(float_tt start, float_tt stop, int nSteps);
	void Next();
	bool maxValue();
	void Reset();
	float_tt operator()();
	ScannedParameter();
	virtual ~ScannedParameter();
	float_tt _start;
	float_tt _stop;
	float_tt _stepSize;
	float_tt _value;
	int _nSteps;
	int _curStep;

};

} /* namespace QSTEM */
#endif /* ScannedParameter_HPP_ */
