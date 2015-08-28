/*
 * ScannedParameter.cpp
 *
 *  Created on: Aug 26, 2015
 *      Author: wenxuan
 */

#include "ScannedParameter.hpp"

namespace QSTEM {
ScannedParameter::ScannedParameter(float_tt start, float_tt stop, int nSteps){
	if (nSteps != 0){
		_start = start;
		_stop = stop;
		_nSteps = nSteps;
		_curStep = 0;
		_value = _start;
		_stepSize = (_stop - _start)/nSteps;
	}else{
		_start = start;
		_stop = start;
		_nSteps = 0;
		_curStep = 0;
		_value = _start;
		_stepSize = 0;
	}
}
ScannedParameter::ScannedParameter(float_tt start) {
	_start = start;
	_stop = start;
	_nSteps = 0;
	_curStep = 0;
	_value = _start;
	_stepSize = 0;
}

ScannedParameter::ScannedParameter() {
}

ScannedParameter::~ScannedParameter() {
}

void ScannedParameter::Next(){
	_curStep += 1;
	_value += _stepSize;
}

bool ScannedParameter::maxValue(){
	return (_curStep == _nSteps);
}

float_tt  ScannedParameter::operator()(){
	return _value;
}

void  ScannedParameter::Reset(){
	_value = _start;
	_curStep = 0;
}

} /* namespace QSTEM */
