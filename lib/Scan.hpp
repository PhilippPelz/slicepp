/*
 * Scan.hpp
 *
 *  Created on: Jul 22, 2015
 *      Author: wenxuan
 */

#ifndef SCAN_HPP_
#define SCAN_HPP_
#include <tuple>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstdlib>
#include "experiments/base.hpp"

namespace slicepp{
class Scan{
public:
	Scan(const ConfigPtr& c);
	std::vector<std::pair<int, int>> GetScanPositions();
	~Scan();

protected:
	int _xStart, _yStart, _xStep, _yStep, _scanType, _probeX, _probeY;
	ConfigPtr _c;

};
typedef boost::shared_ptr<Scan> ScanPtr;
}
#endif /* SCAN_HPP_ */
