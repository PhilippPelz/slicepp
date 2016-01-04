/*
 * Scan.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: wenxuan
 */

#include "Scan.hpp"
namespace slicepp{
Scan::Scan(const ConfigPtr& c) {
	_xStart = c->Scan->xPos;
	_yStart = c->Scan->yPos;
	_xStep = c->Scan->xStep;
	_yStep = c->Scan->yStep;
	_probeX = c->Wave->nx;
	_probeY = c->Wave->ny;
	_scanType = c->Scan->scanType;
	_c = ConfigPtr(c);
}
Scan::~Scan() {
}
std::vector<std::pair<int, int>> Scan::GetScanPositions() {
	std::vector<std::pair<int, int>> positions =
			std::vector<std::pair<int, int>>();
	if (_scanType == 1) {
		for (int i = _xStart; i < _c->Model->n[0]; i += _xStep) {
			for (int j = _yStart; j < _c->Model->n[1]; j += _yStep) {
				if ((i + _probeX <= _c->Model->n[0]) && (j + _probeY <= _c->Model->n[1])) {
					positions.push_back(std::make_pair(i, j));
				}
			}
		}
	}
	return positions;
}
}



