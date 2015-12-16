/*
 QSTEM - image simulation for TEM/STEM/CBED
 Copyright (C) 2000-2010  Christoph Koch
 Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "configs.hpp"
#include <stdlib.h>
#include <string.h>

StructureConfig* StructureConfig_clone(const StructureConfig* cloneMe){
	StructureConfig* p = (StructureConfig*)malloc(sizeof(StructureConfig));
	memcpy((void*)p,(const void*) cloneMe,sizeof(StructureConfig));
	return p;
}

StructureConfig* StructureConfig_new(){
	StructureConfig* p =  (StructureConfig*)malloc(sizeof(StructureConfig));
	memset((void*)p,0,sizeof(StructureConfig));
	return p;
}
ModelConfig* ModelConfig_new(){
	ModelConfig* p =  (ModelConfig*)malloc(sizeof(ModelConfig));
	memset((void*)p,0,sizeof(ModelConfig));
	return p;
}
WaveConfig* WaveConfig_new(){
	WaveConfig* p =  (WaveConfig*)malloc(sizeof(WaveConfig));
	memset((void*)p,0,sizeof(WaveConfig));
	return p;
}
OutputConfig* OutputConfig_new(){
	OutputConfig* p =  (OutputConfig*)malloc(sizeof(OutputConfig));
	memset((void*)p,0,sizeof(OutputConfig));
	return p;
}
DetectorConfig* DetectorConfig_new(){
	DetectorConfig* p =  (DetectorConfig*)malloc(sizeof(DetectorConfig));
	memset((void*)p,0,sizeof(DetectorConfig));
	return p;
}
ScanConfig* ScanConfig_new(){
	ScanConfig* p =  (ScanConfig*)malloc(sizeof(ScanConfig));
	memset((void*)p,0,sizeof(ScanConfig));
	return p;
}
c_Config* c_Config_new(){
	c_Config* p =  (c_Config*)malloc(sizeof(c_Config));
	memset((void*)p,0,sizeof(c_Config));
	return p;
}

