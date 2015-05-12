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
#define VIB_IMAGE_TEST

#ifndef _WIN32
#define UNIX
#endif
/* #define USE_FFT_POT */
// for memory leak checking in windows.  Should not affect speed of release builds.
#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>

#ifdef _WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#include <string.h>
#include <omp.h>
#include <time.h>

#include "Bootstrapper.hpp"
#include "config_IO/read_qsc.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include <iostream>

using boost::property_tree::ptree;
using namespace QSTEM;
using boost::property_tree::info_parser::read_info;

void usage()
{
	printf("usage: stem [input file='stem.dat']\n\n");
}



int main(int argc, char *argv[])
{
	Bootstrapper b(argc,argv);
	b.Initialize();
	auto expt = b.GetExperiment();
	expt->Run();

//	std::cout << "Press key to finish";
//	char name[256];
//	std::cin.getline (name,256);

#if _DEBUG
	_CrtDumpMemoryLeaks();
#endif
	return 0;
}
