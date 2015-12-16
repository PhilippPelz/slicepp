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
#ifndef _WIN32
#define UNIX
#endif
#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>

#ifdef _WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#include "Bootstrapper.hpp"
#include <string.h>
#include <omp.h>
#include <time.h>
#include <iostream>

using boost::property_tree::ptree;
using namespace QSTEM;
using boost::property_tree::info_parser::read_info;

int main(int argc, char *argv[])
{
	Bootstrapper b(argc,argv);
	b.Initialize();
	auto expt = b.GetExperiment();
	expt->Run();

#if _DEBUG
	_CrtDumpMemoryLeaks();
#endif
	return 0;
}
