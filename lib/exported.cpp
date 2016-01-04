/*
 * exported.cpp
 *
 *  Created on: Dec 16, 2015
 *      Author: philipp
 */

#include "exported.hpp"
#include "Bootstrapper.hpp"

void run_simulation(c_Config* conf){
    slicepp::Bootstrapper b;
    b.Initialize(conf);
    auto e = b.GetExperiment();
    e->Run();
}
