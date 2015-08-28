#!/bin/sh
cd ./doxyfiles
doxygen Doxyfile_Config
doxygen Doxyfile_Data
doxygen Doxyfile_Detector
doxygen Doxyfile_Structure
doxygen Doxyfile_Potential
doxygen Doxyfile_Wave
doxygen Doxyfile_Experiment
doxygen Doxyfile_Overview
cd ..
make html 
