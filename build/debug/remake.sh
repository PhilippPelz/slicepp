#!/bin/sh
rm CMakeCache.txt
rm -rf CMakeFiles
cmake ../.. -DCMAKE_BUILD_TYPE=Debug -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DArrayFire_DIR=/opt/arrayfiredebug
make -j8
#/usr/local/share/ArrayFire/cmake
