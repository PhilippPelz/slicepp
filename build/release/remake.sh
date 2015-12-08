#!/bin/sh
rm CMakeCache.txt
rm -rf CMakeFiles
cmake ../.. -DCMAKE_BUILD_TYPE=Release -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DArrayFire_DIR=/usr/local/share/ArrayFire/cmake
make -j8
