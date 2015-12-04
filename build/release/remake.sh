#!/bin/sh
rm CMakeCache.txt
cmake ../.. -DCMAKE_BUILD_TYPE=Release -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda
make -j10
