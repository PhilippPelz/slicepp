#!/bin/sh
rm CMakeCache.txt
cmake ../.. -DCMAKE_BUILD_TYPE=Release
make -j2
