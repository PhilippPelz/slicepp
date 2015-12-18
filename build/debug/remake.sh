#!/bin/sh
rm CMakeCache.txt
rm -rf CMakeFiles
rm -rf bin
rm -rf lib
rm -rf exe
rm -rf hdfboost
rm -rf sliceth

OLDPATH=$PATH
if [[ $(echo $PATH | grep anaconda) ]]; then
    export PATH=$(echo $PATH | tr ':' '\n' | grep -v "anaconda/bin" | grep -v "anaconda/lib" | grep -v "anaconda/include" | uniq | tr '\n' ':')
fi

cmake ../.. -DCMAKE_BUILD_TYPE=Debug -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DArrayFire_DIR=/usr/local/share/ArrayFire/cmake -DHDF5_ROOT=/usr/local
make -j8

export PATH=$OLDPATH
