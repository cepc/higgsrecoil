#!/usr/bin/env bash 

rm -fr build
mkdir build
cd build
#HFcmake  # avoid alias in the script. Need to update in the next version. 
cmake -C /afs/ihep.ac.cn/soft/common/gcc/v01-17-05/ILCSoft.cmake .. 
make install
cd ..


