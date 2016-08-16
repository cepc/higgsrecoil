#!/usr/bin/env bash 

shopt -s expand_aliases
source /besfs/groups/higgs/Software/v01-17-05_slc6/init_ilcsoft.sh
rm -fr build
mkdir build
cd build
HFcmake
make install
cd ..


