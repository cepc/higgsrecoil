#!/usr/bin/env bash

# Main driver to submit jobs 
# Author SHI Xin <shixin@ihep.ac.cn>
# Created [2016-08-16 Tue 08:29] 

usage() {
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "0.1.1"    "Run on signal samples" 
    printf "\nDATE\n"
    printf "\n\t%-5s\n" "AUGUST 2016"     
}


if [[ $# -eq 0 ]]; then
    usage
fi


option=$1

case $option in 
    0.1.1) echo "Running on signal sample..."
	   unset MARLIN_DLL
	   export MARLIN_DLL=./lib/libHiggsRecoilMass.so
	   mkdir result 
	   Marlin steer/signal.steer 
	   ;;
   
esac

