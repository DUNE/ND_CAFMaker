#!/bin/bash

TOPDIR=${PWD}

FORCE=no
if [[ $# == 1 && x$1 == x-f ]]; then FORCE=yes; fi

# set up software
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cmake v3_22_2
setup gcc v9_3_0
setup pycurl
setup ifdhc
setup geant4 v4_11_0_p01c -q e20:debug
setup dk2nugenie   v01_10_01k -q debug:e20
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup genie_phyopt v3_04_00 -q dkcharmtau
setup jobsub_client
setup eigen v3_3_5
setup duneanaobj v03_06_01b -q e20:prof
setup srproxy v00.44 -q py3913
setup hdf5 v1_12_0b -q e20:prof
setup h5cpp v1_10_4_6c
setup fhiclcpp v4_15_03 -q debug:e20
setup edepsim v3_2_0c -q debug:e20

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# Just use the edep-sim UPS product, don't clone master branch off repos!
# Get edep-sim and build it
#if [ $FORCE == yes ]; then rm -rf edep-sim; fi
#git clone https://github.com/ClarkMcGrew/edep-sim.git
#cd edep-sim
#sed -i 's/add_definitions(-DEDEPSIM_FORCE_PRIVATE_FIELDS)//g' */CMakeLists.txt
#source setup.sh
#source build/edep-build.sh

cd ${TOPDIR}

## Update Jan 5 2022:
## DIRT-II task force is currently working on deciding on GENIE 3 models to work with for the ND TDR.
## Since nusystematics' internal fhicl-cpp replacement conflicts with the official fhicl-cpp UPS product,
## for the moment, we're disabling the nusystematics interface.
## It will be re-enabled and cleaned up with help from the DIRT-II group once the model & uncertainties are settled.

# Get nusystematics and built it "artless"
# The ART-dependent version can be ups setup but for ND we need this special build
#
# v00_04_01 is the last version to depend on genie v2 (specifically v2_12_10d)
#if [ $FORCE == yes ]; then rm -rf nusystematics; fi
#git clone ssh://p-nusystematics@cdcvs.fnal.gov/cvs/projects/nusystematics -b v00_04_01
#mkdir nusystematics/build
#cd nusystematics/build
#cmake ../ -DUSEART=0 -DLIBXML2_LIB=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/lib/ -DLIBXML2_INC=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/include/libxml2 -DPYTHIA6=/cvmfs/larsoft.opensciencegrid.org/products/pythia/v6_4_28i/Linux64bit+2.6-2.12-gcc640-prof/lib
#make -j systematicstools # force this to build first
#make -j
#make -j install

cd ${TOPDIR}

# Get Geometric efficiency library and build
if [ $FORCE == yes ]; then rm -rf DUNE_ND_GeoEff; fi
git clone --recurse-submodules https://github.com/cvilelasbu/DUNE_ND_GeoEff.git
cd DUNE_ND_GeoEff
cmake -DPYTHON_EXECUTABLE:FILEPATH=`which python` .
make -j pyGeoEff

cd ${TOPDIR}

# Add nusystematics to the paths
#export LD_LIBRARY_PATH=${TOPDIR}/nusystematics/build/Linux/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=${TOPDIR}/nusystematics/build/nusystematics/artless:$LD_LIBRARY_PATH
#export FHICL_FILE_PATH=${TOPDIR}/nusystematics/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${TOPDIR}/DUNE_ND_GeoEff/lib/

# make tarballs of edep-sim and nusystematics for grid jobs
#tar -zcf edep-sim.tar.gz edep-sim
#tar -zcf nusystematics.tar.gz nusystematics
#tar -zcf DUNE_ND_GeoEff.tar.gz DUNE_ND_GeoEff
