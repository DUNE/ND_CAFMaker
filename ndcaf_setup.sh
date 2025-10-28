#!/bin/bash

BUILD_QUAL="prof"
if [ "$1" = "debug" ]; then
	BUILD_QUAL="debug"
elif [ "$1" = "prof" ]; then
	true
else
	echo
	echo "WARNING: build qualifier was unspecified (options: 'prof' or 'debug').  Assuming 'prof' by default!"
	echo
fi

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# some dependencies use CMake to build
setup cmake v3_22_2

# used for downloading packges (see ndcaf_setup.sh)
setup pycurl
setup curl

# used for building the ND_GeomEff dependency
setup eigen v23_08_01_66e8f 

# note in the following that the "eXX" qualifier
# implies a particular compiler and GCC toolchain,
# and needs to be the same throughout where used

# edep-sim has both ROOT and Geant4 dependencies.
# the upstream ND_Production containers are currently setting up
# a customized edep-sim version, but it boils down to 3.2.0
setup edepsim v3_2_0e -q e26:${BUILD_QUAL}

# The base GENIE version is 3.04.00
setup dk2nugenie   v01_10_01p -q e26:${BUILD_QUAL}
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup genie_phyopt v3_04_00 -q dkcharmtau

# duneanaobj contains the 'StandardRecord' specification of the CAF format
setup duneanaobj v03_11_00 -q e26:${BUILD_QUAL}

# direct dependencies of CAFMaker itself
setup hdf5 v1_12_2b  -q e26:${BUILD_QUAL}
setup fhiclcpp v4_18_04 -q e26:${BUILD_QUAL}

export LD_LIBRARY_PATH=$CURL_ROOT/lib:$LD_LIBRARY_PATH

#edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# shut up ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

# nusystematics paths
#export NUSYST=${PWD}/nusystematics
#export LD_LIBRARY_PATH=${NUSYST}/build/Linux/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=${NUSYST}/build/nusystematics/artless:$LD_LIBRARY_PATH
#export FHICL_FILE_PATH=${NUSYST}/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${PWD}/DUNE_ND_GeoEff/lib/

# duneananobj needs to be in the libs too
export LD_LIBRARY_PATH=${DUNEANAOBJ_LIB}:$LD_LIBRARY_PATH
# finally, add our lib & bin to the paths
mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$mydir/lib:$LD_LIBRARY_PATH
export PATH=$mydir/bin:$PATH

# our FCL needs to be findable too
export FHICL_FILE_PATH="$FHICL_FILE_PATH:$mydir/cfg"
