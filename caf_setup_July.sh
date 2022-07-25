

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

#recalling enviromental variables

export MAIN_FOLDER=${PWD}
export CUSTOM_UPS="/dune/app/users/$USER/ups"
export PRODUCTS="$PRODUCTS:$CUSTOM_UPS"
DUNEANAOBJ_SRC=/dune/app/users/ldinoto/dunesoft2/duneanaobj/
DUNEANAOBJ_BUILD=/dune/app/users/ldinoto/duneanaobj-build2/

setup cmake v3_9_0
setup gcc v6_4_0
setup pycurl
setup ifdhc
setup dk2nugenie   v01_06_01f -q debug:e15
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup geant4 v4_10_3_p01b -q e15:debug
setup jobsub_client
setup eigen v3_3_5
#setup duneanaobj v01_01_01 -q e15:gv1:debug
setup hdf5 v1_10_2a -q e15
setup fhiclcpp v4_06_08 -q debug:e15

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# shut up ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

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

echo "I keep you to the duneanaobj folder for building our own duneanaobj"

cd $DUNEANAOBJ_BUILD
unsetup cmake
setup cmake v3_17_3
setup root v6_12_06a -q e15:debug

. $DUNEANAOBJ_SRC/ups/setup_for_development -d e15:gv1
buildtool -I "$CUSTOM_UPS" -bti

cp /dune/app/users/jwolcott/dunesoft/duneanaobj/ups/duneanaobj.table $CUSTOM_UPS/duneanaobj/v01_01_01/ups/

setup duneanaobj v01_01_01 -q debug:e15:gv1

cd $MAIN_FOLDER
unsetup cmake
setup cmake v3_9_0


