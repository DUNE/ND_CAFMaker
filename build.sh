TOPDIR=${PWD}

# set up software
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cmake v3_9_0
setup gcc v6_4_0
setup pycurl
setup ifdhc
setup dk2nugenie   v01_06_01f -q debug:e15
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup geant4 v4_10_3_p01b -q e15:prof
setup jobsub_client
setup eigen v3_3_5
setup python v2_7_3

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# Get edep-sim and build it
git clone https://github.com/ClarkMcGrew/edep-sim.git
cd edep-sim
source setup.sh
source build/edep-build.sh

cd ${TOPDIR}

# Get nusystematics and built it "artless"
# The ART-dependent version can be ups setup but for ND we need this special build
git clone ssh://p-nusystematics@cdcvs.fnal.gov/cvs/projects/nusystematics
mkdir nusystematics/build
cd nusystematics/build
cmake ../ -DUSEART=0 -DLIBXML2_LIB=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/lib/ -DLIBXML2_INC=/cvmfs/larsoft.opensciencegrid.org/products/libxml2/v2_9_5/Linux64bit+2.6-2.12-prof/include/libxml2 -DPYTHIA6=/cvmfs/larsoft.opensciencegrid.org/products/pythia/v6_4_28i/Linux64bit+2.6-2.12-gcc640-prof/lib
make systematicstools # force this to build first
make
make install

cd ${TOPDIR}

# Get Geometric efficiency library and build
git clone --recurse-submodules https://github.com/cvilelasbu/DUNE_ND_GeoEff.git
cd DUNE_ND_GeoEff
cmake -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_DIR/bin/python .
make

cd ${TOPDIR}

# Add nusystematics to the paths
export LD_LIBRARY_PATH=${TOPDIR}/nusystematics/build/Linux/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${TOPDIR}/nusystematics/build/nusystematics/artless:$LD_LIBRARY_PATH
export FHICL_FILE_PATH=${TOPDIR}/nusystematics/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${TOPDIR}/DUNE_ND_GeoEff/lib/

# make tarballs of edep-sim and nusystematics for grid jobs
tar -zcf edep-sim.tar.gz edep-sim
tar -zcf nusystematics.tar.gz nusystematics
tar -zcf DUNE_ND_GeoEff.tar.gz DUNE_ND_GeoEff
