#! /usr/bin/env bash

HORN=$1
FIRST=$2
NPOT=$3
OFFAXIS=$4
FLUX=$5
TEST=$6

TIME_START=`date +%s`

if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

MODE="neutrino"
RHC=""
if [ "${HORN}" = "RHC" ]; then
MODE="antineutrino"
RHC=" --rhc"
fi

if [ "${NPOT}" = "" ]; then
echo "POT not specified, using 1E16"
NPOT=1E16
fi

if [ "${FIRST}" = "" ]; then
echo "First run number not specified, using 0"
FIRST=0
fi

if [ "${OFFAXIS}" = "" ]; then
echo "Off axis position not specified, assuming on-axis"
OFFAXIS=0
fi

if [ "${FLUX}" = "" ]; then
echo "Flux not specified, using dk2nu"
FLUX="dk2nu"
fi

if [ "${FLUX}" = "dk2nu" ]; then
FLUXOPT="--dk2nu"
FLUXDIR="/pnfs/dune/persistent/users/ljf26/fluxfiles/g4lbne/v3r5p4/QGSP_BERT"
fi
if [ "${FLUX}" = "gsimple" ]; then
FLUXOPT=""
FLUXDIR="/pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d"
fi

CP="ifdh cp"
if [ "${TEST}" = "test" ]; then
echo "Test mode"
PROCESS=0
mkdir -p test
cd test
fi

# ifdhc doen't have a mkdir -p equivalent, which is fine 
# as long as you always remember to include this convenient function in your scripts
ifdh_mkdir_p() {
    local dir=$1
    local force=$2
    if [ `ifdh ls $dir 0 $force | wc -l` -gt 0 ] 
    then
        : # we're done
    else
        ifdh_mkdir_p `dirname $dir` $force
        ifdh mkdir $dir $force
    fi
}

OADIR="${OFFAXIS}m"
if [ "${HORN}" = "RHC" ]; then
OADIR="${OFFAXIS}mRHC"
fi

RUNNO=$((${PROCESS}+${FIRST}))
RNDSEED=$((1000000*${OFFAXIS}+${RUNNO}+1000000))

# 5E16 is about 15000 events on-axis which runs in ~6 hours
NEVENTS="-e ${NPOT}"      # -n XXXX number of events, -e XE16 for POT

RDIR=$((${RUNNO} / 1000))
if [ ${RUNNO} -lt 10000 ]; then
RDIR=0$((${RUNNO} / 1000))
fi

GEOMETRY="MPD_SPY_LAr"
TOPVOL="volArgonCubeActive"

TARDIR="/pnfs/dune/persistent/users/LBL_TDR/sw_tarballs"
OUTDIR="/pnfs/dune/persistent/users/marshalc/nd_offaxis/v7"

# Don't try over and over again to copy a file when it isn't going to work
export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

##################################################

# Setup UPS and required products
echo "Setting up software"
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup ifdhc
setup dk2nugenie   v01_06_01f -q debug:e15
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup geant4 v4_10_3_p01b -q e15:prof

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

##################################################
# Get the binaries & other files that are needed
${CP} ${TARDIR}/sim_inputs.tar.gz sim_inputs.tar.gz
${CP} ${TARDIR}/edep-sim.tar.gz edep-sim.tar.gz
${CP} ${TARDIR}/nusystematics.tar.gz nusystematics.tar.gz
${CP} ${TARDIR}/nusyst_inputs.tar.gz nusyst_inputs.tar.gz
${CP} ${TARDIR}/DUNE_ND_GeoEff.tar.gz DUNE_ND_GeoEff.tar.gz


tar -xzf sim_inputs.tar.gz
tar -xzf edep-sim.tar.gz
tar -xzf nusystematics.tar.gz
tar -xzf nusyst_inputs.tar.gz
tar -xzf DUNE_ND_GeoEff.tar.gz
mv sim_inputs/* ${PWD}

# Get flux files to local node
# dk2nu files: /pnfs/dune/persistent/users/ljf26/fluxfiles/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/flux/dk2nu
# gsimple files: /pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d/OptimizedEngineeredNov2017/neutrino/
chmod +x copy_dune_ndtf_flux
./copy_dune_ndtf_flux --top ${FLUXDIR} --outpath local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=300 ${FLUXOPT}

# GENIE for some reason doesn't recognize *.dk2nu.root as dk2nu format, but it works if dk2nu is at the front?
if [ "${FLUX}" = "dk2nu" ]; then
cd local_flux_files
for f in *.dk2nu.root
do
    mv "$f" "dk2nu_$f"
done
cd ..
fi

##################################################

# Modify GNuMIFlux.xml to the specified off-axis position
sed -i "s/<beampos> ( 0.0, 0.05387, 6.66 )/<beampos> ( ${OFFAXIS}, 0.05387, 6.66 )/g" GNuMIFlux.xml

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

# Run GENIE
echo "Running gevgen"
TIME_GENIE=`date +%s`
gevgen_fnal \
    -f local_flux_files/${FLUX}*,DUNEND \
    -g ${GEOMETRY}.gdml \
    -t ${TOPVOL} \
    -L cm -D g_cm3 \
    ${NEVENTS} \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC
    #-m ${GEOMETRY}.${TOPVOL}.maxpl.xml \

##################################################

# Convert the genie output to rootracker

export LD_LIBRARY_PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/lib:${LD_LIBRARY_PATH}
export PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/bin:${PATH}

echo "Running gntpc"
TIME_ROOTRACKER=`date +%s`
${CP} ${MODE}.${RNDSEED}.ghep.root input_file.ghep.root
gntpc -i input_file.ghep.root -f rootracker \
      --event-record-print-level 0 \
      --message-thresholds Messenger_production.xml

# edep-sim wants number of events, but we are doing POT so the files will be slightly different
# get the number of events from the GENIE files to pass it to edep-sim
NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | genie -l -b input_file.ghep.root 2>/dev/null  | tail -1)

##################################################

# Run edep-sim
echo "Running edep-sim with ${NPER} events."
TIME_EDEPSIM=`date +%s`
edep-sim \
    -C \
    -g ${GEOMETRY}.gdml \
    -o edep.${RNDSEED}.root \
    -e ${NPER} \
    dune-nd.mac

##################################################

# The MakeProject in dumpTree.py won't work if edep-sim is in the library path for reasons unknown
# This is a hack, please avert your eyes if you don't want to see my garbage
unset LD_LIBRARY_PATH
setup ifdhc
setup dk2nugenie   v01_06_01f -q debug:e15
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup geant4 v4_10_3_p01b -q e15:prof

# Add nusystematics to the paths
export LD_LIBRARY_PATH=${PWD}/nusystematics/build/Linux/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${PWD}/nusystematics/build/nusystematics/artless:${LD_LIBRARY_PATH}
export FHICL_FILE_PATH=${PWD}/nusystematics/nusystematics/fcl:${FHICL_FILE_PATH}

# add pyGeoEff to pythonpath, and libgeoEff to LD_LIBRARY_PATH
export PYTHONPATH=${PWD}/DUNE_ND_GeoEff/lib/:${PYTHONPATH}
export LD_LIBRARY_PATH=${PWD}/DUNE_ND_GeoEff/lib:${LD_LIBRARY_PATH}

# Run dumpTree to make a root file, you can start reading again if you averted your eyes before
TIME_DUMPTREE=`date +%s`
python dumpTree.py --infile edep.${RNDSEED}.root --outfile ${HORN}.${RNDSEED}.root --seed ${RNDSEED}

# Run CAFMaker
TIME_CAFMAKER=`date +%s`
./makeCAF --infile ${HORN}.${RNDSEED}.root --gfile ${MODE}.${RNDSEED}.ghep.root --outfile ${HORN}.${RNDSEED}.CAF.root --fhicl ./fhicl.fcl --seed ${RNDSEED} ${RHC} --oa ${OFFAXIS}

##################################################
# Copy the output files
echo "It's copy time, here are the files that I have:"
TIME_COPY=`date +%s`

ifdh_mkdir_p ${OUTDIR}/genie/${OADIR}/${RDIR}
ifdh_mkdir_p ${OUTDIR}/dump/${OADIR}/${RDIR}
ifdh_mkdir_p ${OUTDIR}/CAF/${OADIR}/${RDIR}
ifdh_mkdir_p ${OUTDIR}/edep/${OADIR}/${RDIR}

# GENIE, this is usually small and good idea to save
${CP} ${MODE}.${RNDSEED}.ghep.root ${OUTDIR}/genie/${OADIR}/${RDIR}/${HORN}.${RNDSEED}.ghep.root

# G4/edep-sim file is HUGE and probably we can't save it
#${CP} edep.${RNDSEED}.root ${OUTDIR}/edep/${OADIR}/${RDIR}/${HORN}.${RNDSEED}.edep.root

# "dump" file is useful for various analyses
${CP} ${HORN}.${RNDSEED}.root ${OUTDIR}/dump/${OADIR}/${RDIR}/${HORN}.${RNDSEED}.dump.root

${CP} ${HORN}.${RNDSEED}.CAF.root ${OUTDIR}/CAF/${OADIR}/${RDIR}/${HORN}.${RNDSEED}.CAF.root

TIME_END=`date +%s`
# Print out a single thing that says the time of each step
TIME_S=$((${TIME_GENIE}-${TIME_START}))
TIME_G=$((${TIME_ROOTRACKER}-${TIME_GENIE}))
TIME_R=$((${TIME_EDEPSIM}-${TIME_ROOTRACKER}))
TIME_E=$((${TIME_DUMPTREE}-${TIME_EDEPSIM}))
TIME_D=$((${TIME_CAFMAKER}-${TIME_DUMPTREE}))
TIME_M=$((${TIME_COPY}-${TIME_CAFMAKER}))
TIME_C=$((${TIME_END}-${TIME_COPY}))
echo "Start-up time: ${TIME_S}"
echo "gevgen time: ${TIME_G}"
echo "gntpc time: ${TIME_R}"
echo "edep-sim time: ${TIME_E}"
echo "dumpTree time: ${TIME_D}"
echo "makeCAF time: ${TIME_M}"
echo "Copy time: ${TIME_C}"
