rm edep-sim.tar.gz nusystematics.tar.gz nusyst_inputs.tar.gz DUNE_ND_GeoEff.tar.gz sim_inputs.tar.gz
rm -rf sim_inputs

tar -zcf edep-sim.tar.gz edep-sim
tar -zcf nusystematics.tar.gz nusystematics
tar -zcf nusyst_inputs.tar.gz nusyst_inputs
tar -zcf DUNE_ND_GeoEff.tar.gz DUNE_ND_GeoEff

mkdir sim_inputs
cp copy_dune_ndtf_flux sim_inputs
cp Messenger_production.xml sim_inputs
cp GNuMIFlux.xml sim_inputs
cp MPD_SPY_LAr.gdml sim_inputs
cp MPD_SPY_LAr.volArgonCubeActive.maxpl.xml sim_inputs
cp dune-nd.mac sim_inputs
cp fhicl.fcl sim_inputs
cp makeCAF sim_inputs
cp dumpTree.py sim_inputs

tar -zcf sim_inputs.tar.gz sim_inputs

TARDIR="/pnfs/dune/persistent/users/$USER/cafmaker/tarballs"
mkdir -p $TARDIR
rm -f ${TARDIR}/*.gz

cp edep-sim.tar.gz nusystematics.tar.gz nusyst_inputs.tar.gz DUNE_ND_GeoEff.tar.gz sim_inputs.tar.gz ${TARDIR}

