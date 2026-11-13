WORKSPACE="$PWD" apptainer shell --bind /cvmfs --bind "$WORKSPACE:$WORKSPACE" /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest



cd /path/to/ND_CAFMaker
source ndcaf_setup_deps.sh prof

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
cmake --install build --prefix install
source install/bin/ndcaf_setup.sh prof
command -v makeCAF || echo 'makeCAF not in PATH'
makeCAF --help




makeCAF --fcl=cafmaker_0000965.fcl -n 1

/exp/dune/app/users/abooth/Postdoc/nd-sim-reco/ND_Production/run-cafmaker

 scp -r /path/to/yourfolder/. user@login:/path/known

 scp -r [folder] [user]@perlmutter.nersc.gov:/global/cfs/cdirs/dune/www/data/emir_test/ND_CAFMaker_CI