# ND_CAFMaker
Code for running ND parameterized reconstruction and making CAFs

source build.sh in this directory initially.
This will get external dependencies edep-sim and nusystematics and compile them locally. 

Set up software at fermilab with ndcaf_setup.sh.

Running is done with run_everything.sh which can be called with
jobsub_submit  --group dune --role=Analysis -N 100 --OS=SL7 --expected-lifetime=12h --memory=2000MB file://run_everything.sh HORN FIRSTRUN POT OFFAXIS FLUX
where HORN is FHC or RHC
FIRSTRUN is the first run number you want
POT is the POT per job, suggest ~1e16 for a ~6hr job using dk2nu
OFFAXIS is the off-axis position of the ND, in meters
FLUX is either gsimple or dk2nu; gsimple is faster for on-axis running, but dk2nu is recommended for off-axis

This script will get tarballs of pre-compiled binaries from the /pnfs/dune/persistent/users/LBL_TDR area. 
If you want to update the code, it is recommended that you make new tarballs, put them in /pnfs, and modify 
the grid script to use your new ones.
