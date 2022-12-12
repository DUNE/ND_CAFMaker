# ND CAFMaker
`ND_CAFMaker` takes input `edep-sim`, `GENIE`, and reconstructed objects from the DUNE ND and combines them into the Common Analysis Format ("CAF").

# Setup
The package has a number of dependencies, all accessible through the `ups` framework at Fermilab and github. Simply run
```
./build.sh
```
to build the dependencies of `DUNE_ND_GeoEff` and `edep-sim`. Then finally
```
source ndcaf_setup.sh
```
and your environment should be set up. Depending on your role (developer or user), you may have to set up different `duneanaobj` versions, perhaps even your own custom UPS product.

## Inputs
The package is controlled by `fhicl` config files, found in the `cfg` directory. The `cfg/ndcafmakerjob.fcl` shows the basic setup.

For the minimal test setup:
* Provide an `InputGHEPFile`, which contains the GENIE truth information in `GHEP` format
* Provide an `OutputFile`, where your file will be saved

Extending upon the minimal test setup you can:
* Provide a `NDLArRecoFile`, which contains the output of the ND LAr reconstruction
* Provide a `TMSRecoFile`, which contains the output of the TMS reconstruction
* Provide a `SANDRecoFile`, which contains the output of the SAND reconstruction

To add variables and inspect what is set and how, check `src/Params.h`.


## Building
Once you've set up your environment, it's just a matter of typing 
```
make
```
in the `ND_CAFMaker` folder, which goes through and builds the objects, library and single `makeCAF` executable.

# Running
There is currently only one main executable, `makeCAF`, which is controlled entirely by the input `fhicl` file. To run with the `fhicl` file that was specified under `Inputs`, do
```
./makeCAF --fcl=path_to_your_fhicl_file.fcl
```

You can also override some of the `fhicl` inputs, which are specified by typing 
```
./makeCAF --help
```
and should output
```
Usage: ./makeCAF [options] <driver.fcl> [options]

General options:
  -h [ --help ]          print this help message

FCL overrides (for quick tests; edit your .fcl for regular usage):
  -g [ --ghep ] arg      input GENIE .ghep file
  -o [ --out ] arg       output CAF file
  --startevt arg         event number to start at
  -n [ --numevts ] arg   total number of events to process (-1 means 'all')
```

# Output tree and event format (this section needs expanding)

The output contains a number of different `TTree` `ROOT` objects. `cafTree` contains the information from the reconstruction and some truth information from the `edep-sim` detector simulation, and `genieEvt` contains the true `GENIE` information from the neutrino interaction simulation.

# Contact
* Jeremy Wolcott (jwolcott@fnal.gov), most certainly the lead author
* Chris Marshall (chris.marshall@rochester.edu), originally wrote much of the package before Jeremy's significant update
* Cris Vilela (c.vilela@cern.ch), also originally wrote much of the package before Jeremy's significant update
* Clarence Wret (c.wret@rochester.edu), tagged on TMS reconstruction and general updates

# Old instructions (semi-valid still)
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
