# ND CAFMaker
`ND_CAFMaker` takes input `edep-sim`, `GENIE`, and reconstructed objects from the DUNE ND and combines them into the Common Analysis Format ("CAF").

# Setup
The package has a number of dependencies, all accessible through the `ups` framework at Fermilab and github. 
The first time you use it'll you'll need to build the dependencies.  Simply run
```
./build.sh
```

For subsequent use, 
```
source ndcaf_setup.sh
```
is sufficient to set up your environment for the common use cases.
**If, however, you need to set up a local version of `duneanaobj` for testing purposes, you must edit `ndcaf_setup.sh` before sourcing it.**

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

# Output tree and event format

The output contains a number of different `TTree` `ROOT` objects. `cafTree` contains the information from the reconstruction and some truth information from the `edep-sim` detector simulation, and `genieEvt` contains the true `GENIE` information from the neutrino interaction simulation.

The object format for the `StandardRecord` objects is defined in [`duneanaobj`](https://github.com/DUNE/duneanaobj).
See the README there for instructions on how to set it up and modify it.

## Using a locally installed `duneanaobj` for testing

`duneanaobj`'s [README](https://github.com/DUNE/duneanaobj/blob/master/README.md) explains how to construct a custom UPS product that can be used in place of the version used by default by `ndcaf_setup.sh`.
If you have such a custom version that you wish to use, edit `ndcaf_setup.sh` to add your custom UPS directory to `$PRODUCTS` and to set up the version name for your working copy.
Consult the README for more details.

# Contact
* Jeremy Wolcott (jwolcott@fnal.gov) current lead author / maintainer

### Other major contributors
* Chris Marshall (chris.marshall@rochester.edu)
* Clarence Wret (clarence.wret@physics.ox.ac.uk)
* Cris Vilela (c.vilela@cern.ch)
* Noë Roy (noe.roylamoureux@gmail.com)
* Sindhujha Kumaran (s.kumaran@uci.edu)

# Licensing
Copyright © 2023 FERMI NATIONAL ACCELERATOR LABORATORY for the benefit
of the DUNE Collaboration.

This repository, and all software contained within, is licensed under
the Apache License, Version 2.0 (the "License"); you may not use its code
except in compliance with the License. You may obtain a copy of
the License at
     [http://www.apache.org/licenses/LICENSE-2.0];
The license is also available in the LICENSE.md file in this repository.

Copyright is granted to FERMI NATIONAL ACCELERATOR LABORATORY on behalf
of the Deep Underground Neutrino Experiment (DUNE). Unless required by
applicable law or agreed to in writing, software distributed under the
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for
the specific language governing permissions and limitations under the
License.
