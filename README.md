# ND CAFMaker
`ND_CAFMaker` takes input `edep-sim`, `GENIE`, and reconstructed objects from the DUNE ND and combines them into the Common Analysis Format ("CAF").

# Setup
The package has a number of dependencies, all accessible through the `ups` framework at Fermilab and github. 
The first time you use it'll you'll need to build the dependencies.  Simply run
```
./build_deps.sh
```

For subsequent use, two scripts are provided:

- **`ndcaf_setup_deps.sh`** — sets up the UPS dependencies (ROOT, GENIE, edep-sim,
  duneanaobj, fhiclcpp, HDF5, ...). Source this **before** running `cmake`.
- **`ndcaf_setup.sh`** — sets up the runtime environment from an installed tree
  (puts `makeCAF` on `PATH`, `libND_CAFMaker.so` on `LD_LIBRARY_PATH`, and
  `<prefix>/cfg` on `FHICL_FILE_PATH`). It internally sources `ndcaf_setup_deps.sh`,
  so users running an installed `makeCAF` only need to source `ndcaf_setup.sh`.

Both scripts take an optional build qualifier (`prof` or `debug`).
(If you're hacking on the CAFMaker, you probably want `debug`.  If not, `prof` will run faster.)

**If, however, you need to set up a local version of `duneanaobj` for testing purposes, you must edit `ndcaf_setup_deps.sh` before sourcing it.**

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

### Build environment

The build requires the UPS dependencies provided by the Fermilab software stack.
The recommended way to get a compatible environment is via the Singularity image:

```
singularity shell --bind /cvmfs /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
```

Once inside the container, source the dependency setup script before configuring:

```
source ndcaf_setup_deps.sh {prof|debug}
```

#### Building with Spack instead of UPS

`CMakeLists.txt` and `cmake/FindUPSPackage.cmake` also support locating
dependencies without UPS: any dependency whose UPS `*_FQ_DIR`/`*_LIBRARY`
environment variable is unset falls back to a normal CMake search that
honors `CMAKE_PREFIX_PATH`. This lets the project be built inside a Spack
environment (e.g. on AlmaLinux 9) that provides `root`, `edep-sim`,
`duneanaobj`, `fhiclcpp`, `hdf5`, `log4cpp`, `gsl`, `lhapdf`, `genie`,
`libxml2`, `pythia6`, and `curl` on `CMAKE_PREFIX_PATH`/`PATH`, instead of
sourcing `ndcaf_setup_deps.sh`. A Spack package definition for `ND_CAFMaker`
itself lives in the [`DUNE/dune_spack`](https://github.com/DUNE/dune_spack)
repository, not here.

### CMake configuration

Once you've set up your environment, configure and build with CMake.
From the project root folder:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug # or Release
cmake --build build 
```

The build products are installed under `build/` by default. The `makeCAF` executable will be at `build/bin/makeCAF`.

### CMake build options

| Option | Default | Description |
|--------|---------|-------------|
| `ENABLE_TMS` | `ON` | Enable TMS reconstruction branch filler |
| `CMAKE_BUILD_TYPE` | — | Set to `Debug` or `Release` (overrides the `-g -O2` defaults) |

Example: disable TMS and build the test executable:

```
cmake -S . -B build -DENABLE_TMS=OFF
cmake --build build
```

### Installing

To install the library and executables to a custom prefix:

```
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build --target install
```

This will place:
- `libND_CAFMaker.so` → `<prefix>/lib/`
- `makeCAF` (and optionally `testHDF`) → `<prefix>/bin/`
- `ndcaf_setup.sh`, `ndcaf_setup_deps.sh` → `<prefix>/bin/`
- the `cfg/` fcl files → `<prefix>/cfg/`

If no `CMAKE_INSTALL_PREFIX` is given, the default system prefix (`/usr/local`) is used.

The installed `makeCAF` is built with `$ORIGIN/../lib` (plus the CVMFS paths
of all UPS dependencies) in its `RUNPATH`, so it loads `libND_CAFMaker.so`
without any `LD_LIBRARY_PATH` manipulation. Sourcing `ndcaf_setup.sh` is still
needed to set up the GCC/UPS runtime environment and to point
`FHICL_FILE_PATH` at the installed `cfg/` directory.

# Running
After install, source the runtime setup once per shell to get `makeCAF`
on your `PATH`:
```
source <prefix>/bin/ndcaf_setup.sh {prof|debug}
makeCAF --fcl=path_to_your_fhicl_file.fcl
```

You can also override some of the `fhicl` inputs, which are specified by typing 
```
makeCAF --help
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

## Known Limitations

ROOT's `TTreeFormula` class has an intrinsic limitation: it supports only 2 levels of variable-sized collections. This prevents correct scanning of complex `duneanaobj` fields using `TTree::Scan` in 'structured' CAFs (those that directly embed `StandardRecord` objects).  For example, scanning fields like `rec.nd.lar.tracks.truthOverlap` may produce incomplete or incorrect output, with some entries not properly displayed.

Note: the 'flat' CAFs also emitted by the CAFMaker do not suffer from this limitation.

For more details and examples, see [GitHub Issue #142](https://github.com/DUNE/ND_CAFMaker/issues/142).

## Using a locally installed `duneanaobj` for testing

`duneanaobj`'s [README](https://github.com/DUNE/duneanaobj/blob/master/README.md) explains how to construct a custom UPS product that can be used in place of the version used by default by `ndcaf_setup_deps.sh`.
If you have such a custom version that you wish to use, edit `ndcaf_setup_deps.sh` to add your custom UPS directory to `$PRODUCTS` and to set up the version name for your working copy.
Consult the README for more details.

# Contact
* Jeremy Wolcott (jwolcott@fnal.gov) current lead author / maintainer

### Other major contributors
* Chris Marshall (chris.marshall@rochester.edu)
* Clarence Wret (clarence.wret@physics.ox.ac.uk)
* Cris Vilela (c.vilela@cern.ch)
* Noë Roy (noe.roylamoureux@gmail.com)
* Sindhujha Kumaran (s.kumaran@uci.edu)

---------------------------------

# Licensing
Copyright © 2023 FERMI NATIONAL ACCELERATOR LABORATORY for the benefit
of the DUNE Collaboration.

This repository, and all software contained within, is licensed under
the Apache License, Version 2.0 (the "License"); you may not use its code
except in compliance with the License. You may obtain a copy of
the License at
[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0);
the license is also available in the LICENSE.md file in this repository.

Copyright is granted to FERMI NATIONAL ACCELERATOR LABORATORY on behalf
of the Deep Underground Neutrino Experiment (DUNE). Unless required by
applicable law or agreed to in writing, software distributed under the
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for
the specific language governing permissions and limitations under the
License.
