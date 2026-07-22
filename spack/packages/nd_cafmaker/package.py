# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
# NOTE: This file is a proposed contribution to the DUNE/dune_spack repository
# (https://github.com/DUNE/dune_spack), and is expected to eventually live at
# packages/nd_cafmaker/package.py over there. It is kept here, in the
# ND_CAFMaker repository itself, only because this sandbox has no write
# access to DUNE/dune_spack. Copy it verbatim to
# <dune_spack>/packages/nd_cafmaker/package.py to add the package.

from spack.package import *


class NdCafmaker(CMakePackage):
    """ND_CAFMaker -- combines edep-sim, GENIE, and reconstructed objects
    from the DUNE Near Detector into the Common Analysis Format (CAF)."""

    homepage = "https://github.com/DUNE/ND_CAFMaker"
    url = "https://github.com/DUNE/ND_CAFMaker/archive/refs/tags/v5.1.1.tar.gz"
    git = "https://github.com/DUNE/ND_CAFMaker.git"

    maintainers("alexbooth92")

    license("MIT")

    version("main", branch="main")

    # Only v5.1.0 and later build with CMake (earlier releases used a
    # hand-rolled GNU Makefile and are not supported by this package).
    version("5.1.1", sha256="7ced72723ff2059d10e6f985ca705e958b1a9a95bc76f1a345600c39edf1f0f5")
    version("5.1.0", sha256="fde758ba763e9b1da0e6f70260f626c0b87a8a699beffe1e1cbde01ecaf9cb8d")

    variant(
        "cxxstd",
        default="17",
        values=("17",),
        multi=False,
        description="Use the specified C++ standard when building.",
    )
    variant("tms", default=True, description="Enable TMS reconstruction branch filler")
    variant("testexe", default=False, description="Build the testHDF test executable")
    variant("sand", default=False, description="Enable SANDReco branch filler support")

    depends_on("cxx", type="build")
    depends_on("cmake@3.20:", type="build")

    depends_on("root +mathmore +geom +pythia6")
    depends_on("edep-sim")
    depends_on("duneanaobj")
    depends_on("fhiclcpp")
    depends_on("boost +program_options")
    depends_on("hdf5 +cxx")
    depends_on("log4cpp")
    depends_on("gsl")
    depends_on("lhapdf")
    depends_on("nlohmann-json")
    depends_on("py-srproxy")
    depends_on("libxml2")
    depends_on("pythia6")
    depends_on("curl")

    # GENIE and its data tables (genie_xsec, genie_phyopt in the UPS build).
    # This is the trickiest dependency to translate: ND_CAFMaker's CMake
    # integration shells out to `genie-config --libs` rather than consuming
    # a CMake config package, so no extra Spack-side wiring is required
    # beyond having `genie-config` on PATH. The exact `genie` spec/version
    # this should pin to needs to be confirmed against the
    # dunesw-10_20_03d01-justin-01_06_01-prototype Spack environment
    # (available via /cvmfs/dune.opensciencegrid.org/spack/setup-env.sh),
    # which was not reachable from the sandbox that authored this file.
    depends_on("genie")

    # SANDReco is optional and auto-detected upstream via SANDRECO_INC/LIB
    # env vars; gate it behind a variant instead so Spack builds are
    # reproducible regardless of ambient environment state.
    depends_on("sandreco", when="+sand")

    def cmake_args(self):
        args = [
            self.define_from_variant("CMAKE_CXX_STANDARD", "cxxstd"),
            self.define_from_variant("ENABLE_TMS", "tms"),
            self.define_from_variant("ENABLE_TESTEXE", "testexe"),
        ]
        if self.spec.satisfies("+sand"):
            spec = self.spec["sandreco"]
            args.append(self.define("SANDRECO_INC", spec.prefix.include))
            args.append(self.define("SANDRECO_LIB", spec.prefix.lib))
        return args

    def setup_run_environment(self, env):
        env.append_path("FHICL_FILE_PATH", "{0}/cfg".format(self.prefix))

    def setup_dependent_run_environment(self, env, dependent_spec):
        env.append_path("FHICL_FILE_PATH", "{0}/cfg".format(self.prefix))
