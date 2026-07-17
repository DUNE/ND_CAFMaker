# Spack packaging (draft)

This directory holds a draft Spack `package.py` for `ND_CAFMaker`, intended
to be contributed to the [`DUNE/dune_spack`](https://github.com/DUNE/dune_spack)
repository at `packages/nd_cafmaker/package.py`. It is **not** consumed by
this repository's own build (see the top-level `README.md` for how to build
`ND_CAFMaker` itself, either via UPS/CVMFS or a Spack-provided environment).

It is kept here only because packaging work usually happens in `dune_spack`
directly, but this file is provided so it can be reviewed and copied over by
someone with write access to that repository.

## Open items before submitting to `dune_spack`

- Confirm the exact `genie` package spec (name/version/variants) against the
  `dunesw-10_20_03d01-justin-01_06_01-prototype` Spack environment
  (`/cvmfs/dune.opensciencegrid.org/spack/setup-env.sh`).
- Confirm `fhiclcpp`, `duneanaobj`, `edep-sim`, `py-srproxy`, `sandreco`, and
  `log4cpp` package names/specs match what's actually available in
  `dune_spack` (some already exist there; verify others need to be added).
- Validate `spack install nd_cafmaker` end-to-end on AlmaLinux 9, including
  the `+tms`/`~tms` and `+sand`/`~sand` variants.
