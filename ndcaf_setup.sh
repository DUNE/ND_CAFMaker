#!/bin/bash

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
	echo "ERROR: setup.sh must be sourced, not executed." >&2
	exit 1
fi

_setup_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Pull in UPS dependencies first (forwards the build qualifier).
source "${_setup_dir}/ndcaf_setup_deps.sh" "$@"

# Put the installed ND_CAFMaker on the user's environment. Paths are derived
# from this script's location, so any install prefix works.
export PATH="${_setup_dir}:${PATH}"
export LD_LIBRARY_PATH="${_setup_dir}/install/lib:${LD_LIBRARY_PATH}"
export FHICL_FILE_PATH="${_setup_dir}/install/cfg:${FHICL_FILE_PATH}"

unset _setup_dir
