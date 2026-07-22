# cmake-format: off
# ── find_ups_package ───────────────────────────────────────────────────────────
# Creates an INTERFACE IMPORTED target for a dependency that ships neither a
# CMake config package nor a discoverable pkg-config file.
#
# Two discovery modes are supported:
#   - UPS mode (SL7/CVMFS): the package is located via an environment
#     variable set by ndcaf_setup.sh (e.g. HDF5_LIB), searched with
#     NO_DEFAULT_PATH so nothing from the system is picked up by accident.
#   - Standard/Spack mode: when the relevant env var is *not* set, library
#     lookups fall back to a normal find_library() search, which honors
#     CMAKE_PREFIX_PATH (as populated by Spack's dependency environment)
#     and standard system paths. Include-dir lookups (INC_VAR/INC_SUFFIX),
#     however, are simply skipped in this mode: INC_SUFFIX is a directory
#     suffix rather than a header filename, so it can't be passed to
#     find_path(). Callers relying on INC_VAR/INC_SUFFIX in Spack mode must
#     ensure the include dir is already visible to the compiler (e.g. via
#     CMAKE_PREFIX_PATH/CPATH, or because the dependency ships a proper
#     CMake config target instead of going through find_ups_package()).
#
# Usage:
#   find_ups_package(
#     TARGET_NAME <target>      # e.g. deps::hdf5
#     [INC_VAR    <env-var>]    # include dir env var, e.g. HDF5_INC
#     [INC_SUFFIX <subdir>]     # appended to inc dir, e.g. EDepSim
#     [LIB_VAR    <env-var>]    # lib dir env var, e.g. HDF5_LIB
#     [LIBS       <lib1> ...]   # library names without lib prefix / .so
#     [REQUIRED]
#   )
# cmake-format: on

function(find_ups_package)
  cmake_parse_arguments(UPS "REQUIRED" "TARGET_NAME;INC_VAR;INC_SUFFIX;LIB_VAR"
                        "LIBS" ${ARGN})

  if(NOT UPS_TARGET_NAME)
    message(FATAL_ERROR "find_ups_package: TARGET_NAME is required")
  endif()

  if(UPS_INC_VAR)
    if(DEFINED ENV{${UPS_INC_VAR}})
      set(_inc_dir "$ENV{${UPS_INC_VAR}}")
      if(UPS_INC_SUFFIX)
        set(_inc_dir "${_inc_dir}/${UPS_INC_SUFFIX}")
      endif()
    else()
      # Not running under UPS (e.g. Spack build): UPS_INC_SUFFIX is a
      # directory suffix, not a header filename, so it can't be used with
      # find_path(). Leave _inc_dir unset here; Spack-provided packages are
      # expected to expose their include dir via CMAKE_PREFIX_PATH/CPATH
      # (or a proper CMake config), so the compiler picks it up without
      # needing an explicit INTERFACE_INCLUDE_DIRECTORIES entry.
    endif()
  endif()

  set(_libs "")
  if(UPS_LIB_VAR AND UPS_LIBS)
    foreach(_lib IN LISTS UPS_LIBS)
      if(DEFINED ENV{${UPS_LIB_VAR}})
        find_library(
          _found_${_lib}
          NAMES ${_lib}
          PATHS "$ENV{${UPS_LIB_VAR}}"
          NO_DEFAULT_PATH)
      else()
        # Not running under UPS (e.g. Spack build): search
        # CMAKE_PREFIX_PATH and system locations instead.
        find_library(_found_${_lib} NAMES ${_lib})
      endif()
      if(_found_${_lib})
        list(APPEND _libs "${_found_${_lib}}")
      elseif(UPS_REQUIRED)
        if(DEFINED ENV{${UPS_LIB_VAR}})
          message(
            FATAL_ERROR
              "find_ups_package: '${_lib}' not found under $ENV{${UPS_LIB_VAR}}"
          )
        else()
          message(
            FATAL_ERROR
              "find_ups_package: '${_lib}' not found (set the ${UPS_LIB_VAR} environment variable, or populate CMAKE_PREFIX_PATH)"
          )
        endif()
      endif()
    endforeach()
  endif()

  add_library(${UPS_TARGET_NAME} INTERFACE IMPORTED)
  if(_inc_dir)
    set_target_properties(
      ${UPS_TARGET_NAME} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_inc_dir}")
  endif()
  if(_libs)
    set_target_properties(${UPS_TARGET_NAME} PROPERTIES INTERFACE_LINK_LIBRARIES
                                                        "${_libs}")
  endif()
endfunction()
