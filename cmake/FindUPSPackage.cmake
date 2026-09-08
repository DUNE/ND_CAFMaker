# cmake-format: off
# ── find_ups_package ───────────────────────────────────────────────────────────
# Creates an INTERFACE IMPORTED target from a UPS package located via
# environment variables set by ndcaf_setup.sh (NOT from system paths).
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
    set(_inc_dir "$ENV{${UPS_INC_VAR}}")
    if(UPS_INC_SUFFIX)
      set(_inc_dir "${_inc_dir}/${UPS_INC_SUFFIX}")
    endif()
  endif()

  set(_libs "")
  if(UPS_LIB_VAR AND UPS_LIBS)
    foreach(_lib IN LISTS UPS_LIBS)
      find_library(
        _found_${_lib}
        NAMES ${_lib}
        PATHS "$ENV{${UPS_LIB_VAR}}"
        NO_DEFAULT_PATH)
      if(_found_${_lib})
        list(APPEND _libs "${_found_${_lib}}")
      elseif(UPS_REQUIRED)
        message(
          FATAL_ERROR
            "find_ups_package: '${_lib}' not found in $ENV{${UPS_LIB_VAR}}")
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
