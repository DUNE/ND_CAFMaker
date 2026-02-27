cmake_minimum_required(VERSION 3.20)

# Derive version from the latest git tag (e.g. v4.9.0 → 4.9.0). Falls back to
# 0.0.0 when building outside a git repository.
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --match "v*" --abbrev=0
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE _git_tag
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET
    RESULT_VARIABLE _git_rc)
  if(_git_rc EQUAL 0)
    string(REGEX REPLACE "^v" "" _version "${_git_tag}")
  else()
    set(_version "0.0.0")
  endif()
else()
  set(_version "0.0.0")
endif()

project(
  ND_CAFMaker
  VERSION ${_version}
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_TMS "Enable TMS reconstruction branch filler" ON)
option(ENABLE_TESTEXE "Build the testHDF test executable" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_options(-Wall)
add_compile_definitions(NO_ART)

# cmake-format: off
# ── UPS package helper ─────────────────────────────────────────────────────────
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

# ROOT — imported targets from ROOTConfig.cmake (ships with ROOT >= 6)
find_package(ROOT REQUIRED COMPONENTS Core MathMore Geom EGPythia6 GenVector)

add_library(deps::root_glibs INTERFACE IMPORTED)
set_target_properties(
  deps::root_glibs
  PROPERTIES
    INTERFACE_LINK_LIBRARIES
    "ROOT::Core;ROOT::MathMore;ROOT::Geom;ROOT::EGPythia6;ROOT::GenVector")

# GENIE — uses genie-config --libs
find_program(GENIE_CONFIG genie-config HINTS "$ENV{GENIE_INC}/../bin" REQUIRED)

execute_process(
  COMMAND ${GENIE_CONFIG} --libs
  OUTPUT_VARIABLE _genie_libs
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(deps::genie INTERFACE IMPORTED)

set_target_properties(
  deps::genie PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "$ENV{GENIE_INC}/GENIE"
                         INTERFACE_LINK_LIBRARIES "${_genie_libs}")

# HDF5 — library name differs between debug and prof builds
if("$ENV{HDF5_FQ_DIR}" MATCHES "debug")
  set(_hdf5_lib hdf5_cpp_debug)
else()
  set(_hdf5_lib hdf5_cpp)
endif()

find_ups_package(
  TARGET_NAME deps::hdf5
  INC_VAR HDF5_INC
  LIB_VAR HDF5_LIB
  LIBS ${_hdf5_lib}
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::log4cpp
  INC_VAR LOG4CPP_INC
  LIB_VAR LOG4CPP_LIB
  LIBS log4cpp
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::edepsim
  INC_VAR EDEPSIM_INC
  INC_SUFFIX EDepSim
  LIB_VAR EDEPSIM_LIB
  LIBS edepsim
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::boost_po
  INC_VAR BOOST_INC
  LIB_VAR BOOST_LIB
  LIBS boost_program_options
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::fhiclcpp
  INC_VAR FHICLCPP_INC
  LIB_VAR FHICLCPP_LIB
  LIBS fhiclcpp fhiclcpp_types
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::duneanaobj
  INC_VAR DUNEANAOBJ_INC
  LIB_VAR DUNEANAOBJ_LIB
  LIBS duneanaobj_StandardRecord duneanaobj_StandardRecordFlat
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::tbb
  LIB_VAR TBB_LIB
  LIBS tbb
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::gsl
  LIB_VAR GSL_LIB
  LIBS gsl gslcblas
  REQUIRED)

find_ups_package(
  TARGET_NAME deps::lhapdf
  LIB_VAR LHAPDF_LIB
  LIBS LHAPDF
  REQUIRED)

find_ups_package(TARGET_NAME deps::nlohmann INC_VAR NLOHMANN_JSON_INC)

find_ups_package(TARGET_NAME deps::srproxy INC_VAR SRPROXY_INC)

# cetlib + cetlib_except — two env-var pairs merged into one target
find_library(
  CETLIB_LIB_PATH
  NAMES cetlib
  PATHS "$ENV{CETLIB_LIB}"
  NO_DEFAULT_PATH)

find_library(
  CETLIB_EXCEPT_LIB_PATH
  NAMES cetlib_except
  PATHS "$ENV{CETLIB_EXCEPT_LIB}"
  NO_DEFAULT_PATH)

add_library(deps::cetlib INTERFACE IMPORTED)

set_target_properties(
  deps::cetlib
  PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
             "$ENV{CETLIB_INC};$ENV{CETLIB_EXCEPT_INC}"
             INTERFACE_LINK_LIBRARIES
             "${CETLIB_LIB_PATH};${CETLIB_EXCEPT_LIB_PATH}")

# libxml2 — LIBXML2_FQ_DIR points to product root, lib is in its lib/ subdir
find_library(
  LIBXML2_LIB_PATH
  NAMES xml2
  PATHS "$ENV{LIBXML2_FQ_DIR}/lib"
  NO_DEFAULT_PATH)

add_library(deps::libxml2 INTERFACE IMPORTED)

set_target_properties(deps::libxml2 PROPERTIES INTERFACE_LINK_LIBRARIES
                                               "${LIBXML2_LIB_PATH}")

# Pythia6 — PYTHIA6_LIBRARY points directly to the lib dir
find_library(
  PYTHIA6_LIB_PATH
  NAMES Pythia6
  PATHS "$ENV{PYTHIA6_LIBRARY}"
  NO_DEFAULT_PATH)

add_library(deps::pythia6 INTERFACE IMPORTED)

set_target_properties(deps::pythia6 PROPERTIES INTERFACE_LINK_LIBRARIES
                                               "${PYTHIA6_LIB_PATH}")

# curl — system library
find_library(
  CURL_LIB_PATH
  NAMES curl
  PATHS /usr/lib64
  NO_DEFAULT_PATH)

add_library(deps::curl INTERFACE IMPORTED)

set_target_properties(
  deps::curl PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "/usr/include/curl"
                        INTERFACE_LINK_LIBRARIES "${CURL_LIB_PATH}")

# SANDReco — optional, auto-detected from environment
if(DEFINED ENV{SANDRECO_INC} OR DEFINED ENV{SANDRECO_LIB})
  set(SAND_ENABLED TRUE)

  find_library(
    SANDRECO_LIB_PATH
    NAMES Struct
    PATHS "$ENV{SANDRECO_LIB}"
    NO_DEFAULT_PATH)

  add_library(deps::sand INTERFACE IMPORTED)

  set_target_properties(
    deps::sand PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "$ENV{SANDRECO_INC}"
                          INTERFACE_LINK_LIBRARIES "${SANDRECO_LIB_PATH}")

  message(STATUS "SANDReco: enabled (SANDRECO_INC=$ENV{SANDRECO_INC})")
else()
  set(SAND_ENABLED FALSE)

  add_library(deps::sand INTERFACE IMPORTED) # no-op target for unconditional
                                             # linking
  message(STATUS "SANDReco: disabled (SANDRECO_INC/LIB not set)")
endif()

add_subdirectory(src)
