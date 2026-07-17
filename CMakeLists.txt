cmake_minimum_required(VERSION 3.20)
project(ND_CAFMaker LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 5.1.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_TMS "Enable TMS reconstruction branch filler" ON)
option(ENABLE_TESTEXE "Build the testHDF test executable" OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address,undefined")

string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address,undefined")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -D_GLIBCXX_SANITIZE_STD_ALLOCATOR")
endif()

add_compile_definitions(NO_ART)

# Make the installed makeCAF find libND_CAFMaker.so without LD_LIBRARY_PATH:
# embed $ORIGIN/../lib as an RPATH. CMAKE_INSTALL_RPATH_USE_LINK_PATH also
# bakes in the link-time dependency paths (CVMFS UPS libs) so the executable
# runs even from a shell that hasn't sourced ndcaf_setup.sh.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

include(cmake/FindUPSPackage.cmake)

# Versions are pinned by ndcaf_setup_deps.sh; env vars (and CMAKE_PREFIX_PATH for
# the cetmodules-built packages) are set by UPS at runtime.

# ROOT — always via ROOTConfig.cmake (ships with ROOT >= 6)
find_package(ROOT 6 REQUIRED COMPONENTS Core MathMore Geom EGPythia6 GenVector)

# UPS ships sqlite as libsqlite3_ups.so (not libsqlite3.so) and the SL7 base
# image only has libsqlite3.so.0. Pre-populate the FindSQLite3 cache so
# fhiclcpp -> cetlib's find_dependency(SQLite3) finds the UPS install.
# Under Spack (or any other environment providing a normal libsqlite3.so via
# CMAKE_PREFIX_PATH), SQLITE_FQ_DIR is unset and the standard find_dependency(SQLite3)
# lookup triggered by fhiclcpp/cetlib is left to run unmodified.
if(DEFINED ENV{SQLITE_FQ_DIR})
  find_library(SQLite3_LIBRARY
    NAMES sqlite3_ups
    PATHS "$ENV{SQLITE_FQ_DIR}/lib"
    NO_DEFAULT_PATH)
  find_path(SQLite3_INCLUDE_DIR
    NAMES sqlite3.h
    PATHS "$ENV{SQLITE_FQ_DIR}/include"
    NO_DEFAULT_PATH)
endif()

# Force TBB to be picked up via TBB's own CMake config (it ships TBB::tbb).
# Without this, fhiclcpp -> find_dependency(TBB) loads Geant4's FindTBB.cmake
# module instead, which only sets variables and never creates the target.
find_package(TBB CONFIG REQUIRED)

# Modern UPS products built with cetmodules ship CMake configs under
# ${PRODUCT}_FQ_DIR/lib/<product>/cmake/. ndcaf_setup_deps.sh prepends those
# directories to CMAKE_PREFIX_PATH, so find_package(... CONFIG) works.
find_package(duneanaobj CONFIG REQUIRED)
find_package(EDepSim CONFIG REQUIRED)
find_package(fhiclcpp CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS program_options)

find_program(GENIE_CONFIG genie-config HINTS "$ENV{GENIE_FQ_DIR}/bin"
                                             REQUIRED)

# GENIE ships no CMake config package and no usable UPS *_LIB variable —
# integration goes through the legacy `genie-config --libs` script, which
# returns raw linker flags ("-L<path> -lGFwReg -lGFwAlg ..."). We wrap that
# into an INTERFACE IMPORTED target so consumers can just link deps::genie
# uniformly alongside the other deps::* targets.
execute_process(
  COMMAND ${GENIE_CONFIG} --libs
  OUTPUT_VARIABLE _genie_libs
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(DEFINED ENV{GENIE_FQ_DIR})
  set(_genie_inc_dir "$ENV{GENIE_FQ_DIR}/include/GENIE")
else()
  # Not running under UPS (e.g. Spack build): derive the include dir from
  # the location of genie-config itself (<prefix>/bin/genie-config).
  get_filename_component(_genie_bin_dir "${GENIE_CONFIG}" DIRECTORY)
  get_filename_component(_genie_prefix "${_genie_bin_dir}" DIRECTORY)
  set(_genie_inc_dir "${_genie_prefix}/include/GENIE")
endif()
add_library(deps::genie INTERFACE IMPORTED)
set_target_properties(
  deps::genie PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_genie_inc_dir}"
                         INTERFACE_LINK_LIBRARIES "${_genie_libs}")

# Side effect of the raw-flag integration above: CMake never sees an absolute
# path to any GENIE library, only "-l<name>" tokens. CMAKE_INSTALL_RPATH_USE_LINK_PATH
# only walks paths it knows about, so the GENIE lib dir would be missing from
# the installed binary's RPATH and ld.so would fail to resolve libGFwMsg & co.
# at runtime (unless LD_LIBRARY_PATH is set). Extract the "-L<path>" token from
# genie-config's output and append it to the install RPATH by hand.
string(REGEX MATCH "-L([^ ]+)" _genie_ldir_match "${_genie_libs}")
if(CMAKE_MATCH_1)
  list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_MATCH_1}")
endif()

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

if(DEFINED ENV{LIBXML2_FQ_DIR})
  find_library(
    LIBXML2_LIB_PATH
    NAMES xml2
    PATHS "$ENV{LIBXML2_FQ_DIR}/lib"
    NO_DEFAULT_PATH)
else()
  find_library(LIBXML2_LIB_PATH NAMES xml2)
endif()
add_library(deps::libxml2 INTERFACE IMPORTED)
set_target_properties(deps::libxml2 PROPERTIES INTERFACE_LINK_LIBRARIES
                                               "${LIBXML2_LIB_PATH}")

if(DEFINED ENV{PYTHIA6_LIBRARY})
  find_library(
    PYTHIA6_LIB_PATH
    NAMES Pythia6
    PATHS "$ENV{PYTHIA6_LIBRARY}"
    NO_DEFAULT_PATH)
else()
  find_library(PYTHIA6_LIB_PATH NAMES Pythia6)
endif()
add_library(deps::pythia6 INTERFACE IMPORTED)
set_target_properties(deps::pythia6 PROPERTIES INTERFACE_LINK_LIBRARIES
                                               "${PYTHIA6_LIB_PATH}")

# curl — no dedicated UPS product/env var; find_library()/find_path() search
# default system locations (which include /usr/lib64 on RHEL-family systems)
# as well as CMAKE_PREFIX_PATH, so this works both on the SL7/UPS base image
# and in a Spack environment.
find_library(CURL_LIB_PATH NAMES curl)
find_path(CURL_INC_PATH NAMES curl/curl.h)
if(NOT CURL_LIB_PATH)
  message(FATAL_ERROR "curl library not found. Install curl (e.g. the 'libcurl-devel'/'libcurl-dev' package) or set CMAKE_PREFIX_PATH to its location.")
endif()
if(NOT CURL_INC_PATH)
  message(FATAL_ERROR "curl/curl.h not found. Install curl development headers (e.g. the 'libcurl-devel'/'libcurl-dev' package) or set CMAKE_PREFIX_PATH to their location.")
endif()
add_library(deps::curl INTERFACE IMPORTED)
set_target_properties(
  deps::curl PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CURL_INC_PATH}"
                        INTERFACE_LINK_LIBRARIES "${CURL_LIB_PATH}")

# SANDReco — optional, auto-detected from environment (both modes)
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
  message(STATUS "SANDReco: disabled (SANDRECO_INC/LIB not set)")
endif()

add_subdirectory(src)

# Install the setup scripts next to the makeCAF executable so they can
# locate the install tree via $BASH_SOURCE. Install fcl/cfg files alongside.
install(PROGRAMS ndcaf_setup.sh ndcaf_setup_deps.sh DESTINATION bin)
install(DIRECTORY cfg/ DESTINATION cfg)
