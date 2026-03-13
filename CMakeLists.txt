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

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

string(
  APPEND
  CMAKE_CXX_FLAGS
  " -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion"
  " -Wshadow -Wimplicit-fallthrough -Wextra-semi -Wold-style-cast"
  " -fno-omit-frame-pointer")

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address,undefined")

string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address,undefined")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -D_GLIBCXX_SANITIZE_STD_ALLOCATOR")
endif()

add_compile_definitions(NO_ART)

include(cmake/FindUPSPackage.cmake)

# ROOT — always via ROOTConfig.cmake (ships with ROOT >= 6)
find_package(ROOT 6 REQUIRED COMPONENTS Core MathMore Geom EGPythia6 GenVector)
add_library(deps::root_glibs INTERFACE IMPORTED)
set_target_properties(
  deps::root_glibs
  PROPERTIES
    INTERFACE_LINK_LIBRARIES
    "ROOT::Core;ROOT::MathMore;ROOT::Geom;ROOT::EGPythia6;ROOT::GenVector")

# Versions are pinned by ndcaf_setup.sh; env vars are set by UPS at runtime.

find_program(GENIE_CONFIG genie-config HINTS "$ENV{GENIE_INC}/../bin"
                                             REQUIRED)

execute_process(
  COMMAND ${GENIE_CONFIG} --libs
  OUTPUT_VARIABLE _genie_libs
  OUTPUT_STRIP_TRAILING_WHITESPACE)
add_library(deps::genie INTERFACE IMPORTED)
set_target_properties(
  deps::genie PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "$ENV{GENIE_INC}/GENIE"
                         INTERFACE_LINK_LIBRARIES "${_genie_libs}")

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
  LIB_VAR EDEPSIM_LIB
  LIBS edepsim
  INC_SUFFIX EDepSim
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

find_library(
  LIBXML2_LIB_PATH
  NAMES xml2
  PATHS "$ENV{LIBXML2_FQ_DIR}/lib"
  NO_DEFAULT_PATH)
add_library(deps::libxml2 INTERFACE IMPORTED)
set_target_properties(deps::libxml2 PROPERTIES INTERFACE_LINK_LIBRARIES
                                               "${LIBXML2_LIB_PATH}")

find_library(
  PYTHIA6_LIB_PATH
  NAMES Pythia6
  PATHS "$ENV{PYTHIA6_LIBRARY}"
  NO_DEFAULT_PATH)
add_library(deps::pythia6 INTERFACE IMPORTED)
set_target_properties(deps::pythia6 PROPERTIES INTERFACE_LINK_LIBRARIES
                                               "${PYTHIA6_LIB_PATH}")

find_library(
  CURL_LIB_PATH
  NAMES curl
  PATHS /usr/lib64
  NO_DEFAULT_PATH)
add_library(deps::curl INTERFACE IMPORTED)
set_target_properties(
  deps::curl PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "/usr/include/curl"
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
  add_library(deps::sand INTERFACE IMPORTED) # no-op target for unconditional
                                             # linking
  message(STATUS "SANDReco: disabled (SANDRECO_INC/LIB not set)")
endif()

add_subdirectory(src)
