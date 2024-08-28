cmake_minimum_required(VERSION 3.10) # Probably works with older CMake, but idk?

# Set the project name and version.  The project name will be used to
# import this into other CMAKE files.
project(NDCAFMaker VERSION 4.6.4)

message("DUNE Near-Detector CAF Maker -- ${VERSION}")

# Define the options that can be set in the cache, or on the cmake
# command line.
set(CMAKE_BUILD_TYPE Release)

# Check to see if this is MACOS, not that we try to deal with that. Y'all on your own
if(APPLE)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif (APPLE)

# Make sure that ROOT is available.  ROOT is absolutely required.
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
find_package(Geant4 10.0 REQUIRED)
find_package(duneanaobj REQUIRED)
#find_package(ifdhc REQUIRED)
find_package(hdf5 REQUIRED)
#find_package(hf5cpp REQUIRED)
#find_package(fhicl-cpp REQUIRED)
#find_package(tbb REQUIRED)
#find_package(gsl REQUIRED)
#find_package(pythia6 REQUIRED)
find_package(libxml2 REQUIRED)
#find_package(dk2nugenie REQUIRED)
find_package(log4cpp REQUIRED)

# Look for GEANT4.  If it's found, then build the simulation.
# Otherwise, this can be skipped, but only the i/o library is built.

add_subdirectory(src)

# Install some extra files.  These are in subdirectories, but are not
# really part of the package.  They are still useful for users to have
# available as examples.

install(DIRECTORY cfg scripts DESTINATION share/NDCAFMaker
  FILES_MATCHING
  PATTERN README*
  PATTERN *.fcl
  PATTERN *.py)

# install(DIRECTORY test DESTINATION share/NDCAFMaker
#   FILES_MATCHING
#   PATTERN README*
#   PATTERN *.sh
#   PATTERN *.C
#   PATTERN *.py)

# install(DIRECTORY tools DESTINATION share/NDCAFMaker
#   FILES_MATCHING
#   PATTERN README*
#   PATTERN *.sh
#   PATTERN *.C
#   PATTERN *.py)

#############################################################
#
# Prepare the package so that it can be used with the find_package interface.
#
#############################################################

# Include module with function 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Build the targets description so that the package can be configured
# using find_package.
install(EXPORT NDCAFMakerTargets
  NAMESPACE NDCAFMaker::
  DESTINATION lib/cmake/NDCAFMaker)

# Write the 'NDCAFMakerConfigVersion.cmake' file which can be used to
# check if a version meets the requested properties.
write_basic_package_version_file(
  NDCAFMakerConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)

# Write the 'NDCAFMakerConfig.cmake' file so that a user package can
# access this with find_package.
configure_package_config_file(
  PackageConfig.cmake.in
  NDCAFMakerConfig.cmake
  PATH_VARS CMAKE_INSTALL_PREFIX
  INSTALL_DESTINATION lib/cmake/NDCAFMaker)

# Install the config files.
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/NDCAFMakerConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/NDCAFMakerConfigVersion.cmake
  DESTINATION lib/cmake/NDCAFMaker )
