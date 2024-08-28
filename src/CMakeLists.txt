message("ND_CAFMaker in src/")

# Suck up all of the .cc files for the source.  This isn't good CMAKE
# practice, but it makes maintaining this file easier.
# file(GLOB_RECURSE source *.cc)
# file(GLOB_RECURSE includes NDCAFMaker*.hh)
file(GLOB cafmaker_source makeCAF.C)
file(GLOB cafmaker_includes *.h)

set(source ${cafmaker_source})
set(includes ${cafmaker_includes})

#set(GENIE_INCLUDES $ENV{GENIE_INC})
include_directories(${GENIE_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/reco ${CMAKE_CURRENT_SOURCE_DIR}/truth ${CMAKE_CURRENT_SOURCE_DIR}/util)

# Where the macro files can be found.
#add_definitions(-DDETSIM_INSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\")

add_executable(makeCAF makeCAF.C)
 
target_link_libraries(makeCAF LINK_PUBLIC cafmaker)
install(TARGETS makeCAF RUNTIME DESTINATION bin)

# Build the library.
add_library(cafmaker SHARED ${source})

target_include_directories(cafmaker PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/NDCAFMaker>")

target_link_libraries(cafmaker PUBLIC
  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# Install the G4 macro files used to control the MC.
#install(FILES cafmaker-defaults-1.0.mac DESTINATION lib/NDCAFMaker)

# Install the library for edep-sim
install(TARGETS cafmaker
  EXPORT NDCAFMakerTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/NDCAFMaker )

# Install the header files.
install(FILES ${includes} DESTINATION include/NDCAFMaker)
