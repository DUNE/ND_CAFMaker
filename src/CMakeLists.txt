# Shared library libND_CAFMaker.so
set(LIB_SOURCES
    CAF.cxx
    beam/IFBeam.cxx
    reco/DLP_h5_classes.cxx
    reco/MINERvARecoBranchFiller.cxx
    reco/MLNDLArRecoBranchFiller.cxx
    reco/NDLArDLPH5DatasetReader.cxx
    reco/NDLArMINERvAMatchRecoFiller.cxx
    reco/NDLArTMSMatchRecoFiller.cxx
    reco/NDLArTMSUniqueMatchRecoFiller.cxx
    reco/PandoraLArRecoNDBranchFiller.cxx
    reco/SANDRecoBranchFiller.cxx # always compiled; ENABLE_SAND gates behaviout
    reco/TMSRecoBranchFiller.cxx
    reco/readH5/DatasetBuffer.cxx
    reco/readH5/H5DataView.cxx
    reco/readH5/IH5Viewer.cxx
    truth/FillTruth.cxx
    util/FloatMath.cxx
    util/GENIEBannerBypass.cxx
    util/GENIEQuiet.cxx
    util/IFBeamUtils.cxx
    util/Loggable.cxx
    util/Logger.cxx
    util/Progress.cxx)

add_library(ND_CAFMaker SHARED ${LIB_SOURCES})

# Expose src/ itself so #include "CAF.h" etc. work for consumers
target_include_directories(ND_CAFMaker PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# Enable SAND code paths only when the library was detected at configure time
if(SAND_ENABLED)
  target_compile_definitions(ND_CAFMaker PRIVATE ENABLE_SAND)
endif()

target_link_libraries(
  ND_CAFMaker
  PUBLIC deps::log4cpp
         deps::tbb
         deps::libxml2
         deps::hdf5
         deps::pythia6
         deps::duneanaobj
         deps::gsl
         deps::lhapdf
         deps::genie
         deps::root_glibs
         deps::boost_po
         deps::cetlib
         deps::fhiclcpp
         deps::edepsim
         deps::nlohmann
         deps::srproxy
         deps::curl
         deps::sand)

# Executable makeCAF
set_source_files_properties(makeCAF.C PROPERTIES LANGUAGE CXX)
add_executable(makeCAF makeCAF.C)
target_include_directories(makeCAF PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(makeCAF PRIVATE ND_CAFMaker)

# Executable testHDF (optional, off by default)
if(ENABLE_TESTEXE)
  set_source_files_properties(testHDF.C PROPERTIES LANGUAGE CXX)
  add_executable(testHDF testHDF.C)
  target_include_directories(testHDF PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(testHDF PRIVATE ND_CAFMaker)
endif()

# Install rules
install(TARGETS ND_CAFMaker LIBRARY DESTINATION lib)

install(TARGETS makeCAF RUNTIME DESTINATION bin)

if(ENABLE_TESTEXE)
  install(TARGETS testHDF RUNTIME DESTINATION bin)
endif()
