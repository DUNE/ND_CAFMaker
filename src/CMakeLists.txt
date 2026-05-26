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
    reco/SANDRecoBranchFiller.cxx
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

target_include_directories(ND_CAFMaker PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

if(SAND_ENABLED)
  target_compile_definitions(ND_CAFMaker PRIVATE ENABLE_SAND)
endif()

target_link_libraries(
  ND_CAFMaker
  PUBLIC deps::log4cpp
         deps::libxml2
         deps::hdf5
         deps::pythia6
         deps::gsl
         deps::lhapdf
         deps::genie
         deps::nlohmann
         deps::srproxy
         deps::curl
         ROOT::Core
         ROOT::MathMore
         ROOT::Geom
         ROOT::EGPythia6
         ROOT::GenVector
         Boost::program_options
         fhiclcpp::fhiclcpp
         fhiclcpp::types
         EDepSim::edepsim_io
         duneanaobj::StandardRecord
         duneanaobj::StandardRecordFlat)

# SANDReco is optional — link only when SAND_ENABLED at configure time.
if(SAND_ENABLED)
  target_link_libraries(ND_CAFMaker PUBLIC deps::sand)
endif()

set_source_files_properties(makeCAF.C PROPERTIES LANGUAGE CXX)
add_executable(makeCAF makeCAF.C)
target_include_directories(makeCAF PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(makeCAF PRIVATE ND_CAFMaker)

# WARNING: testHDF is currently broken. testHDF.C includes
# reco/NDLArProductFiller.h and reco/NDLArSummaryH5DatasetReader.h, which no
# longer exist in the tree, so enabling ENABLE_TESTEXE will fail at compile
# time. Kept behind the flag (rather than deleted) until someone decides
# whether to port the test to the current H5 readers or drop it entirely.
# if(ENABLE_TESTEXE)
#   set_source_files_properties(testHDF.C PROPERTIES LANGUAGE CXX)
#   add_executable(testHDF testHDF.C)
#   target_include_directories(testHDF PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
#   target_link_libraries(testHDF PRIVATE ND_CAFMaker)
# endif()

install(TARGETS ND_CAFMaker LIBRARY DESTINATION lib)

install(TARGETS makeCAF RUNTIME DESTINATION bin)

if(ENABLE_TESTEXE)
  install(TARGETS testHDF RUNTIME DESTINATION bin)
endif()
