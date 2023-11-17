#ifndef CAF_h
#define CAF_h

#include "TFile.h"
#include "TTree.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/Flat/FwdDeclare.h"

// fixme: this is a do-nothing replacement for nusystematics stuff until it's re-enabled
//#include "nusystematics/artless/response_helper.hh"
namespace nusyst
{
  typedef std::string response_helper;
}

class CAF {

public:
  CAF(const std::string &filename, const std::string &rw_fhicl_filename, bool makeFlatCAF, bool storeGENIE);
  ~CAF() = default;
  void fill();
  void fillPOT();
  void write();
  void Print();
  void setToBS();

  // Make ntuple variables public so they can be set from other file
  caf::StandardRecord sr;
  caf::SRGlobal srglobal;

  // Event-by-event geometric efficiency throw results
  std::vector< std::vector < std::vector < uint64_t > > > * geoEffThrowResults;

  // meta
  double pot;
  int meta_run, meta_subrun;
  int version;

  TFile * cafFile;
  TTree * cafSR;
  TTree * cafSRGlobal;
  TTree * cafMVA;
  TTree * cafPOT;

  // store the GENIE record as a branch, if requested
  genie::NtpMCEventRecord * mcrec = nullptr;
  TTree                   * genie = nullptr;

  /// Callback function that can be used to store a GENIE event in the GENIE tree,
  /// if client can't manipulate the `genie` tree above directly
  ///
  /// \param   evtIn  Memory location where the GENIE record to be copied is
  /// \return         Index in the TTree of the output where the record was copied
  int StoreGENIEEvent(const genie::NtpMCEventRecord *evtIn);

  TFile * flatCAFFile                             = nullptr;
  TTree * flatCAFTree                             = nullptr;
  flat::Flat<caf::StandardRecord>* flatCAFRecord  = nullptr;

  nusyst::response_helper rh;
};

#endif
