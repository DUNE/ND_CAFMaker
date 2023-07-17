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
  CAF(const std::string &filename, const std::string &rw_fhicl_filename, bool makeFlatCAF);
  ~CAF() = default;
  void fill();
  void fillPOT();
  void write();
  void Print();
  void setToBS();

  // Make ntuple variables public so they can be set from other file
  caf::StandardRecord sr;
  caf::SRGlobal srglobal;

  // store the GENIE record as a branch
  genie::NtpMCEventRecord * mcrec;

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
  TTree * genie;

  TFile * flatCAFFile                             = nullptr;
  TTree * flatCAFTree                             = nullptr;
  flat::Flat<caf::StandardRecord>* flatCAFRecord  = nullptr;

  nusyst::response_helper rh;
};

#endif
