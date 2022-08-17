#ifndef CAF_h
#define CAF_h

#include "TFile.h"
#include "TTree.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"

// fixme: this is a do-nothing replacement for nusystematics stuff until it's re-enabled
//#include "nusystematics/artless/response_helper.hh"
namespace nusyst
{
  typedef std::string response_helper;
}

class CAF {

public:
  CAF(  const std::string& filename, const std::string& rw_fhicl_filename );
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

  nusyst::response_helper rh;
};

#endif

