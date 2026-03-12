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
  CAF(const std::string &filename, const std::string &sandFilename, const std::string &rw_fhicl_filename, bool makeFlatCAF, bool storeGENIE);
  ~CAF() = default;
  void fill();
  void fillPOT();
  void write();
  void Print();
  void setToBS();

  // funzione per settare meta da SANDRecoFile
  void SetMetaFromRecoFile(const std::string &filename);

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

  enum class InteractionType {
      kUnknownMode = -100,
      kQE = 1,
      kSingleKaon = 2,
      kDIS = 3,
      kRes = 4,
      kCoh = 5,
      kDiffractive = 6,
      kNuElectronElastic = 7,
      kInvMuonDecay = 8,
      kAMNuGamma = 9,
      kMEC = 10,
      kCohElastic = 11,
      kInverseBetaDecay = 12,
      kGlashowResonance = 13,
      kIMDAnnihilation = 14,
      kPhotonCoh = 15,
      kPhotonRes = 16,
      kDarkMatterElastic = 101,
      kDarkMatterDIS = 102,
      kDarkMatterElectron = 103
  };

  enum class SelectionType {
      kGRAINcontained,
      kECalcontained,
      kTrackercontained
  };

  enum class RecoType {
      kNone,
      kRealReco,
      kSmearing
  };
  InteractionType interaction_type;
  SelectionType selection_type;
  RecoType ecal_reco_type;
  RecoType grain_reco_type;
  RecoType tracker_reco_type;
};

#endif