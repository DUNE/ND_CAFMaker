/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi, S. Lanzi

#ifndef ND_CAFMAKER_SANDRECOBRANCHFILLER_H
#define ND_CAFMAKER_SANDRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class TFile;
class TTree;

namespace caf {
class StandardRecord;
}

namespace cafmaker {

enum class SANDRecoVersion : std::int8_t { legacy, experimental };

template<SANDRecoVersion V>
class SANDRecoBranchFiller;

template<>
class SANDRecoBranchFiller<SANDRecoVersion::experimental> : public IRecoBranchFiller {
 public:
  explicit SANDRecoBranchFiller(const std::string& cafFilename, std::string treeName = "cafTree");
  ~SANDRecoBranchFiller();

  std::deque<Trigger> GetTriggers(int triggerType = -1, bool beamOnly = false) const override;

  RecoFillerType FillerType() const override {
    return RecoFillerType::BaseReco;
  }

 private:
  void _FillRecoBranches(const Trigger&, caf::StandardRecord&, const Params&,
                         const TruthMatcher*) const override;

  std::unique_ptr<TFile> fFile;
  TTree* fCAFTree{nullptr}; // owned by fFile
  caf::StandardRecord* fSR{nullptr};
  std::unordered_map<long int, long int> fEvtToEntry;
  std::vector<Trigger> fTriggers;
};

#ifdef ENABLE_SANDRECO_LEGACY

template<>
class SANDRecoBranchFiller<SANDRecoVersion::legacy> : public IRecoBranchFiller {
 public:
  explicit SANDRecoBranchFiller(const std::string& recoFilename);
  ~SANDRecoBranchFiller();

  std::deque<Trigger> GetTriggers(int triggerType = -1, bool beamOnly = false) const override;

  RecoFillerType FillerType() const override {
    return RecoFillerType::BaseReco;
  }

 private:
  void _FillRecoBranches(const Trigger&, caf::StandardRecord&, const Params&,
                         const TruthMatcher*) const override;

  void FillECalClusters(const TruthMatcher*, caf::StandardRecord&, std::vector<cluster>&) const;
  void FillTracks(const TruthMatcher*, caf::StandardRecord&, std::vector<track>&) const;

  std::unique_ptr<TFile> fFile;
  TTree* fRecoTree{nullptr}; // owned by fFile
  struct event* fEvent{nullptr};
  std::vector<Trigger> fTriggers;
};

#endif // ENABLE_SANDRECO_LEGACY

} // namespace cafmaker

#endif // ND_CAFMAKER_SANDRECOBRANCHFILLER_H
