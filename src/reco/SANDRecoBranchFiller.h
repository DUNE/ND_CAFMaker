/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi
/// \date    Apr. 2022
///

#ifndef ND_CAFMAKER_SANDRECOBRANCHFILLER_H
#define ND_CAFMAKER_SANDRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"
#include <deque>
#include <string>

#ifdef ENABLE_SAND
#include "struct.h"
#endif

namespace cafmaker
{
  class SANDRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
  public:
    SANDRecoBranchFiller(const std::string &SANDRecoFilename);

    std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

    RecoFillerType FillerType() const override { return RecoFillerType::BaseReco; }

  private:
#ifdef ENABLE_SAND
    void _FillRecoBranches(const Trigger &trigger,
                           caf::StandardRecord &sr,
                           const cafmaker::Params &par,
                           const TruthMatcher *truthMatcher) const override;

    void FillECalClusters(const TruthMatcher * truthMatch,
                          caf::StandardRecord &sr, std::vector<cluster> &cl) const;
    void FillTracks(const TruthMatcher * truthMatch,
                    caf::StandardRecord &sr, std::vector<track> &tr) const;

    TFile* fSANDRecoFile;
    TTree* NDSANDRecoTree;
    TTree* NDSANDEventTree;
    struct event* fEvent;

    mutable std::vector<cafmaker::Trigger> fTriggers;
    mutable decltype(fTriggers)::const_iterator fLastTriggerReqd;

#else
    void _FillRecoBranches(const Trigger &trigger, caf::StandardRecord &sr, const cafmaker::Params &par, const TruthMatcher *truthMatcher) const override;
    // These members exist even when ENABLE_SAND is false
    // because ROOT needs them to generate the dictionary correctly.
    TFile* fSANDRecoFile = nullptr;
    TTree* NDSANDRecoTree = nullptr;
    TTree* NDSANDEventTree = nullptr;
    struct event* fEvent = nullptr;

    mutable std::vector<cafmaker::Trigger> fTriggers;
    mutable decltype(fTriggers)::const_iterator fLastTriggerReqd;
#endif
  };
}

#endif