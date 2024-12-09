/// Fill SAND reco branches using SAND reco data
///
/// \author  L. Di Noto, reworked by M. Vicenzi
/// \date    Apr. 2022
///

#ifndef ND_CAFMAKER_SANDRECOBRANCHFILLER_H
#define ND_CAFMAKER_SANDRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

class TFile;
class TTree;

struct event;

namespace cafmaker
{
  class SANDRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      SANDRecoBranchFiller(const std::string &SANDRecoFilename);
      
      std::deque<Trigger> GetTriggers(int triggerType) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::BaseReco; }


    private:
      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;


      void FillECalClusters(const TruthMatcher * truthMatch,
                            caf::StandardRecord &sr) const;

      void FillInteractions(const TruthMatcher * truthMatch,
                            caf::StandardRecord &sr) const;

      TFile* fSANDRecoFile;
      TTree* NDSANDRecoTree;
      TTree* NDSANDEventTree;
      
      //struct event* fEvent;

      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;    ///< the last trigger requested using _FillRecoBranches()
  
  };

}

#endif //ND_CAFMAKER_SANDRECOBRANCHFILLER_H
