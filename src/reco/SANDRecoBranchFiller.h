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
#include "duneanaobj/StandardRecord/SRTruthBranch.h"
#include "duneanaobj/StandardRecord/SRTrueInteraction.h"
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

      void FillGRAINTracks(const TruthMatcher * truthMatch,
                            caf::StandardRecord &sr, std::vector<track_grain_lens> &trackgrain) const;

      //new -- FillSANDTracks 
      //Jeremy's suggestion 
      
      void FillSANDTracks(const TruthMatcher * truthMatch,
                            caf::StandardRecord &sr, std::vector<track> &tr, std::vector<track_grain_lens> &trackgrain) const;

    
      //new -- neutrino level info
      void FillRecoNu(caf::StandardRecord &sr, const event* evt, int event_num) const; 
      //new -- SRRecoParticle 
      void FillRecoParticle(caf::StandardRecord &sr, const event* evt, int event_num) const; 
      //new -- SRTrueParticle
      void FillTruePrimaryParticles(caf::StandardRecord &sr, const event* evt, int event_num) const;
      void FillTrueSecondaryParticles(caf::StandardRecord &sr, const event* evt, int event_num) const;
      //new -- true neutrino 
      void FillTrueNu(caf::StandardRecord &sr, const event* evt, int event_num) const;

      TFile* fSANDRecoFile;
      TTree* NDSANDRecoTree;
      TTree* NDSANDEventTree;
      
      //struct event* fEvent;
      //new 
      event* evt = nullptr;

      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;
 #else
       void _FillRecoBranches(const Trigger &trigger, caf::StandardRecord &sr, const cafmaker::Params &par, const TruthMatcher *truthMatcher) const override;
       TFile* fSANDRecoFile;
       TTree* NDSANDRecoTree;
       struct event* fEvent;
       mutable std::vector<cafmaker::Trigger> fTriggers;
       mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;
#endif
  };
}

#endif // ND_CAFMAKER_SANDRECOBRANCHFILLER_H
