/// \file NDGArRecoBranchFiller.h
/// Fill NDGAr reco branches using NDGAr reco data
///
/// \author  F. Martinez Lopez <f.martinezlopez@qmul.ac.uk>
/// \date    Oct. 2022
///

#ifndef ND_CAFMAKER_NDGArRECOBRANCHFILLER_H
#define ND_CAFMAKER_NDGArRECOBRANCHFILLER_H

#include "IRecoBranchFiller.h"

#include <vector>
#include <cmath>

class TFile;
class TTree;

namespace caf
{
  class SRTrueInteraction;
  class SRTrueParticle;
}

namespace cafmaker
{
  class NDGArRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      NDGArRecoBranchFiller(const std::string &NDGArRecoFilename);
      std::deque<Trigger> GetTriggers(int triggerType) const override;

    private:
      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

      void FillTrueInteraction(caf::SRTrueInteraction & srTrueInt) const;

      void FillInteractions(const TruthMatcher * truthMatch,
                            caf::StandardRecord &sr) const;
      
      void FillTruth(caf::SRTrueParticle & srTruePart,
                     size_t iTrack) const;
          
      void FillParticles(const TruthMatcher * truthMatch,
                         caf::StandardRecord &sr) const;

      void FillTracks(const TruthMatcher * truthMatch,
                      caf::StandardRecord &sr) const;

      void FillClusters(caf::StandardRecord &sr) const;

      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;    ///< the last trigger requested using _FillRecoBranches()

      TFile* fNDGArRecoFile;
      TTree* NDGArRecoTree;
      
      //int                fRun;
      //int                fSubRun;
      int                fEvent;

      std::vector<int> *      fGPartIdx=0;

      std::vector<int> *      fMCNuPDG=0;
      std::vector<float_t> *  fMCNuPX=0;
      std::vector<float_t> *  fMCNuPY=0;
      std::vector<float_t> *  fMCNuPZ=0;

      std::vector<float_t> *  fMCVertX=0;
      std::vector<float_t> *  fMCVertY=0;
      std::vector<float_t> *  fMCVertZ=0;

      // MC Particles info
      std::vector<int> * fMCTrkID=0;
      std::vector<int> * fMCPDG=0;
      std::vector<int> * fMCMotherIndex=0;
      std::vector<int> * fMCPDGMother=0;

      std::vector<float_t> * fMCStartX=0;
      std::vector<float_t> * fMCStartY=0;
      std::vector<float_t> * fMCStartZ=0;
      std::vector<float_t> * fMCStartPX=0;
      std::vector<float_t> * fMCStartPY=0;
      std::vector<float_t> * fMCStartPZ=0;

      //std::vector<float_t> * fMCEndX=0;
      //std::vector<float_t> * fMCEndY=0;
      //std::vector<float_t> * fMCEndZ=0;
      //std::vector<float_t> * fMCEndPX=0;
      //std::vector<float_t> * fMCEndPY=0;
      //std::vector<float_t> * fMCEndPZ=0;

      // Track-related info
      std::vector<float_t> * fTrackStartX=0;
      std::vector<float_t> * fTrackStartY=0;
      std::vector<float_t> * fTrackStartZ=0;
      std::vector<float_t> * fTrackEndX=0;
      std::vector<float_t> * fTrackEndY=0;
      std::vector<float_t> * fTrackEndZ=0;

      std::vector<float_t> * fTrackStartPX=0;
      std::vector<float_t> * fTrackStartPY=0;
      std::vector<float_t> * fTrackStartPZ=0;
      std::vector<float_t> * fTrackEndPX=0;
      std::vector<float_t> * fTrackEndPY=0;
      std::vector<float_t> * fTrackEndPZ=0;

      std::vector<float_t> * fTrackLenF=0;
      std::vector<float_t> * fTrackLenB=0;
      //std::vector<float_t> * fTrackPF=0;
      //std::vector<float_t> * fTrackPB=0;
      std::vector<float_t> * fTrackAvgIonF=0;
      std::vector<float_t> * fTrackAvgIonB=0;

      std::vector<int> * fTrackIDNumber=0;
      std::vector<int> * fTrackNClusters=0;

      std::vector<int>     * fTrackPIDF=0;
      std::vector<float_t> * fTrackPIDProbF=0;
      std::vector<int>     * fTrackPIDB=0;
      std::vector<float_t> * fTrackPIDProbB=0;

      std::vector<int>     * fTrackMCindex=0;
      std::vector<float_t> * fTrackMCfrac=0;

      // ECAL-related info
      std::vector<float_t> * fECALClusterX=0;
      std::vector<float_t> * fECALClusterY=0;
      std::vector<float_t> * fECALClusterZ=0;

      std::vector<int> * fECALClusterIDNumber=0;

      std::vector<float_t> * fECALClusterEnergy=0;
      std::vector<int>     * fECALClusterNhits=0;

      std::vector<int>     * fECALClusterMCindex=0;
      std::vector<float_t> * fECALClusterMCfrac=0;

      // ECAL-track associations
      std::vector<int> * fECALAssn_ClusterID=0;
      std::vector<int> * fECALAssn_TrackID=0;

      // MuID-related info
      std::vector<float_t> * fMuIDClusterX=0;
      std::vector<float_t> * fMuIDClusterY=0;
      std::vector<float_t> * fMuIDClusterZ=0;

      std::vector<int> * fMuIDClusterIDNumber=0;

      std::vector<float_t> * fMuIDClusterEnergy=0;
      std::vector<int>     * fMuIDClusterNhits=0;

      std::vector<int>     * fMuIDClusterMCindex=0;
      std::vector<float_t> * fMuIDClusterMCfrac=0;

      // MuID-track associations
      //std::vector<int> * fMuIDAssn_ClusterID=0;
      //std::vector<int> * fMuIDAssn_TrackID=0;

      // Reco particle info
      std::vector<float_t> * fRecoParticleMomentum=0;
      std::vector<int> * fRecoParticleNHitsECAL=0;
      std::vector<float_t> * fRecoParticleMuonScore=0;
      std::vector<int> * fRecoParticlePID=0;

      // Reco interaction info
      std::vector<float_t> * fRecoNuEnergy=0;

  };

}

#endif //ND_CAFMAKER_NDGArRECOBRANCHFILLER_H
