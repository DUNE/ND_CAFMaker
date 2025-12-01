/// Fill Pandora LArRecoND branches
///
/// \author  John Back <J.J.Back@warwick.ac.uk>
/// \date    Aug 2024

#ifndef ND_CAFMAKER_PandoraLArRecoNDBranchFiller_H
#define ND_CAFMAKER_PandoraLArRecoNDBranchFiller_H

#include <vector>

// The virtual base class
#include "reco/IRecoBranchFiller.h"
#include "truth/FillTruth.h"

// ROOT headers
#include "TFile.h"
#include "TTree.h"

// duneanaobj
#include "duneanaobj/StandardRecord/StandardRecord.h"

namespace cafmaker
{

  class PandoraLArRecoNDBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      PandoraLArRecoNDBranchFiller(const std::string &pandoraLArRecoNDFilename,
           const float LArDensity = 1.3973);

      std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      bool IsBeamTrigger(int triggerType) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::BaseReco; }

    private:
      void _FillRecoBranches(const Trigger &trigger,
           caf::StandardRecord &sr,
           const cafmaker::Params &par,
           const TruthMatcher *truthMatch = nullptr) const override;

      void FillTracks(caf::StandardRecord &sr, const int nClusters, const std::vector<int> &uniqueSliceIDs,
          std::vector<caf::SRInteraction> &nuInteractions, const TruthMatcher *truthMatch) const;
      void FillShowers(caf::StandardRecord &sr, const int nClusters, const std::vector<int> &uniqueSliceIDs,
           std::vector<caf::SRInteraction> &nuInteractions, const TruthMatcher *truthMatch) const;
      
      std::unique_ptr<TFile> m_LArRecoNDFile;
      std::unique_ptr<TTree> m_LArRecoNDTree;

      int m_eventId;
      int m_run;
      int m_subRun;
      int m_unixTime;
      int m_unixTimeUsec;
      int m_startTime;
      int m_triggerType;
      std::vector<int> *m_isShowerVect = nullptr;
      std::vector<int> *m_sliceIdVect = nullptr;
      std::vector<float> *m_startXVect = nullptr;
      std::vector<float> *m_startYVect = nullptr;
      std::vector<float> *m_startZVect = nullptr;
      std::vector<float> *m_endXVect = nullptr;
      std::vector<float> *m_endYVect = nullptr;
      std::vector<float> *m_endZVect = nullptr;
      std::vector<float> *m_dirXVect = nullptr;
      std::vector<float> *m_dirYVect = nullptr;
      std::vector<float> *m_dirZVect = nullptr;
      std::vector<float> *m_energyVect = nullptr;
      std::vector<int> *m_n3DHitsVect = nullptr;
      std::vector<long> *m_mcNuIdVect = nullptr;
      std::vector<long> *m_mcLocalIdVect = nullptr;
      std::vector<int> *m_isPrimaryVect = nullptr;
      std::vector<float> *m_completenessVect = nullptr;
      std::vector<float> *m_nuVtxXVect = nullptr;
      std::vector<float> *m_nuVtxYVect = nullptr;
      std::vector<float> *m_nuVtxZVect = nullptr;
      std::vector<int> *m_isRecoPrimaryVect = nullptr;
      std::vector<int> *m_recoPDGVect = nullptr;
      // TRACK
      std::vector<float> *m_trackScoreVect = nullptr;
      std::vector<float> *m_trkfitPID_Mu = nullptr;
      std::vector<float> *m_trkfitPID_Pro = nullptr;
      std::vector<int> *m_trkfitPID_NDF = nullptr;
      std::vector<float> *m_trkfitIsContained = nullptr; // TODO : branch not yet in PandoraOuterface
      std::vector<float> *m_trkfitLength = nullptr;
      std::vector<float> *m_trkfitKEFromLengthMuon = nullptr;
      std::vector<float> *m_trkfitKEFromLengthProton = nullptr;
      std::vector<float> *m_trkfitPFromLengthMuon = nullptr;
      std::vector<float> *m_trkfitPFromLengthProton = nullptr;
      std::vector<float> *m_trkfitStartX = nullptr;
      std::vector<float> *m_trkfitStartY = nullptr;
      std::vector<float> *m_trkfitStartZ = nullptr;
      std::vector<float> *m_trkfitEndX = nullptr;
      std::vector<float> *m_trkfitEndY = nullptr;
      std::vector<float> *m_trkfitEndZ = nullptr;
      std::vector<float> *m_trkfitStartDirX = nullptr;
      std::vector<float> *m_trkfitStartDirY = nullptr;
      std::vector<float> *m_trkfitStartDirZ = nullptr;
      // SHOWER
      std::vector<float> *m_shwrfitLength = nullptr;
      std::vector<float> *m_shwrfitCentroidX = nullptr;
      std::vector<float> *m_shwrfitCentroidY = nullptr;
      std::vector<float> *m_shwrfitCentroidZ = nullptr;
      std::vector<float> *m_shwrfitStartX = nullptr;
      std::vector<float> *m_shwrfitStartY = nullptr;
      std::vector<float> *m_shwrfitStartZ = nullptr;
      std::vector<float> *m_shwrfitDirX = nullptr;
      std::vector<float> *m_shwrfitDirY = nullptr;
      std::vector<float> *m_shwrfitDirZ = nullptr;
      std::vector<int> *m_shwrSliceId = nullptr;
      std::vector<int> *m_shwrClusterId = nullptr;
      std::vector<float> *m_startTrkSlidingFitDirX = nullptr;
      std::vector<float> *m_startTrkSlidingFitDirY = nullptr;
      std::vector<float> *m_startTrkSlidingFitDirZ = nullptr;
      std::vector<float> *m_shwrStartPointsX = nullptr;
      std::vector<float> *m_shwrStartPointsY = nullptr;
      std::vector<float> *m_shwrStartPointsZ = nullptr;
      std::vector<float> *m_shwrdEdx = nullptr;
      std::vector<float> *m_shwrTotalE = nullptr;
      std::vector<int> *m_shwrStartPointsRecoId = nullptr;
      std::vector<float> *m_minProjection = nullptr;
      std::vector<float> *m_medianQ = nullptr;
      std::vector<int> *m_chargePerHit = nullptr;
      std::vector<float> *m_chargePerHitStartPoints = nullptr;
      std::vector<int> *m_pitchValue = nullptr;
      std::vector<float> *m_shwrPCAX = nullptr;
      std::vector<float> *m_shwrPCAY = nullptr;
      std::vector<float> *m_shwrPCAZ = nullptr;
      std::vector<float> *m_shwrStartHitPositionX = nullptr;
      std::vector<float> *m_shwrStartHitPositionY = nullptr;
      std::vector<float> *m_shwrStartHitPositionZ = nullptr; 

      mutable std::vector<cafmaker::Trigger> m_Triggers;
      mutable decltype(m_Triggers)::const_iterator  m_LastTriggerReqd; ///< the last trigger requested using _FillRecoBranches
      mutable std::map<int, int> fEntryMap; //Map of the filtered trigger entries stored in the caf file
      const float m_LArDensity;
      const float m_TrackShowerCut = 0.5; ///< threshold on the trackScore variable to decide if a reco particle is track(>=0.5) or shower(<0.5)
  };

}
#endif //ND_CAFMAKER_PandoraLArRecoNDBranchFiller_H
