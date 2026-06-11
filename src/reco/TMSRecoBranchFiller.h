/// Fill TMS reco branches from TMS reco output.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

// Adapted from MINERvA version by Liam O'Sullivan <liam.osullivan@uni-mainz.de>

#ifndef ND_CAFMAKER_TMSRECOBRANCHFILLER_H
#define ND_CAFMAKER_TMSRECOBRANCHFILLER_H

#include <cstddef>
#include <deque>
#include <iostream>
#include <map>
#include <vector>

// The virtual base class
#include "IRecoBranchFiller.h"

// File handlers from ROOT
#include "TFile.h"
#include "TTree.h"

// The duneanaobj includes
#include "duneanaobj/StandardRecord/StandardRecord.h"

namespace cafmaker
{
  class TMSRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      TMSRecoBranchFiller(const std::string & tmsRecoFilename);

      std::deque<Trigger> GetTriggers(int triggerType, bool beamOnly) const override;

      RecoFillerType FillerType() const override { return RecoFillerType::BaseReco; }

      ~TMSRecoBranchFiller();

      struct TimingSummary
      {
        double totalMs = 0.0;
        std::size_t calls = 0;
      };

      struct TimingStats
      {
        bool enabled = false;
        TimingSummary loadTruthSpillEntry;
        TimingSummary buildTruthSpillEntryMap;
        TimingSummary resolveTrueInteractionIDFromVertexIndex;
        TimingSummary resolveRecoTrackTruthParticleIndex;
        TimingSummary findTruthSpillParticleIndex;
        TimingSummary resolvePrimaryTruthParticleIndex;
        TimingSummary resolveRecoTrackInteractionID;
        TimingSummary getTrueInteraction;
        TimingSummary getTrueParticle;
        TimingSummary findSRTrueInteractionIndex;
        TimingSummary findSRTrueParticleIndex;
        TimingSummary fillRecoBranches;
        TimingSummary fillInteractions;
        TimingSummary getTriggers;

        std::size_t buildTruthSpillEntryMapEntries = 0;
        std::size_t loadTruthSpillEntryMisses = 0;
        std::size_t resolvePrimaryParentSteps = 0;
        std::size_t resolvePrimaryLoopFallbacks = 0;
        std::size_t resolvePrimaryMissingParentFallbacks = 0;
        std::size_t findTruthSpillParticleIndexEntriesScanned = 0;
        std::size_t findTruthSpillParticleIndexHits = 0;
        std::size_t resolveRecoTrackInteractionVertexScans = 0;
        std::size_t resolveRecoTrackInteractionExactMatches = 0;
        std::size_t resolveRecoTrackInteractionNearestFallbacks = 0;
        std::size_t resolveRecoTrackInteractionMultiMatches = 0;
        std::size_t fillRecoBranchesTracksProcessed = 0;
        std::size_t fillRecoBranchesSpillsProcessed = 0;
        std::size_t getTriggersEntriesScanned = 0;
        std::size_t getTriggersCreated = 0;
      };

    private:
      void LoadTruthSpillEntry(int spillNo) const;
      unsigned long ResolveTrueInteractionIDFromVertexIndex(const TruthMatcher * truthMatch, int trueVtxIdx) const;
      int ResolveRecoTrackTruthParticleIndex(int recoTrackIdx) const;
      int FindTruthSpillParticleIndex(int vertexId, int trackId) const;
      int ResolvePrimaryTruthParticleIndex(int particleIdx, int recoTrackIdx) const;
      unsigned long ResolveRecoTrackInteractionID(const TruthMatcher * truthMatch, int recoTrackIdx) const;
      void FillInteractions(const TruthMatcher * truthMatch, caf::StandardRecord &sr) const;

      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

      TFile *fTMSRecoFile;
      TTree *TMSRecoTree;
      TTree *TMSTrueTree;
      TTree *TMSTrueSpill;
      TTree *TMSLCTree;

      static constexpr int kTMSMaxTracks = 100;
      static constexpr int kTMSMaxLineHits = 200;
      static constexpr int kTMSMaxTrueParticles = 20000;
      static constexpr int kTMSMaxTrueVertices = 5000;

      // Save the branches that we're reading in
      int   _RunNo;                      ///< Run Number
      int   _nLines;                     ///< Number of Hough Lines reconstructed
      int   _EventNo;                    ///< Event Number
      int   _SliceNo;                    ///< (Time) Slide Number
      int   _SpillNo;                    ///< Spill Number
      int   _nTracks;                    ///< Number of tracks in Interaction
      int   _nHitsInTrack[kTMSMaxTracks];           ///< Number of hits in reco. track
      int   _TrackCharge[kTMSMaxTracks];            ///< Reconstructed charge of track
      float _TrackLength[kTMSMaxTracks];            ///< Length [cm] of the reco. track
      float _TrackArealDensity[kTMSMaxTracks];      ///< Areal Density [g/cm^2] traversed by the reco. track
      float _TrackMomentum[kTMSMaxTracks];          ///< Reco. momentum of the track [MeV]
      float _TrackTotalEnergy[kTMSMaxTracks];       ///< Total reco. track energy [MeV]
      float _TrackEnergyDeposit[kTMSMaxTracks];     ///< Visible energy of reco. track [MeV]
      float _TrackStartPos[kTMSMaxTracks][4];       ///< Reco. start position of the track (x,y,z,t)
      float _TrackEndPos[kTMSMaxTracks][4];         ///< Reco. end position of the track (x,y,z,t)
      float _TrackStartDirection[kTMSMaxTracks][3]; ///< Reco. track direction vector at start (x,y,z)
      float _TrackEndDirection[kTMSMaxTracks][3];   ///< Reco. track direction vector at end (x,y,z)
      float _Occupancy[kTMSMaxTracks];              ///< Fraction of true energy deposits included in the reco. track

      float _TMSStartTime[kTMSMaxTracks];
      float _TrackTime[kTMSMaxTracks];

      float _DirectionX_Downstream[kTMSMaxTracks];
      float _DirectionZ_Downstream[kTMSMaxTracks];
      float _DirectionX_Upstream[kTMSMaxTracks];
      float _DirectionZ_Upstream[kTMSMaxTracks];

      float _TrackHitPos[kTMSMaxTracks][kTMSMaxLineHits][3];     ///< Kalman filtered reco. track hit positions (x,y,z)
      float _TrackRecoHitPos[kTMSMaxTracks][kTMSMaxLineHits][4]; ///< Reco. track hit positions (x,y,z,t)

      // Truth_Info branches used for reco-track -> truth lookup
      int _RecoTruePartId[kTMSMaxTracks];    ///< RecoTrackPrimaryParticleIndex
      int _RecoTruePartIdSec[kTMSMaxTracks]; ///< RecoTrackSecondaryParticleIndex

      // Truth_Spill branches used as the authoritative truth representation
      mutable int _TruthSpillSpillNo;
      mutable int _TruthSpillRunNo;
      mutable int _TruthSpillNTrueParticles;
      mutable int _TruthSpillTrueVtxN;
      mutable int _TruthSpillParticleVertexID[kTMSMaxTrueParticles];
      mutable int _TruthSpillParent[kTMSMaxTrueParticles];
      mutable int _TruthSpillTrackID[kTMSMaxTrueParticles];
      mutable float _TruthSpillBirthPosition[kTMSMaxTrueParticles][4];
      mutable int _TruthSpillTrueVtxID[kTMSMaxTrueVertices];
      mutable float _TruthSpillTrueVtxX[kTMSMaxTrueVertices];
      mutable float _TruthSpillTrueVtxY[kTMSMaxTrueVertices];
      mutable float _TruthSpillTrueVtxZ[kTMSMaxTrueVertices];

      bool is_data;
      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;    ///< the last trigger requested using _FillRecoBranches()
      mutable std::map<int, Long64_t> fTruthSpillEntryBySpillNo;
      mutable TimingStats fTiming;

  };

}
#endif //ND_CAFMAKER_TMSRECOBRANCHFILLER_H
