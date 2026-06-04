/// Fill TMS reco branches from TMS reco output.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

// Adapted from MINERvA version by Liam O'Sullivan <liam.osullivan@uni-mainz.de>

#ifndef ND_CAFMAKER_TMSRECOBRANCHFILLER_H
#define ND_CAFMAKER_TMSRECOBRANCHFILLER_H

#include <iostream>
#include <deque>
#include <map>

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

    private:
      void LoadTruthSpillEntry(int spillNo) const;
      unsigned long ResolveTrueInteractionIDFromVertexIndex(const TruthMatcher * truthMatch, int trueVtxIdx) const;
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

      // Save the branches that we're reading in
      int   _RunNo;                      ///< Run Number
      int   _nLines;                     ///< Number of Hough Lines reconstructed
      int   _EventNo;                    ///< Event Number
      int   _SliceNo;                    ///< (Time) Slide Number
      int   _SpillNo;                    ///< Spill Number
      int   _nTracks;                    ///< Number of tracks in Interaction
      int   _nHitsInTrack[10];           ///< Numebr of Hits in reco. track
      int   _TrackCharge[10];            ///< Reconstructed Charge of track
      float _TrackLength[10];            ///< Length [cm] of the reco. track
      float _TrackArealDensity[10];      ///< Areal Density [g/cm^2] traversed by the reco. track
      float _TrackMomentum[10];          ///< Reco. momentum of the track [MeV]
      float _TrackTotalEnergy[10];       ///< Total reco. track energy [MeV]
      float _TrackEnergyDeposit[10];     ///< Visible energy of reco. track [MeV]
      float _TrackStartPos[10][3];       ///< Reco. start position of the track (x,y,z)
      float _TrackEndPos[10][3];         ///< Reco. end position of the track (x,y,z)
      float _TrackStartDirection[10][3]; ///< Reco. track direction vector at start (x,y,z)
      float _TrackEndDirection[10][3];   ///< Reco. track direction vector at start (x,y,z)
      float _Occupancy[10];              ///< Fraction of true energy deposits included in the reco. track

      float _TMSStartTime[10];
      float _TrackTime[10];

      float _DirectionX_Downstream[10];
      float _DirectionZ_Downstream[10];
      float _DirectionX_Upstream[10];
      float _DirectionZ_Upstream[10];

      // [100][200][4] needs to match TMS reco output (check TMS Reco file if in doubt)
      float _TrackHitPos[100][200][4];     ///< Array of reco. track hit positions (x,y,z,t)
      float _TrackRecoHitPos[100][200][4]; ///< Array of Kalman filtered reco. track hit positions (x,y,z,t)

      // Truth_Info branches used for reco-track -> truth lookup
      int _RecoTruePartId[5000];    ///< RecoTrackPrimaryParticleIndex
      int _RecoTruePartIdSec[5000]; ///< Secondary 

      // Truth_Spill branches used as the authoritative truth representation
      mutable int _TruthSpillSpillNo;
      mutable int _TruthSpillRunNo;
      mutable int _TruthSpillNTrueParticles;
      mutable int _TruthSpillTrueVtxN;
      mutable int _TruthSpillParticleVertexID[5000];
      mutable int _TruthSpillParent[5000];
      mutable float _TruthSpillBirthPosition[5000][4];
      mutable int _TruthSpillTrueVtxID[505];
      mutable float _TruthSpillTrueVtxX[505];
      mutable float _TruthSpillTrueVtxY[505];
      mutable float _TruthSpillTrueVtxZ[505];

      bool is_data;
      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;    ///< the last trigger requested using _FillRecoBranches()
      mutable std::map<int, Long64_t> fTruthSpillEntryBySpillNo;

  };

}
#endif //ND_CAFMAKER_TMSRECOBRANCHFILLER_H
