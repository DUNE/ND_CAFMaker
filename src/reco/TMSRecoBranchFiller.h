/// Fill TMS reco branches from TMS reco output.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

// Adapted from MINERvA version by Liam O'Sullivan <liam.osullivan@uni-mainz.de>

#ifndef ND_CAFMAKER_TMSRECOBRANCHFILLER_H
#define ND_CAFMAKER_TMSRECOBRANCHFILLER_H

#include <iostream>
#include <deque>

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
      void _FillRecoBranches(const Trigger &trigger,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par,
                             const TruthMatcher *truthMatcher) const override;

      TFile *fTMSRecoFile;
      TTree *TMSRecoTree;

      // Save the branches that we're reading in
      int _nLines;
      int _EventNo;
      int _SliceNo;
      int _SpillNo;

      int _nTracks;
      int _nHitsInTrack[10];
      float _TrackLength[10];
      int _TrackCharge[10];
      float _TrackMomentum[10];
      float _TrackTotalEnergy[10];
      float _TrackEnergyDeposit[10];
      float _TrackStartPos[10][3];
      float _TrackEndPos[10][3];
      float _TrackStartDirection[10][3];
      float _TrackEndDirection[10][3];
      float _Occupancy[10];

      float _DirectionX_Downstream[10];
      float _DirectionZ_Downstream[10];
      float _DirectionX_Upstream[10];
      float _DirectionZ_Upstream[10];

      // [100][200][4] needs to match TMS reco output (check file if in doubt)
      float _TrackHitPos[100][200][4];
      float _TrackRecoHitPos[100][200][4];

      bool is_data;
      mutable std::vector<cafmaker::Trigger> fTriggers;
      mutable decltype(fTriggers)::const_iterator  fLastTriggerReqd;    ///< the last trigger requested using _FillRecoBranches()

  };

}
#endif //ND_CAFMAKER_TMSRECOBRANCHFILLER_H
