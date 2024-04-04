/// Fill TMS reco branches from TMS reco output.
///
/// \author  J. Wolcott <jwolcott@fnal.gov> & F. Akbar <fakbar@ur.rochester.edu>
/// \date    Nov. 2021

#ifndef ND_CAFMAKER_TMSRECOBRANCHFILLER_H
#define ND_CAFMAKER_TMSRECOBRANCHFILLER_H

#include <iostream>

// The virtual base class
#include "IRecoBranchFiller.h"

// File handlers from ROOT
#include "TFile.h"
#include "TTree.h"

// The duneanaobj includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
//#include "duneanaobj/StandardRecord/SRTrack.h"

namespace cafmaker
{
  class TMSRecoBranchFiller : public cafmaker::IRecoBranchFiller
  {
    public:
      TMSRecoBranchFiller(const std::string & tmsRecoFilename);
      ~TMSRecoBranchFiller() {
        delete TMSRecoTree;
        fTMSRecoFile->Close();
        delete fTMSRecoFile;
        TMSRecoTree = NULL;
        fTMSRecoFile = NULL;
      }

    private:
      void _FillRecoBranches(std::size_t evtIdx,
                             caf::StandardRecord &sr,
                             const cafmaker::Params &par) const override;

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
      float _TrackCharge[10];
      float _TotalTrackEnergy[10];
      float _TrackEnergyDeposit[10];
      float _TrackStartPos[10][3];
      float _TrackEndPos[10][3];
      float _TrackDirection[10][3];
      float _Occupancy[10];

      float _DirectionX_Downstream[10];
      float _DirectionZ_Downstream[10];
      float _DirectionX_Upstream[10];
      float _DirectionZ_Upstream[10];

      // [10][200][2] needs to match TMS reco output (check file if in doubt)
      float _TrackHitPos[10][200][2];
      float _TrackRecoHitPos[10][200][2];
  };

}
#endif //ND_CAFMAKER_TMSRECOBRANCHFILLER_H
