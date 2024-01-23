#include "TMSRecoBranchFiller.h"
#include "truth/FillTruth.h"


namespace cafmaker
{
  TMSRecoBranchFiller::TMSRecoBranchFiller(const std::string &tmsRecoFilename)
    : IRecoBranchFiller("TMS")
  {
    fTMSRecoFile = new TFile(tmsRecoFilename.c_str(), "READ");
    name = std::string("TMS");

    if (!fTMSRecoFile->IsZombie()) {
      SetConfigured(true);

      // Save pointer to input tree
      TMSRecoTree = dynamic_cast<TTree*>(fTMSRecoFile->Get("Reco_Tree"));
      if (!TMSRecoTree) {
        std::cerr << "Did not find TMS reco tree Line_Candidates in input file " << tmsRecoFilename << std::endl;
        std::cerr << "Are you sure this is a TMS reco file?" << std::endl;
        throw;
      }

      TMSRecoTree->SetBranchAddress("EventNo",               &_EventNo);
      TMSRecoTree->SetBranchAddress("SliceNo",               &_SliceNo);
      TMSRecoTree->SetBranchAddress("SpillNo",               &_SpillNo);
      TMSRecoTree->SetBranchAddress("nTracks",               &_nTracks);
      TMSRecoTree->SetBranchAddress("nHits",                 _nHitsInTrack);
      TMSRecoTree->SetBranchAddress("Length",                _TrackLength);
      TMSRecoTree->SetBranchAddress("kCharge",               _TrackCharge);
      TMSRecoTree->SetBranchAddress("Energy",                _TrackEnergy);
      TMSRecoTree->SetBranchAddress("EnergyDeposit",         _TrackEnergyDeposit);
      TMSRecoTree->SetBranchAddress("Occupancy",             _Occupancy);

      //TMSRecoTree->SetBranchAddress("HitPos",           _TrackHitPos); // TODO: how get hits from da vectur??
      //TMSRecoTree->SetBranchAddress("RecoHitPos",       _TrackRecoHitPos);
      TMSRecoTree->SetBranchAddress("StartPos",              _TrackStartPos);
      TMSRecoTree->SetBranchAddress("EndPos",                _TrackEndPos);
      TMSRecoTree->SetBranchAddress("Direction",             _TrackDirection);
    } else {
      fTMSRecoFile = NULL;
      TMSRecoTree  = NULL;
      std::cerr << "The TMS reco file you provided: " << tmsRecoFilename 
                << " appears to be a Zombie 🧟" << std::endl;
      throw;
    }
  }

  // ---------------------------------------------------------------------------

  // here we copy all the TMS reco into the SRTMS branch of the StandardRecord object.
  void TMSRecoBranchFiller::_FillRecoBranches(const Trigger &trigger,
                                              caf::StandardRecord &sr,
                                              const cafmaker::Params &par,
                                              const TruthMatcher *truthMatcher) const
  {
#ifndef DISABLE_TMS
    // Get nth entry from tree
    TMSRecoTree->GetEntry(evtIdx);

    // First set number of tracks
    sr.nd.tms.nixn = _nTracks;
    sr.nd.tms.ixn.resize(sr.nd.tms.nixn);

    //sr.nd.tms.nixn = _nLines;
    //sr.SRTMSInt.ntracks= _nLines;

    // Fill in the track info 
    for (int i = 0; i < _nTracks; ++i) {
      /* Currently we only really care about catching (anti-)muons from LAr,
       * so we assume that each track is probably a(n) (anti-)muon originating
       * from a single interaction upstream in LAr.
       *     Liam     */

      sr.nd.tms.ixn[i].ntracks = 1; // One reco track per interaction (ixn)
      sr.nd.tms.ixn[i].tracks.resize(sr.nd.tms.ixn[i].ntracks);

      // Save first and last hit in track
      // TMS Reco info is saved in mm whereas CAFs use CM as default -> do conversion here
      sr.nd.tms.ixn[i].tracks[0].start = caf::SRVector3D(_TrackStartPos[i][0]/10., _TrackStartPos[i][1]/10., _TrackStartPos[i][1]/10.);
      sr.nd.tms.ixn[i].tracks[0].end   = caf::SRVector3D(_TrackEndPos[i][0]/10., _TrackEndPos[i][1]/10., _TrackEndPos[i][1]/10.);

      // Track info
      sr.nd.tms.ixn[i].tracks[0].len_gcm2  = _TrackLength[i]/10.;
      sr.nd.tms.ixn[i].tracks[0].qual      = _Occupancy[i];
      sr.nd.tms.ixn[i].tracks[0].E         = _TrackEnergy[i];

      // Get the directions
      // TODO: At present tracks are completely straight objects, so dir is the same here for both
      sr.nd.tms.ixn[i].tracks[0].dir     = caf::SRVector3D(_TrackDirection[i][0], _TrackDirection[i][1] , _TrackDirection[i][2]);
      sr.nd.tms.ixn[i].tracks[0].enddir  = caf::SRVector3D(_TrackDirection[i][0], _TrackDirection[i][1] , _TrackDirection[i][2]);
    }
  }

#endif  // DISABLE_TMS

} // end namespace
