#include "TMSRecoBranchFiller.h"

namespace cafmaker
{
  TMSRecoBranchFiller::TMSRecoBranchFiller(const std::string &tmsRecoFilename)
  {
    fTMSRecoFile = new TFile(tmsRecoFilename.c_str(), "READ");
    name = std::string("TMS");
    if (!fTMSRecoFile->IsZombie()) {
      SetConfigured(true);
      // Save pointer to input tree
      TMSRecoTree = dynamic_cast<TTree*>(fTMSRecoFile->Get("Line_Candidates"));
      if (!TMSRecoTree) {
        std::cerr << "Did not find TMS reco tree Line_Candidates in input file " << tmsRecoFilename << std::endl;
        std::cerr << "Are you sure this is a TMS reco file?" << std::endl;
        throw;
      }
      TMSRecoTree->SetBranchAddress("nLines",               &_nLines);
      TMSRecoTree->SetBranchAddress("nHitsInTrack",          _nHitsInTrack);
      TMSRecoTree->SetBranchAddress("TrackLength",           _TrackLength);
      TMSRecoTree->SetBranchAddress("TotalTrackEnergy",      _TotalTrackEnergy);
      TMSRecoTree->SetBranchAddress("Occupancy",             _Occupancy);

      TMSRecoTree->SetBranchAddress("TrackHitPos",           _TrackHitPos);
      TMSRecoTree->SetBranchAddress("DirectionX_Upstream",   _DirectionX_Upstream);
      TMSRecoTree->SetBranchAddress("DirectionZ_Upstream",   _DirectionZ_Upstream);
      TMSRecoTree->SetBranchAddress("DirectionX_Downstream", _DirectionX_Downstream);
      TMSRecoTree->SetBranchAddress("DirectionZ_Downstream", _DirectionZ_Downstream);
    } else {
      fTMSRecoFile = NULL;
      TMSRecoTree  = NULL;
      std::cerr << "The TMS reco file you provided: " << tmsRecoFilename 
                << " appears to be a Zombie 🧟<" << std::endl;
      throw;
    }
  }

  // ---------------------------------------------------------------------------

  // here we copy all the TMS reco into the SRTMS branch of the StandardRecord object.
  void TMSRecoBranchFiller::_FillRecoBranches(std::size_t evtIdx,
                                              caf::StandardRecord &sr,
                                              const cafmaker::Params &par) const
  {
    // Get nth entry from tree
    TMSRecoTree->GetEntry(evtIdx);

    // First set number of tracks
    sr.nd.tms.nixn = _nLines;
    sr.nd.tms.ixn.resize(sr.nd.tms.nixn);

    //sr.nd.tms.nixn = _nLines;
    //sr.SRTMSInt.ntracks= _nLines;

    // Fill in the track info 
    for (int i = 0; i < _nLines; ++i) {
      /* Currently we only really care about catching (anti-)muons from LAr,
       * so we assume that each track is probably a(n) (anti-)muon originating
       * from a single interaction upstream in LAr.
       *     Liam     */

      sr.nd.tms.ixn[i].ntracks = 1; // One reco track per interaction (ixn)
      sr.nd.tms.ixn[i].tracks.resize(sr.nd.tms.ixn[i].ntracks);

      double prevz = -9E10;
      // Hit info
      // Loop might be unnecessary if you really want optimisation...
      for (int j = 0; j < _nHitsInTrack[i]; ++j) {
        double z = _TrackHitPos[i][j][0];

        // Should all be ordered in z, check this
        // TODO: Clarence says something further down the chain needs z ordering, investigate
        if (z < prevz - 81) { // TODO: There was previously no -81 here, but allow one layer missmatch. Trackfinding a bit wonky?
          // TODO: Flag these tracks somehow?
          std::cerr << "hits in z not ordered " << (prevz - z) << std::endl;
          throw;
        }
        prevz = z;
      }

      // Save first and last hit in track
      // TMS Reco info is saved in mm whereas CAFs use CM as default -> do conversion here
      sr.nd.tms.ixn[i].tracks[0].start = caf::SRVector3D(_TrackHitPos[i][0][1]/10., -999, _TrackHitPos[i][0][0]/10.);
      sr.nd.tms.ixn[i].tracks[0].end   = caf::SRVector3D(_TrackHitPos[i][_nHitsInTrack[i]-1][1]/10., -999, _TrackHitPos[i][_nHitsInTrack[i]-1][0]/10.);

      // Track info
      sr.nd.tms.ixn[i].tracks[0].len_gcm2  = _TrackLength[i]/10.;
      sr.nd.tms.ixn[i].tracks[0].qual      = _Occupancy[i];
      sr.nd.tms.ixn[i].tracks[0].E         = _TotalTrackEnergy[i];

      // Get the directions
      sr.nd.tms.ixn[i].tracks[0].dir     = caf::SRVector3D(_DirectionX_Upstream[i], -999, _DirectionZ_Upstream[i]);
      sr.nd.tms.ixn[i].tracks[0].enddir  = caf::SRVector3D(_DirectionX_Downstream[i], -999, _DirectionZ_Downstream[i]);
    }
  }

} // end namespace
